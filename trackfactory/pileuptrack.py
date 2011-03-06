'''
Created on Sep 29, 2010

@author: mkiyer
'''
import re
import logging
import collections
import operator
import tables
import numpy as np

import pysam

from sequel2.io.cbedgraph import array_to_bedgraph
from sequel2.io.sam import bam_pileup_to_array, bam_pileup_to_array_multihit, bam_to_array, bam_to_array_chunked
from sequel2.track.base import parse_interval
from sequel2.track.arraytrack import ArrayTrack

PILEUP_TRACK_DTYPE = np.float64

class PileupTrack(ArrayTrack):
    
    def __init__(self, hdf_group):
        super(PileupTrack, self).__init__(hdf_group, PILEUP_TRACK_DTYPE)

    def get_num_reads(self):
        return self.hdf_group._v_attrs.reads
    def get_unmapped_reads(self):
        return self.hdf_group._v_attrs.unmapped_reads
    def get_qcfail_reads(self):
        return self.hdf_group._v_attrs.qcfail_reads
    def get_mapped_reads(self):
        return self.get_num_reads() - self.get_unmapped_reads()
    def get_duplicate_reads(self):
        return self.hdf_group._v_attrs.duplicate_reads
    def get_fragment_length(self):        
        return self.hdf_group._v_attrs.fragment_length if hasattr(self.hdf_group._v_attrs, "fragment_length") else -1
    def get_read_length(self):
        return self.hdf_group._v_attrs.read_length
    def get_read_length_distribution(self):
        return self.hdf_group._v_attrs.read_length_distribution
    def get_mismatches_per_cycle(self):
        return self.hdf_group._v_attrs.cycle_mismatches
    def get_mismatches_per_read(self):
        return self.hdf_group._v_attrs.read_mismatches

    def coverage(self, interval, norm=True):
        ref, start, end = parse_interval(interval)
        assert ref in self.hdf_group
        norm_length = max(self.get_fragment_length(), self.get_read_length())        
        ca = self.hdf_group._f_getChild(ref)
        return ca[start:end] * (1.0e6 / (self.get_mapped_reads() * norm_length))

    def count(self, interval):
        ref, start, end = parse_interval(interval)        
        assert ref in self.hdf_group
        group_attrs = self.hdf_group._v_attrs        
        ca = self.hdf_group._f_getChild(ref)
        return np.sum(ca[start:end]) / float(group_attrs.read_length)
        
    def rpkm(self, interval):
        ref, start, end = parse_interval(interval)        
        assert ref in self.hdf_group
        group_attrs = self.hdf_group._v_attrs        
        ca = self.hdf_group._f_getChild(ref)
        return np.sum(ca[start:end]) * (1.0e9 / (group_attrs.reads * group_attrs.read_length * (end - start)))

    def to_bedgraph(self, fileh, span=1, norm=False, mirror=False):
        span = max(1, span)
        factor = -1.0 if mirror else 1.0
        norm_length = max(self.get_fragment_length(), self.get_read_length())
        if norm:
            factor = factor * (1.0e9  / (self.get_mapped_reads() * norm_length))
        for ref in self._get_references():
            ca = self._get_array(ref)
            array_to_bedgraph(ref, ca, fileh, span=span, factor=factor)
    
    def insert_sam(self, filename, genome=False, fragment_length=-1, rmdup=False,
                   multimap=False):
        bamfh = pysam.Samfile(filename, "rb")
        # count reads and get other info from BAM file
        logging.debug("%s: Computing track statistics" % (filename))
        hits = 0
        reads = 0.0
        duplicate_reads = 0.0
        unmapped_reads = 0.0
        mate_unmapped_reads = 0.0
        qcfail_reads = 0.0
        read_lengths = collections.defaultdict(lambda: 0.0)
        multihit_counts = collections.defaultdict(lambda: 0.0)
        cycle_mismatch_counts = collections.defaultdict(lambda: 0.0)
        read_mismatch_counts = collections.defaultdict(lambda: 0.0)
        mismatch_re = re.compile(r"(\d+)([agtcAGTC]*)")

        for read in bamfh:
            hits += 1
            if read.mate_is_unmapped:
                mate_unmapped_reads += 1
            if read.is_unmapped or read.is_qcfail:
                unmapped_reads += 1
                num_read_hits = 0
                if read.is_qcfail:
                    qcfail_reads += 1
            else:            
                if multimap:
                    num_read_hits = read.opt('NH')
                else:
                    num_read_hits = 1
                weighted_cov = 1.0 / num_read_hits
                # count reads
                reads += weighted_cov
                if read.is_duplicate:
                    duplicate_reads += weighted_cov
                read_lengths[read.rlen] += weighted_cov
                # mismatch counts along read
                mm_iter = mismatch_re.finditer(read.opt('MD'))
                offset = 0
                for mismatch in mm_iter:
                    skip, mmbase = mismatch.groups()
                    offset += int(skip)
                    if mmbase:
                        pos = (read.rlen - offset - 1) if read.is_reverse else 0
                        cycle_mismatch_counts[pos] += weighted_cov
                # total mismatches per read
                read_mismatch_counts[read.opt('NM')] += weighted_cov
            multihit_counts[num_read_hits] += 1

        # number of reads
        parent_group = self.hdf_group
        parent_group._v_attrs.hits = hits
        logging.info("%s: Number of alignment hits=%d" % (filename,hits))
        parent_group._v_attrs.reads = reads
        logging.info("%s: Number of reads=%d" % (filename,int(reads)))
        parent_group._v_attrs.duplicate_reads = duplicate_reads
        logging.info("%s: Number of qc fail reads=%d" % (filename,int(qcfail_reads)))
        parent_group._v_attrs.qcfail_reads = qcfail_reads           
        logging.info("%s: Number of read flagged as duplicates=%d" % (filename,int(duplicate_reads)))
        parent_group._v_attrs.unmapped_reads = unmapped_reads
        logging.info("%s: Number of mapped reads=%d" % (filename,int(reads-unmapped_reads)))
        logging.info("%s: Number of unmapped reads=%d" % (filename,unmapped_reads))
        parent_group._v_attrs.mate_unmapped_reads = mate_unmapped_reads
        logging.info("%s: Number of reads with unmapped mate=%d" % (filename,mate_unmapped_reads))
        parent_group._v_attrs.read_length_variable = True if len(read_lengths) > 1 else False
        # set to the most abundant read length in the case of variable read length runs
        parent_group._v_attrs.read_length = int(sorted(read_lengths.items(), key=operator.itemgetter(1), reverse=True)[0][0])
        logging.info("%s: Most common read length=%d" % (filename, parent_group._v_attrs.read_length))
        parent_group._v_attrs.read_length_distribution = dict(read_lengths)
        logging.info("%s: Number of different read lengths=%d" % (filename, len(read_lengths)))
        parent_group._v_attrs.fragment_length = fragment_length
        logging.info("%s: Fragment length (specified as input)=%d" % (filename, parent_group._v_attrs.fragment_length))

        max_multihits = max(multihit_counts.keys())
        arr = np.zeros(max_multihits, dtype=np.int)
        for i in xrange(max_multihits):
            if i in multihit_counts:
                arr[i] = int(multihit_counts[i])
        parent_group._v_attrs.multihits = list(arr)
        logging.info("%s: Multimapping read counts=%s" % (filename, parent_group._v_attrs.multihits))

        max_rlen = max(read_lengths.keys())
        arr = np.zeros(max_rlen, dtype=np.int)
        for i in xrange(max_rlen):
            if i in cycle_mismatch_counts:
                arr[i] = int(cycle_mismatch_counts[i])
        parent_group._v_attrs.cycle_mismatches = list(arr)
        logging.info("%s: Mismatches per cycle=%s" % (filename, parent_group._v_attrs.cycle_mismatches))

        max_mismatches = max(read_mismatch_counts.keys())
        arr = np.zeros(max_mismatches, dtype=np.int)
        for i in xrange(max_mismatches):
            if i in read_mismatch_counts:
                arr[i] = int(read_mismatch_counts[i])
        parent_group._v_attrs.read_mismatches = list(arr)
        logging.info("%s: Mismatches per read=%s" % (filename, parent_group._v_attrs.read_mismatches))

        # insert data
        for ref,length in zip(self._get_references(), self._get_lengths()):
            if ref not in bamfh.references:
                logging.debug("Skipping reference %s in BAM file that is not part of the track" % (ref))
                continue
            logging.debug("Adding file %s ref %s to database" % (filename, ref))
            ca = self.hdf_group._f_getChild(ref)
            # insert the reads from the BAM file into the array
            if genome:
                bam_to_array_chunked(bamfh.fetch(ref), ca, fragment_length=fragment_length, rmdup=rmdup)
            elif multimap:
                bam_pileup_to_array_multihit(bamfh.pileup(ref), ca, dtype=ca.atom.dtype)
            else:
                bam_pileup_to_array(bamfh.pileup(ref), ca, dtype=ca.atom.dtype)
        bamfh.close()
