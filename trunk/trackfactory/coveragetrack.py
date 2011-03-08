'''
.. module:: arraytrack
   :synopsis: track for storing 1-D arrays
.. moduleauthor:: Matthew Iyer <mkiyer@umich.edu>
.. versionadded:: 0.0.1

Created on Mar 7, 2011

@author: mkiyer
'''
import tables
import logging
import numpy as np
from track import TrackError, parse_interval
from arraytrack import ArrayTrack
from io.sam import BamCoverageIterator
from io.interval import write_interval_data_to_array
from io.bedgraph import array_to_bedgraph

COVERAGE_TRACK_DTYPE = np.float32

NUM_READS_ATTR = "num_reads"
TOTAL_COV_ATTR = "total_cov"
NORM_RLEN_ATTR = "norm_rlen"

class CoverageTrack(ArrayTrack):
    
    def __init__(self, hdf_group):
        super(CoverageTrack, self).__init__(hdf_group, COVERAGE_TRACK_DTYPE)

    def get_total_coverage(self):
        return self.hdf_group._v_attrs[TOTAL_COV_ATTR] 

    def count(self, interval):
        ref, start, end = parse_interval(interval)        
        assert ref in self.hdf_group
        group_attrs = self.hdf_group._v_attrs        
        ca = self.hdf_group._f_getChild(ref)
        return np.sum(ca[start:end]) / float(group_attrs.read_length)
        
    def coverage(self, interval, norm=True):
        ref, start, end = parse_interval(interval)
        assert ref in self.hdf_group
        ca = self.hdf_group._f_getChild(ref)
        return ca[start:end] * (1.0e6 / self.get_total_coverage())

    def rpkm(self, interval):
        ref, start, end = parse_interval(interval)        
        assert ref in self.hdf_group
        ca = self.hdf_group._f_getChild(ref)
        return np.sum(ca[start:end]) * (1.0e9 / (self.get_total_coverage() * (end - start)))

    def tobedgraph(self, interval, fileh, span=1, norm=False, mirror=False):
        span = max(1, span)
        factor = -1.0 if mirror else 1.0        
        if norm:
            factor = factor * (1.0e6  / (self.get_total_coverage()))
        ref, start, end = parse_interval(interval)
        arr = self._get_array(ref)
        # use the iterrows() method to iterate over the array
        array_iter = arr.iterrows(start=start, stop=end)
        array_to_bedgraph(ref, array_iter, fileh, factor=factor)
#            array_to_bedgraph(ref, ca, fileh, span=span, factor=factor)

    def frombam(self, bamfh, 
                norm_rlen=False,
                num_hits_tag=None,
                hit_prob_tag=None,
                max_multimaps=None,
                keep_dup=True,
                keep_qcfail=True):
        # keep aggregated statistics
        num_reads = 0
        total_cov = 0.0
        for rname in self.get_rnames():
            if rname not in bamfh.references:
                logging.debug("Reference %s not found in BAM header" % (rname))
                continue
            logging.debug("Adding ref '%s'" % (rname))            
            intervalcoviter = BamCoverageIterator(bamfh.fetch(rname),
                                                  norm_rlen, 
                                                  num_hits_tag,
                                                  hit_prob_tag,
                                                  max_multimaps=max_multimaps, 
                                                  keep_dup=keep_dup,
                                                  keep_qcfail=keep_qcfail)
            arr = self._get_array(rname)
            write_interval_data_to_array(intervalcoviter, arr, 
                                         chunksize=(self.h5_chunksize << 4))
            # store coverage statistics to allow calculations
            bamstats = intervalcoviter.stats
            arr._v_attrs[NUM_READS_ATTR] = bamstats.num_reads
            arr._v_attrs[TOTAL_COV_ATTR] = bamstats.total_cov
            num_reads += bamstats.num_reads
            total_cov += bamstats.total_cov
            logging.debug("\tProcessed '%d' valid reads" % (bamstats.num_reads))
            logging.debug("\tTotal coverage '%f'" % (bamstats.total_cov))
        # store statistics
        self.hdf_group._v_attrs[NUM_READS_ATTR] = num_reads
        self.hdf_group._v_attrs[TOTAL_COV_ATTR] = total_cov
        self.hdf_group._v_attrs[NORM_RLEN_ATTR] = norm_rlen
        #parent_group._v_attrs.read_length = int(sorted(read_lengths.items(), key=operator.itemgetter(1), reverse=True)[0][0])

