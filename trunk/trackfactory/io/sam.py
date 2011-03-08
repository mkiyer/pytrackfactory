'''
Created on Mar 6, 2011

@author: mkiyer
'''
import array
import numpy as np

CIGAR_M = 0 #match  Alignment match (can be a sequence match or mismatch)
CIGAR_I = 1 #insertion  Insertion to the reference
CIGAR_D = 2 #deletion  Deletion from the reference
CIGAR_N = 3 #skip  Skipped region from the reference
CIGAR_S = 4 #softclip  Soft clip on the read (clipped sequence present in <seq>)
CIGAR_H = 5 #hardclip  Hard clip on the read (clipped sequence NOT present in <seq>)
CIGAR_P = 6 #padding  Padding (silent deletion from the padded reference sequence)

def get_genomic_intervals(read):
    intervals = []
    rseq = read.seq
    qseq = array.array('c')
    qstart = 0
    astart = read.pos
    aend = astart
    for op,length in read.cigar:
        if (op == CIGAR_D):
            aend += length
        elif (op == CIGAR_I) or (op == CIGAR_S):
            qstart += length
        elif (op == CIGAR_M):            
            qseq.fromstring(rseq[qstart:qstart + length])
            qstart += length
            aend += length
        elif (op == CIGAR_N):
            assert len(qseq) == (aend - astart)
            if aend > astart:
                intervals.append((astart, aend, qseq))
            astart = aend + length
            aend = astart
            qseq = array.array('c')
    if aend > astart:
        assert len(qseq) == (aend - astart)
        intervals.append((astart, aend, qseq))
    assert aend == read.aend
    return intervals


#
#def profile_bam(bamfh):
#    for read in bamfh:
#        hits += 1
#        if read.mate_is_unmapped:
#            mate_unmapped_reads += 1
#        if read.is_unmapped or read.is_qcfail:
#            unmapped_reads += 1
#            num_read_hits = 0
#            if read.is_qcfail:
#                qcfail_reads += 1
#        else:            
#            if multimap:
#                num_read_hits = read.opt('NH')
#            else:
#                num_read_hits = 1
#            weighted_cov = 1.0 / num_read_hits
#            # count reads
#            reads += weighted_cov
#            if read.is_duplicate:
#                duplicate_reads += weighted_cov
#            read_lengths[read.rlen] += weighted_cov
#            # mismatch counts along read
#            mm_iter = mismatch_re.finditer(read.opt('MD'))
#            offset = 0
#            for mismatch in mm_iter:
#                skip, mmbase = mismatch.groups()
#                offset += int(skip)
#                if mmbase:
#                    pos = (read.rlen - offset - 1) if read.is_reverse else 0
#                    cycle_mismatch_counts[pos] += weighted_cov
#            # total mismatches per read
#            read_mismatch_counts[read.opt('NM')] += weighted_cov
#        multihit_counts[num_read_hits] += 1


class BamCoverageStatistics:
    def __init__(self):
        self.num_reads = 0
        self.total_cov = 0.0
    
class BamCoverageIterator:
    """
    :param norm_rlen: whether to normalize pileup coverage by read length
    :param num_hits_tag: samtools tag that contains number of hits per read
    :param hit_prob_tag: (optional)
    :param keep_dup: whether to omit duplicate reads
    :param ignore qcfail: whether to consider marked as QC fail
    """
    def __init__(self, read_iterator, norm_rlen=False, 
                 num_hits_tag=None, hit_prob_tag=None,
                 max_multimaps=None, keep_dup=True,
                 keep_qcfail=True):
        self.read_iterator = read_iterator
        self.norm_rlen = norm_rlen
        self.num_hits_tag = num_hits_tag
        self.hit_prob_tag = hit_prob_tag
        self.max_multimaps = max_multimaps
        self.keep_dup = keep_dup
        self.keep_qcfail = keep_qcfail
        self.stats = BamCoverageStatistics()
        self.intervals = []
        self.done = False

    def __iter__(self):        
        return self
    
    def next(self):
        while True:
            if len(self.intervals) > 0:
                return self.intervals.pop(0)
            if self.done:
                raise StopIteration
            try:
                self.get_next_read()
            except StopIteration:
                self.done = True

    def get_next_read(self):
        while len(self.intervals) == 0:
            read = self.read_iterator.next()
            #print 'READ', read
            #print 'CIGAR', read.cigar
            # check if read is usable
            if (not self.keep_dup) and read.is_duplicate:
                continue
            if (not self.keep_qcfail) and read.is_qcfail:
                continue
            # compute read coverage
            cov = 1.0
            if self.hit_prob_tag is not None:
                cov *= read.opt(self.hit_prob_tag)
            elif self.num_hits_tag is not None:
                nh = read.opt(self.num_hits_tag)
                if (self.max_multimaps is not None) and (nh > self.max_multimaps):
                    continue
                cov /= nh
            self.stats.num_reads += 1
            # find genomic intervals of read alignment
            for start, end, seq in get_genomic_intervals(read):
                #print 'START', start, 'END', end, 'SEQ', seq        
                if self.norm_rlen:
                    cov /= (end - start)
                #print start, end, read.is_reverse, cov, seq
                self.stats.total_cov += (end - start) * cov
                self.intervals.append((start, end, read.is_reverse, cov, seq))

def calc_cov(pileupcolumn, norm_rlen, nh_tag, prob_tag, mmap_max):
    cov = 0
    for pileupread in pileupcolumn.pileups:
        alignment = pileupread.alignment
        readcov = 1.0
        if prob_tag is not None:
            readcov *= alignment.opt(prob_tag)
        if nh_tag is not None:
            nh = alignment.opt(nh_tag)
            if nh > mmap_max:
                continue
            readcov /= nh
        if norm_rlen:
            # TODO: should actually use the genomic footprint of
            # the read rather than the read length, because indels
            # and other events can lead to extra coverage added
            readcov /= alignment.rlen
        cov += readcov
    return cov

def bam_pileup_to_array(pileup_iterator, arr,
                        norm_rlen=False,
                        num_hits_tag=None,
                        hit_prob_tag=None,
                        max_multimaps=None):
    """stores read pileups in a chromosome-sized array
    
    :param norm_rlen: whether to normalize pileup coverage by read length
    :param num_hits_tag: samtools tag that contains number of hits per read
    :param hit_prob_tag: (optional)
    """
    dtype = arr.dtype
    chunk_arr = array.array(dtype.char)
    start = None
    end = None
    for pileupcolumn in pileup_iterator:
        pos = pileupcolumn.pos        
        if start == None:
            start, end = pos, pos
        elif end != pos:
            arr[start:end] = np.array(chunk_arr, dtype=dtype)
            del chunk_arr
            chunk_arr = array.array(dtype.char)
            start, end = pos, pos
        end += 1
        cov = calc_cov(pileupcolumn, norm_rlen, num_hits_tag, hit_prob_tag, 
                       max_multimaps)
        chunk_arr.append(cov)
    # flush remaining array elements
    if len(chunk_arr) > 0:
        arr[start:(start + len(chunk_arr))] = np.array(chunk_arr, dtype=dtype)

def calc_cov_intervals(read_iterator,
                       norm_rlen=False,
                       num_hits_tag=None,
                       hit_prob_tag=None,
                       max_multimaps=None):
    """
    :param norm_rlen: whether to normalize pileup coverage by read length
    :param num_hits_tag: samtools tag that contains number of hits per read
    :param hit_prob_tag: (optional)
    """
    for read in read_iterator:
        strand = "-" if read.is_reverse else "+"
        # compute read coverage
        cov = 1.0
        if hit_prob_tag is not None:
            cov *= read.opt(hit_prob_tag)
        elif num_hits_tag is not None:
            nh = read.opt(num_hits_tag)
            if (max_multimaps is not None) and (nh > max_multimaps):
                continue
            cov /= nh
        # find genomic intervals of read alignment
        for start,end,seq in get_genomic_intervals(read):        
            if norm_rlen:
                cov /= (end - start)
            yield start, end, strand, cov, seq

