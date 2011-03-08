'''
Created on Mar 6, 2011

@author: mkiyer
'''
import array
import collections

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

class BamCoverageStatistics:
    def __init__(self):
        self.num_reads = 0
        self.total_cov = 0.0
        self.read_lengths = collections.defaultdict(lambda: 0.0)
        
    
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
            self.stats.read_lengths[read.rlen] += 1
            # find genomic intervals of read alignment
            for start, end, seq in get_genomic_intervals(read):
                #print 'START', start, 'END', end, 'SEQ', seq        
                if self.norm_rlen:
                    cov /= (end - start)
                #print start, end, read.is_reverse, cov, seq
                self.stats.total_cov += (end - start) * cov                
                self.intervals.append((start, end, read.is_reverse, cov, seq))
