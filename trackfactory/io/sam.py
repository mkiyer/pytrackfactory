'''
Created on Mar 6, 2011

@author: mkiyer
'''
import array
import collections
import operator

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

class BamCoverageStatistics(object):
    NUM_READS_ATTR = "num_reads"
    STRAND_READS_ATTR = "strand_reads"
    TOTAL_COV_ATTR = "total_cov"
    STRAND_COV_ATTR = "strand_cov"
    RLEN_ATTR = "rlen"
    RLEN_DIST_ATTR = "rlen_dist"
    
    def __init__(self):
        self.num_reads = 0
        self.strand_reads = [0, 0]
        self.total_cov = 0.0
        self.strand_cov = [0.0, 0.0]
        self.rlen_dict = collections.defaultdict(lambda: 0.0)
    
    def add(self, read, cov):
        strand = int(read.is_reverse)
        self.num_reads += 1
        self.strand_reads[strand] += 1
        self.rlen_dict[read.rlen] += 1
        self.total_cov += cov
        self.strand_cov[strand] += cov

    def update(self, other):
        self.num_reads += other.num_reads
        self.total_cov += other.total_cov
        for rlen in other.rlen_dict:
            if rlen not in self.rlen_dict:
                self.rlen_dict[rlen] = 0
            self.rlen_dict[rlen] += other.rlen_dict[rlen]

    def fromhdf(self, hdf_group):
        if self.NUM_READS_ATTR not in hdf_group._v_attrs:
            return
        self.num_reads = hdf_group._v_attrs[self.NUM_READS_ATTR]
        self.total_cov = hdf_group._v_attrs[self.TOTAL_COV_ATTR]
        self.strand_reads = hdf_group._v_attrs[self.STRAND_READS_ATTR]
        self.strand_cov = hdf_group._v_attrs[self.STRAND_COV_ATTR]
        rlen_dict = hdf_group._v_attrs[self.RLEN_DIST_ATTR]
        for rlen in rlen_dict:
            self.rlen_dict[rlen] += rlen_dict[rlen]    
        # most common read length
        rlen = int(sorted(self.rlen_dict.items(), 
                          key=operator.itemgetter(1), 
                          reverse=True)[0][0])
        self.read_length = rlen

    def tohdf(self, hdf_group):
        hdf_group._v_attrs[self.NUM_READS_ATTR] = self.num_reads
        hdf_group._v_attrs[self.TOTAL_COV_ATTR] = self.total_cov
        hdf_group._v_attrs[self.RLEN_DIST_ATTR] = dict(self.rlen_dict) 
        hdf_group._v_attrs[self.STRAND_READS_ATTR] = self.strand_reads
        hdf_group._v_attrs[self.STRAND_COV_ATTR] = self.strand_cov
        
    
    
class BamCoverageIterator:
    """
    :param norm_rlen: whether to normalize pileup coverage by read length
    :param num_hits_tag: samtools tag that contains number of hits per read
    :param hit_prob_tag: (optional)
    :param keep_dup: whether to omit duplicate reads
    :param ignore qcfail: whether to consider marked as QC fail
    """
    def __init__(self, bamfh, norm_rlen=False, 
                 num_hits_tag=None, hit_prob_tag=None,
                 max_multimaps=None, keep_dup=True,
                 keep_qcfail=True):
        self.bamfh = bamfh
        self.read_iterator = bamfh.fetch()
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
            # get reference name
            rname = self.bamfh.getrname(read.tid)
            # find genomic intervals of read alignment
            total_cov = 0.0
            for start, end, seq in get_genomic_intervals(read):
                #print 'START', start, 'END', end, 'SEQ', seq
                interval_cov = cov
                if self.norm_rlen:
                    interval_cov /= (end - start)
                total_cov += (end - start) * interval_cov
                #print start, end, read.is_reverse, cov, seq
                self.intervals.append((rname, start, end, read.is_reverse, interval_cov, seq))
            self.stats.add_read(read, total_cov)
