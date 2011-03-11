'''
Created on Mar 9, 2011

@author: mkiyer
'''
import numpy as np
import tables
import logging
import pysam

from track import Track, TrackError, parse_interval
from intervaltrack import IntervalTrack, get_base_dtype_fields, \
    REF_COL_NAME, START_COL_NAME, END_COL_NAME
from vectortrack import VectorTrack
from io.sam import BamCoverageIterator, BamCoverageStatistics

JUNCTION_GROUP = "junctions"
COVERAGE_GROUP = "coverage"

junction_dtype = (get_base_dtype_fields() +
                  [('id', '<u4'),
                   ('seqdata_name', 'a36'),
                   ('strand', 'a4'),
                   ('reads', '<i4'),
                   ('left_coverage', '<u4'),
                   ('right_coverage', '<u4')])
junction_expectedrows = 1e5

def tophat_bed_to_juncs(source, line_iter):
    '''
    chr5    180666535       180668532       JUNC00002848    36      -       180666535       180668532       255,0,0 2       47,41   0,1956
    '''
    juncarr = np.empty(1, dtype=junction_dtype)
    junc = juncarr[0]
    for line in line_iter:
        if not line:
            continue
        if line.startswith("#") or line.startswith("track"):
            continue
        fields = line.strip().split('\t')
        # parse junction
        start = int(fields[1])
        end = int(fields[2])
        left_cov, right_cov = map(int, fields[10].split(','))
        ignore, junc_length = map(int, fields[11].split(','))
        assert ignore == 0
        assert start + junc_length + right_cov == end
        junc_start = start + left_cov
        junc_end = start + junc_length
        #print '\t'.join(map(str, [fields[0], junc_start, junc_end, fields[3], fields[4], fields[5]]))
        # store junction in track
        junc['id'] = 0
        junc['seqdata_name'] = source
        junc[REF_COL_NAME] = fields[0]
        junc[START_COL_NAME] = junc_start
        junc[END_COL_NAME] = junc_end
        junc['reads'] = int(fields[4])
        junc['strand'] = fields[5]        
        junc['left_coverage'] = left_cov
        junc['right_coverage'] = right_cov
        yield junc

class RnaseqTrack(Track):
    '''contains both genomic coverage data (array track) and splice
    junction data (interval track)
    '''
    def __init__(self, hdf_group, stranded=False):
        super(RnaseqTrack, self).__init__(hdf_group)        
        h5file = self._get_hdf_file()
        if JUNCTION_GROUP not in self.hdf_group:
            junc_group = h5file.createGroup(self.hdf_group, JUNCTION_GROUP)
        else:
            junc_group = self.hdf_group._f_getChild(JUNCTION_GROUP)            
        self.junc_track = IntervalTrack(junc_group, junction_dtype, 
                                        junction_expectedrows)
        if COVERAGE_GROUP not in self.hdf_group:
            cov_group = h5file.createGroup(self.hdf_group, COVERAGE_GROUP)
        else:
            cov_group = self.hdf_group._f_getChild(COVERAGE_GROUP)
        self.cov_track = VectorTrack(cov_group)

    def get_junction_track(self):
        junc_group = self.hdf_group._f_getChild(JUNCTION_GROUP)
        return IntervalTrack(junc_group)
    def get_coverage_track(self):
        cov_group = self.hdf_group._f_getChild(COVERAGE_GROUP)
        return VectorTrack(cov_group)

    def fromtophat(self, accepted_hits_bam, junctions_bed):
        # insert splice junction track
        track_name = self.hdf_group._v_name
        rnames = set(self.get_rnames())
        logging.info("[RnaseqTrack] adding junctions")        
        junc_iter = tophat_bed_to_juncs(track_name, open(junctions_bed))        
        for junc in junc_iter:
            if junc[REF_COL_NAME] not in rnames:
                logging.debug('Skipping junc %s' % str(junc))
            else:
                junc['id'] = self.junc_track.num_intervals
                self.junc_track.add(junc)
        self.junc_track.index(persist=True)
        # insert coverage track
        logging.info("creating coverage track")        
        bamfh = pysam.Samfile(accepted_hits_bam, "rb")        
        #cmdline = bamfh.header["PG"][0]["CL"]        
        #re.search(r'--max-multihits(?:\s+|=)(\d+)', cmdline)
        intervalcoviter = BamCoverageIterator(bamfh,
                                              norm_rlen=True,
                                              num_hits_tag="NH",
                                              hit_prob_tag=None,
                                              max_multimaps=None,
                                              keep_dup=True,
                                              keep_qcfail=False)
        self.cov_track.fromintervals(intervalcoviter)
        # store coverage statistics to allow normalization calculations
        stats = intervalcoviter.stats
        logging.debug("\tProcessed '%d' valid reads" % (stats.num_reads))
        logging.debug("\tTotal coverage '%f'" % (stats.total_cov))
        bamfh.close()
