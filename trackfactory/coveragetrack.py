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
import operator
import numpy as np
from track import TrackError, parse_interval
from arraytrack import ArrayTrack
from io.sam import BamCoverageIterator, BamCoverageStatistics
from io.interval import write_interval_data_to_array
from io.cbedgraph import array_to_bedgraph
#from io.bedgraph import array_to_bedgraph

COVERAGE_TRACK_DTYPE = np.float32
NORM_RLEN_ATTR = "norm_rlen"

class CoverageTrack(ArrayTrack):
    
    def __init__(self, hdf_group):
        super(CoverageTrack, self).__init__(hdf_group, COVERAGE_TRACK_DTYPE)
        self.stats = BamCoverageStatistics()
        self.stats.fromhdf(self.hdf_group)

    def get_total_coverage(self):        
        return self.hdf_group._v_attrs[self.stats.total_cov]

    def count(self, interval):
        ref, start, end = parse_interval(interval)        
        assert ref in self.hdf_group
        ca = self.hdf_group._f_getChild(ref)
        return np.sum(ca[start:end])
        
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
        if start is None: start = 0
        if end is None: end = -1
        arr = self._get_array(ref)
        array_to_bedgraph(ref, arr, fileh, 
                          start=start, end=end, factor=factor, span=span,
                          chunksize=self.h5_chunksize)

    def frombam(self, bamfh, 
                norm_rlen=False,
                num_hits_tag=None,
                hit_prob_tag=None,
                max_multimaps=None,
                keep_dup=True,
                keep_qcfail=True):
        rname_array_dict = self._get_arrays()
        # keep aggregated statistics
        intervalcoviter = BamCoverageIterator(bamfh,
                                              norm_rlen, 
                                              num_hits_tag,
                                              hit_prob_tag,
                                              max_multimaps=max_multimaps, 
                                              keep_dup=keep_dup,
                                              keep_qcfail=keep_qcfail)
        write_interval_data_to_array(intervalcoviter, 
                                     rname_array_dict, 
                                     dtype=self._get_dtype(), 
                                     chunksize=(self.h5_chunksize << 4))
        # store coverage statistics to allow calculations
        self.stats = intervalcoviter.stats
        logging.debug("\tProcessed '%d' valid reads" % (self.stats.num_reads))
        logging.debug("\tTotal coverage '%f'" % (self.stats.total_cov))
        # store parameters used tabulate coverage
        self.hdf_group._v_attrs[NORM_RLEN_ATTR] = norm_rlen
        # store statistics
        self.stats.tohdf(self.hdf_group)
#
#        self.stats = bamstats
#        
#        for rname in self.get_rnames():
#            if rname not in bamfh.references:
#                logging.debug("Reference %s not found in BAM header" % (rname))
#                continue
#            logging.debug("Adding ref '%s'" % (rname))            
#            arr = self._get_array(rname)
#            write_interval_data_to_array(intervalcoviter, 
#                                         rname_array_dict, 
#                                         dtype=self._get_dtype(), 
#                                         chunksize=(self.h5_chunksize << 4))
#            # store coverage statistics to allow calculations
#            refstats = intervalcoviter.stats
#            refstats.tohdf(arr)
#            bamstats.update(refstats)
#            logging.debug("\tProcessed '%d' valid reads" % (refstats.num_reads))
#            logging.debug("\tTotal coverage '%f'" % (refstats.total_cov))
#        # store parameters used tabulate coverage
#        self.hdf_group._v_attrs[NORM_RLEN_ATTR] = norm_rlen
#        # store statistics
#        bamstats.tohdf(self.hdf_group)
#        self.stats = bamstats
