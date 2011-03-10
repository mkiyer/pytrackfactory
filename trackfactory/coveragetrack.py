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
from track import TrackError, parse_interval, NO_STRAND, POS_STRAND, NEG_STRAND
from arraytrack import ArrayTrack
from io.interval import write_interval_data_to_array, \
    write_interval_data_to_stranded_array
from io.bedgraph import array_to_bedgraph

COVERAGE_DTYPE = np.float32
NORM_RLEN_ATTR = "norm_rlen"
TOTAL_COV_ATTR = "total_cov"
NUM_FEATURES_ATTR = "num_features"

class CoverageTrack(ArrayTrack):
    
    def __init__(self, hdf_group):
        super(CoverageTrack, self).__init__(hdf_group, COVERAGE_DTYPE)
        self._init_attrs()

    def _init_attrs(self):
        # store parameters used tabulate coverage
        if TOTAL_COV_ATTR not in self.hdf_group._v_attrs:
            self.hdf_group._v_attrs[TOTAL_COV_ATTR] = 0
            self.hdf_group._v_attrs[NUM_FEATURES_ATTR] = 0
        self.total_cov = self._get_total_cov()
        self.num_features = self._get_num_features()

    def _get_total_cov(self):
        return self.hdf_group._v_attrs[TOTAL_COV_ATTR]
    def _get_num_features(self):
        return self.hdf_group._v_attrs[NUM_FEATURES_ATTR]

    def count(self, interval, strand=None):
        arr, start, end, strand = self._parse_interval(interval)
        self._check_bounds(arr, start, end)
        return np.sum(arr[start:end])

    def coverage(self, interval, multiplier=1.0e6):
        arr, start, end, strand = self._parse_interval(interval)
        self._check_bounds(arr, start, end)
        return arr[start:end] * (multiplier / self.total_cov)

    def density(self, interval, multiplier=1.0e9):
        arr, start, end, strand = self._parse_interval(interval)
        self._check_bounds(arr, start, end)
        return np.sum(arr[start:end]) * (multiplier / (self.total_cov * (end - start)))

    def tobedgraph(self, interval, fileh, span=1, factor=1.0, 
                   norm=False, multiplier=1.0e6):
        if norm:
            factor *= (multiplier  / (self.total_cov))
        ArrayTrack.tobedgraph(self, interval, fileh, span, factor)

    def fromintervals(self, interval_iter, channel=0):
        if channel is None: channel = 0
        rname_array_dict = self._get_arrays()
        intervals, total_cov = \
            write_interval_data_to_array(interval_iter, 
                                         rname_array_dict, 
                                         dtype=self._get_dtype(),
                                         channel=channel,                                     
                                         chunksize=(self.h5_chunksize << 4))
        self.hdf_group._v_attrs[TOTAL_COV_ATTR] = total_cov
        self.total_cov = total_cov
        self.hdf_group._v_attrs[NUM_FEATURES_ATTR] = intervals
        self.num_features = intervals        

class StrandedCoverageTrack(CoverageTrack):
    
    def __init__(self, hdf_group):
        ArrayTrack.__init__(self, hdf_group, COVERAGE_DTYPE, channels=2)
        self._init_attrs()

    def _count(self, arr, start, end, strand):
        if strand == NO_STRAND:
            return arr[start:end].sum()
        else:
            return arr[start:end,strand].sum()
        
    def count(self, interval):
        arr, start, end, strand = self._parse_interval(interval)
        self._check_bounds(arr, start, end)
        return self._count(arr, start, end, strand)

    def coverage(self, interval, multiplier=1.0e6):
        arr, start, end, strand = self._parse_interval(interval)
        self._check_bounds(arr, start, end)
        if strand == NO_STRAND:
            data = arr[start:end].sum(axis=1)
        else:
            data = arr[start:end,strand]
        return data * (multiplier / self.total_cov)

    def density(self, interval, multiplier=1.0e9):
        arr, start, end, strand = self._parse_interval(interval)
        self._check_bounds(arr, start, end)
        count = self._count(arr, start, end, strand)
        return count * (multiplier / (self.total_cov * (end - start)))

    def tobedgraph(self, interval, fileh, span=1, factor=1.0, 
                   norm=False, multiplier=1.0e6):
        if span < 1: span = 1
        ref, start, end, strand = parse_interval(interval)
        if ref is None: 
            rnames = self.get_rnames()
        else:
            rnames = [ref]
        if start is None: start = 0
        if end is None: end = -1
        if strand == NO_STRAND:
            channels = None
        else:
            channels = strand            
        if norm:
            factor *= (multiplier / (self.total_cov))
        for rname in rnames:
            array_to_bedgraph(rname, self._get_array(rname), fileh, 
                              start=start, end=end, factor=factor, span=span,
                              chunksize=self.h5_chunksize, channels=channels)

    def fromintervals(self, interval_iter):
        rname_array_dict = self._get_arrays()
        intervals, total_cov = \
            write_interval_data_to_stranded_array(interval_iter, 
                                                  rname_array_dict, 
                                                  dtype=self._get_dtype(),
                                                  chunksize=(self.h5_chunksize << 4))
        self.hdf_group._v_attrs[TOTAL_COV_ATTR] = total_cov
        self.total_cov = total_cov
        self.hdf_group._v_attrs[NUM_FEATURES_ATTR] = intervals
        self.num_features = intervals