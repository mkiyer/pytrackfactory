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
from track import TrackError, NO_STRAND, POS_STRAND, NEG_STRAND
from arraytrack import ArrayTrack
from io.interval import write_interval_data_to_array
#from io.interval import write_interval_data_to_array, \
#    write_interval_data_to_stranded_array, \
#    write_interval_data_to_stranded_allele_array
from io.bedgraph import array_to_bedgraph

_vector_dtypes = {"i": np.int32,
                  "f": np.float32}
def check_vector_dtype(dtype):
    if dtype not in _vector_dtypes:
        raise TrackError("'dtype' must be one of '%s'" % 
                         (str(_vector_dtypes.keys())))
    dtype = np.dtype(_vector_dtypes[dtype])
    return dtype

NORM_RLEN_ATTR = "norm_rlen"
TOTAL_COUNT_BY_REF = "count_by_ref"
TOTAL_COUNT = "total_count"
NUM_FEATURES_BY_REF = "num_features_by_ref"
NUM_FEATURES = "num_features"

class VectorTrack(ArrayTrack):
    
    def __init__(self, hdf_group, dtype="f", channels=1):
        dtype = check_vector_dtype(dtype)
        super(VectorTrack, self).__init__(hdf_group, dtype, channels)
        self._init_attrs()

    def _init_attrs(self):
        # store parameters used tabulate coverage
        if TOTAL_COUNT_BY_REF not in self.hdf_group._v_attrs:
            self.hdf_group._v_attrs[TOTAL_COUNT_BY_REF] = {}
            self.hdf_group._v_attrs[TOTAL_COUNT] = 0
            self.hdf_group._v_attrs[NUM_FEATURES_BY_REF] = {}
            self.hdf_group._v_attrs[NUM_FEATURES] = 0
        self.total = self._get_total_count()
        self.num_features = self._get_num_features()

    def _get_total_count(self):
        return self.hdf_group._v_attrs[TOTAL_COUNT]
    def _get_num_features(self):
        return self.hdf_group._v_attrs[NUM_FEATURES]
    
    def _set_count_attrs(self, total_dict, num_features_dict):
        self.total = sum(total_dict.values())
        self.hdf_group._v_attrs[TOTAL_COUNT_BY_REF] = dict(total_dict)
        self.hdf_group._v_attrs[TOTAL_COUNT] = self.total
        self.num_features = sum(num_features_dict.values())        
        self.hdf_group._v_attrs[NUM_FEATURES_BY_REF] = dict(num_features_dict)
        self.hdf_group._v_attrs[NUM_FEATURES] = self.num_features
    
    def _select_channels(self, strand=None, alleles=None, channel=0):
        return (channel,)

    def count(self, interval, channel=0):
        ref, start, end, strand = self._parse_interval(interval)
        arr = self._get_array(ref)
        self._check_bounds(arr, start, end)
        channels = self._select_channels(strand, channel=channel)
        return arr[start:end,channels].sum()

    def coverage(self, interval, multiplier=1.0e6, channel=0):
        ref, start, end, strand = self._parse_interval(interval)
        arr = self._get_array(ref)
        self._check_bounds(arr, start, end)
        channels = self._select_channels(strand, channel=channel)
        data = arr[start:end,channels].sum(axis=1)
        return data * (multiplier / self.total)

    def density(self, interval, multiplier=1.0e9, channel=0):
        ref, start, end, strand = self._parse_interval(interval)
        arr = self._get_array(ref)
        self._check_bounds(arr, start, end)
        channels = self._select_channels(strand, channel=channel)
        count = arr[start:end,channels].sum()
        return count * (multiplier / (self.total * (end - start)))

    def tobedgraph(self, interval, fileh, span=1, factor=1.0,
                   norm=False, multiplier=1.0e6, channel=0):
        ref, start, end, strand = self._parse_interval(interval)
        if ref is None: 
            rnames = self.get_rnames()
        else:
            rnames = [ref]
        if start is None: start = 0
        if end is None: end = -1
        if span < 1: span = 1
        if norm: factor *= (multiplier  / (self.total))
        channels = self._select_channels(strand, channel)
        for rname in rnames:
            array_to_bedgraph(rname, self._get_array(rname), fileh, 
                              start=start, end=end, factor=factor, span=span,
                              chunksize=self.h5_chunksize,
                              channels=channels)

    def fromintervals(self, interval_iter, channel=0):
        rname_array_dict = self._get_arrays()
        num_features_dict, total_dict = \
            write_interval_data_to_array(interval_iter, 
                                         rname_array_dict, 
                                         dtype=self._get_dtype(),
                                         chunksize=(self.h5_chunksize << 4),
                                         mode="channel",
                                         channel=channel)
        self._set_count_attrs(total_dict, num_features_dict)                                         

class StrandedVectorTrack(VectorTrack):
    
    def __init__(self, hdf_group, dtype="f"):
        dtype = check_vector_dtype(dtype)
        ArrayTrack.__init__(self, hdf_group, dtype, channels=2)
        self._init_attrs()

    def _select_channels(self, strand, alleles=None, channel=0):
        if strand == NO_STRAND: return (0,1)
        return (strand,)

    def fromintervals(self, interval_iter):
        rname_array_dict = self._get_arrays()
        num_features_dict, total_dict = \
            write_interval_data_to_array(interval_iter, 
                                         rname_array_dict, 
                                         dtype=self._get_dtype(),
                                         chunksize=(self.h5_chunksize << 4),
                                         mode="strand")
        self._set_count_attrs(total_dict, num_features_dict)                                         

class StrandedAlleleVectorTrack(StrandedVectorTrack):

    def __init__(self, hdf_group, dtype="f"):
        dtype = check_vector_dtype(dtype)
        ArrayTrack.__init__(self, hdf_group, dtype, channels=8)
        self._init_attrs()

    def _select_channels(self, strand, alleles=None, channel=0):
        if strand == NO_STRAND: return (0,1,2,3,4,5,6,7)
        elif strand == POS_STRAND: return (0,1,2,3)
        else: return (4,5,6,7)

    def fromintervals(self, interval_iter):
        rname_array_dict = self._get_arrays()
        num_features_dict, total_dict = \
            write_interval_data_to_array(interval_iter, 
                                         rname_array_dict, 
                                         dtype=self._get_dtype(), 
                                         chunksize=(self.h5_chunksize << 2),
                                         mode="allele")
        self._set_count_attrs(total_dict, num_features_dict)
