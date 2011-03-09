'''
.. module:: arraytrack
   :synopsis: track for storing 1-D arrays
.. moduleauthor:: Matthew Iyer <mkiyer@umich.edu>
.. versionadded:: 0.0.1

Created on Sep 14, 2010

@author: mkiyer
'''
import logging
import tables
import numpy as np
from track import Track, TrackError, parse_interval
from io.interval import write_interval_data_to_array
from io.cbedgraph import array_to_bedgraph

DTYPE_ATTR = 'dtype'

class ArrayTrack(Track):
    h5_complevel = 1
    h5_complib = 'lzo'
    h5_filters = tables.Filters(complevel=h5_complevel, 
                                complib=h5_complib,
                                shuffle=True,
                                fletcher32=False)    
    h5_chunksize = 1<<18 # 256kb chunks
    
    def __init__(self, hdf_group, dtype=None):
        super(ArrayTrack, self).__init__(hdf_group)
        if not self._check_arrays():
            if (dtype is None):
                raise TrackError("Must pass valid numpy dtype to ArrayTrack (got '%s')" % (str(dtype)))
            self._create_arrays(np.dtype(dtype))
            hdf_group._v_attrs[DTYPE_ATTR] = np.dtype(dtype)

    def _check_arrays(self):
        # check if arrays exist
        hasarrays = True
        for rname in self.get_rnames():
            if rname in self.hdf_group:
                continue
            hasarrays = False
            break
        return hasarrays

    def _create_arrays(self, dtype):        
        # compute chunks sizes tuple (find nearest power of 2)        
        chunksize = 2**int(np.log2(self.h5_chunksize/float(dtype.itemsize)))
        atom = tables.Atom.from_dtype(np.dtype(dtype))
        h5file = self._get_hdf_file()
        # create a separate array for each reference
        for ref, length in self.get_refs():
            if ref in self.hdf_group:
                continue
            # adjust chunksize relative to length
            if chunksize > length:
                chunksize = length
            chunkshape = (chunksize,)
            h5file.createCArray(self.hdf_group, ref, 
                                atom=atom,
                                shape=(length,),
                                filters=self.h5_filters,
                                chunkshape=chunkshape)

    def _get_dtype(self):
        return self.hdf_group._v_attrs[DTYPE_ATTR]
    
    def _get_array(self, ref):
        return self.hdf_group._f_getChild(ref)
    
    def _get_arrays(self):
        myarrays = {}
        for rname in self.get_rnames():
            myarrays[rname] = self._get_array(rname)            
        return myarrays

    def _check_bounds(self, arr, start, end):
        if ((start < 0) or (end > arr.shape[0])):
            raise TrackError("index out of bounds '(%d,%d)'" % (start,end))

    def __getitem__(self, key):
        arr, start, end = self._parse_interval(key)
        return arr[start:end]
    
    def __setitem__(self, key, value):
        arr, start, end = self._parse_interval(key)
        arr[start:end] = value

    def tobedgraph(self, interval, fileh, span=1, norm=False, mirror=False):
        span = max(1, span)
        factor = -1.0 if mirror else 1.0        
        if norm:
            factor = factor * (1.0e6  / (self.get_total_coverage()))
        ref, start, end = parse_interval(interval)
        if ref is None: 
            rnames = self.get_rnames()
        else:
            rnames = [ref]
        if start is None: start = 0
        if end is None: end = -1
        for rname in rnames:
            array_to_bedgraph(rname, self._get_array(rname), fileh, 
                              start=start, end=end, factor=factor, span=span,
                              chunksize=self.h5_chunksize)

    def fromintervals(self, interval_iter):
        dtype = np.dtype(self._get_dtype())
        rname_array_dict = self._get_arrays()
        write_interval_data_to_array(interval_iter, 
                                     rname_array_dict, 
                                     dtype=self._get_dtype(), 
                                     chunksize=(self.h5_chunksize << 4))

