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
from track import Track, TrackError
from io.interval import write_interval_data_to_array

DTYPE_ATTR = 'dtype'
CHANNELS_ATTR = 'channels'

class ArrayTrack(Track):
    h5_complevel = 1
    h5_complib = 'lzo'
    h5_filters = tables.Filters(complevel=h5_complevel, 
                                complib=h5_complib,
                                shuffle=True,
                                fletcher32=False)    
    h5_chunksize = 1<<18 # 256kb chunks
    
    def __init__(self, hdf_group, dtype=None, channels=1):
        super(ArrayTrack, self).__init__(hdf_group)
        if not self._check_arrays():
            if (dtype is None):
                raise TrackError("Must pass valid numpy dtype to ArrayTrack (got '%s')" % (str(dtype)))
            self._create_arrays(np.dtype(dtype), channels)

    def _check_arrays(self):
        # check if arrays exist
        hasarrays = True
        for rname in self.get_rnames():
            if rname in self.hdf_group:
                continue
            hasarrays = False
            break
        return hasarrays

    def _create_arrays(self, dtype, ndim):        
        # compute chunks sizes tuple (find nearest power of 2)        
        chunksize = 2**int(np.log2(self.h5_chunksize/(ndim*float(dtype.itemsize))))
        atom = tables.Atom.from_dtype(np.dtype(dtype))
        h5file = self._get_hdf_file()
        # create a separate array for each reference
        for ref, length in self.get_refs():
            if ref in self.hdf_group:
                continue
            # adjust chunksize relative to length
            if chunksize > length:
                chunksize = length
            # each chunk stores all dimensions (otherwise, users
            # should just create multiple tracks)
            chunkshape = (chunksize,ndim)
            h5file.createCArray(self.hdf_group, ref, 
                                atom=atom,
                                shape=(length,ndim),
                                filters=self.h5_filters,
                                chunkshape=chunkshape)
        # save attributes
        self.hdf_group._v_attrs[DTYPE_ATTR] = dtype
        self.hdf_group._v_attrs[CHANNELS_ATTR] = ndim

    def _get_dtype(self):
        return self.hdf_group._v_attrs[DTYPE_ATTR]
    def _get_nchannels(self):
        return self.hdf_group._v_attrs[CHANNELS_ATTR]    
    def _get_array(self, ref):
        return self.hdf_group._f_getChild(ref)    
    def _get_arrays(self):
        myarrays = {}
        for rname in self.get_rnames():
            myarrays[rname] = self._get_array(rname)            
        return myarrays

    def _check_bounds(self, arr, start, end, channel=None):
        if channel is not None:
            if isinstance(channel, slice):
                pts = (channel.start, channel.stop)
            else:
                pts = (channel,)
            for x in pts:
                if ((x < 0) or (x > arr.shape[1])):
                    raise TrackError("dimension out of bounds '%d'" % (channel))        
        if ((start < 0) or (end > arr.shape[0])):
            raise TrackError("index out of bounds '(%d,%d)'" % (start,end))

    def __getitem__(self, key):
        ref, start, end, strand = self._parse_interval(key)
        arr = self._get_array(ref)
        return arr[start:end]
    
    def __setitem__(self, key, value):
        ref, start, end, strand = self._parse_interval(key)
        arr = self._get_array(ref)
        arr[start:end] = value
                                       