'''
.. module:: bitarraytrack
   :synopsis: track for storing bitarrays
.. moduleauthor:: Matthew Iyer <mkiyer@umich.edu>
.. versionadded:: 0.0.1

Created on Mar 5, 2011

@author: mkiyer
'''
import tables
import numpy as np
from bitarray import bitarray

from track import Track, TrackError

DTYPE_ATTR = 'dtype'
BPB_ATTR = 'bpb'

def translate_bpb(pos, bpb):
    byte, rem = divmod(pos, 8)
    offset, bit = divmod((rem * bpb), 8)
    byte = (byte * bpb) + offset
    return byte,bit

def frame_bpb(start, end, bpb):
    # convert from genomic position to byte/bit position
    # and align to nearest 3bpb frame
    startbyte, startbit = translate_bpb(start, bpb)
    startframebyte = startbyte - (startbyte % bpb)
    startframebit = startbit + 8 * (startbyte - startframebyte)
    # convert from genomic position to byte/bit position
    # end frame offset based on start frame
    endbyte, endbit = translate_bpb(end, bpb)
    endframebit = endbit + 8 * (endbyte - startframebyte)
    endframebyte = endbyte - (endbyte % bpb) + bpb
    return startframebyte, endframebyte, startframebit, endframebit

class BitArrayTrack(Track):
    """Stores genome annotation data in 'bitarray' format to efficiently
    pack data into an array of bytes.  Useful for data that can be 
    represented using less than or equal to 4 bits per base.
    
    :param bpb: number of bits required to store a value 
    """
    h5_complevel = 1
    h5_complib = 'lzo'
    h5_filters = tables.Filters(complevel=h5_complevel, 
                                complib=h5_complib,
                                shuffle=True,
                                fletcher32=False)    
    h5_chunksize = 1<<18 # 256kb chunks
        
    def __init__(self, hdf_group, bpb=4):
        super(BitArrayTrack, self).__init__(hdf_group)
        if not self._check_arrays():
            self._create_arrays(np.dtype(np.byte), bpb)
            hdf_group._v_attrs[DTYPE_ATTR] = np.dtype(np.byte)
            hdf_group._v_attrs[BPB_ATTR] = bpb

    def _bpb(self):
        return self.hdf_group._v_attrs[BPB_ATTR]
    bpb = property(_bpb, None, None, 
                   "bits per base used internally to store data")

    def _check_arrays(self):
        # check if arrays exist
        hasarrays = True
        for rname in self.get_rnames():
            if rname in self.hdf_group:
                continue
            hasarrays = False
            break
        return hasarrays

    def _create_arrays(self, dtype, bpb):        
        # compute chunks sizes tuple (find nearest power of 2)        
        chunksize = 2**int(np.log2(self.h5_chunksize/float(dtype.itemsize)))
        atom = tables.Atom.from_dtype(np.dtype(dtype))
        h5file = self._get_hdf_file()
        # create a separate array for each reference
        for ref, length in self.get_refs():
            if ref in self.hdf_group:
                continue
            # adjust length to reflect compression by 
            # bpb parameter
            bytes,bits = translate_bpb(length, bpb)
            bytes = bytes - (bytes % bpb) + bpb
            # adjust chunksize relative to length in bytes
            if chunksize > bytes:
                chunksize = bytes
            chunkshape = (chunksize,)
            h5file.createCArray(self.hdf_group, ref, 
                                atom=atom,
                                shape=(bytes,),
                                filters=self.h5_filters,
                                chunkshape=chunkshape)

    def _read(self, arr, start, end):
        bpb = self.bpb
        # find 'in frame' start/end positions
        startbyte, endbyte, startbit, endbit = frame_bpb(start, end, bpb)
        # load 'frame' as a bitarray
        frameba = bitarray()
        frameba.fromstring(arr[startbyte:endbyte].tostring())
        # grab the relevant bits from the 'frame' and decode them
        b = frameba[startbit:endbit]
        return b

    def _write(self, arr, start, end, b):
        bpb = self.bpb
        # find 'in frame' start/end positions
        startbyte, endbyte, startbit, endbit = frame_bpb(start, end, bpb) 
        # load 'frame' as a bitarray
        frameba = bitarray()
        frameba.fromstring(arr[startbyte:endbyte].tostring())
        # write the sequence into the bitarray
        frameba[startbit:endbit] = b
        # compress back into the numpy array
        arr[startbyte:endbyte] = np.fromstring(frameba.tostring(), dtype=np.byte)

    def _get_array(self, ref):
        return self.hdf_group._f_getChild(ref)

    def __getitem__(self, key):
        ref, start, end, strand = self._parse_interval(key)
        arr = self._get_array(ref)
        return self._read(arr, start, end)

    def __setitem__(self, key, value):
        ref, start, end, strand = self._parse_interval(key)
        arr = self._get_array(ref)
        self._write(arr, start, end, value)
