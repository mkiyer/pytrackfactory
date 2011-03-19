'''
Created on Mar 10, 2011

@author: mkiyer
'''
import numpy as np
cimport numpy as np

DEFAULT_CHUNKSIZE = (1 << 20)
FORMAT_STRING = "%s\t%d\t%d\t"
FORMAT_STRING_INT = FORMAT_STRING + "%d"
FORMAT_STRING_FLOAT = FORMAT_STRING + "%.6f"

cdef _array_to_bedgraph(char * ref, int start, 
                        np.ndarray arr, 
                        object fileh, 
                        char * fmt_string,
                        float factor=1.0):
    cdef int end = start
    cdef float cov = 0
    cdef float n
    cdef int i    
    for i in xrange(0, arr.shape[0]):
        n = arr[i]
        if (cov != n):
            if (start != end) and (cov != 0):
                print >>fileh, fmt_string % (ref, start, end, factor * cov)
            start = end
            cov = n
        end += 1
    if (start != end):
        print >>fileh, fmt_string % (ref, start, end, factor * cov)

cdef _array_to_bedgraph_span(char * ref, int offset, 
                             np.ndarray arr, 
                             object fileh, 
                             char * fmt_string,
                             float factor=1.0,
                             int span=1):                            
    cdef float cov = 0
    cdef int i = 0
    cdef int start = 0
    cdef int length = arr.shape[0]
    cdef np.ndarray vals
    cdef float n
    if span < 1: span = 1
    if span > length: span = length    
    while (i + span) <= length:
        vals = arr[i:i+span]
        #if len((<object>vals).shape) > 1:
        #    vals = vals.sum(axis=1)
        n = np.average(vals) if np.any(vals) else 0
        if (cov != n):       
            if (start != i) and (cov != 0):
                print >>fileh, fmt_string % (ref, offset + start, offset + i, factor * cov)
            cov = n
            start = i
        i += span
    # cleanup current interval
    if (start != i):
        print >>fileh, fmt_string % (ref, offset + start, offset + i, factor * cov)
    # cleanup final interval
    if (i < length):
        vals = arr[i:length]
        cov = np.average(vals) if np.any(vals) else 0
        print >>fileh, fmt_string % (ref, offset + i, offset + length, factor * cov)


def array_to_bedgraph(char * ref, object arr, object fileh, int start=0, 
                      int end=-1, float factor=1.0, int span=1,  
                      int chunksize=DEFAULT_CHUNKSIZE, 
                      object channels=None):
    cdef char * fmt_string
    cdef np.ndarray chunk_arr
   
    if arr.dtype.kind == "i":
        fmt_string = FORMAT_STRING_INT
    elif arr.dtype.kind == "f":
        fmt_string = FORMAT_STRING_FLOAT
    else:
        raise "invalid array dtype '%s'" % (arr.dtype.kind)
    if end < start:
        end = arr.shape[0]
    if span < 1:
        span = 1
    if span > (end - start):
        span = (end - start)
    if chunksize > (end - start):
        chunksize = (end - start)
    if channels is None:
        channels = np.s_[0:arr.shape[1]]
    while start < (end - chunksize):
        # extract chunk from array 
        chunk_arr = arr[start:start+chunksize,channels].sum(axis=1)
        # write chunk
        if span == 1:
            _array_to_bedgraph(ref, start, chunk_arr, fileh, fmt_string, factor)
        else:
            _array_to_bedgraph_span(ref, start, chunk_arr, fileh, 
                                    fmt_string, factor, span)            
        start += chunksize

    if start < end:
        # extract chunk from array 
        chunk_arr = arr[start:end,channels].sum(axis=1)
        if span == 1:
            _array_to_bedgraph(ref, start, chunk_arr, fileh, fmt_string, 
                               factor)
        else:
            _array_to_bedgraph_span(ref, start, chunk_arr, fileh, 
                                    fmt_string, factor, span)