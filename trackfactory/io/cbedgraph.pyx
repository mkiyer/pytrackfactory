'''
Created on Mar 10, 2011

@author: mkiyer
'''
import numpy as np

cimport numpy as np

DEFAULT_CHUNKSIZE = (1 << 20)
PRINT_FORMAT_STRING = "%s\t%d\t%d\t%.2f"

def _array_to_bedgraph(char * ref, int start, 
                       np.ndarray arr, 
                       object channels, 
                       object fileh, 
                       float factor=1.0):
    cdef int end = start
    cdef float cov = 0
    cdef float n    
    for i in xrange(0, arr.shape[0]):
        n = arr[i,channels].sum()
        if (cov != n):
            if (start != end) and (cov != 0):
                print >>fileh, PRINT_FORMAT_STRING % (ref, start, end, factor * cov)
            start = end
            cov = n
        end += 1
    if (start != end):
        print >>fileh, PRINT_FORMAT_STRING % (ref, start, end, factor * cov)

def _array_to_bedgraph_span(char * ref, int offset, 
                            np.ndarray arr, 
                            object channels, 
                            object fileh, 
                            int span, 
                            float factor=1.0):
    cdef float cov = 0
    cdef int i = 0
    cdef int start = 0
    cdef int length = arr.shape[0]
    cdef np.ndarray vals
    cdef float n
    if span < 1: span = 1
    if span > length: span = length    
    while (i + span) <= length:
        vals = arr[i:i+span,channels]
        if len((<object>vals).shape) > 1:
            vals = vals.sum(axis=1)
        n = np.average(vals) if np.any(vals) else 0
        if (cov != n):       
            if (start != i) and (cov != 0):
                print >>fileh, PRINT_FORMAT_STRING % (ref, offset + start, offset + i, factor * cov)
            cov = n
            start = i
        i += span
    # cleanup current interval
    if (start != i):
        print >>fileh, PRINT_FORMAT_STRING % (ref, offset + start, offset + i, factor * cov)
    # cleanup final interval
    if (i < length):
        vals = arr[i:length]
        cov = np.average(vals) if np.any(vals) else 0
        print >>fileh, PRINT_FORMAT_STRING % (ref, offset + i, offset + length, factor * cov)


def array_to_bedgraph(char * ref, object arr, object fileh, int start=0, 
                      int end=-1, int span=1, float factor=1.0, 
                      int chunksize=DEFAULT_CHUNKSIZE, 
                      object channels=None):
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
        if span == 1:
            _array_to_bedgraph(ref, start, arr[start:start+chunksize], channels, fileh, factor)
        else:
            _array_to_bedgraph_span(ref, start, arr[start:start+chunksize], channels, fileh, span, factor)            
        start += chunksize

    if start < end:
        if span == 1:
            _array_to_bedgraph(ref, start, arr[start:end], channels, fileh, factor)
        else:
            _array_to_bedgraph_span(ref, start, arr[start:start+chunksize], channels, fileh, span, factor)


    
    