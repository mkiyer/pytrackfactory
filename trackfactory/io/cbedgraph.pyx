'''
Created on Sep 30, 2010

@author: mkiyer
'''
import numpy as np

cimport numpy as np
cimport libc.stdio as stdio
cimport cython

np.import_array()

cdef extern from "Python.h":
    stdio.FILE* PyFile_AsFile(object)

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t
DEFAULT_CHUNKSIZE = int(1e7)
cdef char * PRINT_FORMAT_STRING = "%s\t%d\t%d\t%.2f\n"

def _array_to_bedgraph(char * ref,
                       int start,
                       np.ndarray[DTYPE_t, ndim=1, cast=True] arr,
                       object fileh,
                       float factor):
    cdef int i = 0
    cdef int end = start
    cdef stdio.FILE * fp
    cdef float cov = 0
    cdef float n = 0
    fp = PyFile_AsFile(fileh)

    for i in xrange(0, arr.shape[0]):
        n = float(arr[i])
        if (cov != n):
            if (start != end) and (cov > 0):
                stdio.fprintf(fp, PRINT_FORMAT_STRING, ref, start, end, factor * cov)
            start = end
            cov = n
        end += 1
    if (start != end):
        stdio.fprintf(fp, PRINT_FORMAT_STRING, ref, start, end, factor * cov)

def _array_to_bedgraph_span(char * ref,
                            int offset,
                            np.ndarray[DTYPE_t, ndim=1, cast=True] arr,
                            object fileh,
                            int span,
                            float factor):
    cdef int length = arr.shape[0]
    cdef stdio.FILE * fp
    fp = PyFile_AsFile(fileh)
    cdef float n = 0
    cdef float cov = 0
    cdef int start = 0
    cdef int i = 0

    while (i + span) <= length:
        vals = arr[i:i+span]
        n = np.average(vals) if np.any(vals) else 0
        if (cov != n):       
            if (start != i) and (cov > 0):
                stdio.fprintf(fp, PRINT_FORMAT_STRING, ref, offset + start, offset + i, factor * cov)
            cov = n
            start = i
        i += span
    # cleanup current interval
    if (start != i):
        stdio.fprintf(fp, PRINT_FORMAT_STRING, ref, offset + start, offset + i, factor * cov)
    # cleanup final interval
    if (i < length):
        vals = arr[i:length]
        cov = np.average(vals) if np.any(vals) else 0
        stdio.fprintf(fp, PRINT_FORMAT_STRING, ref, offset + i, offset + length, factor * cov)

def array_to_bedgraph(char * ref, object arr, object fileh, int span=1, float factor=1.0):
    cdef int chunksize
    cdef int start = 0
    span = max(1, span)
    span = min(arr.shape[0], span)
    chunksize = min(arr.shape[0], DEFAULT_CHUNKSIZE)

    while start < (arr.shape[0] - chunksize):
        if span == 1:
            _array_to_bedgraph(ref, start, arr[start:start+chunksize], fileh, factor)
        else:
            _array_to_bedgraph_span(ref, start, arr[start:start+chunksize], fileh, span, factor)            
        start += chunksize

    if start < arr.shape[0]:
        if span == 1:
            _array_to_bedgraph(ref, start, arr[start:arr.shape[0]], fileh, factor)
        else:
            _array_to_bedgraph_span(ref, start, arr[start:start+chunksize], fileh, span, factor)            
