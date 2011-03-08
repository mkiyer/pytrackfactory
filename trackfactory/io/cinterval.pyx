'''
Created on Mar 6, 2011

@author: mkiyer
'''
import numpy as np
cimport numpy as np
cimport cython
cimport libc.stdio as stdio

cdef class IntervalToArrayChunks:
    cdef object chunk_chrom
    cdef int chunk_start
    cdef int chunk_end
    cdef int max_end
    cdef int chunksize
    cdef object arr
    cdef object interval_iter
    cdef bint eof
    cdef list arraystack
    
    def __init__( self, interval_iter, arr):
        self.interval_iter = interval_iter
        self.arr = arr
        if arr.shape[0] <= 0:
            raise "invalid zero-length array"
        self.chunksize = arr.shape[0]
        self.chunk_chrom = None
        self.chunk_start = 0
        self.chunk_end = self.chunksize
        self.max_end = 0
        self.eof = False
        self.arraystack = []

    def __iter__( self ):
        return self

    def __next__( self ):
        cdef object oldarr
        cdef object oldchrom
        cdef int oldstart
        cdef int oldend
        if self.eof:
            raise StopIteration()        
        # read an interval
        while True:
            # if there are arrays to return, return them here
            if len(self.arraystack) > 0:
                return self.arraystack.pop(0)
            # attempt to read next interval
            try:
                chrom, start, end, strand, value = self.interval_iter.next()
            except StopIteration:
                # set EOF flag
                self.eof = True
                # return final array if necessary
                if (self.max_end > self.chunk_start):
                    oldend = self.max_end
                    self.max_end = self.chunk_start
                    return self.chunk_chrom, self.chunk_start, oldend, self.arr[:(oldend - self.chunk_start)]
            # process interval
            if (self.chunk_chrom != chrom) or (self.chunk_end <= start):
                # see if need to return anything
                if self.chunk_chrom is not None:                    
                    # copy state
                    oldarr = self.arr[:(self.max_end - self.chunk_start)].copy()
                    oldchrom = self.chunk_chrom
                    oldstart = self.chunk_start
                    #oldend = self.chunk_end
                    oldend = self.max_end
                    self.arraystack.append((oldchrom, oldstart, oldend, oldarr))
                # reset state
                self.arr[:] = 0
                self.chunk_chrom = chrom
                self.chunk_start = start
                self.chunk_end = start + self.chunksize
                self.max_end = self.chunk_start
            # deal with intervals larger than the chunk size
            while (end > self.chunk_end):
                # fill up rest of current chunk with value
                self.arr[start-self.chunk_start:self.chunksize] = np.cast[self.arr.dtype](value)
                # save values to return
                oldarr = self.arr[:self.chunksize].copy()
                oldchrom = self.chunk_chrom
                oldstart = self.chunk_start
                oldend = self.chunk_end
                self.arraystack.append((oldchrom, oldstart, oldend, oldarr))            
                # set new state values
                self.arr[:] = 0
                self.chunk_start = self.chunk_end
                self.chunk_end = self.chunk_start + self.chunksize
                start = self.chunk_start 
                self.max_end = self.chunk_start
            # fill array with data
            if end > self.max_end:
                self.max_end = end
            self.arr[start-self.chunk_start:end-self.chunk_start] = np.cast[self.arr.dtype](value)

