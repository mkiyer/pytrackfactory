'''
Created on Mar 6, 2011

@author: mkiyer
'''
import numpy as np
cimport numpy as np
cimport cython
cimport libc.stdio as stdio

def write_interval_data_to_array(interval_iter, rname_array_dict, dtype, 
                                 channel, chunkshape, mode):
    if mode == "channel":
        write_func = write_channel
    elif mode == "strand":
        write_func = write_strand
    elif mode == "allel":
        write_func = write_strand_allele
    # check params
    chunksize = chunkshape[0]
    num_channels = chunkshape[1]
    if chunksize <= 0:
        chunksize = 1
    # initialize array chunk
    arr = None
    chunk_arr = np.zeros(chunkshape, dtype=dtype)
    chunk_chrom = None
    chunk_start = 0
    chunk_end = 0
    dirty = False
    # keep track of statistics
    intervals = 0
    total_cov = 0
    # parse intervals
    for interval in interval_iter:
        chrom, start, end, strand, value, seq = interval[:6]
        if chrom not in rname_array_dict:
            continue
        # stats
        intervals += 1
        total_cov += ((end - start) * value)
        # check if the new interval is outside the current chunk
        if ((chunk_chrom != chrom) or (start < chunk_start) or 
            (start >= (chunk_start + chunksize))):
            if dirty:
                # write chunk
                arr[chunk_start:chunk_end] += chunk_arr[:chunk_end-chunk_start]
                dirty = False
            # reset chunk
            chunk_arr[:] = 0
            chunk_chrom = chrom
            chunk_start = start
            chunk_end = chunk_start
            arr = rname_array_dict[chunk_chrom]
        # ensure value is compatible with array
        value = np.cast[dtype](value)     
        # deal with intervals larger than the chunk size
        if end > (chunk_start + chunksize):
            # fill up rest of current chunk with value
            write_func(chunk_arr, (start-chunk_start), chunksize, channel, 
                       strand, seq[:chunk_start+chunksize-start], value)
            # do one big write of the rest of the data
            write_func(arr, (chunk_start+chunksize), end, channel, 
                       strand, seq[chunk_start+chunksize-start:], value)
            # update chunk end
            chunk_end = chunk_start + chunksize
        else:
            # fill array with data
            if end > chunk_end:
                chunk_end = end
            write_func(chunk_arr, (start-chunk_start), (end-chunk_start), 
                       channel, strand, seq, value)
        dirty = True
    # write final chunk
    if dirty:        
        arr[chunk_start:chunk_end] += chunk_arr[:chunk_end-chunk_start]
    return intervals, total_cov

cdef class Interval:
    """
    Basic interval class
    
    Based on bx-python Interval class by James Taylor
    """
    cdef public int start, end, strand
    cdef public object value, chrom

    def __init__(self, int start, int end, object value=None, object chrom=None, object strand=None ):
        assert start <= end, "start must be less than end"
        self.start  = start
        self.end   = end
        self.value = value
        self.chrom = chrom
        self.strand = strand

    def __repr__(self):
        fstr = "Interval(%d, %d" % (self.start, self.end)
        if not self.value is None:
            fstr += ", value=" + str(self.value)
        fstr += ")"
        return fstr

    def __cmp__(self, other):
        return cmp( self.start, other.start ) or cmp( self.end, other.end )
    
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

