'''
Created on Mar 6, 2011

@author: mkiyer
'''
import numpy as np

def intervals_to_array(interval_iter, dtype, chunksize=(1 << 20)):
    chunksize = max(1, chunksize)
    dtype = np.dtype(dtype)
    arr = np.zeros(chunksize, dtype=dtype)
    chunk_chrom = None
    chunk_start = 0
    chunk_end = chunksize
    max_end = 0
    for chrom, start, end, strand, value in interval_iter:
        #print 'interval', chrom, start, end, strand, value
        #print 'state   ', chunk_chrom, chunk_start, chunk_end, max_end, arr
        if (chunk_chrom == chrom) and (start < chunk_start):
            raise ValueError("(start < chunk_start) intervals must be sorted")
        if (chunk_chrom != chrom) or (chunk_end <= start):
            # reset the chunk
            if chunk_chrom is not None:
                #print "reset"
                yield chunk_chrom, chunk_start, max_end, arr[:(max_end-chunk_start)]
            arr[:] = 0
            chunk_chrom = chrom
            chunk_start = start
            chunk_end = start + chunksize
        # deal with intervals larger than the chunk size
        while end > chunk_end:
            # fill up rest of current chunk with value and return it
            arr[start-chunk_start:chunksize] = np.cast[dtype](value)
            #print 'returning   ', chunk_chrom, chunk_start, chunk_end, arr[:chunksize]
            yield chunk_chrom, chunk_start, chunk_end, arr[:chunksize]
            # advance interval and chunk
            arr[:] = 0
            chunk_start = chunk_end
            chunk_end = chunk_start + chunksize
            start = chunk_start 
            max_end = 0
        # fill array with data
        if end > max_end:
            max_end = end
        arr[start-chunk_start:end-chunk_start] = np.cast[dtype](value)
        #print 'newstate', chunk_chrom, chunk_start, chunk_end, max_end, arr        
    # return the last chunk 
    yield chunk_chrom, chunk_start, max_end, arr[:(max_end-chunk_start)]
