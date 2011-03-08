'''
Created on Mar 6, 2011

@author: mkiyer
'''
import numpy as np

def write_interval_data_to_array(interval_iter, arr, chunksize):
    # check params
    if chunksize <= 0:
        chunksize = 1
    # initialize array chunk
    chunk_arr = np.zeros(chunksize, dtype=arr.dtype)    
    chunk_start = 0
    chunk_end = 0
    dirty = False
    # parse intervals
    for start, end, strand, value, seq in interval_iter:
        #print "coord", start, end, value, "chunk", chunk_start, chunk_end
        if (start < chunk_start) or (start >= (chunk_start + chunksize)):
            if dirty:
                # write chunk
                #print "writing", chunk_start, chunk_end
                arr[chunk_start:chunk_end] += chunk_arr[:chunk_end-chunk_start]
                dirty = False
            # reset chunk
            chunk_arr[:] = 0
            chunk_start = start
            chunk_end = chunk_start
        # ensure value is compatible with array
        value = np.cast[arr.dtype](value)     
        # deal with intervals larger than the chunk size
        if end > (chunk_start + chunksize):
            # fill up rest of current chunk with value
            chunk_arr[start-chunk_start:chunksize] += value
            # do one big write of the rest of the data
            arr[chunk_start+chunksize:end] += value
            # update chunk end
            chunk_end = chunk_start + chunksize
        else:
            # fill array with data
            if end > chunk_end:
                chunk_end = end
            chunk_arr[start-chunk_start:end-chunk_start] += value
        dirty = True
        #print 'newstate', chunk_start, chunk_end
        #debug += 1
        #if debug % 100000 == 0:
        #    print "writing", debug, chunk_start, chunk_end
    # write final chunk
    if dirty:
        arr[chunk_start:chunk_end] += chunk_arr[:chunk_end-chunk_start]


def intervals_to_array(interval_iter, dtype, chunksize=(1 << 20)):
    # check params
    if chunksize <= 0:
        chunksize = 1
    # must be a valid dtype
    dtype = np.dtype(dtype)
    # initialize array chunk
    arr = np.zeros(chunksize, dtype=dtype)
    chunk_chrom = None
    chunk_start = 0
    chunk_end = chunksize
    max_end = 0
    # parse intervals
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
            max_end = chunk_start
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
            max_end = chunk_start
        # fill array with data
        if end > max_end:
            max_end = end
        arr[start-chunk_start:end-chunk_start] = np.cast[dtype](value)
        #print 'newstate', chunk_chrom, chunk_start, chunk_end, max_end, arr        
    # return the last chunk 
    yield chunk_chrom, chunk_start, max_end, arr[:(max_end-chunk_start)]
