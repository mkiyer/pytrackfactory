'''
Created on Mar 6, 2011

@author: mkiyer
'''
import numpy as np

from trackfactory.track import POS_STRAND, NEG_STRAND, NO_STRAND

def write_interval_data_to_array(interval_iter, rname_array_dict, dtype, 
                                 channel, chunksize):
    # check params
    if chunksize <= 0:
        chunksize = 1
    # initialize array chunk
    arr = None
    chunk_arr = np.zeros(chunksize, dtype=dtype)
    chunk_chrom = None
    chunk_start = 0
    chunk_end = 0
    dirty = False
    # keep track of statistics
    intervals = 0
    total_cov = 0
    # parse intervals
    for interval in interval_iter:
        chrom, start, end, strand, value = interval[:5]
        if chrom not in rname_array_dict:
            continue
        # stats
        intervals += 1
        total_cov += ((end - start) * value)
        if ((chunk_chrom != chrom) or (start < chunk_start) or 
            (start >= (chunk_start + chunksize))):
            if dirty:
                # write chunk
                arr[chunk_start:chunk_end,channel] += chunk_arr[:chunk_end-chunk_start]
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
            chunk_arr[start-chunk_start:chunksize] += value
            # do one big write of the rest of the data
            arr[chunk_start+chunksize:end,channel] += value
            # update chunk end
            chunk_end = chunk_start + chunksize
        else:
            # fill array with data
            if end > chunk_end:
                chunk_end = end
            chunk_arr[start-chunk_start:end-chunk_start] += value
        dirty = True
    # write final chunk
    if dirty:        
        arr[chunk_start:chunk_end,channel] += chunk_arr[:chunk_end-chunk_start]
    return intervals, total_cov

def write_interval_data_to_stranded_array(interval_iter, rname_array_dict, dtype, chunksize):
    # check params
    if chunksize <= 0:
        chunksize = 1
    # initialize array chunk
    arr = None
    chunk_arr = np.zeros((chunksize,2), dtype=dtype)
    chunk_chrom = None
    chunk_start = 0
    chunk_end = 0
    dirty = False
    # keep track of statistics
    intervals = 0
    total_cov = 0
    # parse intervals
    for interval in interval_iter:
        chrom, start, end, strand, value = interval[:5]
        if chrom not in rname_array_dict:
            continue
        # stats
        intervals += 1
        total_cov += ((end - start) * value)
        # handle strands
        if strand == NO_STRAND:
            strands = (POS_STRAND, NEG_STRAND)
        else:
            strands = (strand,)
        value /= float(len(strands))
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
            for s in strands:
                chunk_arr[start-chunk_start:chunksize,s] += value
                # do one big write of the rest of the data
                arr[chunk_start+chunksize:end,s] += value
            # update chunk end
            chunk_end = chunk_start + chunksize
        else:
            # fill array with data
            if end > chunk_end:
                chunk_end = end                
            for s in strands:
                chunk_arr[start-chunk_start:end-chunk_start,s] += value
        dirty = True
    # write final chunk
    if dirty:        
        arr[chunk_start:chunk_end] += chunk_arr[:chunk_end-chunk_start]
    return intervals, total_cov