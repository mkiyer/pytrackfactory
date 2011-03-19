'''
Created on Mar 6, 2011

@author: mkiyer
'''
import collections
import numpy as np

def write_interval_data_to_array(interval_iter, rname_array_dict, dtype, 
                                 chunksize, num_channels, channel_dict):
    # check params
    if chunksize <= 0:
        chunksize = 1
    chunkshape = (chunksize,num_channels)
    channel_range = range(num_channels)
    # initialize array chunk
    arr = None
    chunk_arr = np.zeros(chunkshape, dtype=dtype)
    chunk_chrom = None
    chunk_start = 0
    chunk_end = 0
    dirty = False
    # keep track of statistics
    interval_dict = collections.defaultdict(lambda: 0)
    total_cov_dict = collections.defaultdict(lambda: 0)
    # parse intervals
    for ival in interval_iter:
        chrom = ival.ref
        start = ival.start
        end = ival.end
        value = ival.value
        if chrom not in rname_array_dict:
            continue
        # ensure value is compatible with array
        value = np.cast[dtype](value)
        # stats
        interval_dict[chrom] += 1
        total_cov_dict[chrom] += ((end - start) * value)
        # check if the new interval is outside the current chunk
        if ((chunk_chrom != chrom) or (start < chunk_start) or 
            (start >= (chunk_start + chunksize))):
            if dirty:
                # write chunk
                arr[chunk_start:chunk_end,:] += chunk_arr[:chunk_end-chunk_start]
                dirty = False
            # reset chunk
            chunk_arr[:] = 0
            chunk_chrom = chrom
            chunk_start = start
            chunk_end = chunk_start
            arr = rname_array_dict[chunk_chrom]

        # deal with intervals larger than the chunk size
        if end > (chunk_start + chunksize):
            # fill up rest of current chunk with value
            begin_size = chunk_start+chunksize-start
            ival.write_array(chunk_arr, (start-chunk_start), chunksize, 0, 
                             begin_size, channel_dict) 
            # do one big write of the rest of the data
            ival.write_array(arr, (chunk_start+chunksize), end,
                             begin_size, end-start, channel_dict)
            # update chunk end
            chunk_end = chunk_start + chunksize
        else:
            if end > chunk_end:
                chunk_end = end                
            # fill array with data
            ival.write_array(chunk_arr, (start-chunk_start), (end-chunk_start), 
                             0, end-start, channel_dict)
        dirty = True
    # write final chunk
    if dirty:        
        arr[chunk_start:chunk_end,:] += chunk_arr[:chunk_end-chunk_start]
    return interval_dict, total_cov_dict

