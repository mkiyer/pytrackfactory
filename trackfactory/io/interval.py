'''
Created on Mar 6, 2011

@author: mkiyer
'''
import collections
import logging
import numpy as np

from cinterval import SequenceInterval, Interval
from trackfactory.track import TrackError

def write_array(ival, istart, iend, dst, start, end, channel_dict):
    channels = channel_dict[(ival.readnum, ival.strand, None)]
    dst[start:end,channels] += (ival.value / float(len(channels)))

def write_allele_array(ival, istart, iend, dst, start, end, 
                       channel_dict):
    for i in xrange(istart, iend):
        dna = ival.seq[i]
        channels = channel_dict[(ival.readnum, ival.strand, dna)]
        dst[start+i,channels] += (ival.value / float(len(channels)))

def write_interval_data_to_array(interval_iter, rname_array_dict, dtype, 
                                 chunksize, num_channels, channel_dict,
                                 allele_writer=False):
    # check params
    if chunksize <= 0:
        chunksize = 1
    chunkshape = (chunksize,num_channels)
    if allele_writer:
        write_func = write_allele_array
    else:
        write_func = write_array
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
    # debugging
    debug_every = 1000000
    debug_next = debug_every
    debug_count = 0
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
            write_func(ival, 0, begin_size, chunk_arr, (start-chunk_start), 
                       chunksize, channel_dict) 
            # do one big write of the rest of the data
            write_func(ival, begin_size, end-start, arr, 
                       (chunk_start+chunksize), end, channel_dict)
            # update chunk end
            chunk_end = chunk_start + chunksize
        else:
            if end > chunk_end:
                chunk_end = end                
            # fill array with data
            write_func(ival, 0, end-start, chunk_arr, (start-chunk_start), 
                       (end-chunk_start), channel_dict)
        dirty = True
        # debugging
        if debug_count == debug_next:
            logging.debug("\tWrote %d intervals" % (debug_count))
            debug_next += debug_every
        debug_count += 1
    # write final chunk
    if dirty:        
        arr[chunk_start:chunk_end,:] += chunk_arr[:chunk_end-chunk_start]
    return interval_dict, total_cov_dict
