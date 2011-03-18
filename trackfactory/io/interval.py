'''
Created on Mar 6, 2011

@author: mkiyer
'''
import collections
import numpy as np

from trackfactory.track import POS_STRAND, NEG_STRAND, NO_STRAND

def write_channel(dst, start, end, interval, istart, iend, 
                  channels=(0,)):
    dst[start:end,channels] += interval.value

def write_strand(dst, start, end, interval, istart, iend, 
                 channels=(0,)):
    strand = interval.strand
    if strand == NO_STRAND: 
        channels = (0, 1)
    else: 
        channels = (strand,)
    dst[start:end,channels] += (interval.value / float(len(channels)))

_channel_lookup_dict = {(POS_STRAND, "A"): (0,),
                        (POS_STRAND, "G"): (1,),
                        (POS_STRAND, "C"): (2,),                      
                        (POS_STRAND, "T"): (3,),
                        (POS_STRAND, "N"): (0,1,2,3),
                        (NEG_STRAND, "A"): (4,),
                        (NEG_STRAND, "G"): (5,),
                        (NEG_STRAND, "C"): (6,),                      
                        (NEG_STRAND, "T"): (7,),
                        (NEG_STRAND, "N"): (4,5,6,7),
                        (NO_STRAND, "A"): (0,4),
                        (NO_STRAND, "G"): (1,5),
                        (NO_STRAND, "C"): (2,6),                    
                        (NO_STRAND, "T"): (3,7),
                        (NO_STRAND, "N"): (0,1,2,3,4,5,6,7)}

def write_strand_allele(dst, start, end, interval, istart, iend, 
                        channels=(0,)):
    for i in xrange(istart, iend):
        base = interval.seq[i]
        channels = _channel_lookup_dict[(interval.strand, base)]
        dst[start+i,channels] += (interval.value / float(len(channels)))

def write_interval_data_to_array(interval_iter, rname_array_dict, dtype, 
                                 chunksize, mode, channel=0):
    # check params
    if chunksize <= 0:
        chunksize = 1
    if mode == "channel":
        write_func = write_channel
        chunkshape = (chunksize,1)
        channels = (channel,)
    elif mode == "strand":
        write_func = write_strand
        chunkshape = (chunksize,2)
        channels = (0,1)
    elif mode == "allele":
        write_func = write_strand_allele
        chunkshape = (chunksize,8)
        channels = (0,1,2,3,4,5,6,7)
    else:
        raise "invalid mode '%s'" % (mode)
    # initialize array chunk
    arr = None
    chunk_arr = np.zeros(chunkshape, dtype=dtype)
    chunk_chrom = None
    chunk_start = 0
    chunk_end = 0
    dirty = False
    # keep track of statistics
    intervals = collections.defaultdict(lambda: 0)
    total_cov = collections.defaultdict(lambda: 0)
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
        intervals[chrom] += 1
        total_cov[chrom] += ((end - start) * value)
        # check if the new interval is outside the current chunk
        if ((chunk_chrom != chrom) or (start < chunk_start) or 
            (start >= (chunk_start + chunksize))):
            if dirty:
                # write chunk
                arr[chunk_start:chunk_end,channels] += chunk_arr[:chunk_end-chunk_start]
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
            write_func(chunk_arr, (start-chunk_start), chunksize, 
                       ival, 0, begin_size) 
            # do one big write of the rest of the data
            write_func(arr, (chunk_start+chunksize), end, ival,
                       begin_size, end-start) 
            # update chunk end
            chunk_end = chunk_start + chunksize
        else:
            if end > chunk_end:
                chunk_end = end                
            # fill array with data
            write_func(chunk_arr, (start-chunk_start), (end-chunk_start), 
                       ival, 0, end-start)
        dirty = True
    # write final chunk
    if dirty:        
        arr[chunk_start:chunk_end,channels] += chunk_arr[:chunk_end-chunk_start]
    return intervals, total_cov
