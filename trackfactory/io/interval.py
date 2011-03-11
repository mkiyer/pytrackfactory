'''
Created on Mar 6, 2011

@author: mkiyer
'''
import numpy as np

from trackfactory.track import POS_STRAND, NEG_STRAND, NO_STRAND

def write_channel(dst, start, end, channel, strand, seq, value):
    dst[start:end,channel] += value

def write_strand(dst, start, end, channel, strand, seq, value):
    if strand == NO_STRAND: 
        channels = (0, 1)
    else: 
        channels = (strand,)
    dst[start:end,channels] += (value / float(len(channels)))

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
        # check if the new interval is outside the current chunk
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


_base_channel_dict = {'a': 0,
                      'g': 1,
                      'c': 2,
                      't': 3,
                      'A': 0,
                      'G': 1,
                      'C': 2,
                      'T': 3}

def make_allele_array(start, end, strand, value, seq, dtype):
    valarray = np.zeros((end-start,8), dtype=dtype)
    for pos in xrange(end - start):
        base = seq[pos]
        channel = _base_channel_dict.get(base, None)
        if (channel is None):
            if strand == NO_STRAND:
                valarray[pos,:] = value / 8.0
            elif strand == POS_STRAND:
                valarray[pos,:4] = value / 4.0
            else:
                valarray[pos,4:8] = value / 4.0                   
        else:
            if strand == NO_STRAND:
                valarray[pos,channel] = value / 2.0
                valarray[pos,channel+4] = value / 2.0                    
            elif strand == POS_STRAND:
                valarray[pos,channel] = value
            else:
                valarray[pos,channel+4] = value
    return valarray
 
_channel_lookup_dict = {(POS_STRAND, "A"): (0,),
                        (POS_STRAND, "G"): (1,),
                        (POS_STRAND, "C"): (2,),                      
                        (POS_STRAND, "T"): (3,),
                        (POS_STRAND, "N"): (0,1,2,3),
                        (NEG_STRAND, "A"): (4,),
                        (NEG_STRAND, "G"): (5,),
                        (NEG_STRAND, "C"): (6,),                      
                        (NEG_STRAND, "T"): (7,),
                        (POS_STRAND, "N"): (4,5,6,7),
                        (NO_STRAND, "A"): (0,4),
                        (NO_STRAND, "G"): (1,5),
                        (NO_STRAND, "C"): (2,6),                    
                        (NO_STRAND, "T"): (3,7),
                        (NO_STRAND, "N"): (0,1,2,3,4,5,6,7)}

def write_strand_allele(dst, start, end, channel, strand, seq, value):
    for i,base in enumerate(seq):
        channels = _channel_lookup_dict[base]
        dst[start+i,channels] += (value / float(len(channels)))

def write_interval_data_to_stranded_allele_array(interval_iter, rname_array_dict, dtype, chunksize):
    # check params
    if chunksize <= 0:
        chunksize = 1
    # initialize array chunk
    arr = None
    chunk_arr = np.zeros((chunksize,8), dtype=dtype)
    chunk_chrom = None
    chunk_start = 0
    chunk_end = 0
    dirty = False
    # keep track of statistics
    intervals = 0
    total_cov = 0
    # parse intervals
    for interval in interval_iter:
        chrom, start, end, strand, value, seq = interval
        if chrom not in rname_array_dict:
            continue
        # stats
        intervals += 1
        total_cov += ((end - start) * value)
        # handle alleles
        valarray = make_allele_array(start, end, strand, value, seq, dtype)
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
        # deal with intervals larger than the chunk size
        if end > (chunk_start + chunksize):
            # fill up rest of current chunk with value
            # (uses broadcasting to get all alleles)
            begin_size = chunk_start+chunksize-start
            #begin_size = (chunksize - (start - chunk_start))
            chunk_arr[start-chunk_start:chunksize] += valarray[:begin_size]
            # do one big write of the rest of the data
            arr[chunk_start+chunksize:end] += valarray[begin_size:]
            # update chunk end
            chunk_end = chunk_start + chunksize
        else:
            if end > chunk_end:
                chunk_end = end                
            # fill array with data
            # (uses broadcasting to add alleles)
            chunk_arr[start-chunk_start:end-chunk_start] += valarray
        dirty = True
    # write final chunk
    if dirty:        
        arr[chunk_start:chunk_end] += chunk_arr[:chunk_end-chunk_start]
    return intervals, total_cov
