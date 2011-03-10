'''
Created on Sep 29, 2010

@author: mkiyer
'''
import numpy as np
import collections

DEFAULT_CHUNKSIZE = (1 << 20)
PRINT_FORMAT_STRING = "%s\t%d\t%d\t%.2f"

def _array_to_bedgraph(ref, start, arr, channels, fileh, factor=1.0):
    end = start
    cov = 0
    for i in xrange(0, arr.shape[0]):
        n = arr[i,channels].sum()
        if (cov != n):
            if (start != end) and (cov != 0):
                print >>fileh, PRINT_FORMAT_STRING % (ref, start, end, factor * cov)
            start = end
            cov = n
        end += 1
    if (start != end):
        print >>fileh, PRINT_FORMAT_STRING % (ref, start, end, factor * cov)

def _array_to_bedgraph_span(ref, offset, arr, channels, fileh, span, factor=1.0):
    length = arr.shape[0]
    if span < 1: span = 1
    if span > length: span = length    
    cov = 0
    i = 0
    start = 0
    while (i + span) <= length:
        vals = arr[i:i+span,channels]
        if len(vals.shape) > 1:
            vals = vals.sum(axis=1)
        n = np.average(vals) if np.any(vals) else 0
        if (cov != n):       
            if (start != i) and (cov != 0):
                print >>fileh, PRINT_FORMAT_STRING % (ref, offset + start, offset + i, factor * cov)
            cov = n
            start = i
        i += span
    # cleanup current interval
    if (start != i):
        print >>fileh, PRINT_FORMAT_STRING % (ref, offset + start, offset + i, factor * cov)
    # cleanup final interval
    if (i < length):
        vals = arr[i:length]
        cov = np.average(vals) if np.any(vals) else 0
        print >>fileh, PRINT_FORMAT_STRING % (ref, offset + i, offset + length, factor * cov)

def array_to_bedgraph(ref, arr, fileh, start=0, end=-1, span=1, factor=1.0, 
                      chunksize=DEFAULT_CHUNKSIZE, channels=None):
    if end < start:
        end = arr.shape[0]
    if span < 1:
        span = 1
    if span > (end - start):
        span = (end - start)
    if chunksize > (end - start):
        chunksize = (end - start)
    if channels is None:
        channels = np.s_[0:arr.shape[1]]
    while start < (end - chunksize):
        if span == 1:
            _array_to_bedgraph(ref, start, arr[start:start+chunksize], channels, fileh, factor)
        else:
            _array_to_bedgraph_span(ref, start, arr[start:start+chunksize], channels, fileh, span, factor)            
        start += chunksize

    if start < end:
        if span == 1:
            _array_to_bedgraph(ref, start, arr[start:end], channels, fileh, factor)
        else:
            _array_to_bedgraph_span(ref, start, arr[start:start+chunksize], channels, fileh, span, factor)


def bedgraph_to_array(fileh, ref=None):
    intervals = collections.defaultdict(lambda: [])
    maxend = collections.defaultdict(lambda: 0)
    for line in fileh:
        fields = line.strip().split('\t')
        ref, start, end, cov = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
        intervals[ref].append((start, end, cov))
        maxend[ref] = max(maxend[ref], end)
    covarrays = {}
    for ref, values in intervals.iteritems():
        covarrays[ref] = np.zeros(maxend[ref], dtype=np.float)
        for start, end, cov in values:  
            covarrays[ref][start:end] += cov
    return covarrays

    
    