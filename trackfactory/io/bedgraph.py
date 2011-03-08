'''
Created on Sep 29, 2010

@author: mkiyer
'''
import numpy as np
import collections

def array_to_bedgraph(ref, array_iter, fileh, factor=1.0):
    start, end, cov = 0, 0, 0
    for n in array_iter:
        if (start != end) and (cov != n):
            if cov > 0:
                print >>fileh, "%s\t%d\t%d\t%.2f" % (ref, start, end, factor * cov)
            start = end
            cov = n
        end += 1
    if (start != end):
        print >>fileh, "%s\t%d\t%d\t%.2f" % (ref, start, end, factor * cov)

#def _array_to_bedgraph(ref, arr, fileh, factor=1.0):
#    start, end, cov = 0, 0, 0
#    for n in arr:
#        if (start != end) and (cov != n):
#            if cov > 0:
#                print >>fileh, "%s\t%d\t%d\t%.2f" % (ref, start, end, factor * cov)
#            start = end
#            cov = n
#        end += 1
#    if (start != end):
#        print >>fileh, "%s\t%d\t%d\t%.2f" % (ref, start, end, factor * cov)
#
#def _array_to_bedgraph_span(ref, arr, fileh, span, factor=1.0):
#    length = arr.shape[0]
#    span = max(1, span)
#    span = min(length, span)
#    wrote_something = False
#    start = 0
#    while (start + span) <= length:
#        vals = arr[start:start+span]
#        if np.any(vals):                
#            print >>fileh, "%s\t%d\t%d\t%.2f" % (ref, start, start+span, factor * np.average(vals))
#            wrote_something = True    
#        start += span
#    # cleanup final interval
#    if start < length:
#        vals = arr[start:length]
#        if np.any(vals):                
#            print >>fileh, "%s\t%d\t%d\t%.2f" % (ref, start, length, factor * np.average(vals))
#            wrote_something = True    
#    # make sure to write something
#    if not wrote_something:
#        print >>fileh, "%s\t%d\t%d\t%.2f" % (ref, 0, length, 0)
#
#def array_to_bedgraph(ref, arr, fileh, span=1, factor=1.0):
#    if span == 1:
#        _array_to_bedgraph(ref, arr, fileh, factor)
#    else:
#        _array_to_bedgraph_span(ref, arr, fileh, span, factor)
#
#def bedgraph_to_array(fileh, ref=None):
#    intervals = collections.defaultdict(lambda: [])
#    maxend = collections.defaultdict(lambda: 0)
#    for line in fileh:
#        fields = line.strip().split('\t')
#        ref, start, end, cov = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
#        intervals[ref].append((start, end, cov))
#        maxend[ref] = max(maxend[ref], end)
#    covarrays = {}
#    for ref, values in intervals.iteritems():
#        covarrays[ref] = np.zeros(maxend[ref], dtype=np.float)
#        for start, end, cov in values:  
#            covarrays[ref][start:end] += cov
#    return covarrays

    
    