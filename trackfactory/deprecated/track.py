'''
Created on Mar 9, 2011

@author: mkiyer
'''
#def parse_interval(interval):
#    """parses a genomic interval specifier in either the string
#    format `<ref:start-end>` or tuple (ref, start, end)
#    
#    :param interval: interval specifier
#    :type interval: str or tuple
#    
#    >>> parse_interval("chr1:100-200")
#    ("chr1", 100, 200)
#    >>> parse_interval("chr1:1,000-20,000")
#    ("chr1", 1000, 20000)
#    >>> parse_interval("chrX")
#    ("chrX", None, None)
#    >>> parse_interval("chrX:5")
#    ("chrX", 5, 6)
#    >>> parse_interval(("chr2", 100, 300))
#    ("chr2", 100, 300) 
#    """
#    if interval is None:
#        return None, None, None
#    if isinstance(interval, basestring):
#        # remove commas
#        k = interval.replace(',','')            
#        # split chrom on ':' character
#        if k.find(':') != -1:
#            ref, startend = k.split(':')                
#            # separate start/end on '-' character
#            if k.find('-') != -1:
#                start, end = map(int, startend.split('-'))
#            else:
#                start = int(startend)
#                end = start + 1
#        else:
#            ref = k
#            start = None
#            end = None
#        return ref, start, end
#    else:
#        if len(interval) == 1:
#            return interval[0], None, None
#        elif len(interval) == 2:
#            return interval[0], interval[1], interval[1]+1
#        elif len(interval) >= 3:
#            return interval[:3]