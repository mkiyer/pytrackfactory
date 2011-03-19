'''
Created on Mar 6, 2011

@author: mkiyer
'''
import numpy as np
cimport numpy as np
cimport cython
cimport libc.stdio as stdio

from cinterval import Interval, SequenceInterval, BedInterval

# constants
DEF POS_STRAND = 0
DEF NEG_STRAND = 1
DEF NO_STRAND = 2

_strand_str_to_int = {"+": POS_STRAND,
                      "-": NEG_STRAND,
                      ".": NO_STRAND}
_strand_int_to_str = dict([(v,k) for k,v in _strand_str_to_int.items()])
def strand_str_to_int(strand):
    return _strand_str_to_int.get(strand, NO_STRAND)
def strand_int_to_str(strand):
    return _strand_int_to_str.get(strand, NO_STRAND)


cdef class Interval:
    """
    Basic interval class
    
    Based on bx-python Interval class by James Taylor
    """
    def __init__(self, object ref, int start, int end, int strand=NO_STRAND,
                 object value=None, int readnum=-1):                 
        assert start <= end, "start must be less than end"
        self.start = start
        self.end = end
        self.value = value
        self.ref = ref
        self.strand = strand
        self.readnum = readnum

    def write_array(self, object dst, int start, int end, int istart, 
                    int iend, dict channel_dict):
        channels = channel_dict[(self.readnum, self.strand, None)]
        dst[start:end,channels] += (self.value / float(len(channels)))

    def __repr__(self):
        fstr = "Interval(%d, %d" % (self.start, self.end)
        if not self.value is None:
            fstr += ", value=" + str(self.value)
        fstr += ")"
        return fstr

cdef class SequenceInterval(Interval):
    """Interval with sequence information
    """
    def __init__(self, object ref, int start, int end, int strand=NO_STRAND,
                 object value=None, object seq=None, int readnum=-1):
        Interval.__init__(self, ref, start, end, strand, value, readnum)
        self.seq = seq

    def write_array(self, object dst, int start, int end, int istart, 
                    int iend, dict channel_dict):
        for i in xrange(istart, iend):
            dna = self.seq[i]
            channels = channel_dict[(self.readnum, self.strand, dna)]
            dst[start+i,channels] += (self.value / float(len(channels)))

    def __repr__(self):
        fstr = "Interval(%d, %d" % (self.start, self.end)
        if not self.value is None:
            fstr += ", value=" + str(self.value)
        if not self.seq is None:
            fstr += ", seq=" + str(self.seq)
        fstr += ")"
        return fstr

cdef class BedInterval(Interval):
    """Interval with name information
    """
    def __init__(self, object ref, int start, int end, int strand=NO_STRAND,
                 object value=None, int readnum=-1, object name=None):
        Interval.__init__(self, ref, start, end, strand, value, readnum)
        self.name = name

    def __repr__(self):
        fstr = "Interval(%d, %d" % (self.start, self.end)
        if not self.value is None:
            fstr += ", value=" + str(self.value)
        if not self.name is None:
            fstr += ", name=" + str(self.name)
        fstr += ")"
        return fstr

