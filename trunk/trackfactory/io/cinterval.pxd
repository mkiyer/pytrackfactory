'''
Created on Mar 11, 2011

@author: mkiyer
'''

cdef class Interval:
    """Basic interval class

    Based on bx-python Interval class by James Taylor
    """
    cdef public object ref
    cdef public int start
    cdef public int end 
    cdef public int strand
    cdef public object value
    cdef public int readnum

cdef class BedInterval(Interval):
    cdef public object name

cdef class SequenceInterval(Interval):
    cdef public object seq
