"""
.. module:: cwiggle
   :synopsis: parser for wiggle format
.. moduleauthor:: Matthew Iyer <mkiyer@umich.edu>
.. versionadded:: 0.0.1

Adapted from James Taylor's bx-python library

Support for scores in the `wiggle`_ file format used by the UCSC Genome 
Browser.

The positions in the wiggle format are 1-relative, however,
the positions returned match the BED/interval format which is zero-based, half-open.

.. _wiggle: http://genome.ucsc.edu/goldenPath/help/wiggle.html

Created on Mar 5, 2011
"""
from cinterval cimport Interval
from cinterval import strand_str_to_int

# constants
DEF POS_STRAND = 0
DEF NEG_STRAND = 1
DEF NO_STRAND = 2

cdef enum linemode:
    MODE_BED
    MODE_FIXED
    MODE_VARIABLE

def parse_header( line ):
    return dict( [ field.split( '=' ) for field in line.split()[1:] ] )

cdef class WiggleReader:
    """
    Iterator yielding chrom, start, end, strand, value.
    Values are zero-based, half-open.
    Regions which lack a score are ignored.
    """
    cdef object file
    cdef object current_chrom
    cdef long current_pos
    cdef long current_step
    cdef long current_span
    cdef linemode mode

    def __init__( self, file ):
        self.file = file
        self.current_chrom = None
        self.current_pos = -1
        self.current_step = -1
        self.current_span = -1
        self.mode = MODE_BED

    def __iter__( self ):
        return self

    def __next__( self ):
        while True:
            line = self.file.readline()
            if not line:
                raise StopIteration()
            if line.isspace():
                continue    
            if line[0] == "#":
                continue
            if line.startswith( "track" ) or line.startswith( "browser" ):
                continue
            if line.startswith( "variableStep" ):
                header = parse_header( line )
                self.current_chrom = header['chrom']
                self.current_pos = -1
                self.current_step = -1
                if 'span' in header:
                    self.current_span = int( header['span'] )
                else:
                    self.current_span = 1
                self.mode = MODE_VARIABLE
                continue
            elif line.startswith( "fixedStep" ):
                header = parse_header( line )
                self.current_chrom = header['chrom']
                self.current_pos = int( header['start'] ) - 1
                self.current_step = int( header['step'] )
                if 'span' in header:
                    self.current_span = int( header['span'] )
                else:
                    self.current_span = 1
                self.mode = MODE_FIXED
                continue
            elif self.mode == MODE_BED:
                fields = line.split()
                if len( fields ) > 3:
                    if len( fields ) > 5:
                        # chrom, start, end, strand, value.
                        return Interval(ref=fields[0], start=int(fields[1]), 
                                        end=int(fields[2]), 
                                        strand=strand_str_to_int(fields[5]), 
                                        value=float(fields[3]))
                    else:
                        return Interval(ref=fields[0], start=int(fields[1]), 
                                        end=int(fields[2]), strand=POS_STRAND, 
                                        value=float(fields[3]))
            elif self.mode == MODE_VARIABLE: 
                fields = line.split()
                try:
                    pos = int( fields[0] ) - 1
                    val = float( fields[1] )
                except ValueError:
                    continue
                return Interval(self.current_chrom, pos, 
                                pos + self.current_span, strand=POS_STRAND, 
                                value=val)
            elif self.mode == MODE_FIXED:
                fields = line.split()
                try:
                    val = float( fields[0] )
                except ValueError:
                    continue
                return Interval(self.current_chrom, self.current_pos, 
                                self.current_pos + self.current_span, 
                                strand=POS_STRAND, value=val)
                self.current_pos += self.current_step
            else:
                raise "Unexpected input line: %s" % line.strip()
