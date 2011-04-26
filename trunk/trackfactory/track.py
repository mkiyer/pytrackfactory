"""
.. module:: track
   :synopsis: base classes and definitions
.. moduleauthor:: Matthew Iyer <mkiyer@umich.edu>
.. versionadded:: 0.0.1
"""
import subprocess
import pysam
import re

# constants
POS_STRAND = 0
NEG_STRAND = 1
NO_STRAND = 2
_strand_str_to_int = {"+": POS_STRAND,
                      "-": NEG_STRAND,
                      ".": NO_STRAND}
_strand_int_to_str = dict([(v,k) for k,v in _strand_str_to_int.items()])
def strand_str_to_int(strand):
    return _strand_str_to_int.get(strand, NO_STRAND)
def strand_int_to_str(strand):
    return _strand_int_to_str.get(strand, NO_STRAND)

REF_TRACK_NAME = 'references'
REF_NODE_PATH = "/" + REF_TRACK_NAME
TRACK_CLASS_ATTR = 'track_class'

_interval_string_re = re.compile(r'^(?P<ref>\w+)(?:\[(?P<strand>[+-.])\])?(?::(?P<start>[\d,]+)(?:-(?P<end>[\d,]+))?)?')

def get_ref_length(tbl, ref):
    length = None
    for r in tbl.where('name == ref'):
        length = r['length']
        break
    return length

class TrackError(Exception):
    """
    Exception class
    """
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return repr(self.msg)

class Track(object):
    """
    Base class for all Tracks.
    
    Contains methods to access the references and lengths that the 
    TrackFactory was created with.
    """
    def __init__(self, hdf_group):
        self.hdf_group = hdf_group

    def get_name(self):
        return self.hdf_group._v_name
    def get_type(self):
        return self.hdf_group._v_attrs[TRACK_CLASS_ATTR]
    def _get_hdf_file(self):
        return self.hdf_group._v_file
    def _get_references_node(self):
        return self.hdf_group._v_file.getNode(REF_NODE_PATH)
    def _parse_interval(self, key):
        ref, start, end, strand = parse_interval(key)
        if ref is None:
            return None, None, None, strand
        if not self.has_ref(ref): 
            raise TrackError("unknown reference '%s'" % (ref))
        if start is None:
            start = 0
        if end is None:
            end = self.get_ref_length(ref)        
        return ref, start, end, strand

    def has_ref(self, ref):
        res = False
        tbl = self._get_references_node()
        for r in tbl.where('name == ref'):
            return True
        return False

    def get_ref_length(self, ref):
        """return the size of the reference `ref`"""
        tbl = self._get_references_node()
        return get_ref_length(tbl, ref)

    def get_rnames(self):
        """generator function that returns list of references names from the 
        :mod:`TrackFactory` used to create this track
        
        :returns: iterator of references names (ordered by reference id)
        """
        for row in self._get_references_node():
            yield row['name']
    
    def get_lengths(self):
        """generator function that returns list of reference lengths from the 
        :mod:`TrackFactory` used to create this track
        
        :returns: iterator lengths (ordered by reference id)
        """
        for row in self._get_references_node():
            yield row['name']

    def get_refs(self):
        """generator function that returns list of (reference,length) tuples
        from the :mod:`TrackFactory` used to create this track
        
        :returns: iterator of (reference,length) tuples (ordered by ref id)
        """        
        for row in self._get_references_node():
            yield row.fetch_all_fields()


def _parse_interval_string(interval_string):
    """ref[strand]:start-end 
    
    `strand` can be '+', '-', or '.'    
    """
    m = _interval_string_re.match(interval_string)
    if m is None:
        ref, start, end, strand = None, None, None, NO_STRAND
        # if interval string is just a single character,
        # try to extract something useful
        if len(interval_string) == 1:
            strand = strand_str_to_int(interval_string[0])
    else:
        ref = m.group('ref')
        start = m.group('start')
        end = m.group('end')
        strand = m.group('strand')
        if start is not None: 
            start = int(start.replace(",", ""))
        if end is not None: 
            end = int(end.replace(",", ""))
        if strand is None: 
            strand = NO_STRAND
        else:
            strand = strand_str_to_int(strand)
    return ref, start, end, strand

def _parse_interval_tuple(interval):
    if interval is None:
        return None, None, None, NO_STRAND
    if len(interval) == 1:
        return interval[0], None, None, NO_STRAND
    elif len(interval) == 2:
        return interval[0], interval[1], interval[1]+1, NO_STRAND
    elif len(interval) == 3:
        return interval[0], interval[1], interval[2], NO_STRAND
    else:
        if isinstance(interval[3], basestring):
            strand = strand_str_to_int(interval[3])
        else:
            strand = interval[3]
        assert (strand >= POS_STRAND) and (strand <= NO_STRAND)
        return interval[0], interval[1], interval[2], strand

def parse_interval(interval):
    """parses a genomic interval specifier in either the string
    format `<ref[strand]:start-end>` or tuple (ref, start, end, strand)
    
    :param interval: interval specifier
    :type interval: str or tuple
    
    >>> parse_interval("chr1:100-200")
    ("chr1", 100, 200, ".")
    >>> parse_interval("chr1:1,000-20,000")
    ("chr1", 1000, 20000)
    >>> parse_interval("chr1[-]:1,000-20,000")
    ("chr1", 1000, 20000, "-")
    >>> parse_interval("chrX")
    ("chrX", None, None, ".")
    >>> parse_interval("chrX:5")
    ("chrX", 5, 6, ".")
    >>> parse_interval(("chr2", 100, 300))
    ("chr2", 100, 300, ".")
    """
    if isinstance(interval, basestring):
        interval = _parse_interval_string(interval)
    interval = _parse_interval_tuple(interval)
    return interval

def get_refs_from_bowtie_index(bowtie_index, 
                               bowtie_inspect_bin='bowtie-inspect'):
    args = [bowtie_inspect_bin, '-s', bowtie_index]    
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    output = p.communicate()[0]
    refs = []
    for line in output.split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        if fields[0].startswith('Sequence'):
            refs.append((fields[1], int(fields[2])))
    return refs

def get_refs_from_bam(bamfile):
    f = pysam.Samfile(bamfile, "rb")
    refs = list(zip(f.references, f.lengths))
    f.close()
    return refs

def get_refs_from_sam(samfile):
    f = pysam.Samfile(samfile, "r")
    refs = list(zip(f.references, f.lengths))
    f.close()
    return refs
