"""
.. module:: track
   :synopsis: base classes and definitions
.. moduleauthor:: Matthew Iyer <mkiyer@umich.edu>
.. versionadded:: 0.0.1
"""
import subprocess
import pysam

# constants
REF_TRACK_NAME = 'references'
REF_NODE_PATH = "/" + REF_TRACK_NAME
TRACK_CLASS_ATTR = 'track_class'

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

    def _get_hdf_file(self):
        return self.hdf_group._v_file
    def _get_references_node(self):
        return self.hdf_group._v_file.getNode(REF_NODE_PATH)
    def _parse_interval(self, key):
        ref, start, end = parse_interval(key)
        if (ref not in self.hdf_group):
            raise TrackError("unknown reference '%s'" % (ref))
        if end is None:
            end = self.get_ref_length(ref)        
        ca = self.hdf_group._f_getChild(ref)
        return ca, start, end

    def get_type(self):
        return self.hdf_group._v_attrs[TRACK_CLASS_ATTR]

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

def parse_interval(interval):
    """parses a genomic interval specifier in either the string
    format `<ref:start-end>` or tuple (ref, start, end)
    
    :param interval: interval specifier
    :type interval: str or tuple
    
    >>> parse_interval("chr1:100-200")
    ("chr1", 100, 200)
    >>> parse_interval("chr1:1,000-20,000")
    ("chr1", 1000, 20000)
    >>> parse_interval("chrX")
    ("chrX", None, None)
    >>> parse_interval("chrX:5")
    ("chrX", 5, 6)
    >>> parse_interval(("chr2", 100, 300))
    ("chr2", 100, 300) 
    """
    if interval is None:
        return None, None, None
    if isinstance(interval, basestring):
        # remove commas
        k = interval.replace(',','')            
        # split chrom on ':' character
        if k.find(':') != -1:
            ref, startend = k.split(':')                
            # separate start/end on '-' character
            if k.find('-') != -1:
                start, end = map(int, startend.split('-'))
            else:
                start = int(startend)
                end = start + 1
        else:
            ref = k
            start = None
            end = None
        return ref, start, end
    else:
        if len(interval) == 1:
            return interval[0], None, None
        elif len(interval) == 2:
            return interval[0], interval[1], interval[1]+1
        elif len(interval) >= 3:
            return interval[:3]

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
