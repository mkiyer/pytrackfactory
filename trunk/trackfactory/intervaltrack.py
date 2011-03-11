'''
.. module:: intervaltrack
   :synopsis: track for storing genomic intervals
.. moduleauthor:: Matthew Iyer <mkiyer@umich.edu>
.. versionadded:: 0.0.1

Created on Sep 14, 2010

The IntervalTrack was derived from the IntervalTree implementation from the
bx-python project by James Taylor, Ian Schenk, Brent Pedersen, and others.
Here, it has been adapted to the TrackFactory interfaces via HDF5 and numpy.

:Authors: James Taylor (james@jamestaylor.org),
          Ian Schenk (ian.schenck@gmail.com),
          Brent Pedersen (bpederse@gmail.com),
          Matthew Iyer (matthew.iyer@gmail.com)
'''
# Historical note:
#    This module original contained an implementation based on sorted endpoints
#    and a binary search, using an idea from Scott Schwartz and Piotr Berman.
#    Later an interval tree implementation was implemented by Ian for Galaxy's
#    join tool (see `bx.intervals.operations.quicksect.py`). This was then
#    converted to Cython by Brent, who also added support for
#    upstream/downstream/neighbor queries. This was modified by James to
#    handle half-open intervals strictly, to maintain sort order, and to
#    implement the same interface as the original Intersecter.
#
#   Matthew Iyer then adapted the algorithm to create a persistent version of the
#   IntervalTree using HDF5 and numpy
import collections
import operator
import numpy as np
import tables
import logging

from track import Track, TrackError, parse_interval
from lib.intervaltree import IntervalTree

REF_COL_NAME = "ref"
START_COL_NAME = "start"
END_COL_NAME = "end"

INTERVAL_TABLE = "interval_data"
INTERVAL_INDEX_GROUP = "interval_trees"
ROW_ATTR = 'row'
ROOT_ID_ATTR = 'root_id'
DTYPE_ATTR = 'dtype'

def get_base_dtype_fields():
    '''Get field names and numpy types for fields *required* to support the 
    :class:`IntervalTrack`.  For use when deriving subclasses.
    :returns: list of (name,type) tuples
    '''    
    return [(REF_COL_NAME, 'a32'),
            (START_COL_NAME, '<i4'),
            (END_COL_NAME, '<i4')]
    
def get_bed_dtype_fields():
    fields = get_base_dtype_fields()
    fields.extend([("name", "a31"),
                   ("strand", 'a1'),
                   ("score", np.int)])
    return fields

class IntervalTrack(Track):
    '''Associates genomic intervals with custom metadata. This is a base 
    class that implements the indexing and lookup functionality but does
    little else.  Users will want to use one of the derived classes 
    (or derive their own) to achieve the benefit of this track.
    
    For performance reasons, it is best if you have some idea of the number 
    of intervals you expect to add.  The default is 10,000, but this is 
    arbitrary.  You can pass this to the 'length' parameter when creating
    the track.
    
    The `dtype` parameter should contain a :class:numpy.dtype describing 
    the metadata to be stored at each interval.    
    '''
    # pytables configuration
    h5_complevel = 0
    h5_complib = None
    #h5_complevel = 1
    #h5_complib = 'lzo'
    h5_filters = tables.Filters(complevel=h5_complevel, 
                                complib=h5_complib,
                                shuffle=True,
                                fletcher32=False)    
    h5_expectedrows = 1e5
    h5_chunksize = 1<<18 # 256kb chunks

    class _Index(object):
        def __init__(self, tree=None, dirty=False):
            self.tree = tree
            self.dirty = dirty

    def __init__(self, hdf_group, dtype=None, expectedrows=None):        
        super(IntervalTrack, self).__init__(hdf_group)
        if not INTERVAL_TABLE in self.hdf_group:
            if dtype is None:
                dtype = get_base_dtype_fields()
            dtype = np.dtype(dtype)
            self._init_table(dtype, expectedrows)
            hdf_group._v_attrs[DTYPE_ATTR] = dtype            
        self._init_index()

    def _get_dtype(self):
        return self.hdf_group._v_attrs[DTYPE_ATTR]
    def _get_interval_table(self):
        return self.hdf_group._f_getChild(INTERVAL_TABLE)

    def _init_table(self, dtype, expectedrows):
        if expectedrows is None:
            expectedrows = self.h5_expectedrows
        # create interval table
        h5file = self._get_hdf_file()
        h5file.createTable(self.hdf_group, 
                           name=INTERVAL_TABLE, 
                           description=dtype,
                           filters=self.h5_filters,
                           expectedrows=expectedrows)

    def _init_index(self):
        if not INTERVAL_INDEX_GROUP in self.hdf_group:
            h5file = self._get_hdf_file()
            h5file.createGroup(self.hdf_group, INTERVAL_INDEX_GROUP)
        parentgroup = self.hdf_group._f_getChild(INTERVAL_INDEX_GROUP)
        self.indexes = {}
        for rname in self.get_rnames():
            idx = self._Index()
            # load index if it exists
            if rname in parentgroup:
                tree = IntervalTree.fromhdf(parentgroup._f_getChild(rname))
                idx.tree = tree
            self.indexes[rname] = idx

    def _get_ref_count_dict(self):
        refdict = collections.defaultdict(lambda: 0)
        for row in self.hdf_group._f_getChild(INTERVAL_TABLE):
            refdict[row[REF_COL_NAME]] += 1
        return refdict
    
    def index(self, persist=True):
        """(re)builds the interval tree index into the interval list you
        have stored.
        
        .. note: interval tree performance will not be obtained unless 
        the user *explicitly* calls this method.  This is typically 
        called after a set of intervals has been added.
        
        :param persist: permanently store the index on disk
        :type persist: bool
        """
        refcountdict = None
        tbl = self.hdf_group._f_getChild(INTERVAL_TABLE)
        tree_group = self.hdf_group._f_getChild(INTERVAL_INDEX_GROUP)
        for rname, idx in self.indexes.iteritems():
            # see if index needs to be generated
            if (idx.tree is None) or idx.dirty:
                # count number of items for each reference
                if refcountdict is None:
                    refcountdict = self._get_ref_count_dict()
                # iterate and insert all intervals
                tree = IntervalTree(refcountdict[rname])            
                cur_id = 1
                for row in tbl.where('ref == rname'):
                    start = row[START_COL_NAME]
                    end = row[END_COL_NAME]
                    tree.insert(cur_id, start, end, row.nrow)
                    cur_id += 1
                # update index
                idx.tree = tree
                idx.dirty = False
            # save for future use
            if persist:
                tree.tohdf(tree_group, rname)
    
    def __getitem__(self, key):
        ref, start, end, strand = parse_interval(key)
        if start is None:
            start = 0
        if end is None:
            end = self.get_ref_length(ref)
        return self.intersect(ref, start, end)

    def __iter__(self):
        for row in self.hdf_group[INTERVAL_TABLE]:
            yield row

    def _get_num_intervals(self):
        return self.hdf_group._f_getChild(INTERVAL_TABLE).nrows
    num_intervals = property(_get_num_intervals, None, None, 
                             "get number of intervals in track")

    def intersect(self, rname, start, end):
        """Finds all intervals that have at least 1-bp of overlap
        with the interval provided
        """
        print rname, start, end
        tbl = self.hdf_group._f_getChild(INTERVAL_TABLE)
        idx = self.indexes[rname]
        if (idx.tree is None) or idx.dirty:
            expr = '((ref == myref) & (start < myend) & (end > mystart))'
            mappings = {'myref': rname, 'mystart': start, 'myend': end}
            hits = [r.fetch_all_fields() for r in tbl.where(expr, mappings)]
        else:
            row_ids = idx.tree.intersect(start, end)
            hits = tbl[row_ids]
        return hits

    def before(self, ref, position, num_intervals=1, max_dist=2500):
        """Find `num_intervals` intervals that lie before `position` and 
        are no further than `max_dist` distance away
        """
        tbl = self.hdf_group.interval_data
        idx = self.indexes[ref]
        if (idx.tree is None) or idx.dirty:
            expr = ('((ref == myref) & (end <= position) & '
                    '(position < end + max_dist))')
            mappings = {'myref': ref, 'position': position, 'max_dist': max_dist}
            results = [(r['end'],r['id']) for r in tbl.where(expr, mappings)]
            results.sort(key=operator.itemgetter(0), reverse=True)
            row_ids = [x[1] for x in results[:num_intervals]]
            hits = tbl[row_ids]           
        else:
            row_ids = idx.tree.left(position, num_intervals, max_dist)
            hits = tbl[row_ids]
        return hits

    def after(self, ref, position, num_intervals=1, max_dist=2500):
        """Find `num_intervals` intervals that lie after `position` and are 
        no further than `max_dist` distance away
        """
        tbl = self.hdf_group.interval_data
        idx = self.indexes[ref]
        if (idx.tree is None) or idx.dirty:
            expr = ('((ref == myref) & (position < start) & '
                    '(start <= (position + max_dist)))')
            mappings = {'myref': ref, 'position': position, 'max_dist': max_dist}
            results = [(r['start'],r['id']) for r in tbl.where(expr, mappings)]
            results.sort(key=operator.itemgetter(0))
            row_ids = [x[1] for x in results[:num_intervals]]
            hits = tbl[row_ids]
        else:
            row_ids = idx.tree.right(position, num_intervals, max_dist)
            hits = tbl[row_ids]
        return hits

    def add(self, value):
        """add an interval to the track
        
        :param value: interval and metadata to add
        :type value: `numpy.void`, dictionary, or object types accepted
        
        Example::
          # create a dummy track
          import pytrackfactory as ptf
          tf = ptf.TrackFactory("foo.tf", "w", refs=(("chr1", 10000)))          
          t = tf.create_track("a", IntervalTrack, length=10)
          
          # add via dictionary
          t.add({'ref': 'chr1', 'start': 0, 'end': 1000, 'id': 100})
          
          # add via numpy void type
          import numpy as np        
          rec = np.zeros(1, dtype=get_base_dtype_fields())[0]        
          rec['ref'] = 'chr1'        
          rec['start'] = 0        
          rec['end'] = 1000        
          rec['id'] = 100        
          t.add(rec)

          # add via user-defined class        
          class MyInterval(object):        
              def __init__(self, ref, start, end, id):        
                  self.ref = None
                  self.start = 0
                  self.end = 0
                  self.id = 0        
          x = MyInterval('chr1', 0, 1000, 100)
          t.add(x)
        
        Note that in all three cases, the attributes or fields *ref*, 
        *start*, *end*, and *id* were required.
        """
        if (isinstance(value, np.void) or 
            isinstance(value, dict)):
            getfunc = value.__getitem__            
        else:
            getfunc = lambda a: getattr(value, a)
        # add interval to intervals table
        tbl = self.hdf_group._f_getChild(INTERVAL_TABLE)
        row = tbl.row
        for colname in tbl.colnames:
            row[colname] = getfunc(colname)
        row.append()
        tbl.flush()
        # set index to dirty
        rname = getfunc(REF_COL_NAME)
        self.indexes[rname].dirty = True  

    def fromintervals(self, interval_iter, index=True):
        tbl = self.hdf_group._f_getChild(INTERVAL_TABLE)
        fields = self._get_dtype().names
        row = tbl.row
        flush_every = 1000
        rows = 0
        logging.info("[%s] inserting intervals" % (self.__class__.__name__))
        for interval in interval_iter:
            for i,colname in enumerate(fields):
                row[colname] = interval[i]
            row.append()
            rows += 1
            if rows % (flush_every) == 0:
                tbl.flush()
        tbl.flush()
        if index:
            logging.info("[%s] building index" % (self.__class__.__name__))
            self.index(persist=True)
