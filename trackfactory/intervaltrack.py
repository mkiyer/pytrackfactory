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

from track import Track, TrackError, parse_interval
import lib.npintervaltree as intervaltree

REF_COL_NAME = "ref"
START_COL_NAME = "start"
END_COL_NAME = "end"
ID_COL_NAME = "id"
INTERVAL_DATA_TABLE = "interval_data"
INTERVAL_INDEX_GROUP = "interval_trees"
ROW_ATTR = 'row'
ROOT_ID_ATTR = 'root_id'

def get_base_dtype_fields():
    '''Get field names and numpy types for fields *required* to support the 
    :class:`IntervalTrack`.  For use when deriving subclasses.
    :returns: list of (name,type) tuples
    '''    
    return [(REF_COL_NAME, 'a32'),
            (START_COL_NAME, '<i4'),
            (END_COL_NAME, '<i4'),
            (ID_COL_NAME, '<u4')]

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
        def __init__(self, tree=None, root_id=0, dirty=True):
            self.tree = tree
            self.root_id = root_id
            self.dirty = dirty

    def __init__(self, hdf_group, dtype=None, length=None):        
        super(IntervalTrack, self).__init__(hdf_group)
        if not INTERVAL_DATA_TABLE in self.hdf_group:
            if dtype is None:
                dtype = np.dtype(get_base_dtype_fields())
            self._init_interval_table(dtype, length)
        self._init_trees()

    def __getitem__(self, key):
        ref, start, end = parse_interval(key)
        return self.intersect(ref, start, end)

    def __iter__(self):
        for row in self.hdf_group[INTERVAL_DATA_TABLE]:
            yield row

    def _get_num_intervals(self):
        return self.hdf_group._f_getChild(INTERVAL_DATA_TABLE).nrows
    num_intervals = property(_get_num_intervals, None, None, "get number of intervals in track")

    def _init_interval_table(self, dtype, length):
        if length is None:
            length = self.h5_expectedrows
        # create interval table
        h5file = self._get_hdf_file()
        h5file.createTable(self.hdf_group, 
                           name=INTERVAL_DATA_TABLE, 
                           description=dtype,
                           filters=self.h5_filters,
                           expectedrows=length)

    def _init_trees(self):
        h5file = self._get_hdf_file()
        if not INTERVAL_INDEX_GROUP in self.hdf_group:            
            h5file.createGroup(self.hdf_group, INTERVAL_INDEX_GROUP)
        self.tree_indexes = {}
        for rname in self.get_rnames():
            self.tree_indexes[rname] = self._load_index(rname)
    
    def _get_tree(self, rname):
        return self.hdf_group.interval_trees._f_getChild(rname)
    
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
        tbl = self.hdf_group._f_getChild(INTERVAL_DATA_TABLE)
        row = tbl.row
        # assign interval a unique primary key
        for colname in tbl.colnames:
            row[colname] = getfunc(colname)
        row.append()
        tbl.flush()
        # set index to dirty
        rname = getfunc(REF_COL_NAME)    
        self.tree_indexes[rname].dirty = True

    def _save_index(self, ref, tree, root_id):
        h5file = self._get_hdf_file()
        parentgroup = self.hdf_group._f_getChild(INTERVAL_INDEX_GROUP)
        if not ref in parentgroup:
            tbl = h5file.createTable(parentgroup, ref, tree)
            tbl.attrs[ROOT_ID_ATTR] = root_id

    def _load_index(self, ref):
        parentgroup = self.hdf_group._f_getChild(INTERVAL_INDEX_GROUP)
        if not ref in parentgroup:
            return self._Index()
        tbl = self._get_tree(ref)
        tree = tbl.read()
        root_id = tbl.attrs[ROOT_ID_ATTR]
        return self._Index(tree, root_id, False)

    def _get_ref_count_dict(self):
        refdict = collections.defaultdict(lambda: 0)
        for row in self.hdf_group._f_getChild(INTERVAL_DATA_TABLE):
            refdict[row[REF_COL_NAME]] += 1
        return refdict

    def _create_index(self, myref, nrows):                      
        # create index from scratch
        tbl = self.hdf_group._f_getChild(INTERVAL_DATA_TABLE)
        tree, root_id = intervaltree.init_tree(nrows)
        # iterate and insert all intervals
        cur_id = 1
        for row in tbl.where('ref == myref'):
            start = row[START_COL_NAME]
            end = row[END_COL_NAME]
            root_id = intervaltree.insert(tree, root_id, cur_id, start, end, row.nrow)
            cur_id += 1
        idx = self._Index(tree, root_id, False)
        return idx
    
    def index(self, persist=True):
        """(re)builds the interval tree index into the interval list you
        have stored.
        
        .. note: interval tree performance will not be obtained unless the user *explicitly* calls this method.  This is typically called after a set of intervals has been added.
        
        :param persist: permanently store the index on disk
        :type persist: bool
        """
        refcountdict = None
        for rname in self.get_rnames():
            # see if in-memory index exists
            idx = self.tree_indexes[rname]
            if (idx.tree is None) or idx.dirty:
                # count number of items for each reference
                if refcountdict is None:
                    refcountdict = self._get_ref_count_dict()
                # create index from scratch
                idx = self._create_index(rname, refcountdict[rname])
                # save for future use
                if persist:
                    self._save_index(rname, idx.tree, idx.root_id)
                self.tree_indexes[rname] = idx

    def intersect(self, rname, start, end):
        """Finds all intervals that have at least 1-bp of overlap
        with the interval provided
        """
        tbl = self.hdf_group._f_getChild(INTERVAL_DATA_TABLE)
        idx = self.tree_indexes[rname]
        if (idx.tree is None) or idx.dirty:
            hits = [r.fetch_all_fields() for r in tbl.where('((ref == myref) & (start < myend) & (end > mystart))', {'myref': rname, 'mystart': start, 'myend': end})]
            return hits
        else:
            row_ids = intervaltree.intersect(idx.tree, idx.root_id, start, end)
            return tbl[row_ids]

    def before(self, ref, position, num_intervals=1, max_dist=2500):
        """Find `num_intervals` intervals that lie before `position` and 
        are no further than `max_dist` distance away
        """
        tbl = self.hdf_group.interval_data
        idx = self.tree_indexes[ref]
        if (idx.tree is None) or idx.dirty:
            results = [(r['end'],r['id']) for r in tbl.where('((ref == myref) & (end <= position) & (position < end + max_dist))', {'myref': ref, 'position': position, 'max_dist': max_dist})]
            results.sort(key=operator.itemgetter(0), reverse=True)
            row_ids = [x[1] for x in results[:num_intervals]]
            return tbl[row_ids]           
        else:
            row_ids = intervaltree.left(idx.tree, idx.root_id, position, num_intervals, max_dist)
            return tbl[row_ids]

    def after(self, ref, position, num_intervals=1, max_dist=2500):
        """Find `num_intervals` intervals that lie after `position` and are 
        no further than `max_dist` distance away
        """
        tbl = self.hdf_group.interval_data
        idx = self.tree_indexes[ref]
        if (idx.tree is None) or idx.dirty:
            results = [(r['start'],r['id']) for r in tbl.where('((ref == myref) & (position < start) & (start <= (position + max_dist)))', {'myref': ref, 'position': position, 'max_dist': max_dist})]
            results.sort(key=operator.itemgetter(0))
            row_ids = [x[1] for x in results[:num_intervals]]
            return tbl[row_ids]
        else:
            row_ids = intervaltree.right(idx.tree, idx.root_id, position, num_intervals, max_dist)
            return tbl[row_ids]
