'''
.. module:: test_intervaltrack
   :synopsis: unit tests
.. moduleauthor:: Matthew Iyer <mkiyer@umich.edu>
.. versionadded:: 0.0.1

Created on September 8, 2010
@author: mkiyer
'''
import unittest
import tempfile
import os
import collections
import numpy as np

from trackfactory.track import TrackError
from trackfactory.trackfactory import TrackFactory
from trackfactory.intervaltrack import IntervalTrack, get_base_dtype_fields

def mktemp(prefix, suffix):
    fh, filename = tempfile.mkstemp(suffix=suffix, prefix=prefix)
    os.close(fh)
    return filename

class Interval(object):
    __slots__ = ('ref', 'start', 'end', 'id')
    def __init__(self, ref, start, end, id):
        self.ref = ref
        self.start = start
        self.end = end
        self.id = id

class TestIntervalTrack(unittest.TestCase):

    def setUp(self):
        filename = mktemp(prefix="tmp", suffix=".h5")
        self.filename = filename
        self.refs = (('gene1', 1200), ('gene2', 200))
        self.tf = TrackFactory(self.filename, 'w', refs=self.refs) 
        self.interval_dtype = get_base_dtype_fields() + [('id', 'u4')]
        self.t = self.tf.create_track('intervals1', IntervalTrack,
                                      dtype=self.interval_dtype)
        
    def tearDown(self):
        self.tf.close()
        if os.path.exists(self.filename):
            os.remove(self.filename)

    def test_save_load_index(self):
        """create a table, add intervals, save the index, and 
        ensure the index persists when reloading the interval table
        """
        mytrack = self.t
        ref = 'gene1'
        intervals = []
        for i in xrange(100):
            start = np.random.randint(0, 1150)
            end = start + np.random.randint(1, 50)
            interval = Interval(ref, start, end, i)
            intervals.append(interval)
            mytrack.add(interval)
        # create and save index
        mytrack.index(persist=True)
        # save the index table
        tblcopy = mytrack.indexes['gene1'].tree.tbl.copy()
        tblroot = mytrack.indexes['gene1'].tree.root_id
        # close and reopen the file
        self.tf.close()
        self.tf = TrackFactory(self.filename, 'r', refs=self.refs)
        mytrack = self.tf.get_track('intervals1')
        # compare indexes
        self.assertTrue(np.all(tblcopy == mytrack.indexes['gene1'].tree.tbl))
        self.assertEqual(tblroot, mytrack.indexes['gene1'].tree.root_id)

    def testNonOverlapping(self):
        mytrack = self.t
        # insert some non-overlapping intervals
        intervals = []
        for i,start in enumerate(xrange(0, 900, 50)):
            ref = 'gene1'
            end = start + 10
            interval = Interval(ref, start, end, i)
            intervals.append(interval)
            mytrack.add(interval)        
        # check intervals
        for i,interval in enumerate(intervals):
            hits = mytrack.intersect(interval.ref, interval.start, interval.end)
            self.assertTrue([x['id'] for x in hits] == [i])
        # create interval tree index
        mytrack.index()
        for i,interval in enumerate(intervals):
            hits = mytrack.intersect(interval.ref, interval.start, interval.end)
            self.assertTrue([x['id'] for x in hits] == [i])
 
    def testOverlappingIntervals(self):
        mytrack = self.t
        # insert overlapping intervals
        intervals = []
        for i,start in enumerate(xrange(0, 900, 50)):
            ref = 'gene1'
            end = start + 60
            interval = Interval(ref, start, end, i)
            intervals.append(interval)
            mytrack.add(interval)
        # check intervals        
        self.assertTrue([x['id'] for x in mytrack.intersect(intervals[0].ref, intervals[0].start, intervals[0].end)] == [0,1])
        for i in xrange(1, len(intervals) - 1):
            interval = intervals[i]
            self.assertTrue([x['id'] for x in mytrack.intersect(interval.ref, interval.start, interval.end)] == [i-1, i, i+1])
        self.assertTrue([x['id'] for x in mytrack.intersect(intervals[-1].ref, intervals[-1].start, intervals[-1].end)] == [len(intervals) - 2, len(intervals)-1])
        # create interval tree index
        mytrack.index()
        self.assertTrue([x['id'] for x in mytrack.intersect(intervals[0].ref, intervals[0].start, intervals[0].end)] == [0,1])
        for i in xrange(1, len(intervals) - 1):
            interval = intervals[i]
            self.assertTrue([x['id'] for x in mytrack.intersect(interval.ref, interval.start, interval.end)] == [i-1, i, i+1])
        self.assertTrue([x['id'] for x in mytrack.intersect(intervals[-1].ref, intervals[-1].start, intervals[-1].end)] == [len(intervals) - 2, len(intervals)-1])

    def testBeforeAfter(self):
        mytrack = self.t
        # insert overlapping intervals
        intervals = []
        for i,start in enumerate(xrange(0, 1000, 50)):
            ref = 'gene1'
            end = start + 10
            interval = Interval(ref, start, end, i)            
            intervals.append(interval)
            mytrack.add(interval)
        # check intervals
        mytrack.index()
        for i,interval in enumerate(intervals):
            res = mytrack.before('gene1', interval.start, num_intervals=len(intervals), max_dist=1200)
            if [x['id'] for x in res] != range(i-1, -1, -1):
                print res, range(i-1, -1, -1)
            self.assertTrue([x['id'] for x in res] == range(i-1, -1, -1))
            res = mytrack.after('gene1', interval.end, num_intervals=len(intervals), max_dist=1200)
            if [x['id'] for x in res] != range(i+1, len(intervals)):
                print res, range(i-1, -1, -1)
            self.assertTrue([x['id'] for x in res] == range(i+1, len(intervals)))

    def _boundary_checks(self, mytrack):
        # test left boundary
        self.assertTrue(len(mytrack.before('gene1', 199, 1, 2000)) == 0) 
        self.assertTrue(len(mytrack.before('gene1', 200, 1, 2000)) == 1) 
        self.assertTrue(len(mytrack.before('gene1', 201, 1, 2000)) == 1) 
        self.assertTrue(len(mytrack.before('gene1', 202, 1, 2000)) == 1) 
        # test right boundary
        self.assertTrue(len(mytrack.after('gene1', 101, 1, 2000)) == 0) 
        self.assertTrue(len(mytrack.after('gene1', 100, 1, 2000)) == 0) 
        self.assertTrue(len(mytrack.after('gene1', 99, 1, 2000)) == 1) 
        self.assertTrue(len(mytrack.after('gene1', 98, 1, 2000)) == 1) 
        # test left max dist
        self.assertTrue(len(mytrack.before('gene1', 200, 1, 1)) == 1) 
        self.assertTrue(len(mytrack.before('gene1', 201, 1, 1)) == 0) 
        self.assertTrue(len(mytrack.before('gene1', 202, 1, 1)) == 0) 
        self.assertTrue(len(mytrack.before('gene1', 300, 1, 102)) == 1) 
        self.assertTrue(len(mytrack.before('gene1', 300, 1, 101)) == 1) 
        self.assertTrue(len(mytrack.before('gene1', 300, 1, 100)) == 0) 
        self.assertTrue(len(mytrack.before('gene1', 300, 1, 99)) == 0) 
        # test right max dist
        self.assertTrue(len(mytrack.after('gene1', 99, 1, 1)) == 1) 
        self.assertTrue(len(mytrack.after('gene1', 98, 1, 1)) == 0) 
        self.assertTrue(len(mytrack.after('gene1', 98, 1, 2)) == 1) 
        self.assertTrue(len(mytrack.after('gene1', 0, 1, 100)) == 1) 
        self.assertTrue(len(mytrack.after('gene1', 0, 1, 99)) == 0)

    def testBeforeAfterBoundaries(self):
        mytrack = self.t
        # insert overlapping intervals
        mytrack.add(Interval('gene1', 100, 200, 0))
        self._boundary_checks(mytrack)
        mytrack.index()
        self._boundary_checks(mytrack)
 
    def _intersect_checks(self, mytrack):
        self.assertTrue(len(mytrack.intersect('gene1', 0, 100)) == 0)
        self.assertTrue(len(mytrack.intersect('gene1', 0, 101)) == 1)
        self.assertTrue(len(mytrack.intersect('gene1', 200, 210)) == 0)
        self.assertTrue(len(mytrack.intersect('gene1', 199, 210)) == 1)

    def testIntersectBoundaries(self):
        mytrack = self.t
        # insert overlapping intervals
        mytrack.add(Interval('gene1', 100, 200, 0))
        self._intersect_checks(mytrack)
        mytrack.index()
        self._intersect_checks(mytrack)

    def testMultipleReferences(self):
        mytrack = self.t
        # add intervals to different references
        ref = 'gene1'        
        for i in xrange(0, 10):
            mytrack.add(Interval(ref, i, i+10, i))
            if ref == 'gene1':
                ref = 'gene2'
            else:
                ref = 'gene1'
        ids = [r['id'] for r in mytrack.intersect('gene1', 0, 20)]
        self.assertTrue(set(ids) == set(range(0, 10, 2)))
        ids = [r['id'] for r in mytrack.intersect('gene2', 0, 20)]
        self.assertTrue(set(ids) == set(range(1, 10, 2)))

    def testAdd(self):
        mytrack = self.t
        myrow = np.zeros(1, dtype=self.interval_dtype)[0]
        intervals = []
        for i,start in enumerate(xrange(0, 1000, 50)):
            myrow['ref'] = 'gene1'
            myrow['start'] = start
            myrow['end'] = start + 10
            myrow['id'] = i
            ref = 'gene1'
            end = start + 10
            interval = Interval(ref, start, end, i)
            mytrack.add(myrow)
            intervals.append(interval)
        for i,interval in enumerate(intervals):
            res = mytrack[(interval.ref, interval.start, interval.end)]
            self.assertTrue(len(res) == 1)
            hit = res[0]
            self.assertTrue((interval.ref, interval.start, interval.end, interval.id) ==
                            (hit['ref'], hit['start'], hit['end'], hit['id']))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
