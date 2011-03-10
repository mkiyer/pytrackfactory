'''
Created on Mar 9, 2011

@author: mkiyer
'''
import unittest
import tempfile
import os
import numpy as np

from ..track import TrackError
from ..trackfactory import TrackFactory
from ..arraytrack import ArrayTrack

def mktemp(prefix, suffix):
    fh, filename = tempfile.mkstemp(suffix=suffix, prefix=prefix)
    os.close(fh)
    return filename

def random_intervals(n, length, isize_max, dtype):
    intervals = []
    correct = np.zeros(100000, dtype)
    for i in xrange(1000):
        start = np.random.randint(0, 100000-isize_max)
        end = start + np.random.randint(0, isize_max)
        intervals.append(('gene1', start, end, ".", i))
        correct[start:end] += i
    return intervals, correct
    
class TestReferences(unittest.TestCase):
    def setUp(self):
        filename = mktemp(prefix="tmp", suffix=".h5")
        self.filename = filename
        self.refs = (('gene1', 100000), ('gene2', 10))
        self.tf = TrackFactory(self.filename, 'w', refs=self.refs) 

    def tearDown(self):
        self.tf.close()
        if os.path.exists(self.filename):
            os.remove(self.filename)

    def test_fromintervals(self):
        dtype = "i4"
        intervals1, correct1 = random_intervals(1000, 100000, 500, dtype)
        # try one channel array
        t = self.tf.create_track("a", ArrayTrack, dtype, channels=1)
        t.fromintervals(iter(intervals1), channel=0)
        self.assertTrue(np.all(t["gene1"][:,0] == correct1))
        # try three channel array
        intervals2, correct2 = random_intervals(1000, 100000, 500, dtype)
        intervals3, correct3 = random_intervals(1000, 100000, 500, dtype)
        t = self.tf.create_track("b", ArrayTrack, dtype, channels=3)
        t.fromintervals(iter(intervals1), channel=0)
        t.fromintervals(iter(intervals2), channel=1)
        t.fromintervals(iter(intervals3), channel=2)
        self.assertTrue(np.all(t["gene1"][:,0] == correct1))
        self.assertFalse(np.all(t["gene1"][:,0] == correct2))
        self.assertTrue(np.all(t["gene1"][:,1] == correct2))
        self.assertFalse(np.all(t["gene1"][:,1] == correct3))
        self.assertTrue(np.all(t["gene1"][:,2] == correct3))
        self.assertFalse(np.all(t["gene1"][:,2] == correct1))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()