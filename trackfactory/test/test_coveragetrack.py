'''
Created on Mar 9, 2011

@author: mkiyer
'''
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
from ..coveragetrack import StrandedCoverageTrack

def mktemp(prefix, suffix):
    fh, filename = tempfile.mkstemp(suffix=suffix, prefix=prefix)
    os.close(fh)
    return filename

def random_intervals(n, length, isize_max, dtype):
    intervals = []
    correct = np.zeros(length, dtype)
    for i in xrange(1000):
        start = np.random.randint(0, length-isize_max)
        end = start + np.random.randint(1, isize_max)
        intervals.append(('gene1', start, end, 0, i))
        correct[start:end] += i
    return intervals, correct

def random_stranded_intervals(n, length, isize_max, dtype):
    intervals = []
    correct = np.zeros((length,2), dtype)
    for i in xrange(n):
        start = np.random.randint(0, length-isize_max)
        end = start + np.random.randint(1, isize_max)
        intervals.append(('gene1', start, end, 0, i+1))
        intervals.append(('gene1', start, end, 1, -i))
        correct[start:end] += [i+1,-i]
    return intervals, correct
    
class TestCoverageTrack(unittest.TestCase):
    def setUp(self):
        filename = mktemp(prefix="tmp", suffix=".h5")
        self.filename = filename
        self.length = 100000
        self.isize_max = 500
        self.refs = (('gene1', self.length), ('gene2', 10))
        self.tf = TrackFactory(self.filename, 'w', refs=self.refs) 

    def tearDown(self):
        self.tf.close()
        if os.path.exists(self.filename):
            os.remove(self.filename)

    def test_stranded_intervals(self):
        dtype = "i4"
        intervals1, correct1 = random_stranded_intervals(1000, self.length, self.isize_max, dtype)
        total_cov = correct1.sum()
        t = self.tf.create_track("a", StrandedCoverageTrack)
        # test loading from intervals
        t.fromintervals(iter(intervals1))
        self.assertTrue(np.all(t["gene1"] == correct1))
        # test count function
        intervals2, correct2 = random_intervals(1000, self.length, self.isize_max, dtype)
        for interval in intervals2:
            ref, start, end, strand, val = interval
            # check plus strand
            mycount = t.count((ref, start, end, 0, val))
            mycov = t.coverage((ref, start, end, 0, val), multiplier=1.0)
            mydens = t.density((ref, start, end, 0, val), multiplier=1.0)
            correctcount = correct1[start:end,0].sum()
            correctcov = correct1[start:end,0] / float(total_cov)
            correctdens = correctcount / float(total_cov * (end - start))
            self.assertAlmostEqual(mycount, correctcount)
            self.assertTrue(np.allclose(mycov, correctcov, atol=1e-4))
            self.assertTrue(np.allclose(mydens, correctdens, atol=1e-4))            
            # check minus strand
            mycount = t.count((ref, start, end, 1, val))
            mycov = t.coverage((ref, start, end, 1, val), multiplier=1.0)
            mydens = t.density((ref, start, end, 1, val), multiplier=1.0)
            correctcount = correct1[start:end,1].sum()
            correctcov = correct1[start:end,1] / float(total_cov)
            correctdens = correctcount / float(total_cov * (end - start))
            self.assertAlmostEqual(mycount, correctcount)
            self.assertTrue(np.allclose(mycov, correctcov, atol=1e-4))
            self.assertTrue(np.allclose(mydens, correctdens, atol=1e-4))            
            # check both strands
            mycount = t.count((ref, start, end, 2, val))
            mycov = t.coverage((ref, start, end, 2, val), multiplier=1.0)
            mydens = t.density((ref, start, end, 2, val), multiplier=1.0)
            correctcount = correct1[start:end].sum()
            correctcov = correct1[start:end].sum(axis=1) / float(total_cov)
            correctdens = correctcount / float(total_cov * (end - start))
            self.assertAlmostEqual(mycount, correctcount)
            self.assertTrue(np.allclose(mycov, correctcov, atol=1e-4))
            self.assertTrue(np.allclose(mydens, correctdens, atol=1e-4))            



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()