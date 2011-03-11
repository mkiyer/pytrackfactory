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
import random

from trackfactory.track import TrackError, POS_STRAND, NEG_STRAND, NO_STRAND
from trackfactory import TrackFactory
from trackfactory.vectortrack import VectorTrack, StrandedVectorTrack, StrandedAlleleVectorTrack
from trackfactory.io.cinterval import Interval, SequenceInterval

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
        intervals.append(Interval('gene1', start, end, POS_STRAND, i))
        correct[start:end] += i
    return intervals, correct

def random_stranded_intervals(n, length, isize_max, dtype):
    intervals = []
    correct = np.zeros((length,2), dtype)
    for i in xrange(n):
        start = np.random.randint(0, length-isize_max)
        end = start + np.random.randint(1, isize_max)
        intervals.append(Interval('gene1', start, end, POS_STRAND, i+1))
        intervals.append(Interval('gene1', start, end, NEG_STRAND, -i))
        correct[start:end] += [i+1,-i]
    return intervals, correct

def random_stranded_allele_intervals(n, length, isize_max, dtype):
    base_channel_dict = {"A": 0,
                         "G": 1,
                         "C": 2,
                         "T": 3}
    intervals = []
    correct = np.zeros((length,8), dtype)
    for i in xrange(n):
        start = np.random.randint(0, length-isize_max)
        end = start + np.random.randint(1, isize_max)        
        seq = []
        for x in xrange(start, end):            
            base = random.choice("ATGCN")
            channel = base_channel_dict.get(base, None)
            if channel is None:
                correct[x,0:4] += 0.5
                correct[x,4:8] -= 0.25
            else:
                correct[x,channel] += 2
                correct[x,channel+4] -= 1
            seq.append(base)
        seq = ''.join(seq)        
        intervals.append(SequenceInterval('gene1', start, end, POS_STRAND, 2, 
                                          seq=seq))
        intervals.append(SequenceInterval('gene1', start, end, NEG_STRAND, -1, 
                                          seq=seq))
    return intervals, correct

class TestVectorTrack(unittest.TestCase):
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

    def test_fromintervals(self):
        dtype = "i"
        intervals1, correct1 = random_intervals(1000, self.length, 
                                                self.isize_max, dtype)
        # try one channel array
        t = self.tf.create_track("a", VectorTrack, dtype, channels=1)
        t.fromintervals(iter(intervals1), channel=0)
        self.assertTrue(np.all(t["gene1"][:,0] == correct1))
        # try three channel array
        intervals2, correct2 = random_intervals(1000, self.length, 
                                                self.isize_max, dtype)
        intervals3, correct3 = random_intervals(1000, self.length, 
                                                self.isize_max, dtype)
        t = self.tf.create_track("b", VectorTrack, dtype, channels=3)
        t.fromintervals(iter(intervals1), channel=0)
        t.fromintervals(iter(intervals2), channel=1)
        t.fromintervals(iter(intervals3), channel=2)
        self.assertTrue(np.all(t["gene1"][:,0] == correct1))
        self.assertFalse(np.all(t["gene1"][:,0] == correct2))
        self.assertTrue(np.all(t["gene1"][:,1] == correct2))
        self.assertFalse(np.all(t["gene1"][:,1] == correct3))
        self.assertTrue(np.all(t["gene1"][:,2] == correct3))
        self.assertFalse(np.all(t["gene1"][:,2] == correct1))

    def test_stranded_intervals(self):
        """testing allocating coverage to both strands"""
        dtype = "i4"
        intervals1, correct1 = random_stranded_intervals(100, self.length, self.isize_max, dtype)
        total_cov = correct1.sum()
        t = self.tf.create_track("a", StrandedVectorTrack)
        # test loading from intervals
        t.fromintervals(iter(intervals1))
        self.assertTrue(np.all(t["gene1"] == correct1))
        # test count function
        intervals2, correct2 = random_intervals(1000, self.length, self.isize_max, dtype)
        for ival in intervals2:
            ref = ival.ref
            start = ival.start
            end = ival.end
            strand = ival.strand
            val = ival.value 
            # check plus strand
            # count
            mycount = t.count((ref, start, end, POS_STRAND, val))
            correctcount = correct1[start:end,0].sum()
            self.assertAlmostEqual(mycount, correctcount)
            # coverage
            mycov = t.coverage((ref, start, end, POS_STRAND, val), multiplier=1.0)
            correctcov = correct1[start:end,0] / float(total_cov)
            self.assertTrue(np.allclose(mycov, correctcov, atol=1e-4))
            # density
            mydens = t.density((ref, start, end, POS_STRAND, val), multiplier=1.0)
            correctdens = correctcount / float(total_cov * (end - start))            
            self.assertTrue(np.allclose(mydens, correctdens, atol=1e-4))            
            # check minus strand
            # count
            mycount = t.count((ref, start, end, NEG_STRAND, val))
            correctcount = correct1[start:end,1].sum()
            self.assertAlmostEqual(mycount, correctcount)
            # coverage
            mycov = t.coverage((ref, start, end, NEG_STRAND, val), multiplier=1.0)
            correctcov = correct1[start:end,1] / float(total_cov)
            self.assertTrue(np.allclose(mycov, correctcov, atol=1e-4))
            # density
            mydens = t.density((ref, start, end, NEG_STRAND, val), multiplier=1.0)
            correctdens = correctcount / float(total_cov * (end - start))
            self.assertTrue(np.allclose(mydens, correctdens, atol=1e-4))            
            # check both strands
            mycount = t.count((ref, start, end, NO_STRAND, val))
            correctcount = correct1[start:end].sum()
            self.assertAlmostEqual(mycount, correctcount)
            # cov
            mycov = t.coverage((ref, start, end, NO_STRAND, val), multiplier=1.0)
            correctcov = correct1[start:end].sum(axis=1) / float(total_cov)
            self.assertTrue(np.allclose(mycov, correctcov, atol=1e-4))
            # density
            mydens = t.density((ref, start, end, NO_STRAND, val), multiplier=1.0)
            correctdens = correctcount / float(total_cov * (end - start))
            self.assertTrue(np.allclose(mydens, correctdens, atol=1e-4))            

    def test_stranded_allele_intervals(self):
        """testing coverage with allele frequencies"""
        dtype = "f"
        intervals1, correct1 = \
            random_stranded_allele_intervals(100, self.length, 
                                             self.isize_max, dtype)
        total_cov = correct1.sum()
        t = self.tf.create_track("a", StrandedAlleleVectorTrack)
        # test loading from intervals
        t.fromintervals(iter(intervals1))
        self.assertTrue(np.all(t["gene1"] == correct1))        
        # test count function
        intervals2, correct2 = random_intervals(10, self.length, self.isize_max, dtype)
        for ival in intervals2:
            ref = ival.ref
            start = ival.start
            end = ival.end
            strand = ival.strand
            val = ival.value
            # check plus strand
            mycount = t.count((ref, start, end, POS_STRAND, val))
            correctcount = correct1[start:end,0:4].sum()
            self.assertAlmostEqual(mycount, correctcount)
            mycov = t.coverage((ref, start, end, POS_STRAND, val), multiplier=1.0)
            correctcov = correct1[start:end,0:4].sum(axis=1) / float(total_cov)
            self.assertTrue(np.allclose(mycov, correctcov, atol=1e-4))
            mydens = t.density((ref, start, end, POS_STRAND, val), multiplier=1.0)
            correctdens = correctcount / float(total_cov * (end - start))
            self.assertTrue(np.allclose(mydens, correctdens, atol=1e-4))
            # check minus strand
            mycount = t.count((ref, start, end, NEG_STRAND, val))
            correctcount = correct1[start:end,4:8].sum()
            self.assertAlmostEqual(mycount, correctcount)
            mycov = t.coverage((ref, start, end, NEG_STRAND, val), multiplier=1.0)
            correctcov = correct1[start:end,4:8].sum(axis=1) / float(total_cov)
            self.assertTrue(np.allclose(mycov, correctcov, atol=1e-4))
            mydens = t.density((ref, start, end, NEG_STRAND, val), multiplier=1.0)
            correctdens = correctcount / float(total_cov * (end - start))
            self.assertTrue(np.allclose(mydens, correctdens, atol=1e-4))            
            # check both strands
            mycount = t.count((ref, start, end, NO_STRAND, val))
            correctcount = correct1[start:end].sum()
            self.assertAlmostEqual(mycount, correctcount)
            mycov = t.coverage((ref, start, end, NO_STRAND, val), multiplier=1.0)
            correctcov = correct1[start:end].sum(axis=1) / float(total_cov)
            self.assertTrue(np.allclose(mycov, correctcov, atol=1e-4))
            mydens = t.density((ref, start, end, NO_STRAND, val), multiplier=1.0)
            correctdens = correctcount / float(total_cov * (end - start))
            self.assertTrue(np.allclose(mydens, correctdens, atol=1e-4))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()