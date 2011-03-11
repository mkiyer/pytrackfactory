'''
Created on Sep 29, 2010

@author: mkiyer
'''
import unittest
import numpy as np
import tempfile
import os
import StringIO

from trackfactory.io.cbedgraph import array_to_bedgraph
#from trackfactory.io.bedgraph import array_to_bedgraph
from trackfactory.io.bedgraph import bedgraph_to_array

def interval(ref, start, end, cov):
    return "%s\t%d\t%d\t%.2f\n" % (ref, start, end, cov)

def make_temp():
    fd, filename = tempfile.mkstemp(suffix='.bedgraph', prefix='tmp')
    os.close(fd)
    fh = open(filename, "w+")
    return fh, filename

def getvalue(fh):
    fh.seek(0)
    return bedgraph_to_array(fh)
    #return ''.join(fh.readlines())

class TestBedGraph(unittest.TestCase):

    def test_array_dtypes(self):
        # float
        a = np.array([1.25, 2.01, 3.98, 10.1, 4.02], dtype=np.float32)
        a.resize((a.shape[0],1))
        output, filename = make_temp()
        array_to_bedgraph('chr1', a, output)
        contents = getvalue(output)
        self.assertTrue(np.array_equal(contents['chr1'], a[:,0]))
        os.remove(filename)
        # int
        b = np.array([1.25, 2.01, 3.98, 10.1, 4.02], dtype=np.int32)
        b.resize((b.shape[0],1))
        bcorrect = np.array([1, 2, 3, 10, 4], dtype=np.int32)
        bcorrect.resize((b.shape[0],1))
        output, filename = make_temp()
        array_to_bedgraph('chr1', b, output)
        contents = getvalue(output)
        self.assertTrue(np.array_equal(contents['chr1'], bcorrect[:,0]))
        os.remove(filename)

    def test_array_to_bedgraph_zeros(self):
        a = np.array([0, 0, 0, 0, 0])
        a.resize((a.shape[0],1))
        output, filename = make_temp()
        array_to_bedgraph('chr1', a, output)
        contents = getvalue(output)
        self.assertTrue(np.array_equal(contents['chr1'], a[:,0]))
        os.remove(filename)        

    def test_array_to_bedgraph_jumps(self):
        a = np.array([0, 0, 1, 1, 0, 0, 0, 0, 10, 10, 0, 0], dtype=np.float)
        a.resize((a.shape[0],1))
        output, filename = make_temp()
        array_to_bedgraph('chr1', a, output)
        contents = getvalue(output)
        self.assertTrue(np.array_equal(contents['chr1'], a[:,0]))
        os.remove(filename)

    def test_array_to_bedgraph_consecutive(self):
        a = np.arange(0, 100, dtype=np.float)
        a.resize((a.shape[0],1))
        output, filename = make_temp()
        array_to_bedgraph('chr1', a, output)
        contents = getvalue(output)
        self.assertTrue(np.array_equal(contents['chr1'], a[:,0]))
        os.remove(filename)

    def test_array_to_bedgraph_sparse(self):
        a = np.zeros((10000,1), dtype=np.float)
        a[6655] = 100
        output, filename = make_temp()
        array_to_bedgraph('chr1', a, output)
        contents = getvalue(output)   
        self.assertTrue(np.array_equal(contents['chr1'], a[:,0]))
        os.remove(filename)

    def test_array_to_bedgraph_sparse_spans(self):
        a = np.zeros((10000,1), dtype=np.float)
        a[6655] = 100
        # check that spans work too
        output, filename = make_temp()
        array_to_bedgraph('chr1', a, output, span=10)
        contents = getvalue(output)
        correct = np.zeros(10000, dtype=np.float)
        correct[6650:6660] = 10
        self.assertTrue(np.array_equal(contents['chr1'], correct))
        os.remove(filename)

    def test_array_to_bedgraph_span_zeros(self):
        a = np.zeros((10000,1), dtype=np.float)
        output, filename = make_temp()
        array_to_bedgraph('chr1', a, output, span=1)
        contents = getvalue(output)        
        self.assertTrue(np.array_equal(contents['chr1'], a[:,0]))
        os.remove(filename)

    def test_array_to_bedgraph_span_jumps(self):
        a = np.array([0, 0, 2, 2, 0, 0, 4, 4, 0, 0, 0, 0, 10, 10, 10, 10], dtype=np.float)
        a.resize((a.shape[0],1))
        # span of 2
        output, filename = make_temp()
        array_to_bedgraph('chr1', a, output, span=2)
        contents = getvalue(output)
        self.assertTrue(np.array_equal(contents['chr1'], a[:,0]))
        os.remove(filename)
        # span of 4
        output, filename = make_temp()
        array_to_bedgraph('chr1', a, output, span=4)
        contents = getvalue(output)
        correct = np.array([1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 10, 10, 10, 10])
        self.assertTrue(np.array_equal(contents['chr1'], correct))
        os.remove(filename)

    def test_array_to_bedgraph_span_uneven(self):
        a = np.array([0, 0, 2, 2, 0, 0, 4, 4, 0, 0, 0, 0, 10, 10, 10, 10, 999], dtype=np.float)
        a.resize((a.shape[0],1))
        # span of 2
        output, filename = make_temp()
        array_to_bedgraph('chr1', a, output, span=2)
        contents = getvalue(output)
        self.assertTrue(np.array_equal(contents['chr1'], a[:,0]))
        os.remove(filename)
        # span of 4
        a = np.array([0, 0, 2, 2, 0, 0, 4, 4, 0, 0, 0, 0, 10, 10, 10, 10, 999, 998, 997])
        a.resize((a.shape[0],1))
        output, filename = make_temp()
        array_to_bedgraph('chr1', a, output, span=4)
        contents = getvalue(output)
        correct = np.array([1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 10, 10, 10, 10, 998, 998, 998])
        self.assertTrue(np.array_equal(contents['chr1'], correct))
        os.remove(filename)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()