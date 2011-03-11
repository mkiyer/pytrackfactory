'''
Created on Mar 9, 2011

@author: mkiyer
'''
import unittest
import tempfile
import os
import numpy as np

from trackfactory.track import TrackError
from trackfactory.trackfactory import TrackFactory
from trackfactory.arraytrack import ArrayTrack
from trackfactory.io.cinterval import Interval

def mktemp(prefix, suffix):
    fh, filename = tempfile.mkstemp(suffix=suffix, prefix=prefix)
    os.close(fh)
    return filename
    
class TestArrayTrack(unittest.TestCase):
    def setUp(self):
        filename = mktemp(prefix="tmp", suffix=".h5")
        self.filename = filename
        self.refs = (('gene1', 100000), ('gene2', 10))
        self.tf = TrackFactory(self.filename, 'w', refs=self.refs) 

    def tearDown(self):
        self.tf.close()
        if os.path.exists(self.filename):
            os.remove(self.filename)



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()