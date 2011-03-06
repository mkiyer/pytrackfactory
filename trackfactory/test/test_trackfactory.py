'''
.. module:: test_trackfactory
   :synopsis: unit tests
.. moduleauthor:: Matthew Iyer <mkiyer@umich.edu>
.. versionadded:: 0.0.1

Created on September 8, 2010
@author: mkiyer
'''
import unittest
import tempfile
import os

from ..track import TrackError
from ..trackfactory import TrackFactory
from ..arraytrack import ArrayTrack

def mktemp(prefix, suffix):
    fh, filename = tempfile.mkstemp(suffix=suffix, prefix=prefix)
    os.close(fh)
    return filename

class TestReferences(unittest.TestCase):
    def setUp(self):
        filename = mktemp(prefix="tmp", suffix=".h5")
        self.filename = filename

    def tearDown(self):
        if os.path.exists(self.filename):
            os.remove(self.filename)

    def test_references(self):
        refs = (('gene1', 100000), ('gene2', 200))
        tf = TrackFactory(self.filename, 'w', refs=refs) 
        rnames = list(tf.get_rnames())
        lengths = list(tf.get_lengths())
        tfrefs = list(tf.get_refs())
        for i,ref in enumerate(refs):
            self.assertTrue(rnames[i] == ref[0])
            self.assertTrue(lengths[i] == ref[1])
            self.assertEqual(ref[0], tfrefs[i][0])
            self.assertEqual(ref[1], tfrefs[i][1])
        # test getting sizes of specific references
        self.assertEqual(tf.get_ref_length('gene1'), 100000)
        self.assertEqual(tf.get_ref_length('gene2'), 200)        
        tf.close()

    def test_track_factory(self):
        refs = (('gene1', 100000), ('gene2', 10))
        tf = TrackFactory(self.filename, 'w', refs=refs) 
        t = tf.create_track("a", ArrayTrack, "i4")
        t['gene2'] = range(10)        
        for i in range(10):
            self.assertEqual(t[('gene2', i)], i)        
        # assert cannot create same track twice
        self.assertRaises(TrackError, tf.create_track, "a", ArrayTrack, "i4")
        self.assertTrue(tf.has_track("a"))
        # delete and recreate track
        tf.delete_track("a")
        t = tf.create_track("a", ArrayTrack, "i4")
        t['gene2'] = range(10)
        for i in range(10):
            self.assertEqual(t[('gene2', i)], i)        
        tf.close()
        # test persistence of data
        tf = TrackFactory(self.filename)
        t = tf.get_track("a")
        for i in range(10):
            self.assertEqual(t[('gene2', i)], i)
        tf.close() 


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()