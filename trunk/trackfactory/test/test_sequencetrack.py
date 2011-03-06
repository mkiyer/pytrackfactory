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
import itertools
import collections
import numpy as np

from trackfactory.track import TrackError
from trackfactory.trackfactory import TrackFactory
from trackfactory.sequencetrack import SequenceTrack

def mktemp(prefix, suffix):
    fh, filename = tempfile.mkstemp(suffix=suffix, prefix=prefix)
    os.close(fh)
    return filename

class TestSequenceTrack(unittest.TestCase):

    def setUp(self):
        filename = mktemp(prefix="tmp", suffix=".h5")
        self.filename = filename
        self.refs = (('chr1', 100), ('chr2', 200))
        self.tf = TrackFactory(self.filename, 'w', refs=self.refs) 

    def tearDown(self):
        self.tf.close()
        if os.path.exists(self.filename):
            os.remove(self.filename)

    def test_bpb(self):
        for bpb in (2,3,4):
            tname = 't%d' % (bpb)
            t = self.tf.create_track(tname, SequenceTrack, bpb=bpb)
            # write a cycle of bases one at a time and then
            # check them
            dna_iter = itertools.cycle('AATTTGCCTGC')
            for i in xrange(100):
                t[('chr1', i)] = dna_iter.next()
            dna_iter = itertools.cycle('AATTTGCCTGC')
            for i in xrange(100):
                self.assertEqual(t[('chr1', i)], dna_iter.next())
            # write chunks of bases and check them
            dna_iter = itertools.repeat('CACATGTAGAGCT')
            for i in xrange(0, 195, 13):
                t[('chr2', i, i+13)] = dna_iter.next()            
            for i in xrange(0, 195, 13):            
                self.assertEqual(t[('chr2', i, i+13)], dna_iter.next()) 
            # overwrite a range of bases without affecting others
            t[('chr2', 99, 102)] = 'TTT'
            self.assertEqual(t[('chr2', 94,107)], 'ATGTATTTCTCAC')

    def test_bases(self):
        t2 = self.tf.create_track("t2", SequenceTrack, bpb=2)
        t3 = self.tf.create_track("t3", SequenceTrack, bpb=3)
        t4 = self.tf.create_track("t4", SequenceTrack, bpb=4)
        # test default bases
        self.assertEqual(t2[('chr1', 35, 45)], 'AAAAAAAAAA')
        self.assertEqual(t3[('chr1', 35, 45)], 'NNNNNNNNNN')
        self.assertEqual(t4[('chr1', 35, 45)], 'NNNNNNNNNN')
        # test storing lower case bases and 'N's
        seq = 'NattgcgcNN'
        t2[('chr1', 35, 45)] = seq
        t3[('chr1', 35, 45)] = seq
        t4[('chr1', 35, 45)] = seq
        self.assertEqual(t2[('chr1', 35, 45)], seq.upper().replace('N', 'A'))
        self.assertEqual(t3[('chr1', 35, 45)], seq.upper())
        self.assertEqual(t4[('chr1', 35, 45)], seq)

    def test_fromfasta(self):
        # read/write to chromosome
        seqlist = [">chr2",
                   "ATGCAGTGAC",
                   "GTGACGAGAG",
                   "TGTAGAGAGA",
                   "GTGATGTATG",
                   "GGGGGGGGGG",
                   "GGGCCCCCGC",
                   "CCCCCCCCCC"
                   "CCCCCCCCCC",
                   "ATGAAAAAAG",
                   "TTGCCAAAAA",
                   "AAAAAAAAAA",
                   "AAAAAAAAAA",
                   "AAAGTGATTT",
                   "TTTTTTGCCG",
                   "AGCGCGCGCG"]
        fullseq = ''.join(seqlist[1:])
        t2 = self.tf.create_track("t2", SequenceTrack, bpb=2)
        t2.fromfasta(iter(seqlist))
        self.assertEqual(fullseq, t2[('chr2', 0, len(fullseq))])
        # insert at different places
        for startpos in (0, 10, 20, 30, 40):
            seqlist[0] = ">chr2:%d" % (startpos)
            t2.fromfasta(iter(seqlist))
            self.assertEqual(fullseq, t2[('chr2', startpos, startpos + len(fullseq))])
            
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
