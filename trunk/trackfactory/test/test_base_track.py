'''
Created on Mar 18, 2011

@author: mkiyer
'''
import unittest
import tempfile
import os
import numpy as np

from trackfactory.track import TrackError, POS_STRAND, NEG_STRAND, NO_STRAND
from trackfactory.trackfactory import TrackFactory
from trackfactory.io.cinterval import Interval
from trackfactory.track import _parse_interval_string, _parse_interval_tuple

    
class TestTrack(unittest.TestCase):

    def test_parse_interval(self):
        # test intervals with only strand information
        self.assertEqual(_parse_interval_string("+"),
                         (None, None, None, POS_STRAND))
        self.assertEqual(_parse_interval_string("-"),
                         (None, None, None, NEG_STRAND))
        self.assertEqual(_parse_interval_string("."),
                         (None, None, None, NO_STRAND))
        # test intervals with only reference info
        self.assertEqual(_parse_interval_string("1"),
                         ("1", None, None, NO_STRAND))
        self.assertEqual(_parse_interval_string("chr1"),
                         ("chr1", None, None, NO_STRAND))
        self.assertEqual(_parse_interval_string("asdfasdf"),
                         ("asdfasdf", None, None, NO_STRAND))
        # test intervals with reference and strand
        self.assertEqual(_parse_interval_string("chrX[+]"),
                         ("chrX", None, None, POS_STRAND))
        self.assertEqual(_parse_interval_string("1[-]"),
                         ("1", None, None, NEG_STRAND))
        self.assertEqual(_parse_interval_string("gene_1234[.]"),
                         ("gene_1234", None, None, NO_STRAND))
    


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()