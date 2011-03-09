'''
Created on Mar 8, 2011

@author: mkiyer
'''
import unittest

import trackfactory.lib.cintervaltree as itree

class TestIntervalTree(unittest.TestCase):

    def test_init(self):    
        tree, root = itree.init_tree(100)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()