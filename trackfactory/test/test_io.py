'''
Created on Mar 5, 2011

@author: mkiyer
'''
import unittest
import itertools
import numpy as np

from trackfactory.io.sequence import parse_fasta_as_chunks
from trackfactory.io.interval import write_interval_data_to_array
from trackfactory.io.cinterval import Interval

class TestIO(unittest.TestCase):

    def test_fasta_chunks(self):
        """Testing chunked fasta parser"""
        seqlist = [">hello",
                   "ATGCAGTGACGTGACGAGAGTGTAGAGAGAGTGATGTATG",
                   "GGGGGGGGGGGGGCCCCCGCCCCCCCCCCCCCCCCCCCCC",
                   "ATGAAAAAAGTTGCC",
                   "AAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGA",
                   "TTTTTTTTTGCCGAGCGCGCGCGCGCG"]
        fullseq = ''.join(seqlist[1:])
        for chunksize in xrange(0, 50):
            newseq = ''
            linenum = 0
            for tag,start,end,seq in parse_fasta_as_chunks(iter(seqlist), 
                                                           chunksize=0):            
                #print tag, start, end, seq
                self.assertEqual(tag, "hello")
                self.assertEqual(fullseq[start:end], seq)
                newseq += seq
            self.assertEqual(fullseq, newseq)
        # test multiple sequences
        seqlist = [">0",
                   "ATGCAGTGACGTGACGAGAGTGTAGAGAGAGTGATGTATG",
                   "GGGGGGGGGGGGGCCCCCGCCCCCCCCCCCCCCCCCCCCC",
                   ">1",
                   "ATGAAAAAAGTTGCC",
                   "AAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGA",
                   "TTTTTTTTTGCCGAGCGCGCGCGCGCG",
                   ">2 asdf asdf asdf",
                   "acgtacgataAGTACAgacNANACNtgatNACNagt",
                   "NNNNNNNNNNNNNNtgtgacaTGTGTACACAGTGT"]
        fullseqs = [''.join(seqlist[1:3]),
                    ''.join(seqlist[4:7]),
                    ''.join(seqlist[8:10])]
        for chunksize in xrange(0, 100):
            for tag,start,end,seq in parse_fasta_as_chunks(iter(seqlist), 
                                                           chunksize=0):
                seqnum = int(tag.split(None, 1)[0])
                self.assertEqual(fullseqs[seqnum][start:end], seq)
        
    def test_interval_to_array(self):        
        """Testing interval to array code"""
        ref = "chr1"        
        endpos = 1000        
        chunkstep = 10
        dtype = "i"
        # test different interval sizes
        for intervalsize in xrange(1,100,10):
            fullarr = np.zeros((endpos+intervalsize,1), dtype=dtype)
            intervals = []
            val = 0
            for start in xrange(0, endpos, intervalsize):
                intervals.append(Interval(ref, start, start+intervalsize, "+", val))
                fullarr[start:start+intervalsize,0] = val
                val += 2
            # test different chunk sizes
            for chunksize in xrange(chunkstep, endpos, chunkstep):
                testarr = np.zeros((endpos+intervalsize,1), dtype=dtype)
                write_interval_data_to_array(iter(intervals), 
                                             {"chr1": testarr}, 
                                             dtype=dtype, 
                                             channel=0,
                                             chunksize=chunksize,
                                             mode="channel")
                self.assertTrue(np.all(testarr == fullarr))   
        #
        # test intervals on different chromosomes
        #
        refiter = itertools.cycle(itertools.chain(itertools.repeat("chr1", 3),
                                                  itertools.repeat("chr2", 3) ,                                       
                                                  itertools.repeat("chr3", 3)))                                        
        for intervalsize in (1,10,50,100):
            #print "intervalsize", intervalsize
            fullarr = {"chr1":np.zeros((endpos+intervalsize,1), dtype=dtype),
                       "chr2":np.zeros((endpos+intervalsize,1), dtype=dtype),
                       "chr3":np.zeros((endpos+intervalsize,1), dtype=dtype)}                       
            intervals = []
            val = 0
            for start in xrange(0, endpos, intervalsize):
                ref = refiter.next()
                intervals.append(Interval(ref, start, start+intervalsize, "+", val))
                fullarr[ref][start:start+intervalsize] = val
                val += 2
            # test different chunk sizes
            for chunksize in (1, 8, 16, 32, 64, 128, 256, 512, 1024):
                testarr = {"chr1":np.zeros((endpos+intervalsize,1), dtype=dtype),
                           "chr2":np.zeros((endpos+intervalsize,1), dtype=dtype),
                           "chr3":np.zeros((endpos+intervalsize,1), dtype=dtype)}
                write_interval_data_to_array(iter(intervals),
                                             testarr, 
                                             dtype=dtype, 
                                             channel=0,
                                             chunksize=chunksize,
                                             mode="channel")
                for chrom in testarr:
                    self.assertTrue(np.all(testarr[chrom] == fullarr[chrom]))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()