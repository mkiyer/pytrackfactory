'''
Created on Mar 5, 2011

@author: mkiyer
'''
import unittest
import itertools
import numpy as np

from trackfactory.io.sequence import parse_fasta_as_chunks
from trackfactory.io.intervals import intervals_to_array
from trackfactory.io.cintervals import IntervalToArrayChunks

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
            fullarr = np.zeros(endpos+intervalsize, dtype=dtype)
            intervals = []
            val = 0
            for start in xrange(0, endpos, intervalsize):
                intervals.append((ref, start, start+intervalsize, "+", val))
                fullarr[start:start+intervalsize] = val
                val += 2
            # test different chunk sizes
            testarr = np.zeros(endpos+intervalsize, dtype=dtype)
            for chunksize in xrange(chunkstep, endpos, chunkstep):
                for chrom, start, end, arr in intervals_to_array(iter(intervals), "i", chunksize):
                    #print chrom, start, end, arr, fullarr[start:end]
                    testarr[start:end] = arr                    
                    self.assertTrue(np.all(arr == fullarr[start:end]))
            self.assertTrue(np.all(testarr == fullarr))
        # test intervals on different chromosomes
        refiter = itertools.cycle(itertools.chain(itertools.repeat("chr1", 3),
                                                  itertools.repeat("chr2", 3) ,                                       
                                                  itertools.repeat("chr3", 3)))                                        
        for intervalsize in (1,10,50,100):
            #print "intervalsize", intervalsize
            fullarr = {"chr1":np.zeros(endpos+intervalsize, dtype=dtype),
                       "chr2":np.zeros(endpos+intervalsize, dtype=dtype),
                       "chr3":np.zeros(endpos+intervalsize, dtype=dtype)}                       
            intervals = []
            val = 0
            for start in xrange(0, endpos, intervalsize):
                ref = refiter.next()
                intervals.append((ref, start, start+intervalsize, "+", val))
                fullarr[ref][start:start+intervalsize] = val
                val += 2
            # test different chunk sizes
            testarr = {"chr1":np.zeros(endpos+intervalsize, dtype=dtype),
                       "chr2":np.zeros(endpos+intervalsize, dtype=dtype),
                       "chr3":np.zeros(endpos+intervalsize, dtype=dtype)}
            for chunksize in (1, 8, 16, 32, 64, 128, 256, 512, 1024):
                #print "chunksize", chunksize
                for chrom, start, end, arr in intervals_to_array(iter(intervals), "i", chunksize):
                    #print chrom, start, end, arr, fullarr[chrom][start:end]
                    testarr[chrom][start:end] = arr                    
                    self.assertTrue(np.all(arr == fullarr[chrom][start:end]))
            for chrom in testarr:
                self.assertTrue(np.all(testarr[chrom] == fullarr[chrom]))
        
    def test_cinterval_to_array(self):        
        """Testing interval to array code"""
        ref = "chr1"        
        endpos = 1000        
        chunkstep = 10
        dtype = "i"
        # test different interval sizes
        for intervalsize in xrange(1,100,10):
            fullarr = np.zeros(endpos+intervalsize, dtype=dtype)
            intervals = []
            val = 0
            for start in xrange(0, endpos, intervalsize):
                intervals.append((ref, start, start+intervalsize, "+", val))
                fullarr[start:start+intervalsize] = val
                val += 2
            # test different chunk sizes
            testarr = np.zeros(endpos+intervalsize, dtype=dtype)
            for chunksize in xrange(chunkstep, endpos, chunkstep):
                buf = np.zeros(chunksize, dtype="i")
                array_iter = IntervalToArrayChunks(iter(intervals), buf)
                for chrom, start, end, arr in array_iter:
                    #print chrom, start, end, arr, fullarr[start:end]
                    testarr[start:end] = arr
                    self.assertTrue(np.all(arr == fullarr[start:end]))
            self.assertTrue(np.all(testarr == fullarr))
        # test intervals on different chromosomes
        refiter = itertools.cycle(itertools.chain(itertools.repeat("chr1", 3),
                                                  itertools.repeat("chr2", 3) ,                                       
                                                  itertools.repeat("chr3", 3)))                                        
        for intervalsize in (1,10,50,100):
            #print "intervalsize", intervalsize
            fullarr = {"chr1":np.zeros(endpos+intervalsize, dtype=dtype),
                       "chr2":np.zeros(endpos+intervalsize, dtype=dtype),
                       "chr3":np.zeros(endpos+intervalsize, dtype=dtype)}                       
            intervals = []
            val = 0
            for start in xrange(0, endpos, intervalsize):
                ref = refiter.next()
                intervals.append((ref, start, start+intervalsize, "+", val))
                fullarr[ref][start:start+intervalsize] = val
                val += 2
            # test different chunk sizes
            testarr = {"chr1":np.zeros(endpos+intervalsize, dtype=dtype),
                       "chr2":np.zeros(endpos+intervalsize, dtype=dtype),
                       "chr3":np.zeros(endpos+intervalsize, dtype=dtype)}
            for chunksize in (1, 8, 16, 32, 64, 128, 256, 512, 1024):
                buf = np.zeros(chunksize, dtype="i")
                array_iter = IntervalToArrayChunks(iter(intervals), buf)
                #print "chunksize", chunksize
                for chrom, start, end, arr in array_iter:
                    #print chrom, start, end, arr, fullarr[chrom][start:end]
                    testarr[chrom][start:end] = arr                    
                    self.assertTrue(np.all(arr == fullarr[chrom][start:end]))
            for chrom in testarr:
                self.assertTrue(np.all(testarr[chrom] == fullarr[chrom]))    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()