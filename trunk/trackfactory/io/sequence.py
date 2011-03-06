'''
Created on Sep 8, 2010

@author: mkiyer
'''
from cStringIO import StringIO
from string import maketrans

_dna_compl_translate = maketrans('acgtACGT', 'tgcaTGCA')

def complement(s):
    return s.translate(_dna_compl_translate)
def reverse_complement(s):
    return s.translate(_dna_compl_translate)[::-1]

def parse_fasta_as_chunks(line_iter, chunksize=(1 << 20)):
    """
    Parses a multi-line sequence file in FASTA format and 
    return chunks of sequence whenever the parser reads
    more than `chunksize` bases of the file.  Chunks are not
    guaranteed to be fixed in size.
    
    :returns: iterator of (tag, start, end, sequence) tuples
    """
    # initialize state
    chunksize = max(1, chunksize)
    tag = None
    start = 0
    end = 0    
    seqbuf = StringIO()
    try:
        while True:
            line = line_iter.next().strip()           
            if line.startswith('>'):
                seq = seqbuf.getvalue()
                if len(seq) > 0:                
                    yield tag, start, end, seq 
                # clear sequence buffer
                seqbuf.close()
                seqbuf = StringIO()
                start = end
                # get new tag
                tag = line[1:]
                start = 0
                end = 0
            else:
                if (end - start) >= chunksize:
                    seq = seqbuf.getvalue()
                    if len(seq) > 0:                
                        yield tag, start, end, seq
                    # clear sequence buffer
                    seqbuf.close()
                    seqbuf = StringIO()
                    start = end
                seqbuf.write(line)        
                end += len(line)
    except StopIteration:
        pass
    seq = seqbuf.getvalue()
    if len(seq) > 0:
        yield tag, start, end, seq
