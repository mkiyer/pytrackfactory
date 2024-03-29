'''
.. module:: sequencetrack
   :synopsis: track for storing 1-D arrays
.. moduleauthor:: Matthew Iyer <mkiyer@umich.edu>
.. versionadded:: 0.0.1

Created on Sep 14, 2010

@author: mkiyer
'''
import logging
#import collections
#import operator
from bitarray import bitarray

from track import TrackError, NEG_STRAND
from bitarraytrack import BitArrayTrack
from io.sequence import parse_fasta_as_chunks, dna_reverse_complement

_encode_dict_2bpb = {'N': bitarray('00'),
                     'n': bitarray('00'),
                     'A': bitarray('00'),
                     'G': bitarray('01'),
                     'C': bitarray('10'),
                     'T': bitarray('11'),
                     'a': bitarray('00'),
                     'g': bitarray('01'),
                     'c': bitarray('10'),
                     't': bitarray('11')}
_decode_dict_2bpb = {'A': bitarray('00'),
                     'G': bitarray('01'),
                     'C': bitarray('10'),
                     'T': bitarray('11')}

_encode_dict_3bpb = {'N': bitarray('000'),
                     'n': bitarray('000'),
                     'A': bitarray('100'),
                     'G': bitarray('101'),
                     'C': bitarray('110'),
                     'T': bitarray('111'),
                     'a': bitarray('100'),
                     'g': bitarray('101'),
                     'c': bitarray('110'),
                     't': bitarray('111')}
_decode_dict_3bpb = {'N': bitarray('000'),
                     'A': bitarray('100'),
                     'G': bitarray('101'),
                     'C': bitarray('110'),
                     'T': bitarray('111')}

_encode_dict_4bpb = {'N': bitarray('0000'),
                     'n': bitarray('0100'),
                     'A': bitarray('1000'),
                     'G': bitarray('1001'),
                     'C': bitarray('1010'),
                     'T': bitarray('1011'),
                     'a': bitarray('1100'),
                     'g': bitarray('1101'),
                     'c': bitarray('1110'),
                     't': bitarray('1111')}

_encoders = {2: _encode_dict_2bpb,
             3: _encode_dict_3bpb,
             4: _encode_dict_4bpb}
_decoders = {2: _decode_dict_2bpb,
             3: _decode_dict_3bpb,
             4: _encode_dict_4bpb}

class SequenceTrack(BitArrayTrack):
    """Stores DNA/RNA sequence information in compact bitarrays
    
    Three formats are currently available via the `bpb` parameter
    
    2 bits per base: Stores only 'A', 'C', 'G', and 'T'.  
    'N' bases are converted to 'A'.
    
    3 bits per base: Stores 'N', 'A', 'C', 'G', and 'T'.  
    
    4 bits per base (default): Stores 'N', 'actg', and 'ACGT' such that 
    lower/upper case base information is retained.
    
    :param bpb: 2, 3, or 4 
    """        
    def __init__(self, hdf_group, bpb=4):
        super(SequenceTrack, self).__init__(hdf_group, bpb)

    def _read(self, arr, start, end, strand):
        # get bitarray
        b = BitArrayTrack._read(self, arr, start, end)
        # decode to sequence
        codec_dict = _decoders[self.bpb]
        seq = ''.join(b.decode(codec_dict))
        if strand == NEG_STRAND:
            seq = dna_reverse_complement(seq)
        return seq

    def _write(self, arr, start, end, strand, value):
        if strand == NEG_STRAND:
            value = dna_reverse_complement(value)
        # encode to bitarray
        b = bitarray()
        b.encode(_encoders[self.bpb], value)        
        # write
        BitArrayTrack._write(self, arr, start, end, b)

    def __getitem__(self, key):
        ref, start, end, strand = self._parse_interval(key)
        arr = self._get_array(ref)
        return self._read(arr, start, end, strand)

    def __setitem__(self, key, value):
        ref, start, end, strand = self._parse_interval(key)
        arr = self._get_array(ref)
        self._write(arr, start, end, strand, value)

    def fromfasta(self, line_iter, ref=None, start=None, end=None,
                  split_tag=True):
        """Write a sequence (fasta) file into the track
        
        If no fasta header is provided, the user must specify the `ref`, 
        `start`, and `end` arguments to the function.  Otherwise, the
        fasta header must be formatted in `>ref[strand]:start-end` string 
        format and will be parsed accordingly.
        
        You can use this function to import entire chromosomes (or genomes)
        as inserting into the track is done in a efficient, chunked fashion
        """
        for tag, tag_start, tag_end, seq in parse_fasta_as_chunks(line_iter):
            logging.debug("%s:%d-%d" % (tag, tag_start, tag_end))
            if ref is None:
                if split_tag:
                    tag = tag.split(None,1)[0]
                try:
                    ref, arr_start, arr_end, strand = \
                        self._parse_interval(tag)
                except TrackError as e:                    
                    logging.warning("[SequenceTrack] %s" % (str(e)))
                    continue
            else:
                ref, arr_start, arr_end, strand = self._parse_interval((ref, start, end))
            write_start = arr_start + tag_start
            write_end = arr_start + tag_end
            arr = self._get_array(ref)
            self._write(arr, write_start, write_end, strand, seq)
