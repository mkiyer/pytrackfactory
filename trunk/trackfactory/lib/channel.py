'''
Created on Mar 18, 2011

@author: mkiyer
'''
from trackfactory.track import NO_STRAND, POS_STRAND, NEG_STRAND

_allele_dict = {"A": (0,),
                "G": (1,),
                "C": (2,),
                "T": (3,),
                "N": (0,1,2,3),
                None: (0,1,2,3)}

def _select_pe_channel(pe, strand, allele):
    if ((pe is None) or (pe < 0)): return (0,1)
    return (pe,)
def _select_strand_channel(pe, strand, allele):
    if strand == NO_STRAND: return (0,1)
    return (strand,)
def _select_allele_channel(pe, strand, allele):
    if allele is None: return (0,1,2,3)
    return _allele_dict[allele]
def _select_pe_strand_channel(pe, strand, allele):
    pe_channels = _select_pe_channel(pe, strand, allele)
    strand_channels = _select_strand_channel(pe, strand, allele)
    channels = []
    for x in pe_channels:
        for y in strand_channels:
            channels.append((y << 1) | (x))
    return tuple(sorted(channels))
def _select_pe_allele_channel(pe, strand, allele):
    pe_channels = _select_pe_channel(pe, strand, allele)
    allele_channels = _select_allele_channel(pe, strand, allele)
    channels = []
    for x in pe_channels:
        for y in allele_channels:
            channels.append((y << 1) | (x))
    return tuple(sorted(channels))
def _select_strand_allele_channel(pe, strand, allele):
    strand_channels = _select_strand_channel(pe, strand, allele)
    allele_channels = _select_allele_channel(pe, strand, allele)
    channels = []
    for x in strand_channels:
        for y in allele_channels:
            channels.append((y << 1) | (x))
    return tuple(sorted(channels))
def _select_pe_strand_allele_channel(pe, strand, allele):
    pe_strand_channels = _select_pe_strand_channel(pe, strand, allele)
    allele_channels = _select_allele_channel(pe, strand, allele)
    channels = []
    for x in pe_strand_channels:
        for y in allele_channels:
            channels.append((y << 2) | (x))
    return tuple(sorted(channels))

_channel_dict = {(False, False, False): lambda pe,strand,allele: (0,),
                 (True, False, False): _select_pe_channel,
                 (False, True, False): _select_strand_channel,
                 (False, False, True): _select_allele_channel,
                 (True, True, False): _select_pe_strand_channel,
                 (True, False, True): _select_pe_allele_channel,
                 (False, True, True): _select_strand_allele_channel,
                 (True, True, True): _select_pe_strand_allele_channel}

def get_channel_dict(is_pe=False, is_strand=False, is_allele=False):
    selector = _channel_dict[(is_pe, is_strand, is_allele)]
    d = {}
    for read in (None,-1,0,1):
        for strand in (NO_STRAND, POS_STRAND, NEG_STRAND):
            for allele in _allele_dict.keys():
                d[(read, strand, allele)] = selector(read, strand, allele)
    return d




