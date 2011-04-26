'''
Created on Mar 8, 2011

@author: mkiyer
'''
from cinterval import BedInterval, strand_str_to_int

def parse_bed6(line_iter, sep="\t"):
    for line in line_iter:
        if not line:
            continue
        line = line.strip()
        if not line:
            continue
        if line[0] == "#":
            continue
        if line.startswith("track") or line.startswith("browser"):
            continue
        fields = line.split(sep)
        ref = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3]        
        score = int(fields[4])
        strand = strand_str_to_int(fields[5])
        yield BedInterval(ref, start, end, strand=strand, value=score,
                          name=name) 


class BedGene(object):
    __slots__ = ('chrom', 'tx_start', 'tx_end', 'name', 'score', 'strand',
                 'cds_start', 'cds_end', 'exon_count', 'exons', 'introns')
    
    @staticmethod
    def parse_line(line):
        if line is None:
            return None
        line = line.strip()
        if line.startswith('#'):
            return None
        if line.startswith('track'):
            return None
        fields = line.split('\t')
        g = BedGene()
        g.chrom = fields[0]
        g.tx_start = int(fields[1])
        g.tx_end = int(fields[2])
        g.name = fields[3]
        g.score = fields[4]
        g.strand = fields[5]
        g.cds_start = int(fields[6])
        g.cds_end = int(fields[7])
        g.exon_count = int(fields[9])
        block_sizes = map(int, fields[10].split(',')[:-1])
        block_starts = map(int, fields[11].split(',')[:-1])        
        exon_starts = [(g.tx_start + start) for start in block_starts]        
        exon_ends = [(start + size) for start, size in zip(exon_starts, block_sizes)]
        g.exons = zip(exon_starts, exon_ends)
        g.introns = zip(exon_ends, exon_starts[1:])        
        return g
    
    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            g = BedGene.parse_line(line)
            if g is None:
                continue
            yield g
