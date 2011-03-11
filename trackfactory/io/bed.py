'''
Created on Mar 8, 2011

@author: mkiyer
'''
from io.cinterval import BedInterval, strand_str_to_int

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


#def parse_bed12_line(line):
#    if line is None:
#        return None
#    line = line.strip()
#    if line.startswith('#'):
#        logging.debug("skipping comment line: %s" % (line))
#        return None
#    if line.startswith('track'):
#        logging.debug("skipping track header line: %s"  % (line))
#        return None
#    thisfields = line.split('\t')
#    # first six fields are required
#    g = BEDGene()
#    g.chrom = thisfields[0]
#    g.tx_start = int(thisfields[1])
#    g.tx_end = int(thisfields[2])
#    g.name = thisfields[3]
#    if len(thisfields) <= 4:
#        g.score = 0
#        g.strand = '.'
#    else:
#        g.score = thisfields[4]
#        g.strand = thisfields[5]        
#    if len(thisfields) <= 6:
#        g.cds_start = g.tx_start
#        g.cds_end = g.tx_end
#        g.exon_count = 1
#        g.exons = [(g.tx_start, g.tx_end)]
#        g.introns = []
#    else:
#        g.cds_start = int(thisfields[6])
#        g.cds_end = int(thisfields[7])
#        g.exon_count = int(thisfields[9])
#        block_sizes = map(int, thisfields[10].split(',')[:-1])
#        block_starts = map(int, thisfields[11].split(',')[:-1])        
#        g.exon_starts = [(g.tx_start + start) for start in block_starts]        
#        g.exon_ends = [(start + size) for start, size in zip(g.exon_starts, block_sizes)]
#        g.exons = zip(g.exon_starts, g.exon_ends)
#        g.introns = zip(g.exon_ends, g.exon_starts[1:])        
#    return g