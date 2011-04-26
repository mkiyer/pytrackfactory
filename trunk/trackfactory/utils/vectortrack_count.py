'''
Created on Mar 20, 2011

@author: mkiyer
'''
import argparse
import logging
import collections

from trackfactory.io.bed import BedGene
from trackfactory import TrackFactory
from trackfactory.track import POS_STRAND, NEG_STRAND

from bx.intervals.intersection import Interval, IntervalTree

def filter_strand_conflicts(bedfile):
    trees = collections.defaultdict(lambda: IntervalTree())    
    genes = []     
    for g in BedGene.parse(open(bedfile)):
        for start, end in g.exons:   
            trees[g.chrom].insert_interval(Interval(start, end, strand=g.strand, chrom=g.chrom))
        genes.append(g)
    exon_count = 0
    new_exon_count = 0
    for g in genes:
        new_exons = []
        for start, end in g.exons:
            exon_count += 1
            ambig = False
            for hit in trees[g.chrom].find(start,end):
                if hit.strand != g.strand:
                    ambig = True
                    logging.debug("Filtering ambiguous exon %s:%d-%d" % (g.chrom, start, end))
                    break
            if not ambig:
                new_exon_count += 1
                new_exons.append((start,end))
        g.exons = new_exons
    logging.info("%d exons" % (exon_count))
    logging.info("%d unambiguous exons" % (new_exon_count))         
    return genes

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--stranded", dest="stranded", action="store_true", default=False)
    parser.add_argument("--ambiguous", dest="ambiguous", action="store_true", default=False)
    parser.add_argument("--aliases", dest="alias_file", default=None)
    parser.add_argument("bed")
    parser.add_argument("track_files", nargs="+")
    options = parser.parse_args()
    
    alias_dict = {}
    alias_header = []
    if options.alias_file is not None:
        for line in open(options.alias_file):
            if line.startswith("#"):
                alias_header = line.strip()[1:].split('\t')
                continue    
            fields = line.strip().split('\t')
            alias_dict[fields[0]] = fields[1:]
    
    header_fields = alias_header + ["gene_name", "gene_interval", "gene_length"]
    tracks = []
    for track_path in options.track_files:
        name, path = track_path.split("@")
        file_path, h5_path = path.split(":")
        tf = TrackFactory(file_path, "r")
        t = tf.get_track(h5_path)
        tracks.append((name, tf, t, set(t.get_rnames())))
        if options.stranded:
            header_fields.append("%s_sense" % name)
            header_fields.append("%s_antisense" % name)
        else:
            header_fields.append(name)
    # output header
    print '\t'.join(map(str, header_fields))

    # read genes
    if options.ambiguous:
        genes = list(BedGene.parse(open(options.bed)))
    else:
        genes = filter_strand_conflicts(options.bed)
    # get counts
    for g in genes:
        alias_fields = alias_dict.get(g.name, ["None"] * len(alias_header))
        fields = ([g.name] + alias_fields +
                  ["%s[%s]:%d-%d" % (g.chrom, g.strand, g.tx_start, g.tx_end),
                   sum((end-start) for start,end in g.exons)])
        sense_strand = NEG_STRAND if g.strand == "+" else POS_STRAND
        antisense_strand = int(not sense_strand)
        rname_found = False
        for name, tf, t, rnames in tracks:
            if g.chrom not in rnames:
                continue
            rname_found = True        
            if options.stranded:
                sense_count = 0
                antisense_count = 0
                for start, end in g.exons:
                    sense_count += t.count((g.chrom, start, end, sense_strand))
                    antisense_count += t.count((g.chrom, start, end, antisense_strand))
                fields.append(sense_count)
                fields.append(antisense_count)
            else:
                count = 0
                for start, end in g.exons:
                    count += t.count((g.chrom, start, end))
                fields.append(count)
        if rname_found:
            print '\t'.join(map(str, fields))

    for name,tf,t,rnames in tracks:
        tf.close()
    
if __name__ == '__main__':
    main()
