'''
Created on Mar 10, 2011

@author: mkiyer
'''
from trackfactory.io.feature import parse_gtf_features
import sys

chrom_xref = {}
for x in xrange(1, 23):
    chrom_xref[str(x)] = "chr%d" % (x)
chrom_xref["MT"] = "chrM"
chrom_xref["X"] = "chrX"
chrom_xref["Y"] = "chrY"

for f in parse_gtf_features(open(sys.argv[1])):
    f.source = "Ensembl"
    if f.ref not in chrom_xref:
        print >>sys.stderr, "SKIPPING", f.ref
        continue
    f.ref = chrom_xref[f.ref]
    print str(f)
    