'''
Created on Nov 2, 2010

@author: mkiyer
'''
GTF_EMPTY_FIELD = '.'
GTF_ATTR_SEP = ';'
GTF_ATTR_TAGVALUE_SEP = ' '

class Feature(object):
    '''
    1. seqname - The name of the sequence. Must be a chromosome or scaffold.
    2. source - The program that generated this feature.
    3. feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
    4. start - The starting position of the feature in the sequence. The first base is numbered 1.
    5. end - The ending position of the feature (inclusive).
    6. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
    7. strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
    8. phase - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.

    chr1    Cufflinks       transcript      136546  137059  1000    .       .       gene_id "VCAP_SHEZH2.657699"; transcript_id "VCAP_SHEZH2.657699.1"; FPKM "100.7219943204"; frac "1.000000"; conf_lo "80.649925"; conf_hi "120.794064"; cov "2.198209";
    '''
    __slots__ = ('ref', 'source', 'feature_type', 'start', 'end', 
                 'score', 'strand', 'gene_id', 'transcript_id',
                 'exon_number', 'aliases', 'attrs')

    def __str__(self):
        line = [self.ref,
                self.source,
                self.feature_type,
                str(self.start),
                str(self.end),
                str(self.score),
                str(self.strand),
                '.'] #self.phase]
        attr_str = ' '.join('%s "%s";' % (k, v) for (k, v) in self.attrs.iteritems())
        line.append(attr_str)
        return '\t'.join(line)

    @staticmethod
    def fromgtf(line):
        # read the GTF line
        fields = line.strip().split('\t')
        f = Feature()  
        f.ref = fields[0]
        f.source = fields[1]        
        f.feature_type = fields[2]
        f.start = int(fields[3])
        f.end = int(fields[4])
        f.score = 0 if (fields[5] == '.') else float(fields[5])
        f.strand = fields[6]
        attrs = {}
        if fields[8] != GTF_EMPTY_FIELD:
            attr_strings = fields[8].split(GTF_ATTR_SEP)
            for a in attr_strings:
                a = a.strip()
                if len(a) == 0:
                    continue
                tag, value = a.split(GTF_ATTR_TAGVALUE_SEP, 1)
                # remove quotes
                value = value.split('"')[1]
                attrs[tag] = value
                # apply parsing function
                #if (attr_defs != None) and (tag in attr_defs) and (attr_defs[tag] != None):
                #    value = attr_defs[tag](value)
        f.attrs = attrs
        f.gene_id = attrs['gene_id'] 
        f.transcript_id = attrs['transcript_id']
        f.exon_number = attrs['exon_number'] 
        f.aliases = [('gene_name', attrs['gene_name'])]
        #for alias_id in alias_attrs:
        #    if alias_id in attrs:
        #        f.aliases.append((alias_id, attrs[alias_id]))
        return f

    @staticmethod
    def frombed(self, line, attr_defs=None):
        pass



def parse_gtf_features(line_iter):
    for line in line_iter:
        if not line:
            continue
        if not line.strip():
            continue            
        if line.startswith("#"):
            continue
        fields = line.strip().split('\t')
        # only allow 'exon' type
        if fields[2] != "exon":
            continue
        yield Feature.fromgtf(line)
