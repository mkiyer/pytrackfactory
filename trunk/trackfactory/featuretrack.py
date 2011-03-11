'''
Created on Sep 28, 2010

@author: mkiyer
'''
import collections
import tables
import numpy as np
import logging

from track import TrackError
from intervaltrack import IntervalTrack, get_base_dtype_fields, \
    REF_COL_NAME, START_COL_NAME, END_COL_NAME

EXON_TYPE = "exon"
GENE_TYPE = "gene"
TRANSCRIPT_TYPE = "transcript"

FEATURE_TABLE = "features"
FEATURE_ALIAS_TABLE = "aliases"
FEATURE_ASSOC_TABLE = "feature_assoc"
ALIAS_ASSOC_TABLE = "alias_assoc"

# the numpy dtype used to represent features 
feature_dtype = (get_base_dtype_fields() +
                 [('id', '<u4'),
                  ('name', 'a32'),
                  ('source', 'a24'), 
                  ('feature_type', 'a23'), 
                  ('strand', 'a1'),
                  ('score', '<i4')])
feature_expectedrows = 1e5

feature_assoc_dtype = [('parent_id', '<u4'),
                       ('child_id', '<u4'),
                       ('value', '<u4')]
feature_alias_dtype = [('id', '<u4'),
                       ('source', 'a16'),
                       ('name', 'a128')]
alias_assoc_dtype = [('alias_id', '<u4'),
                     ('feature_id', '<u4')]

class FeatureTrack(IntervalTrack):

    def __init__(self, hdf_group, expectedrows=None):
        if expectedrows is None:
            expectedrows = feature_expectedrows
        super(FeatureTrack, self).__init__(hdf_group, dtype=feature_dtype,
                                           expectedrows=expectedrows)
        self._init_tables(expectedrows)
        self.alias_to_feature_index = None
        self.feature_to_alias_index = None
        self.parent_id_index = None
        self.child_id_index = None

    def _init_tables(self, expectedrows):
        h5file = self._get_hdf_file()
        if not FEATURE_ASSOC_TABLE in self.hdf_group:
            h5file.createTable(self.hdf_group, FEATURE_ASSOC_TABLE,
                               np.dtype(feature_assoc_dtype),
                               expectedrows=expectedrows)            
        if not FEATURE_ALIAS_TABLE in self.hdf_group:
            h5file.createTable(self.hdf_group, FEATURE_ALIAS_TABLE, 
                               np.dtype(feature_alias_dtype), 
                               expectedrows=expectedrows) 
        if not ALIAS_ASSOC_TABLE in self.hdf_group:
            h5file.createTable(self.hdf_group, ALIAS_ASSOC_TABLE,
                               np.dtype(alias_assoc_dtype),
                               expectedrows=expectedrows)
    
    def index(self, persist=True):
        super(FeatureTrack, self).index(persist)
        self._index_aliases()
        self._index_assoc()
    
    def _index_aliases(self):
        # automatically index the aliases for fast lookup       
        alias_to_feature_index = collections.defaultdict(lambda: set())
        feature_to_alias_index = collections.defaultdict(lambda: set())
        for row in self.hdf_group.alias_assoc:
            alias_to_feature_index[row['alias_id']].add(row['feature_id'])
            feature_to_alias_index[row['feature_id']].add(row['alias_id'])
        self.alias_to_feature_index = alias_to_feature_index
        self.feature_to_alias_index = feature_to_alias_index
        return alias_to_feature_index, feature_to_alias_index

    def _index_assoc(self):
        # automatically index the feature associations for fast lookup 
        parent_ids = collections.defaultdict(lambda: set())
        child_ids = collections.defaultdict(lambda: set())
        for row in self.hdf_group.feature_assoc:            
            parent_ids[row['child_id']].add((row['parent_id'], row['value']))
            child_ids[row['parent_id']].add((row['child_id'], row['value']))
        self.child_id_index = child_ids
        self.parent_id_index = parent_ids
        return parent_ids, child_ids

    def get_feature(self, feature_id):
        tbl = self._get_interval_table()
        return tbl[feature_id]
    
    def get_feature_aliases(self, my_feature_id):
        if (self.feature_to_alias_index is not None):
            return [self.hdf_group.aliases[alias_id] for alias_id in self.feature_to_alias_index[my_feature_id]]   
        return [self.hdf_group.aliases[r['alias_id']] for r in self.hdf_group.alias_assoc.where('feature_id == my_feature_id')]
    def get_alias_feature_ids(self, alias_name):
        alias_ids = [r['id'] for r in self.hdf_group.aliases.where('name == alias_name')]
        feature_ids = set()
        for my_alias_id in alias_ids:            
            feature_ids.update([r['feature_id'] for r in self.hdf_group.alias_assoc.where('alias_id == my_alias_id')])
        return feature_ids
    def get_parents(self, feature_id):
        tbl = self._get_interval_table()
        if (self.parent_id_index is not None) and (feature_id in self.parent_id_index):
            return [(tbl[parent_id],value) for parent_id,value in self.parent_id_index[feature_id]]            
        return [(tbl[r['parent_id']],r['value']) for r in self.hdf_group.feature_assoc.where('child_id == feature_id')]
    def get_parent_ids(self, feature_id):
        if (self.parent_id_index is not None) and (feature_id in self.parent_id_index):
            return [(parent_id,value) for parent_id,value in self.parent_id_index[feature_id]]            
        return [(r['parent_id'],r['value']) for r in self.hdf_group.feature_assoc.where('child_id == feature_id')]
    def get_children(self, feature_id):
        tbl = self._get_interval_table()
        if (self.child_id_index is not None) and (feature_id in self.child_id_index):
            return [(tbl[child_id],value) for child_id,value in self.child_id_index[feature_id]]            
        return [(tbl[r['child_id']],r['value']) for r in self.hdf_group.feature_assoc.where('parent_id == feature_id')]
    def get_child_ids(self, feature_id):
        if (self.child_id_index is not None) and (feature_id in self.child_id_index):
            return [(child_id,value) for child_id,value in self.child_id_index[feature_id]]            
        return [(r['child_id'],r['value']) for r in self.hdf_group.feature_assoc.where('parent_id == feature_id')]

    def _add_alias(self, feature_id, alias_source, alias_name):
        # add alias
        tbl = self.hdf_group.aliases
        row = tbl.row
        alias_id = tbl.nrows
        row['id'] = alias_id
        row['name'] = alias_name
        row['source'] = alias_source
        row.append()
        tbl.flush()
        # add association
        tbl = self.hdf_group.alias_assoc
        row = tbl.row
        row['alias_id'] = alias_id
        row['feature_id'] = feature_id
        row.append()
        tbl.flush()

    def _add_assoc(self, parent_id, child_id, value):
        tbl = self.hdf_group.feature_assoc
        row = tbl.row
        row['parent_id'] = parent_id
        row['child_id'] = child_id
        row['value'] = value
        row.append()
        tbl.flush()

#    def import_alias_file(self, fileh, source):
#        for line in fileh:
#            print line
#            fields = line.strip().split('\t')
#            known_alias, new_alias = fields
#            feature_ids = self.get_alias_feature_ids(known_alias)
#            for id in feature_ids:
#                print id, source, known_alias, new_alias                
#                self._add_alias(id, source, new_alias)
    
    def fromfeatures(self, feature_iter, source=None): 
        num_lines = 0
        debug_every = 10000
        debug_count = debug_every
        rownum = 0
        interval_exon_map = {}
        tx_id_map = {}
        gene_id_map = {}
        if source is None:
            source = "default"
        
        tbl = self.hdf_group.interval_data
        for f in feature_iter:
            # check if reference is valid
            try:
                ref, start, end, strand = self._parse_interval((f.ref, f.start, f.end, f.strand))
            except TrackError as e:
                logging.warning("[FeatureTrack] %s" % (str(e)))
                continue
            # map unique intervals to exons
            interval = (ref, start, end, strand)
            if interval in interval_exon_map:
                exon_id = interval_exon_map[interval]
            else:
                # add exon
                exon_id = rownum
                rownum += 1
                row = tbl.row
                row[REF_COL_NAME] = ref
                row[START_COL_NAME] = start
                row[END_COL_NAME] = end
                row['id'] = exon_id
                row['name'] = 'EXON%s' % exon_id
                row['source'] = f.source
                row['feature_type'] = EXON_TYPE
                row['score'] = f.score
                row['strand'] = strand
                row.append()
                interval_exon_map[interval] = exon_id
            # map exons to transcripts
            if f.transcript_id in tx_id_map:
                tx_id = tx_id_map[f.transcript_id]
                tx_start = tbl[tx_id][START_COL_NAME]
                tx_end = tbl[tx_id][END_COL_NAME]
                # update transcript
                if start < tx_start:
                    tx_start = start
                if end > tx_end:
                    tx_end = end
                tbl[tx_id][START_COL_NAME] = tx_start
                tbl[tx_id][END_COL_NAME] = tx_end             
                tbl[tx_id]['score'] += f.score                
            else:
                # create transcript id
                tx_id = rownum
                rownum += 1
                tx_start = start
                tx_end = end
                # add feature
                row = tbl.row
                row[REF_COL_NAME] = ref
                row[START_COL_NAME] = start
                row[END_COL_NAME] = end
                row['id'] = tx_id
                row['name'] = f.transcript_id
                row['source'] = f.source
                row['feature_type'] = TRANSCRIPT_TYPE
                row['score'] = f.score
                row['strand'] = strand
                row.append()
                tx_id_map[f.transcript_id] = tx_id      
            # map transcripts to genes
            if f.gene_id in gene_id_map:
                gene_id = gene_id_map[f.gene_id]
                gene_start = tbl[gene_id][START_COL_NAME]
                gene_end = tbl[gene_id][END_COL_NAME]
                # update gene
                if start < gene_start:
                    tbl[gene_id][START_COL_NAME] = start
                if end > gene_end:
                    tbl[gene_id][END_COL_NAME] = end                    
                tbl[gene_id]['score'] += f.score                    
            else:
                # create gene id
                gene_id = rownum
                rownum += 1
                # add feature
                row = tbl.row
                row[REF_COL_NAME] = ref
                row[START_COL_NAME] = tx_start
                row[END_COL_NAME] = tx_end
                row['id'] = gene_id
                row['name'] = f.gene_id
                row['source'] = f.source
                row['feature_type'] = GENE_TYPE
                row['score'] = f.score
                row['strand'] = strand
                row.append()
                gene_id_map[f.gene_id] = gene_id
            # add associations
            self._add_assoc(tx_id, exon_id, f.exon_number)
            self._add_assoc(gene_id, tx_id, 0)
            # add gene aliases
            for alias_id,alias_value in f.aliases:
                self._add_alias(gene_id, alias_id, alias_value)
            # done, flush table
            tbl.flush()
            # debugging
            num_lines += 1
            if num_lines == debug_count:
                logging.debug("Read %d lines" % (num_lines))
                debug_count += debug_every
        tbl.flush()
