'''
Created on Sep 28, 2010

@author: mkiyer
'''
import collections
import tables
import numpy as np
import logging

from sequel2.track.intervaltrack import IntervalTrack, REF_COL_NAME, START_COL_NAME, END_COL_NAME
from sequel2.io.bed import parse_bed12_file

FEATURE_TABLE = "features"
FEATURE_ALIAS_TABLE = "aliases"
FEATURE_ASSOC_TABLE = "feature_assoc"
ALIAS_ASSOC_TABLE = "alias_assoc"

EXON_TYPE = "exon"
GENE_TYPE = "gene"
TRANSCRIPT_TYPE = "transcript"

# the numpy dtype used to represent interval nodes
feature_dtype = np.dtype([('id', '<u4'),
                          (REF_COL_NAME, 'a25'), 
                          ('source', 'a25'), 
                          ('feature_type', 'a25'), 
                          (START_COL_NAME, '<i4'), 
                          (END_COL_NAME, '<i4'), 
                          ('score', '<i4'),
                          ('strand', 'a1'),
                          ('phase', '<i4')])
feature_assoc_dtype = np.dtype([('parent_id', '<u4'),
                                ('child_id', '<u4'),
                                ('value', '<u4')])
feature_alias_dtype = np.dtype([('id', '<u4'),
                                ('source', 'a16'),
                                ('name', 'a128')])
alias_assoc_dtype = np.dtype([('alias_id', '<u4'),
                              ('feature_id', '<u4')])

class FeatureTrack(IntervalTrack):
    _expectedrows = 1e5

    def __init__(self, hdf_group):
        super(FeatureTrack, self).__init__(hdf_group, descriptor=feature_dtype)
        self._init_tables()
        self.alias_to_feature_index = None
        self.feature_to_alias_index = None
        self.parent_id_index = None
        self.child_id_index = None

    def _init_tables(self):
        h5file = self._get_hdf_file()
        if not FEATURE_ALIAS_TABLE in self.hdf_group:
            h5file.createTable(self.hdf_group, FEATURE_ALIAS_TABLE, 
                               feature_alias_dtype, 
                               expectedrows=self._expectedrows) 
        if not FEATURE_ASSOC_TABLE in self.hdf_group:
            h5file.createTable(self.hdf_group, FEATURE_ASSOC_TABLE,
                               feature_assoc_dtype,
                               expectedrows=self._expectedrows)
        if not ALIAS_ASSOC_TABLE in self.hdf_group:
            h5file.createTable(self.hdf_group, ALIAS_ASSOC_TABLE,
                               alias_assoc_dtype,
                               expectedrows=self._expectedrows)
    
    def index_aliases(self):
        # automatically index the aliases for fast lookup       
        alias_to_feature_index = collections.defaultdict(lambda: set())
        feature_to_alias_index = collections.defaultdict(lambda: set())
        for row in self.hdf_group.alias_assoc:
            alias_to_feature_index[row['alias_id']].add(row['feature_id'])
            feature_to_alias_index[row['feature_id']].add(row['alias_id'])
        self.alias_to_feature_index = alias_to_feature_index
        self.feature_to_alias_index = feature_to_alias_index
        return alias_to_feature_index, feature_to_alias_index

    def index_assoc(self):
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
        return self.hdf_group.interval_data[feature_id]
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
        if (self.parent_id_index is not None) and (feature_id in self.parent_id_index):
            return [(self.hdf_group.interval_data[parent_id],value) for parent_id,value in self.parent_id_index[feature_id]]            
        return [(self.hdf_group.interval_data[r['parent_id']],r['value']) for r in self.hdf_group.feature_assoc.where('child_id == feature_id')]
    def get_parent_ids(self, feature_id):
        if (self.parent_id_index is not None) and (feature_id in self.parent_id_index):
            return [(parent_id,value) for parent_id,value in self.parent_id_index[feature_id]]            
        return [(r['parent_id'],r['value']) for r in self.hdf_group.feature_assoc.where('child_id == feature_id')]
    def get_children(self, feature_id):
        if (self.child_id_index is not None) and (feature_id in self.child_id_index):
            return [(self.hdf_group.interval_data[child_id],value) for child_id,value in self.child_id_index[feature_id]]            
        return [(self.hdf_group.interval_data[r['child_id']],r['value']) for r in self.hdf_group.feature_assoc.where('parent_id == feature_id')]
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
        
    def import_alias_file(self, fileh, source):
        for line in fileh:
            print line
            fields = line.strip().split('\t')
            known_alias, new_alias = fields
            feature_ids = self.get_alias_feature_ids(known_alias)
            for id in feature_ids:
                print id, source, known_alias, new_alias                
                self._add_alias(id, source, new_alias)

    def import_bed_file(self, fileh, source):
        num_lines = 0
        debug_every = 100
        debug_count = debug_every
        rownum = 0
        exon_id_map = {}

        tbl = self.hdf_group.interval_data
        references = set(self._get_references())
        for g in parse_bed12_file(fileh):            
            # debugging
            num_lines += 1
            if num_lines == debug_count:
                logging.debug("Read %d lines" % (num_lines))
                debug_count += debug_every
            # check if reference is valid
            if g.chrom not in references:
                logging.warning("Reference %s not found in table" % g.chrom)
                continue
            # add feature
            row = tbl.row
            tx_id = rownum
            rownum += 1            
            row['id'] = tx_id
            row[REF_COL_NAME] = g.chrom
            row['source'] = source
            row['feature_type'] = TRANSCRIPT_TYPE
            row[START_COL_NAME] = g.tx_start
            row[END_COL_NAME] = g.tx_end
            row['score'] = g.score
            row['strand'] = g.strand
            row['phase'] = 0
            row.append()
            # add alias
            self._add_alias(tx_id, source, g.name)
            # add transcript exons 
            for exon_num,exon in enumerate(g.exons):
                start, end = exon
                interval = (g.chrom, start, end, g.strand)
                if interval in exon_id_map:
                    exon_id = exon_id_map[interval]
                else:
                    exon_id = rownum
                    rownum += 1
                    row['id'] = exon_id
                    row[REF_COL_NAME] = g.chrom
                    row['source'] = source
                    row['feature_type'] = EXON_TYPE
                    row[START_COL_NAME] = start
                    row[END_COL_NAME] = end
                    row['score'] = g.score
                    row['strand'] = g.strand
                    row['phase'] = 0
                    row.append()
                exon_id_map[interval] = exon_id
                # add association between transcript and exon
                self._add_assoc(tx_id, exon_id, exon_num)
        tbl.flush()
