'''
Created on Oct 20, 2010

@author: mkiyer
'''
import numpy as np
import logging
import multiprocessing
import tables
import tempfile
import os
import operator
import collections

from track import Track, TrackError
from trackfactory import TrackFactory
from featuretrack import FeatureTrack, EXON_TYPE, TRANSCRIPT_TYPE, GENE_TYPE

def make_temp(prefix, suffix, dir):
    fd, tmpfilename = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=dir)
    os.close(fd)
    return tmpfilename

FEATURE_TRACK_LINK = "feature_track_link"
SEQDATA_TABLE = "seqdata"
PROPS_TABLE = "props"
EXPRMAT_ARRAY = "exprmat"

seqdata_dtype = np.dtype([('index', '<u4'),
                          ('group', 'a64'),
                          ('experiment', 'a64'),
                          ('library', 'a64'),
                          ('seqdata', 'a64'),
                          ('num_sequences', '<u8'),
                          ('num_mapped_reads', '<u8'),
                          ('num_duplicate_reads', '<u8')])

props_dtype = np.dtype([('index', '<u4'), 
                        ('name', 'a64'), 
                        ('value', 'a256')])

def fetch_counts(vector_track_file, vector_track_name, feature_track_file, feature_track_name):
    logging.debug("Vector track file: %s" % vector_track_file)
    logging.debug("Vector track name: %s" % vector_track_name)
    logging.debug("Feature track file: %s" % feature_track_file)
    logging.debug("Feature track name: %s" % feature_track_name)
    # open files and get tracks
    vector_tf = TrackFactory(vector_track_file, "r")
    vectortrack = vector_tf.get_track(vector_track_name)
    feature_tf = TrackFactory(feature_track_file, "r")
    featuretrack = feature_tf.get_track(feature_track_name)
    # allocate array to store temporary expression    
    arr = np.zeros(featuretrack.num_intervals, dtype=np.float)
    # index parent/child associations for fast lookup
    logging.debug("%s: indexing features" % (vector_track_file))
    parent_id_dict, child_id_dict = featuretrack.index_assoc()
    logging.debug("%s: counting regions" % (vector_track_file))
    debug_count = 0
    debug_every = 10000
    debug_next = debug_every
    for feature in featuretrack:
        if feature['feature_type'] == EXON_TYPE:
            interval = (feature['ref'], feature['start'], feature['end'])
            count = vectortrack.count(interval)
            arr[feature['id']] = count
            #logging.debug("Exon id=%d count=%f" % (feature['id'], count))
            for parent_id,exon_num in parent_id_dict[feature['id']]:
                arr[parent_id] += count
                #logging.debug("Transcript id=%d count=%f" % (parent_id, arr[parent_id]))
        debug_count += 1
        if debug_count == debug_next:
            logging.debug("%s: Finished counting %d features" % (vector_track_file, debug_count))
            debug_next += debug_every
    feature_tf.close()
    vector_tf.close()
    return arr

def fetch_counts_worker(arglist):
    index = arglist[0]
    arr = fetch_counts(*arglist[1:])
    return index, arr

class ExpressionTrack(Track):            
    _expected_seqdata_rows = 1000
    _expr_data_dtype = np.float

    default_complevel = 1
    default_complib = 'lzo'
    default_filters = tables.Filters(complevel=default_complevel, 
                                     complib=default_complib,
                                     shuffle=True)
    
    def __init__(self, hdf_group, feature_track=None):
        super(ExpressionTrack, self).__init__(hdf_group)
        h5file = self._get_hdf_file()
        # link to feature track        
        if not FEATURE_TRACK_LINK in self.hdf_group:
            if feature_track is None:
                raise TrackError("'feature_track' param invalid, got '%s'" % (feature_track))
            h5file.createSoftLink(self.hdf_group, FEATURE_TRACK_LINK, feature_track.hdf_group)
        # create table for sequence data
        if not SEQDATA_TABLE in self.hdf_group:
            h5file.createTable(self.hdf_group, SEQDATA_TABLE, seqdata_dtype,
                               expectedrows=self._expected_seqdata_rows)
        # create table for sequence data properties
        if not PROPS_TABLE in self.hdf_group:
            h5file.createTable(self.hdf_group, PROPS_TABLE, props_dtype,
                               expectedrows=self._expected_seqdata_rows)
        # create expression matrix array
        if not EXPRMAT_ARRAY in self.hdf_group:
            if feature_track is None:
                raise TrackError("'feature_track' param invalid, got '%s'" % (feature_track))
            shape = (feature_track.num_intervals,0)            
            atom = tables.Atom.from_dtype(np.dtype(self._expr_data_dtype))
            h5file.createEArray(self.hdf_group, EXPRMAT_ARRAY,
                                atom=atom,
                                shape=shape,
                                filters=self.default_filters,
                                chunkshape=None)

    def _get_feature_track(self):
        return FeatureTrack(self.hdf_group.feature_track_link())

    def add_vector_tracks(self, seqdata_line_iter, processes=1):
        # copy the feature track file so it can be opened by 
        # multiple processes in read-only mode
        h5file = self._get_hdf_file()
        ftrack = self._get_feature_track()
        ftrack_name = ftrack.get_name()
        logging.debug("Creating temporary copy of file '%s'" % (h5file.filename))
        # open temp file and copy feature track
        ftrack_file = make_temp(suffix=".h5", prefix="ftrack", dir=os.path.dirname(h5file.filename))        
        f2 = tables.openFile(ftrack_file, 'w')
        new_ftrack = ftrack._f_copy(f2.root, ftrack_name, recursive=True)
        f2.close()
        # setup list of counting tasks
        exprmat = self.hdf_group.exprmat
        tasklist = []
        # make field index lists
        header_fields = seqdata_line_iter.next().strip().split('\t')
        seqdata_fields = set(seqdata_dtype.names).difference(set(['index']))       
        seqdata_field_indexes = [header_fields.index(f) for f in seqdata_fields]
        data_path_index = header_fields.index('track_path')         
        prop_field_indexes = set(range(len(header_fields)))
        prop_field_indexes.difference_update(seqdata_field_indexes + [data_path_index])        
        # setup individual tasks
        for line in seqdata_line_iter:
            fields = line.strip().split('\t')
            # get and validate data path
            data_path = fields[data_path_index]
            file_path, h5_path = data_path.split(":")
            if not os.path.isfile("file_path"):
                logging.error("file path '%s' not found" % (file_path))
                continue
            # add to metadata
            tbl = self.hdf_group.seqdata
            row = tbl.row
            index = tbl.nrows
            row['index'] = index
            for ind in seqdata_field_indexes:
                row[header_fields[ind]] = fields[ind]
            row.append()
            tbl.flush()        
            # add to properties
            tbl = self.hdf_group.props
            row = tbl.row
            for ind in prop_field_indexes:                
                row['index'] = index
                row['name'] = header_fields[ind]
                row['value'] = fields[ind]
                row.append()
            tbl.flush()        
            # add to earray
            exprmat.append(np.zeros((exprmat.shape[0],1), dtype=self._expr_data_dtype))
            # add to list of tasks
            logging.debug("Adding task index=%d data_path=%s track=%s feature_file=%s feature_track=%s" % 
                          (index, data_path, h5_path, ftrack_file, ftrack_name))
            tasklist.append((index, data_path, h5_path, ftrack_file, ftrack_name))
        logging.debug("Number of sequence data files: %d" % len(tasklist))
        logging.debug("Number of features: %d" % exprmat.shape[0])
        # kickoff jobs on multiple processors to count reads on each of the genes        
        pool = multiprocessing.Pool(processes=processes)
        result_iter = pool.imap_unordered(fetch_counts_worker, tasklist)
        for index,arr in result_iter:
            exprmat[:,index] = arr
            logging.debug("Finished task index=%d" % (index))
        pool.close()
        pool.join()
        # remove temp feature track 
        os.remove(ftrack_file)

    def _get_feature_length(self, feature):
        featuretrack = self._get_feature_track()
        if feature['feature_type'] == EXON_TYPE:
            return feature['end'] - feature['start']
        elif feature['feature_type'] == TRANSCRIPT_TYPE:
            length = 0
            for child in featuretrack.get_children(feature['id']):
                length += child['end'] - child['start']
            return length
        else:
            # TODO: need a recursive function for gene level
            return 1
    
    def group_matrix(self, group_by="library"):
        seqdata_tbl = self.hdf_group.seqdata
        # build lookup from grouped to original seqdata
        group_names = []
        group_metadata = []
        group_props = []
        group_to_seqdata_dict = collections.defaultdict(lambda: [])
        seqdata_to_group_dict = {}
        current_index = 0
        for row in seqdata_tbl:
            if row[group_by] not in group_names:
                group_index = current_index
                # add to group names
                group_names.append(row[group_by])                
                group_metadata.append(collections.defaultdict(lambda: set()))
                group_props.append(collections.defaultdict(lambda: set()))
                current_index += 1
            else:
                # get group index based on group name
                group_index = group_names.index(row[group_by])
            # add to group metadata
            for fieldname in seqdata_dtype.names:
                group_metadata[group_index][fieldname].add(row[fieldname])
            # add to group properties            
            for prop_row in self.hdf_group.props.where('index == rowindex', {'rowindex': row['index']}):
                group_props[group_index][prop_row['name']].add(prop_row['value'])
            # keep a mapping between new and old indexes
            group_to_seqdata_dict[group_index].append(row['index'])
            seqdata_to_group_dict[row['index']] = group_index
        # join columns
        numeric_columns = ['num_sequences', 'num_mapped_reads', 'num_duplicate_reads']
        for ind in xrange(len(group_names)):
            metadict = group_metadata[ind]
            for col in metadict.keys():
                if col in numeric_columns:
                    metadict[col] = sum(metadict[col])
                else:
                    metadict[col] = '||'.join(map(str, metadict[col]))
            propdict = group_props[ind]
            for col in propdict.keys():
                propdict[col] = '||'.join(map(str, propdict[col]))            
        # make new expression matrix
        num_groups = len(group_names)        
        seqdata_exprmat = self.hdf_group.exprmat
        grouped_exprmat = np.zeros((seqdata_exprmat.shape[0],num_groups), dtype=self._expr_data_dtype)
        # populate the matrix with summed columns
        for group_index,seqdata_indexes in group_to_seqdata_dict.iteritems():
            for seqdata_index in seqdata_indexes:
                grouped_exprmat[:,group_index] += seqdata_exprmat[:,seqdata_index]         
        return group_names, group_metadata, group_props, grouped_exprmat

    def write_tabular_text(self, outfh, index=True, rpkm=False, group_by="experiment"):
        featuretrack = self._get_feature_track()
        if index:
            featuretrack.index_assoc()
            featuretrack.index_aliases()
        # group seqdata according to grouping key
        group_names, group_metadata, group_props, exprmat = self.group_matrix(group_by)
        num_cols = len(group_names)
        num_seq_arr = np.array([m['num_sequences'] for m in group_metadata])
        num_seq_norm = 1.0e6 / num_seq_arr

        # write primary header and metadata
        header = ["id", "feature_type", "aliases", "genomic_coords", "feature_size", "exon_num"]
        print >>outfh, '\t'.join(header + group_names)
        metadata_keys = set()
        for m in group_metadata:
            metadata_keys.update(m.keys())
        for fieldname in metadata_keys:
            padding = [''] * (len(header) - 1)                        
            print >>outfh, '\t'.join(padding + [fieldname] + list(map(str, [m[fieldname] for m in group_metadata]))) 

        # write property metadata
        prop_keys = set()
        for p in group_props:
            prop_keys.update(p.keys())
        for prop_name in prop_keys:
            padding = [''] * (len(header) - 1)                        
            print >>outfh, '\t'.join(padding + [prop_name] + [p[prop_name] for p in group_props])

        def write_feature(f, name, length, exon_num=-1, rpkm=False):
            cols = [f['id'], 
                    f['feature_type'],
                    name,
                    "%s:%d-%d_%s" % (f['ref'], f['start'], f['end'], f['strand']),
                    length,
                    exon_num]
            exprs = exprmat[f['id'],:]
            if rpkm:
                exprs = (1e3 * exprs / length) * num_seq_norm
            cols.extend(["%.3f" % x for x in exprs])
            return cols

        # index features for fast lookup
        for feature in featuretrack:
            if feature['feature_type'] == TRANSCRIPT_TYPE:     
                alias_string = '|'.join(["%s=%s" % (a['source'],a['name']) for a in featuretrack.get_feature_aliases(feature['id'])])
                tx_length = 0
                # output children first
                lines = []
                for child,exon_num in featuretrack.get_children(feature['id']):
                    exon_length = child['end'] - child['start']
                    tx_length += exon_length 
                    lines.append(write_feature(child, alias_string, exon_length, exon_num, rpkm=rpkm))
                # output parent
                lines.append(write_feature(feature, alias_string, tx_length, len(lines), rpkm=rpkm))
                lines = sorted(lines, key=operator.itemgetter(5))
                for fields in lines:
                    print >>outfh, '\t'.join(map(str, fields))
