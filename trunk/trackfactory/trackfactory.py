'''
.. module:: trackfactory
   :synopsis: defines TrackFactory class
.. moduleauthor:: Matthew Iyer <mkiyer@umich.edu>
.. versionadded:: 0.0.1
'''
import os
import numpy as np
import tables

from track import TrackError, REF_TRACK_NAME, TRACK_CLASS_ATTR, get_ref_length
from arraytrack import ArrayTrack
from coveragetrack import CoverageTrack
from bitarraytrack import BitArrayTrack
from sequencetrack import SequenceTrack
from intervaltrack import IntervalTrack
from rnaseqtrack import RnaseqTrack
#from sequel2.track.featuretrack import FeatureTrack
#from sequel2.track.expressiontrack import ExpressionTrack

def _open_hdf_file(filename, mode='r'):
    '''opens an hdf file
    
    mode has several options
    ------------------------
    'r'  Read-only; no data can be modified.
    'w'  Write; a new file is created (an existing file with the same name would be deleted).
    'a'  Append; an existing file is opened for reading and writing, and if the file does not exist it is created.
    'r+' It is similar to 'a', but the file must already exist.
    
    :returns: :class:`tables.File` object
    '''
    if (mode != 'w') and os.path.exists(filename):
        if not tables.isHDF5File(filename):
            raise IOError('file %s is not an HDF5 file' % filename)
        if not tables.isPyTablesFile(filename):
            raise IOError('file %s is not a PyTables file' % filename)
    return tables.openFile(filename, mode)

def _close_hdf_file(h5file):
    '''closes an hdf5 file
    '''
    h5file.close()

class TrackFactory(object):
    '''
    Open/Create a TrackFactory object
    
    filename: trackfactory file
    mode: read/write mode for file
    refs (create only): list (reference,length) tuples
    
    modes
    'r'  Read-only; no data can be modified.
    'w'  Write; a new file is created (an existing file with the same name would be deleted).
    'a'  Append; an existing file is opened for reading and writing, and if the file does not exist it is created.
    'r+' It is similar to 'a', but the file must already exist.
    '''    
    # constants
    track_classes = [ArrayTrack, 
                     CoverageTrack,
                     BitArrayTrack,
                     SequenceTrack, 
                     IntervalTrack,
                     RnaseqTrack]
    track_name_class_dict = dict((cls.__name__, cls) for cls in track_classes)

    def __init__(self, filename, mode="r", refs=None):
        self.h5file = _open_hdf_file(filename, mode)
        if ("/" + REF_TRACK_NAME) not in self.h5file:
            if refs is None:
                raise TrackError("TrackFactory: Must specify (ref,length) tuples in 'refs' argument")
            if mode == "r":
                raise TrackError("TrackFactory: Cannot create track factory when mode='r' (read only)")
            self._create_reference_track(refs)

    def _get_hdf_file(self):
        return self.h5file._v_file

    def create_track(self, track_name, track_class, *args, **kwargs):
        '''Add new :class:`Track` to the :class:`TrackFactory`
        
        :param track_name: track name
        :type track_name: str
        :param track_class: type of track to create
        :type track_class: :class:`ArrayTrack`, :class:`IntervalTrack`, etc
        :param args: args passed to Track creation method
        :param kwargs: kwargs passed to Track creation method        
        :returns: :class:`Track` object
        :raises: :class:`TrackError` if :class:`Track` already exists 
        '''
        # create subgroup
        if ('/' + track_name) in self.h5file:
            raise TrackError('track %s already exists' % track_name)
        group = self.h5file.createGroup(self.h5file.root, track_name)        
        # the subgroup needs to remember the track type
        setattr(group._v_attrs, TRACK_CLASS_ATTR, track_class.__name__)
        return track_class(group, *args, **kwargs)

    def get_track(self, track_name):
        """Get existing :class:Track
        
        :param track_name: track name
        :type track_name: str
        :returns:  :class:`Track` object
        :raises: :class:`TrackError` if 'track_name' not found 
        """
        if ('/' + track_name) not in self.h5file:
            raise TrackError('track %s does not exist' % track_name)        
        group = self.h5file.getNode(self.h5file.root, track_name)        
        track_class_name = getattr(group._v_attrs, TRACK_CLASS_ATTR)
        return self.track_name_class_dict[track_class_name](group)

    def has_track(self, name):
        '''returns `True` if `name` is a valid :class:`Track` in the 
        :class:`TrackFactory`, and `False` otherwise
        '''        
        return ('/' + name) in self.h5file
    
    def delete_track(self, track_name):
        """Delete a track from the TrackFactory
        
        :param track_name: track name to delete
        :type track_name: str
        :raises: TrackError if 'track_name' not found

        .. note: this is a permanent operation
        """        
        if not ('/' + track_name) in self.h5file:
            raise TrackError('error trying to delete track %s that does not exist' % track_name)
        self.h5file.removeNode(self.h5file.root, name=track_name, recursive=True)

    def close(self):
        self.h5file.close()
        self.h5file = None

    # reference track code
    ref_track_dtype = np.dtype([('name', 'a32'),
                                ('length', np.uint64)])    
    def _create_reference_track(self, refs):
        tbl = self.h5file.createTable(self.h5file.root, 
                                      REF_TRACK_NAME,
                                      self.ref_track_dtype,
                                      expectedrows=len(refs))
        for ref,length in refs:
            row = tbl.row
            row['name'] = ref
            row['length'] = length
            row.append()
        tbl.flush()

    def get_ref_length(self, ref):
        """return the size of the reference `ref`"""
        tbl = self.h5file.root._f_getChild(REF_TRACK_NAME)
        return get_ref_length(tbl, ref)

    def get_rnames(self):
        """generator function that returns list of references names from the 
        :mod:`TrackFactory` used to create this track
        
        :returns: iterator of references names (ordered by reference id)
        """
        for row in self.h5file.root._f_getChild(REF_TRACK_NAME):
            yield row['name']

    def get_lengths(self):
        """generator function that returns list of reference lengths from the 
        :mod:`TrackFactory` used to create this track
        
        :returns: array of lengths (ordered by reference id)
        """
        for row in self.h5file.root._f_getChild(REF_TRACK_NAME):
            yield row['length']

    def get_refs(self):
        """generator function that returns list of (reference,length) tuples
        from the :mod:`TrackFactory` used to create this track
        
        :returns: array of (reference,length) tuples (ordered by ref id)
        """        
        for row in self.h5file.root._f_getChild(REF_TRACK_NAME):
            yield row.fetch_all_fields()
