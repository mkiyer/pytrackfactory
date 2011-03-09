'''
Created on Mar 8, 2011

@author: mkiyer
'''
import numpy as np
import operator

# constant HDF attribute to store root id
ROOT_ID_ATTR = "root_id"

# the numpy dtype used to represent interval nodes
_interval_node_dtype = np.dtype([('id', '<u4'),
                                 ('cleft_id', '<u4'), 
                                 ('cright_id', '<u4'), 
                                 ('croot_id', '<u4'), 
                                 ('start', '<i4'), 
                                 ('end', '<i4'), 
                                 ('minstart', '<i4'), 
                                 ('maxend', '<i4'), 
                                 ('minend', '<i4'),
                                 ('priority', '<f8'),
                                 ('value', np.int)])

# set an index for EMPTY nodes
EMPTY = 0
# constant for prioritizing nodes
NLOG = -1.0 / np.log(0.5)

def imax2(a, b):
    if b > a: return b
    return a

def imax3(a, b, c):
    if b > a:
        if c > b:
            return c
        return b
    if a > c:
        return a
    return c

def imin3(a, b, c):
    if b < a:
        if c < b:
            return c
        return b
    if a < c:
        return a
    return c

def imin2(a, b):
    if b < a: return b
    return a

def _init_node(tbl, id, start, end, value):
    node = tbl[id]
    node['id'] = id
    node['cleft_id'] = EMPTY
    node['cright_id'] = EMPTY
    node['croot_id'] = EMPTY
    node['start'] = start
    node['end'] = end
    node['maxend'] = end
    node['minstart'] = start
    node['minend'] = end
    # set priority according to binomial distribution
    node['priority'] = np.ceil(NLOG * np.log(-1.0/(np.random.random() - 1)))
    node['value'] = value
    return id

def _set_ends(node, cleft, cright):
    if ((node['cright_id'] != EMPTY) and (node['cleft_id'] != EMPTY)):
        node['maxend'] = imax3(node['end'], cright['maxend'], cleft['maxend'])
        node['minend'] = imin3(node['end'], cright['minend'], cleft['minend'])
        node['minstart'] = imin3(node['start'], cright['minstart'], cleft['minstart'])
    elif (node['cright_id'] != EMPTY):
        node['maxend'] = imax2(node['end'], cright['maxend'])
        node['minend'] = imin2(node['end'], cright['minend'])
        node['minstart'] = imin2(node['start'], cright['minstart'])
    elif (node['cleft_id'] != EMPTY):
        node['maxend'] = imax2(node['end'], cleft['maxend'])
        node['minend'] = imin2(node['end'], cleft['minend'])
        node['minstart'] = imin2(node['start'], cleft['minstart'])


def _insert(tbl, croot_id, id, start, end, value):
    """
    Insert a new IntervalNode into the tree of which the 'node' specified is
    currently the root. The return value is the new root of the tree (which
    may or may not be this node!)
    """
    node = tbl[croot_id]
    # If starts are the same, decide which to add interval to based on
    # end, thus maintaining sortedness relative to start/end        
    decision_endpoint = start
    if start == node['start']:
        decision_endpoint = end
    
    if decision_endpoint > node['start']:
        # insert to cright tree
        if node['cright_id'] != EMPTY:
            cright_id = _insert(tbl, node['cright_id'], id, start, end, value)                
        else:
            cright_id = _init_node(tbl, id, start, end, value)
        # modify the 'cright_id' column        
        node['cright_id'] = cright_id
        # rebalance tree
        if node['priority'] < tbl[cright_id]['priority']:
            croot_id = _rotate_left(tbl, node['id'])
    else:
        # insert to cright tree
        if node['cleft_id'] != EMPTY:
            cleft_id = _insert(tbl, node['cleft_id'], id, start, end, value)                
        else:
            cleft_id = _init_node(tbl, id, start, end, value)     
        # modify the 'cleft_id' column
        node['cleft_id'] = cleft_id
        # rebalance tree
        if node['priority'] < tbl[cleft_id]['priority']:
            croot_id = _rotate_right(tbl, node['id'])
    # set new croot
    croot = tbl[croot_id]
    _set_ends(croot, tbl[croot['cleft_id']], tbl[croot['cright_id']])
    tbl[node['cleft_id']]['croot_id'] = croot_id
    tbl[node['cright_id']]['croot_id'] = croot_id
    return croot_id

def _rotate_right(tbl, id):
    node = tbl[id]
    croot = tbl[node['cleft_id']]
    node['cleft_id'] = croot['cright_id']
    croot['cright_id'] = node['id']
    _set_ends(node, tbl[node['cleft_id']], tbl[node['cright_id']])
    return croot['id']

def _rotate_left(tbl, id):
    node = tbl[id]
    croot = tbl[node['cright_id']]
    node['cright_id'] = croot['cleft_id']
    croot['cleft_id'] = node['id']
    _set_ends(node, tbl[node['cleft_id']], tbl[node['cright_id']])
    return croot['id']

def _intersect(tbl, id, start, end, results):
    node = tbl[id]
    # Left subtree
    if ((node['cleft_id'] != EMPTY) and 
        (tbl[node['cleft_id']]['maxend'] > start)):
        _intersect(tbl, node['cleft_id'], start, end, results)
    # This interval
    if ((node['end'] > start) and (node['start'] < end)):
        results.append(node['value'])
    # Right subtree
    if ((node['cright_id'] != EMPTY) and 
        (node['start'] < end)):
        _intersect(tbl, node['cright_id'], start, end, results)

def _traverse(tbl, id, func):
    node = tbl[id]
    if node['cleft_id'] != EMPTY:
        _traverse(tbl, tbl[node['cleft_id']], func)
    func(node['value'])
    if node['cright_id'] != EMPTY:
        _traverse(tbl, tbl[node['cright_id']], func)

def _seek_left(tbl, id, position, results, n, max_dist):   
    node = tbl[id]         
    # we know we can bail in these 2 cases.
    if position >= node['maxend'] + max_dist: return
    if position < node['minstart']: return

    # the ordering of these 3 blocks makes it so the results are
    # ordered nearest to furthest from the query position
    if node['cright_id'] != EMPTY:
        _seek_left(tbl, node['cright_id'], position, results, n, max_dist)

    if node['end'] <= position < node['end'] + max_dist:
        results.append((node['end'], node['value']))

    if node['cleft_id'] != EMPTY:
        _seek_left(tbl, node['cleft_id'], position, results, n, max_dist)

def _seek_right(tbl, id, position, results, n, max_dist):
    node = tbl[id]
    # we know we can bail in these 2 cases.
    if position >= node['maxend']: return        
    if position + max_dist < node['minstart']: return
    # the ordering of these 3 blocks makes it so the results are
    # ordered nearest to furthest from the query position
    if node['cleft_id'] != EMPTY:
        _seek_right(tbl, node['cleft_id'], position, results, n, max_dist)

    if position < node['start'] <= position + max_dist:        
        results.append((node['start'], node['value']))

    if node['cright_id'] != EMPTY:
        _seek_right(tbl, node['cright_id'], position, results, n, max_dist)

class IntervalTree(object):

    def __init__(self, max_intervals):
        self.tbl = np.zeros(max_intervals+1, dtype=_interval_node_dtype)
        self.root_id = _init_node(self.tbl, EMPTY, 0, 0, 0)

    def insert(self, cur_id, start, end, value):
        if self.root_id == EMPTY:
            # create the first node and set to root
            self.root_id = _init_node(self.tbl, cur_id, start, end, value)
        else:
            self.root_id = _insert(self.tbl, self.root_id, cur_id, start, 
                                   end, value)

    def intersect(self, start, end):
        """
        given an id, start, and end, return a list of features
        falling within that range
        """
        results = []
        _intersect(self.tbl, self.root_id, start, end, results)
        return results

    def left(self, position, n=1, max_dist=2500):       
        """
        find n features with a start > than `position`
        f: a Interval object (or anything with an `end` attribute)
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        results = []
        _seek_left(self.tbl, self.root_id, position, results, n, max_dist)
        results.sort(key=operator.itemgetter(0), reverse=True)
        return [x[1] for x in results[:n]]
    
    def right(self, position, n=1, max_dist=2500):        
        """
        find n features with a end < than position
        f: a Interval object (or anything with a `start` attribute)
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        results = []
        _seek_right(self.tbl, self.root_id, position, results, n, max_dist)
        results.sort(key=operator.itemgetter(0))
        return [x[1] for x in results[:n]]

    def traverse(self, func):
        _traverse(self.tbl, self.root_id, func)

    def tohdf(self, parentgroup, name):
        h5file = parentgroup._v_file
        if name in parentgroup:
            h5file.removeNode(parentgroup, name, recursive=False)
        tbl = h5file.createTable(parentgroup, name, self.tbl)
        tbl.attrs[ROOT_ID_ATTR] = self.root_id

    @staticmethod    
    def fromhdf(hdf_tbl):    
        tree = IntervalTree(hdf_tbl.nrows-1)        
        tree.tbl[:] = hdf_tbl.read()
        tree.root_id = hdf_tbl.attrs[ROOT_ID_ATTR]
        return tree

