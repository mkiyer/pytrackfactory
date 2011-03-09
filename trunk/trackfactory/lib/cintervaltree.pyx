'''
Created on Sep 15, 2010

@author: mkiyer
'''
import operator
import numpy as np

cimport numpy as np
cimport cython

#np.import_array()

cdef extern from "stdlib.h":
    int ceil(float f)
    float log(float f)
    int RAND_MAX
    int rand()
    int strlen(char *)
    int iabs(int)

cdef float nlog = -1.0 / log(0.5)

cdef inline int imax2(int a, int b):
    if b > a: return b
    return a

cdef inline int imax3(int a, int b, int c):
    if b > a:
        if c > b:
            return c
        return b
    if a > c:
        return a
    return c

cdef inline int imin3(int a, int b, int c):
    if b < a:
        if c < b:
            return c
        return b
    if a < c:
        return a
    return c

cdef inline int imin2(int a, int b):
    if b < a: return b
    return a

# the C structure used to hold interval node information
cdef packed struct IntervalNode:
    np.uint32_t id
    np.uint32_t cleft_id
    np.uint32_t cright_id
    np.uint32_t croot_id
    np.int32_t start
    np.int32_t end
    np.int32_t minstart
    np.int32_t maxend
    np.int32_t minend
    np.float64_t priority
    np.uint32_t value

# the numpy dtype used to represent interval nodes
interval_node_dtype = np.dtype([('id', '<u4'),
                                ('cleft_id', '<u4'), 
                                ('cright_id', '<u4'), 
                                ('croot_id', '<u4'), 
                                ('start', '<i4'), 
                                ('end', '<i4'), 
                                ('minstart', '<i4'), 
                                ('maxend', '<i4'), 
                                ('minend', '<i4'),
                                ('priority', '<f8'),
                                ('value', '<u4')])

# set an index for EMPTY nodes
cdef int EMPTY = 0

cdef inline void _init_pnode(IntervalNode* node, int id, int start, int end, int value):
    node.id = id
    node.cleft_id = EMPTY
    node.cright_id = EMPTY
    node.croot_id = EMPTY
    node.start = start
    node.end = end
    node.maxend = end
    node.minstart = start
    node.minend = end
    # set priority according to binomial distribution
    node.priority = ceil(nlog * log(-1.0/(1.0 * rand()/RAND_MAX - 1)))
    node.value = value
    
def _init_node(np.ndarray[IntervalNode, ndim=1] tbl, int id, int start, int end, int value):
    _init_pnode(&tbl[id], id, start, end, value)
    return id

cdef inline void _set_ends(IntervalNode* p_node, IntervalNode* cleft, IntervalNode* cright):
    if ((p_node.cright_id != EMPTY) and (p_node.cleft_id != EMPTY)):
        p_node.maxend = imax3(p_node.end, cright.maxend, cleft.maxend)
        p_node.minend = imin3(p_node.end, cright.minend, cleft.minend)
        p_node.minstart = imin3(p_node.start, cright.minstart, cleft.minstart)
    elif (p_node.cright_id != EMPTY):
        p_node.maxend = imax2(p_node.end, cright.maxend)
        p_node.minend = imin2(p_node.end, cright.minend)
        p_node.minstart = imin2(p_node.start, cright.minstart)
    elif (p_node.cleft_id != EMPTY):
        p_node.maxend = imax2(p_node.end, cleft.maxend)
        p_node.minend = imin2(p_node.end, cleft.minend)
        p_node.minstart = imin2(p_node.start, cleft.minstart)

def _insert(np.ndarray[IntervalNode, ndim=1] tbl, int croot_id, int id, int start, int end, int value):
    """
    Insert a new IntervalNode into the tree of which the 'node' specified is
    currently the root. The return value is the new root of the tree (which
    may or may not be this node!)
    """
    cdef IntervalNode* p_node = &(tbl[croot_id])
    cdef int cright_id, cleft_id

    # If starts are the same, decide which to add interval to based on
    # end, thus maintaining sortedness relative to start/end        
    cdef int decision_endpoint = start
    if start == p_node.start:
        decision_endpoint = end
    
    if decision_endpoint > p_node.start:
        # insert to cright tree
        if p_node.cright_id != EMPTY:
            cright_id = _insert(tbl, p_node.cright_id, id, start, end, value)                
        else:
            cright_id = _init_node(tbl, id, start, end, value)
        # modify the 'cright_id' column        
        p_node.cright_id = cright_id
        # rebalance tree
        if p_node.priority < tbl[cright_id].priority:
            croot_id = _rotate_left(tbl, p_node.id)
    else:
        # insert to cright tree
        if p_node.cleft_id != EMPTY:
            cleft_id = _insert(tbl, p_node.cleft_id, id, start, end, value)                
        else:
            cleft_id = _init_node(tbl, id, start, end, value)     
        # modify the 'cleft_id' column
        p_node.cleft_id = cleft_id
        # rebalance tree
        if p_node.priority < tbl[cleft_id].priority:
            croot_id = _rotate_right(tbl, p_node.id)
    # set new croot
    cdef IntervalNode *croot = &tbl[croot_id]
    _set_ends(croot, &tbl[croot.cleft_id], &tbl[croot.cright_id])
    tbl[p_node.cleft_id].croot_id = croot_id
    tbl[p_node.cright_id].croot_id = croot_id
    return croot_id

def _rotate_right(np.ndarray[IntervalNode, ndim=1] tbl, int id):
    cdef IntervalNode* node = &tbl[id]
    cdef IntervalNode* croot = &tbl[node.cleft_id]
    node.cleft_id = croot.cright_id
    croot.cright_id = node.id
    _set_ends(node, &tbl[node.cleft_id], &tbl[node.cright_id])
    return croot.id

def _rotate_left(np.ndarray[IntervalNode, ndim=1] tbl, int id):
    cdef IntervalNode* node = &tbl[id]
    cdef IntervalNode* croot = &tbl[node.cright_id]
    node.cright_id = croot.cleft_id
    croot.cleft_id = node.id
    _set_ends(node, &tbl[node.cleft_id], &tbl[node.cright_id])
    return croot.id

def insert(np.ndarray[IntervalNode, ndim=1] tbl, int root_id, int cur_id, int start, int end, int value):
    if root_id == EMPTY:
        # create the first node and set to root
        root_id = _init_node(tbl, cur_id, start, end, value)
    else:
        root_id = _insert(tbl, root_id, cur_id, start, end, value)
    return root_id

def init_tree(num_rows):    
    cdef np.ndarray[IntervalNode, ndim=1] tbl = np.zeros(num_rows+1, dtype=interval_node_dtype)
    #tbl = np.zeros(num_rows+1, dtype=interval_node_dtype)
    return tbl, EMPTY
    _init_node(tbl, EMPTY, 0, 0, 0)
    return tbl, EMPTY

def _intersect(np.ndarray[IntervalNode, ndim=1] tbl, 
               int id, int start, int end, list results):
    cdef IntervalNode* node = &tbl[id]
    # Left subtree
    if ((node.cleft_id != EMPTY) and 
        (tbl[node.cleft_id].maxend > start)):
        _intersect(tbl, node.cleft_id, start, end, results)
    # This interval
    if ((node.end > start) and (node.start < end)):
        results.append(node.value)
    # Right subtree
    if ((node.cright_id != EMPTY) and 
        (node.start < end)):
        _intersect(tbl, node.cright_id, start, end, results)

def intersect(np.ndarray[IntervalNode, ndim=1] tbl, int id, int start, int end):
    """
    given an id, start, and end, return a list of features
    falling within that range
    """
    cdef list results = []
    _intersect(tbl, id, start, end, results)
    return results

def traverse(np.ndarray[IntervalNode, ndim=1] tbl, int id, object func):
    cdef IntervalNode* node = &tbl[id]
    if node.cleft_id != EMPTY:
        traverse(tbl, tbl[node.cleft_id], func)
    func(node.value)
    if node.cright_id != EMPTY:
        traverse(tbl, tbl[node.cright_id], func)

def _seek_left(np.ndarray[IntervalNode, ndim=1] tbl, int id, int position, 
               list results, int n, int max_dist):   
    cdef IntervalNode* node = &tbl[id]         
    # we know we can bail in these 2 cases.
    if position >= node.maxend + max_dist: return
    if position < node.minstart: return

    # the ordering of these 3 blocks makes it so the results are
    # ordered nearest to furthest from the query position
    if node.cright_id != EMPTY:
        _seek_left(tbl, node.cright_id, position, results, n, max_dist)

    if node.end <= position < node.end + max_dist:
        results.append((node.end, node.value))

    if node.cleft_id != EMPTY:
        _seek_left(tbl, node.cleft_id, position, results, n, max_dist)

def _seek_right(np.ndarray[IntervalNode, ndim=1] tbl, int id, int position, 
                list results, int n, int max_dist):
    cdef IntervalNode* node = &tbl[id]
    # we know we can bail in these 2 cases.
    if position >= node.maxend: return        
    if position + max_dist < node.minstart: return

    # the ordering of these 3 blocks makes it so the results are
    # ordered nearest to furthest from the query position
    if node.cleft_id != EMPTY:
        _seek_right(tbl, node.cleft_id, position, results, n, max_dist)

    if position < node.start <= position + max_dist:        
        results.append((node.start, node.value))

    if node.cright_id != EMPTY:
        _seek_right(tbl, node.cright_id, position, results, n, max_dist)

def left(np.ndarray[IntervalNode, ndim=1] tbl, int id, int position, 
         int n=1, int max_dist=2500):       
    """
    find n features with a start > than `position`
    f: a Interval object (or anything with an `end` attribute)
    n: the number of features to return
    max_dist: the maximum distance to look before giving up.
    """
    cdef list results = []
    _seek_left(tbl, id, position, results, n, max_dist)
    #if len(results) == n: return results
    results.sort(key=operator.itemgetter(0), reverse=True)
    return [x[1] for x in results[:n]]

def right(np.ndarray[IntervalNode, ndim=1] tbl, int id, int position, 
          int n=1, int max_dist=2500):        
    """
    find n features with a end < than position
    f: a Interval object (or anything with a `start` attribute)
    n: the number of features to return
    max_dist: the maximum distance to look before giving up.
    """
    cdef list results = []
    # use end + 1 becuase .right() assumes strictly right-of
    _seek_right(tbl, id, position, results, n, max_dist)
    #if len(results) == n: return results
    results.sort(key=operator.itemgetter(0))
    return [x[1] for x in results[:n]]
