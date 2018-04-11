 
import numpy as np
import sys
import cPickle
import itertools

try:
    from scipy.spatial import cKDTree
    SCIPY_FLAG = True
except:
    SCIPY_FLAG = False

################################################################################
## GLOBAL VARIABLES
################################################################################
DELTA = 10**-7
def intersect1d(tup_of_arr):
    """ intersects multiple 1d arrays.
    Args:
        tup_of_arr - tuple of arrays that are intersected
    Returns:
        rslt_arr   - intersected array
    """
    if len(tup_of_arr)<2:
        print 'ERROR: NOT ENOUGH ARRAYS TO INTERSECT'
    if 0 in map(len,tup_of_arr):
        return np.array([],dtype=int)
    rslt_arr = tup_of_arr[0]
    for arr_i in tup_of_arr[1:]:
        rslt_arr = np.intersect1d( rslt_arr, arr_i)
    return rslt_arr
# READ AND WRITE INTO ABQ
class AbstractMesh(object):
    """    
    """
    def __init__(self,dim):        
        self.dim = dim
        self.nds = np.ndarray([0,self.dim])
        self.nds_sets= {}
        self.ele_sets= {}
    def add_ele(self):
        return
    def add_nds_ele(self,nds,ele,ele_set_nam=None):
        nb_ele = len(self.ele)
        if np.shape(self.ele)[1]!=np.shape(ele)[1]:
            print 'ERROR: ADDED ELEMENTS ARE NOT THE SAME'
        if ele_set_nam!=None and isinstance(ele_set_nam,str):
            if ele_set_nam in self.ele_sets.keys():
                self.ele_sets[ele_set_nam] = \
                    np.concatenate( (self.ele_sets[ele_set_nam],
                                     np.arange(nb_ele,nb_ele+len(ele))))
            else:
                self.ele_sets.update({ele_set_nam:np.arange(nb_ele,nb_ele+len(ele))})            
        self.ele = np.vstack( (self.ele,ele+len(self.nds)) )        
        self.nds = np.vstack( (self.nds,nds) )    
    def merge_nds(self,delta=10**-8):
        self.nds, idx_cor = merge_nds(self.nds,delta=delta)
        self.ele          = idx_cor[self.ele]
        for key_i, nd_set_i in self.nds_sets.iteritems():
            nd_set_i = idx_cor[nd_set_i]        
    def rm_extra_nds_ele(self):
        self.nds,self.ele = rm_extra_nodes(self.nds,self.ele)
        self.ele          = unique_rows(np.sort(self.ele,axis=1))        
    def order_ele_set(self,merge_incl_sets=False):
        if len(self.ele_sets)<=1:
            print 'WARNING: NO ELEMENT-SETS TO ORDER'
        else:
            for key in self.ele_sets.keys():
                self.ele_sets[key+'old'] = self.ele_sets.pop(key)
            set_nams  = self.ele_sets.keys()
            set_len_s = -1*np.array([ len(ele_set)
                                      for ele_set in self.ele_sets.values()])
            ids = set_len_s.argsort()
            if merge_incl_sets:
                self.ele_sets['MAT_0'] = self.ele_sets.pop(set_nams[ids[0]])
                self.ele_sets.update({'MAT_1':[]})
                for id_i in ids[1:]:
                    self.ele_sets['MAT_1']+= self.ele_sets.pop(set_nams[id_i])
            else:
                for i,id_i in enumerate(ids):
                    old_key = set_nams[id_i]
                    self.ele_sets['MAT'+str(i)] = self.ele_sets.pop(old_key) 
class Read_abq_file(AbstractMesh):
    '''
    Read the input file by going through it twice:
     - First time is to check how many nodes, elements and sets there are.
     - Second time is to read the nodes coordinates, element connectivities and
       element sets.
    The sets are stored in a long one-dimensional array, where each entry
    represents the element-set and it's index the element number.
    '''
    def __init__(self,f_name):
        AbstractMesh.__init__(self)
        if f_name.endswith('.inp'):
            f_name = f_name[:-4]
        with open(f_name+'.inp', "rU") as f:
            nb_nds = 0
            nb_ele = 0
            nb_set = 0
            self.nds      = []
            self.ele      = []
            self.ele_sets = {}
            while True:
                line = f.readline()
                if line=='':
                    break
                line = line.upper().strip()
                if line.startswith('**') or line.startswith('\n') or line=='':
                    continue
                if line.startswith('*'):
                    nds_flag = False
                    ele_flag = False
                    set_flag = False                
                    keywords = line.split(',')                
                    if   keywords[0].upper().startswith('*NODE') \
                         and not keywords[0].upper().startswith('*NODE OUTPUT'):
                        nds_flag = True
                        continue
                    elif keywords[0].upper().startswith('*ELEMENT'):
                        ele_flag = True
                        matching = [s for s in keywords if "ELSET" in s]
                        if len(matching)>0:
                            set_name      = matching[0].split('=')[-1]
                            self.ele_sets[set_name]=[]
                            nb_set       += 1
                            set_flag      = True
                        continue
                    elif keywords[0].upper().startswith('*ELSET'):
                        set_flag = True
                        matching = [s for s in keywords if "ELSET" in s]
                        if len(matching)>0:
                            set_name      = matching[-1].split('=')[-1]                      
                            self.ele_sets[set_name]=[]
                            nb_set       += 1
                        continue    
                if nds_flag:
                    self.nds.append( map(float,line.split(',')[1:]) )
                    nb_nds += 1
                elif ele_flag:
                    line = line.split(',')
                    self.ele.append( map(int,line[1:]) )
                    if set_flag:
                        self.ele_sets[set_name].append(int(line[0]))                
                elif set_flag:
                    self.ele_sets[set_name].extend(map(int,line.strip(',\n').split(',')))
            # ACCOUNT FOR ABAQUS STARTING COUNTING WITH 1
            self.nds = np.array(self.nds)
            self.ele = np.array(self.ele)-1
            for (set_nam,ele_set_i) in self.ele_sets.iteritems():
                self.ele_sets[set_nam] = np.array(self.ele_sets[set_nam])-1

class ABQ_inp(object):
    def __init__(self,name):
        self.dim      = 3
        self.name     = name
        self.parts    = {}
        self.assembly = None
        self.step     = []
        self.materials= {}
        self.equations= ''
        self.sets     = ''
class ABQ_assembly(object):    
    def __init__(self,name):
        self.name          = name
        self.instances     = {}
        self.nds           = []
        self.nds_sets      = {}
        self.srf_nds_based = {}
        self.srf_ele_based = {}
        self.interactions  = []
class ABQ_part(object):
    def __init__(self,name):
        self.name           = name
        self.nds            = []
        self.elements       = []
        self.ele_type_order = []
        self.nds_sets       = {}
        self.ele_sets       = {}
        self.srf_nds_based  = {}
        self.srf_ele_based  = {}
        self.sections       = []
class ABQ_instance(object):
    def __init__(self,part):        
        self.part = ''
        self.sets = {}
class ABQ_step(object):
    def __init__(self):
        self.name = ''
        self.bc   = ''
def read_abq_part(f,abq_inp,part_name):
    ''' read the part information
    '''
    elements       = {}
    ele_type_order = []
    ele_sets       = {}
    solid_sections = []
    while True:
        line = f.readline().upper().replace('\n','')
        if line.startswith('*END PART'):
            break
        if line.startswith('*NODE'):
            nd_ids,nds = read_abq_nodes(f)
            # RENUMBER NODES (IF NODES ARE NOT GIVEN IN ASCENDING ORDER)
            nd_ids        = [0]+nd_ids
            new_nd_labels = np.arange(max(nd_ids)+1)
            for i,nd_id in enumerate(nd_ids):
                new_nd_labels[nd_id] = i            
        if line.startswith('*ELEMENT'):
            el_type         = line.split('TYPE')[-1].split(',')[0].replace('=','').strip()            
            ele_ids,ele = read_abq_elements(f)            
            # RENUMBER ELEMENTS WRT NODE-LABELS
            ele             = new_nd_labels[np.array(ele)]            
            if 'ELSET' in line:
                el_set_nam = line.split('ELSET')[-1].split(',')[0].replace('=','').strip()
                ele_sets.update( {el_set_nam:ele_ids} )
            #
            elements.update( {el_type:{'ids':ele_ids,'conn':(np.array(ele)-1).tolist()}} )
        if line.startswith('*ELSET'):
            el_set_nam  = line.split('ELSET')[-1].split(',')[0].replace('=','').strip()
            if 'GENERATE' in line:
                option = 'GENERATE'
            else:                
                option = None
            ele_ids     = read_abq_elset(f,option=option)
            ele_sets.update( {el_set_nam:ele_ids} )
        if line.startswith('*SOLID SECTION'):
            solid_sections.append( line+'\n'+f.readline())
    part          = ABQ_part(part_name)
    part.nds      = nds
    part.elements = elements
    part.ele_sets = ele_sets
    part.sections = solid_sections
    abq_inp.parts.update( {part_name:part} )
    return abq_inp
def read_abq_material(f,abq_inp,material_name):
    '''read abaqus material section
    '''
    mat_str = ''
    while True:
        pos  = f.tell()
        line = f.readline()
        if line=='':
            break
        line = line.upper()
        if line.startswith('**'):
            continue
        if line.startswith('*STEP') or line=='':
            f.seek(pos)
            break
        if line.startswith('*MATERIAL') or line=='':
            f.seek(pos)
            break
        mat_str +=line
    abq_inp.materials.update({material_name:mat_str})    
    return abq_inp
def read_abq_elset(f,option=None):
    ''' read the elment sets
    Args: 
        f - file object 
    Returns: 
        ele_ids - element ids        
    '''
    ele_ids = []    
    while True:
        pos  = f.tell()        
        line = f.readline().replace(' ','').replace(',\n','')  
        if line.startswith('**') or line.startswith('\n'):
            continue
        if line.upper().startswith('*') or line=='':
            f.seek(pos)
            break
        tmp = map(int,line.split(','))
        if option=='GENERATE':
            ele_ids += range(tmp[0],tmp[1]+1,tmp[2])
        else:
            ele_ids += tmp
    return ele_ids

def read_abq_elements(f):
    ''' read the elements in abq format
    Args: 
        f - file object 
    Returns: 
        ele_ids - element ids
        ele    - nested integer list of element connectivity
    '''
    ele_ids = []
    ele    = []    
    while True:
        pos  = f.tell()
        line = f.readline().replace(' ','')
        if line=='\n' or line.startswith('**'):
            continue
        if line.startswith('**') or line.startswith('\n'):
            continue
        if line.startswith('*') or line=='':
            f.seek(pos)
            break
        line = line.split(',')
        ele_ids.append( int(line[0]) )
        ele.append( map(int,line[1:]) )
    return ele_ids,ele

def read_abq_nodes(f):
    ''' read the nodes in abq format
    Args: 
        f - file object 
    Returns: 
        nds - nested list of nodal coordinates
        dim - integer of nodal dimension
    '''
    nds    = []
    nd_ids = []
    while True:
        pos  = f.tell()
        line = f.readline().replace(' ','')
        if line.startswith('**') or line.startswith('\n'):
            continue
        if line.startswith('*') or line=='':
            f.seek(pos)
            break
        line_split = line.split(',')
        nd_ids.append(int(line_split[0]))
        nds.append( map(float,line_split[1:]) )        
    return nd_ids,nds
def read_abq_inp(f_name):
    '''
    '''    
    if f_name.endswith('.inp'):
        f_name = f_name[:-4] 
    abq_model = ABQ_inp(f_name)
    #
    with open(f_name+'.inp', "rU") as f:
        while True:
            line = f.readline()
            if line=='':
                break
            line = line.upper().replace('\n','')
            if line.startswith('**'): # COMMENT
                continue
            if line.startswith('*PART'):
                part_name   = line.split('NAME=')[-1].replace(' ','')
                abq_model   = read_abq_part(f,abq_model,part_name)
            if line.startswith('*MATERIAL'):
                material_name = line.split('NAME=')[-1].replace(' ','')
                abq_model     = read_abq_material(f,abq_model,material_name)
    abq_model.dim = len(abq_model.parts[abq_model.parts.keys()[0]].nds[0])
    return abq_model
def mesh2abq(f_name,nds,ele,ele_sets={},solid_sections=[]):
    ''' writes an abaqus *.inp file frome the nodes,ele, and sets array
    Args: 
       - f_name -> FILE NAME
       - nodes  -> ARRAY OF NODES
       - ele    -> ELEMENT CONNECTIVITY ARRAY
       - sets   -> ARRAY OF ELEMENTS SETS
    Returns:
       -> f_name.inp
    '''
    #
    f_name = f_name.split('.')[0]
    nds = np.array(nds)
    ele = np.array(ele)
    with open(f_name+'.inp','w') as inpFile:
        inpFile.write('*Part, name='+f_name+' \n \n')
        # NODE-SECTION
        inpFile.write('*Node \n')
        if len(nds[0])==2: w_str = '%d, %16.14F, %16.14F \n'
        else:                w_str = '%d, %16.14F, %16.14F, %16.14F \n'
        for i,nd_i in enumerate(nds):
            inpFile.write( w_str % ((i+1,)+tuple(nd_i)) )
        inpFile.write('\n')
        # ELEMENT-SECTION
        if   len(ele[0])==2:
            inpFile.write('*ELEMENT, TYPE=T3D2 \n')
            w_str = '%d, %d, %d \n'
        elif len(ele[0])==3 and len(nds[0])==2: 
            inpFile.write('*ELEMENT, TYPE=CPE3 \n')
            w_str = '%d, %d, %d, %d \n'
        elif len(ele[0])==3 and len(nds[0])==3: # SURFACE-ELEMENTS
            inpFile.write('*ELEMENT, TYPE=STRI3 \n')
            w_str = '%d, %d, %d, %d \n'
        elif len(ele[0])==4 and len(nds[0])==3:
            inpFile.write('*ELEMENT, TYPE=C3D4 \n')
            w_str = '%d, %d, %d, %d, %d \n'
        elif len(ele[0])==4 and len(nds[0])==2:
            inpFile.write('*ELEMENT, TYPE=CPE4 \n')
            w_str = '%d, %d, %d, %d, %d \n'
        elif len(ele[0])==6:
            inpFile.write('*ELEMENT, TYPE=CPE6 \n')
            w_str = '%d, %d, %d, %d, %d, %d, %d \n'
        elif len(ele[0])==10:
            inpFile.write('*ELEMENT, TYPE=C3D10 \n')
            w_str = '%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d \n'
        elif len(ele[0])==8:
            inpFile.write('*ELEMENT, TYPE=C3D8 \n')
            w_str = '%d, %d, %d, %d, %d, %d, %d, %d, %d \n'        
        for i,ele_i in enumerate(ele):
            inpFile.write( w_str % ((i+1,)+tuple(ele_i+1)) ) 
        inpFile.write('\n')
        # SET-SECTION
        for (set_nam,ele_set) in ele_sets.iteritems():            
            inpFile.write('*ELSET, ELSET='+set_nam+'\n')
            ele_set        = np.array(ele_set)+1 # ABAQUS STARTS COUNTING WITH 1
            len_ele_set10  = int(np.floor(len(ele_set)/10))
            ele_set10      = ele_set[:len_ele_set10*10]
            ele_set_remain = ele_set[len_ele_set10*10:]
            for eles_i in ele_set10.reshape((len_ele_set10,10)):
                inpFile.write(' %d, %d, %d, %d, %d, %d, %d, %d, %d, %d,\n'% tuple(eles_i))
            final_set_str  = ' '
            for ele_j in ele_set_remain:
                final_set_str += '%d, ' % ele_j
            final_set_str  = final_set_str[:-2]+'\n' # REMOVE THE LAST ,
            inpFile.write(final_set_str)
        inpFile.write(' \n')
        if solid_sections==[]:
            solid_sections = ele_sets.iterkeys()
        # ADD SECTION ACCORDING TO THE SETS
        for set_nam in solid_sections:
            inpFile.write( '*SOLID SECTION, ELSET='+set_nam
                          +', MATERIAL='+set_nam+' \n')
            inpFile.write(',\n \n')            
        # END THE PART
        inpFile.write('*End Part \n \n')        
    return
# MESH OPERATION ROUTINES
def msh_conv_tet_l2q(nds,ele):
    """
This function converts a given linear tetra mesh to a quadratic tetra mesh
Input:
    nds,    A numpy array consisting of three colums, respresenting
             the x, y and z locations of each node. (floats)
    ele,    A numpy array consisting of four colums, representing
             the nodes each element is connected to. (integers)

Output:
    new_nds,    A numpy array consisting of three colums, representing
             the x, y and z locations of each node. The nodes at the
             beginning of this list are the same as the nodes given
             as input. (floats)
    new_ele,   A numpy array consisting of ten colums, representing the
             the ten nodes each element is connected to. The order is
             such that it matches the ABAQUS C3D10 element convention:
             [(1), (2), (3), (4),(1-2),(2-3),(1-3),(1-4),(2-4),(3-4)]
             (integers)

Internal function variables:
    tetra_num,      The number of tetrahedrons in the linear mesh (int)
    edge_data,      Contains all edge information of the original mesh.
                    Edge_data consists solely out of integers.
                    The dimensions are 4 x 6*tetra_num, where each row
                    represents an edge in the original mesh.
                      Column one and two are the global indices of the
                    nodes connected to the endpoints of the edge.
                      Colum three is a representation of the location
                    of the node in the local framework.
                      Column four is an integer indicating the element
                    the edge is connected to.
    uniq,           A boolean list with size 1 x 6*tetra_num. If uniq[i]
                    is True, it means that edge edge_data[i] is unique.
                    The amount of True elements in this lists equals the
                    number of unique edges in the mesh, which equals the
                    number of nodes that should be added to the mesh.

Working princple of the function:
First, the array edge_data is created. In this array, each row represents
one edge. Edge_data is carefully sorted. If the first two elements in a
specific row are the same as the first two elements in the next row, the
edge is not unique. New nodes are created for each unique edge, after
which they are carefully assigned to the elements.
    """
    nds     = np.array(nds)
    ele     = np.array(ele)
    dim     = len(nds[0])
    ele_num = len(ele)    
    # Initialize edge_data
    if dim == 3:
        edge_data = np.ndarray(shape=(6*ele_num, 4), dtype=int)
        # All (hardcoded) possible edges on a linear tetrahedron
        edges     = np.array([[0, 1], [0, 2], [0, 3],
                              [1, 2], [1, 3], [2, 3]])
        # Ordering of the edges according to ABAQUS convention
        ele_vert = np.array([4, 6, 7, 5, 8, 9])
    elif dim == 2:
        edge_data = np.ndarray(shape=(3*ele_num, 4), dtype=int)
        # All (hardcoded) possible edges on a linear triangle element
        edges = np.array([[0, 1], [1, 2], [2, 0]])
        # Ordering of the edges according to ABAQUS convention
        ele_vert = np.array([3, 4, 5])

    # Sort the nodes that are at the endpoints of the edge before
    # storing them
    edge_vertices = np.reshape(ele[:, edges],
                               (len(edges)*ele_num, 2))
    edge_data[:, [0, 1]] = np.sort(edge_vertices)
    edge_data[:, 2]      = np.tile(ele_vert, ele_num)
    edge_data[:, 3]      = np.repeat(range(ele_num), len(edges))

    # sort edge_data array
    edge_data = edge_data[np.lexsort((edge_data[:, 1],
                                      edge_data[:, 0]))]

    # globInd are the first two colums of edge_data after it is sorted
    globInd = edge_data[:, [0, 1]]
    # Find unique rows
    tmp  = (globInd[1:]-globInd[:-1]).any(axis=1)
    # The first row is always unique, this one is added here
    uniq = np.concatenate(([True], tmp))
    ndsQ = np.ndarray((np.sum(uniq), dim), dtype=float)
    ab   = edge_data[uniq][:, [0, 1]]
    ndsQ = (nds[ab[:, 0]] + nds[ab[:, 1]])/2.0

    if   dim == 3:
        eleQ = np.ndarray((len(ele), 6), dtype=int)
    elif dim == 2:
        eleQ= np.ndarray((len(ele), 3), dtype=int)
    new_ele = np.hstack((ele, eleQ))
    # Assign nodes to elements
    nd_counter = len(nds)-1
    for idx, x in enumerate(uniq):
        if x:
            nd_counter += 1
        new_ele[edge_data[idx, 3]][edge_data[idx, 2]] = nd_counter
    # Merge original node set with the newly created nodes
    new_nds = np.concatenate((nds, ndsQ), axis=0)
    return new_nds,new_ele








def conv_lin_abq_inp_msh2_quad_abq_inp(input_file_name,out_put_file_name):
    mesh = read_abq_inp(input_file_name)
    part_name_RVE = mesh.parts.keys()[0]
    part = mesh.parts[mesh.parts.keys()[0]]
    nds, ele = msh_conv_tet_l2q(mesh.parts[ part_name_RVE ].nds,
                                mesh.parts[ part_name_RVE ].elements['C3D4']['conn'])
    # take care of abq-counting
    for key in part.ele_sets.keys():
        part.ele_sets[key] = list(np.array(part.ele_sets[key])-1)    
    mesh2abq(out_put_file_name, nds, ele,ele_sets=mesh.parts[ part_name_RVE ].ele_sets)
    return
def gen_node_set_on_RVE(abq_model,RVE_dim=[2.5,0.92,2.5]):
    ''' Generate node sets in correspondece with the RVE.
    Args:
    Returns:
        abq_model - update abq model
    '''
    dim = abq_model.dim
    # GENERATE NODE SETS
    if len(abq_model.parts.keys())!=1:
        print 'ERROR: ONLY ONE PART IS MANAGABLE WITH THIS CODE'
        return None
    part_name = abq_model.parts.keys()[0]
    #
    nds       = np.array(abq_model.parts[part_name].nds)
    nds_sets  = gen_node_sets(nds,RVE_dim)
    abq_model.parts[part_name].nds_sets.update(nds_sets)    
    return abq_model
def gen_node_sets(nds,RVE_dim):
    ''' generates a dictionary of node sets of nodes on the corners / edges and
    faces of a cuboid
    Args:
        nds       - list of nodal coordinates
        RVE_dim   - array of RVE dimension
    Returns:
        node_sets - dictionary of nodal sets
    '''
    dim       = len(RVE_dim)
    nds       = np.array(nds)
    nds_sets  = {}
    if dim==1:
        planes    = ['X0',            'X1'            ]
        dim_idxs  = [   0,              0             ]
        shift     = [ 0.0,     RVE_dim[0]             ]
    elif dim==2:
        planes    = ['X0','Y0',       'X1',       'Y1']
        dim_idxs  = [   0,   1,          0,          1]
        shift     = [ 0.0, 0.0, RVE_dim[0], RVE_dim[1]]
    else:
        planes    = ['X0','Y0','Z0',      'X1',        'Y1',       'Z1']
        dim_idxs  = [   0,   1,   2,         0,           1,          2]
        shift     = [ 0.0, 0.0, 0.0, RVE_dim[0], RVE_dim[1], RVE_dim[2]]    
    for coo,dim_idx,sh in zip(planes,dim_idxs,shift):
        if len(nds)==0:
            n_ids_imin,n_ids_ipls = np.array([]),np.array([])
        else:
            n_ids_imin = np.where( nds[:,dim_idx] < DELTA+sh )[0]
            n_ids_ipls = np.where( nds[:,dim_idx] >-DELTA+sh )[0]
        nds_sets.update( {coo+'pls':n_ids_ipls})
        nds_sets.update( {coo+'min':n_ids_imin})
    ############################################################################
    # GENERATE SETS BY 1D-integer array operations
    
    nds_sets.update({'ALL':np.arange(len(nds))})    
    # ONE HALFSPACE
    all_nd_sets = []
    space_reduction_1 = [''.join(hs) for hs in itertools.combinations(planes,1)]
    for hs in space_reduction_1:        
        nds_tmp = intersect1d((nds_sets[hs+'pls'],nds_sets[hs+'min']))        
        nds_sets.update({hs:nds_tmp})
        all_nd_sets += nds_tmp.tolist()
    nds_sets.update({'BND':np.unique(all_nd_sets)})
    nds_sets.update({'NO_BND':np.setdiff1d(nds_sets['ALL'],nds_sets['BND'])})
    if dim==1:
        return nds_sets
    
    # TWO HALFSPACES
    all_nd_sets = []
    space_reduction_2 = [sorted(hss) for hss in itertools.combinations(planes,2)]
    for hss in space_reduction_2: # hss = HALFSPACES
        # CHECK FOR HALFSPACES WITH THE SAME COORDINATE LETTER ('X0','X1')
        if len(set([dimension[0] for dimension in hss]))!= len(set(hss)):
            continue        
        nds_tmp   = intersect1d((nds_sets[hss[0]],
                                 nds_sets[hss[1]]))        
        nds_sets.update({''.join(hss):nds_tmp})
        all_nd_sets += nds_tmp.tolist()
    
    if dim==3:
        dim_reduction_key = 'EDGES'
    elif dim==2:
        dim_reduction_key = 'CORNERS'
    
    nds_sets.update({dim_reduction_key:np.unique(all_nd_sets)})
    for hs in space_reduction_1:        
        nds_sets.update({hs+'_NO_BND':np.setdiff1d(nds_sets[hs], 
                                                   nds_sets[dim_reduction_key])})    
    
    if dim==2:
        return nds_sets

    # THREE HALFSPACES
    all_nd_sets = []
    space_reduction_3 = [sorted(hss)  for hss in itertools.combinations(planes,3)]    
    for hss in space_reduction_3:
        if len(set([dimension[0] for dimension in hss]))!= len(set(hss)):
            continue                
        corner_name = ''.join(hss)
        nds_tmp     = intersect1d((nds_sets[hss[0]],
                                   nds_sets[hss[1]],
                                   nds_sets[hss[2]]))
        nds_sets.update({corner_name:nds_tmp})
        all_nd_sets += nds_tmp.tolist()
        
    nds_sets.update({'CORNERS':np.unique(all_nd_sets)})
    for hss in space_reduction_2:
        if len(set([dimension[0] for dimension in hss]))!= len(set(hss)):
            continue  
        hss_name = ''.join(hss)
        nds_sets.update({hss_name+'_NO_BND':np.setdiff1d(nds_sets[hss_name], 
                                                         nds_sets['CORNERS'])})    
    return nds_sets
def gen_PBC(abq_model,RVE_dim=np.array([2.5,0.92,2.5])):
    ''' application of periodic boundary conditions to a abqaqus model. the 
    abaqus model must be placed with its edges alligned with the COS axis, one 
    corner placed at the COS origin. The cuboid RVE extends towards positive
    COS directions.

                      7----------- Y1Z1 ---------------6
                     / |                              /|  
                    /  |                             / |  
                   /   |                            /  |  
                 X0Z1  |                         X1Z1  |  
                 /    X0Y1                        /   X1Y1 
                /      |                         /     |  
               /       |                        /      |  
              4------------- Y0Z1 -------------5       |  
              |        |                       |       |  
              |        3--------- Y1Z0 --------|-------2
              |       /                        |       /  
              |      /                         |      /   
            X0Y0    /                        X1Y0    /    
              |   X0Z0                         |   X1Z0    
              |   /                            |   /      
              |  /                             |  /       
    ^  /      | /                              | /        
   z| /y      |/                               |/         
    |/--->    0----------- Y0Z0 ---------------1
     x                     
    CORNERS: 0: X0Y0Z0 1: X1Y0Z0 2: X1Y1Z0 3: X0Y1Z0
            4: X0Y0Z1 5: X1Y0Z1 6: X1Y1Z1 7: X0Y1Z1
    EDGES:   x-direction: -  
    Y0Z0 (FRONT-DOWNSIDE); Y1Z0 (BACK-DOWNSIDE); Y1Z1 (BACK-TOP); Y0Z1 (FRONT-TOP)
            y-direction: |  
    X0Y0 (LEFT-FRONT)    ; X1Y0 (FRONT-RIGHT)  ; X1Y1 (RIGHT-BACK); X0Y1 (BACK-LEFT)
            z-direction: /  
    X0Z0 (LEFT-DOWN)     ; X1Z0 (RIGHT-DOWN)   ; X1Z1 (RIGHT-TOP) ; X0Z1 (LEFT-TOP)
    AREAS:      Y0 (FRONT); Y1 (BACK); Z1 (TOP); Z0 (DOWNSIDE); X0 (LEFT); X1 (RIGHT)
    '''
    dim      = abq_model.dim
    part_name= abq_model.parts.keys()[0]
    nds      = np.array(abq_model.parts[part_name].nds)
    nds_sets = abq_model.parts[part_name].nds_sets
    pairs    = get_coupling_pairs(nds,nds_sets,dim)
    # FORMULATE EQUATIONS
    dof_arr = [1,2]
    if dim==3:
        dof_arr +=[3]
    # ASSEMBLY
    abq_model.assembly = ABQ_assembly('MAIN_ASSEMBLY')
    abq_model.assembly.instances.update( {'RVE':abq_model.parts[part_name]} )
    abq_model          = gen_ref_nodes(abq_model,RVE_dim)
    # GENERATE EQUATIONS
    # EQUATIONS FOR PLANES / EDGES IN 2D
    abq_model = gen_hs_eqs(abq_model,pairs,dim,inst_name='RVE')
    # GENERATE EQUATIONS ON EDGES
    if dim==3:
        abq_model = gen_edge_eqs(abq_model,pairs,inst_name='RVE')
    # EQUATIONS FOR THE CORNERS
    abq_model = gen_corner_eqs(abq_model,nds_sets,dim,inst_name='RVE')
    abq_model.equations += '**'
    return abq_model
def gen_ref_nodes(abq_model,RVE_dim):
    ''' create reference nodes in the assembly
    Args:
        abq_model - object of abaqus model        
    Returns:
        abq_model - updated abq_model
    '''
    REF_sets = ['X_REF_PNT','Y_REF_PNT']
    if abq_model.dim==2:
        abq_model.assembly.nds = (np.array([[2.0,0.5],
                                            [0.5,2.0]])*RVE_dim).tolist()
    else:
        abq_model.assembly.nds = (np.array([[2.0,0.5,0.5],
                                            [0.5,2.0,0.5],
                                            [0.5,0.5,2.0]])*RVE_dim).tolist()
        abq_model.assembly.nds_sets.update({ 'Z_REF_PNT': [2]})
        REF_sets += ['Z_REF_PNT']
    abq_model.assembly.nds_sets.update({ 'X_REF_PNT': [0]})
    abq_model.assembly.nds_sets.update({ 'Y_REF_PNT': [1]})    
    return abq_model
def gen_hs_eqs(abq_model,pairs,dim,inst_name='RVE'):
    ''' adds equations to abq_model instance between halfspace boundaries
    of oposing sites of the cuboid or rectangle.
    
    u_dof1 - u_dof0 = u_REF_dof
    ##
    nd_set_str = '*Nset, nset=NODE_%d, instance='+inst_name+'\n %d\n'
    eq_str = ['*Equation\n',
              '3\n',
              'NSET1',',','DOF',', 1\n',
              'NSET0',',','DOF',',-1\n',
              'REF'  ,',','DOF',',-1\n']
    ##
    Args:
        abq_model - object of abaqus model
        pairs     - dictionary of pairs
    Returns:
        abq_model - updated abq_model
    '''    
    coos   = ['X','Y']
    dofs   = ['1','2']
    if dim ==3:
        coos += ['Z']
        dofs += ['3']
    REF_pnt_str = '_REF_PNT'        
    for coo in coos:
        if pairs[coo]==[]:
            continue
        for pair in pairs[coo]:
            abq_model,nd_set_names = gen_nd_set_abq(abq_model,pair,inst_name=inst_name)
            for dof in dofs:
                abq_model = gen_eq_by_ndset(abq_model,
                                            nd_set_names[1],nd_set_names[0],
                                            [coo.upper()+REF_pnt_str],dof)
    return abq_model
def gen_nd_set_abq(abq_model,nd_ids,inst_name='RVE',set_names=None):
    ''' generate nodeset for abaqus
    '''
    # CHECK IF nd_ids is single number
    if not hasattr(nd_ids, "__len__"):
        nd_ids = [nd_ids]
    nd_set_names  = []
    # +1 to account for abaqus counting style
    nd_ids = [nd_id+1 for nd_id in nd_ids]
    if set_names==None:
        if inst_name==None:
            nd_set_str    = '*Nset, nset=NODE_%d\n %d\n'            
        else:
            nd_set_str    = '*Nset, nset='+inst_name+'_NODE_%d, instance='+inst_name+'\n %d\n'
        for nd_id in nd_ids:
            abq_model.sets += nd_set_str % (nd_id,nd_id)
            if inst_name==None:
                nd_set_names.append( 'NODE_'+str(nd_id) )
            else:
                nd_set_names.append( inst_name+'_NODE_'+str(nd_id) )
    else:
        if not type(set_names) is list:
            set_names    = [set_names]
        if inst_name==None:
            nd_set_str     = '*Nset, nset=set_name\n %d\n'
        else:
            nd_set_str     = '*Nset, nset=set_name, instance='+inst_name+'\n %d\n'
        for nd_id,set_name in zip(nd_ids,set_names):
            abq_model.sets += (nd_set_str % (nd_id)).replace('set_name',set_name)
            nd_set_names.append( set_name )
    return (abq_model,nd_set_names)
def gen_eq_by_ndset(abq_model,nd_set_nam0,nd_set_nam1,ref_pnts,dof):
    ''' adds equation to abq_model-object in the following form:
    
        *Equation
        3
        nd_set_nam0,2, 1
        nd_set_nam1,2,-1
        Y_REF_PNT  ,2,-1
        
    Args:
        nd_set_nam0 - string of first nodeset (this one is canceled aus)
        nd_set_nam1 - string of second nodeset
        ref_pnts    - list of ref_pnts
    Returns:
        abq_model   - updated abq_model-object
    '''    
    eq_str = ['*Equation\n',
              str(len(ref_pnts)+2)+'\n',
              nd_set_nam0,',',dof,', 1\n',
              nd_set_nam1,',',dof,',-1\n']
    for ref_pnt in ref_pnts:
        eq_str += ref_pnt,',',dof,',-1\n'
    abq_model.equations += ''.join(eq_str)
    return abq_model
def gen_edge_eqs(abq_model,pairs,inst_name='RVE'):
    ''' generate PBC constraint equtaions for the edges
    
    pairs = [u^{X0Y0} u^{X1Y0} u^{X1Y1} u^{X0Y1}]
    
    u^{X1Y0}-u^{X0Y0} = u^{Rx}
    u^{X1Y1}-u^{X0Y0} = u^{Rx}+u^{Ry}
    u^{X0Y1}-u^{X0Y0} = u^{Ry}
    '''
    eq_str = ['*Equation\n',
              '3\n',
              'NSET1',',','DOF',', 1\n',
              'NSET0',',','DOF',',-1\n',
              'REF'  ,',','DOF',',-1\n']
    coos       = ['X','Y','Z']
    dofs       = ['1','2','3']
    REF_pnt_str= '_REF_PNT'
    for coo in coos:
        if pairs[coo+'_dir']==[]:
            continue        
        for pair in pairs[coo+'_dir']:
            abq_model,nd_set_names = gen_nd_set_abq(abq_model,pair,inst_name=inst_name)
            tmp = ['X','Y','Z']
            tmp.remove(coo)
            for dof in dofs:
                # FIRST SLAVE EDGE
                abq_model = gen_eq_by_ndset( abq_model,
                                             nd_set_names[1],nd_set_names[0],
                                             [tmp[0].upper()+REF_pnt_str],dof)
                # SECOND SLAVE EDGE (DIAGONAL)
                REF_pnts = [tmp[0].upper()+REF_pnt_str,
                            tmp[1].upper()+REF_pnt_str]
                abq_model = gen_eq_by_ndset( abq_model,
                                             nd_set_names[2],nd_set_names[0],
                                             REF_pnts,dof)
                # THIRD SLAVE EDGE
                abq_model = gen_eq_by_ndset( abq_model,
                                             nd_set_names[3],nd_set_names[0],
                                             [tmp[1].upper()+REF_pnt_str],dof)
    return abq_model
def gen_corner_eqs(abq_model,nds_sets,dim,inst_name='RVE'):
    ''' define equations for corners
    '''
    eq_str = ['*Equation\n',
              '3\n',
              'NSET1',',','DOF',', 1\n',
              'NSET0',',','DOF',',-1\n',
              'REF'  ,',','DOF',',-1\n']
    REF_pnt_str = '_REF_PNT'
    if list(nds_sets['CORNERS'])==[]:
        return abq_model
    if dim==2:
        coos  = ['x','y']
        dofs  = ['1','2']
        pair  = np.concatenate([nds_sets['X0Y0'],
                                nds_sets['X1Y0'],
                                nds_sets['X1Y1'],
                                nds_sets['X0Y1']])
        abq_model,nd_set_names = gen_nd_set_abq(abq_model,pair,inst_name=inst_name)
        for dof in dofs:
            # FIRST SLAVE EDGE
            gen_eq_by_ndset( abq_model,
                             nd_set_names[1],nd_set_names[0],['X_REF_PNT'],dof)
            # SECOND SLAVE EDGE (DIAGONAL)            
            REF_pnts = ['X_REF_PNT','Y_REF_PNT']
            gen_eq_by_ndset( abq_model,
                             nd_set_names[2],nd_set_names[0],REF_pnts,dof)            
            # THIRD SLAVE EDGE
            gen_eq_by_ndset( abq_model,
                             nd_set_names[3],nd_set_names[0],['Y_REF_PNT'],dof)
    else:
        # GENERATE SETS
        nd_ids = np.concatenate([nds_sets['X0Y0Z0'],
                                 nds_sets['X1Y0Z0'],
                                 nds_sets['X1Y1Z0'],
                                 nds_sets['X0Y1Z0'],
                                 nds_sets['X0Y0Z1'],
                                 nds_sets['X1Y0Z1'],
                                 nds_sets['X1Y1Z1'],
                                 nds_sets['X0Y1Z1']] )
        abq_model,nd_set_names = gen_nd_set_abq(abq_model,nd_ids,inst_name=inst_name)
        for dof in ['1','2','3']:
            abq_inp = gen_eq_by_ndset(abq_model,
                                      nd_set_names[1],nd_set_names[0],
                                      ['X_REF_PNT'],dof)
            abq_inp = gen_eq_by_ndset(abq_model,
                                      nd_set_names[3],nd_set_names[0],
                                      ['Y_REF_PNT'],dof)
            abq_inp = gen_eq_by_ndset(abq_model,
                                      nd_set_names[4],nd_set_names[0],
                                      ['Z_REF_PNT'],dof)
            #
            abq_inp = gen_eq_by_ndset(abq_model,
                                      nd_set_names[2],nd_set_names[0],
                                      ['X_REF_PNT','Y_REF_PNT'],dof)
            abq_inp = gen_eq_by_ndset(abq_model,
                                      nd_set_names[5],nd_set_names[0],
                                      ['X_REF_PNT','Z_REF_PNT'],dof)
            abq_inp = gen_eq_by_ndset(abq_model,
                                      nd_set_names[7],nd_set_names[0],
                                      ['Y_REF_PNT','Z_REF_PNT'],dof)
            #
            abq_inp = gen_eq_by_ndset(abq_model,
                                      nd_set_names[6],nd_set_names[0],
                                      ['X_REF_PNT','Y_REF_PNT','Z_REF_PNT'],dof)            
    return abq_model
def gen_bc_PBC(abq_model,bc,RVE_dim):
    '''generate boundary conditions for homogenization
    Args:
        abq_model - object that contains the data of a abaqus model
        bc        - dictionary with boundary conditions
        RVE_dim   - np.array with boundary conditions
    Returns:
        bc_string - string of boundary conditions
    '''
    bc_string = ''
    dim = abq_model.dim
    if dim==2:
        ref_pnts = ['X_REF_PNT','Y_REF_PNT']
    else:
        ref_pnts = ['X_REF_PNT','Y_REF_PNT','Z_REF_PNT']    
    # FIX ONE NODE
    part_name = 'RVE'
    if dim==2:
        corner_flag = bool(len(abq_model.parts[part_name].nds_sets['X0Y0']))
        corner_key  = 'X0Y0'
    else:
        corner_flag = bool(len(abq_model.parts[part_name].nds_sets['X0Y0Z0']))
        corner_key  = 'X0Y0Z0'
    if corner_flag:
        abq_model,set_name = gen_nd_set_abq(abq_model,
                                abq_model.parts[part_name].nds_sets[corner_key],
                                inst_name='RVE',
                                set_names='FIXED_NODE')
    else:
        # PICK A RANDOM NODE TO BE CONSTRAINED
        all_nd_ids = []
        for ele_by_type in abq_model.parts['RVE'].elements.itervalues():
            all_nd_ids += sum(ele_by_type['conn'],[])
        nd_id = np.setdiff1d(all_nd_ids,
                             abq_model.parts['RVE'].nds_sets['BND'])[0]
        abq_model,set_name = gen_nd_set_abq(abq_model,
                                [nd_id],
                                inst_name='RVE',
                                set_names='FIXED_NODE')
    bc_string += '*Boundary\n'
    for i in np.arange(dim)+1:
        bc_string += set_name[0]+','+str(i)+','+str(i)+'\n'
    ########
    # APPLY BC TO REFERENCE NODES
    if 'F' in bc.keys():
        H    = bc['F']-np.eye(dim)
        u_arr= np.dot(H,np.eye(dim)*RVE_dim) # [u_Rx,u_Ry,u_Rz]
    elif 'load_dir' in bc.keys():        
        H     = np.eye(dim)*bc['u']
        u_arr = np.dot(H,np.eye(dim)*RVE_dim)
        if   bc['load_dir']=='X':
            u_arr[1,1] = -10.1
            if dim==3: u_arr[2,2] = -10.1
        elif bc['load_dir']=='Y':
            u_arr[0,0] = -10.1
            if dim==3: u_arr[2,2] = -10.1
        elif bc['load_dir']=='Z':
            u_arr[0,0] = -10.1
            u_arr[1,1] = -10.1
    # GENERATE BC-STRING
    for i,ref_pnt in enumerate(ref_pnts):
        bc_string += '*Boundary\n'
        for j in np.arange(dim)+1:
            if u_arr[j-1,i]==-10.1:
                continue            
            tmp_str    = ref_pnt+','+str(j)+','+str(j)+',%20.16F\n' % u_arr[j-1,i]
            bc_string +=  tmp_str        
    return bc_string
def gen_step(abq_model,name='Model-1_step',nl_geo=True,step_type='static',
             step_para='0.1, 1., 1e-05, 1.\n',
             bc='',
             RVE_dim = np.array([1.,1.,1.]),
             output  = \
'''*Output, field, number interval=10
*Node Output
RT, UT
*Element Output, directions=YES
IVOL, LE, S, SDV
*Output, history, number interval=10
*Energy Output
ALLAE, ALLIE
                          '''):
    ''' generates the step
    Args:
        abq_model - obj containing the abaqus model information
        name      - string of step name
        nl_geo    - bool of nonlinear geometry flag
        step_type - string of step type
        step_para - string of step parameter
        bc        - string for boundary conditions
        RVE_dim   - np.array of RVE dimensions
        outpout   - string of output requests 
    Returns:
    '''
    if nl_geo:
        abq_model.step.append('*Step, name='+name+', nlgeom=YES\n')
    else:
        abq_model.step.append('*Step, name='+name+', nlgeom=NO\n')
    abq_model.step.append('*Static\n')
    abq_model.step.append(step_para)
    # SET BOUNDARY CONDITIONS
    abq_model.step.append( bc )
    # OUTPUT
    abq_model.step.append(output)
    abq_model.step.append('\n*END STEP')
    return abq_model
def get_coupling_pairs(nds,nds_sets,dim):
    ''' find node pairs for coupling; the node pairs have to be situated on 
    opposite sites of the RVE
    Args:
        nds      - list of all nodes
        nds_sets - dictionary of node ids that form sets
    Returns:
        pairs - dictionary of pairs
    '''
    pairs = {}
    if dim==2:
        coos   = ['X','Y']
        coo_ids= [[1],[0]]        
    else:
        coos   = ['X','Y','Z']
        coo_ids= [[1,2],[0,2],[0,1]]
    # FACES / EDGES IN 2D
    for coo,coo_id in zip(coos,coo_ids):
        if len(nds_sets[coo+'0_NO_BND'])!=len(nds_sets[coo+'1_NO_BND']):
            print len(nds_sets[coo+'0_NO_BND'])
            print len(nds_sets[coo+'1_NO_BND'])
            print 'ERROR: NUMBER OF NODES FOR '+coo+'0 AND '+coo+'1 ARE NOT MATCHING.'
            return None
        if SCIPY_FLAG:
            nd_id_pairs = get_pairs(nds[:,coo_id],
                                    nds_sets[coo+'0_NO_BND'],
                                    nds_sets[coo+'1_NO_BND'])
        else:
            nd_id_pairs = get_pairs_noKDTree(nds[:,coo_id],
                                             nds_sets[coo+'0_NO_BND'],
                                             nds_sets[coo+'1_NO_BND'])            
        pairs.update({coo:nd_id_pairs})
    # EDGES IN 3D
    if dim==3:
        # EDGES
        x_dir_edge = np.vstack( ( 
                  nds_sets['Y0Z0_NO_BND'][np.argsort(nds[nds_sets['Y0Z0_NO_BND']][:,0])],
                  nds_sets['Y1Z0_NO_BND'][np.argsort(nds[nds_sets['Y1Z0_NO_BND']][:,0])],
                  nds_sets['Y1Z1_NO_BND'][np.argsort(nds[nds_sets['Y1Z1_NO_BND']][:,0])],
                  nds_sets['Y0Z1_NO_BND'][np.argsort(nds[nds_sets['Y0Z1_NO_BND']][:,0])] ) ).T
        y_dir_edge = np.vstack( ( 
                  nds_sets['X0Z0_NO_BND'][np.argsort(nds[nds_sets['X0Z0_NO_BND']][:,1])],
                  nds_sets['X1Z0_NO_BND'][np.argsort(nds[nds_sets['X1Z0_NO_BND']][:,1])],
                  nds_sets['X1Z1_NO_BND'][np.argsort(nds[nds_sets['X1Z1_NO_BND']][:,1])],
                  nds_sets['X0Z1_NO_BND'][np.argsort(nds[nds_sets['X0Z1_NO_BND']][:,1])] ) ).T
        z_dir_edge = np.vstack( ( 
                  nds_sets['X0Y0_NO_BND'][np.argsort(nds[nds_sets['X0Y0_NO_BND']][:,2])],
                  nds_sets['X1Y0_NO_BND'][np.argsort(nds[nds_sets['X1Y0_NO_BND']][:,2])],
                  nds_sets['X1Y1_NO_BND'][np.argsort(nds[nds_sets['X1Y1_NO_BND']][:,2])],
                  nds_sets['X0Y1_NO_BND'][np.argsort(nds[nds_sets['X0Y1_NO_BND']][:,2])] ) ).T
        pairs.update( {'X_dir':x_dir_edge,
                       'Y_dir':y_dir_edge,
                       'Z_dir':z_dir_edge} )
    return pairs
def get_pairs(nds,nd_ids_A,nd_ids_B):
    ''' FIND NODAL PAIRS
    
    Args:
        nds      - list of lists with nodal coordinates
        nd_ids_A - node ids of set A
        nd_ids_B - node ids of set B
    Returns:
        pairs    - list of list of node pair ids
    '''
    ###
    nds   = np.array(nds)
    tree  = cKDTree( nds[nd_ids_A] )
    pairs = []
    for nd,id_B in zip(nds[nd_ids_B],nd_ids_B):
        ids = tree.query_ball_point(nd, DELTA)
        if len(ids)!=1:
            print 'ERROR:'
        pairs.append( [nd_ids_A[ids[0]],id_B] )
    return pairs
################################################################################
##PERIODIC BOUNDARY CONDITIONS
def write_abq_model(abq_model,f_name='TEST'):
    '''
    Args:
        abq_model - object storing all the model information
        f_name    - string of 
    Returns:
        inp_file ready for JOB submission
    '''
    float_str = '%20.16F , %20.16F'
    if abq_model.dim == 3:
        float_str += ' , %20.16F'        
    with open(f_name+'.inp','w') as f:
        f.write('**\n** PARTS\n**\n')
        for part_name,part in abq_model.parts.iteritems():
            f.write('*PART,NAME='+part_name+'\n')
            f.write('*NODE\n')
            for i,nd in enumerate(part.nds):
                f.write( ('%10d, '+float_str+'\n') % ((i+1,)+tuple(nd)) )            
            # WRITE ELEMENTS
            for ele_type,ele_info in part.elements.iteritems():
                f.write('*ELEMENT, TYPE='+ele_type+'\n')
                ele_ids = ele_info['ids']
                eles    = ele_info['conn']
                ele_str = ''.join(['%d,' for i in xrange(len(eles[0]))])[:-1]
                eles    = (np.array(eles)+1).tolist()
                for ele_id,ele in zip(ele_ids,eles):
                    f.write( ('%d,'+ele_str+'\n')%((ele_id,)+tuple(ele)))
            # WRITE NODE_SETS
            for set_nam,nd_ids in part.nds_sets.iteritems():
                if len(nd_ids)==0:
                    continue                 
                f.write('*NSET, NSET='+set_nam+'\n')
                nd_ids = (np.array(nd_ids)+1).tolist()
                nd_ids = [nd_ids[i:i+10] for i in xrange(0, len(nd_ids), 10)]
                for nd_id in nd_ids:
                    f.write( ','.join(map(str, nd_id))+'\n' )                     
            # WRITE ELE SETS
            for set_nam,ele_ids in part.ele_sets.iteritems():
                if len(ele_ids)==0:
                    continue                
                f.write('*ELSET, ELSET='+set_nam+'\n')                
                ele_ids = [ele_ids[i:i+10] for i in xrange(0, len(ele_ids), 10)]
                for ele_id in ele_ids:
                    f.write( ','.join(map(str, ele_id))+'\n' )
            # WRITE SURFACES
            for surf_name,data_line in part.srf_nds_based.iteritems():
                f.write('*SURFACE, TYPE=NODE, NAME='+surf_name+'\n')
                f.write(data_line+'\n')
            for surf_name,data_line in part.srf_ele_based.iteritems():
                f.write('*SURFACE, TYPE=ELEMENT, NAME='+surf_name+'\n')
                f.write(data_line+'\n')
            # WRITE SECTIONS
            f.write( ''.join(part.sections) )
            f.write('**\n')
            f.write('*END PART\n')
        f.write('**\n')
        # WRITE ASSEMBLY
        if abq_model.assembly!=None:
            f.write('*ASSEMBLY,NAME='+abq_model.assembly.name+'\n')
            for inst_nam,inst in abq_model.assembly.instances.iteritems():
                f.write('*INSTANCE, NAME='+inst_nam+', PART='+inst.name+'\n')
                f.write('*END INSTANCE\n')   
            f.write('**\n')
            # WRITE NODES ON ASSEMBLY:
            for i,nd in enumerate(abq_model.assembly.nds):
                f.write('*NODE\n')
                f.write( ('%10d, '+float_str+'\n') % ((i+1,)+tuple(nd)) )
            f.write('**\n')            
            # WRITE NODE SETS ON ASSEMBLY
            for set_nam,nd_ids in abq_model.assembly.nds_sets.iteritems():
                if len(nd_ids)==0:
                    continue
                f.write('*NSET,NSET='+set_nam+'\n')
                nd_ids = (np.array(nd_ids)+1).tolist() # ABAQUS COUNTING STARTS WITH 1
                nd_ids = [nd_ids[i:i+10] for i in xrange(0, len(nd_ids), 10)]
                for nd_id in nd_ids:
                    f.write( ','.join(map(str, nd_id))+'\n' )
            f.write('**\n')            
            # WRITE OTHER-SETS
            f.write(''.join(abq_model.sets))
            # WRITE EQUATIONS
            f.write(''.join(abq_model.equations))
            f.write('**\n')
            # WRITE SURFACES            
            for srf_nam,datalines in abq_model.assembly.srf_ele_based.iteritems():
                f.write('*SURFACE, TYPE=ELEMENT, NAME='+srf_nam+'\n')
                f.write(datalines)
            for srf_nam,datalines in abq_model.assembly.srf_nds_based.iteritems():
                f.write('*SURFACE, TYPE=NODE, NAME='+srf_nam+'\n')
                f.write(datalines)
            # WRITE INTERACTIONS
            f.write(''.join(abq_model.assembly.interactions))
            f.write('*END ASSEMBLY\n')
            f.write('**\n')
        # MATERIALS
        if abq_model.materials!={}:
            for mat_name,mat_str in abq_model.materials.iteritems():
                f.write('*MATERIAL, NAME='+mat_name+'\n')
                f.write(''.join(mat_str))        
        # STEPS
        if abq_model.step!=[]:
            f.write(''.join(abq_model.step))            
    return f_name
def pbc_abq(abq_model,bc,RVE_dim=[1,1,1],f_name='TEST'):
    '''
    Args:
        abq_model
        bc
        RVE_dim
        hom_type
        f_name
    Returns:
        inp_job_name     
    '''
    abq_model_loc = cPickle.loads(cPickle.dumps(abq_model, -1))
    #abq_model_loc = copy.deepcopy(abq_model)
    # RENAME RVE-SOURCE PART
    abq_model_loc.parts['RVE']      = abq_model_loc.parts.pop(abq_model.parts.keys()[0])
    abq_model_loc.parts['RVE'].name = 'RVE'    
    # GENERATE NODE_SETS
    abq_model_loc = gen_node_set_on_RVE(abq_model_loc,RVE_dim=RVE_dim)    
    abq_model_loc = gen_PBC(abq_model_loc,RVE_dim=RVE_dim) 
    bc_string     = gen_bc_PBC(abq_model_loc,bc,RVE_dim)
    if 'F' in bc.keys():
        bc_type_str = '_'+'PBC'+'_F'
    else:
        bc_type_str = '_'+'PBC'+'_'+bc['load_dir']
    # STEP
    abq_model_loc = gen_step(abq_model_loc,
                             name='QUASI_STATIC_LOADING',
                             nl_geo=False,
                             step_type='static',
                             step_para='0.1, 1., 1e-05, 1.\n',
                             bc=bc_string,
                             RVE_dim = RVE_dim)
    # WRITE
    inp_job_name = write_abq_model(abq_model_loc,f_name=f_name+bc_type_str+'_JOB') 
    return inp_job_name
def apply_PBC(inp_file,u,epsilon,isotropy_flag=True,E_Matrix=3800.,nu_Matrix=0.38,
                                                    E_Yarn  =10000.,nu_Yarn=0.15):
    # BOUNDARY CONDITIONS
    F_eps11   = np.array([[epsilon+1.0,0.0,0.0],
                          [0.0,1.0 ,0.0],
                          [0.0,0.0 ,1.0]])
    F_eps22   = np.array([[1.0,0.0,0.0],
                          [0.0,epsilon+1.0 ,0.0],
                          [0.0,0.0 ,1.0]])
    F_eps33   = np.array([[1.0,0.0,0.0],
                          [0.0,1.0 ,0.0],
                          [0.0,0.0 ,epsilon+1.0]])
    F_eps12   = np.array([[1.0,2.0*epsilon,0.0],
                          [0.0,1.0 ,0.0],
                          [0.0,0.0 ,1.0]])
    F_eps13   = np.array([[1.0,0.0,2.0*epsilon],
                          [0.0,1.0 ,0.0],
                          [0.0,0.0 ,1.0]])    
    F_eps23   = np.array([[1.0,0.0,0.0],
                          [0.0,1.0 ,2.0*epsilon],
                          [0.0,0.0 ,1.0]])
    # READ-ABQ-INPUT-FILE
    abq_model = read_abq_inp(inp_file)    
    # GET RVE-SIZE
    RVE_dim = np.max(abq_model.parts[abq_model.parts.keys()[0]].nds,axis=0)
    
    # UPDATE MATERIAL
    if isotropy_flag:
        abq_model.materials.update({'MATRIX':'*ELASTIC\n'+str(E_Matrix)+', '+str(nu_Matrix)+'\n'})
        abq_model.materials.update({'YARN':'*ELASTIC\n'+str(E_Yarn)+', '+str(nu_Yarn)+'\n'})
    # GENERATE SIX STRAIN-DEFORMATION INP-FILES
    bc= {'F':F_eps11}
    pbc_abq(abq_model,bc,RVE_dim=RVE_dim,f_name=inp_file+'_eps11')
    bc= {'F':F_eps22}
    pbc_abq(abq_model,bc,RVE_dim=RVE_dim,f_name=inp_file+'_eps22')
    bc= {'F':F_eps33}
    pbc_abq(abq_model,bc,RVE_dim=RVE_dim,f_name=inp_file+'_eps33')
    bc= {'F':F_eps12}
    pbc_abq(abq_model,bc,RVE_dim=RVE_dim,f_name=inp_file+'_eps12')
    bc= {'F':F_eps23}
    pbc_abq(abq_model,bc,RVE_dim=RVE_dim,f_name=inp_file+'_eps23')
    bc= {'F':F_eps13}
    pbc_abq(abq_model,bc,RVE_dim=RVE_dim,f_name=inp_file+'_eps13')    
    # TENSILE MODI
    bc= {'u':u,'load_dir':'X'}
    pbc_abq(abq_model,bc,RVE_dim=RVE_dim,f_name=inp_file+'_X')
    bc= {'u':u,'load_dir':'Y'}
    pbc_abq(abq_model,bc,RVE_dim=RVE_dim,f_name=inp_file+'_Y')
    bc= {'u':u,'load_dir':'Z'}
    pbc_abq(abq_model,bc,RVE_dim=RVE_dim,f_name=inp_file+'_Z')
    return

