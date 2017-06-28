import numpy as np
from cvxopt import matrix, spmatrix, solvers, spdiag, mul, div, sparse
from util import find_basis
import rank_nullspace as rn
import networkx as nx


def constraints(graph, rm_redundant = False):
    """Construct constraints for the UE link flow
    
    Parameters
    ----------
    graph: graph object
    rm_redundant: if True, remove redundant constraints for the ue solver
    
    Return value
    ------------
    Aeq, beq: equality constraints Aeq*x = beq
    """
    C, ind = nodelink_incidence(graph, rm_redundant)
    ds = [get_demands(graph, ind, id) for id,node in graph.nodes.items() if len(node.endODs) > 0]
    import ipdb; ipdb.set_trace() # BREAKPOINT
    p = len(ds)
    m,n = C.size
    Aeq, beq = spmatrix([], [], [], (p*m,p*n)), matrix(ds)
    for k in range(p): Aeq[k*m:(k+1)*m, k*n:(k+1)*n] = C
    return Aeq, beq


def nodelink_incidence(graph, rm_redundant = False):
    """
    get node-link incidence matrix

    Parameters
    ----------
    graph: graph object
    rm_redundant: if True, remove redundant constraints for the ue solver

    Return value
    ------------
    C: matrix of incidence node-link
    ind: indices of a basis formed by the rows of C
    """
    m, n = graph.numnodes, graph.numlinks
    entries, I, J = [], [], []
    for id1,node in graph.nodes.items():
        for id2,link in node.inlinks.items(): entries.append(1.0); I.append(id1-1); J.append(graph.indlinks[id2])
        for id2,link in node.outlinks.items(): entries.append(-1.0); I.append(id1-1); J.append(graph.indlinks[id2])
    C = spmatrix(entries, I, J, (m,n))
    if rm_redundant:
        M = matrix(C)
        r = rn.rank(M)
        if r < m:
            print 'Remove {} redundant constraint(s)'.format(m-r)
            ind = find_basis(M.trans())
            return C[ind,:], ind
    return C, range(m)


def linkpath_incidence(graph):
    """
    get link-path incidence matrix

    Parameters
    ----------
    graph: graph object

    Return value
    ------------
    C: matrix of incidence link-path
    """
    nnode, nlink = graph.numnodes, graph.numlinks
    npair = len(graph.ODs.keys())

    G = nx.DiGraph()
    G.add_nodes_from(graph.nodes.keys())
    G.add_edges_from([(key[0],key[1]) for key in graph.links.keys()])
    indcol = -1
    # C = np.ones((nlink, npair*nlink))
    entries, I, J = [], [], []
    for OD in graph.ODs.itervalues():
        for nodes_on_path in nx.all_simple_paths(G, OD.o, OD.d):
            graph.add_path_from_nodes(nodes_on_path)
            indcol += 1
            for u,v,route in graph.links.keys():
                indrow = graph.indlinks[(u,v,route)]
                for i in xrange(len(nodes_on_path)-1):
                    if np.array_equal([u,v], [nodes_on_path[i], nodes_on_path[i+1]]):
                        entries.append(1.0); I.append(indrow), J.append(indcol)
    C = spmatrix(entries, I, J, (nlink,indcol+1))
    return C


def pathOD_incidence(graph, haspath=False):
    """
    get path-OD incidence matrix

    Parameters
    ----------
    graph: graph object
    haspath: if True, simple has already been set before

    Return value
    ------------
    C: matrix of incidence path-OD
    """
    nnode, nlink = graph.numnodes, graph.numlinks
    npair = len(graph.ODs.keys())

    G = nx.DiGraph()
    G.add_nodes_from(graph.nodes.keys())
    G.add_edges_from([(key[0],key[1]) for key in graph.links.keys()])
    if not haspath:
        for OD in graph.ODs.itervalues():
            for nodes_on_path in nx.all_simple_paths(G, OD.o, OD.d):
                graph.add_path_from_nodes(nodes_on_path)
    npath = len(graph.paths)
    entries, I, J = [], [], []
    for OD in graph.ODs.itervalues():
        u = OD.o; v=OD.d
        indcol = graph.indods[(u,v)]
        for path in graph.paths.iterkeys():
            if (u == path[0]) and (v == path[-1]):
                indrow = graph.indpaths[path]
                entries.append(1.0); I.append(indrow), J.append(indcol)
    C = spmatrix(entries, I, J, (npath, npair))
    return C


def get_demands(graph, ind, node_id):
    """
    get demands for all OD pairs sharing the same destination
    
    Parameters
    ----------
    graph: graph object
    ind: indices of a basis formed by the rows of C
    node_id: id of the destination node
    """
    d = matrix(0.0, (graph.numnodes,1))
    for OD in graph.nodes[node_id].endODs.values():
        d[node_id-1] += OD.flow
        d[OD.o-1] = -OD.flow
    return d[ind]



def objective_poly(x, z, ks, p, w_obs=0.0, obs=None, l_obs=None, w_gap=1.0):
    """Objective function of UE program with polynomial delay functions
    f(x) = sum_i f_i(l_i) (+ 0.5*w_obs*||l[obs]-l_obs||^2)
    f_i(u) = sum_{k=1}^degree ks[i,k] u^k
    with l = sum_w x_w
    
    Parameters
    ----------
    x,z: variables for the F(x,z) function for cvxopt.solvers.cp
    ks: matrix of size (n,degree)
    p: number of w's
    w_obs: weight on the observation residual
    obs: indices of the observed links
    l_obs: observations
    """
    n, d = ks.size
    if x is None: return 0, matrix(1.0/p, (p*n,1))
    l = matrix(0.0, (n,1))
    for k in range(p): l += x[k*n:(k+1)*n]
    f, Df, H = 0.0, matrix(0.0, (1,n)), matrix(0.0, (n,1))
    for i in range(n):
        tmp = matrix(np.power(l[i],range(d+1)))
        f += ks[i,:] * tmp[1:]
        Df[i] = ks[i,:] * mul(tmp[:-1], matrix(range(1,d+1)))
        H[i] = ks[i,1:] * mul(tmp[:-2], matrix(range(2,d+1)), matrix(range(1,d)))
    if w_gap != 1.0: f, Df, H = w_gap*f, w_gap*Df, w_gap*H
    
    if w_obs > 0.0:
        num_obs, e = len(obs), l[obs]-l_obs
        f += 0.5*w_obs*e.T*e
        Df += w_obs*spmatrix(e, [0]*num_obs, obs, (1,n))
        H += spmatrix([w_obs]*num_obs, obs, [0]*num_obs, (n,1))
    
    Df = matrix([[Df]]*p)
    if z is None: return f, Df
    return f, Df, matrix([[spdiag(z[0] * H)]*p]*p)


def objective_hyper(x, z, ks, p):
    """Objective function of UE program with hyperbolic delay functions
    f(x) = sum_i f_i(v_i) with v = sum_w x_w
    f_i(u) = ks[i,0]*u - ks[i,1]*log(ks[i,2]-u)
    
    Parameters
    ----------
    x,z: variables for the F(x,z) function for cvxopt.solvers.cp
    ks: matrix of size (n,3) 
    p: number of destinations
    (we use multiple-sources single-sink node-arc formulation)
    """
    n = ks.size[0]
    if x is None: return 0, matrix(1.0/p, (p*n,1))
    l = matrix(0.0, (n,1))
    for k in range(p): l += x[k*n:(k+1)*n]
    f, Df, H = 0.0, matrix(0.0, (1,n)), matrix(0.0, (n,1))
    for i in range(n):
        tmp = 1.0/(ks[i,2]-l[i])
        f += ks[i,0]*l[i] - ks[i,1]*np.log(max(ks[i,2]-l[i], 1e-13))
        Df[i] = ks[i,0] + ks[i,1]*tmp
        H[i] = ks[i,1]*tmp**2
    Df = matrix([[Df]]*p)
    if z is None: return f, Df
    return f, Df, sparse([[spdiag(z[0] * H)]*p]*p)


def objective_hyper_SO(x, z, ks, p):
    """Objective function of SO program with hyperbolic delay functions
    f(x) = \sum_i f_i(v_i) with v = sum_w x_w
    f_i(u) = ks[i,0]*u + ks[i,1]*u/(ks[i,2]-u)
    
    Parameters
    ----------
    x,z: variables for the F(x,z) function for cvxopt.solvers.cp
    ks: matrix of size (n,3) where ks[i,j] is the j-th parameter of the delay on link i
    p: number of destinations
    (we use multiple-sources single-sink node-arc formulation)
    """
    n = ks.size[0]
    if x is None: return 0, matrix(1.0/p, (p*n,1))
    l = matrix(0.0, (n,1))
    for k in range(p): l += x[k*n:(k+1)*n]
    f, Df, H = 0.0, matrix(0.0, (1,n)), matrix(0.0, (n,1))
    for i in range(n):
        tmp = 1.0/(ks[i,2]-l[i])
        f += ks[i,0]*l[i] + ks[i,1]*l[i]*tmp
        Df[i] = ks[i,0] + ks[i,1]*tmp + ks[i,1]*l[i]*tmp**2
        H[i] = 2.0*ks[i,1]*tmp**2 + 2.0*ks[i,1]*l[i]*tmp**3
    Df = matrix([[Df]]*p)
    if z is None: return f, Df
    return f, Df, matrix([[spdiag(z[0] * H)]*p]*p)


def get_data(graph):
    """Get data for the ue solver"""
    ## TODO deprecated
    Aeq, beq = constraints(graph)
    ffdelays = graph.get_ffdelays()
    parameters, type = graph.get_parameters()
    return Aeq, beq, ffdelays, parameters, type


