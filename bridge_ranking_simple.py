import os
import sys
import psycopg2
import numpy as np
import scipy.stats as stats
import pyNBI.bridge as pybridge
import pyNBI.traffic as pytraffic

import pyDUE.generate_graph as g
import pyDUE.ue_solver as ue
from pyNataf.nataf import natafcurve
from pyNBI.risk import compute_cost
from cvxopt import matrix, mul

import itertools

import time
import datetime

theta = matrix([1.0, 0.0, 0.0, 0.0])
delaytype = 'Polynomial'
bridge_db = np.array(
    [['      1', 1, 1, 1., 1.,  5, 6, 7, 10.0e-3, [(1, 2, 1)]],
    ['      2', 1, 1, 1., 1., 5, 6, 6, 1.0e-3, [(1, 3, 1)]],
    ['      3', 1, 1, 1., 1., 5, 6, 6, 1.0e-3, [(2, 4, 1)]],
    ['      4', 1, 1, 1., 1., 5, 6, 7, 10.0e-3, [(3, 4, 1)]]], dtype=object)
pmatrix = np.load('pmatrix.npy')
graph0 = g.braess_paradox()
nlink = len(graph0.links)
# capacity drop
cap_drop_array = np.ones(np.asarray(bridge_db, dtype=object).shape[0])*0.1
# time of interest
t = 50
# get current cs distribution
cs_dist = pytraffic.condition_distribution(t, bridge_db, pmatrix)
# number of smps
nsmp = int(1e4)
# initial capacity without failed bridges
all_capacity = np.zeros(nlink)
for link, link_indx in graph0.indlinks.iteritems():
    all_capacity[link_indx] = graph0.links[link].capacity
# initial delay
res0 = ue.solver_fw(graph0, full=True)
#delay0 = res0[1][0,0]
#length_vector = np.zeros(len(graph0.links.keys()))
#for link_key, link_indx in graph0.indlinks.iteritems():
    #length_vector[link_indx] = graph0.links[link_key].length
#distance0  = (res0[0].T * matrix(length_vector))[0,0]
#cost0 = compute_cost(bridge_db, delay0, distance0, t)
cost0 = delay0
res_bench = ue.solver(graph0)
# create bookkeeping dict
bookkeeping = {}
# correlation
correlation = None
#correlation = np.array([[1.,0.5, 0., 0.], [0.5, 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]])
#correlation = np.array([[1.,0.9, 0., 0.], [0.9, 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]])
# nataf
nataf=None
#popt = np.load('nataf_popt.npy')
#def nataf(x):
    #return natafcurve(x,*popt)


def loop_over_bridges(bridge_indx):
    indx, smp = pytraffic.delay_samples(nsmp, graph0, cost0, all_capacity, t, bridge_indx,
            bridge_db, cs_dist, cap_drop_array, theta, delaytype, correlation, nataf, bookkeeping=bookkeeping)
    return indx, smp

if __name__ == '__main__':
    # calculation
    start_delta_time = time.time()
    print 'CALC:'
    res = map(loop_over_bridges, np.arange(bridge_db.shape[0]))
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))

    bridge_indx = np.asarray(res, dtype=object)[:,0].astype('int')
    bridge_risk_data = np.vstack(np.asarray(res, dtype=object)[:,1]).T

    # postprocessing
    import matplotlib.pyplot as plt
    plt.rc('font', family='serif', size=12)
    plt.rc('text', usetex=True)

    fig, ax = plt.subplots(1,1)
    ax.boxplot(bridge_risk_data, showmeans=True)
    plt.xlabel('Bridge index')
    plt.ylabel('Risk of bridge failure (time unit)')

    # save data
    import shelve
    dir_name = os.path.join(os.path.abspath('./'), 'figures',
        'ranking_simple '+str(datetime.datetime.now()).replace(':', '-'))
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    plt.savefig(os.path.join(dir_name,'bridge_ranking_simple.eps'))
    filename=os.path.join(dir_name,'data_shelve.out')
    my_shelf = shelve.open(filename,'n') # 'n' for new
    for key in dir():
        try:
            my_shelf[key] = globals()[key]
        #except TypeError:
        except:
            #
            # __builtins__, my_shelf, and imported modules can not be shelved.
            #
            if not key.startswith("_"):
                print('ERROR shelving: {0}'.format(key))
    my_shelf.close()
    # to restore workspace, uncommon the follows
    #my_shelf = shelve.open(filename)
    #for key in my_shelf:
        #globals()[key]=my_shelf[key]
    #my_shelf.close()

    plt.ion()
    plt.show()
