"""
Created on Thu Jun 04 12:19:41 2015

@author: cedavidyang
"""
__author__ = 'cedavidyang'

import os
import sys
import shelve

import numpy as np
import scipy.io as sio
from cvxopt import matrix, mul, div

import pyNBI.traffic as pytraffic
import pyDUE.ue_solver as ue
from pyNBI.risk import bridge_cost, social_cost

import time
import datetime

# to restore workspace import global variables
filename = os.path.join(os.path.abspath('./'), 'Data', 'Python', 'metadata.out')
my_shelf = shelve.open(filename, 'r')
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()

# time of interest
t = 10
# get current cs distribution and socialcost0
cs_dist = pytraffic.condition_distribution(t, bridge_db, pmatrix)
cost0 = social_cost(delay0, distance0, t)
# number of smps
nsmp = int(5)

def extract_data(graph, res):
    fr = res[0]
    nlink = len(fr)
    lfdata = np.empty((nlink,2), dtype=object)
    for lnkkey in graph.links.iterkeys():
        indx = graph.indlinks[lnkkey]
        lfdata[indx, 0] = lnkkey
        lfdata[indx, 1] = fr[indx]

    totalDelay = res[1][0,0]
    lengthVector = np.zeros(len(graph.links.keys()))
    for link_key, link_indx in graph.indlinks.iteritems():
        lengthVector[link_indx] = graph.links[link_key].length

    return lfdata, lengthVector


def estimate_cost(res, lengthVector, bridgeIndx=None, t=t, cost0=cost0, bridge_db=bridge_db):
    totalDelay = res[1][0,0]
    totalDistance = (res[0].T * matrix(lengthVector))[0,0]
    socialCost = social_cost(totalDelay, totalDistance, t)
    if bridgeIndx is None:
        bridgeCost = 0.
    else:
        bridgeCost = bridge_cost(bridge_db[bridgeIndx][np.newaxis], t)
    return bridgeIndx, socialCost


def initial_network(graph0=graph0, res0=res0, t=t, cost0=cost0):
    start_delta_time = time.time()
    print 'CALC: Series version'
    res = ue.solver_fw(graph0, full=True, x0=res0[0])
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))

    # flow data
    lfdata, lengthVector = extract_data(graph0, res)

    # cost data
    bridgeCost, socialCost = estimate_cost(res, lengthVector, t=t, cost0=cost0)

    sio.savemat('./figs/la/flow_net0.mat',
            {'lfdata':lfdata, 'cost':socialCost})

    return locals()


def broken_network(bridge_indx, nsmp=nsmp, graph0=graph0, cost0=cost0,
        all_capacity=all_capacity, t=t, bridge_db=bridge_db, cs_dist=cs_dist,
        cap_drop_array=cap_drop_array, theta=theta, delaytype=delaytype,
        correlation=norm_cov, nataf=nataf):

    start_delta_time = time.time()
    print 'CALC: Series version'
    indx, smp, graphs, graphres, bridgeCond = pytraffic.flow_samples(nsmp, graph0, cost0, all_capacity, t, bridge_indx,
            bridge_db, cs_dist, cap_drop_array, theta, delaytype,
            correlation=norm_cov, nataf=nataf, corrcoef=0., x0=res0[0], bookkeeping={})
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))

    riskData = smp

    i = 1
    for graph, res in zip(graphs,graphres):
        # flow data
        lfdata, lengthVector = extract_data(graph, res)
        # cost data
        bridgeCost, socialCost = estimate_cost(res, lengthVector,
                bridgeIndx=bridge_indx, t=t, cost0=cost0, bridge_db=bridge_db)
        sio.savemat('./figs/la/flow_bridge{}_MC{}.mat'.format(bridge_indx, i),
                {'lfdata':lfdata, 'cost':socialCost})
        i += 1

    return locals()


if __name__ == '__main__':
    np.random.seed(1)
    # results = initial_network()
    results = broken_network(0)
