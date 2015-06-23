"""
Created on Thu Jun 04 12:19:41 2015

@author: cedavidyang
"""
__author__ = 'cedavidyang'

import os
import sys
import psycopg2
import numpy as np
import scipy.stats as stats
import pyNBI.bridge as pybridge
import pyNBI.traffic as pytraffic

import pyDUE.generate_graph as g
import pyDUE.ue_solver as ue
from cvxopt import matrix, mul

from multiprocessing import Pool, Manager
import itertools

import time
import datetime

# global variables for parallel computing... stupid multiprocessing in Python

# open databases
conn_gis = psycopg2.connect("dbname='gisdatabase' user='amadeus' host='localhost' password=''")
cur_gis = conn_gis.cursor()
conn_nbi = psycopg2.connect("dbname='nbi' user='amadeus' host='localhost' password=''")
cur_nbi = conn_nbi.cursor()

# retrieve initial condition states of bridges
bridge_db = pytraffic.retrieve_bridge_db(cur_gis, cur_nbi)
# get transition matrix
if os.path.isfile('./pmatrix.npy'):
    pmatrix = np.load('pmatrix.npy')
else:
    pmatrix = pybridge.transition_matrix()

#create graph
theta = matrix([0.0,0.0,0.0,0.15])
delaytype = 'Polynomial'
graph0 = g.test_LA(parameters=theta,delaytype=delaytype)
nlink = len(graph0.links)

# capacity drop
cap_drop_array = np.ones(np.asarray(bridge_db, dtype=object).shape[0])*0.1
# time of interest
t = 50
# get current cs distribution
cs_dist = pytraffic.condition_distribution(t, bridge_db, pmatrix)
# number of smps
nsmp = int(1e2)
# initial capacity without failed bridges
all_capacity = np.zeros(nlink)
for link, link_indx in graph0.indlinks.iteritems():
    all_capacity[link_indx] = graph0.links[link].capacity
# initial delay
res0 = ue.solver_fw(graph0, full=True)
delay0 = res0[1][0,0]
res_bench = ue.solver(graph0)
# create bookkeeping dict
bookkeeping = {}
manager = Manager()
bookkeeping = manager.dict(bookkeeping)

def loop_over_bridges(bridge_indx):
    indx, smp = pytraffic.delay_samples(nsmp, graph0, delay0, all_capacity, bridge_indx,
            bridge_db, cs_dist, cap_drop_array, theta, delaytype, bookkeeping)
    return indx, smp

def main():
    start_delta_time = time.time()
    print 'CALC: Parallel version'
    try:
        pool = Pool(processes = 1)
        res = pool.map_async(loop_over_bridges, np.arange(bridge_db.shape[0])).get(0xFFFF)
        #res = map(loop_over_bridges, np.arange(1))
        #res = pool.map_async(loop_over_bridges,
                #itertools.izip(itertools.repeat(nsmp), itertools.repeat(t),
                    #itertools.repeat(graph), np.arange(bridge_db.shape[0]), itertools.repeat(cs_dist),
                    #itertools.repeat(cap_drop_array), itertools.repeat(theta),
                    #itertools.repeat(delaytype), itertools.repeat(bookkeeping))).get(0xFFFF)
        pool.close()
        pool.join()
    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        pool.terminate()
        pool.join()
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))
    return res

if __name__ == '__main__':
    res = main()
    bridge_indx = np.asarray(res, dtype=object)[:,0].astype('int')
    bridge_risk_data = np.vstack(np.asarray(res, dtype=object)[:,1]).T/3600.

    # postprocessing
    import matplotlib.pyplot as plt
    plt.ion()
    plt.rc('font', family='serif', size=12)
    plt.rc('text', usetex=True)

    fig, ax = plt.subplots(1,1)
    ax.boxplot(bridge_risk_data, showmeans=True)
    plt.xlabel('Bridge index')
    plt.ylabel('Risk of bridge failure (time unit)')
    xtick_label = bridge_db[bridge_indx, 0]
    ax.set_xticklabels(xtick_label, rotation='vertical')
    left = fig.subplotpars.left
    right = fig.subplotpars.right
    top = fig.subplotpars.top
    bottom = fig.subplotpars.bottom
    plt.subplots_adjust(left=left, right=right, top=top+0.07, bottom=bottom+0.07)
