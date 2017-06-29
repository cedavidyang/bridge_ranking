"""
Created on Thu Jun 04 12:19:41 2015

@author: cedavidyang
"""
__author__ = 'cedavidyang'

import os
import sys

import numpy as np
import pyNBI.traffic as pytraffic
from pyNBI.risk import social_cost
import pyDUE.ue_solver as ue

from multiprocessing import Pool, Manager, freeze_support, Queue, Process
import Queue as queue
import itertools

import time
import datetime
import shelve

# global variables for parallel computing... stupid multiprocessing in Python

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

#def loop_over_bridges(bridge_indx, bookkeeping):
def loop_over_bridges(bridge_indx):
    indx, smp = pytraffic.delay_samples(nsmp, graph0, cost0, all_capacity, t, bridge_indx,
            bridge_db, cs_dist, cap_drop_array, theta, delaytype,
            correlation=norm_cov, nataf=nataf, corrcoef=0., x0=res0[0], bookkeeping={})

    return indx, smp

def tmpfunc(bridge_indx,q):
    indx,smp =  loop_over_bridges(bridge_indx)
    q.put((indx,smp))

if __name__ == '__main__':
    freeze_support()
    nprocess = 2
    #manager = Manager()
    #bookkeeping = manager.dict(bookkeeping)

    start_delta_time = time.time()
    print 'CALC: Parallel version'
    try:
        pool = Pool(processes = nprocess)
        #res = pool.map_async(loop_over_bridges, np.arange(bridge_db.shape[0])).get(0xFFFFFFFF)
        #res = pool.map_async(loop_over_bridges, np.arange(10)).get(0xFFFFFFFF)
        #res = pool.map_async(loop_over_bridges, np.zeros(nprocess, dtype=int)).get(0xFFFFFFFF)
        #results = [pool.apply_async(loop_over_bridges, (b,)) for b in np.arange(30)]

        #res = [r.get() for r in results]
        #res = map(loop_over_bridges, np.arange(1))
        #res = pool.map_async(loop_over_bridges,
                #itertools.izip(itertools.repeat(nsmp), itertools.repeat(graph0), itertools.repeat(cost0),
                    #itertools.repeat(all_capacity), np.arange(bridge_db.shape[0]), itertools.repeat(bridge_db),
                    #itertools.repeat(cs_dist), itertools.repeat(cap_drop_array), itertools.repeat(theta),
                    #itertools.repeat(delaytype), itertools.repeat(correlation), itertools.repeat(nataf),
                    #itertools.repeat(bookkeeping))).get(0xFFFF)
        #res = pool.map_async(tmpfunc,itertools.izip(np.arange(bridge_db.shape[0]), itertools.repeat(bookkeeping))).get(0xFFFF)
        pool.close()
        pool.join()

        #q = Queue()
        #for bridge_indx in np.zeros(10, dtype=int):
            #p = Process(target=tmpfunc, args=(bridge_indx,q))
            #p.start()
        #for bridge_indx in np.zeros(10, dtype=int):
            #p.join()
        #res = []
        #while True:
            #try:
                #res.append(q.get(timeout=10))
            #except queue.Empty:
                #break

    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        pool.terminate()
        pool.join()
        sys.exit(1)
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))

    start_delta_time = time.time()
    print 'CALC: Series version'
    res = map(loop_over_bridges, np.arange(1))
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))
    
    start_delta_time = time.time()
    print 'CALC: Series version'
    res0 = ue.solver_fw(graph0, full=True, x0=res0[0])
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))

    bridge_indx = np.asarray(res, dtype=object)[:,0].astype('int')
    bridge_risk_data = np.vstack(np.asarray(res, dtype=object)[:,1]).T
