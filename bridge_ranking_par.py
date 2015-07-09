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
from pyNataf.nataf import natafcurve
from pyNBI.risk import social_cost, bridge_cost
from cvxopt import matrix, mul

from multiprocessing import Pool, Manager, freeze_support, Queue, Process
import Queue as queue
import itertools

import time
import datetime
import shelve

# global variables for parallel computing... stupid multiprocessing in Python

# to restore workspace, uncommon the follows
filename = os.path.join(os.path.abspath('./'), 'Data', 'Python', 'metadata.out')
my_shelf = shelve.open(filename, 'r')
for key in my_shelf:
    #globals()[key]=my_shelf[key]
    bridge_db = my_shelf['bridge_db']
    pmatrix = my_shelf['pmatrix']
    theta = my_shelf['theta']
    delaytype = my_shelf['delaytype']
    graph0 = my_shelf['graph0']
    all_capacity = my_shelf['all_capacity']
    length_vector = my_shelf['length_vector']
    popt = my_shelf['popt']
    res0 = my_shelf['res0']
my_shelf.close()


nlink = len(graph0.links)
cap_drop_array = np.ones(np.asarray(bridge_db, dtype=object).shape[0])*0.1
# time of interest
t = 50
# get current cs distribution
cs_dist = pytraffic.condition_distribution(t, bridge_db, pmatrix)
# number of smps
nsmp = int(1000)
delay0 = res0[1][0,0]
distance0  = (res0[0].T * matrix(length_vector))[0,0]
cost0 = social_cost(delay0, distance0, t)
#res_bench = ue.solver(graph0)
# correlation
corr_length = 8.73
correlation = pybridge.bridge_correlation(bridge_db, corr_length)
#correlation = None
# nataf
#def nataf(x):
    #return natafcurve(x,*popt)
nataf=None
## create bookkeeping dict
#bookkeeping = {}

#def loop_over_bridges(bridge_indx, bookkeeping):
def loop_over_bridges(bridge_indx):
    indx, smp = pytraffic.delay_samples(nsmp, graph0, cost0, all_capacity, t, bridge_indx,
            bridge_db, cs_dist, cap_drop_array, theta, delaytype, correlation, nataf, bookkeeping={})

    return indx, smp

def tmpfunc(bridge_indx,q):
    indx,smp =  loop_over_bridges(bridge_indx)
    q.put((indx,smp))

if __name__ == '__main__':
    freeze_support()
    #manager = Manager()
    #bookkeeping = manager.dict(bookkeeping)

    start_delta_time = time.time()
    print 'CALC: Parallel version'
    try:
        #pool = Pool(processes = 60)
        #res = pool.map_async(loop_over_bridges, np.arange(bridge_db.shape[0])).get(0xFFFFFFFF)
        #res = pool.map_async(loop_over_bridges, np.arange(60)).get(0xFFFFFFFF)
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
        #pool.close()
        #pool.join()

        q = Queue()
        for bridge_indx in np.arange(10):
            print "sub"
            p = Process(target=tmpfunc, args=(bridge_indx,q))
            p.start()
        for bridge_indx in np.arange(10):
            p.join()
        res = []
        while True:
            try:
                res.append(q.get(timeout=1))
            except queue.Empty:
                break

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

    bridge_indx = np.asarray(res, dtype=object)[:,0].astype('int')
    bridge_risk_data = np.vstack(np.asarray(res, dtype=object)[:,1]).T

    # postprocessing
    import matplotlib.pyplot as plt
    plt.ion()
    plt.rc('font', family='serif', size=12)
    #plt.rc('text', usetex=True)

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

    # save data
    import shelve
    dir_name = os.path.join(os.path.abspath('./'), 'figures',
        'ranking_LA '+str(datetime.datetime.now()).replace(':', '-'))
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    plt.savefig(os.path.join(dir_name,'bridge_ranking_LA.eps'))
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
