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
graph = g.test_LA(parameters=theta,delaytype=delaytype)

# capacity drop
cap_drop_array = np.ones(np.asarray(bridge_db, dtype=object).shape[0])*0.1
# get current cs distribution
cs_dist = pytraffic.condition_distribution(t, bridge_db, pmatrix)
# time of interest
time_array = np.arange(0, 110, 10)
# number of smps
nsmp = int(1e4)
# create bookkeeping dict
bookkeeping = {}
manager = Manager()
bookkeeping = manager.dict(bookkeeping)

def loop_over_time(t):
    t, smp = pytraffic.delay_history(nsmp, graph, t, bridge_db, cs_dist, cap_drop_array, theta, delaytype, bookkeeping)
    return t, smp

def main():
    start_delta_time = time.time()
    print 'CALC: Parallel version'
    try:
        pool = Pool(processes = 1)
        res = pool.map_async(loop_over_time, time_array).get(0xFFFF)
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
    time_array = np.asarray(res, dtype=object)[:,0].astype('int')
    total_delay_history = np.vstack(np.asarray(res, dtype=object)[:,1]).T/3600.
    import matplotlib.pyplot as plt
    plt.ion()
    plt.rc('font', family='serif', size=12)
    plt.rc('text', usetex=True)
    fig, ax = plt.subplots(1,1)
    ax.plot(time_array+2014, np.mean(total_delay_history,axis=1), 'o-')
    ax.set_xlim((2010, 2120))
    data_min = np.min(total_delay_history)
    data_max = np.max(total_delay_history)
    ax.set_ylim((data_min, data_max))
    ax.set_xlabel('Time')
    ax.set_ylabel('Total travel time (hour)')

    # draw vertical pdf
    xtick_max = ax.get_xticks()[-1]
    xtick_min = ax.get_xticks()[0]
    #ytick_max = ax.get_yticks()[-1]
    #ytick_min = ax.get_yticks()[0]
    left = fig.subplotpars.left
    right = fig.subplotpars.right
    top = fig.subplotpars.top
    bottom = fig.subplotpars.bottom
    def xdata_location(xdata):
        loc = (right-left)/(xtick_max-xtick_min)*(xdata-xtick_min)+left
        return loc
    for time in time_array[1:-1:2]+2014:
        ax.axvline(x=time, linestyle='--')
        loc = xdata_location(time)
        rect = [loc, bottom, 0.1*(right-left), top-bottom]
        axins = fig.add_axes(rect)
        indx = np.where(time_array+2014==time)[0][0]
        axins.hist(total_delay_history[indx, :], orientation='horizontal', bins=10, histtype='step', color='r')
        axins.set_ylim((data_min, data_max))
        axins.set_axis_off()
    plt.show()
