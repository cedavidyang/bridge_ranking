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
from pyNataf.robust import semidefinitive
from pyNBI.risk import social_cost, bridge_cost
from cvxopt import matrix, mul

from multiprocessing import Pool, Manager, freeze_support
import itertools

import time
import datetime

import shelve

# load existing metadata
filename = os.path.join(os.path.abspath('./'), 'Data', 'Python', 'metadata.out')
my_shelf = shelve.open(filename, 'r')
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()


# update according to new data
# global variables for parallel computing... stupid multiprocessing in Python

# open databases
conn_gis = psycopg2.connect("dbname='gisdatabase' user='amadeus' host='localhost' password='19881229'")
cur_gis = conn_gis.cursor()
conn_nbi = psycopg2.connect("dbname='nbi' user='amadeus' host='localhost' password='19881229'")
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
graph0 = g.LA_county(parameters=theta,delaytype=delaytype, cur_gis = cur_gis)
nlink = len(graph0.links)
conn_gis.close()
conn_nbi.close()

# capacity drop
cap_drop_array = np.ones(np.asarray(bridge_db, dtype=object).shape[0])*0.5

# initial capacity without failed bridges
all_capacity = np.zeros(nlink)
for link, link_indx in graph0.indlinks.iteritems():
    all_capacity[link_indx] = graph0.links[link].capacity

# initial delay
res0 = ue.solver_fw(graph0, full=True, x0=res0[0])
delay0 = res0[1][0,0]
length_vector = np.zeros(len(graph0.links.keys()))
for link_key, link_indx in graph0.indlinks.iteritems():
    length_vector[link_indx] = graph0.links[link_key].length
distance0  = (res0[0].T * matrix(length_vector))[0,0]
#res_bench = ue.solver(graph0)

## create bookkeeping dict
#bookkeeping = {}

## correlation
corr_length = 8.73
correlation = pybridge.bridge_correlation(bridge_db, corr_length)
norm_cov = semidefinitive(correlation, tol=1e-14, deftol=1e-12)
#norm_cov = None

# with nataf
popt = np.load('nataf_popt.npy')
def nataf(x):
    return natafcurve(x,*popt)
nataf=None

if nataf is not None:
    norm_cov = nataf(correlation)
try:
    tmp = np.linalg.cholesky(norm_cov)
except:
    norm_cov = semidefinitive(norm_cov, tol=1e-14, deftol=1e-12)

filename = os.path.join(os.path.abspath('./'), 'Data', 'Python', 'metadata.out')
my_shelf = shelve.open(filename,'n') # 'n' for new
keys = ['bridge_db', 'all_capacity', 'cap_drop_array', 'res0', 'length_vector', 'pmatrix', 'theta',
        'delaytype', 'graph0', 'popt', 'nataf', 'norm_cov', 'corr_length', 'delay0', 'distance0']
for key in keys:
    try:
        my_shelf[key] = globals()[key]
    except:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        if not key.startswith("_"):
            print('ERROR shelving: {0}'.format(key))
my_shelf.close()
