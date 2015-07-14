import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import pyNBI.traffic as pytraffic
from pyNBI.risk import social_cost
from pyDUE.util import distance_on_unit_sphere
from pyNBI.bridge import int_to_degree

from multiprocessing import Pool, Manager, freeze_support, Queue, Process
import Queue as queue
import itertools

import time
import datetime
import shelve

np.random.seed(3)

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
nsmp = int(1000)

corrcoef = 0.
bridge_name = np.asarray(bridge_db, dtype=object)[:,0].astype(str)
bridge_indx = np.where(np.core.defchararray.rfind(bridge_name, '53 1304')==6)[0][0]
bridge_safety_smp, bridge_pfs0 = pytraffic.generate_bridge_safety(cs_dist,
        bridge_indx=None, correlation=norm_cov, nataf=nataf, corrcoef=corrcoef)
bridge_safety_smp, bridge_pfs = pytraffic.generate_bridge_safety(cs_dist,
        bridge_indx, norm_cov, nataf, corrcoef)

lat0 = int_to_degree(bridge_db[bridge_indx,1])
long0 = int_to_degree(bridge_db[bridge_indx,2])

distance = []
pfchange = []
for indx, (name, lat, long, length, width, deck_cs0, super_cs0, sub_cs0, detour, onlink) in enumerate(bridge_db):
    lat = int_to_degree(lat)
    long = int_to_degree(long)
    distance.append(distance_on_unit_sphere(lat,long,lat0, long0)*6371.)
    pf0 = bridge_pfs0[indx][-1]
    pf1 = bridge_pfs[indx][-1]
    if pf1>1: pf1 = 1.
    tmp = abs(pf0-pf1)
    pfchange.append(tmp)

res = np.vstack((distance,pfchange)).T
res = res[res[:,0].argsort()]
plt.rc('font', family='serif', size=12)
plt.rc('text', usetex=True)
plt.plot(res[1:,0], res[1:,1], 'o')
plt.xlabel('Distance to bridge 53 0134 (km)')
plt.ylabel('Change of failure probability, $p_f\'-pf$')
