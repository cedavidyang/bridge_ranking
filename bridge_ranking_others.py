"""
Created on Thu Jun 04 12:19:41 2015

@author: cedavidyang
"""
__author__ = 'cedavidyang'

import os
import sys

import numpy as np
from scipy import stats
from pyNBI.traffic import cs2reliable

import time
import datetime
import shelve

# to restore workspace import global variables
filename = os.path.join(os.path.abspath('./'), 'Data', 'Python', 'metadata.out')
my_shelf = shelve.open(filename, 'r')
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()

# time of interest
t = 0

def bridge_safety(superCs, subCs):
    superCs = np.array(superCs)
    subCs = np.array(subCs)
    Asuper = np.zeros(np.shape(superCs))
    Asuper[superCs<=2] = superCs[superCs<=2]*55.
    Asuper[superCs==3] = superCs[superCs==3]*40.
    Asuper[superCs==4] = superCs[superCs==4]*25.
    Asuper[superCs==5] = superCs[superCs==5]*10.
    Asuper[superCs>=6] = superCs[superCs>=6]*0.

    Asub = np.zeros(np.shape(subCs))
    Asub[subCs<=2] = subCs[subCs<=2]*55.
    Asub[subCs==3] = subCs[subCs==3]*40.
    Asub[subCs==4] = subCs[subCs==4]*25.
    Asub[subCs==5] = subCs[subCs==5]*10.
    Asub[subCs>=6] = subCs[subCs>=6]*0.

    safeindx = (55.-np.minimum(Asuper, Asub))/55.

    return safeindx


def rank_safety(bridge_db):
    superCs = []; subCs = []
    for (name, lat, long, length, width, deck_cs0, super_cs0, sub_cs0, detour, onlink) in bridge_db:
        superCs.append(super_cs0)
        subCs.append(sub_cs0)
    indx = bridge_safety(superCs, subCs)
    return indx


def onebridge_pf(indx, brige_db, t, corrcoef=0.):
    name, lat, long, length, width, deck_cs0, super_cs0, sub_cs0, detour, onlink = bridge_db[indx]
    # super
    beta = cs2reliable(super_cs0)
    super_pf = stats.norm.cdf(-beta)
    # sub
    beta = cs2reliable(sub_cs0)
    sub_pf = stats.norm.cdf(-beta)
    # entire bridge, only super and sub are considered
    bridge_pf = super_pf + sub_pf - (corrcoef*np.sqrt(super_pf*(1-super_pf))*\
            np.sqrt(sub_pf*(1-sub_pf))+super_pf*sub_pf)
    return bridge_pf


def onebridge_directcost(indx, bridge_db, t):
    """compute bridge cost according to Saydam and Frangopol (2011)"""
    DISCOUNT_RATE = 0.00    # discount rate
    REBUILD_COST = 894      # rebuild cost, USD/m2
    total_cost = 0.
    name, lat, long, length, width, deck_cs0, super_cs0, sub_cs0, detour, onlink = bridge_db[indx]
    cost_reb = REBUILD_COST * length*width * (1+DISCOUNT_RATE)**t

    return cost_reb


def onebridge_indirectcost(indx, bridge_db, t, adt):
    """compute social cost according to Saydam and Frangopol (2011)"""
    TRUCK2TRAFFIC = 0.12    # Average Daily Truck Traffic (ADTT) to Average Daily Traffic (ADT)
    TRUCK_COMP = 26.97    # average compensation for truck drivers, USD/h
    CAR_OCCUPY = 1.5        # average vehicle occupancies for cars
    TRUCK_OCCUPY = 1.05     # average vehicle occupancies for trucks
    CAR_WAGE = 22.82        # average wage of car drivers, USD/h
    DISCOUNT_RATE = 0.00    # discount rate
    CAR_RUN_COST = 0.08     # running cost for cars, USD/km
    TRUCK_RUN_COST = 0.375  # running cost for trucks, USD/km
    CARGO_VALUE = 4.        # time value of a cargo

    name, lat, long, length, width, deck_cs0, super_cs0, sub_cs0, detour, onlink = bridge_db[indx]
    freespeed = 30.         # m/s

    price_time = (detour*1e3/freespeed/3600.*(1.-TRUCK2TRAFFIC)*CAR_OCCUPY*CAR_WAGE +
        detour*1e3/freespeed/3600.*TRUCK2TRAFFIC*(TRUCK_COMP*TRUCK_OCCUPY+CARGO_VALUE))*(1+DISCOUNT_RATE)**t
    price_run = (detour*(1.-TRUCK2TRAFFIC)*CAR_RUN_COST +
            detour*TRUCK2TRAFFIC*TRUCK_RUN_COST)*(1+DISCOUNT_RATE)**t
    total_cost = adt*(price_time+price_run)

    return total_cost



def rank_bridge_risk(bridge_db, adt, t):
    risks = []
    for indx in xrange(len(bridge_db)):
        pf = onebridge_pf(indx, bridge_db, t, corrcoef=0.)
        Cb = onebridge_directcost(indx, bridge_db, t)
        Cs = onebridge_indirectcost(indx, bridge_db, t, adt[indx])
        risk = pf*(Cb+Cs)
        risks.append(risk)

    return risks


if __name__ == '__main__':
    import scipy.io as sio
    import csv

    adt = []
    with open('bridge_adt.csv', 'rb') as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            adt.append(float(row[-1]))

    adt = np.array(adt)

    safety = rank_safety(bridge_db)
    risks = rank_bridge_risk(bridge_db, adt, t)

    sio.savemat('rank_others.mat', {'bridge_db':bridge_db, 'safety':safety, 'risks':risks})
