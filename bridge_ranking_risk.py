from pyDUE.generate_graph import braess_paradox
from pyDUE.ue_solver import solver_fw
import pyNBI.traffic as pytraffic

import numpy as np
import copy
from cvxopt import matrix
import time
import datetime

def update_func(graph, capacity, bridge_indx=None):
    theta = matrix([1.0, 0.0, 0.0, 0.0])
    delaytype = 'Polynomial'
    bridge_db = np.array(
        [['      1', 5, 6, 7, 10.0e-3, [(1, 2, 1)]],
        ['      2', 5, 6, 6, 1.0e-3, [(1, 3, 1)]],
        ['      3', 5, 6, 6, 1.0e-3, [(2, 4, 1)]],
        ['      4', 5, 6, 7, 10.0e-3, [(3, 4, 1)]]], dtype=object)
    pmatrix = np.load('pmatrix.npy')
    cap_drop_array = np.ones(np.asarray(bridge_db, dtype=object).shape[0])*0.0
    # time of interest
    t = 50
    # get current cs distribution
    cs_dist = pytraffic.condition_distribution(t, bridge_db, pmatrix)

    bridge_safety_smp, bridge_pfs = pytraffic.generate_bridge_safety(cs_dist)
    # update link input
    bridge_safety_profile = np.asarray(bridge_safety_smp)[:,1].astype('int')
    if bridge_indx is not None:
        bridge_safety_profile[bridge_indx] = 0
    fail_bridges = bridge_db[np.logical_not(bridge_safety_profile.astype(bool))]
    initial_link_cap = pytraffic.get_initial_capacity(graph, capacity, fail_bridges)
    cap_drop_after_fail = cap_drop_array[np.logical_not(bridge_safety_profile.astype(bool))]
    pytraffic.update_links(graph,fail_bridges,initial_link_cap,cap_drop_after_fail,theta,delaytype)

    return bridge_pfs

def main():
    # bridge importance
    graph0 = braess_paradox()
    nlink = len(graph0.links)
    nsmp_array = np.arange(10e3, 11e3, 3e3, dtype='int')
    total_delay_data = []
    # initial capacity
    capacity = np.zeros(nlink)
    for link, link_indx in graph0.indlinks.iteritems():
        capacity[link_indx] = graph0.links[link].capacity
    # delay without failed bridges
    res0 = solver_fw(graph0, full=True)
    delay0 = res0[1][0,0]

    # MC simulation
    total_delay_data=[]
    bridge_risk_data=[]
    for bridge_indx in xrange(4):
        total_delay=[]
        bridge_risk=[]
        for nsmp in nsmp_array:
            for i in xrange(nsmp):
                #update_func(graph, capacity, initial=True)
                graph = copy.deepcopy(graph0)
                bridge_pfs = update_func(graph, capacity, bridge_indx=bridge_indx)
                res = solver_fw(graph, full=True)
                total_delay.append(res[1][0,0])
                bridge_risk.append(bridge_pfs[bridge_indx][-1]*(res[1][0,0]-delay0))
        total_delay_data.append(total_delay)
        bridge_risk_data.append(bridge_risk)
    return total_delay_data, bridge_risk_data


if __name__ == '__main__':
    total_delay_data, bridge_risk_data = main()
    np.save('test_risk_47to0.npy', bridge_risk_data)
    import matplotlib.pyplot as plt
    plt.ion()
    plt.rc('font', family='serif', size=12)
    plt.rc('text', usetex=True)

    ## total delay plot
    #total_delay_data = np.asarray(total_delay_data).T
    #fig, ax = plt.subplots(1,1)
    #ax.boxplot(total_delay_data)
    #plt.xlabel('Bridge index')
    #plt.ylabel('Conditional total travel time (CTTT), time unit')
    #ax.tick_params(axis='x', labelsize=12)
    ##xtick_label = bridge_db[bridge_indx, 0]
    ##ax.set_xticklabels(xtick_label, rotation='vertical')
    ##left = fig.subplotpars.left
    ##right = fig.subplotpars.right
    ##top = fig.subplotpars.top
    ##bottom = fig.subplotpars.bottom
    ##plt.subplots_adjust(left=left, right=right, top=top+0.05, bottom=bottom+0.05)

    # bridge risk plot
    bridge_risk_data = np.asarray(bridge_risk_data).T
    fig, ax = plt.subplots(1,1)
    ax.boxplot(bridge_risk_data, showmeans=True)
    plt.xlabel('Bridge index')
    plt.ylabel('Risk of bridge failure (time unit)')
    ax.tick_params(axis='x', labelsize=12)
