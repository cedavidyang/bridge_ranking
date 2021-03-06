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
import pyNBIdatabase

import generate_graph as g
import ue_solver as ue
import draw_graph as d
from cvxopt import matrix, mul

def retrieve_bridge_db(cur_gis, cur_nbi):
    # get bridge names
    s='select structure_, link_num from network.test_LA_bridges;'
    cur_gis.execute(s)
    bridge_data = cur_gis.fetchall()
    bridge_names = np.asarray(bridge_data,dtype=object)[:,0].flatten()
    # get start node, end node and route number from test_links database
    bridge_links = np.asarray(bridge_data,dtype=object)[:,1].flatten()
    cur_gis.execute(s)
    bridge_data = cur_gis.fetchall()
    links_data = []
    for bridge_link_pair in bridge_links:
        link_pair = []
        for bridge_link in bridge_link_pair:
            s='select fromID, toID, 1 from network.test_LA_links where ID = {};'.format(bridge_link)
            cur_gis.execute(s)
            link_pair.append(cur_gis.fetchall()[0])
        links_data.append(link_pair)
    # get bridge condition ratings
    s='select cast(structure_number_008 as varchar(15)), cast(deck_cond_058 as int), \
    cast(superstructure_cond_059 as int), cast(substructure_cond_060 as int), \
    cast(detour_kilos_019 as real) \
    from nbi2014.ca2014 where ltrim(rtrim(structure_number_008)) in {} and \
    record_type_005a=\'1\';'.format(tuple(bridge_names))
    cur_nbi.execute(s)
    bridge_data = cur_nbi.fetchall()
    bridge_db = []
    for bridge_entry, link_entry in zip(bridge_data, links_data):
        tmp = list(bridge_entry)
        tmp.append(link_entry)
        bridge_db.append(tmp)

    return np.asarray(bridge_db, dtype=object)

def cs2reliable(cs):
    #beta = (4.7-3.0)/(8-2)*(cs-8)+4.7
    #beta = (3.0-1.0)/(8-2)*(cs-8)+3.0
    beta = (4.7-0.0)/(8-2)*(cs-8)+4.7
    return beta

def condition_distribution(year, cs0_data, pmatrix):
    cs_dist = []
    cs = np.arange(8,0,-1)
    if int(year) % 2 !=0:
        print 'illegal year of interest, must be even number'
    else:
        indx = int(year)/2
    for (name, deck_cs0, super_cs0, sub_cs0, dummy, dummy) in cs0_data:
        # create reliability index distribution of deck
        deck_cs0_array = (cs==deck_cs0).astype('float')
        deck_pk = np.dot(np.linalg.matrix_power(pmatrix.T,indx),deck_cs0_array)
        deck_cs_dist = stats.rv_discrete(name='deck_cs_dist', values=(cs,deck_pk))
        # create super
        super_cs0_array = (cs==super_cs0).astype('float')
        super_pk = np.dot(np.linalg.matrix_power(pmatrix.T,indx),super_cs0_array)
        super_cs_dist = stats.rv_discrete(name='super_cs_dist', values=(cs,super_pk))
        # create sub
        sub_cs0_array = (cs==sub_cs0).astype('float')
        sub_pk = np.dot(np.linalg.matrix_power(pmatrix.T,indx),sub_cs0_array)
        sub_cs_dist = stats.rv_discrete(name='sub_cs_dist', values=(cs,sub_pk))

        cs_dist.append( (name, deck_cs_dist, super_cs_dist, sub_cs_dist) )

    return cs_dist

def generate_bridge_safety(cs_dist):
    bridge_smps = []
    for (name, deck_dist, super_dist, sub_dist) in cs_dist:
        # deck
        beta = cs2reliable(deck_dist.rvs(size=1))
        deck_pf = stats.norm.cdf(-beta)
        # super
        beta = cs2reliable(super_dist.rvs(size=1))
        super_pf = stats.norm.cdf(-beta)
        # sub
        beta = cs2reliable(sub_dist.rvs(size=1))
        sub_pf = stats.norm.cdf(-beta)
        # entire bridge, only super and sub are considered and are assumed
        # independent
        bridge_pf = 1. - (1.-super_pf)*(1.-sub_pf)
        smp = stats.bernoulli.rvs(1-bridge_pf)

        bridge_smps.append( (name, smp) )

    return bridge_smps

def assign_traffic(graph, algorithm='FW', output=True):
    # traffic assignment
    l,x = ue.solver(graph, update=True, full=True)
    total_delay = sum([link.delay*link.flow for link in graph.links.values()])
    if output:
        #delay = np.asarray([link.delay for link in graph.links.itervalues()])
        #ffdelay = np.asarray([link.ffdelay for link in graph.links.itervalues()])
        #edge_ratios = delay/ffdelay
        #d.draw(graph)
        d.draw_delays(graph, l)
        #print max(mul(l,graph.get_slopes()))
        print 'cost UE:', sum([link.delay*link.flow for link in graph.links.values()])
    return total_delay

def get_initial_capacity(graph, bridge_db):
    # get initial capacity
    initial_cap=[]
    for bridge in bridge_db:
        on_links = bridge[-1]
        capacity_pair = []
        for link in on_links:
            capacity = graph.links[link].capacity
            capacity_pair.append(capacity)
        initial_cap.append(capacity_pair)
    return initial_cap

def update_links(graph,fail_bridge_db,initial_link_cap,cap_drop_array,parameters,delaytype):
    # update links
    theta = parameters
    degree = len(theta)
    to_update_links = []
    for fail_bridge,ini_caps, cap_drop in zip(fail_bridge_db, initial_link_cap, cap_drop_array):
        bridge_name = fail_bridge[0]
        bridge_detour = fail_bridge[-2]
        on_links = fail_bridge[-1]
        for on_link,ini_cap in zip(on_links,ini_caps):
            length = graph.links[on_link].length
            freespeed = graph.links[on_link].freespeed
            cap_current = graph.links[on_link].capacity
            cap_candidate = ini_cap*(1.-cap_drop)
            capacity = cap_candidate if cap_candidate<cap_current else cap_current
            slope = 1. / capacity
            length += bridge_detour*1e3
            ffdelay=length/freespeed
            coef = [ffdelay*a*b for a,b in zip(theta, np.power(slope, range(1,degree+1)))]
            to_update_links.append((on_link[0], on_link[1], on_link[2], ffdelay, (ffdelay, slope, coef),
                capacity,length,freespeed))
    graph.modify_links_from_lists(to_update_links, delaytype)

def delay_history(time_array):

    # open databases
    conn_gis = psycopg2.connect("dbname='gisdatabase' user='amadeus' host='localhost' password=''")
    cur_gis = conn_gis.cursor()
    conn_nbi = psycopg2.connect("dbname='nbi' user='amadeus' host='localhost' password=''")
    cur_nbi = conn_nbi.cursor()

    # retrieve initial condition states of bridges
    bridge_db = retrieve_bridge_db(cur_gis, cur_nbi)
    # get transition matrix
    if os.path.isfile('./pmatrix.npy'):
        pmatrix = np.load('pmatrix.npy')
    else:
        pmatrix = pyNBIdatabase.transition_matrix()

    #create graph
    theta = matrix([0.0,0.0,0.0,0.15])
    delaytype = 'Polynomial'
    graph = g.test_LA(parameters=theta,delaytype=delaytype)

    # capacity drop
    cap_drop_array = np.ones(np.asarray(bridge_db, dtype=object).shape[0])*0.1

    # loop over time_array
    total_delay_history = []
    bookkeeping = {}
    for yr in time_array:
        # get current cs distribution
        cs_dist = condition_distribution(yr, bridge_db, pmatrix)

        nSmp = int(2e4)
        # start MC
        total_delay_array = []
        for i in xrange(nSmp):
            bridge_safety_smp = generate_bridge_safety(cs_dist)
            # update link input
            bridge_safety_profile = np.asarray(bridge_safety_smp)[:,1].astype('int')
            if tuple(bridge_safety_profile) in bookkeeping.iterkeys():
                total_delay = bookkeeping[tuple(bridge_safety_profile)]
            else:
                fail_bridges = bridge_db[np.logical_not(bridge_safety_profile.astype(bool))]
                initial_link_cap = get_initial_capacity(graph, fail_bridges)
                cap_drop_after_fail = cap_drop_array[np.logical_not(bridge_safety_profile.astype(bool))]
                update_links(graph,fail_bridges,initial_link_cap,cap_drop_after_fail,theta,delaytype)
                total_delay = assign_traffic(graph, algorithm='FW', output=False)
                # save to bookkeeping
                bookkeeping[tuple(bridge_safety_profile)] = total_delay
            # add to total delay samples
            total_delay_array.append(total_delay)
        total_delay_array = np.asarray(total_delay_array)
        total_delay_history.append(total_delay_array)

    return total_delay_history


if __name__ == '__main__':
    time_array = np.arange(0, 110, 10)
    import time
    import datetime
    start_delta_time = time.time()
    total_delay_history = delay_history(time_array)
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))
    #postprocess.display_delay_history(total_delay_history)
    #total_delay_history = np.load('total_delay_history.npy')
    np.save('total_delay_history_47to0.npy', total_delay_history)
    total_delay_history = np.asarray(total_delay_history)/3600.0

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
        axins.hist(total_delay_history[indx, :], orientation='horizontal', bins=100, histtype='step', color='r')
        axins.set_ylim((data_min, data_max))
        axins.set_axis_off()
    plt.show()
