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
import copy

import pyDUE.ue_solver as ue
import pyDUE.draw_graph as d
from pyNataf.robust import semidefinitive
from pyDUE.util import distance_on_unit_sphere
from cvxopt import mul

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
    cast(lat_016 as int), cast(long_017 as int),\
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

def cs2reliable_linear(cs):
    #beta = (4.7-3.0)/(8-2)*(cs-8)+4.7
    #beta = (3.0-0.0)/(8-2)*(cs-8)+3.0
    beta = (4.7-0.0)/(8-2)*(cs-8)+4.7
    return beta

def cs2reliable(cs):
    rho = (4.84-1)/(8.-2.)*(cs-2) + 1.
    beta = (rho-1.)/np.sqrt(rho**2*0.15**2+(2.5*0.15)**2)
    return beta

def condition_distribution(year, bridge_db, pmatrix):
    cs_dist = []
    cs = np.arange(8,0,-1)
    if int(year) % 2 !=0:
        print 'illegal year of interest, must be even number'
    else:
        indx = int(year)/2
    for (name, lat, long, deck_cs0, super_cs0, sub_cs0, dummy, dummy) in bridge_db:
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

def generate_bridge_safety_deprecated(cs_dist):
    bridge_smps = []
    bridge_pfs = []
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
        bridge_pfs.append( (name, bridge_pf[0]) )

    return bridge_smps, bridge_pfs

def generate_bridge_safety(cs_dist, bridge_indx=None, correlation=None, nataf=None, corrcoef=0.):
    bridge_smps = []
    bridge_pfs = []
    if correlation is None:
        correlation = np.eye(len(cs_dist))
    # generate random field
    if nataf is None:
        norm_cov = correlation
    else:
        norm_cov = nataf(correlation)
    norm_cov = semidefinitive(norm_cov, tol=1e-14, deftol=1e-12)
    rv = stats.multivariate_normal(mean=np.zeros(len(cs_dist)), cov=norm_cov, allow_singular=True)
    field_smps = stats.norm.cdf(rv.rvs(size=1))
    # generate pf data
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
        # entire bridge, only super and sub are considered
        bridge_pf = super_pf + sub_pf - (corrcoef*np.sqrt(super_pf*(1-super_pf))*\
                np.sqrt(sub_pf*(1-sub_pf))+super_pf*sub_pf)
        bridge_pfs.append( [name, bridge_pf[0]] )
    # update pf data
    if bridge_indx is not None:
        pfe = bridge_pfs[bridge_indx][1]
        bridge_pfs[bridge_indx][1] = 1.
        for indx, bridge_pf_data in enumerate(bridge_pfs):
            rho = correlation[bridge_indx,indx]
            pf0 = bridge_pfs[indx][1]
            pf1 = (rho*np.sqrt(pfe*(1-pfe))*np.sqrt(pf0*(1-pf0))+pf0*pfe)/pfe
            bridge_pfs[indx][1] = pf1
    # generate bridge state sample
    for bridge_pf, field_smp in zip(bridge_pfs, field_smps):
        smp = bridge_pf[1]<field_smp
        # smp=True: safe; smp=False: fail
        bridge_smps.append( [name, smp] )

    return bridge_smps, bridge_pfs

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

def get_initial_capacity_deprecated(graph, bridge_db):
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

def get_initial_capacity(graph, all_capacity, bridge_db):
    # get initial capacity
    initial_cap=[]
    for bridge in bridge_db:
        on_links = bridge[-1]
        capacity_pair = []
        for link in on_links:
            capacity = all_capacity[graph.indlinks[link]]
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

def delay_samples(nsmp, graph0, delay0, all_capacity, bridge_indx, bridge_db, cs_dist,
        cap_drop_array, theta, delaytype, correlation=None, nataf=None, corrcoef=0., bookkeeping={}):
    # start MC
    bridge_risk_array=[]
    #total_delay_array = []
    for i in xrange(int(nsmp)):
        bridge_safety_smp, bridge_pfs = generate_bridge_safety(cs_dist, bridge_indx,
                correlation, nataf, corrcoef)
        # update link input
        bridge_safety_profile = np.asarray(bridge_safety_smp,dtype=object)[:,1].astype('int')
        bridge_safety_profile[bridge_indx] = 0
        #if tuple(bridge_safety_profile) in bookkeeping.iterkeys():
        #multiprocessing.managers.DictProxy has no iterkeys attribute, see http://bugs.python.org/issue9733
        if tuple(bridge_safety_profile) in iter(bookkeeping.keys()):
            total_delay = bookkeeping[tuple(bridge_safety_profile)]
            bridge_risk = bridge_pfs[bridge_indx][-1]*(total_delay-delay0)
        else:
            graph = copy.deepcopy(graph0)
            fail_bridges = bridge_db[np.logical_not(bridge_safety_profile.astype(bool))]
            initial_link_cap = get_initial_capacity(graph, all_capacity, fail_bridges)
            cap_drop_after_fail = cap_drop_array[np.logical_not(bridge_safety_profile.astype(bool))]
            update_links(graph,fail_bridges,initial_link_cap,cap_drop_after_fail,theta,delaytype)
            res = ue.solver_fw(graph, full=True)
            # save to bookkeeping
            bookkeeping[tuple(bridge_safety_profile)] = res[1][0,0]
            bridge_risk = bridge_pfs[bridge_indx][-1]*(res[1][0,0]-delay0)
        # add to total delay samples and risk samples
        #total_delay_array.append(total_delay)
        bridge_risk_array.append(bridge_risk)
    #total_delay_array = np.asarray(total_delay_array)
    bridge_risk_array = np.asarray(bridge_risk_array)

    return bridge_indx, bridge_risk_array

def delay_history(nsmp, graph, t, bridge_db, cs_dist, cap_drop_array, theta, delaytype, bookkeeping={}):

    # start MC
    total_delay_array = []
    all_capacity = np.zeros(nlink)
    for link, link_indx in graph.indlinks.iteritems():
        all_capacity[link_indx] = graph.links[link].cap
    for i in xrange(nsmp):
        bridge_safety_smp = generate_bridge_safety(cs_dist)
        # update link input
        bridge_safety_profile = np.asarray(bridge_safety_smp)[:,1].astype('int')
        if tuple(bridge_safety_profile) in iter(bookkeeping.keys()):
            total_delay = bookkeeping[tuple(bridge_safety_profile)]
        else:
            fail_bridges = bridge_db[np.logical_not(bridge_safety_profile.astype(bool))]
            initial_link_cap = get_initial_capacity(graph, all_capacity, fail_bridges)
            cap_drop_after_fail = cap_drop_array[np.logical_not(bridge_safety_profile.astype(bool))]
            update_links(graph,fail_bridges,initial_link_cap,cap_drop_after_fail,theta,delaytype)
            total_delay = assign_traffic(graph, algorithm='FW', output=False)
            # save to bookkeeping
            bookkeeping[tuple(bridge_safety_profile)] = total_delay
        # add to total delay samples
        total_delay_array.append(total_delay)
    total_delay_array = np.asarray(total_delay_array)

    return t, total_delay_array

if __name__ == '__main__':
    pass
