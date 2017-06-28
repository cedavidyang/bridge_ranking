# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 16:15:19 2014

@author: hugo
"""
import sys
import numpy as np
import manips as m
from util import distance_on_unit_sphere, closest_node


def Create_ODs_nodes_unique(nodes=None, datapath='', cur_gis=None):
    M=m.Manips(datapath)
    List_TAZ, List_TAZ_ids = M.List_TAZ, M.List_TAZ_ids
    ODs_TAZ = np.asarray(M.Read_TAZ_from_csv())
    TAZ_2_node = Create_dict_TAZ_2_node(List_TAZ, List_TAZ_ids, nodes, cur_gis)
    ODs_nodes_multiple = []
    for od in ODs_TAZ:
        startnode, endnode = TAZ_2_node[od[0]], TAZ_2_node[od[1]]
        ODs_nodes_multiple.append([startnode, endnode, od[2]])
    ODs_nodes = Sum_multiple_ODs(ODs_nodes_multiple)
    return np.asarray(ODs_nodes)

def Create_dict_TAZ_2_node(List_TAZ, List_TAZ_ids, nodes=None, cur_gis=None):
    dict = {}
    for j,i in enumerate(List_TAZ_ids):
        if cur_gis is not None:
            s = ("select (select v.id from network.la_voronoi v "
                        "where ST_contains(v.geom,ST_centroid(t.geom)) limit 1 ) "
                "from taz.ca_taz_2009 t "
                " where t.gid = {};").format(i)
            cur_gis.execute(s)
            try:
                dict[i] = int(cur_gis.fetchall()[0][0])
            except TypeError:
                dict[i] = closest_node(List_TAZ[j][1], List_TAZ[j][2], nodes)
        elif nodes is not None:
            dict[i] = closest_node(List_TAZ[j][1], List_TAZ[j][2], nodes)
        else:
            print "nodes or cur_gis must be given"
            sys.exit(1)
    return dict

def Sum_multiple_ODs(ODs_multiple):
    dict_unique = {}
    ODs_unique = []
    for od in ODs_multiple:
        if od[0]!=od[1]:
            if dict_unique.has_key((od[0], od[1])):
                dict_unique[(od[0], od[1])] += od[2]/2. #We divide by 2 because peak hour lasts 2 hours
            else : dict_unique[(od[0], od[1])] = od[2]/2.
    for od, value in dict_unique.iteritems():
        ODs_unique.append([od[0], od[1], value])
    return ODs_unique

if __name__ == "__main__":
    nodes = np.genfromtxt('Data/Network/CSV/LA_big_box_arterials/nodes_LA_toy.csv', delimiter = ',', skiprows = 1)
    nodes = nodes[:,1:3]
    ODs = Create_ODs_nodes_unique(nodes)

    def Get_key(item):
        return item[2]
    M = m.Manips()
    ODs_sorted = np.asarray(sorted(ODs, key = Get_key))
    temp = ODs_sorted
    ODs_sorted = []

    for od in temp:
        startnode = nodes[od[0]-1]
        endnode = nodes[od[1]-1]
        ODs_sorted.append([od[0], od[1], od[2], distance_on_unit_sphere(nodes[od[0]-1][1], nodes[od[0]-1][0], nodes[od[1]-1][1], nodes[od[1]-1][0]), m.Is_in_I210box(startnode[1], startnode[0], 'medium'), m.Is_in_I210box(endnode[1], endnode[0], 'medium')])
    ODs_sorted = np.asarray(ODs_sorted)

    temp = ODs_sorted
    ODs_sorted  = []

    for od in temp:
        if od[4]+od[5]>1 : ODs_sorted.append(od[0:4])
    ODs_sorted = np.asarray(ODs_sorted)
