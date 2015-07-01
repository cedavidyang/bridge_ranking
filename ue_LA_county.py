'''
Created on May 13, 2015

@author: cedavidyang
'''
import numpy as np
import psycopg2

import pyDUE.ue_solver as ue
import pyDUE.draw_graph as d
from pyDUE.generate_graph import LA_county
from cvxopt import matrix, mul

import time
import datetime

def assign_traffic(algorithm='FW'):
    conn_gis = psycopg2.connect("dbname='gisdatabase' user='amadeus' host='localhost' password='19881229'")
    cur_gis = conn_gis.cursor()
    # generate the graph
    theta = matrix([0.0, 0.0, 0.0, 0.15])
    graph = LA_county(datapath='./', parameters=theta,delaytype='Polynomial', cur_gis=cur_gis)
    d.draw(graph)
    # traffic assignment
    #l,x = ue.solver(graph, update=True, full=True)
    start_delta_time = time.time()
    print 'CALC: traffic assignment'
    l = ue.solver_fw(graph, update=True)
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))
    d.draw_delays(graph, l, width=2)
    delay = np.asarray([link.delay for link in graph.links.itervalues()])
    ffdelay = np.asarray([link.ffdelay for link in graph.links.itervalues()])
    edge_ratios = delay/ffdelay
    #print max(mul(l,graph.get_slopes()))
    #print edge_ratios
    print 'cost UE:', sum([link.delay*link.flow for link in graph.links.values()])
    #l2, x2 = ue.solver(graph, update=True, full=True, SO=True)
    #d.draw_delays(graph, x2)
    ##print l2
    #print max(mul(l2,graph.get_slopes()))
    #print 'cost SO:', sum([link.delay*link.flow for link in graph.links.values()])


if __name__ == '__main__':
    start_delta_time = time.time()
    print 'CALC: test'
    assign_traffic()
    delta_time = time.time() - start_delta_time
    print 'DONE',str(datetime.timedelta(seconds=delta_time))

