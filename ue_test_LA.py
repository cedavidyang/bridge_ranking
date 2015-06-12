'''
Created on May 13, 2015

@author: cedavidyang
'''
import numpy as np

import pyDUE.ue_solver as ue
import pyDUE.draw_graph as d
from pyDUE.generate_graph import test_LA
from cvxopt import matrix, mul

def assign_traffic(algorithm='FW'):
    # generate the graph
    theta = matrix([0.0, 0.0, 0.0, 0.15])
    graph = test_LA(datapath='./', parameters=theta,delaytype='Polynomial')
    d.draw(graph)
    # traffic assignment
    l,x = ue.solver(graph, update=True, full=True)
    d.draw_delays(graph, l)
    delay = np.asarray([link.delay for link in graph.links.itervalues()])
    ffdelay = np.asarray([link.ffdelay for link in graph.links.itervalues()])
    edge_ratios = delay/ffdelay
    print max(mul(l,graph.get_slopes()))
    print 'cost UE:', sum([link.delay*link.flow for link in graph.links.values()])
    #l2, x2 = ue.solver(graph, update=True, full=True, SO=True)
    #d.draw_delays(graph, x2)
    ##print l2
    #print max(mul(l2,graph.get_slopes()))
    #print 'cost SO:', sum([link.delay*link.flow for link in graph.links.values()])


if __name__ == '__main__':
    assign_traffic()
