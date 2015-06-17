'''
Created on May 13, 2015

@author: cedavidyang
'''
import numpy as np
import pyDUE.ue_solver as ue
import pyDUE.draw_graph as d
from pyDUE.generate_graph import braess_paradox
from cvxopt import matrix, mul

def braess():
    # draw the graph
    graph = braess_paradox()
    d.draw(graph)
    # assing delay function
    l,x = ue.solver(graph, update=True, full=True)
    d.draw_delays(graph, x)
    print max(mul(l,graph.get_slopes()))
    print 'cost UE:', sum([link.delay*link.flow for link in graph.links.values()])
    l2, x2 = ue.solver(graph, update=True, full=True, SO=True)
    d.draw_delays(graph, x2)
    #print l2
    print max(mul(l2,graph.get_slopes()))
    print 'cost SO:', sum([link.delay*link.flow for link in graph.links.values()])


if __name__ == '__main__':
    braess()
