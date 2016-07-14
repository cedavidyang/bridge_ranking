import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import shelve

if __name__ == '__main__':
    plt.rc('font', family='serif', size=12)

    filename = os.path.join(os.path.abspath('./'), 'figures', 'ranking_LA 2015-07-13 09-02-58.407000', 'data_shelve.out')
    my_shelf = shelve.open(filename, 'r')
    globals()['nsmp']=my_shelf['nsmp']
    globals()['bridge_risk_data']=my_shelf['bridge_risk_data']
    globals()['bridge_db'] = my_shelf['bridge_db']
    my_shelf.close()

    folders = ['ranking_LA 2015-07-14 09-41-38.529000','ranking_LA 2015-07-15 08-30-07.403000',
             'ranking_LA 2015-07-16 09-53-42.157000', 'ranking_LA 2015-07-17 02-39-39.633756',
             'ranking_LA 2015-07-17 10-28-41.563000', 'ranking_LA 2015-07-17 14-12-51.180460',
             'ranking_LA 2015-07-18 09-03-27.937000', 'ranking_LA 2015-07-20 17-51-01.428000',
             'ranking_LA 2015-07-22 19-41-07.027000', 'ranking_LA 2015-07-26 22-58-54.961000']

    for folder in folders:
        filename = os.path.join(os.path.abspath('./'), 'figures', folder, 'data_shelve.out')
        my_shelf = shelve.open(filename, 'r')
        nsmp += my_shelf['nsmp']
        bridge_risk_data = np.vstack((bridge_risk_data, my_shelf['bridge_risk_data']))
        my_shelf.close()

    def millions(x, pos):
        'The two args are the value and tick position'
        #return '%2.1fM' % (x*1e-6)
        return '{:2.1f}M'.format(x*1e-6)
    formatter = FuncFormatter(millions)

    ## display in subfigures
    #nfig = 5
    #for f in xrange(19):
        #fig, axs = plt.subplots(nfig, 1, sharex=True)
        #for i in xrange(nfig):
            #bridge_indx = i+nfig*f
            ##bridge_indx = i
            #axs[i].yaxis.set_major_formatter(formatter)
            #axs[i].plot(np.arange(nsmp)+1, np.cumsum(bridge_risk_data.T[bridge_indx])/(np.arange(nsmp)+1))
            ##axs[i].locator_params(axis='y', nbins=3)
            #start, end = axs[i].get_ylim()
            #axs[i].yaxis.set_ticks(np.linspace(start, end, 3))
            #axs[i].get_yaxis().set_label_coords(-0.1,0.5)
            #axs[i].set_ylabel(bridge_db[bridge_indx][0].lstrip())

    # display in one figure
    fig, ax = plt.subplots(1, 1, sharex=True)
    ax.yaxis.set_major_formatter(formatter)
    mean_risk = np.mean(bridge_risk_data, axis=0)
    for bridge_indx in mean_risk.argsort()[-10:]:
        ax.plot(np.arange(nsmp)+1, np.cumsum(bridge_risk_data.T[bridge_indx])/(np.arange(nsmp)+1),
                label=bridge_db[bridge_indx][0].lstrip())
        #axs[i].locator_params(axis='y', nbins=3)

    plt.xlabel('Number of samples')
    plt.ylim((0, 100e6))
    plt.ylabel('Risk of bridge failures')
    plt.ion()
    plt.show()
