import os
import sys

import numpy as np
import itertools

FONT_SIZE = 9

import matplotlib as mpl
pgf_with_custom_preamble = {
    "figure.figsize": [3.54, 2.655],
    "figure.subplot.bottom": 0.15,
    "figure.subplot.top": 0.94,
    "figure.subplot.left": 0.13,
    "figure.subplot.right": 0.94,
    "font.family": "serif", # use serif/main font for text elements
    "font.size": FONT_SIZE, # use font size
    "text.usetex": True,    # use inline math for ticks
    'text.latex.unicode': True,
    #"pgf.rcfonts": False,   # don't setup fonts from rc parameters
    #"pgf.preamble": [
        #r'\usepackage{fontspec}',         # load additional packages
        #]
}
mpl.rcParams.update(pgf_with_custom_preamble)
import matplotlib.pyplot as plt
plt.ion()
# from mpldatacursor import datacursor

def CSprobdevelopment():
    # cs evolution
    pmatrix = np.load('pmatrix.npy')
    pmatrix = pmatrix.item()['sub'][:-1,:-1]
    cs0 = np.array([1, 0, 0, 0, 0, 0, 0], dtype='float')
    service_life = np.arange(0, 102, 2)
    cs_array = []
    for indx, year in enumerate(service_life):
        cs = np.dot(np.linalg.matrix_power(pmatrix.T,indx),cs0)
        cs_array.append(cs)
    cs_array = np.array(cs_array)
    #plt.plot(service_life, cs_array)
    cs=7
    for csi,ms in zip(cs_array.T, itertools.cycle('o>^+*')):
        if cs == 7:
            plt.plot(service_life, csi, marker=ms, label='Condition rating $\ge8$')
        elif cs == 1:
            plt.plot(service_life, csi, marker=ms, label='Condition rating $\le2$')
        else:
            plt.plot(service_life, csi, marker=ms, label='Condition rating {}'.format(cs+1))
        cs -= 1
    plt.xlabel('Service time (year)', fontsize=FONT_SIZE)
    plt.ylabel('Probability')
    annotate_text = '{label}'
    datacursor(formatter=annotate_text.format,display='multiple', draggable=True,
            bbox=None, fontsize=FONT_SIZE,
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))


def CSdevelopment(life=100, intl=2, struct='super'):
    # cs evolution
    pmatrix = np.load('pmatrix.npy')
    pmatrix = pmatrix.item()[struct][:-1,:-1]
    cs0 = np.array([1, 0, 0, 0, 0, 0, 0], dtype='float')
    service_life = np.arange(0, life+intl, intl)
    cs_array = []
    for indx, year in enumerate(service_life):
        cs = np.dot(np.linalg.matrix_power(pmatrix.T,indx),cs0)
        cs_array.append(cs)
    cs_array = np.array(cs_array)
    cs=7
    cmap = mpl.cm.get_cmap('Greys')
    cidx = np.linspace(0.1,0.9,7)
    lgdrec = []
    for csi,fc in zip(cs_array.T, cidx):
        if cs == 7:
            cscum = csi
            plt.fill_between(service_life, cscum, y2=0, facecolor=cmap(fc))
            lgdrec.append(plt.Rectangle((0, 0), 1, 1, fc=cmap(fc), label='$\ge8$'))
        elif cs == 1:
            lastcscum,cscum = cscum, cscum+csi
            plt.fill_between(service_life, cscum, y2=lastcscum, facecolor=cmap(fc), label='$\le2$')
            lgdrec.append(plt.Rectangle((0, 0), 1, 1, fc=cmap(fc), label='$\le2$'))
        else:
            lastcscum,cscum = cscum, cscum+csi
            plt.fill_between(service_life, cscum, y2=lastcscum, facecolor=cmap(fc), label='Condition rating {}'.format(cs+1))
            lgdrec.append(plt.Rectangle((0, 0), 1, 1, fc=cmap(fc), label='Condition rating {}'.format(cs+1)))
        cs -= 1
    import scipy.io as sio
    sio.savemat('csPostprocess.mat', {'cs': cs_array, 'life':service_life, 'pmatrix':pmatrix})

    plt.xlabel('Service time (year)', fontsize=FONT_SIZE)
    plt.ylabel('Probability')
    plt.xlim([0, life])
    plt.ylim([0, 1.0])
    plt.legend(lgdrec, ['$\ge8$', '7', '6', '5', '4', '3', '$\le2$'],fontsize=FONT_SIZE, title='Condition Rating', ncol=2).draggable(state=True)


def plotMStobeta():
    fig, ax = plt.subplots()
    msarray = np.arange(1,8)
    msname = ['MS1', 'MS2', 'MS3', 'MS4', 'MS5', 'MS6', 'MS7']
    csname = ['CR $\ge8$', 'CR 7', 'CR 6', 'CR 5', 'CR 4', 'CR 3', 'CR $\le2$']
    Rarray = (4.84-1)/(1.-7)*(msarray-1)+4.84
    betaarray = (Rarray-1)/np.sqrt(Rarray**2*0.15**2+0.375**2)
    bar1 = ax.bar(msarray, betaarray)
    ax.set_xticks(msarray + 0.5)
    ax.set_xticklabels(msname)
    plt.ylim((0,6))
    plt.xlim((0.4,8.1))
    def autolabel(rects, labels=None):
        # attach some text labels
        for rect,lb in zip(rects,labels):
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                lb, ha='center', va='bottom')
    autolabel(bar1, csname)
    plt.ylabel('Reliability index')
    plt.xlabel('Markovian state')


def nataftest():
    from pyNataf.nataf import nataf_transform
    import scipy.stats as stats
    np.random.seed(1)
    xrv = stats.uniform()
    rho1_array = np.arange(0., 1.0, 0.01)
    rho0_array = []
    for rho1 in rho1_array:
        rho0_array.append(nataf_transform(rho1, xrv))
    rho0_array = np.asarray(rho0_array)
    plt.plot(rho0_array, rho1_array, 'b.')
    plt.plot([0.,1.], [0.,1.], 'r-')
    plt.grid()
    plt.xlabel('Original correlation $\\rho$')
    plt.ylabel('Adjusted correlation $\\rho\'$')


def postpfvsdist(year, checkname):
    import shelve
    import pyNBI.traffic as pytraffic
    from pyDUE.util import distance_on_unit_sphere, int_to_degree
    # year of interest
    t = year
    # to restore workspace import global variables
    filename = os.path.join(os.path.abspath('./'), 'Data', 'Python', 'metadata.out')
    my_shelf = shelve.open(filename, 'r')
    for key in my_shelf:
        globals()[key]=my_shelf[key]
    my_shelf.close()

    # get current cs distribution and socialcost0
    cs_dist = pytraffic.condition_distribution(t, bridge_db, pmatrix)
    corrcoef = 0.
    bridge_name = np.asarray(bridge_db, dtype=object)[:,0].astype(str)
    bridge_indx = np.where(np.core.defchararray.rfind(bridge_name, checkname)==6)[0][0]
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

    resunsorted = np.vstack((distance,pfchange)).T
    res = resunsorted[resunsorted[:,0].argsort()]
    # plt.plot(res[1:,0], res[1:,1], 'o')
    # plt.xlabel('Distance to bridge 53 0134 (km)')
    # plt.ylabel('Change of failure probability, $p_f\'-pf$')

    import scipy.io as sio
    sio.savemat('./figs/pfchange_'+checkname[-4:]+'.mat',
            {'res':res, 'unsorted':resunsorted})

    return locals()


def mc_convergence(subfig=False):
    import shelve
    from matplotlib.ticker import FuncFormatter
    # load data
    filename = os.path.join(os.path.abspath('./'), 'figures', 'ranking_LA 2015-07-13 09-02-58.407000', 'data_shelve.out')
    my_shelf = shelve.open(filename, 'r')
    nsmp=my_shelf['nsmp']
    bridge_risk_data=my_shelf['bridge_risk_data']
    bridge_db = my_shelf['bridge_db']
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

    if subfig:
        # display in subfigures
        nfig = 5
        for f in xrange(19):
            fig, axs = plt.subplots(nfig, 1, sharex=True)
            for i in xrange(nfig):
                bridge_indx = i+nfig*f
                #bridge_indx = i
                axs[i].yaxis.set_major_formatter(formatter)
                axs[i].plot(np.arange(nsmp)+1, np.cumsum(bridge_risk_data.T[bridge_indx])/(np.arange(nsmp)+1))
                #axs[i].locator_params(axis='y', nbins=3)
                start, end = axs[i].get_ylim()
                axs[i].yaxis.set_ticks(np.linspace(start, end, 3))
                axs[i].get_yaxis().set_label_coords(-0.1,0.5)
                axs[i].set_ylabel(bridge_db[bridge_indx][0].lstrip())
        plt.xlabel('Number of samples')
    else:
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

if __name__ == '__main__':
    results = postpfvsdist(0, '53 1176')
