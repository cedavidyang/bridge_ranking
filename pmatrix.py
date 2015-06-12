import itertools
import pyNBI

import matplotlib.pyplot as plt
from mpldatacursor import datacursor

def main():
    plt.rc('font', family='serif', size=12)
    plt.rc('text', usetex=True)
    ## debug of get_deck_CS
    #data = pyNBI.get_deck_CS(2014).values()
    #plt.hist(data, bins=np.arange(1,10)-.5)
    #plt.xticks(np.arange(1,9))
    #plt.show()

    ## transition probability of deck_CS
    #step=1
    #time_span = np.arange(1992,2013,step,dtype='int')
    #pmatrix_array = []
    #for yr in time_span:
        #pmatrix = pyNBI.transition_CS(yr, component='deck')
        #pmatrix_array.append(pmatrix)
    #pmatrix_array = np.array(pmatrix_array)
    ## postprocessing
    #for i,ms in zip(np.arange(6,1,-1), itertools.cycle('o>^+*d')):
        #labeltxt = '$a_{{ 7{} }}$'.format(i)
        #plt.plot(time_span, pmatrix_array[:,8-7,8-i], marker=ms, label=labeltxt)
    #plt.xlabel('Year', fontsize=12)
    #plt.ylabel('Transition probability')
    #annotate_text = '{label}'
    #datacursor(formatter=annotate_text.format,display='multiple', draggable=True,
            #bbox=None, fontsize=12,
            #arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))
    ##plt.plot(time_span, pmatrix_array[:,8-7,8-7], 'o-')
    ##plt.show()

    # median transition matrix
    pmatrix = pyNBI.transition_matrix(np.arange(1992,2013,dtype='int'))
    #np.save('pmatrix.npy', pmatrix)

    ## cs evolution
    #pmatrix = np.load('pmatrix.npy')
    #cs0 = np.array([1, 0, 0, 0, 0, 0, 0, 0], dtype='float')
    #service_life = np.arange(0, 102, 2)
    #cs_array = []
    #for indx, year in enumerate(service_life):
        #cs = np.dot(np.linalg.matrix_power(pmatrix.T,indx),cs0)
        #cs_array.append(cs)
    #cs_array = np.array(cs_array)
    ##plt.plot(service_life, cs_array)
    #cs=8
    #for csi,ms in zip(cs_array.T, itertools.cycle('o>^+*')):
        #plt.plot(service_life, csi, marker=ms, label='condition state {}'.format(cs))
        #cs -= 1
    #plt.xlabel('Year', fontsize=12)
    #plt.ylabel('Probability')
    #annotate_text = '{label}'
    #datacursor(formatter=annotate_text.format,display='multiple', draggable=True,
            #bbox=None, fontsize=12,
            #arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))

    ## reliability evolution
    #pmatrix = np.load('pmatrix.npy')
    #cs0 = np.array([1, 0, 0, 0, 0, 0, 0, 0], dtype='float')
    #service_life = np.arange(0, 102, 2)
    #cs_array = []
    #for indx, year in enumerate(service_life):
        #cs = np.dot(np.linalg.matrix_power(pmatrix.T,indx),cs0)
        #cs_array.append(cs)
    #cs_array = np.array(cs_array)
    #cs_mean_array = np.dot(cs_array, np.array([8,7,6,5,4,3,2,1]))
    #cs2_mean_array = np.dot(cs_array, np.array([8,7,6,5,4,3,2,1])**2)
    #cs_std_array = np.sqrt(cs2_mean_array-cs_mean_array**2)
    #beta_mean_array = (4.7-3.0)/(8-2)*(cs_mean_array-8)+4.7
    #beta_std_array = (4.7-3.0)/(8-2)*cs_std_array
    #plt.plot(service_life, beta_mean_array, 'b-', label='$\\mu_{\\beta}$')
    #plt.plot(service_life, beta_mean_array+beta_std_array, 'r--', label='$\\pm \\sigma_{\\beta}$')
    #plt.plot(service_life, beta_mean_array-beta_std_array, 'r--', label='$\\pm \\sigma_{\\beta}$')
    #plt.xlabel('Year', fontsize=12)
    #plt.ylabel('Annual reliability index')
    #annotate_text = '{label}'
    #datacursor(formatter=annotate_text.format,display='multiple', draggable=True,
            #bbox=None, fontsize=12,
            #arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k'))
