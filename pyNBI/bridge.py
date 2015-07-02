'''
Created on May 22, 2015

@author: cedavidyang
'''

__author__ = 'cedavidyang'

import psycopg2
import numpy as np
import matplotlib.pyplot as plt
from pyDUE.util import distance_on_unit_sphere

def select_data(db, query_name, query_condition):
    """ get deck condition state data according to query dictionary """
    # open database
    conn = psycopg2.connect("dbname='nbi' user='amadeus' host='localhost' password=''")
    cur = conn.cursor()
    command = 'SELECT {} FROM {}.{} T WHERE {};'.format(query_name, db['schema'], db['table'], query_condition)
    print command
    cur.execute(command)
    data = cur.fetchall()

    return np.array(data)

def get_deck_CS(year, structure_kind=None):
    db = {'schema':'nbi{}'.format(year), 'table':'ca{}'.format(year)}
    query_name = 'T.STRUCTURE_NUMBER_008, T.DECK_COND_058'
    ## concrete and concrete (continuous) bridges
    #query_condition = '(T.STRUCTURE_KIND_043a=\'1\' or T.STRUCTURE_KIND_043a=\'2\') '
    # select specific bridge types
    query_condition = 'T.STRUCTURE_KIND_043a NOT IN (\'7\', \'8\', \'9\') '
    # built before the database
    query_condition = query_condition + 'and CAST(COALESCE(NULLIF(T.YEAR_BUILT_027,\'    \'),\'0\') AS INT)<1992'
    # remove invalide condition states
    query_condition = query_condition + 'and T.DECK_COND_058 NOT IN (\'N\', \'0\', \'99\')'
    data = select_data(db,query_name, query_condition)
    print '# of entries: {}'.format(data.shape[0])
    name = data[:,0]
    value = data[:,1].flatten()
    name = name[value!=' ']
    value = value[value!=' ']
    value = value.astype('int')
    data_dict = dict(zip(name, value))

    return data_dict

def transition_deck_CS_deprecated(year_i, year_j):
    # data at year i
    db_i = {'schema':'nbi{}'.format(year_i), 'table':'ca{}'.format(year_i)}
    query_name_i = 'T.STRUCTURE_NUMBER_008, T.RECORD_TYPE_005A, T.DECK_COND_058'
    # concrete and concrete (continuous) bridges
    query_condition_i = '(T.STRUCTURE_KIND_043a=\'1\' or T.STRUCTURE_KIND_043a=\'2\') '
    data_i = select_data(db_i,query_name_i, query_condition_i)
    # data at year j
    db_j = {'schema':'nbi{}'.format(year_j), 'table':'ca{}'.format(year_j)}
    query_name_j = 'T.STRUCTURE_NUMBER_008, T.RECORD_TYPE_005A, T.DECK_COND_058'
    # concrete and concrete (continuous) bridges
    query_condition_j = '(T.STRUCTURE_KIND_043a=\'1\' or T.STRUCTURE_KIND_043a=\'2\') '
    data_j = select_data(db_j,query_name_j, query_condition_j)
    p = np.zeros((10,10))
    for i in xrange(10):
        state_i = 9-i
        mask = np.logical_and(data_i[:,-1]==str(state_i),data_i[:,1]=='1')
        num_state_i = np.sum(mask)
        name_state_i = data_i[mask,0]
        for j in range(i, 10):
            state_j = 9-j
            mask = [(entry_name in name_state_i and entry_record=='1')
                    for entry_name, entry_record in zip(data_j[:,0],data_j[:,1])]
            cond_data_j = data_j[mask,:]
            num_state_j = np.sum(cond_data_j[:,-1]==str(state_j))
            #for entry in name_state_i:
                #occur = np.where(data_j[:,0]==entry)[0].size
                #if occur>=2:
                    #print '{} at location {}'.format(entry, np.where(data_j[:,0]==entry)[0])'' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '
            if num_state_i != 0:
                p[i,j] = float(num_state_j)/float(num_state_i)
            else:
                p[i,j] = -1

    return p

def transition_CS(year, component='deck'):
    if component == 'deck':
        query_name = 'T.STRUCTURE_NUMBER_008, T.DATE_OF_INSPECT_090, T.DECK_COND_058'
    elif 'super' in component:
        query_name = 'T.STRUCTURE_NUMBER_008, T.DATE_OF_INSPECT_090, T.SUPERSTRUCTURE_COND_059'
    elif 'sub' in component:
        query_name = 'T.STRUCTURE_NUMBER_008, T.DATE_OF_INSPECT_090, T.SUBSTRUCTURE_COND_060'
    # select specific bridge types
    query_condition = 'T.STRUCTURE_KIND_043a NOT IN (\'7\', \'8\', \'9\') '
    # record type = '1'
    query_condition = query_condition + 'AND T.RECORD_TYPE_005A=\'1\' '
    # data at year
    db0 = {'schema':'nbi{}'.format(year), 'table':'ca{}'.format(year)}
    data0 = select_data(db0,query_name, query_condition)
    # data at year+1
    db1 = {'schema':'nbi{}'.format(year+1), 'table':'ca{}'.format(year+1)}
    data1 = select_data(db1,query_name, query_condition)
    # data at year+2
    db2 = {'schema':'nbi{}'.format(year+2), 'table':'ca{}'.format(year+2)}
    data2 = select_data(db2,query_name, query_condition)
    def next_condition(data0, data1, data2):
        next_cs_array = []
        for name, inspect0, deck in zip(data0[:,0], data0[:,1], data0[:,2]):
            inspect1 = data1[data1[:,0]==name, 1]
            inspect2 = data2[data2[:,0]==name, 1]
            try:
                if inspect1  == inspect0 and inspect2 == inspect0:
                    next_cs = 'na'
                elif inspect1  != inspect0:
                    next_cs = data1[data1[:,0]==name,2][0]
                elif inspect2  != inspect0:
                    next_cs = data2[data2[:,0]==name,2][0]
                else: # for those cases that undergo rehabilitation at year2
                    next_cs = 'na'
            except IndexError:
                #print 'Cannot find bridge \'{}\''.format(name)
                next_cs = 'na'
            next_cs_array.append(next_cs)
        next_cs_array = np.array(next_cs_array)
        # get rid of missing bridges
        data0 = data0[next_cs_array!='na']
        next_cs_array = next_cs_array[next_cs_array!='na']

        return data0, next_cs_array

    data0, next_cs_array = next_condition(data0, data1, data2)
    data0[data0[:,2]>=8,2]=8
    data0[data0[:,2]<=1,2]=1
    next_cs_array[next_cs_array[:]>=8]=8
    next_cs_array[next_cs_array[:]<=1]=1
    p = np.zeros((8,8))
    for i in xrange(8):
        state_i = 8-i
        maski = data0[:,-1]==str(state_i)
        num_state_i = np.sum(maski)
        for j in range(i, 8):
            state_j = 8-j
            maskj = next_cs_array==str(state_j)
            num_state_j = np.sum(np.logical_and(maski, maskj))
            if num_state_i != 0:
                p[i,j] = float(num_state_j)/float(num_state_i)
            else:
                p[i,j] = -1

    return p

def transition_matrix(years, component='deck'):
    pmatrix_array = []
    for yr in years:
        pmatrix = transition_CS(yr, component=component)
        pmatrix_array.append(pmatrix)
    pmatrix_array = np.array(pmatrix_array)
    pmatrix_median = np.zeros((8,8))
    for i in xrange(8):
        for j in range(i, 8):
            pmatrix_median[i,j] = np.median(pmatrix_array[pmatrix_array[:,i,j]!=-1,i,j])
    for i in xrange(8):
        pmatrix_median[i,:] = pmatrix_median[i,:]/np.sum(pmatrix_median[i,:])

    return pmatrix_median

def bridge_correlation(bridge_db, corr_length):
    def int_to_degree(int_value):
        degree = np.floor(int_value/1e6)
        minute = np.floor((int_value - int(1e6)*degree)/1e4)
        second = int_value % int(1e4) / 100.
        return degree+minute/60.+second/3600.
    corr = np.ones((len(bridge_db),len(bridge_db)))
    for i_indx, bridge_i in enumerate(bridge_db):
        lati = int_to_degree(bridge_i[1])
        longi = int_to_degree(bridge_i[2])
        for j_indx in range(i_indx+1, len(bridge_db)):
            bridge_j = bridge_db[j_indx]
            latj = int_to_degree(bridge_j[1])
            longj = int_to_degree(bridge_j[2])
            #distance in km
            dis_ij = distance_on_unit_sphere(lati, longi, latj, longj)*6373.
            rho_ij = np.exp(-dis_ij**2/corr_length**2)
            corr[i_indx,j_indx] = rho_ij
            corr[j_indx,i_indx] = rho_ij

    return corr

if __name__ == '__main__':
    import itertools
    from mpldatacursor import datacursor
    plt.rc('font', family='serif', size=12)
    plt.rc('text', usetex=True)
    ## debug of get_deck_CS
    #data = get_deck_CS(2014).values()
    #plt.hist(data, bins=np.arange(1,10)-.5)
    #plt.xticks(np.arange(1,9))
    #plt.show()

    ## transition probability of deck_CS
    #step=1
    ##time_span = np.arange(1992,1996,step,dtype='int')
    #time_span = np.arange(1992,2013,step,dtype='int')
    #pmatrix_array = []
    #for yr in time_span:
        #pmatrix = transition_CS(yr, component='deck')
        #pmatrix_array.append(pmatrix)
    #pmatrix_array = np.array(pmatrix_array)
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
    pmatrix_deck = transition_matrix(np.arange(1992,2013,dtype='int'), component='deck')
    pmatrix_super = transition_matrix(np.arange(1992,2013,dtype='int'), component='superstructure')
    pmatrix_sub = transition_matrix(np.arange(1992,2013,dtype='int'), component='substructure')
    np.save('pmatrix.npy', {'deck': pmatrix_deck, 'super': pmatrix_super, 'sub': pmatrix_sub})

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
