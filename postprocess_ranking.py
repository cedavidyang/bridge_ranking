import os
import sys
import psycopg2

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import shelve


def retrieve_bridge_db(bridge_names, cur_nbi):
    bridge_names = np.array([name.lstrip() for name in bridge_name])
    # get bridge condition ratings
    s='select cast(structure_number_008 as varchar(15)), \
    cast(lat_016 as int), cast(long_017 as int),\
    cast(structure_len_mt_049 as int), cast(deck_width_mt_052 as int),\
    cast(deck_cond_058 as int),\
    cast(superstructure_cond_059 as int), cast(substructure_cond_060 as int), \
    cast(detour_kilos_019 as real), cast(sufficiency_rating as int) \
    from nbi2014.ca2014 where ltrim(rtrim(structure_number_008)) in {} and \
    record_type_005a=\'1\';'.format(tuple(bridge_names))
    cur_nbi.execute(s)
    bridge_data = cur_nbi.fetchall()
    bridge_db = []
    for bridge_entry in bridge_data:
        tmp = list(bridge_entry)
        bridge_db.append(tmp)

    return np.asarray(bridge_db, dtype=object)

def safety_rating(rating):
    if rating <=2: index = 55
    elif rating == 3: index = 40
    elif rating == 4: index = 25
    elif rating == 5: index = 10
    elif rating >= 6: index = 0
    safety = (55.-index)/55.
    return safety

if __name__ == '__main__':
    plt.rc('font', family='serif', size=12)

    filename = os.path.join(os.path.abspath('./'), 'figures', 'ranking_LA 2015-07-13 09-02-58.407000', 'data_shelve.out')
    my_shelf = shelve.open(filename, 'r')
    globals()['nsmp']=my_shelf['nsmp']
    globals()['bridge_risk_data']=my_shelf['bridge_risk_data']
    globals()['bridge_db'] = my_shelf['bridge_db']
    my_shelf.close()

    filename = os.path.join(os.path.abspath('./'), 'figures', 'ranking_LA 2015-07-14 09-41-38.529000', 'data_shelve.out')
    my_shelf = shelve.open(filename, 'r')
    nsmp += my_shelf['nsmp']
    bridge_risk_data = np.vstack((bridge_risk_data, my_shelf['bridge_risk_data']))
    my_shelf.close()

    mean_risk = np.mean(bridge_risk_data, axis=0)
    bridge_name = np.asarray(bridge_db, dtype=object)[:,0].astype(str)
    rank1 = bridge_name[mean_risk.argsort()][-10:]


    # open databases
    conn_nbi = psycopg2.connect("dbname='nbi' user='amadeus' host='localhost' password='19881229'")
    cur_nbi = conn_nbi.cursor()

    bridge_db = retrieve_bridge_db(bridge_name, cur_nbi)
    nbridges = bridge_db.shape[0]
    safety= np.zeros(nbridges)
    sufficiency= np.zeros(nbridges)
    for indx, (name, lat, long, length, width, deck_cs0, super_cs0, sub_cs0, detour, rating) in enumerate(bridge_db):
        safety[indx] = safety_rating(np.minimum(super_cs0, sub_cs0))
        sufficiency[indx] = rating


    rank2 = bridge_name[safety.argsort()][-10:]
    rank3 = bridge_name[sufficiency.argsort()][-10:]
