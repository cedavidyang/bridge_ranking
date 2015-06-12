import numpy as np
import scipy.stats as stats

FREESPEED = 30.
FFDELAY0=8e3/FREESPEED
CAP0 = 4500.
ALPHA = 0.15
BETA = 4
CAP_DROP = 0.1

def generate_cost_smp(nsmp, link_flow, detour_list, cap_drop_list, reliability_list):
    link_cost_smp = []
    for ismp in xrange(nsmp):
        fail_smp = stats.bernoulli.rvs(stats.norm.cdf(-reliability_list))
        ffdelay = FFDELAY0 + np.sum(detour_list[fail_smp.astype(bool)])/FREESPEED
        try:
            cap = (1. - np.max(cap_drop_list[fail_smp.astype(bool)]))*CAP0
        except ValueError: # all safe, max function return ValueError
            cap = 1.*CAP0
        link_cost = ffdelay*( 1.+ALPHA*(link_flow/cap)**BETA )
        link_cost_smp.append(link_cost)
    return link_cost_smp

def main(link_flow_array):
    size = 20
    nsmp = int(1e5)
    np.random.seed(1)
    detour_list = np.random.randint(3, 5, size=size)*1e3
    cap_drop_list = CAP_DROP*np.ones(size)
    lower = 0.0
    upper = 3.0
    reliability_list = np.random.rand(size)*(upper-lower)+lower
    link_cost_data = []
    for link_flow in link_flow_array:
        link_cost_smp = generate_cost_smp(nsmp, link_flow, detour_list, cap_drop_list, reliability_list)
        link_cost_data.append(link_cost_smp)

    return np.asarray(link_cost_data)

if __name__ == '__main__':
    link_flow_array = np.linspace(0, 9000, num=20)
    #link_cost_data = main(link_flow_array)
    link_cost_data = np.load('link_cost_data.npy')
    import matplotlib.pyplot as plt
    plt.ion()
    fig, ax = plt.subplots(1,1)
    ax.plot(link_flow_array, np.mean(link_cost_data, axis=1), 'b-')
    ax.plot(link_flow_array, FFDELAY0*(1+ALPHA*(link_flow_array/CAP0)**BETA), 'g-')
    data_min = np.min(link_cost_data)
    data_max = np.max(link_cost_data)
    ax.set_ylim((data_min, data_max))
    # draw vertical pdf
    xtick_max = ax.get_xticks()[-1]
    xtick_min = ax.get_xticks()[0]
    #ytick_max = ax.get_yticks()[-1]
    #ytick_min = ax.get_yticks()[0]
    left = fig.subplotpars.left
    right = fig.subplotpars.right
    top = fig.subplotpars.top
    bottom = fig.subplotpars.bottom
    def xdata_location(xdata):
        loc = (right-left)/(xtick_max-xtick_min)*(xdata-xtick_min)+left
        return loc
    for flow in link_flow_array[1:-1:2]:
        ax.axvline(x=flow, linestyle='--')
        loc = xdata_location(flow)
        rect = [loc, bottom, 0.1*(right-left), top-bottom]
        axins = fig.add_axes(rect)
        indx = np.where(link_flow_array==flow)[0][0]
        axins.hist(link_cost_data[indx, :], orientation='horizontal', bins=10, histtype='step', color='r')
        axins.set_ylim((data_min, data_max))
        axins.set_axis_off()
    plt.show()
