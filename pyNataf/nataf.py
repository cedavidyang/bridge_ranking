import numpy as np
import scipy.stats as stats
import scipy.optimize as op
# from pyDUE.util import distance_on_unit_sphere

def nataf_transform(rho1, x1rv, x2rv=None, xtol=1e-6, rtol=1e-5, nsmp=1e4, nconv=1e3, maxiter=100):
    """ rho1: coefficient of correlation of standard normal dist
        x1rv: random variable x1
        x2rv: random variable x2
        return: rho0, coefficient of correlation of x1 and x2 corresponding to rho1
    """
    cov = 1.
    mean = 1.
    int_smp=[]
    rv = stats.multivariate_normal(mean=[0,0], cov=[[1.,rho1],[rho1,1.]], allow_singular=True)
    if x2rv is None:
        x2rv = x1rv
    m1,v1 = x1rv.stats()
    m2,v2 = x2rv.stats()
    i = 0
    while abs(cov) > rtol and abs(mean)>xtol and i<maxiter:
        x1 = x1rv.rvs(size=int(nsmp))
        x2 = x2rv.rvs(size=int(nsmp))
        xcdf = np.vstack((x1rv.cdf(x1), x2rv.cdf(x2))).T
        int_smp0 = (x1-m1)*(x2-m2)*rv.pdf(stats.norm.ppf(xcdf))/\
        (stats.norm.pdf(stats.norm.ppf(x1rv.cdf(x1)))*stats.norm.pdf(stats.norm.ppf(x2rv.cdf(x2))))
        #if xrv.dist.name == 'uniform':
            #m,v = xrv.stats()
            #int_smp0 = (x1-m)*(x2-m)*rv.pdf(stats.norm.ppf(x))/\
            #(stats.norm.pdf(stats.norm.ppf(x1))*stats.norm.pdf(stats.norm.ppf(x2)))
        #elif xrv.dist.name == 'lognorm':
            #m,v = xrv.stats(); std = np.sqrt(v)
            #scale = m/np.sqrt(1+(std/m)**2) # log(scale) = mean of LogX
            #s=np.sqrt(np.log(1+(std/m)**2))
            #int_smp0 = (x1-1.)*(x2-1.)*rv.pdf(1./s*np.log(x/scale))/\
            #(stats.norm.pdf(1./s*np.log(x1/scale))*stats.norm.pdf(1./s*np.log(x2/scale)))
        int_smp = np.append(int_smp, int_smp0)
        means = np.cumsum(int_smp) / (np.arange(int_smp.shape[0])+1)
        mean = np.mean(int_smp)
        cov = np.std(means[-int(nconv):]) / np.mean(means[-int(nconv):])
        i += 1
    rho0 = 1/np.sqrt(v1*v2)*mean

    return rho0

def nataf_transform_inverse(rho0, xtol=1e-5, rtol=1e-4):
    def find_root(rho1):
        if rho1>1.-xtol: rho1=1.-xtol
        elif rho1<-1.+xtol: rho1=-1.+xtol
        return abs(nataf_integral(rho1, xtol=xtol, rtol=rtol)-1./12.*rho0)
    #rho1 = op.brentq(find_root,-1.0, 1.0, rtol=rtol)[0]
    #rho1 = op.newton(find_root, rho0, tol=xtol)
    rho1 = op.minimize(find_root, rho0, bounds=[(-1.+xtol, 1.-xtol)], tol=xtol,
            options={'disp':True}).x[0]
    #rho1 = op.minimize(find_root, 0, method='Nelder-Mead', tol=xtol,
            #options={'disp':True}).x[0]
    #rho1 = op.brute(find_root, ((-1., 1.),), Ns=20, disp=True)
    return rho1


def natafcurve(x, a, b, c):
    rho1 = x*(a*(x-1.)**3+b*(x-1.)**2+c*(x-1.)+1.)
    return rho1

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.ion()
    plt.rc('font', family='serif', size=12)
    plt.rc('text', usetex=True)
    np.random.seed(1)

    ## debug of nataf_transform
    #m=1; std=0.4
    #scale = m/np.sqrt(1+(std/m)**2) # log(scale) = mean of LogX
    #s=np.sqrt(np.log(1+(std/m)**2))
    #xrv = stats.lognorm(s=s, scale=scale)
    #rho1_array = np.arange(-.99, 1.0, 0.01)
    #rho0_array = []
    #for rho1 in rho1_array:
        #rho0_array.append(nataf_transform(rho1, xrv, xtol=1e-5, rtol=1e-4))
    #rho3_array = np.arange(-1., 1.01, 0.01)
    ## analytical solution for lognormal
    #rho4_array = np.log(1+rho3_array*s**2)/np.sqrt(np.log(1+s**2)*np.log(1+s**2))
    #plt.plot(rho0_array, rho1_array, 'b-')
    #plt.plot(rho3_array, rho4_array,'g-')
    #plt.plot([-1.,1.], [-1.,1.], 'r-.')

    # uniform(0,1) to standard normal
    xrv = stats.uniform()
    rho1_array = np.arange(0., 1.0, 0.01)
    rho0_array = []
    for rho1 in rho1_array:
        rho0_array.append(nataf_transform(rho1, xrv))
    rho0_array = np.asarray(rho0_array)
    plt.plot(rho0_array, rho1_array, 'b-')
    plt.plot([0.,1.], [0.,1.], 'r-.')
    plt.grid()
    plt.xlabel('Original correlation $\\rho$')
    plt.ylabel('Adjusted correlation $\\rho\'$')
    ## least square
    #xdata = rho0_array[:-1]
    #ydata = rho1_array[:-1]
    #popt, pcov = op.curve_fit(natafcurve, xdata, ydata)
    #plt.plot(xdata, natafcurve(xdata,*popt))
    #np.save('nataf_popt.npy', popt)
