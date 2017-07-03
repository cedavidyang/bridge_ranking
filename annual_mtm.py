import numpy as np
from scipy.optimize import brentq

def mtm_2yearto1year(pmatrix):
    nstate, nstate = np.shape(pmatrix)
    amatrix = np.zeros((nstate,nstate))
    for i in np.arange(nstate-1, -1, -1):
        for j in np.arange(i, nstate):
            if i == j:
                amatrix[i,j] = np.sqrt(pmatrix[i,j])
            else:
                def sum2year(aij, amatrix=amatrix):
                    amatrix[i,j] = aij
                    pij = 0.
                    for k in np.arange(i,j+1):
                        pij += amatrix[i,k]*amatrix[k,j]
                    return pij
                f = lambda x: sum2year(x) - pmatrix[i,j]
                aij = brentq(f, -1., 2.)
                amatrix[i,j] = aij
    amatrix[amatrix<0] = 0.
    amatrix[amatrix>1] = 1.
    return amatrix

if __name__ == '__main__':
    pdict = np.load('pmatrix.npy')[()]
    psuper = pdict['super']
    asuper = mtm_2yearto1year(psuper)
    psub = pdict['sub']
    asub = mtm_2yearto1year(psub)
