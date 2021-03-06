import numpy as np

def semidefinitive(correlation, weights=None,tol=1e-8, maxiter=1e4, deftol=1e-6):
    """ corrUpdate update sample correlation matrix according to algorithm 3.3 in Higham (2002)
    """
    nV = correlation.shape[0]
    if weights is None:
        weights = np.ones(nV)
    wMatrix = np.diag(weights)

    dS = 0
    newCorrY = np.copy(correlation)
    newCorr = np.eye(nV)

    dXY = np.linalg.norm(newCorr-newCorrY)

    i=0
    while dXY>tol and i<maxiter:
        r = newCorrY - dS
        newCorr = project2s(r, wMatrix, deftol=deftol)
        dS = newCorr - r
        newCorrY = project2u(newCorr, wMatrix)
        dXY = np.linalg.norm(newCorr-newCorrY)
        i += 1

    ## diag=1 and symetric
    #for i in xrange(nV):
        #newCorr[i,i] = 1.
    for i in xrange(nV):
        for j in range(i+1,nV):
            tmp = 0.5*(newCorr[i,j]+newCorr[j,i])
            newCorr[i,j] = tmp
            newCorr[j,i] = tmp

    return newCorr

def project2u(smpCorr, wMatrix):

    theta = np.linalg.solve(np.linalg.inv(wMatrix)*np.linalg.inv(wMatrix),
        np.diag(smpCorr - np.eye(smpCorr.shape[0])) )
    pUcorr = smpCorr - np.dot(np.dot(np.linalg.inv(wMatrix), np.diag(theta)), np.linalg.inv(wMatrix))

    return pUcorr

def project2s(smpCorr, wMatrix, deftol=0.):

    wMatrixSqrt = np.sqrt(wMatrix)
    pSpos = get_corrPos( np.dot(np.dot(wMatrixSqrt,smpCorr),wMatrixSqrt),deftol=deftol)
    pScorr = np.dot(np.dot(np.linalg.inv(wMatrixSqrt),pSpos),np.linalg.inv(wMatrixSqrt))

    return pScorr

def get_corrPos(smpCorr, deftol=0.):

    eigValuePos, eigVector = np.linalg.eig(smpCorr)
    eigValuePos[ eigValuePos<deftol ] = deftol
    eigValueDiagPos = np.diag(eigValuePos)

    corrPos = np.dot(np.dot(eigVector,eigValueDiagPos),eigVector.T)

    return corrPos

if __name__ == '__main__':
    #nV = 4
    #corr = np.ones((nV,nV))*-1.
    #for i in xrange(nV):
        #corr[i,i] = 1.
    #corr1 = semidefinitive(corr,deftol=1e-6)
    #print corr1
    corr = np.load('correlation_after_nataf.npy')
    corr1 = semidefinitive(corr,deftol=1e-6)
