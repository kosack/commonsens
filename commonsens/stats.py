"""
Some functions for calculating statistical significances,
excesses, etc from ON/OFF data.

"""

__all__ = ['signif_lima', 'signif_lima_scalar', 'excess']


import numpy as np
from math import pi

EPSILON = 1.0e-10

def signif_lima_scalar(n_on,n_off,alpha):
    """Li and Ma significance formula, with scalar inputs.
    use the version called signif_lima(), which is vectorized

    :param N_on: number of on-source counts
    :param N_off: nuber of off-source counts
    :param alpha: ratio of on/off exposure

    TODO: could speed this up by making it a true vector function, and
    not using vectorize (would need to change the logic a bit to use
    only numpy ufuncs)

    """
    global EPSILON
    thesum = n_on + n_off + EPSILON
    arg1 = (1.0+alpha)/alpha * ( n_on/thesum)
    arg2 = (1.0+alpha)*(n_off/thesum)

    if (alpha*n_off > n_on):
        return (n_on-alpha*n_off)/np.sqrt(alpha*thesum)
    
    
    if arg1 > 1e-10 and arg2 > 1e-10:
        sigsqr = np.fabs(2.0*(n_on*np.log(arg1) + n_off*np.log(arg2)))
    else:
        if  thesum > 1.0e-5 :
            return (n_on-alpha*n_off)/np.sqrt(alpha*thesum)
        else :
            return 0.0
    if ((n_on - alpha*n_off) > 0.0):
        return np.sqrt(sigsqr)
    else:
        return -np.sqrt(sigsqr)


signif_lima = np.vectorize(signif_lima_scalar, doc=signif_lima_scalar.__doc__)


def excess(N_on, N_off, alpha):
    """ return excess counts given on, off, and exposure ratio """
    return (N_on - alpha*N_off)
