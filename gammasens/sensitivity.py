"""
functions for calculating and plotting sensitivity curves from
ParticleDistributions.
"""

import numpy as np
from matplotlib import pyplot as plt
from math import pi
#import quantities as pq
from astropy import units
from scipy import optimize

from gammasens import inputs
from gammasens import spectra


EPSILON = 1.0e-10
PSTY = {'linestyle':'steps-mid', 'linewidth':2} # default line style    


def solid_angle( theta ):
    """
    Returns the solid angle for a circular region of angular radius
    `theta`
    """
    return (pi * theta**2).to( units.sr )

    # alternate (equivalant calculation):
    # return (1.0-np.cos(theta.to(units.rad)))*2.0*np.pi*units.steradian


@np.vectorize
def signif_lima(n_on,n_off,alpha):
    """
    Li and Ma significance formula 

    :param N_on: number of on-source counts
    :param N_off: nuber of off-source counts
    :param alpha: ratio of on/off exposure

    TODO: could speed this up by making it a true vector function, and not using
    vectorize (would need to change the logic a bit to use only numpy ufuncs)
    """
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


def residual_signif(N_on, N_off, alpha, minsig):
    """
    used my minimization routine to invert the Li and Ma formula
    """
    return minsig - signif_lima( N_on, N_off, alpha)


def excess(N_on, N_off, alpha):
    """ return excess counts given on, off, and exposure ratio """
    return (N_on - alpha*N_off)


def calc_background_rate(gammas, electrons, protons, return_all=False):
    """
    Calculate the background rate, given the input particle distributions

    :param gammas: gamma ray input data
    :type gammas: :class:`~gammasens.inputs.ParticleDistribution`
    :param electrons: electron input data
    :type electrons: :class:`~gammasens.inputs.ParticleDistribution`
    :param protons: proton input data
    :type protons: :class:`~gammasens.inputs.ParticleDistribution`
    """

    # nominal rates in HZ over the theta2 regions used by each
    # particle species note: does FOV have to be used for phi_diffuse
    # if it is much smaller?
    rp_nom = protons.rate_per_solidangle()*solid_angle(protons.phi_diffuse)
    re_nom = electrons.rate_per_solidangle()*solid_angle(electrons.phi_diffuse)
    # now want to normalize to the gamma-ray theta^2 cut (since the
    # protons and electrons were done at different cuts)
    rp = rp_nom * (gammas.thetasqr / (protons.thetasqr+EPSILON*units.deg))
    re = re_nom * (gammas.thetasqr / (electrons.thetasqr+EPSILON*units.deg))

    if return_all:
        return (re+rp), re,rp

    return (re+rp)


def calc_sensitivity_from_distributions( name, gammas, electrons, protons,  
                                         **kwargs ):
    """
    
    Arguments:
    - `name`:
    - `gammas`:
    - `electrons`:
    - `protons`:
    
    """
    
    background_rate = calc_background_rate( gammas, electrons,protons)
    gamma_aeff_reco = gammas.effective_area_reco()
    delta_e = gammas.deltaE

    return calc_sensitivity( name, background_rate, gamma_aeff_reco, 
                             delta_e, **kwargs)

def calc_sensitivity(name, background_rate, gamma_aeff_reco, delta_e,
                     obstime=5*units.h, 
                     num_bg_regions=2, min_signif=5.0, min_events=10.0, 
                     min_sys_pct=5.0, verbose=False):
    """
    calculates a sensitivity curve. Returns a dictionary of energy bin
    ranges, senstivitity measurement, and intermediate values useful
    for debugging and display

    :param name: identifier for plots
    :param background_rate: background rate array, 
                            as given by calc_background_rate()
    :param gamma_aeff_reco: reconstructed effective area of gammas , after all cuts (as output by gammas.effective_area_reco())
    :param delta_e: array of energy bins (not log), as output from gamma.deltaE
    :param obstime: observation time (as a units.Quantity)
    :param num_bg_regions: number of background regions to simulate
    :param min_signif: minimum significance per bin
    :param min_events: minimum number of excess counts per bin
    :param min_sys_pct: minimum background systematic percentage

    """

    # avoid constant warnings about quantity conversions 
    units.quantity.WARN_IMPLICIT_NUMERIC_CONVERSION.set(False)

    if verbose:
        print "CALCULATING SENSITIVITY FOR '{0}':".format(name)
        print "   num_bg_regions: ", num_bg_regions
        print "          object name extends type stime: ", \
            obstime.to(units.h),"(",obstime.to(units.min),")"
        print "       min_signif: ", min_signif
        print "       min_events: ", min_events
        print "      min_sys_pct: ", min_sys_pct,"%"

    
    # now calculate number of BG events in each energy bin:
    N_bg = background_rate * obstime.to(units.s)
    N_bg = N_bg.to(units.count).value

    alpha = np.ones(N_bg.shape) / num_bg_regions
    N_off = N_bg * num_bg_regions
    
    N_off[N_off<EPSILON] = np.nan

    # want to calcualte N_gamma = N_on-alpha*N_off such that:
    #  5sigma = signif_lima( N_on, N_off, alpha)
    #  - if the excess N_gamma = (N_on-alpha*N_off) < 10.0, excess = 10
    #  - if  N_gamma < N_off * 5.0% , then N_gamma = N_off*5%

    # first numerically solve for minimum significance
    #   note factor=0.1 in fsolve seems to give good results (default
    #   is too coarse)
    N_on = optimize.fsolve( residual_signif, np.zeros_like(N_off), 
                            args=(N_off,alpha,min_signif),
                            factor=0.1)

    # apply conditions on minimum excess
    N_on_orig = N_on.copy()
    mask = excess(N_on,N_off,alpha)<min_events
    if any(mask):
        N_on[mask] = min_events + alpha[mask]*N_off[mask]

    # apply conditions on minimum background systematics
    mask = excess(N_on,N_off,alpha)<N_off*min_sys_pct*0.01
    if any(mask):
        N_on[mask] = N_off[mask]*min_sys_pct*0.01 + alpha[mask]*N_off[mask]
           
    
    # calculate sensitivity limit, and conver to proper units 
    N_on[np.isnan(N_off)] = np.nan  #chop off bad values
    sens = excess(N_on,N_off,alpha)*units.count \
           /gamma_aeff_reco/obstime/delta_e

    sens = sens.to("ct cm**-2 s**-1 TeV**-1")

    # if verbose:
    #     print "#logE_lo  logE_hi  Sensitivity",sens.unit
    #     for el,eh,s in zip(gammas.log_e_lo, gammas.log_e_hi,sens.value):
    #         if np.isfinite(s):
    #             print "{0:8.3f} {1:8.3f} {2:10.3g}".format(el,eh,s)

    # return all the output (including intermediate values) in a dict,
    # for later plotting
    return dict( params=dict(obstime=obstime,
                              num_bg_regions = num_bg_regions,
                              min_signif=min_signif,
                              min_events=min_events,
                              min_sys_pct=min_sys_pct),
                 name=name,
                 sensitivity=sens,
                 N_on = N_on,
                 N_on_orig = N_on_orig,
                 N_off = N_off,
                 alpha = alpha )

def plot_effareas( gammas, electrons, protons ):

    plt.semilogy( protons.log_e, protons.effective_area_reco(), 
                  label="p", **PSTY)
    plt.semilogy( electrons.log_e, electrons.effective_area_reco(), 
                  label="e-", **PSTY )
    plt.semilogy( gammas.log_e, gammas.effective_area_reco(), 
                  label=r"$\gamma$", **PSTY)
    plt.ylabel("Effective Area (m)")
    plt.xlabel("log E$_{reco}$/TeV")
    plt.legend(loc="best")


def plot_count_distributions( log_e, sens ):

    plt.semilogy( log_e, sens['N_on'], label="N_on", **PSTY  )
    plt.semilogy( log_e, sens['N_off'], color='r', label="N_off", **PSTY )
    plt.semilogy( log_e, excess(sens['N_on'],sens['N_off'],sens['alpha']), 
                  color='black', label="N_exc", **PSTY )
    plt.ylabel("Counts ({0})".format(sens['params']['obstime']))
    plt.xlabel("Log10(E/TeV)")
    plt.legend(loc='best')


def plot_significances( log_e, sens ):
    """ sens: output dictionary from calc_sensitivity """ 

    plt.scatter( log_e, signif_lima( sens['N_on'], sens['N_off'], 
                                    sens['alpha'] ) )
    plt.scatter( log_e, signif_lima( sens['N_on_orig'], sens['N_off'], 
                                    sens['alpha'] ),
                 color='grey' )
    plt.ylabel("Significance")
    plt.xlabel("Log10(E/TeV)")


def plot_rates( log_e, rate_p, rate_e, sens ):
    plt.semilogy( log_e, rate_p,label="p ", **PSTY)
    plt.semilogy( log_e, rate_e,label="e- ",**PSTY)

    # also overlay the predicted minimum gamma excess rate
    excess_rate = excess(sens['N_on'],sens['N_off'],sens['alpha'])*units.ct \
                  / sens['params']['obstime']

    plt.semilogy( log_e, excess_rate.to(units.ct/units.s),
                  label=r"$\gamma_\mathrm{exc}$",color='g', 
                  drawstyle='steps-mid', linestyle='--'  )


    plt.ylabel("Rate ({0})".format(rate_p.unit.to_string()) )
    plt.xlabel("log E$_{reco}$/TeV")
    plt.legend( loc="best")

def plot_sensitivity(log_e, sens, **kwargs):
    """
    Display the differential sensitivity curve 

    :param sens: sensitivity output dictrionary
    """

    
    sensitivity = sens['sensitivity']
    par = sens['params']

    label=r"{2} {0}, {1} $\sigma$".format(par['obstime'].to(units.h), 
                                             par['min_signif'], sens['name'])

    lines = plt.semilogy( log_e, sensitivity.value, marker="+",
                          label=label,**kwargs )
    plt.ylabel("Sens ({0})".format( sensitivity.unit.to_string() ))
    plt.xlabel("log E$_{reco}$/TeV")
    plt.title("Differential Sensitivity")
    plt.ylim(1e-14, 1e-5)
    return lines

def plot_crab( log_e, sens ):
    """ call after plot_sensitivity() to overplot Crab flux 
    :param sens: output dictionary from calc_sensitivity()
    """

    plt.semilogy( log_e, 
                  spectra.hess_crab_spectrum( 10**(log_e) ), 
                  linestyle="--", color='black',
                  label="100% Crab")
    
    plt.semilogy( log_e, 
                  spectra.hess_crab_spectrum( 10**(log_e),fraction=0.1 ), 
                  linestyle="--", color='gray',
                  label="10% Crab" )
    
    plt.semilogy( log_e, 
                  spectra.hess_crab_spectrum( 10**(log_e),fraction=0.01 ), 
                  linestyle=":", color='gray',
                  label=" 1% Crab" )
    



def plot_sensitivity_crabunits( log_e, sens ):
    """
    display the sensitivity curve as a fraction of the crab flux.
    
    The crab model currently used is only valid in the VHE energy
    range (100 GeV-100TeV)

    :param sens: output dictionary from calc_sensitivity()
    """

    
    sensitivity = sens['sensitivity']
    par = sens['params']

    label=r"{2} {0}, {1} $\sigma$".format(par['obstime'].to(units.h), 
                                             par['min_signif'], sens['name'])
    crabs = spectra.hess_crab_spectrum( 10**log_e )
    plt.plot( log_e, sensitivity/crabs, label=label )
    plt.xlabel("log10(E/TeV)")
    plt.ylabel("Diff Sensitivity (Crab Units)")

if __name__ == '__main__':

    gammas, electrons, protons = inputs.loadAllFromFITS( "test-*.fits" )
    sens = calc_sensitivity( gammas, electrons, protons )


