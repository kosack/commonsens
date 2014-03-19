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

from commonsens import inputs
from commonsens import spectra
from commonsens import stats
from commonsens import config


EPSILON = 1.0e-10
PSTY = {'linestyle':'steps-mid', 'linewidth':2} # default line style    

class SensOutput(dict):
    """
    a dict that allows tab completion in ipython, and access to
    data members using instance.key=value notation
    """
    def __init__(self,**kw):
        dict.__init__(self,kw)
        self.__dict__ = self


def xlabel_energy():
    """ add the x-label for Energy in nice format """
    plt.xlabel(r"$E_{reco} (\mathrm{TeV})$")


def solid_angle( theta ):
    """
    Returns the solid angle for a circular region of angular radius
    `theta`
    """
    return (pi * theta**2).to( units.sr )

    # alternate (equivalant calculation):
    # return (1.0-np.cos(theta.to(units.rad)))*2.0*np.pi*units.steradian


def fill_between(x, y1, y2=0, ax=None, **kwargs):
    """Plot filled region between `y1` and `y2`.

    This function works exactly the same as matplotlib's fill_between, except
    that it also plots a proxy artist (specifically, a rectangle of 0 size)
    so that it can be added it appears on a legend.
    """
    ax = ax if ax is not None else plt.gca()
    ax.fill_between(x, y1, y2, **kwargs)
    p = plt.Rectangle((0, 0), 0, 0, **kwargs)
    ax.add_patch(p)
    return p



def residual_signif(N_on, N_off, alpha, minsig):
    """
    used my minimization routine to invert the Li and Ma formula
    """
    return minsig - stats.signif_lima( N_on, N_off, alpha)




def calc_background_rate(gammas, electrons, protons, return_all=False):
    """
    Calculate the background rate, given the input particle distributions

    :param gammas: gamma ray input data
    :type gammas: :class:`~commonsens.inputs.ParticleDistribution`
    :param electrons: electron input data
    :type electrons: :class:`~commonsens.inputs.ParticleDistribution`
    :param protons: proton input data
    :type protons: :class:`~commonsens.inputs.ParticleDistribution`
    """

    # nominal rates in Hz over the theta2 regions used by each
    # particle species 

    # TODO: does FOV have to be used for phi_diffuse
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


def calc_from_distributions( name, gammas, electrons, protons,  
                                         **kwargs ):
    """Calculate sensitivity from input ParticleDistributions for gammas,
    electrons, and protons.

    Arguments:
    :param `name`:
    :param gammas,electrons,protons: particle input data
    :type gammas,electrons,protons: :class:`~commonsens.inputs.ParticleDistribution`
    :param kwargs:: see commonsens.sensitivity.calc_sensitivity()

    """
    
    background_rate = calc_background_rate( gammas, electrons,protons)
    gamma_aeff_reco = gammas.effective_area_reco()
    delta_e = gammas.delta_e

    return calc_sensitivity( name, background_rate, gamma_aeff_reco, 
                             delta_e, **kwargs)

def calc_sensitivity(name, background_rate, gamma_aeff_reco, delta_e,
                     obstime=5*units.h, 
                     num_bg_regions=2, min_signif=5.0, min_events=10.0, 
                     min_sys_pct=5.0):
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


    isintegral = False
    try:
        test = background_rate.to("s**-1 * TeV")
        isintegral = True
    except units.UnitsError:
        isintegral=False


    if config.verbose:
        print "CALCULATING SENSITIVITY FOR '{0}':".format(name)
        print "   num_bg_regions: ", num_bg_regions
        print "             time: ", \
            obstime.to(units.h),"(",obstime.to(units.min),")"
        print "       min_signif: ", min_signif
        print "       min_events: ", min_events
        print "      min_sys_pct: ", min_sys_pct,"%"
        print "         integral: ", isintegral
    

    # clean the histograms: remove anything below X% of the peak in
    # effective area:
    gamma_aeff_reco = gamma_aeff_reco.copy()
    background_rate = background_rate.copy()
    max_aeff = config.effective_area_fraction_min*np.max(gamma_aeff_reco[np.isfinite(gamma_aeff_reco)])
    badbins = gamma_aeff_reco < max_aeff

    gamma_aeff_reco[badbins] = 0
    background_rate[badbins] = 0

    # now calculate number of BG events in each energy bin:
    N_bg = background_rate * obstime.to(units.s)
    N_bg_unit = N_bg.unit
    N_bg = N_bg.value

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
    mask = stats.excess(N_on,N_off,alpha)<min_events
    if any(mask):
        N_on[mask] = min_events + alpha[mask]*N_off[mask]

    # apply conditions on minimum background systematics
    mask = stats.excess(N_on,N_off,alpha)<N_off*min_sys_pct*0.01
    if any(mask):
        N_on[mask] = N_off[mask]*min_sys_pct*0.01 + alpha[mask]*N_off[mask]
    
    # calculate sensitivity limit, and conver to proper units 
    N_on[np.isnan(N_off)] = np.nan  #chop off bad values
    sens = stats.excess(N_on,N_off,alpha)*N_bg_unit \
           /gamma_aeff_reco/obstime/delta_e

#    sens = sens.to("cm**-2 s**-1 TeV**-1")

    # if verbose:
    #     print "#logE_lo  logE_hi  Sensitivity",sens.unit
    #     for el,eh,s in zip(gammas.log_e_lo, gammas.log_e_hi,sens.value):
    #         if np.isfinite(s):
    #             print "{0:8.3f} {1:8.3f} {2:10.3g}".format(el,eh,s)

    # return all the output (including intermediate values) in a dict,
    # for later plotting
    return SensOutput( params=dict(obstime=obstime,
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

def calc_integral_sensitivity(name, background_rate, 
                              gamma_aeff_reco, delta_e,
                              obstime=5*units.h, 
                              num_bg_regions=2, min_signif=5.0, min_events=10.0, 
                              min_sys_pct=5.0):
    """
    not yet working
    """
    
    bgrate_e = (background_rate).value
    int_bgrate = np.cumsum(bgrate_e[::-1])[::-1]*1/units.s

    aeff_e =(gamma_aeff_reco).value
    int_aeff = np.cumsum( aeff_e[::-1] )[::-1]*units.erg*units.cm**2

    de = delta_e.to(units.TeV).value
    int_delta_e = np.cumsum( de[::-1] )[::-1]*units.TeV


    return calc_sensitivity( name, int_bgrate, int_aeff, int_delta_e,
                             obstime=obstime, num_bg_regions=num_bg_regions,
                             min_signif=min_signif, min_events=min_events,
                             min_sys_pct=min_sys_pct )


def plot_effareas( gammas, electrons, protons ):

    aeff_e = protons.effective_area_reco().to(units.m**2).value
    aeff_p = electrons.effective_area_reco().to(units.m**2).value
    aeff_g = gammas.effective_area_reco().to(units.m**2).value

    plt.loglog( 10**protons.log_e, aeff_p, label="p", **PSTY)
    plt.loglog( 10**electrons.log_e, aeff_e, label="e-", **PSTY )
    plt.loglog( 10**gammas.log_e, aeff_g, label=r"$\gamma$", **PSTY)
    plt.ylabel("Effective Area (m$^2$)")
    xlabel_energy()
    plt.legend(loc="best")


def plot_count_distributions( log_e, sens ):

    plt.loglog( 10**log_e, sens['N_on'], label="N_on", **PSTY  )
    plt.loglog( 10**log_e, sens['N_off'], color='r', label="N_off", **PSTY )
    plt.loglog( 10**log_e, stats.excess(sens['N_on'],sens['N_off'],
                                        sens['alpha']), 
                  color='black', label="N_exc", **PSTY )
    plt.ylabel("Counts ({0})".format(sens['params']['obstime']))
    xlabel_energy()
    plt.legend(loc='best')


def plot_significances( log_e, sens ):
    """ sens: output dictionary from calc_sensitivity """ 

    plt.semilogx()
    plt.scatter( 10**log_e, stats.signif_lima( sens['N_on'], sens['N_off'], 
                                           sens['alpha'] ),label=r"adjusted" )
    plt.scatter( 10**log_e, stats.signif_lima( sens['N_on_orig'], sens['N_off'], 
                                           sens['alpha']),label=r"original" ,
                 color='grey' )
    plt.ylabel("Significance")
    xlabel_energy()
    plt.legend(loc='best')


def plot_rates( log_e, rate_p, rate_e, sens ):
    plt.loglog( 10**log_e, rate_p.to("s^-1").value,label="p ", **PSTY)
    plt.loglog( 10**log_e, rate_e.to("s^-1").value ,label="e- ",**PSTY)

    # also overlay the predicted minimum gamma excess rate
    excess_rate = stats.excess(sens['N_on'],sens['N_off'],
                               sens['alpha']) \
                  / sens['params']['obstime']

    plt.loglog( 10**log_e, excess_rate.to(1/units.s).value,
                  label=r"$\gamma_\mathrm{exc}$",color='g', 
                  drawstyle='steps-mid', linestyle='--'  )


    plt.ylabel("Rate {0}".format(rate_p.unit.to_string(format='latex')) )
    xlabel_energy()
    plt.legend( loc="best")

def plot_sensitivity(log_e, sens, esquared=False, shade_percent=None, **kwargs):
    """
    Display the differential sensitivity curve 

    :param sens: sensitivity output dictionary
    """

    E = 10**log_e * units.TeV
    
    sensitivity = sens['sensitivity'].to(1/units.cm**2/units.s/units.TeV)

    if (esquared):
        sensitivity *= E**2
        sensitivity = sensitivity.to(units.erg/units.cm**2/units.s)

    par = sens['params']
    label=r"{2} {0}, {1} $\sigma$".format(par['obstime'].to(units.h), 
                                          par['min_signif'], sens['name'])

    plt.loglog()

    lines = None
    if shade_percent is not None:
        # do a shaded region:
        fill_between( E.value, 
                      sensitivity.value-sensitivity.value*shade_percent, 
                      sensitivity.value+sensitivity.value*shade_percent, 
                      label=label, **kwargs )
    else :
        lines = plt.loglog( E.value, 
                            sensitivity.value,
                            marker=None,
                            label=label, drawstyle="default",**kwargs )



    plt.ylabel("Sens {0}".format( sensitivity.unit.to_string(format='latex') ))
    xlabel_energy()

    plt.title("Differential Sensitivity")
    plt.ylim(1e-14, 1e-5)
    return lines

def plot_crab( log_e, esquared=False ):
    """ call after plot_sensitivity() to overplot Crab flux 
    """


    crab = spectra.hess_crab_spectrum( 10**(log_e) )
    E = 10**log_e * units.TeV

    if (esquared):
        crab = (crab*E**2).to(units.erg/units.cm**2/units.s) 



    plt.semilogy( E.value, 
                  crab.value * 1.0, 
                  linestyle="--", color='black',
                  label="100% Crab")
    
    plt.semilogy( E.value, 
                  crab.value * 0.1, 
                  linestyle="--", color='gray',
                  label="10% Crab" )
    
    plt.semilogy( E.value, 
                  crab.value * 0.01, 
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
    plt.loglog()
    plt.plot( 10**log_e, (sensitivity/crabs).to("").value, label=label )
    xlabel_energy()
    plt.ylabel("Diff Sensitivity (Crab Units)")

if __name__ == '__main__':

    gammas, electrons, protons = inputs.load_all_from_fits( "test-*.fits" )
    sens = calc_sensitivity( gammas, electrons, protons )


