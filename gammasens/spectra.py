"""norm_energy
Set of particle spectral functions for use when calculating
background and source rates
"""

import numpy as np
from math import pi
from astropy import units

def powerlaw( energy, index, norm, norm_energy=1.0 ):
    return norm*(energy/norm_energy)**(-index)

def exponential_cutoff(energy, cutoff_energy):
    return np.exp(-energy/cutoff_energy)
    

def lognormal( E, center, width ):
        return  1.0/(E*width*np.sqrt(2)*pi) * np.exp( -(np.log(E)-center)**2/(2.0*width**2))


def electron_spectrum(e_true_tev):
    """Cosmic-Ray Electron spectrum CTA version, with Fermi Shoulder, in
    units of :math:`\mathrm{TeV^{-1} s^{-1} m^{-2} sr^{-1}}`

    .. math::
       {dN \over dE dA dt d\Omega} = 

    """

    return units.Quantity(6.85e-5 * e_true_tev**-3.21 + \
        3.18e-3/(e_true_tev*0.776*np.sqrt(2*pi)) * \
        np.exp(-0.5 * (np.log(e_true_tev/0.107)/0.776)**2), 
                          "TeV**-1 s**-1 m**-2 sr**-1")

def electron_spectrum_powerlaw(e_true_tev):
    """ simple power-law cosmic-ray electron spectrum """
    norm = units.Quantity(6.85e-5, "1/(TeV s m**2 sr)") 
    return powerlaw( e_true_tev, index=3.21, norm=norm, norm_energy=1.0 )


def electron_spectrum_fermi(e_true_tev, E_peak_tev=0.107, width=0.776):
    """ cosmic ray electron spectrum including Fermi shoulder 
    units of :math:`\\mathrm{TeV^{-1} s^{-1} m^{-2} sr^{-1}}`
    """

    # including Fermi shoulder: powerlaw plus lognormal

    amp_peak = units.Quantity(3.186e-3,"1/(TeV s m**2 sr)")

    return electron_spectrum_powerlaw(e_true_tev) \
        + amp_peak*lognormal(e_true_tev, 
                             center= np.log(E_peak_tev), 
                             width=width)


def proton_spectrum(e_true_tev):
    """returns cosmic-ray proton spectrum 

    .. math::
        {dN \over dE dA dt d\Omega} = 0.096 \\left({E_t \over TeV}\\right)^{-2.7}
   

    in units of :math::`\\mathrm{TeV^{-1} s^{-1} m^{-2} sr^{-1}}`
    """

    norm = units.Quantity(0.096,"TeV**-1 s**-1 m**-2 sr**-1")
    return powerlaw(e_true_tev, norm=norm, index=2.7,norm_energy=1.0)

def cosmicray_spectrum(e_true_tev):
    """
    takes into account higher-Z CRs, not just protons

    Scale factor comes from Konrad (accounts for the fact that the
    trigger rate at the end is 10% higher than the proton rate due to
    the other species)
    """
    return proton_spectrum(e_true_tev) * 1.1

def hess_crab_spectrum(e_true_tev, fraction=1.0) :
    norm = fraction*units.Quantity(3.76e-11, "ct cm**-2 s**-1 TeV**-1")
    return  powerlaw( e_true_tev, norm=norm,
                      index=2.39, norm_energy=1.0) \
        * exponential_cutoff(e_true_tev, cutoff_energy=14.3)

def hess_binned_crab_spectrum(logEmin, logEmax, fraction=1.0):
    spec = np.zeros_like(logEmin)

    for ii,(logElo, logEhi) in enumerate(zip(logEmin,logEmax)):
        Elo = 10**logElo 
        Ehi = 10**logEhi 
        spec[ii] = integrate.quad( hess_crab_spectrum, Elo,Ehi,
                                   args=(fraction,) )[0]

    return units.Quantity(spec, "ct cm**-2/s**-1")

