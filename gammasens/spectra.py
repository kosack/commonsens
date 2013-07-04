"""
Set of particle spectral functions for use when calculating
background and source rates
"""

import numpy as np
from math import pi
from astropy import units

def electron_spectrum(Etrue_tev):
    """Cosmic-Ray Electron spectrum CTA version, with Fermi Shoulder, in
    units of :math:`\mathrm{TeV^{-1} s^{-1} m^{-2} sr^{-1}}`

    :param Etrue: true energy of electron in TeV 
    """

    return units.Quantity(6.85e-5 * Etrue_tev**-3.21 + \
        3.18e-3/(Etrue_tev*0.776*np.sqrt(2*pi)) * \
        np.exp(-0.5 * (np.log(Etrue_tev/0.107)/0.776)**2), 
                          "TeV**-1 s**-1 m**-2 sr**-1")

def proton_spectrum(Etrue_tev):
    """returns cosmic-ray proton spectrum in units of :math:`\mathrm{TeV^{-1}
    s^{-1} m^{-2} sr^{-1}}`
   
    :param Etrue_tev: true proton energy in TeV
    """

    return units.Quantity(0.096 * Etrue_tev**-2.7,"TeV**-1 s**-1 m**-2 sr**-1")

def cosmicray_spectrum(Etrue_tev):
    """
    takes into account higher-Z CRs, not just protons

    Scale factor comes from Konrad (accounts for the fact that the
    trigger rate at the end is 10% higher than the proton rate due to
    the other species)
    """
    return proton_spectrum(Etrue_tev) * 1.1

def hess_crab_spectrum(Etrue_tev, fraction=1.0) :
    norm = fraction*units.Quantity(3.76e-11, "ct cm**-2 s**-1 TeV**-1")
    return  norm * Etrue_tev**-2.39 * np.exp(-Etrue_tev/14.3)

def hess_binned_crab_spectrum(logEmin, logEmax, fraction=1.0):
    spec = np.zeros_like(logEmin)

    for ii,(logElo, logEhi) in enumerate(zip(logEmin,logEmax)):
        Elo = 10**logElo 
        Ehi = 10**logEhi 
        spec[ii] = integrate.quad( hess_crab_spectrum, Elo,Ehi,
                                   args=(fraction,) )[0]

    return units.Quantity(spec, "ct cm**-2/s**-1")

