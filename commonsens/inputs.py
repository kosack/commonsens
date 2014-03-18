"""This module contains a class (:class:`~commonsens.particledist.ParticleDistribution`)  used to store and calculate information about a simulated particle distribution for a single particle species, as well as some utility functions for loading ParticleDistributions

`ParticleDistributions` of gammas, electrons, and protons are used in
the sensitivity calculation.

ParticleDistributions may be constructed from compatible FITS files
using the from_fits() constructor.  

>>> gammas = GammaDistribution.from_fits( "mysens.fits" )
>>> gammas.plot()  # show some debugging info

you can load gammas,electrons, and protons at once, and set some
default values if you have named your files appropriately (the words
"gamma","electron", and "proton" are substituted for the *). Use:

>>> gammas,electrons,protons = load_all_from_fits( "mysens.fits")



Code Documentation:
===================

"""

__all__ = ['ParticleDistribution', 'load_all_from_fits']

import numpy as np
from matplotlib import pyplot as plt
from math import pi
from scipy.interpolate import RectBivariateSpline, interp1d
from astropy.io import fits
from astropy import units

from commonsens import spectra
from commonsens import config


def extrap(x, xp, yp):
    """np.interp function with linear extrapolation"""
    y = np.interp(x, xp, yp, kind="linear", bounds_error=False)
    y[x < xp[0]] = yp[0] + (x[x<xp[0]]-xp[0]) * (yp[0]-yp[1]) / (xp[0]-xp[1])
    y[x > xp[-1]]= yp[-1] + (x[x>xp[-1]]-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2])
    return y


def normalize_to_prob(x):
    """
    Makes x a probability distribution 
    
    Arguments:
    - `x`: array-like

    >>> np.apply_along_axis( normalize_to_prob, arr=migmatrix, axis=1)

    """
    x[x>0] /= np.sum(x)
    return x

def zero_if_low_stats(x):
    if sum(x>0) < 4:
        x *= 0
    return x


def functional_energy_migration(log_e_true, value, migration_function):
    """very simplistic energy migration (doesn't take into account any
    spread in the energy, just the bias. Returns an array of value as
    a func of log_e_reco (in the same energy bins as log_e_true)
    
    :param log_e_true: array of  log of true energy
    :param value: array of value at log_e_true
    :param migration_function: simple function that returns e_reco given e_true

    """
    
    func = interp1d( log_e_true, value, kind="linear", 
                     bounds_error=False, fill_value=0)
    log_ereco = np.log10( migration_function( 10**log_e_true ) )
    return func(log_ereco)

def matrix_energy_migration( log_e_true, value, log_e_reco, matrix):
    """Migration-Matrix energy migration. Uses a matrix (e_true vs
    e_reco) to migrate from true to reconstructed energy.
    
    :param log_e_true: array N of true energy bin centers
    :param value: array of values for each true energy bin
    :param log_e_reco: array of M reco energy bin centers
    :param matrix: NxM array giving probability of Ereco vs Etrue
    """
    N,M = matrix.shape

    recovalues = np.zeros_like(value)

    for ii in range(N):
        for jj in range(M):
            recovalues[jj] += value[ii] * matrix[ii,jj]

    return recovalues

def _rebin1d(arr, ntimes=1):

    """ helper to rebin a 1D distribution (summing bins) """
    shp = arr.shape[0]//2**ntimes
    sh2 = shp, arr.shape[0]//shp
    return arr.reshape(sh2).sum(1)


def _rebin2d(arr, shape):
    """ helper to rebin a 2D distribution (summing bins) """
    sh = shape[0], arr.shape[0]//shape[0],shape[1], arr.shape[1]//shape[1]
    return arr.reshape(sh).sum(-1).sum(1)

class ParticleDistribution(object):
    """A set of binned data describing a simulated particle distribution
    that has been processed by a gamma-ray analysis code.
   
     Each data member is an array containing the value it is named,
     with energy bins defined by log_e_lo and log_e_hi. The binning
     must be the same for all variables.

     It is assumed that the binning in E_true is the same as the
     binning in E_reco
    """

    def __init__(self, log_e_lo, log_e_hi, name=""):
        """
        - :name: name of particle distribution (e.g. "electrons")
        - :log_e_lo: array of lower energy bin values
        - :log_e_hi: array of upper energy bin values
        """
        self._name = name


        self.log_e_lo = log_e_lo #: lower edge of energy bins
        self.log_e_hi = log_e_hi #: upper edge of energy bins
        self.n_detected = None #: counts detected per E_reco bin
        self.thetasqr = None  #: angular cut used per bin
        self.n_simulated = None #: counts simulated per E_true bin
        self.r_simulated = None #: radius of simulation per E_true bin
        self.phi_diffuse = None #: diffuse simulation cone angle 
        self.e_mig = None #: 2D energy migration matrix (E_true vs E_reco)
        self.r68_psf = None #: 68% containment radius of PSF

        self._e_mig_normalized = None #: normalized 2D energy migration matrix
        self._dnde_true_func = lambda e_true : 0.0*e_true
        self._migration_function = lambda e_true : e_true # default Etrue=Ereco
#        self._energy_migration_method = "matrix"
        self._energy_migration_method = "shifted"


    def set_energy_migration_method(self, method ):
        if method not in ("matrix", "functional"):
            raise ValueError("unknown energy migration method")
        self._energy_migration_method = method
        

    @property
    def e_mig_normalized(self):
        """
        returns the energy migration matrix, normalized such that the
        integral along the true energy is 1.0
        """
        if self._e_mig_normalized == None:
            self._e_mig_normalized = np.apply_along_axis( normalize_to_prob, 
                                                          arr=self.e_mig, 
                                                          axis=1)
        return self._e_mig_normalized


    def _migrate_etrue_to_ereco(self, value ):
        """
        apply energy migration to go from true to reconstructed energy on
        the x-axis
        """
        print self._name," migration via ", self._energy_migration_method

        if self._energy_migration_method == "functional":
            return functional_energy_migration( self.log_e, value,
                                             self._migration_function) 
        elif self._energy_migration_method == "matrix":
            return matrix_energy_migration( self.log_e, value,
                                            self.log_e, self.e_mig_normalized )
        else:
            raise ValueError("Energy Migration Method not implemented")


    @property
    def n_simulated_reco(self):
        return self._migrate_etrue_to_ereco( self.n_simulated )

    @property
    def dnde_reco(self):
        """return :math:`dN/dE_{reco}` using the energy migration function or
        matrix
        """
        return self.dnde_true.unit \
            * self._migrate_etrue_to_ereco( self.dnde_true.value )

    def set_spectrum(self, specfunc):
        """
        Set the differential particle spectrum function (a mathematical
        function of E_true) for this particle species.  

        :param specfunc:  function that returns dNdE(E_true), units of 1/TeV/s/cm^2
        """
        self._dnde_true_func = specfunc
        self.dnde_true = specfunc(10**self.log_e)  #value at bin center

    @property
    def log_e(self):
        """ centers of energy bins """
        return (self.log_e_hi+self.log_e_lo)/2.0

    @property
    def nbins(self):
        """ returns number of bins in energy"""
        return len(self.log_e_hi)


    def effective_area_reco(self, ):
        """
        effective area in reconstructed energy 
        for the analysis that was performed for this
        particle species (including the theta^2 cut)
        """
        # correct for cone angle?
        #n_det = self.n_detected * self.thetasqr/(self.phi_diffuse)
        Aeff = (self.n_detected / self.n_simulated_reco) * pi \
               * self.r_simulated.to(units.m)**2
        return Aeff


    @property
    def delta_e(self):
        """ energy bin widths (not log) """
        return (10**self.log_e_hi - 10**self.log_e_lo)*units.TeV

    def rate_per_solidangle(self, ):
        """
        returns predicted detection rate of this particle species,
        based on the given differential source spectrum (set by set_spectrum)
        """

        E = 10**self.log_e * units.TeV
        deltaE = (10**self.log_e_hi - 10**self.log_e_lo) * units.TeV
        # TODO: should really integrate spectrum for each bin, rather
        # than use rectangle-rule approximation.  Or at least use
        # trapezoid rule... 

        return self.effective_area_reco()*self.dnde_reco*deltaE


    def plot(self):
        """ display some debugging plots for this particle distribution """
        plt.figure(figsize=(10,5))

        # plot the detected and simulated counts (for e_True and
        # e_reco)

        plt.subplot(1,2,1)
        plt.title( self._name )
        plt.loglog( 10**self.log_e, self.n_simulated,color='gray',
                      linewidth=3.0,
                      label="N_sim(E_true)", drawstyle='steps-mid' )
        plt.loglog( 10**self.log_e, self.n_simulated_reco, color='r',
                      linewidth=3.0,
                      label="N_sim(E_reco)", drawstyle='steps-mid' )
        plt.loglog( 10**self.log_e, self.n_detected,color='r', 
                      linewidth=3.0,
                      linestyle="--", label="N_det(E_reco)",
                      drawstyle='steps-mid')
        plt.xlabel("log10(E/TeV)")
        plt.ylabel("N (counts)")
        plt.legend(loc="best", fontsize="small")

        # plot the energy migration matrix + the simple energy
        # migration function + the e_bias given in the data file to
        # see if they compare.

        plt.subplot(1,2,2)
        plt.title("Energy Migration")
        plt.pcolormesh(self.log_e_lo, self.log_e_lo, self.e_mig_normalized.T)
        # plt.errorbar( self.log_e, np.log10(10**self.log_e + self.e_bias), 
        #               self.e_res, color='w', lw=3 )
        plt.plot(self.log_e, 
                 np.log10(self._migration_function(10**self.log_e)),
                 color="red", lw=3, ls="--")
        plt.xlabel("log_e_true")
        plt.ylabel("log_e_reco")
        plt.colorbar()

    @classmethod
    def from_fits(cls, filename, hduname):
        """
        Construct a ParticleDistribution from a FITS file that has a SENS
        extension

        :param name: identifier for distribution (e.g. "electrons")
        :param filename: FITS file with a SENS extension
        """

        print "LOADING: {0}".format(filename)

        sens = fits.open(filename)[hduname]
        log_e_lo = sens.data.field("LOG10_E_LO") 
        log_e_hi = sens.data.field("LOG10_E_HI")

        part = cls( log_e_lo, log_e_hi )
        part.n_detected = sens.data.field("N_detected")
        part.n_simulated = sens.data.field("N_simulated")
        part.r_simulated = sens.data.field("r_simulated") * units.meter
        part.thetasqr = sens.data.field("ThetaSqr") * units.deg
        part.phi_diffuse = sens.data.field("phi_diffuse") * units.deg
#        part.r68_psf = sens.data.field("R68_psf") * units.meter
        part.e_mig =  sens.data.field("E_migration") 
        
        return part


class GammaDistribution(ParticleDistribution):
    def __init__(self, log_e_lo, log_e_hi):
        super(GammaDistribution,self).__init__(log_e_lo, log_e_hi, name="gammas")
        self.set_energy_migration_method(config.gamma_energy_migration_method)
        self.set_spectrum( spectra.hess_crab_spectrum )
        
class ElectronDistribution(ParticleDistribution):
    """ Electron ParticleDistribution """
    def __init__(self, log_e_lo, log_e_hi):
        super(ElectronDistribution,self).__init__(log_e_lo, log_e_hi,name="electrons")
        self.set_spectrum( spectra.electron_spectrum_fermi )
        self.set_energy_migration_method(config.electron_energy_migration_method)

class ProtonDistribution(ParticleDistribution):
    """ Proton ParticleDistribution """
    def __init__(self, log_e_lo, log_e_hi):
        super(ProtonDistribution,self).__init__(log_e_lo, log_e_hi,name="protons")
        self.set_spectrum( spectra.cosmicray_spectrum )
        self._migration_function = lambda e_true : e_true*3.0 
        self.set_energy_migration_method(config.proton_energy_migration_method)

def load_all_from_fits(fitsfile):
    """Load set of particle species (gammas, electrons, hadrons) from a
    FITS files containing GAMMA, ELECTRON, and PROTON HDUs
    
    :returns: (gammas,electrons,protons) 
    """

    gammas = GammaDistribution.from_fits(fitsfile,"GAMMA")
    electrons = ElectronDistribution.from_fits( fitsfile, "ELECTRON")
    protons = ProtonDistribution.from_fits( fitsfile, "PROTON")
    return (gammas,electrons,protons)
    


def load_cta_response(filename):
    """read Bernloehr-style CTA text response file, that is already
    pre-processed (e.g. background rate is already calculated). This
    is included for comparing older CTA sensitivity output to those
    calculated with ParticleDistributions
    
    Assumes the text file has the columns:
    * 0: log10(E) at the bin center (TeV)
    * 1: A_eff for gammas, after all cuts  (m^2)
    * 2-4: r28,r80,Eres (not used)
    * 5: background rate post-cuts, and inside source region (Hz)
    * 6: differential sensitivity pre-calculated, E^2 dN/dE (erg/cm^2/s)
    """
    AA = np.loadtxt(filename)
    log_e = AA[:,0]
    dloge = log_e[1]-log_e[0]

    log_e_lo = log_e - dloge/2.0
    log_e_hi = log_e + dloge/2.0

    sens = AA[:,6]

    return (log_e_lo, log_e_hi, AA[:,1] *units.m**2,
            AA[:,5] * 1/units.s, 
            units.Quantity(AA[:,6],"erg/(cm**2 s)") )
    

def load_cta_response_fits(filename):
    """ read pre-computed Aeff, BgRate, etc from a CTA FITS file """

    table = fits.open(filename)['CTASENS']
    return (table.data.field("LOG10_E_LO"),
            table.data.field("LOG10_E_HI"),
            table.data.field("EffectiveArea") * units.m**2,
            table.data.field("BGRate") * 1.0/units.s,
            units.Quantity(table.data.field("DiffSens"),"erg/(cm**2 s)"))
            
            
    
