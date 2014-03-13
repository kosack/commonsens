"""This module contains a class (:class:`~commonsens.particledist.ParticleDistribution`)  used to store and calculate information about a simulated particle distribution for a single particle species, as well as some utility functions for loading ParticleDistributions

`ParticleDistributions` of gammas, electrons, and protons are used in
the sensitivity calculation.

ParticleDistributions may be constructed from compatible FITS files
using the from_fits() constructor.  They may also be resampled to a new binning.

>>> gammas = ParticleDistribution.from_fits( "gammas", "mysens-gammas.fits" )
>>> gammas.plot()  # show some debugging info

you can load gammas,electrons, and protons at once, and set some
default values if you have named your files appropriately. Use:

>>> gammas,electrons,protons = load_all_from_fits( "MyAnalysis", "mysens-*.fits")


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

from commonsens import spectra, config


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

        self._dnde_true_func = lambda e_true : 0.0*e_true
        self._migration_function = lambda e_true : e_true # default Etrue=Ereco
        self._energy_migration_method = "matrix"

    def _resample(self, newlog_e, values ):
        """ resample values to new bin centers """
        spline =  interp1d( self.log_e, values, kind="linear",
                            bounds_error=False, fill_value=0 )  
        return spline( newlog_e )

    def set_energy_migration_method(self, method ):
        if method not in ("matrix", "functional"):
            raise ValueError("unknown energy migration method")
        self._energy_migration_method = method
        

    def resample(self, log_e_min, log_e_max, nbins):
        """
        returns new ParticleDistribution interpolated using splines to a
        new binning. Be careful to set the min and max in valid ranges
        where there is good statistics per bin in the original
        distribution, otherwise the interpolation will go crazy

        note that this creates new binning from an old one, but does
        not increase statistics. To do that you should rebin first,
        and then resample. Also note that only summable quantities can be
        rebinned, others should be interpolated.

        :param log_e_min: new lowest energy (left edge of first bin)
        :param log_e_max: new highest energy (right edge of last bin)
        :param nbins: new number of bins
        """
        
        therange = np.linspace(log_e_min,log_e_max, nbins+1)
        log_e_lo = therange[0:-1]
        log_e_hi = therange[1:]
        newlog_e = (log_e_hi+log_e_lo)/2.0

        dist = ParticleDistribution( self._name+"_resamp", log_e_lo,log_e_hi)
        newlog_e = dist.log_e

        # kind of a hack:
        # loop over values in the class dictionary (to avoid having to
        # write out each data member). Skip members starting with _ or
        # ending in _mig, and treat the migration matrix separately:
        for key in self.__dict__.keys():
            if  (not key.startswith("_")) and  (not key.endswith("_mig")):
                val = self.__dict__[key]
                if type(val) == units.Quantity:
                    dist.__dict__[key] = self._resample( newlog_e, 
                                                         val.value )*val.unit
                else:
                    dist.__dict__[key] = self._resample( newlog_e, val )

        #resample the migration matrix
        
        miginterp = RectBivariateSpline( self.log_e,self.log_e, self.e_mig )
        dist.e_mig = miginterp( newlog_e,newlog_e )
        
        dist._migration_function = self._migration_function
        dist._dnde_true_func = self._dnde_true_func
        dist._energy_migration_method = self._energy_migration_method

        return dist

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
                                            self.log_e, self.e_mig )
        else:
            raise ValueError("Energy Migration Method not implemented")


    def generateBiasFromMigration(self):
        pass
        # for ii in range(self.nbins):
        #     avg = np.average( self.log_e, weights=self.e_mig[ii] )
        


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
        plt.semilogy( self.log_e, self.n_simulated,color='gray',
                      linewidth=3.0,
                      label="N_sim(E_true)", drawstyle='steps-mid' )
        plt.semilogy( self.log_e, self.n_simulated_reco, color='r',
                      linewidth=3.0,
                      label="N_sim(E_reco)", drawstyle='steps-mid' )
        plt.semilogy( self.log_e, self.n_detected,color='r', 
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
        plt.pcolormesh(self.log_e_lo, self.log_e_lo, self.e_mig.T)
        # plt.errorbar( self.log_e, np.log10(10**self.log_e + self.e_bias), 
        #               self.e_res, color='w', lw=3 )
        plt.plot(self.log_e, 
                 np.log10(self._migration_function(10**self.log_e)),
                 color="red", lw=3, ls="--")
        plt.xlabel("log_e_true")
        plt.ylabel("log_e_reco")
        plt.colorbar()

    @classmethod
    def from_fits(cls, filename):
        """
        Construct a ParticleDistribution from a FITS file that has a SENS
        extension

        :param name: identifier for distribution (e.g. "electrons")
        :param filename: FITS file with a SENS extension
        """

        print "LOADING: {0}".format(filename)

        sens = fits.open(filename)['SENS']
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

        # normalize the migration matrix to be a probability:
        # TODO: which axis? along E_true or E_reco?
        np.apply_along_axis( normalize_to_prob, arr=part.e_mig, axis=1)

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
        self._migration_function = lambda e_true : e_true/3.0 
        self.set_energy_migration_method(config.proton_energy_migration_method)

def load_all_from_fits(filepattern, species = ['gamma','electron','proton']):
    """Load set of particle species (gammas, electrons, hadrons) assuming
    the files are named identically, but with the character "*"
    replaced with "gamma", "electron" or "proton" in the filename.

    :param name: identifier of this data set
    :param filepattern: filename, with the particle type replaced with *
    """

    dists = []

    for particle in species:
        if particle == 'gamma':
            dist = GammaDistribution.from_fits(filepattern.replace("*", particle))
        elif particle == 'electron':
            dist = ElectronDistribution.from_fits(filepattern.replace("*", particle))          
        elif particle == 'proton':
            dist = ProtonDistribution.from_fits(filepattern.replace("*", particle))     
        else: # unknown particle
            dist = ParticleDistribution.from_fits(filepattern.replace("*",particle), name=particle)         


        dists.append(dist)

    return dists
    


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
    
