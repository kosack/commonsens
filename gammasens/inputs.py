"""This module contains a class (:class:`~gammasens.particledist.ParticleDistribution`)  used to store and calculate information about a simulated particle distribution for a single particle species, as well as some utility functions for loading ParticleDistributions

`ParticleDistributions` of gammas, electrons, and protons are used in
the sensitivity calculation.

ParticleDistributions may be constructed from compatible FITS files
using the fromFITS() constructor.  They may also be resampled to a new binning.

>>> gammas = ParticleDistribution.fromFITS( "gammas", "mysens-gammas.fits" )
>>> gammas.plot()  # show some debugging info

you can load gammas,electrons, and protons at once, and set some
default values if you have named your files appropriately. Use:

>>> gammas,electrons,protons = loadAllFromFITS( "MyAnalysis", "mysens-*.fits")


Code Documentation:
===================



"""

__all__ = ['ParticleDistribution','loadAllFromFITS']

import numpy as np
from matplotlib import pyplot as plt
from math import pi
from scipy.interpolate import RectBivariateSpline, InterpolatedUnivariateSpline
from astropy.io import fits
from astropy import units

from gammasens import spectra

def shifted_energy_migration(log_e_true, value, migration_function):
    """very simplistic energy migration (doesn't take into account any
    spread in the energy, just the bias. Returns an array of value as
    a func of log_e_reco (in the same energy bins as log_e_true)
    
    :param log_e_true: array of  log of true energy
    :param value: array of value at log_e_true
    :param migration_function: simple function that returns e_reco given e_true

    """
    
    func = InterpolatedUnivariateSpline( log_e_true, value )
    log_ereco = np.log10( migration_function( 10**log_e_true ) )
    return func(log_ereco)


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

    def __init__(self, name, log_e_lo, log_e_hi):
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
        self.e_res = None #: energy resolution 
        self.e_bias = None #: energy bias
        self.e_mig = None #: 2D energy migration matrix (E_true vs E_reco)
        self.r68_psf = None #: 68% containment radius of PSF

        self._dnde_true_func = lambda e_true : 0.0*e_true
        self._migration_function = lambda e_true : e_true # default Etrue=Ereco
        self._energy_migration_method = "shifted"

    def _resample(self, newlog_e, values ):
        """ resample values to new bin centers """
        spline =  InterpolatedUnivariateSpline( self.log_e, values )
        return spline( newlog_e )


    def getResampledDistribution(self, log_e_min, log_e_max, nbins):
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

    def migrateToEreco(self, value ):
        if self._energy_migration_method == "shifted":
            return shifted_energy_migration( self.log_e, value,
                                             self._migration_function) 
        else:
            raise ValueError("Energy Migration Method not implemented")


    def generateBiasFromMigration(self):
        pass
        # for ii in range(self.nbins):
        #     avg = np.average( self.log_e, weights=self.e_mig[ii] )
        


    @property
    def n_simulated_reco(self):
        return self.migrateToEreco( self.n_simulated )

    @property
    def dnde_reco(self):
        return self.dnde_true.unit * self.migrateToEreco( self.dnde_true.value )

    def setSpectrum(self, specfunc):
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
        return Aeff * units.m**2


    @property
    def deltaE(self):
        """ energy bin widths (not log) """
        return (10**self.log_e_hi - 10**self.log_e_lo)*units.TeV

    def rate_per_solidangle(self, ):
        """
        returns predicted detection rate of this particle species,
        based on the measured flux spectrum
        """

        E = 10**self.log_e * units.TeV
        deltaE = (10**self.log_e_hi - 10**self.log_e_lo) * units.TeV
        # TODO: should really integrate spectrum for each bin, rather
        # than use rectangle-rule approximation.  Or at least use
        # trapezoid rule... 

        return units.count*self.effective_area_reco()*self.dnde_reco*deltaE


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
        plt.legend(loc="best")

        # plot the energy migration matrix + the simple energy
        # migration function + the e_bias given in the data file to
        # see if they compare.

        plt.subplot(1,2,2)
        plt.title("Energy Migration")
        plt.pcolormesh(self.log_e_lo, self.log_e_lo, self.e_mig.T)
        plt.errorbar( self.log_e, np.log10(10**self.log_e + self.e_bias), 
                      self.e_res, color='w', lw=3 )
        plt.plot(self.log_e, 
                 np.log10(self._migration_function(10**self.log_e)),
                 color="red", lw=3, ls="--")
        plt.xlabel("log_e_true")
        plt.ylabel("log_e_reco")

    @classmethod
    def fromFITS(cls, name, filename):
        """
        Construct a ParticleDistribution from a FITS file that has a SENS
        extension

        :param name: identifier for distribution (e.g. "electrons")
        :param filename: FITS file with a SENS extension
        """

        print "LOADING: {0} ({1})".format(filename,name)

        sens = fits.open(filename)['SENS']
        log_e_lo = sens.data.field("LOG10_E_LO") 
        log_e_hi = sens.data.field("LOG10_E_HI")

        part = ParticleDistribution( name , log_e_lo, log_e_hi )
        part.n_detected = sens.data.field("N_detected")
        part.n_simulated = sens.data.field("N_simulated")
        part.r_simulated = sens.data.field("r_simulated") * units.meter
        part.thetasqr = sens.data.field("ThetaSqr") * units.deg
        part.phi_diffuse = sens.data.field("phi_diffuse") * units.deg
        part.r68_psf = sens.data.field("R68_psf") * units.meter
        part.e_bias = sens.data.field("E_bias") * units.TeV
        part.e_res = sens.data.field("e_res")  * units.TeV
        part.e_mig =  sens.data.field("E_migration") 

        return part







def loadAllFromFITS(filepattern, species = ['gamma','electron','proton']):
    """Load set of particle species (gammas, electrons, hadrons) assuming
    the files are named identically, but with the character "*"
    replaced with "gamma", "electron" or "proton" in the filename.

    :param name: identifier of this data set
    :param filepattern: filename, with the particle type replaced with *
    """

    dists = []

    for particle in species:
        dist = ParticleDistribution.fromFITS(particle, 
                                             filepattern.replace("*", particle))
        dists.append(dist)


    # set some defaults:
    dists[1].setSpectrum( spectra.electron_spectrum )
    dists[2].setSpectrum( spectra.cosmicray_spectrum )

    # simplistic proton migration function used in CTA curves (just
    # based on branching ratio)
    dists[2]._migration_function = lambda e_true : e_true*3.0 

    return dists
    
