#!/usr/bin/env python

#
# complete example of calcualting and plotting sensitivity. 
# 
# EXAMPLE OF HOW TO RUN:
#  calcsens.py -l "Dan Zeta Loose" -g dan-zeta-loose-gamma.fits -e dan-zeta-loose-electron.fits -p dan-zeta-loose-proton.fits
#
# OR VIA IPYTHON:
# ipython calcsens.py -- -l "Dan Zeta Loose" -g dan-zeta-loose-gamma.fits -e dan-zeta-loose-electron.fits -p dan-zeta-loose-proton.fits

import commonsens as cs
from matplotlib import pyplot as plt
from astropy import units
from optparse import OptionParser
import sys


if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-g","--gammas", dest="gammas", 
                      help="gamma-ray input file")
    parser.add_option("-e","--electrons", dest="electrons", 
                      help="electron input file")
    parser.add_option("-p","--protons", dest="protons", 
                      help="proton input file")
    parser.add_option("-l","--label", dest="label", 
                      help="name of analysis", default="")

    parser.set_usage("calcsens.py [options] -g <gammafile> -e <electronfile> -p <protonfile> ")

    (options, args) = parser.parse_args()

    if not (options.gammas and options.electrons and options.protons):
        parser.error("wrong arguments. See --help for more info")
        sys.exit(1)
        
    name = options.label

    gammas    = cs.inputs.GammaDistribution.from_fits( options.gammas )
    electrons = cs.inputs.ElectronDistribution.from_fits( options.electrons )
    protons   = cs.inputs.ProtonDistribution.from_fits(options.protons )

    bgrate,rate_p,rate_e = cs.sensitivity.calc_background_rate( gammas, 
                                                                electrons, 
                                                                protons, 
                                                                return_all=True)



    gamma_aeff = gammas.effective_area_reco()
    deltaE = gammas.delta_e


    # make sensitivity plot for several parameters:
    e2 = True
    plt.figure()
    for hours in [0.5,5,50]:
        result = cs.sensitivity.calc_sensitivity( name,
                                                  bgrate,gamma_aeff,deltaE,
                                                  obstime=hours*units.h,)
        cs.sensitivity.plot_sensitivity( gammas.log_e, result,esquared=e2 )

    cs.sensitivity.plot_crab( gammas.log_e, esquared=e2) # overlay Crab contours
    plt.legend(loc="best")
    plt.grid(alpha=0.3)    


    # do the same in crab units:
    plt.figure()
    for hours in [0.5,2,5,50]:
        cs.sensitivity.plot_sensitivity_crabunits( gammas.log_e, \
                            cs.sensitivity.calc_sensitivity(name, 
                                                            bgrate,gamma_aeff,
                                                            deltaE, 
                                                            obstime=hours*units.h))
    plt.legend(loc='best', fontsize='small')
    plt.grid(alpha=0.3)    


    # plot some of the intermediate distributions for debugging:
    if True:
        plt.figure(figsize=(10,9))
        plt.subplot(2,2,1)
        cs.sensitivity.plot_effareas( gammas, electrons, protons )
        plt.subplot(2,2,2)
        cs.sensitivity.plot_count_distributions( gammas.log_e, result )
        plt.subplot(2,2,3)
        cs.sensitivity.plot_significances( gammas.log_e, result )
        plt.subplot(2,2,4)
        cs.sensitivity.plot_rates( gammas.log_e, rate_p, rate_e, result )
        plt.gcf().suptitle("Intermediate Distributions")

    # plot the underlying distributions
    if True:
        gammas.plot()
        electrons.plot()
        protons.plot()

    plt.show()
