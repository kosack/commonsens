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
from commonsens import sensitivity as sens
from matplotlib import pyplot as plt
from astropy import units
from optparse import OptionParser
import sys
import os

if __name__ == '__main__':
    
    parser = OptionParser()
    # parser.add_option("-l","--label", dest="label", 
    #                   help="name of analysis")
    parser.add_option("-b","--batch", dest="batch", 
                      help="batch mode, don't display")
    parser.add_option("-w","--write", dest="write", 
                      help="write out the plots")


    parser.set_usage("calcsens.py [options] <performance file> ")
    (options, args) = parser.parse_args()


    for filename in args:
        
        name, _ = os.path.splitext( os.path.basename( filename ) )
        gammas,electrons,protons = cs.inputs.load_all_from_fits( filename )

        bgrate,rate_p,rate_e = sens.calc_background_rate( gammas,
                                                          electrons,
                                                          protons,
                                                          return_all=True)

        gamma_aeff = gammas.effective_area_reco()
        deltaE = gammas.delta_e


        # make sensitivity plot for several parameters:
        e2 = True
        plt.figure()
        for hours in [0.5,5,50]:
            result = sens.calc_sensitivity( name,
                                            bgrate,gamma_aeff,deltaE,
                                            obstime=hours*units.h,)
            sens.plot_sensitivity( gammas.log_e, result,esquared=e2 )

        sens.plot_crab( gammas.log_e, esquared=e2) # overlay Crab contours
        plt.legend(loc="best")
        plt.grid(alpha=0.3)    
        if (options.write):
            plt.savefig( name+"_sensitivity.pdf")

        # do the same in crab units:
        plt.figure()
        for hours in [0.5,2,5,50]:
            result = sens.calc_sensitivity( name, bgrate,gamma_aeff,
                                            deltaE, obstime=hours*units.h)
            sens.plot_sensitivity_crabunits( gammas.log_e, result )

        plt.legend(loc='best', fontsize='small')
        plt.grid(alpha=0.3)    
        if (options.write):
            plt.savefig( name+"_sensitivity_crab.pdf")

        # plot some of the intermediate distributions for debugging:
        if True:
            plt.figure(figsize=(10,9))
            plt.subplot(2,2,1)
            sens.plot_effareas( gammas, electrons, protons )
            plt.subplot(2,2,2)
            sens.plot_count_distributions( gammas.log_e, result )
            plt.subplot(2,2,3)
            sens.plot_significances( gammas.log_e, result )
            plt.subplot(2,2,4)
            sens.plot_rates( gammas.log_e, rate_p, rate_e, result )
            plt.gcf().suptitle("Intermediate Distributions")
            if (options.write):
                plt.savefig( name+"_intermediate.pdf")

        # plot the underlying distributions
        if True:
            gammas.plot()
            if (options.write):
                plt.savefig( name+"_gammas.pdf")

            electrons.plot()
            if (options.write):
                plt.savefig( name+"_electrons.pdf")

            protons.plot()
            if (options.write):
                plt.savefig( name+"_protons.pdf")

        if not options.batch:
            plt.show()
