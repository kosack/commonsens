import gammasens as gs
from matplotlib import pyplot as plt
from astropy import units

if __name__ == '__main__':
    

    # load the example files
    filepat = "thomas_mono_*_zen020_az180_off0.50.fits"
    gammas,electrons,protons = gs.inputs.load_all_from_fits(filepat)
    bgrate,rate_p,rate_e = gs.sensitivity.calc_background_rate( gammas, 
                                                                electrons, 
                                                                protons, 
                                                                return_all=True)
    gamma_aeff = gammas.effective_area_reco()
    deltaE = gammas.delta_e


    # make sensitivity plot for several parameters:
    plt.figure()
    for hours in [0.5,5,50]:
        result = gs.sensitivity.calc_sensitivity( "thomas",
                                                  bgrate,gamma_aeff,deltaE,
                                                  obstime=hours*units.h)
        gs.sensitivity.plot_sensitivity( gammas.log_e, result )

    gs.sensitivity.plot_crab( gammas.log_e ) # overlay Crab contours
    plt.legend(loc="best")
    plt.grid(alpha=0.3)    


    # do the same in crab units:
    plt.figure()
    for hours in [0.5,2,5,50]:
        gs.sensitivity.plot_sensitivity_crabunits( gammas.log_e, \
                            gs.sensitivity.calc_sensitivity("Thomas", 
                                                            bgrate,gamma_aeff,
                                                            deltaE, 
                                                            obstime=hours*units.h))
    plt.legend(loc='best', fontsize='small')
    plt.grid(alpha=0.3)    


    # plot some of the intermediate distributions for debugging:
    if True:
        plt.figure(figsize=(10,9))
        plt.subplot(2,2,1)
        gs.sensitivity.plot_effareas( gammas, electrons, protons )
        plt.subplot(2,2,2)
        gs.sensitivity.plot_count_distributions( gammas.log_e, result )
        plt.subplot(2,2,3)
        gs.sensitivity.plot_significances( gammas.log_e, result )
        plt.subplot(2,2,4)
        gs.sensitivity.plot_rates( gammas.log_e, rate_p, rate_e, result )
        plt.gcf().suptitle("Intermediate Distributions")

    plt.show()
