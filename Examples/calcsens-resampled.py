from commonsens import inputs
from commonsens.sensitivity import *
from matplotlib import pyplot as plt
from astropy import units

if __name__ == '__main__':
    

    resample = True

    filepat = "thomas_mono_*_zen020_az180_off0.50.fits"
    gammas,electrons,protons = inputs.load_all_from_fits(filepat)


    # also try with resampled distributions:
    rgammas = gammas.resample( -1.7, 0.0, 50 )
    relectrons = electrons.resample( -1.7, 0.0, 50 )
    rprotons = protons.resample( -1.7, 0.0, 50 )


    


    # make sensitivity plot for several parameters:
    plt.figure()
    for hours in [0.5,5,50]:
        out = calc_from_distributions( "thomas",gammas,electrons,protons,
                                obstime=hours*units.h)
        plot_sensitivity( gammas.log_e, out )

        rout = calc_from_distributions( "thomas_r",rgammas,relectrons,rprotons,
                                obstime=hours*units.h)
        plot_sensitivity( rgammas.log_e, rout )

    plot_crab( gammas.log_e )
    plt.legend(loc="best")
    plt.grid(alpha=0.3)    


    # # do the same in crab units:
    # plt.figure()
    # for hours in [0.5,2,5,0]:
    #     plot_sensitivity_crabunits( calc_sensitivity("Thomas", 
    #                                                  gammas,
    #                                                  electrons,
    #                                                  protons,
    #                                                  obstime=hours*units.h))

    #     plot_sensitivity_crabunits( calc_sensitivity("Thomas_r", 
    #                                                  rgammas,
    #                                                  relectrons,
    #                                                  rprotons,
    #                                                  obstime=hours*units.h))
    # plt.legend(loc='best', fontsize='small')
    # plt.grid(alpha=0.3)    



    if True:
        # debugging plots:
        plt.figure(figsize=(8,8))
        plt.subplot(2,2,1)
        plot_effareas( gammas, electrons, protons )
        plot_effareas( rgammas, relectrons, rprotons )
        plt.subplot(2,2,2)
        plot_count_distributions( gammas.log_e, out )
        plot_count_distributions( rgammas.log_e, rout )
        plt.subplot(2,2,3)
        plot_significances( gammas.log_e, out )
        plot_significances( rgammas.log_e, rout )
        plt.subplot(2,2,4)

        plt.show()
