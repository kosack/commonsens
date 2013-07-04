from gammasens import inputs
from gammasens.sensitivity import *
from matplotlib import pyplot as plt
from astropy import units

if __name__ == '__main__':
    

    resample = True

    filepat = "thomas_mono_*_zen020_az180_off0.50.fits"
    gammas,electrons,protons = inputs.loadAllFromFITS("thomas", filepat)


    # also try with resampled distributions:
    rgammas = gammas.getResampledDistribution( -1.4, 0.0, 50 )
    relectrons = electrons.getResampledDistribution( -1.4, 0.0, 50 )
    rprotons = protons.getResampledDistribution( -1.4, 0.0, 50 )



    # make sensitivity plot for several parameters:
    plt.figure()
    for hours in [0.5,5,50]:
        out = calc_sensitivity( "thomas",gammas,electrons,protons,
                                obstime=hours*units.h)
        plot_sensitivity( out )

        rout = calc_sensitivity( "thomas_r",rgammas,relectrons,rprotons,
                                obstime=hours*units.h)
        plot_sensitivity( rout )

    plot_crab( out )
    plt.legend(loc="best")
    plt.grid(alpha=0.3)    


    # do the same in crab units:
    plt.figure()
    for hours in [0.5,2,5,0]:
        plot_sensitivity_crabunits( calc_sensitivity("Thomas", 
                                                     gammas,
                                                     electrons,
                                                     protons,
                                                     obstime=hours*units.h))

        plot_sensitivity_crabunits( calc_sensitivity("Thomas_r", 
                                                     rgammas,
                                                     relectrons,
                                                     rprotons,
                                                     obstime=hours*units.h))
    plt.legend(loc='best', fontsize='small')
    plt.grid(alpha=0.3)    



    if True:
        # debugging plots:
        plt.figure(figsize=(8,8))
        plt.subplot(2,2,1)
        plot_effareas( gammas, electrons, protons )
        plot_effareas( rgammas, relectrons, rprotons )
        plt.subplot(2,2,2)
        plot_count_distributions( out )
        plot_count_distributions( rout )
        plt.subplot(2,2,3)
        plot_significances( out )
        plot_significances( rout )
        plt.subplot(2,2,4)
        plot_rates( out )
        plot_rates( rout )

#    plt.show()
