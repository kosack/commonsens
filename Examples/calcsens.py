import gammasens as gs
from matplotlib import pyplot as plt
from astropy import units

if __name__ == '__main__':
    

    # load the example files
    filepat = "thomas_mono_*_zen020_az180_off0.50.fits"
    gammas,electrons,protons = gs.inputs.loadAllFromFITS(filepat)



    # make sensitivity plot for several parameters:
    plt.figure()
    for hours in [0.5,5,50]:
        out = gs.sensitivity.calc_sensitivity( "thomas",gammas,electrons,protons,
                                               obstime=hours*units.h)
        gs.sensitivity.plot_sensitivity( out )

    gs.sensitivity.plot_crab( out ) # overlay Crab contours
    plt.legend(loc="best")
    plt.grid(alpha=0.3)    


    # do the same in crab units:
    plt.figure()
    for hours in [0.5,2,5,50]:
        gs.sensitivity.plot_sensitivity_crabunits( \
                            gs.sensitivity.calc_sensitivity("Thomas", 
                                                            gammas,
                                                            electrons,
                                                            protons,
                                                            obstime=hours*units.h))
    plt.legend(loc='best', fontsize='small')
    plt.grid(alpha=0.3)    


    # plot some of the intermediate distributions for debugging:
    if True:
        plt.figure(figsize=(10,9))
        plt.subplot(2,2,1)
        gs.sensitivity.plot_effareas( gammas, electrons, protons )
        plt.subplot(2,2,2)
        gs.sensitivity.plot_count_distributions( out )
        plt.subplot(2,2,3)
        gs.sensitivity.plot_significances( out )
        plt.subplot(2,2,4)
        gs.sensitivity.plot_rates( out )
        plt.gcf().suptitle("Intermediate Distributions")

    plt.show()
