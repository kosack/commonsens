from  commonsens import sensitivity, inputs, config
from matplotlib import pyplot as plt
from astropy import units

if __name__ == '__main__':
    

    # load the example files
    config.gamma_energy_migration_method = "matrix"
    config.electron_energy_migration_method = "matrix"
    config.proton_energy_migration_method = "functional"
    

    datasets = {
        "ParisMVA Hybrid":"PerfData/ParisMVA_Hybrid.fits",
        "APR": "PerfData/APR.fits",
        "Loose Zeta": "PerfData/dan_zeta_loose.fits",
        "Tight": "PerfData/tight.fits",
        "Thomas Mono": "PerfData/thomas_mono.fits"       
    }

    e2=True

    for name in datasets:
        gammas,electrons,protons = inputs.load_all_from_fits(datasets[name])

        result=sensitivity.calc_from_distributions( name, 
                                                    gammas,
                                                    electrons,
                                                    protons,
                                                    obstime=50*units.h)
        
        sensitivity.plot_sensitivity( gammas.log_e, result, esquared=e2 )


    sensitivity.plot_crab( gammas.log_e, esquared=e2) # overlay Crab contours
    plt.legend(loc="best")
    plt.grid(alpha=0.3)    
    plt.xlim( 0.02, 100 )

    plt.show()
