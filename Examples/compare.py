from  commonsens import sensitivity, inputs, config
from matplotlib import pyplot as plt
from astropy import units

if __name__ == '__main__':
    

    # load the example files
    config.gamma_energy_migration_method = "matrix"
    config.electron_energy_migration_method = "matrix"
    config.proton_energy_migration_method = "functional"
    

    datasets = {
#        "Thomas":"thomas_mono_*_zen020_az180_off0.50.fits",
        "Markus Old":"markus_mono_*.fits",
        "Markus Tight": "markus_tight_*.fits",
        "Dan Zeta Loose": "dan-zeta-loose-*.fits"
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

    plt.show()
