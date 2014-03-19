from  commonsens import sensitivity, inputs, config
from matplotlib import pyplot as plt
from astropy import units
from astropy import convolution

def zeroLowStats():
    """
    """
    pass


if __name__ == '__main__':
    

    # load the example files
    config.enable_smoothing =False
    config.smooth_parameter=2.0

    e2=True
    bands = True
    
    datasets = {
        "ParisMVA Hybrid":"PerfData/ParisMVA_Hybrid.fits",
        "APR": "PerfData/APR.fits",
        "Loose Zeta": "PerfData/dan_zeta_loose.fits",
        "Mono Tight": "PerfData/Markus_Mono_Tight.fits",
        "Thomas Mono": "PerfData/thomas_mono.fits"       
    }

    colors = list(np.linspace(0,1,len(datasets)))
    plt.figure(figsize=(10,10))

    for name in datasets:
        gammas,electrons,protons = inputs.load_all_from_fits(datasets[name])

        result=sensitivity.calc_from_distributions( name, 
                                                    gammas,
                                                    electrons,
                                                    protons,
                                                    obstime=50*units.h )

        if bands:
            sensitivity.plot_sensitivity( gammas.log_e, result, 
                                          esquared=e2, linewidth=1,
                                          shade_percent=0.5, alpha=0.6, 
                                          color=cm.rainbow(colors.pop()))

        else:
            sensitivity.plot_sensitivity( gammas.log_e, result, 
                                          esquared=e2, linewidth=3,
                                          shade_percent=None, alpha=0.8, 
                                          color=cm.rainbow(colors.pop()))


    # overlay Crab contours
    sensitivity.plot_crab( gammas.log_e, esquared=e2) 


    # overlay CTA for reference
    (log_e_lo, log_e_hi, 
     a_eff_reco, background_rate, 
     precalc_sens) = inputs.load_cta_response_fits("PerfData/aar_south_sensitivity_logcut_ele2_ebin2_s030_b010_f005_r030_meff030_t050.00.fits")

    delta_e = (10**log_e_hi - 10**log_e_lo)*units.TeV
    ctasens = sensitivity.calc_sensitivity( "CTA Prod2 AAR", 
                                            background_rate,
                                            a_eff_reco, delta_e, 
                                            obstime=50*units.h,
                                            num_bg_regions=3,
                                            verbose=True)

    sensitivity.plot_sensitivity( 0.5*(log_e_hi+log_e_lo), ctasens, 
                                  esquared=e2, linewidth=2, linestyle="dashed" )
    
    
    plt.vlines( [0.030,0.050], 1e-14,1e-5, linestyle='dotted', linewidth=2)

    plt.legend(loc="best", ncol=2, frameon=False, fontsize="medium")
    plt.grid(alpha=0.3)    
    plt.xlim( 0.02, 100 )
    plt.show()
