from astropy import units
import gammasens as gs
from matplotlib import pyplot as plt

respdir = "../Data/CTA-Prod1-Bernloehr/"
respfile="kb_E_0.5h_20deg_v3.dat.txt"
#respfile ="kb_E_50h_50deg_v3.dat.txt"

(log_e_lo, log_e_hi, 
 a_eff_reco, background_rate, precalc_sens) = gs.inputs.load_cta_response(respdir+respfile) 

delta_e = (10**(log_e_hi)- 10**(log_e_lo)) * units.TeV
log_e = (log_e_hi + log_e_lo)/2.0


for obstime in [0.5,50]*units.h:

    sens = gs.sensitivity.calc_sensitivity( respfile, 
                                            background_rate,
                                            a_eff_reco, delta_e, 
                                            obstime=obstime,
                                            num_bg_regions=3,
                                            verbose=True)

    gs.sensitivity.plot_sensitivity( log_e, sens )



# plot the sensitivity calculated by this script, and overlay the
# precalculated one (note that the CTA sensitiity is in E^2 units)gs.sensitivity.plot_crab( log_e )
E = 10**log_e * units.TeV
presens = (precalc_sens/E**2).to("1/(cm**2 s TeV)")
plt.plot( log_e, presens, linewidth=4, linestyle="--", color="red", alpha=0.5,
          label="precalculated")
plt.legend()
plt.grid()

# display residuals
plt.figure()
plt.scatter( log_e, sens.sensitivity/units.ct - presens)
plt.title("Residuals between CTA and this calculation")
plt.xlabel("Log10(E_reco)")
plt.grid()






#S = sens['sensitivity']
#plt.plot( log_e, (S*3)/E )  # need to find where error is: get correct
# answer when multiplying by 3/E, so probably a problem in the
# proton spectrum

plt.show()
