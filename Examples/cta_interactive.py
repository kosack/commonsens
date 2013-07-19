"""
A simple example of making an interactive sensitivity plot, with
sliders to adjust the parameters
"""

from gammasens import *
from gammasens.sensitivity import calc_sensitivity,plot_sensitivity_crabunits

from astropy import units
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons


if __name__ == '__main__':

    respdir = "../Data/CTA-Prod1-Bernloehr/"
    respfile="kb_E_0.5h_20deg_v3.dat.txt"
    #respfile ="kb_E_50h_50deg_v3.dat.txt"
    
    (log_e_lo, log_e_hi, 
     a_eff_reco, background_rate, 
     precalc_sens) = inputs.load_cta_response(respdir+respfile) 
    
    delta_e = (10**(log_e_hi)- 10**(log_e_lo)) * units.TeV
    log_e = (log_e_hi + log_e_lo)/2.0
    

    obstime_hrs  = 50.0
    min_sys_pct = 5.0
    min_signif = 5.0
    min_events = 10.0

    ax = subplot(111)
    subplots_adjust(left=0.15, bottom=0.30)
 
    axcolor = 'lightgoldenrodyellow'

    ax_obstime = axes([0.25, 0.11, 0.65, 0.03], axisbg=axcolor)
    slider_obstime= Slider(ax_obstime, 'Time (h)', 0.0, 100.0, valinit=obstime_hrs)

    ax_syspct = axes([0.25, 0.16, 0.65, 0.03], axisbg=axcolor)
    slider_syspct= Slider(ax_syspct, 'Systematics (%)', 0.0, 100.0, 
                          valinit=min_sys_pct)

    ax_minsignif = axes([0.25, 0.06, 0.65, 0.03], axisbg=axcolor)
    slider_minsignif= Slider(ax_minsignif, 'Min Signif', 1.0, 30.0, 
                             valinit=min_signif)

    ax_minevents = axes([0.25, 0.01, 0.65, 0.03], axisbg=axcolor)
    slider_minevents= Slider(ax_minevents, 'Min events', 0.0, 50.0, 
                             valinit=min_events)


    def update(val):
        obstime_hrs = slider_obstime.val
        min_sys_pct = slider_syspct.val
        min_signif = slider_minsignif.val
        min_events = slider_minevents.val
        sca(ax)
        ax.clear()
        out = calc_sensitivity( respfile, 
                                background_rate, a_eff_reco, delta_e,
                                obstime=obstime_hrs*units.h,
                                min_sys_pct=min_sys_pct,
                                min_signif=min_signif,
                                min_events=min_events)
        semilogy()
        plot_sensitivity_crabunits( log_e, out  ) 
        ylim( 0.0, 1.1 )
        grid()
        draw()


    slider_obstime.on_changed(update)
    slider_syspct.on_changed(update)
    slider_minsignif.on_changed(update)
    slider_minevents.on_changed(update)

    update(None)

    show()
    
    
    
