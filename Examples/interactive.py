"""
A simple example of making an interactive sensitivity plot, with
sliders to adjust the parameters
"""

from commonsens import *
from commonsens.sensitivity import calc_from_distributions,plot_sensitivity_crabunits

from astropy import units
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons


if __name__ == '__main__':

#    filepat = "thomas_mono_*_zen020_az180_off0.50.fits"
#    name = "Thomas"
    filepat = "markus_mono_*.fits"
    name = "Markus"    

    filepat = "dan-zeta-loose-*.fits"
    name = "Dan Zeta Loose"    

    gammas,electrons,protons = inputs.load_all_from_fits(filepat)


    obstime_hrs  = 50.0
    min_sys_pct = 5.0
    min_signif = 5.0
    min_events = 10.0

    ax = subplot(111)
    subplots_adjust(left=0.15, bottom=0.30)
 
    axcolor = 'lightgoldenrodyellow'

    ax_obstime = axes([0.25, 0.11, 0.45, 0.03], axisbg=axcolor)
    slider_obstime= Slider(ax_obstime, 'Time (h)', 0.0, 100.0, valinit=obstime_hrs)

    ax_syspct = axes([0.25, 0.16, 0.45, 0.03], axisbg=axcolor)
    slider_syspct= Slider(ax_syspct, 'Systematics (%)', 0.0, 100.0, 
                          valinit=min_sys_pct)

    ax_minsignif = axes([0.25, 0.06, 0.45, 0.03], axisbg=axcolor)
    slider_minsignif= Slider(ax_minsignif, 'Min Signif', 1.0, 30.0, 
                             valinit=min_signif)

    ax_minevents = axes([0.25, 0.01, 0.45, 0.03], axisbg=axcolor)
    slider_minevents= Slider(ax_minevents, 'Min events', 0.0, 50.0, 
                             valinit=min_events)

    ax_method = axes([0.80, 0.01, 0.15, 0.2], axisbg=axcolor)
    method = RadioButtons(ax_method, ('matrix', 'mixed', 'func'))


    def update(val):
        obstime_hrs = slider_obstime.val
        min_sys_pct = slider_syspct.val
        min_signif = slider_minsignif.val
        min_events = slider_minevents.val
        sca(ax)
        ax.clear()
        out = calc_from_distributions( name, 
                                                   gammas, 
                                                   electrons,  
                                                   protons, 
                                                   obstime=obstime_hrs*units.h,
                                                   min_sys_pct=min_sys_pct,
                                                   min_signif=min_signif,
                                                   min_events=min_events)
        loglog()
        plot_sensitivity_crabunits( gammas.log_e, out  ) 
        ylim( 0.01, 10.0 )
        grid()
        draw()

    def update_method(val):
        if val == "matrix":
            gammas.set_energy_migration_method( "matrix")
            electrons.set_energy_migration_method( "matrix")
            protons.set_energy_migration_method( "matrix")
        if val == "mixed":
            gammas.set_energy_migration_method( "matrix")
            electrons.set_energy_migration_method( "matrix")
            protons.set_energy_migration_method( "functional")
        if val == "func":
            gammas.set_energy_migration_method( "functional")
            electrons.set_energy_migration_method( "functional")
            protons.set_energy_migration_method( "functional")
        update(0)

    method.on_clicked(update_method)



    slider_obstime.on_changed(update)
    slider_syspct.on_changed(update)
    slider_minsignif.on_changed(update)
    slider_minevents.on_changed(update)

    update(None)

    show()
    
    
    
