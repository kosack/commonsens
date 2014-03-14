#
# merge 3 fits files for gammas,electrons,protons into a single one
# with multiple HDUs
#

try:
    from astropy.io import fits 
except:
    import pyfits as fits

import sys

if __name__ == '__main__':
    
    gammas = sys.argv.pop(1)
    electrons = sys.argv.pop(1)
    protons = sys.argv.pop(1)
    output = sys.argv.pop(1)


    print "gammas",gammas
    print "electrons",electrons
    print "protons",protons
    print " ===> ",output

    ghdu = fits.open(gammas)['SENS']
    ehdu = fits.open(electrons)['SENS']
    phdu = fits.open(protons)['SENS']

    ghdu.name="GAMMA"
    ehdu.name="ELECTRON"
    phdu.name="PROTON"

    hdus = fits.HDUList( hdus=[fits.PrimaryHDU(),ghdu,ehdu,phdu] )
    hdus.writeto( output )
