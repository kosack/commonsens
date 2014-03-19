try:
    from astropy.io import fits 
except:
    import pyfits as fits
import ROOT
from ROOT import TFile, TH1D,TH2F
import sys
import os
import numpy as np
import fitshistogram

def get_hist_1d(rootfile,name):
    """

    Arguments:
    - `rootfile`:
    - `name`:
    """

    print "GET",name

    hist = rootfile.Get( name )
    xaxis = hist.GetXaxis()
    nbins = xaxis.GetNbins()

    # convert to numpy array:

    Ebin_lo = [ xaxis.GetBinLowEdge(ii+1) for ii in range(nbins)   ]
    Ebin_hi = [ xaxis.GetBinUpEdge(ii+1) for ii in range(nbins)   ]
    Ebin = [ xaxis.GetBinCenter(ii+1) for ii in range(nbins)   ]
    value = [hist.GetBinContent(ii+1) for ii in range(nbins) ]


    return (np.array(Ebin_lo), np.array(Ebin_hi),
            np.array(Ebin),np.array(value))

def get_hist_2d(rootfile,name):
    """

    Arguments:
    - `rootfile`:
    - `name`:
    """

    hist = rootfile.Get( name )
    xaxis = hist.GetXaxis()
    yaxis = hist.GetYaxis()
    nbinsx = xaxis.GetNbins()
    nbinsy = yaxis.GetNbins()

    # convert to numpy array:

    xEbin_lo = [ xaxis.GetBinLowEdge(ii+1) for ii in range(nbinsx)   ]
    xEbin_hi = [ xaxis.GetBinUpEdge(ii+1) for ii in range(nbinsx)   ]
    xEbin = [ xaxis.GetBinCenter(ii+1) for ii in range(nbinsx)   ]

    yEbin_lo = [ xaxis.GetBinLowEdge(ii+1) for ii in range(nbinsy)   ]
    yEbin_hi = [ xaxis.GetBinUpEdge(ii+1) for ii in range(nbinsy)   ]
    yEbin = [ xaxis.GetBinCenter(ii+1) for ii in range(nbinsy)   ]

    value = [[hist.GetBinContent(ii+1,jj+1) for ii in range(nbinsx)] \
             for jj in range(nbinsy)]

    return (np.array(xEbin_lo), np.array(xEbin_hi), np.array(xEbin),
            np.array(yEbin_lo), np.array(yEbin_hi), np.array(yEbin),
            value)


if __name__ == '__main__':

    hists_to_get = ['N_detected',
                    'N_simulated',
                    'r_simulated',
                    'ThetaSqr',
                    'phi_diffuse',
                    'R68_psf']
    #                    'E_res',
    #                    'E_bias' ]


    for rootfilename in  sys.argv:
        print rootfilename

        rootfile = TFile(rootfilename)


        E = None
        cols = []

        # get the 1D columns:
        for histname in hists_to_get:
            lo,hi,cen, val = get_hist_1d( rootfile, histname )

            if len(cols) == 0:
                cols.append(fits.Column( name="LOG10_E_LO", format="D",array=lo ))
                cols.append(fits.Column( name="LOG10_E_HI", format="D",array=hi ))

            cols.append( fits.Column( name=histname, format="D", array=val ))

        #insert the migration matrix:
        try:
            etlo,ethi,et,erlo,erhi,er,mat = get_hist_2d( rootfile, "E_migration" )
        except:
            etlo,ethi,et,erlo,erhi,er,mat = get_hist_2d( rootfile, "EMig" )

        cols.append( fits.Column( name="E_migration",
                                    format="{0}D".format(len(mat[0])),
                                    array=mat))

        tbhdu = fits.new_table( cols )
        tbhdu.name="SENS"



        nx = len(et)
        ny = len(er)
        mig = fitshistogram.Histogram(bins=[nx,ny],\
                                      range=[ [etlo[0],ethi[nx-1]] ,
                                              [erlo[0],erhi[ny-1]] ],\
                        axisNames=['logE_t','logE_r'], name="MIGRATION")
        mig.hist=np.array(mat)

        outputfilename = os.path.splitext(os.path.basename(rootfilename))[0]+".fits"
        print "WRITING:",outputfilename
        hdus = fits.HDUList( hdus=[fits.PrimaryHDU(),tbhdu,mig.asFITS()] )
        hdus.writeto(outputfilename, clobber=True)
