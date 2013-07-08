Data Format
===========

The inputs are FITS tables with the following format  ::

        Name               Type       Dimensions
        ----               ----       ----------
 HDU 1   Primary Array      Null Array                               
 HDU 2   SENS               BinTable    11 cols x 20 rows            

   Col  Name             Format[Units](Range)      Comment
     1 LOG10_E_LO         D     log10TeV           
     2 LOG10_E_HI         D     log10TeV       
     3 N_detected         D     counts             vs E_reco
     4 N_simulated        D     counts             vs E_true
     5 r_simulated        D     deg                vs E_true
     6 ThetaSqr           D     deg                vs E_reco
     7 phi_diffuse        D     deg                vs E_true
     8 R68_psf            D     deg                vs E_reco
     9 E_res              D                    
    10 E_bias             D                    
    11 E_migration        20D                      
