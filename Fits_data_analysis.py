# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 12:50:55 2023

@author: Nikolaj Lange Dons
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from matplotlib.colors import LogNorm

import pandas as pd
import matplotlib.pyplot as plt

# Set up matplotlib
import matplotlib.pyplot as plt
#%matplotlib inline

# Load file
# open it with memmap=True to prevent RAM storage issues.
hdu_list = fits.open('COSMOS2020_flux_mag_zphot.fits')
hdu_list.info()

# Print column names
print(hdu_list[1].columns)


# Convert data in to an astropy.table - Because it's more convenient
evt_data = hdu_list[1].data

# Display table
evt_data

# Panda table

# Remove noise and unreliable values (-99, >300)
evt_data = evt_data[evt_data != -99]
#evt_data.remove_rows([]


# If statement to remove rows with higher values than 300


# m = evt_data > 300
# evt_data = evt_data(m)

# Defin central wavelength (Cwl) in Ångstrom[Å] from used bands from cosmos2020
Cwl_CFHT_u = 3709
Cwl_CFHT_ustar = 3858
Cwl_HSC_g = 4847
Cwl_HSC_r = 6219
Cwl_HSC_i = 7699
Cwl_HSC_z = 8894
Cwl_HSC_y = 9761
Cwl_UVISTA_Y = 10216
Cwl_UVISTA_J = 12525
Cwl_UVISTA_H = 16466
Cwl_UVISTA_Ks = 21557
Cwl_SC_IB427 = 4266
Cwl_SC_IB464 = 4635
Cwl_SC_IA484 = 4851
Cwl_SC_IB505 = 5064
Cwl_SC_IA527 = 5261
Cwl_SC_IB574 = 5766
Cwl_SC_IA624 = 6232
Cwl_SC_IA679 = 6780
Cwl_SC_IB709 = 7073
Cwl_SC_IA738 = 7361
Cwl_SC_IA767 = 7694
Cwl_SC_IB827 = 8243
Cwl_SC_NB711 = 7121
Cwl_SC_NB816 = 8150
Cwl_IRAC_CH1 = 35686
Cwl_IRAC_CH2 = 45067
Cwl_SPLASH_CH2 = 57788 # Tror der er sket en fejl | Det er nok IRAC_CH3
Cwl_SPLASH_CH3 = 79958 # Tror der er sket en fejl | Det er nok IRAC_CH4
Cwl_GALEX_FUV = 1526
Cwl_GALEX_NUV = 2307
Cwl_F814W = 8333

Lambda = np.array([Cwl_CFHT_u,Cwl_CFHT_ustar,Cwl_HSC_g,Cwl_HSC_r,Cwl_HSC_i,Cwl_HSC_z,Cwl_HSC_y,Cwl_UVISTA_Y,
          Cwl_UVISTA_J,Cwl_UVISTA_H,Cwl_UVISTA_Ks,Cwl_SC_IB427,Cwl_SC_IB464,Cwl_SC_IA484,Cwl_SC_IB505,
          Cwl_SC_IA527,Cwl_SC_IB574,Cwl_SC_IA624,Cwl_SC_IA679,Cwl_SC_IB709,Cwl_SC_IA738,Cwl_SC_IA767,
          Cwl_SC_IB827,Cwl_SC_NB711,Cwl_SC_NB816,Cwl_IRAC_CH1,Cwl_IRAC_CH2,Cwl_SPLASH_CH2,Cwl_SPLASH_CH3,
          Cwl_GALEX_FUV,Cwl_GALEX_NUV,Cwl_F814W])


#### Plot 1 scatter(wl,flux) of 10 galaxies

G10_flux = np.array([evt_data['CFHT_u_FLUX'][1:10],evt_data['CFHT_g_FLUX'][1:10],evt_data['CFHT_u_FLUX'][1:10]])

G10_fluxerr = np.array([evt_data['GALEX_FUV_FLUXERR'][1:10],evt_data['GALEX_NUV_MAGERR'][1:10],evt_data['CFHT_u_FLUXERR'][1:10],evt_data['CFHT_g_FLUXERR'][1:10],evt_data['CFHT_r_FLUXERR'][1:10],
                        evt_data['CFHT_i_FLUXERR'][1:10],evt_data['CFHT_z_FLUXERR'][1:10],evt_data['CFHT_y_FLUXERR'][1:10],
                        evt_data['CFHT_y_FLUXERR'][1:10],])

#plt.scatter(a, b)
  
c = [1, 3, 2, 1]
  
#plt.errorbar(a, b, yerr=c, fmt="o")
#plt.show()



# Create array with flux from bands in order from cosmos2020
Band_flux = 1

# Close FITS file so it won't use up excess memory
hdu_list.close()