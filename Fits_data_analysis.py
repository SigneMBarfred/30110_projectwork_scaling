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
Cwl_GALEX_FUV = 1526
Cwl_GALEX_NUV = 2307
Cwl_CFHT_u = 3709
Cwl_CFHT_ustar = 3858
Cwl_SC_IB427 = 4266
Cwl_SC_IB464 = 4635
Cwl_HSC_g = 4847
Cwl_SC_IA484 = 4851
Cwl_SC_IB505 = 5064
Cwl_SC_IA527 = 5261
Cwl_SC_IB574 = 5766
Cwl_HSC_r = 6219
Cwl_SC_IA624 = 6232
Cwl_SC_IA679 = 6780
Cwl_SC_IB709 = 7073
Cwl_SC_NB711 = 7121
Cwl_SC_IA738 = 7361
Cwl_SC_IA767 = 7694
Cwl_HSC_i = 7699
Cwl_SC_NB816 = 8150
Cwl_SC_IB827 = 8243
Cwl_F814W = 8333
Cwl_HSC_z = 8894
Cwl_HSC_y = 9761
Cwl_UVISTA_Y = 10216
Cwl_UVISTA_J = 12525
Cwl_UVISTA_H = 16466
Cwl_UVISTA_Ks = 21557
Cwl_IRAC_CH1 = 35686
Cwl_IRAC_CH2 = 45067
Cwl_SPLASH_CH2 = 57788 # Tror der er sket en fejl | Det er nok IRAC_CH3
Cwl_SPLASH_CH3 = 79958 # Tror der er sket en fejl | Det er nok IRAC_CH4

Lambda = np.array([Cwl_CFHT_u,Cwl_CFHT_ustar,Cwl_HSC_g,Cwl_HSC_r,Cwl_HSC_i,Cwl_HSC_z,Cwl_HSC_y,Cwl_UVISTA_Y,
          Cwl_UVISTA_J,Cwl_UVISTA_H,Cwl_UVISTA_Ks,Cwl_SC_IB427,Cwl_SC_IB464,Cwl_SC_IA484,Cwl_SC_IB505,
          Cwl_SC_IA527,Cwl_SC_IB574,Cwl_SC_IA624,Cwl_SC_IA679,Cwl_SC_IB709,Cwl_SC_IA738,Cwl_SC_IA767,
          Cwl_SC_IB827,Cwl_SC_NB711,Cwl_SC_NB816,Cwl_IRAC_CH1,Cwl_IRAC_CH2,Cwl_SPLASH_CH2,Cwl_SPLASH_CH3,
          Cwl_GALEX_FUV,Cwl_GALEX_NUV,Cwl_F814W])

def selection_sort(Lambda):
    for i in range(len(Lambda)):
        swap = i + np.argmin(Lambda[i:])
        (Lambda[i], Lambda[swap]) = (Lambda[swap], Lambda[i])
    return Lambda
Lambda = selection_sort(Lambda)

    
#### Plot 1 scatter(wl,flux) of 10 galaxies



G10_flux = np.matrix([evt_data['GALEX_FUV_FLUX'][1:10],evt_data['GALEX_NUV_FLUX'][1:10],evt_data['CFHT_u_FLUX'][1:10],
                     evt_data['CFHT_ustar_FLUX'][1:10],evt_data['SC_IB427_FLUX'][1:10],evt_data['SC_IB464_FLUX'][1:10],
                     evt_data['HSC_g_FLUX'][1:10],evt_data['SC_IA484_FLUX'][1:10],evt_data['SC_IB505_FLUX'][1:10],
                     evt_data['SC_IA527_FLUX'][1:10],evt_data['SC_IB574_FLUX'][1:10],evt_data['HSC_r_FLUX'][1:10],
                     evt_data['SC_IA624_FLUX'][1:10],evt_data['SC_IA679_FLUX'][1:10],evt_data['SC_IB709_FLUX'][1:10],
                     evt_data['SC_NB711_FLUX'][1:10],evt_data['SC_IA738_FLUX'][1:10],evt_data['SC_IA767_FLUX'][1:10],
                     evt_data['HSC_i_FLUX'][1:10],evt_data['SC_NB816_FLUX'][1:10],evt_data['SC_IB827_FLUX'][1:10],
                     evt_data['F814W_FLUX'][1:10],evt_data['HSC_z_FLUX'][1:10],evt_data['HSC_y_FLUX'][1:10],
                     evt_data['UVISTA_Y_FLUX'][1:10],evt_data['UVISTA_J_FLUX'][1:10],evt_data['UVISTA_H_FLUX'][1:10],
                     evt_data['UVISTA_Ks_FLUX'][1:10],evt_data['IRAC_CH1_FLUX'][1:10],evt_data['IRAC_CH2_FLUX'][1:10],
                     evt_data['SPLASH_CH2_FLUX'][1:10],evt_data['SPLASH_CH3_FLUX'][1:10]])


G10_fluxerr = np.matrix([evt_data['GALEX_FUV_FLUXERR'][1:10],evt_data['GALEX_NUV_MAGERR'][1:10],
                        evt_data['CFHT_u_FLUXERR'][1:10],evt_data['CFHT_ustar_FLUXERR'][1:10],
                        evt_data['SC_IB427_FLUXERR'][1:10],evt_data['SC_IB464_FLUXERR'][1:10],
                        evt_data['HSC_g_FLUXERR'][1:10],evt_data['SC_IA484_FLUXERR'][1:10],
                        evt_data['SC_IB505_FLUXERR'][1:10],evt_data['SC_IA527_FLUXERR'][1:10],
                        evt_data['SC_IB574_FLUXERR'][1:10],evt_data['HSC_r_FLUXERR'][1:10],
                        evt_data['SC_IA624_FLUXERR'][1:10],evt_data['SC_IA679_FLUXERR'][1:10],
                        evt_data['SC_IB709_FLUXERR'][1:10],evt_data['SC_NB711_FLUXERR'][1:10],
                        evt_data['SC_IA738_FLUXERR'][1:10],evt_data['SC_IA767_FLUXERR'][1:10],
                        evt_data['HSC_i_FLUXERR'][1:10],evt_data['SC_NB816_FLUXERR'][1:10],
                        evt_data['SC_IB827_FLUXERR'][1:10],evt_data['F814W_FLUXERR'][1:10],evt_data['HSC_z_FLUXERR'][1:10],
                        evt_data['UVISTA_Y_FLUXERR'][1:10], evt_data['UVISTA_J_FLUXERR'][1:10],
                        evt_data['UVISTA_H_FLUXERR'][1:10], evt_data['UVISTA_Ks_FLUXERR'][1:10],
                        evt_data['IRAC_CH1_FLUXERR'][1:10], evt_data['IRAC_CH2_FLUXERR'][1:10],
                        evt_data['SPLASH_CH4_FLUXERR'][1:10], evt_data['SPLASH_CH4_FLUXERR'][1:10]])
# SPLASH_CH4 Burde nok fjernes der mangler en fluxerr


plt.scatter(Lambda, G10_flux)

plt.errorbar(Lambda, G10_flux, yerr=G10_fluxerr, fmt="o")
plt.show()


# Create array with flux from bands in order from cosmos2020
Band_flux = 1

# Close FITS file so it won't use up excess memory
hdu_list.close()