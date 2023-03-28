# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 12:50:55 2023

@author: Nikolaj Lange Dons
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#loading data as recarrays
hdu = fits.open('COSMOS2020_flux_mag_zphot.fits')
data = hdu[1].data


# Remove noise and unreliable values
for colname in data.names:
    data = data[data[colname] >= -98]
    data = data[data[colname] <= 300]

# Define central wavelength (Cwl) in Ångstrom[Å] from used bands/filters from cosmos2020
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

Lambda = np.array([Cwl_CFHT_u,Cwl_CFHT_ustar,Cwl_HSC_g,Cwl_HSC_r,Cwl_HSC_i,Cwl_HSC_z,Cwl_HSC_y,Cwl_UVISTA_Y,
          Cwl_UVISTA_J,Cwl_UVISTA_H,Cwl_UVISTA_Ks,Cwl_SC_IB427,Cwl_SC_IB464,Cwl_SC_IA484,Cwl_SC_IB505,
          Cwl_SC_IA527,Cwl_SC_IB574,Cwl_SC_IA624,Cwl_SC_IA679,Cwl_SC_IB709,Cwl_SC_IA738,Cwl_SC_IA767,
          Cwl_SC_IB827,Cwl_SC_NB711,Cwl_SC_NB816,Cwl_IRAC_CH1,Cwl_IRAC_CH2,Cwl_GALEX_FUV,Cwl_GALEX_NUV,Cwl_F814W])

#sorts lambda array from lowest to highest
def selection_sort(Lambda):
    for i in range(len(Lambda)):
        swap = i + np.argmin(Lambda[i:])
        (Lambda[i], Lambda[swap]) = (Lambda[swap], Lambda[i])
    return Lambda
Lambda = selection_sort(Lambda)

#converts lambda from ångstrøm to micrometer
Lambda = Lambda*10**-4

#matrix with 10 columns of galaxies and 30 rows with filters
G10_flux = np.array([data['GALEX_FUV_FLUX'][0:9],data['GALEX_NUV_FLUX'][0:9],
                     data['CFHT_u_FLUX'][0:9],data['CFHT_ustar_FLUX'][0:9],
                     data['SC_IB427_FLUX'][0:9],data['SC_IB464_FLUX'][0:9],
                     data['HSC_g_FLUX'][0:9],data['SC_IA484_FLUX'][0:9],
                     data['SC_IB505_FLUX'][0:9],data['SC_IA527_FLUX'][0:9],
                     data['SC_IB574_FLUX'][0:9],data['HSC_r_FLUX'][0:9],
                     data['SC_IA624_FLUX'][0:9],data['SC_IA679_FLUX'][0:9],
                     data['SC_IB709_FLUX'][0:9],data['SC_NB711_FLUX'][0:9],
                     data['SC_IA738_FLUX'][0:9],data['SC_IA767_FLUX'][0:9],
                     data['HSC_i_FLUX'][0:9],data['SC_NB816_FLUX'][0:9],
                     data['SC_IB827_FLUX'][0:9],data['F814W_FLUX'][0:9],
                     data['HSC_z_FLUX'][0:9],data['HSC_y_FLUX'][0:9],
                     data['UVISTA_Y_FLUX'][0:9],data['UVISTA_J_FLUX'][0:9],
                     data['UVISTA_H_FLUX'][0:9],data['UVISTA_Ks_FLUX'][0:9],
                     data['IRAC_CH1_FLUX'][0:9],data['IRAC_CH2_FLUX'][0:9]])

G10_fluxerr = np.array([data['GALEX_FUV_FLUXERR'][0:9],data['GALEX_NUV_FLUXERR'][0:9],
                        data['CFHT_u_FLUXERR'][0:9],data['CFHT_ustar_FLUXERR'][0:9],
                        data['SC_IB427_FLUXERR'][0:9],data['SC_IB464_FLUXERR'][0:9],
                        data['HSC_g_FLUXERR'][0:9],data['SC_IA484_FLUXERR'][0:9],
                        data['SC_IB505_FLUXERR'][0:9],data['SC_IA527_FLUXERR'][0:9],
                        data['SC_IB574_FLUXERR'][0:9],data['HSC_r_FLUXERR'][0:9],
                        data['SC_IA624_FLUXERR'][0:9],data['SC_IA679_FLUXERR'][0:9],
                        data['SC_IB709_FLUXERR'][0:9],data['SC_NB711_FLUXERR'][0:9],
                        data['SC_IA738_FLUXERR'][0:9],data['SC_IA767_FLUXERR'][0:9],
                        data['HSC_i_FLUXERR'][0:9],data['SC_NB816_FLUXERR'][0:9],
                        data['SC_IB827_FLUXERR'][0:9],data['F814W_FLUXERR'][0:9],
                        data['HSC_z_FLUXERR'][0:9],data['HSC_y_FLUXERR'][0:9],
                        data['UVISTA_Y_FLUXERR'][0:9], data['UVISTA_J_FLUXERR'][0:9],
                        data['UVISTA_H_FLUXERR'][0:9], data['UVISTA_Ks_FLUXERR'][0:9],
                        data['IRAC_CH1_FLUXERR'][0:9], data['IRAC_CH2_FLUXERR'][0:9]])

#change units from janskys to W/m^2
# G10_flux = G10_flux*10**-26
# G10_fluxerr = G10_fluxerr*10**-26

#initialize layout of galaxy 1
fig, ax = plt.subplots()

#add scatterplot and error bars of G1
ax.errorbar(Lambda, G10_flux[:,0], yerr = G10_fluxerr[:,0],fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "k", markersize = 2)

#set logarithmic scale on the x and y variable
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(10**-1, 10**1)
ax.set_ylim(10**-1, 10**2)
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Flux density (W/m2)")
ax.set_title("Flux Density against Wavelength of Galaxy 1")

#initialize layout of galaxy 2
fig, ax = plt.subplots()

#add scatterplot and error bars of G2
ax.errorbar(Lambda, G10_flux[:,1], yerr = G10_fluxerr[:,1],fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "k", markersize = 2)

#set logarithmic scale on the x and y variable
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(10**-1, 10**1)
ax.set_ylim(10**-1, 10**2)
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Flux density (W/m2)")
ax.set_title("Flux Density against Wavelength of Galaxy 2")

# Close FITS file so it won't use up excess memory
hdu.close()
