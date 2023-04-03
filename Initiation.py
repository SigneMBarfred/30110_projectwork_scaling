# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 19:53:48 2023

@author: Nikolaj Lange Dons
"""

from astropy.table import Table
import numpy as np
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt

#loading data as recarrays
hdu = fits.open('COSMOS2020_best.fits',memmap=True)
data = hdu[1].data

# create a sample table
t = Table(data)  # convert to masked table

# convert to numpy array
tnew = Table.to_pandas(t)
SED = pd.DataFrame(tnew).to_numpy()

# Remove noise and unreliable values
SED = SED[SED[:,tnew.columns.get_loc('mass')] > 10**8]
SED = SED[SED[:,tnew.columns.get_loc('mass')] < 10**13]
SED = SED[SED[:,tnew.columns.get_loc('z_phot')] > 0]
SED = SED[SED[:,(tnew.columns.get_loc('UVISTA_Ks_FLUX')-1)] != -99]
SED[SED == -99] = np.nan


# Define central wavelength (Cwl) in Ångstrom[Å] from used bands/filters from cosmos2020
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
Cwl_HSC_z = 8894
Cwl_HSC_y = 9761
Cwl_UVISTA_Y = 10216
Cwl_UVISTA_J = 12525
Cwl_UVISTA_H = 16466
Cwl_UVISTA_Ks = 21557
Cwl_IRAC_CH1 = 35686
Cwl_IRAC_CH2 = 45067

Lambda = np.array([Cwl_GALEX_NUV, Cwl_CFHT_u, Cwl_CFHT_ustar, Cwl_SC_IB427, Cwl_SC_IB464, Cwl_HSC_g,
     Cwl_SC_IA484,Cwl_SC_IB505, Cwl_SC_IA527, Cwl_SC_IB574, Cwl_HSC_r, Cwl_SC_IA624,
     Cwl_SC_IA679, Cwl_SC_IB709, Cwl_SC_NB711, Cwl_SC_IA738, Cwl_SC_IA767, Cwl_HSC_i,
     Cwl_SC_NB816, Cwl_SC_IB827, Cwl_HSC_z, Cwl_HSC_y, Cwl_UVISTA_Y, Cwl_UVISTA_J, 
     Cwl_UVISTA_H, Cwl_UVISTA_Ks, Cwl_IRAC_CH1, Cwl_IRAC_CH2])

# #converts lambda from ångstrøm to micrometer
Lambda = Lambda*10**-4

#matrix with 10 columns of galaxies and 30 rows with filters
G10_flux = np.array([[SED[:,tnew.columns.get_loc('GALEX_NUV_FLUXERR')]],[SED[:,tnew.columns.get_loc('CFHT_u_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('CFHT_ustar_FLUX')]],[SED[:,tnew.columns.get_loc('SC_IB427_FLUX')]],
        [SED[:,tnew.columns.get_loc('SC_IB464_FLUX')]],[SED[:,tnew.columns.get_loc('HSC_g_FLUX')]],
        [SED[:,tnew.columns.get_loc('SC_IA484_FLUX')]],[SED[:,tnew.columns.get_loc('SC_IB505_FLUX')]],
        [SED[:,tnew.columns.get_loc('SC_IA527_FLUX')]],[SED[:,tnew.columns.get_loc('SC_IB574_FLUX')]],
        [SED[:,tnew.columns.get_loc('HSC_r_FLUX')]],[SED[:,tnew.columns.get_loc('SC_IA624_FLUX')]],
        [SED[:,tnew.columns.get_loc('SC_IA679_FLUX')]],[SED[:,tnew.columns.get_loc('SC_IB709_FLUX')]],
        [SED[:,tnew.columns.get_loc('SC_NB711_FLUX')]],[SED[:,tnew.columns.get_loc('SC_IA738_FLUX')]],
        [SED[:,tnew.columns.get_loc('SC_IA767_FLUX')]],[SED[:,tnew.columns.get_loc('HSC_i_FLUX')]],
        [SED[:,tnew.columns.get_loc('SC_NB816_FLUX')]],[SED[:,tnew.columns.get_loc('SC_IB827_FLUX')]],
        [SED[:,tnew.columns.get_loc('HSC_z_FLUX')]],[SED[:,tnew.columns.get_loc('HSC_y_FLUX')]],
        [SED[:,tnew.columns.get_loc('UVISTA_Y_FLUX')]],[SED[:,tnew.columns.get_loc('UVISTA_J_FLUX')]],
        [SED[:,tnew.columns.get_loc('UVISTA_H_FLUX')]],[SED[:,tnew.columns.get_loc('UVISTA_Ks_FLUX')]],
        [SED[:,tnew.columns.get_loc('IRAC_CH1_FLUX')]],[SED[:,tnew.columns.get_loc('IRAC_CH2_FLUX')]]])

G10_flux = G10_flux.reshape(len(SED[:,1]),1,len(Lambda))


G10_flux = G10_flux/np.array([SED[:,tnew.columns.get_loc('GALEX_NUV_FLUXERR')]])

G10_fluxerr = np.array([[SED[:,tnew.columns.get_loc('GALEX_NUV_FLUXERR')]],[SED[:,tnew.columns.get_loc('CFHT_u_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('CFHT_ustar_FLUXERR')]],[SED[:,tnew.columns.get_loc('SC_IB427_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('SC_IB464_FLUXERR')]],[SED[:,tnew.columns.get_loc('HSC_g_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('SC_IA484_FLUXERR')]],[SED[:,tnew.columns.get_loc('SC_IB505_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('SC_IA527_FLUXERR')]],[SED[:,tnew.columns.get_loc('SC_IB574_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('HSC_r_FLUXERR')]],[SED[:,tnew.columns.get_loc('SC_IA624_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('SC_IA679_FLUXERR')]],[SED[:,tnew.columns.get_loc('SC_IB709_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('SC_NB711_FLUXERR')]],[SED[:,tnew.columns.get_loc('SC_IA738_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('SC_IA767_FLUXERR')]],[SED[:,tnew.columns.get_loc('HSC_i_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('SC_NB816_FLUXERR')]],[SED[:,tnew.columns.get_loc('SC_IB827_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('HSC_z_FLUXERR')]],[SED[:,tnew.columns.get_loc('HSC_y_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('UVISTA_Y_FLUXERR')]], [SED[:,tnew.columns.get_loc('UVISTA_J_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('UVISTA_H_FLUXERR')]], [SED[:,tnew.columns.get_loc('UVISTA_Ks_FLUXERR')]],
        [SED[:,tnew.columns.get_loc('IRAC_CH1_FLUXERR')]], [SED[:,tnew.columns.get_loc('IRAC_CH2_FLUXERR')]]])

G10_fluxerr = G10_fluxerr.reshape(len(SED[:,1]),1,len(Lambda))


# #initialize layout of galaxy 1
fig, ax = plt.subplots()

#add scatterplot and error bars of G1
ax.errorbar(Lambda, G10_flux[:,0,0], yerr = G10_fluxerr[:,0,0],fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "k", markersize = 2)

#set logarithmic scale on the x and y variable
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(10**-1, 10**1)
ax.set_ylim(10**-6, 10**2)
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Flux density")
ax.set_title("Flux Density against Wavelength of Galaxy 1")

#initialize layout of galaxy 2
fig, ax = plt.subplots()

#add scatterplot and error bars of G2
ax.errorbar(Lambda, G10_flux[:,0,1], yerr = G10_fluxerr[:,0,1],fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "k", markersize = 2)


#set logarithmic scale on the x and y variable
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(10**-1, 10**1)
ax.set_ylim(10**-6, 10**2)
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Flux density")
ax.set_title("Flux Density against Wavelength of Galaxy 2")

#initialize layout of G3
fig, ax = plt.subplots()

#add scatterplot and error bars of G3
ax.errorbar(Lambda, G10_flux[:,0,4], yerr = G10_fluxerr[:,0,4],fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "k", markersize = 2)

#set logarithmic scale on the x and y variable
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(10**-1, 10**1)
ax.set_ylim(10**-6, 10**2)
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Flux density")
ax.set_title("Flux Density against Wavelength of Galaxy 3")

# Close FITS file so it won't use up excess memory
hdu.close()

