# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 10:47:56 2023

@author: Nikolaj Lange Dons
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Column
import matplotlib.pyplot as plt
from astropy.table import Table
import itertools
import csv
import os

# Set directory
os.chdir('C://Users//Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet\Dokumenter/4 semester/Fagprojekt/Data behandling')

#loading data as recarrays
hdu = fits.open('filtered_SED.fits',memmap=True)
data = hdu[1].data

# Filter for U, V and J band
Uband = np.genfromtxt('CFHT_CFH12k.B.dat') #skip_header=1, skip_footer=1, names=True, dtype=None, delimiter=' '
Vband = np.genfromtxt('CFHT_CFH12k.R.dat')
Jband = np.genfromtxt('2MASS_2MASS.J.dat')

# Load fits models
Fit1 = pd.read_csv ('fits_models.csv')
Fit = pd.DataFrame(Fit1)

##############################
# Command room # Choose everything here

# Choose redshift interval
lower_z = 0.1
upper_z = 0.3

# Choose galaxy from catalouge (max = 511005)
galaxy = 1
Gal = galaxy+1

# Choose galaxy fit model (Choose between 1:11)
Fit_model1 = 1
Fit_model2 = 2
Fit_model3 = 3

# maximum of filter # Choose value to fit under SED's
flux_min = 3*10**-2;

##############################

# Format the different filters to fit in the plot
Uband[:,1] = 10**Uband[:,1]
Vband[:,1] = 10**Vband[:,1]
Jband[:,1] = 10**Jband[:,1]

Uband_max = np.max(Uband[:,1]);
Vband_max = np.max(Vband[:,1]);
Jband_max = np.max(Jband[:,1]);

Uband[:,1] = Uband[:,1] / Uband_max #Normalize
Vband[:,1] = Vband[:,1] / Vband_max
Jband[:,1] = Jband[:,1] / Jband_max

# Match wl to correct microns
Uband[:,0] = Uband[:,0]/10000
Vband[:,0] = Vband[:,0]/10000
Jband[:,0] = Jband[:,0]/10000

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

Lambda = np.array([Cwl_CFHT_u,Cwl_CFHT_ustar,Cwl_HSC_g,Cwl_HSC_r,Cwl_HSC_i,Cwl_HSC_z,Cwl_HSC_y,Cwl_UVISTA_Y,
          Cwl_UVISTA_J,Cwl_UVISTA_H,Cwl_UVISTA_Ks,Cwl_SC_IB427,Cwl_SC_IB464,Cwl_SC_IA484,Cwl_SC_IB505,
          Cwl_SC_IA527,Cwl_SC_IB574,Cwl_SC_IA624,Cwl_SC_IA679,Cwl_SC_IB709,Cwl_SC_IA738,Cwl_SC_IA767,
          Cwl_SC_IB827,Cwl_SC_NB711,Cwl_SC_NB816,Cwl_IRAC_CH1,Cwl_IRAC_CH2,Cwl_GALEX_NUV])

def selection_sort(Lambda):
    for i in range(len(Lambda)):
        swap = i + np.argmin(Lambda[i:])
        (Lambda[i], Lambda[swap]) = (Lambda[swap], Lambda[i])
    return Lambda
Lambda = selection_sort(Lambda)

#converts lambda from ångstrøm to micrometer
Lambda = Lambda*10**-4

#matrix with flux from all filters
Gflux = np.array([data['GALEX_NUV_FLUX'],data['CFHT_u_FLUX'],
                     data['CFHT_ustar_FLUX'],data['SC_IB427_FLUX'],
                     data['SC_IB464_FLUX'],data['HSC_g_FLUX'],
                     data['SC_IA484_FLUX'],data['SC_IB505_FLUX'],
                     data['SC_IA527_FLUX'],data['SC_IB574_FLUX'],
                     data['HSC_r_FLUX'],data['SC_IA624_FLUX'],
                     data['SC_IA679_FLUX'],data['SC_IB709_FLUX'],
                     data['SC_NB711_FLUX'],data['SC_IA738_FLUX'],
                     data['SC_IA767_FLUX'],data['HSC_i_FLUX'],
                     data['SC_NB816_FLUX'],data['SC_IB827_FLUX'],
                     data['HSC_z_FLUX'],data['HSC_y_FLUX'],
                     data['UVISTA_Y_FLUX'],data['UVISTA_J_FLUX'],
                     data['UVISTA_H_FLUX'],data['UVISTA_Ks_FLUX'],
                     data['IRAC_CH1_FLUX'],data['IRAC_CH2_FLUX']])
Gflux.reshape(28, 511006)

# Normalize around Ks_filter
Gflux = Gflux/data['UVISTA_Ks_FLUX']

# Do the same for flux_err
Gfluxerr = np.array([data['GALEX_NUV_FLUXERR'],data['CFHT_u_FLUXERR'],
                        data['CFHT_ustar_FLUXERR'],data['SC_IB427_FLUXERR'],
                        data['SC_IB464_FLUXERR'],data['HSC_g_FLUXERR'],
                        data['SC_IA484_FLUXERR'],data['SC_IB505_FLUXERR'],
                        data['SC_IA527_FLUXERR'],data['SC_IB574_FLUXERR'],
                        data['HSC_r_FLUXERR'],data['SC_IA624_FLUXERR'],
                        data['SC_IA679_FLUXERR'],data['SC_IB709_FLUXERR'],
                        data['SC_NB711_FLUXERR'],data['SC_IA738_FLUXERR'],
                        data['SC_IA767_FLUXERR'],data['HSC_i_FLUXERR'],
                        data['SC_NB816_FLUXERR'],data['SC_IB827_FLUXERR'],
                        data['HSC_z_FLUXERR'],data['HSC_y_FLUXERR'],
                        data['UVISTA_Y_FLUXERR'], data['UVISTA_J_FLUXERR'],
                        data['UVISTA_H_FLUXERR'], data['UVISTA_Ks_FLUXERR'],
                        data['IRAC_CH1_FLUXERR'], data['IRAC_CH2_FLUXERR']])
Gfluxerr.reshape(28, 511006)

Gfluxerr = Gfluxerr/data['UVISTA_Ks_FLUX']

# Normalize fit model to K band
Fit.iloc[:,[1,2,3,4,5,6,7,8,9,10,11,12]] = Fit.iloc[:,[1,2,3,4,5,6,7,8,9,10,11,12]]/data['UVISTA_Ks_FLUX'][0,galaxy]


#######################################

#initialize layout of plot1 
select_cat1_flux = np.logical_and((data['z_phot'] > lower_z), (data['z_phot'] < upper_z))

Gflux_cat1 = Gflux[:,0,select_cat1_flux[0,:]]
Gfluxerr_cat1 = Gfluxerr[:,0,select_cat1_flux[0,:]]

# Format fits model to z = (upper_z+lower_z)/2 # mean of interval
cat1_z = (upper_z+lower_z)/2
wl_cat1 = Fit['wl']+Fit['wl']*cat1_z

fig1, ax = plt.subplots()

# Add scatterplot and error bars of choosen galaxy from catalouge
Galaxy = ax.errorbar(Lambda, Gflux_cat1[:,galaxy], yerr = Gfluxerr_cat1[:,galaxy],fmt = "o", alpha = 0.8,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "c", markersize = 2)


# Add fits model from our first catalouge
Fitting1_txt = "Fit model nr. %d"
Fitting1 = ax.plot(Fit['wl'], Fit.iloc[:,Fit_model1], 'forestgreen', alpha = 0.6, label=(Fitting1_txt%Fit_model1))

# Fitting2_txt = "Fit model nr. %d"
# Fitting2 = ax.plot(Fit['wl'], Fit.iloc[:,Fit_model2], 'limegreen', alpha = 0.6, label=(Fitting2_txt%Fit_model2))

# Fitting3_txt = "Fit model nr. %d"
# Fitting3 = ax.plot(Fit['wl'], Fit.iloc[:,Fit_model3], 'darkgreen', alpha = 0.6, label=(Fitting3_txt%Fit_model3))

# Add Uband, Vband and Jband 
U1 = ax.plot(Uband[:,0],Uband[:,1]*flux_min, 'b', alpha = 0.5) # label="U band"
V1 = ax.plot(Vband[:,0],Vband[:,1]*flux_min, 'k', alpha = 0.5) # label="V band"
J1 = ax.plot(Jband[:,0],Jband[:,1]*flux_min, 'r', alpha = 0.5) # label="J band"
Band = []
Band.append([U1, V1, J1])

#set logarithmic scale on the x and y variable
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(10**-1, 10**1)
ax.set_ylim(10**-2, 10**1)
ax.set_xlabel("Wavelength (µm)")
ax.set_ylabel("Flux density")
fig1_txt = "Flux Density for Galaxy %d at %.1f < z < %.1f"
ax.set_title(fig1_txt%(Gal,lower_z,upper_z))
legend1_txt = "Catalouge galaxy nr. %d"
# ax.legend((Galaxy),('Catalouge galaxy nr. '), scatterpoints=1, 
#           loc='lower right', ncol=1, fontsize=8)
leg = ax.legend(loc='lower right', ncol=1);

#########################

# #initialize layout of plot3 containing 3 galaxies [red, yellow, blue] in 0.9 < z < 1.1
# select_cat2_flux = np.logical_and((data['z_phot'] > 0.9), (data['z_phot'] < 1.1))

# Gflux_cat2 = Gflux[:,0,select_cat2_flux[0,:]]
# Gfluxerr_cat2 = Gfluxerr[:,0,select_cat2_flux[0,:]]

# fig3, ax = plt.subplots()

# #add scatterplot and error bars of G7
# G7 = ax.errorbar(Lambda, Gflux_cat2[:,6], yerr = Gfluxerr_cat2[:,6],fmt = "o", alpha = 0.7,
#             elinewidth = 0.7, capsize = 5, markeredgecolor = "r", markersize = 2)

# # #add scatterplot and error bars of G8
# G8 = ax.errorbar(Lambda, Gflux_cat2[:,7], yerr = Gfluxerr_cat2[:,7],fmt = "o", alpha = 0.7,
#             elinewidth = 0.7, capsize = 5, markeredgecolor = "y", markersize = 2)

# # #add scatterplot and error bars of G9
# G9 = ax.errorbar(Lambda, Gflux_cat2[:,8], yerr = Gfluxerr_cat2[:,8],fmt = "o", alpha = 0.7,
#             elinewidth = 0.7, capsize = 5, markeredgecolor = "b", markersize = 2)

# # Add Uband, Vband and Jband 
# U3 = ax.plot(Uband[:,0],Uband[:,1]*flux_min, 'b', alpha = 0.5)
# V3 = ax.plot(Vband[:,0],Vband[:,1]*flux_min, 'k', alpha = 0.5)
# J3 = ax.plot(Jband[:,0],Jband[:,1]*flux_min, 'r', alpha = 0.5)

# #set logarithmic scale on the x and y variable
# ax.set_xscale("log")
# ax.set_yscale("log")
# ax.set_xlim(10**-1, 10**1)
# ax.set_ylim(10**-2, 10**1)
# ax.set_xlabel("Wavelength (µm)")
# ax.set_ylabel("Flux density")
# ax.legend((G7, G8, G9, U3, V3, J3),
#            ('Galaxy 7', 'Galaxy 8', 'Galaxy 9', 'U band', 'V band', 'J band'),
#            scatterpoints=1,
#            loc='lower right',
#            ncol=1,
#            fontsize=8)
# ax.set_title("Flux Density against Wavelength for Galaxy 7,8,9 at 0.9 < z < 1.1")

# #########################

# #initialize layout of plot4 containing 3 galaxies [red, yellow, blue] in 1.9 < z < 2.1
# select_cat3_flux = np.logical_and((data['z_phot'] > 1.9), (data['z_phot'] < 2.1))

# Gflux_cat3 = Gflux[:,0,select_cat3_flux[0,:]]
# Gfluxerr_cat3 = Gfluxerr[:,0,select_cat3_flux[0,:]]

# fig4, ax = plt.subplots()

# #add scatterplot and error bars of G10
# G10 = ax.errorbar(Lambda, Gflux_cat3[:,9], yerr = Gfluxerr_cat3[:,9],fmt = "o", alpha = 0.7,
#             elinewidth = 0.7, capsize = 5, markeredgecolor = "r", markersize = 2)

# # #add scatterplot and error bars of G11
# G11 = ax.errorbar(Lambda, Gflux_cat3[:,10], yerr = Gfluxerr_cat3[:,10],fmt = "o", alpha = 0.7,
#             elinewidth = 0.7, capsize = 5, markeredgecolor = "y", markersize = 2)

# # #add scatterplot and error bars of G12
# G12 = ax.errorbar(Lambda, Gflux_cat3[:,11], yerr = Gfluxerr_cat3[:,11],fmt = "o", alpha = 0.7,
#             elinewidth = 0.7, capsize = 5, markeredgecolor = "b", markersize = 2)

# # Add Uband, Vband and Jband 
# U4 = ax.plot(Uband[:,0],Uband[:,1]*flux_min, 'b', alpha = 0.5)
# V4 = ax.plot(Vband[:,0],Vband[:,1]*flux_min, 'k', alpha = 0.5)
# J4 = ax.plot(Jband[:,0],Jband[:,1]*flux_min, 'r', alpha = 0.5)

# #set logarithmic scale on the x and y variable
# ax.set_xscale("log")
# ax.set_yscale("log")
# ax.set_xlim(10**-1, 10**1)
# ax.set_ylim(10**-2, 10**1)
# ax.set_xlabel("Wavelength (µm)")
# ax.set_ylabel("Flux density")
# ax.legend((G10, G11, G12, U4, V4, J4),
#            ('Galaxy 10', 'Galaxy 11', 'Galaxy 12', 'U band', 'V band', 'J band'),
#            scatterpoints=1,
#            loc='lower right',
#            ncol=1,
#            fontsize=8)
# ax.set_title("Flux Density against Wavelength for Galaxy 10,11,12 at 1.9 < z < 2.1")

# Close FITS file so it won't use up excess memory
hdu.close()
