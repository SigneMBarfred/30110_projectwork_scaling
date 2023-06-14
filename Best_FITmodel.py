# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 16:21:53 2023

@author: Nikolaj Lange Dons & Signe Barfred
"""

import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt

#loading data as recarrays
hdu = fits.open('filtered_SED.fits',memmap=True)
data = hdu[1].data

# Filter for U, V and J band
Uband = np.genfromtxt('CFHT_CFH12k.B.dat') #skip_header=1, skip_footer=1, names=True, dtype=None, delimiter=' '
Vband = np.genfromtxt('CFHT_CFH12k.R.dat')
Jband = np.genfromtxt('2MASS_2MASS.J.dat')

# Load fits models
Fit1 = pd.read_csv ('fits_models.csv')
df = pd.DataFrame(Fit1)
Fit = df.to_numpy()

##############################
# Command room # Choose variables here

# Choose redshift interval
lower_z = 2
upper_z = 2.5

# Choose galaxy from catalogue (max = 511005)
galaxy1 = 14
galaxy2 = 6427
galaxy3 = 3

Gal1 = galaxy1+1
Gal2 = galaxy2+1
Gal3 = galaxy3+1

# Choose galaxy fit model (Choose between 1:11)
Fit_model1 = 2
Fit_model2 = 3
Fit_model3 = 3

# maximum of filter # Choose value to fit under SED's
flux_min = 5*10**-2;

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

#sort filters from lowest central wavelength to largest central wavelength
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
Fit_Ks_Index = np.logical_and((Fit[:,0] > 21457*10**-4), (Fit[:,0] < 21657*10**-4))
Fit_Norm = Fit[Fit_Ks_Index,Fit_model1]
Norm_faktor = np.mean(Fit_Norm)

Fit[:,1:12]= Fit[:,1:12]/Norm_faktor

#######################################

#initialize layout of plot1 
#
select_cat1_flux = np.logical_and((data['z_phot'] > lower_z), (data['z_phot'] < upper_z))

Gflux_cat1 = Gflux[:,0,select_cat1_flux[0,:]]
Gfluxerr_cat1 = Gfluxerr[:,0,select_cat1_flux[0,:]]

# Determine redshift for the galaxy chosen
z = data['z_phot']
z = z[:,select_cat1_flux[0,:]]
Galaxy_z = z[:,galaxy1]
z_txt = "z = %0.2f"

# Format fits model to z = (upper_z+lower_z)/2 # mean of interval
cat1_z = (upper_z+lower_z)/2

## we shift the template model to fit the galaxy's redshift
#1 + z = fit_shifted / fit_z0
#fit_shifted = fit_z0 + fit_z0 * z
wl_cat1 = Fit[:,0]+Fit[:,0]*Galaxy_z

ub_cat1 = Uband[:,0]+Uband[:,0]*Galaxy_z
vb_cat1 = Vband[:,0]+Vband[:,0]*Galaxy_z
jb_cat1 = Jband[:,0]+Jband[:,0]*Galaxy_z

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 6), gridspec_kw={'height_ratios': [3, 1]})
# Add scatterplot and error bars of choosen galaxy from catalouge
Galaxy_txt = "Galaxy nr. %d"
Galaxy1 = ax1.errorbar(Lambda, Gflux_cat1[:,galaxy1], yerr = Gfluxerr_cat1[:,galaxy1],fmt = "o", alpha = 0.8,
            elinewidth = 0.7, capsize = 5, mfc ='black', ecolor = 'black', markeredgecolor = "black", markersize = 3.5, label=Galaxy_txt%galaxy1)
# Add fits model from our first catalouge
Fitting1_txt = "Fit model nr. %d"
Fitting1 = ax1.plot(wl_cat1, Fit[:,Fit_model1], 'y', color = 'blue', lw = 0.75, alpha = 0.75, label=(Fitting1_txt%Fit_model1))

Fitting2_txt = "Fit model nr. %d"
Fitting2 = ax1.plot(wl_cat1, Fit[:,Fit_model2], 'y', color = 'cornflowerblue', lw = 0.95, alpha = 0.75, label=(Fitting2_txt%Fit_model2))

#colors of templates are set to match the colorcoding we did when illustrating templates in matlab

# Add Uband, Vband and Jband and scale to the flux
U1 = ax2.plot(ub_cat1,Uband[:,1]*flux_min, 'b', alpha = 0.5) # label="U band"
V1 = ax2.plot(vb_cat1,Vband[:,1]*flux_min, 'k', alpha = 0.5) # label="V band"
J1 = ax2.plot(jb_cat1,Jband[:,1]*flux_min, 'r', alpha = 0.5) # label="J band"
ax2.fill_between(ub_cat1,Uband[:,1]*flux_min, 'b', alpha = 0.3, label ='U-Band')
ax2.fill_between(vb_cat1,Vband[:,1]*flux_min, color = 'gray', alpha = 0.3, label = 'V-Band')
ax2.fill_between(jb_cat1,Jband[:,1]*flux_min, color = 'red', alpha = 0.3, label = 'J-Band')
Band = []
Band.append([U1, V1, J1])

#set logarithmic scale on the x and y variable
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlim(3*10**-1, 10**1)
ax1.set_ylim(0.3*10**-2, 0.55*10**1)
ax2.set_ylim(0.003, 0.07)

ax1.text(0.05,0.95,z_txt%Galaxy_z,ha='left',va='top',transform=ax1.transAxes,fontsize=10)
ax2.set_xlabel("Wavelength [µm] at "+ z_txt%Galaxy_z)
ax1.set_ylabel("Flux density [µJy]")
ax2.set_ylabel("")
ax2.set_yticklabels([])
ax2.legend()

ax1.set_title("Photometric flux data & recommended SED-template \n",fontsize = 14)
legend1_txt = "Catalouge galaxy nr. %d"
leg = ax1.legend(loc='lower right', ncol=1);


plt.show()
####################
# Residual plot
# How good is the fit
Index1 = np.logical_and((Fit[:,0] > Lambda[0]-300*10**-6), (Fit[:,0] < Lambda[0]+300*10**-6))
Index2 = np.logical_and((Fit[:,0] > Lambda[1]-300*10**-6), (Fit[:,0] < Lambda[1]+300*10**-6))
Index3 = np.logical_and((Fit[:,0] > Lambda[2]-300*10**-6), (Fit[:,0] < Lambda[2]+300*10**-6))
Index4 = np.logical_and((Fit[:,0] > Lambda[3]-300*10**-6), (Fit[:,0] < Lambda[3]+300*10**-6))
Index5 = np.logical_and((Fit[:,0] > Lambda[4]-300*10**-6), (Fit[:,0] < Lambda[4]+300*10**-6))
Index6 = np.logical_and((Fit[:,0] > Lambda[5]-300*10**-6), (Fit[:,0] < Lambda[5]+300*10**-6))
Index7 = np.logical_and((Fit[:,0] > Lambda[6]-300*10**-6), (Fit[:,0] < Lambda[6]+300*10**-6))
Index8 = np.logical_and((Fit[:,0] > Lambda[7]-300*10**-6), (Fit[:,0] < Lambda[7]+300*10**-6))
Index9 = np.logical_and((Fit[:,0] > Lambda[8]-300*10**-6), (Fit[:,0] < Lambda[8]+300*10**-6))
Index10 = np.logical_and((Fit[:,0] > Lambda[9]-300*10**-6), (Fit[:,0] < Lambda[9]+300*10**-6))
Index11 = np.logical_and((Fit[:,0] > Lambda[10]-300*10**-6), (Fit[:,0] < Lambda[10]+300*10**-6))
Index12 = np.logical_and((Fit[:,0] > Lambda[11]-300*10**-6), (Fit[:,0] < Lambda[11]+300*10**-6))
Index13 = np.logical_and((Fit[:,0] > Lambda[12]-300*10**-6), (Fit[:,0] < Lambda[12]+300*10**-6))
Index14 = np.logical_and((Fit[:,0] > Lambda[13]-300*10**-6), (Fit[:,0] < Lambda[13]+300*10**-6))
Index15 = np.logical_and((Fit[:,0] > Lambda[14]-500*10**-6), (Fit[:,0] < Lambda[14]+500*10**-6))
Index16 = np.logical_and((Fit[:,0] > Lambda[15]-500*10**-6), (Fit[:,0] < Lambda[15]+500*10**-6))
Index17 = np.logical_and((Fit[:,0] > Lambda[16]-500*10**-6), (Fit[:,0] < Lambda[16]+500*10**-6))
Index18 = np.logical_and((Fit[:,0] > Lambda[17]-1000*10**-6), (Fit[:,0] < Lambda[17]+1000*10**-6))
Index19 = np.logical_and((Fit[:,0] > Lambda[18]-1000*10**-6), (Fit[:,0] < Lambda[18]+1000*10**-6))
Index20 = np.logical_and((Fit[:,0] > Lambda[19]-1000*10**-6), (Fit[:,0] < Lambda[19]+1000*10**-6))
Index21 = np.logical_and((Fit[:,0] > Lambda[20]-1000*10**-6), (Fit[:,0] < Lambda[20]+1000*10**-6))
Index22 = np.logical_and((Fit[:,0] > Lambda[21]-1000*10**-6), (Fit[:,0] < Lambda[21]+1000*10**-6))
Index23 = np.logical_and((Fit[:,0] > Lambda[22]-1500*10**-6), (Fit[:,0] < Lambda[22]+1500*10**-6))
Index24 = np.logical_and((Fit[:,0] > Lambda[23]-1500*10**-6), (Fit[:,0] < Lambda[23]+1500*10**-6))
Index25 = np.logical_and((Fit[:,0] > Lambda[24]-2000*10**-6), (Fit[:,0] < Lambda[24]+2000*10**-6))
Index26 = np.logical_and((Fit[:,0] > Lambda[25]-1500*10**-6), (Fit[:,0] < Lambda[25]+1500*10**-6))
Index27 = np.logical_and((Fit[:,0] > Lambda[26]-1500*10**-6), (Fit[:,0] < Lambda[26]+1500*10**-6))
Index28 = np.logical_and((Fit[:,0] > Lambda[27]-1000000*10**-6), (Fit[:,0] < Lambda[27]+1000000*10**-6))

idx = np.array([[Index1],[Index2],[Index3],[Index4],[Index5],[Index6],[Index7],[Index8],[Index9],[Index10],
       [Index11],[Index12],[Index13],[Index14],[Index15],[Index16],[Index17],[Index18],[Index19],[Index20],
       [Index21],[Index22],[Index23],[Index24],[Index25],[Index26],[Index27],[Index28]])

mylist  = [[np.mean(Fit[idx[j,0,:],i]) for i in range(1,13)] for j in range(28)]
myarray = np.array(mylist) # (row , colname) = (flux from filter, fit models)

Fit_res = [myarray[:,i] - Gflux_cat1[:,galaxy1] for i in range(12)]
res = np.array(np.absolute(Fit_res)) # (row , colname) = (fit models, residual from filters)

Final_res = [np.mean(res[i,1]) for i in range(12)]

def selection_sort(Final_res):
    for i in range(len(Final_res)):
        swap = i + np.argmin(Final_res[i:])
        (Final_res[i], Final_res[swap]) = (Final_res[swap], Final_res[i])
    return Final_res
Final_res = selection_sort(Final_res)

# Ide Vise de 3 bedtse modeller, men array skal behodle det orginale index, så vi
# ved hvilken fit model der er tale om

best_fit = np.argmin(Final_res)+1
print(best_fit)



# Close FITS file so it won't use up excess memory
hdu.close()
