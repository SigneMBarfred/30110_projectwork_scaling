# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 12:21:59 2023

@author: Nikolaj Lange Dons
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 10:03:50 2023

@author: Nikolaj Lange Dons, Rasmus Bruun, Signe Barfred
"""


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#loading data as recarrays
hdu = fits.open('filtered_SED.fits',memmap=True)
data = hdu[1].data

selection_all = np.logical_and((data['z_phot'] > 0), (data['z_phot'] < 4), (data['CFHT_u_MAG'] < 50))

U = 23.9 - 2.5*np.log10(data['restU'][0,:])
V = 23.9 - 2.5*np.log10(data['restV'][0,:])
J = 23.9 - 2.5*np.log10(data['restJ'][0,:])
mass = data['mass'][0,:]

UV = U[selection_all[0,:]] - V[selection_all[0,:]]
VJ = V[selection_all[0,:]] - J[selection_all[0,:]]

UV.reshape(476363)
VJ.reshape(476363)

UVz = U - V
VJz = V - J

UVz.reshape(511006)
VJz.reshape(511006)

test0 = "%d galaxies"
int0 = len(UV)
print(test0%int0)

######################

# Select QGs and SFGs

# iduv1 = UV >= 1.3

# iduv2 = UV <= 3

# idvj1 = VJ <= 1.6

# idvj2 = VJ >= -1

# iduvj = UV >= 0.875*VJ+0.6

# IDS = iduv1&idvj1&iduvj&idvj2&iduv2

# NQ = sum(IDS)

# NSF = sum(~IDS)

# Grate = NSF/NQ

# print("Number of quiessent galaxies: %d"%NQ)
# print("Number of starforming galaxies: %d"%NSF)
# print("Number of starforming over quiessent galaxies: %.2f"%Grate)

######################

# UVJ diagram for galaxies at all redshifts

fig0, ax = plt.subplots(figsize = (8, 6))
hb0 = ax.hexbin(VJ, UV, vmax = 3000, cmap = "hot", gridsize = (173,100), mincnt = 0)
cb0 = fig0.colorbar(hb0, ax = ax)
cb0.set_label("Number of galaxies")
ax.set_xlim(0,2.5)
ax.set_ylim(0,2.5)
ax.set_xlabel("V - J [AB mag]", fontsize=12)
ax.set_ylabel("U - V [AB mag]", fontsize=12)
ax.set_title("UVJ diagram for galaxies at all redshifts", fontsize=18)

######################

# # UVJ diagram for cat1
# # 0.1 < z < 0.3

# selection_cat1 = np.logical_and((data['z_phot'] > 0.1), (data['z_phot'] < 0.3), (data['CFHT_u_MAG'] < 50))

# UV_cat1 = U[selection_cat1[0,:]] - V[selection_cat1[0,:]]
# VJ_cat1 = V[selection_cat1[0,:]] - J[selection_cat1[0,:]]

# fig1, ax = plt.subplots(figsize = (8, 6))
# hb1 = ax.hexbin(VJ_cat1, UV_cat1, vmax = 50, cmap = "hot", gridsize = (173,100), mincnt = 0)
# cb1 = fig1.colorbar(hb1, ax = ax)
# cb1.set_label("Number of galaxies")
# ax.plot([-4, 0.693], [1.3, 1.3], color = "g", lw = 5) #draws the horizontal line
# ax.plot([1.6, 1.6], [2.098, 4], color = "g", lw = 5) #draws the vertical line
# ax.plot([0.693, 1.6], [1.3, 2.098], color = "g", lw = 5) #draws the tilted line
# ax.set_xlim(0,2.5)
# ax.set_ylim(0,2.5)
# ax.set_xlabel("V - J [AB mag]", fontsize=12)
# ax.set_ylabel("U - V [AB mag]", fontsize=12)
# ax.set_title("UVJ diagram for 0.1 < z < 0.3", fontsize=18)

# ######################

# # UVJ diagram for cat2
# # 0.9 < z < 1.1

# selection_cat2 = np.logical_and((data['z_phot'] > 0.9), (data['z_phot'] < 1.1), (data['CFHT_u_MAG'] < 50))

# UV_cat2 = U[selection_cat2[0,:]] - V[selection_cat2[0,:]]
# VJ_cat2 = V[selection_cat2[0,:]] - J[selection_cat2[0,:]]

# fig2, ax = plt.subplots(figsize = (8, 6))
# hb2 = ax.hexbin(VJ_cat2, UV_cat2, vmax = 50, cmap = "hot", gridsize = (173,100), mincnt = 0)
# cb2 = fig2.colorbar(hb2, ax = ax)
# cb2.set_label("Number of galaxies")
# ax.plot([-4.0, 0.81], [1.3, 1.3], color = "g", lw = 5) #draws the horizontal line
# ax.plot([1.6, 1.6], [2.0, 4], color = "g", lw = 5) #draws the vertical line
# ax.plot([0.81, 1.6], [1.3, 2.0], color = "g", lw = 5) #draws the tilted line
# ax.set_xlim(0,2.5)
# ax.set_ylim(0,2.5)
# ax.set_xlabel("V - J [AB mag]", fontsize=12)
# ax.set_ylabel("U - V [AB mag]", fontsize=12)
# ax.set_title("UVJ diagram for 0.9 < z < 1.1", fontsize=18)

# ######################

# # UVJ diagram for cat3
# # 1.9 < z < 2.1

# selection_cat3 = np.logical_and((data['z_phot'] > 1.9), (data['z_phot'] < 2.1), (data['CFHT_u_MAG'] < 50))

# UV_cat3 = U[selection_cat3[0,:]] - V[selection_cat3[0,:]]
# VJ_cat3 = V[selection_cat3[0,:]] - J[selection_cat3[0,:]]

# fig3, ax = plt.subplots(figsize = (8, 6))
# hb3 = ax.hexbin(VJ_cat3, UV_cat3, vmax = 50, cmap = "hot", gridsize = (173,100), mincnt = 0)
# cb3 = fig3.colorbar(hb3, ax = ax)
# cb3.set_label("Number of galaxies")
# ax.plot([-4.0, 0.92], [1.3, 1.3], color = "g", lw = 5) #draws the horizontal line
# ax.plot([1.6, 1.6], [1.9, 4], color = "g", lw = 5) #draws the vertical line
# ax.plot([0.92, 1.6], [1.3, 1.9], color = "g", lw = 5) #draws the tilted line
# ax.set_xlim(0,2.5)
# ax.set_ylim(0,2.5)
# ax.set_xlabel("V - J [AB mag]", fontsize=12)
# ax.set_ylabel("U - V [AB mag]", fontsize=12)
# ax.set_title("UVJ diagram for 1.9 < z < 2.1", fontsize=18)

##################
# Q/SF
SFR = data['sfr'][0,:]


# cat1: 0 < z < 0.5
idUVJ_cat1 = np.logical_and((data['z_phot'] > 0), (data['z_phot'] < 0.5), (data['CFHT_u_MAG'] < 50))
iduv1_cat1 = UVz[idUVJ_cat1[0,:]] >= 1.3
iduv2_cat1 = UVz[idUVJ_cat1[0,:]] <= 3
idvj1_cat1 = VJz[idUVJ_cat1[0,:]] <= 1.6
idvj2_cat1 = VJz[idUVJ_cat1[0,:]] >= -1
iduvj_cat1 = UVz[idUVJ_cat1[0,:]] >= 0.88*VJz[idUVJ_cat1[0,:]]+0.69
IDS_cat1 = iduv1_cat1 & idvj1_cat1 & iduvj_cat1 & idvj2_cat1 & iduv2_cat1
#For the SFR plots we need to sort out outliers so we create additional ID
id_valid_sfr_cat1 = np.logical_and(np.log10(SFR[idUVJ_cat1[0,:]])<4,np.log10(SFR[idUVJ_cat1[0,:]])>-4)
id_sf_and_valid_sfr_cat1 = ~IDS_cat1 & id_valid_sfr_cat1
nQ_cat1 = sum(IDS_cat1)
NSF_cat1 = sum(~IDS_cat1)
cat1_total = nQ_cat1 + NSF_cat1
# Mass
m1 = mass[idUVJ_cat1[0,:]]
NQ_m1 = np.average(m1[IDS_cat1])
nSF_m1 = np.average(m1[~IDS_cat1])

## sub category with smaller z interval , 0<z<0.1
idUVJ_cat1_1 = np.logical_and((data['z_phot'] > 0), (data['z_phot'] < 0.1), (data['CFHT_u_MAG'] < 50))
iduv1_cat1_1 = UVz[idUVJ_cat1_1[0,:]] >= 1.3
iduv2_cat1_1 = UVz[idUVJ_cat1_1[0,:]] <= 3
idvj1_cat1_1 = VJz[idUVJ_cat1_1[0,:]] <= 1.6
idvj2_cat1_1 = VJz[idUVJ_cat1_1[0,:]] >= -1
iduvj_cat1_1 = UVz[idUVJ_cat1_1[0,:]] >= 0.88*VJz[idUVJ_cat1_1[0,:]]+0.69
IDS_cat1_1 = iduv1_cat1_1 & idvj1_cat1_1 & iduvj_cat1_1 & idvj2_cat1_1 & iduv2_cat1_1
nQ_cat1_1 = sum(IDS_cat1_1)
NSF_cat1_1 = sum(~IDS_cat1_1)
cat1_1_total = nQ_cat1_1 + NSF_cat1_1

# cat2: 0.25 < z < 0.75

#vælger datainterval for rødforskydning og frasorterer ugyldige magnituder
idUVJ_cat2_lower = np.logical_and((data['z_phot'] > 0.25), (data['z_phot'] < 0.5), (data['CFHT_u_MAG'] < 50)) 
idUVJ_cat2_upper = np.logical_and((data['z_phot'] > 0.5), (data['z_phot'] < 0.75), (data['CFHT_u_MAG'] < 50))

iduv1_cat2_lower = UVz[idUVJ_cat2_lower[0,:]] >= 1.3
iduv2_cat2_lower = UVz[idUVJ_cat2_lower[0,:]] <= 3
idvj1_cat2_lower = VJz[idUVJ_cat2_lower[0,:]] <= 1.6
idvj2_cat2_lower = VJz[idUVJ_cat2_lower[0,:]] >= -1
iduvj_cat2_lower = UVz[idUVJ_cat2_lower[0,:]] >= 0.88*VJz[idUVJ_cat2_lower[0,:]]+0.69
IDS_cat2_lower = iduv1_cat2_lower & idvj1_cat2_lower & iduvj_cat2_lower & idvj2_cat2_lower & iduv2_cat2_lower

iduv1_cat2_upper = UVz[idUVJ_cat2_upper[0,:]] >= 1.3
iduv2_cat2_upper = UVz[idUVJ_cat2_upper[0,:]] <= 3
idvj1_cat2_upper = VJz[idUVJ_cat2_upper[0,:]] <= 1.6
idvj2_cat2_upper = VJz[idUVJ_cat2_upper[0,:]] >= -1
iduvj_cat2_upper = UVz[idUVJ_cat2_upper[0,:]] >= 0.88*VJz[idUVJ_cat2_upper[0,:]]+0.59
IDS_cat2_upper = iduv1_cat2_upper & idvj1_cat2_upper & iduvj_cat2_upper & idvj2_cat2_upper & iduv2_cat2_upper

nQ_cat2 = sum(IDS_cat2_lower)+sum(IDS_cat2_upper)
NSF_cat2 = sum(~IDS_cat2_lower)+sum(~IDS_cat2_upper)
cat2_total = nQ_cat2 + NSF_cat2

# cat3: 0.5 < z < 1.0
idUVJ_cat3 = np.logical_and((data['z_phot'] > 0.5), (data['z_phot'] < 1), (data['CFHT_u_MAG'] < 50))
iduv1_cat3 = UVz[idUVJ_cat3[0,:]] >= 1.3
iduv2_cat3 = UVz[idUVJ_cat3[0,:]] <= 3
idvj1_cat3 = VJz[idUVJ_cat3[0,:]] <= 1.6
idvj2_cat3 = VJz[idUVJ_cat3[0,:]] >= -1
iduvj_cat3 = UVz[idUVJ_cat3[0,:]] >= 0.88*VJz[idUVJ_cat3[0,:]]+0.59
IDS_cat3 = iduv1_cat3 & idvj1_cat3 & iduvj_cat3 & idvj2_cat3 & iduv2_cat3
nQ_cat3 = sum(IDS_cat3)
NSF_cat3 = sum(~IDS_cat3)
cat3_total = nQ_cat3 + NSF_cat3
# Mass
m3 = mass[idUVJ_cat3[0,:]]
NQ_m3 = np.average(m3[IDS_cat3])
nSF_m3 = np.average(m3[~IDS_cat3])

# cat4: 0.75 < z < 1.25
idUVJ_cat4_lower = np.logical_and((data['z_phot'] > 0.75), (data['z_phot'] < 1), (data['CFHT_u_MAG'] < 50))
idUVJ_cat4_upper = np.logical_and((data['z_phot'] > 0.5), (data['z_phot'] < 1.25), (data['CFHT_u_MAG'] < 50))

iduv1_cat4_lower = UVz[idUVJ_cat4_lower[0,:]] >= 1.3
iduv2_cat4_lower = UVz[idUVJ_cat4_lower[0,:]] <= 3
idvj1_cat4_lower = VJz[idUVJ_cat4_lower[0,:]] <= 1.6
idvj2_cat4_lower = VJz[idUVJ_cat4_lower[0,:]] >= -1
iduvj_cat4_lower = UVz[idUVJ_cat4_lower[0,:]] >= 0.88*VJz[idUVJ_cat4_lower[0,:]]+0.59
IDS_cat4_lower = iduv1_cat4_lower & idvj1_cat4_lower & iduvj_cat4_lower & idvj2_cat4_lower & iduv2_cat4_lower

iduv1_cat4_upper = UVz[idUVJ_cat4_upper[0,:]] >= 1.3
iduv2_cat4_upper = UVz[idUVJ_cat4_upper[0,:]] <= 3
idvj1_cat4_upper = VJz[idUVJ_cat4_upper[0,:]] <= 1.6
idvj2_cat4_upper = VJz[idUVJ_cat4_upper[0,:]] >= -1
iduvj_cat4_upper = UVz[idUVJ_cat4_upper[0,:]] >= 0.88*VJz[idUVJ_cat4_upper[0,:]]+0.49
IDS_cat4_upper = iduv1_cat4_upper & idvj1_cat4_upper & iduvj_cat4_upper & idvj2_cat4_upper & iduv2_cat4_upper

nQ_cat4 = sum(IDS_cat4_lower)+sum(IDS_cat4_upper)
NSF_cat4 = sum(~IDS_cat4_lower)+sum(~IDS_cat4_upper)
cat4_total = nQ_cat4 + NSF_cat4

# cat5: 1.0 < z < 1.5
idUVJ_cat5 = np.logical_and((data['z_phot'] > 1), (data['z_phot'] < 1.5), (data['CFHT_u_MAG'] < 50))
iduv1_cat5 = UVz[idUVJ_cat5[0,:]] >= 1.3
iduv2_cat5 = UVz[idUVJ_cat5[0,:]] <= 3
idvj1_cat5 = VJz[idUVJ_cat5[0,:]] <= 1.6
idvj2_cat5 = VJz[idUVJ_cat5[0,:]] >= -1
iduvj_cat5 = UVz[idUVJ_cat5[0,:]] >= 0.88*VJz[idUVJ_cat5[0,:]]+0.49
IDS_cat5 = iduv1_cat5 & idvj1_cat5 & iduvj_cat5 & idvj2_cat5 & iduv2_cat5
NQ_cat5 = sum(IDS_cat5)
NSF_cat5 = sum(~IDS_cat5)
cat5_total = NQ_cat5 + NSF_cat5
# Mass
m5 = mass[idUVJ_cat5[0,:]]
NQ_m5 = np.average(m5[IDS_cat5])
nSF_m5 = np.average(m5[~IDS_cat5])

# cat6: 1.25 < z < 1.75
idUVJ_cat6 = np.logical_and((data['z_phot'] > 1.25), (data['z_phot'] < 1.75), (data['CFHT_u_MAG'] < 50))
iduv1_cat6 = UVz[idUVJ_cat6[0,:]] >= 1.3
iduv2_cat6 = UVz[idUVJ_cat6[0,:]] <= 3
idvj1_cat6 = VJz[idUVJ_cat6[0,:]] <= 1.6
idvj2_cat6 = VJz[idUVJ_cat6[0,:]] >= -1
iduvj_cat6 = UVz[idUVJ_cat6[0,:]] >= 0.88*VJz[idUVJ_cat6[0,:]]+0.49
IDS_cat6 = iduv1_cat6 & idvj1_cat6 & iduvj_cat6 & idvj2_cat6 & iduv2_cat6
NQ_cat6 = sum(IDS_cat6)
NSF_cat6 = sum(~IDS_cat6)
cat6_total = NQ_cat6 + NSF_cat6

# cat7: 1.5 < z < 2
idUVJ_cat7 = np.logical_and((data['z_phot'] > 1.5), (data['z_phot'] < 2), (data['CFHT_u_MAG'] < 50))
iduv1_cat7 = UVz[idUVJ_cat7[0,:]] >= 1.3
iduv2_cat7 = UVz[idUVJ_cat7[0,:]] <= 3
idvj1_cat7 = VJz[idUVJ_cat7[0,:]] <= 1.6
idvj2_cat7 = VJz[idUVJ_cat7[0,:]] >= -1
iduvj_cat7 = UVz[idUVJ_cat7[0,:]] >= 0.88*VJz[idUVJ_cat7[0,:]]+0.49
IDS_cat7 = iduv1_cat7 & idvj1_cat7 & iduvj_cat7 & idvj2_cat7 & iduv2_cat7
NQ_cat7 = sum(IDS_cat7)
NSF_cat7 = sum(~IDS_cat7)
cat7_total = NQ_cat7 + NSF_cat7
# Mass
m7 = mass[idUVJ_cat7[0,:]]
NQ_m7 = np.average(m7[IDS_cat7])
nSF_m7 = np.average(m7[~IDS_cat7])

# cat8: 1.75 < z < 2.25
idUVJ_cat8_lower = np.logical_and((data['z_phot'] > 1.75), (data['z_phot'] < 2), (data['CFHT_u_MAG'] < 50))
idUVJ_cat8_upper = np.logical_and((data['z_phot'] > 2), (data['z_phot'] < 2.25), (data['CFHT_u_MAG'] < 50))

iduv1_cat8_lower = UVz[idUVJ_cat8_lower[0,:]] >= 1.3
iduv2_cat8_lower = UVz[idUVJ_cat8_lower[0,:]] <= 3
idvj1_cat8_lower = VJz[idUVJ_cat8_lower[0,:]] <= 1.6
idvj2_cat8_lower = VJz[idUVJ_cat8_lower[0,:]] >= -1
iduvj_cat8_lower = UVz[idUVJ_cat8_lower[0,:]] >= 0.88*VJz[idUVJ_cat8_lower[0,:]]+0.49
IDS_cat8_lower = iduv1_cat8_lower & idvj1_cat8_lower & iduvj_cat8_lower & idvj2_cat8_lower & iduv2_cat8_lower

iduv1_cat8_upper = UVz[idUVJ_cat8_upper[0,:]] >= 1.3
iduv2_cat8_upper = UVz[idUVJ_cat8_upper[0,:]] <= 3
idvj1_cat8_upper = VJz[idUVJ_cat8_upper[0,:]] <= 1.6
idvj2_cat8_upper = VJz[idUVJ_cat8_upper[0,:]] >= -1
iduvj_cat8_upper = UVz[idUVJ_cat8_upper[0,:]] >= 0.88*VJz[idUVJ_cat8_upper[0,:]]+0.39
IDS_cat8_upper = iduv1_cat8_upper & idvj1_cat8_upper & iduvj_cat8_upper & idvj2_cat8_upper & iduv2_cat8_upper

NQ_cat8 = sum(IDS_cat8_lower)+sum(IDS_cat8_upper)
NSF_cat8 = sum(~IDS_cat8_lower)+sum(~IDS_cat8_upper)
cat8_total = NQ_cat8 + NSF_cat8

# cat9: 2 < z < 2.5
idUVJ_cat9 = np.logical_and((data['z_phot'] > 2), (data['z_phot'] < 2.5), (data['CFHT_u_MAG'] < 50))
iduv1_cat9 = UVz[idUVJ_cat9[0,:]] >= 1.3
iduv2_cat9 = UVz[idUVJ_cat9[0,:]] <= 3
idvj1_cat9 = VJz[idUVJ_cat9[0,:]] <= 1.6
idvj2_cat9 = VJz[idUVJ_cat9[0,:]] >= -1
iduvj_cat9 = UVz[idUVJ_cat9[0,:]] >= 0.88*VJz[idUVJ_cat9[0,:]]+0.39
IDS_cat9 = iduv1_cat9 & idvj1_cat9 & iduvj_cat9 & idvj2_cat9 & iduv2_cat9
NQ_cat9 = sum(IDS_cat9)
NSF_cat9 = sum(~IDS_cat9)
cat9_total = NQ_cat9 + NSF_cat9
# Mass
m9 = mass[idUVJ_cat9[0,:]]
NQ_m9 = np.average(m9[IDS_cat9])
nSF_m9 = np.average(m9[~IDS_cat9])


# cat10: 2.25 < z < 2.75
idUVJ_cat10 = np.logical_and((data['z_phot'] > 2.25), (data['z_phot'] < 2.75), (data['CFHT_u_MAG'] < 50))
iduv1_cat10 = UVz[idUVJ_cat10[0,:]] >= 1.3
iduv2_cat10 = UVz[idUVJ_cat10[0,:]] <= 3
idvj1_cat10 = VJz[idUVJ_cat10[0,:]] <= 1.6
idvj2_cat10 = VJz[idUVJ_cat10[0,:]] >= -1
iduvj_cat10 = UVz[idUVJ_cat10[0,:]] >= 0.88*VJz[idUVJ_cat10[0,:]]+0.39
IDS_cat10 = iduv1_cat10 & idvj1_cat10 & iduvj_cat10 & idvj2_cat10 & iduv2_cat10
nQ_cat10 = sum(IDS_cat10)
NSF_cat10 = sum(~IDS_cat10)
cat10_total = nQ_cat10 + NSF_cat10

# cat11: 2.5 < z < 3
idUVJ_cat11 = np.logical_and((data['z_phot'] > 2.5), (data['z_phot'] < 3), (data['CFHT_u_MAG'] < 50))
iduv1_cat11 = UVz[idUVJ_cat11[0,:]] >= 1.3
iduv2_cat11 = UVz[idUVJ_cat11[0,:]] <= 3
idvj1_cat11 = VJz[idUVJ_cat11[0,:]] <= 1.6
idvj2_cat11 = VJz[idUVJ_cat11[0,:]] >= -1
iduvj_cat11 = UVz[idUVJ_cat11[0,:]] >= 0.88 * VJz[idUVJ_cat11[0,:]] + 0.39
IDS_cat11 = iduv1_cat11 & idvj1_cat11 & iduvj_cat11 & idvj2_cat11 & iduv2_cat11
nQ_cat11 = sum(IDS_cat11)
NSF_cat11 = sum(~IDS_cat11)
cat11_total = nQ_cat11 + NSF_cat11
# Mass
m11 = mass[idUVJ_cat11[0,:]]
NQ_m11 = np.average(m11[IDS_cat11])
nSF_m11 = np.average(m11[~IDS_cat11])

# cat12: 2.75 < z < 3.25
idUVJ_cat12 = np.logical_and((data['z_phot'] > 2.75), (data['z_phot'] < 3.25), (data['CFHT_u_MAG'] < 50))
iduv1_cat12 = UVz[idUVJ_cat12[0,:]] >= 1.3
iduv2_cat12 = UVz[idUVJ_cat12[0,:]] <= 3
idvj1_cat12 = VJz[idUVJ_cat12[0,:]] <= 1.6
idvj2_cat12 = VJz[idUVJ_cat12[0,:]] >= -1
iduvj_cat12 = UVz[idUVJ_cat12[0,:]] >= 0.88 * VJz[idUVJ_cat12[0,:]] + 0.39
IDS_cat12 = iduv1_cat12 & idvj1_cat12 & iduvj_cat12 & idvj2_cat12 & iduv2_cat12
nQ_cat12 = sum(IDS_cat12)
NSF_cat12 = sum(~IDS_cat12)
cat12_total = nQ_cat12 + NSF_cat12

# cat13: 3 < z < 3.5
idUVJ_cat13 = np.logical_and((data['z_phot'] > 3), (data['z_phot'] < 3.5), (data['CFHT_u_MAG'] < 50))
iduv1_cat13 = UVz[idUVJ_cat13[0,:]] >= 1.3
iduv2_cat13 = UVz[idUVJ_cat13[0,:]] <= 3
idvj1_cat13 = VJz[idUVJ_cat13[0,:]] <= 1.6
idvj2_cat13 = VJz[idUVJ_cat13[0,:]] >= -1
iduvj_cat13 = UVz[idUVJ_cat13[0,:]] >= 0.88 * VJz[idUVJ_cat13[0,:]] + 0.39
IDS_cat13 = iduv1_cat13 & idvj1_cat13 & iduvj_cat13 & idvj2_cat13 & iduv2_cat13
nQ_cat13 = sum(IDS_cat13)
NSF_cat13 = sum(~IDS_cat13)
cat13_total = nQ_cat13 + NSF_cat13
# Mass
m13 = mass[idUVJ_cat13[0,:]]

# cat14: 3.25 < z < 3.75
idUVJ_cat14 = np.logical_and((data['z_phot'] > 3.25), (data['z_phot'] < 3.75), (data['CFHT_u_MAG'] < 50))
iduv1_cat14 = UVz[idUVJ_cat14[0,:]] >= 1.3
iduv2_cat14 = UVz[idUVJ_cat14[0,:]] <= 3
idvj1_cat14 = VJz[idUVJ_cat14[0,:]] <= 1.6
idvj2_cat14 = VJz[idUVJ_cat14[0,:]] >= -1
iduvj_cat14 = UVz[idUVJ_cat14[0,:]] >= 0.88 * VJz[idUVJ_cat14[0,:]] + 0.39
IDS_cat14 = iduv1_cat14 & idvj1_cat14 & iduvj_cat14 & idvj2_cat14 & iduv2_cat14
nQ_cat14 = sum(IDS_cat14)
NSF_cat14 = sum(~IDS_cat14)
cat14_total = nQ_cat14 + NSF_cat14

# cat15: 3.5 < z < 4
idUVJ_cat15 = np.logical_and((data['z_phot'] > 3.5), (data['z_phot'] < 4), (data['CFHT_u_MAG'] < 50))
iduv1_cat15 = UVz[idUVJ_cat15[0,:]] >= 1.3
iduv2_cat15 = UVz[idUVJ_cat15[0,:]] <= 3
idvj1_cat15 = VJz[idUVJ_cat15[0,:]] <= 1.6
idvj2_cat15 = VJz[idUVJ_cat15[0,:]] >= -1
iduvj_cat15 = UVz[idUVJ_cat15[0,:]] >= 0.88 * VJz[idUVJ_cat15[0,:]] + 0.39
IDS_cat15 = iduv1_cat15 & idvj1_cat15 & iduvj_cat15 & idvj2_cat15 & iduv2_cat15
nQ_cat15 = sum(IDS_cat15)
NSF_cat15 = sum(~IDS_cat15)
cat15_total = nQ_cat15 + NSF_cat15
#mass
m15 = mass[idUVJ_cat15[0,:]]


# Create a matrix containing intervals and no. of galaxies in each interval
QSF = np.array([[0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75], 
         [nQ_cat1,nQ_cat2,nQ_cat3,nQ_cat4,NQ_cat5,NQ_cat6,NQ_cat7,NQ_cat8,NQ_cat9,
          nQ_cat10,nQ_cat11,nQ_cat12,nQ_cat13,nQ_cat14,nQ_cat15], 
         [NSF_cat1,NSF_cat2,NSF_cat3,NSF_cat4,NSF_cat5,NSF_cat6,NSF_cat7,NSF_cat8,
          NSF_cat9,NSF_cat10,NSF_cat11,NSF_cat12,NSF_cat13,NSF_cat14,NSF_cat15],
         [cat1_total,cat2_total,cat3_total,cat4_total,cat5_total,cat6_total,
          cat7_total,cat8_total,cat9_total,cat10_total,cat11_total,cat12_total,
          cat13_total,cat14_total,cat15_total]])
# QSF = [z, num quiesent, num starforming, total]

#initialize layout of plot for nQ/nSF
fig1, ax = plt.subplots()

#add scatterplot for cat1
SF = ax.plot(QSF[0,:],QSF[2,:]/QSF[3,:], 'b', marker='o', label="nSF")
Q = ax.plot(QSF[0,:],QSF[1,:]/QSF[3,:],'r', marker='o', label="nQ")
plt.legend()

ax.set_xlim(0,4)
ax.set_xlabel("Redshift z")
ax.set_ylabel("number of galaxies")
ax.set_title("Number of galaxies at redshift z")
########################
# Quiessent boundary

def intersect(a, b, c, d):
    t = ((a[0] - c[0]) * (c[1] - d[1]) - (a[1] - c[1]) * (c[0] - d[0])) / ((a[0] - b[0]) * (c[1] - d[1]) - (a[1] - b[1]) * (c[0] - d[0]))
    u = ((a[0] - c[0]) * (a[1] - b[1]) - (a[1] - c[1]) * (a[0] - b[0])) / ((a[0] - b[0]) * (c[1] - d[1]) - (a[1] - b[1]) * (c[0] - d[0]))

    # check if line actually intersect
    if (0 <= t and t <= 1 and 0 <= u and u <= 1):
        return [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
    else: 
        return False
# Function for intersection between two line


# 0 < z < 0.5

# Horizontal line 1
A = (0, 1.3)
B = (2.5, 1.3)
# cross line 1
C1 = (np.nanmin(VJ),0.88*(np.nanmin(VJ))+0.69)
D1 = (np.nanmax(VJ),0.88*(np.nanmax(VJ))+0.69)
# vertical line 1
E = (1.6, 0)
F = (1.6, 2.5)

P11 = intersect(A, B, C1, D1)
P12 = intersect(E, F, C1, D1)

# 0.5 < z < 1
# cross line 2
C2 = (np.nanmin(VJ),0.88*(np.nanmin(VJ))+0.59)
D2 = (np.nanmax(VJ),0.88*(np.nanmax(VJ))+0.59)

P21 = intersect(A, B, C2, D2)
P22 = intersect(E, F, C2, D2)

# 1 < z < 2
# cross line 3
C3 = (np.nanmin(VJ),0.88*(np.nanmin(VJ))+0.49)
D3 = (np.nanmax(VJ),0.88*(np.nanmax(VJ))+0.49)

P31 = intersect(A, B, C3, D3)
P32 = intersect(E, F, C3, D3)

# Filled polygon
Qpoly_x = [A[0],P31[0],P32[0],P32[0],A[0],A[0]]
Qpoly_y = [A[1],P31[1],P32[1],2.5,2.5,A[1]]

SFpoly_x = [A[0],2.5,2.5,P32[0],P32[0],P31[0],A[0]]
SFpoly_y = [A[0],0,2.5,2.5,P32[1],P31[1],A[1]]

########################
#UVJ diagram at each redshift

# UVJ diagram for cat1,2,3 og 4
idcat1 = np.logical_and((data['z_phot'] > 0), (data['z_phot'] < 0.5))
idcat2 = np.logical_and((data['z_phot'] > 0.25), (data['z_phot'] < 0.75))
idcat3 = np.logical_and((data['z_phot'] > 0.5), (data['z_phot'] < 1))
idcat4 = np.logical_and((data['z_phot'] > 0.75), (data['z_phot'] < 1.25))

UVcat1 = U[idcat1[0,:]] - V[idcat1[0,:]]
VJcat1 = V[idcat1[0,:]] - J[idcat1[0,:]]
UVcat2 = U[idcat2[0,:]] - V[idcat2[0,:]]
VJcat2 = V[idcat2[0,:]] - J[idcat2[0,:]]
UVcat3 = U[idcat3[0,:]] - V[idcat3[0,:]]
VJcat3 = V[idcat3[0,:]] - J[idcat3[0,:]]
UVcat4 = U[idcat4[0,:]] - V[idcat4[0,:]]
VJcat4 = V[idcat4[0,:]] - J[idcat4[0,:]]

figure, ax = plt.subplots(2, 2, sharex=(True), sharey=(True), figsize=(12,10) )
plt.tight_layout()

# cat1
hb1 = ax[0,0].hexbin(VJcat1, UVcat1, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
cb1 = figure.colorbar(hb1, ax = ax)
cb1.set_label("Number of galaxies", fontsize=15)
ax[0,0].plot((A[0],P11[0]), (A[1],P11[1]), 'black', linestyle="--")
ax[0,0].plot((P11[0],P12[0]), (P11[1],P12[1]), 'black', linestyle="--")
ax[0,0].plot((P12[0],F[0]), (P12[1],F[1]), 'black', linestyle="--")
ax[0,0].set_xlim(0,2.5)
ax[0,0].set_ylim(0,2.5)
ax[0,0].set_ylabel("$(U-V)_{rest}$", fontsize=12)
ax[0,0].text(0.25,0.925,"0 < z < 0.5",ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)
ax[0,0].text(0.25,0.55,"QGs: %d" %(nQ_cat1),color='red',ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)
ax[0,0].text(0.3,0.45,"SFGs: %d" %(NSF_cat1),color='blue',ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)

# cat2
hb2 = ax[0,1].hexbin(VJcat2, UVcat2, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[0,1].plot((A[0],P21[0]), (A[1],P21[1]), 'black', linestyle="--")
ax[0,1].plot((P21[0],P22[0]), (P21[1],P22[1]), 'black', linestyle="--")
ax[0,1].plot((P22[0],F[0]), (P22[1],F[1]), 'black', linestyle="--")
ax[0,1].set_xlim(0,2.5)
ax[0,1].set_ylim(0,2.5)
ax[0,1].text(0.35,0.925,"0.25 < z <0.75",ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)
ax[0,1].text(0.3,0.65,"QGs: %d" %(nQ_cat2),color='red',ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)
ax[0,1].text(0.27,0.45,"SFGs: %d" %(NSF_cat2),color='blue',ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)

# cat3
hb3 = ax[1,0].hexbin(VJcat3, UVcat3, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,0].plot((A[0],P21[0]), (A[1],P21[1]), 'black', linestyle="--")
ax[1,0].plot((P21[0],P22[0]), (P21[1],P22[1]), 'black', linestyle="--")
ax[1,0].plot((P22[0],F[0]), (P22[1],F[1]), 'black', linestyle="--")
ax[1,0].set_xlim(0,2.5)
ax[1,0].set_ylim(0,2.5)
ax[1,0].set_xlabel("$(V-J)_{rest}$", fontsize=12)
ax[1,0].set_ylabel("$(U-V)_{rest}$", fontsize=12)
ax[1,0].text(0.25,0.925,"0.5 < z < 1",ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)
ax[1,0].text(0.25,0.65,"QGs: %d" %(nQ_cat3),color='red',ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)
ax[1,0].text(0.25,0.45,"SFGs: %d" %(NSF_cat3),color='blue',ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)

# cat4
hb4 = ax[1,1].hexbin(VJcat4, UVcat4, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,1].plot((A[0],P21[0]), (A[1],P21[1]), 'black', linestyle="--")
ax[1,1].plot((P21[0],P22[0]), (P21[1],P22[1]), 'black', linestyle="--")
ax[1,1].plot((P22[0],F[0]), (P22[1],F[1]), 'black', linestyle="--")
ax[1,1].set_xlim(0,2.5)
ax[1,1].set_ylim(0,2.5)
ax[1,1].set_xlabel("$(V-J)_{rest}$", fontsize=12)
ax[1,1].text(0.35,0.925,"0.75 < z < 1.25",ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)
ax[1,1].text(0.3,0.55,"QGs: %dd" %(nQ_cat4),color='red',ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)
ax[1,1].text(0.3,0.45,"SFGs: %d" %(NSF_cat4),color='blue',ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)

plt.show()

# UVJ diagram for cat5,6,7 og 8
idcat5 = np.logical_and((data['z_phot'] > 1), (data['z_phot'] < 1.5))
idcat6 = np.logical_and((data['z_phot'] > 1.25), (data['z_phot'] < 1.75))
idcat7 = np.logical_and((data['z_phot'] > 1.5), (data['z_phot'] < 2))
idcat8 = np.logical_and((data['z_phot'] > 1.75), (data['z_phot'] < 2.25))

UVcat5 = U[idcat5[0,:]] - V[idcat5[0,:]]
VJcat5 = V[idcat5[0,:]] - J[idcat5[0,:]]
UVcat6 = U[idcat6[0,:]] - V[idcat6[0,:]]
VJcat6 = V[idcat6[0,:]] - J[idcat6[0,:]]
UVcat7 = U[idcat7[0,:]] - V[idcat7[0,:]]
VJcat7 = V[idcat7[0,:]] - J[idcat7[0,:]]
UVcat8 = U[idcat8[0,:]] - V[idcat8[0,:]]
VJcat8 = V[idcat8[0,:]] - J[idcat8[0,:]]

figure, ax = plt.subplots(2, 2, sharex=(True), sharey=(True), figsize=(12,10) )
plt.tight_layout()

# cat5
hb5 = ax[0,0].hexbin(VJcat5, UVcat5, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
cb5 = figure.colorbar(hb5, ax = ax)
cb5.set_label("Number of galaxies", fontsize=15)
ax[0,0].plot((A[0],P31[0]), (A[1],P31[1]), 'black', linestyle="--")
ax[0,0].plot((P31[0],P32[0]), (P31[1],P32[1]), 'black', linestyle="--")
ax[0,0].plot((P32[0],F[0]), (P32[1],F[1]), 'black', linestyle="--")
ax[0,0].set_xlim(0,2.5)
ax[0,0].set_ylim(0,2.5)
ax[0,0].set_ylabel("$(U-V)_{rest}$", fontsize=12)
ax[0,0].text(0.25,0.925,"1 < z < 1.5",ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)
ax[0,0].text(0.25,0.55,"QGs: %d"%(NQ_cat5),color='red',ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)
ax[0,0].text(0.3,0.45,"SFGs: %d"%(NSF_cat5),color='blue',ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)

# cat6
hb6 = ax[0,1].hexbin(VJcat6, UVcat6, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[0,1].plot((A[0],P31[0]), (A[1],P31[1]), 'black', linestyle="--")
ax[0,1].plot((P31[0],P32[0]), (P31[1],P32[1]), 'black', linestyle="--")
ax[0,1].plot((P32[0],F[0]), (P32[1],F[1]), 'black', linestyle="--")
ax[0,1].set_xlim(0,2.5)
ax[0,1].set_ylim(0,2.5)
ax[0,1].text(0.35,0.925,"1.25 < z <1.75",ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)
ax[0,1].text(0.25,0.55,"QGs: %d"%(NQ_cat6),color='red',ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)
ax[0,1].text(0.3,0.45,"SFGs: %d"%(NSF_cat6),color='blue',ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)

# cat7
hb7 = ax[1,0].hexbin(VJcat7, UVcat7, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,0].plot((A[0],P31[0]), (A[1],P31[1]), 'black', linestyle="--")
ax[1,0].plot((P31[0],P32[0]), (P31[1],P32[1]), 'black', linestyle="--")
ax[1,0].plot((P32[0],F[0]), (P32[1],F[1]), 'black', linestyle="--")
ax[1,0].set_xlim(0,2.5)
ax[1,0].set_ylim(0,2.5)
ax[1,0].set_xlabel("$(V-J)_{rest}$", fontsize=12)
ax[1,0].set_ylabel("$(U-V)_{rest}$", fontsize=12)
ax[1,0].text(0.25,0.925,"1.5 < z < 2",ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)
ax[1,0].text(0.25,0.55,"QGs: %d"%(NQ_cat7),color='red',ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)
ax[1,0].text(0.3,0.45,"SFGs: %d"%(NSF_cat7),color='blue',ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)

# cat8
hb8 = ax[1,1].hexbin(VJcat8, UVcat8, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,1].plot((A[0],P31[0]), (A[1],P31[1]), 'black', linestyle="--")
ax[1,1].plot((P31[0],P32[0]), (P31[1],P32[1]), 'black', linestyle="--")
ax[1,1].plot((P32[0],F[0]), (P32[1],F[1]), 'black', linestyle="--")
ax[1,1].set_xlim(0,2.5)
ax[1,1].set_ylim(0,2.5)
ax[1,1].set_xlabel("$(V-J)_{rest}$", fontsize=12)
ax[1,1].text(0.35,0.925,"1.75 < z < 2.25",ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)
ax[1,1].text(0.25,0.55,"QGs: %d"%(NQ_cat8),color='red',ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)
ax[1,1].text(0.3,0.45,"SFGs: %d"%(NSF_cat8),color='blue',ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)
plt.show()

# UVJ diagram for cat9,10,11 og 12
idcat9 = np.logical_and((data['z_phot'] > 2), (data['z_phot'] < 2.5))
idcat10 = np.logical_and((data['z_phot'] > 2.25), (data['z_phot'] < 2.75))
idcat11 = np.logical_and((data['z_phot'] > 2.5), (data['z_phot'] < 3))
idcat12 = np.logical_and((data['z_phot'] > 2.75), (data['z_phot'] < 3.25))

UVcat9 = U[idcat9[0,:]] - V[idcat9[0,:]]
VJcat9 = V[idcat9[0,:]] - J[idcat9[0,:]]
UVcat10 = U[idcat10[0,:]] - V[idcat10[0,:]]
VJcat10 = V[idcat10[0,:]] - J[idcat10[0,:]]
UVcat11 = U[idcat11[0,:]] - V[idcat11[0,:]]
VJcat11 = V[idcat11[0,:]] - J[idcat11[0,:]]
UVcat12 = U[idcat12[0,:]] - V[idcat12[0,:]]
VJcat12 = V[idcat12[0,:]] - J[idcat12[0,:]]

figure, ax = plt.subplots(2, 2, sharex=(True), sharey=(True), figsize=(12,10) )
plt.tight_layout()

# cat9
hb9 = ax[0,0].hexbin(VJcat9, UVcat9, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
cb9 = figure.colorbar(hb9, ax = ax)
cb9.set_label("Number of galaxies", fontsize=15)
ax[0,0].plot((A[0],P31[0]), (A[1],P31[1]), 'black', linestyle="--")
ax[0,0].plot((P31[0],P32[0]), (P31[1],P32[1]), 'black', linestyle="--")
ax[0,0].plot((P32[0],F[0]), (P32[1],F[1]), 'black', linestyle="--")
ax[0,0].set_xlim(0,2.5)
ax[0,0].set_ylim(0,2.5)
ax[0,0].set_ylabel("$(U-V)_{rest}$", fontsize=12)
ax[0,0].text(0.25,0.925,"2 < z < 2.5",ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)
ax[0,0].text(0.25,0.55,"QGs: %d"%(NQ_cat9),color='red',ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)
ax[0,0].text(0.3,0.45,"SFGs: %d"%(NSF_cat9),color='blue',ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)

# cat10
hb10 = ax[0,1].hexbin(VJcat10, UVcat10, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[0,1].plot((A[0],P31[0]), (A[1],P31[1]), 'black', linestyle="--")
ax[0,1].plot((P31[0],P32[0]), (P31[1],P32[1]), 'black', linestyle="--")
ax[0,1].plot((P32[0],F[0]), (P32[1],F[1]), 'black', linestyle="--")
ax[0,1].set_xlim(0,2.5)
ax[0,1].set_ylim(0,2.5)
ax[0,1].text(0.35,0.925,"2.25 < z <2.75",ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)
ax[0,1].text(0.25,0.55,"QGs: %d"%(nQ_cat10),color='red',ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)
ax[0,1].text(0.3,0.45,"SFGs: %d"%(NSF_cat10),color='blue',ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)

# cat11
hb11 = ax[1,0].hexbin(VJcat11, UVcat11, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,0].plot((A[0],P31[0]), (A[1],P31[1]), 'black', linestyle="--")
ax[1,0].plot((P31[0],P32[0]), (P31[1],P32[1]), 'black', linestyle="--")
ax[1,0].plot((P32[0],F[0]), (P32[1],F[1]), 'black', linestyle="--")
ax[1,0].set_xlim(0,2.5)
ax[1,0].set_ylim(0,2.5)
ax[1,0].set_xlabel("$(V-J)_{rest}$", fontsize=12)
ax[1,0].set_ylabel("$(U-V)_{rest}$", fontsize=12)
ax[1,0].text(0.25,0.925,"2.5 < z < 3",ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)
ax[1,0].text(0.25,0.55,"QGs: %d"%(nQ_cat11),color='red',ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)
ax[1,0].text(0.3,0.45,"SFGs: %d"%(NSF_cat11),color='blue',ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)

# cat12
hb12 = ax[1,1].hexbin(VJcat12, UVcat12, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,1].plot((A[0],P31[0]), (A[1],P31[1]), 'black', linestyle="--")
ax[1,1].plot((P31[0],P32[0]), (P31[1],P32[1]), 'black', linestyle="--")
ax[1,1].plot((P32[0],F[0]), (P32[1],F[1]), 'black', linestyle="--")
ax[1,1].set_xlim(0,2.5)
ax[1,1].set_ylim(0,2.5)
ax[1,1].set_xlabel("$(V-J)_{rest}$", fontsize=12)
ax[1,1].text(0.35,0.925,"2.75 < z < 3.25",ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)
ax[1,1].text(0.25,0.55,"QGs: %d"%(nQ_cat12),color='red',ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)
ax[1,1].text(0.3,0.45,"SFGs: %d"%(NSF_cat12),color='blue',ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)
plt.show()

# UVJ diagram for cat13,14,15 og blank
idcat13 = np.logical_and((data['z_phot'] > 2), (data['z_phot'] < 2.5))
idcat14 = np.logical_and((data['z_phot'] > 2.25), (data['z_phot'] < 2.75))
idcat15 = np.logical_and((data['z_phot'] > 2.5), (data['z_phot'] < 3))

UVcat13 = U[idcat13[0,:]] - V[idcat13[0,:]]
VJcat13 = V[idcat13[0,:]] - J[idcat13[0,:]]
UVcat14 = U[idcat14[0,:]] - V[idcat14[0,:]]
VJcat14 = V[idcat14[0,:]] - J[idcat14[0,:]]
UVcat15 = U[idcat15[0,:]] - V[idcat15[0,:]]
VJcat15 = V[idcat15[0,:]] - J[idcat15[0,:]]

figure, ax = plt.subplots(2, 2, sharex=(True), sharey=(True), figsize=(12,10) )
plt.tight_layout()

# cat13
hb13 = ax[0,0].hexbin(VJcat13, UVcat13, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
cb13 = figure.colorbar(hb13, ax = ax)
cb13.set_label("Number of galaxies", fontsize=15)
ax[0,0].plot((A[0],P31[0]), (A[1],P31[1]), 'black', linestyle="--")
ax[0,0].plot((P31[0],P32[0]), (P31[1],P32[1]), 'black', linestyle="--")
ax[0,0].plot((P32[0],F[0]), (P32[1],F[1]), 'black', linestyle="--")
ax[0,0].set_xlim(0,2.5)
ax[0,0].set_ylim(0,2.5)
ax[0,0].set_ylabel("$(U-V)_{rest}$", fontsize=12)
ax[0,0].text(0.25,0.925,"3 < z < 3.5",ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)
ax[0,0].text(0.25,0.55,"QGs: %d"%(nQ_cat13),color='red',ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)
ax[0,0].text(0.3,0.45,"SFGs: %d"%(NSF_cat13),color='blue',ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)

# cat14
hb14 = ax[0,1].hexbin(VJcat14, UVcat14, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[0,1].plot((A[0],P31[0]), (A[1],P31[1]), 'black', linestyle="--")
ax[0,1].plot((P31[0],P32[0]), (P31[1],P32[1]), 'black', linestyle="--")
ax[0,1].plot((P32[0],F[0]), (P32[1],F[1]), 'black', linestyle="--")
ax[0,1].set_xlim(0,2.5)
ax[0,1].set_ylim(0,2.5)
ax[0,1].text(0.35,0.925,"3.25 < z <3.75",ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)
ax[0,1].text(0.25,0.55,"QGs: %d"%(nQ_cat14),color='red',ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)
ax[0,1].text(0.3,0.45,"SFGs: %d"%(NSF_cat14),color='blue',ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)

# cat15
hb11 = ax[1,0].hexbin(VJcat15, UVcat15, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,0].plot((A[0],P31[0]), (A[1],P31[1]), 'black', linestyle="--")
ax[1,0].plot((P31[0],P32[0]), (P31[1],P32[1]), 'black', linestyle="--")
ax[1,0].plot((P32[0],F[0]), (P32[1],F[1]), 'black', linestyle="--")
ax[1,0].set_xlim(0,2.5)
ax[1,0].set_ylim(0,2.5)
ax[1,0].set_xlabel("$(V-J)_{rest}$", fontsize=12)
ax[1,0].set_ylabel("$(U-V)_{rest}$", fontsize=12)
ax[1,0].text(0.25,0.925,"3.5 < z < 4",ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)
ax[1,0].text(0.25,0.55,"QGs: %d"%(nQ_cat15),color='red',ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)
ax[1,0].text(0.3,0.45,"SFGs: %d"%(NSF_cat15),color='blue',ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)

# fig1_txt = "Flux Density for 3 galaxies at %d < z < %d"
# ax[1,0].text("Flux Density for 3 galaxies at %d < z < %d"%(A[0],A[1]))

# cat16
# hb16 = ax[1,1].hexbin(VJ, UV, vmax = 20, cmap = "binary", mincnt = 0, gridsize=(866,500))
ax[1,1].fill(Qpoly_x,Qpoly_y,facecolor='lightsalmon')
ax[1,1].fill(SFpoly_x,SFpoly_y,facecolor='lightblue')
ax[1,1].set_xlim(0,2.5)
ax[1,1].set_ylim(0,2.5)
ax[1,1].set_xlabel("$(V-J)_{rest}$", fontsize=12)
ax[1,1].text(0.45,0.8,"Quiescent",ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=20)
ax[1,1].text(0.65,0.25,"Star-Forming",ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=20)
plt.show()


#############################

## Mass

# 0 < z < 0.5
m1_cat1 = np.log10(mass[idUVJ_cat1[0,:]])
mass_axe = np.array([8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13])

id_SF_bin1_cat1 = np.logical_and((m1_cat1[~IDS_cat1] > 7.75), (m1_cat1[~IDS_cat1] < 8.25))
id_Q_bin1_cat1 = np.logical_and((m1_cat1[IDS_cat1] > 7.75), (m1_cat1[IDS_cat1] < 8.25))
id_SF_bin2_cat1 = np.logical_and((m1_cat1[~IDS_cat1] > 8.25), (m1_cat1[~IDS_cat1] < 8.75))
id_Q_bin2_cat1 = np.logical_and((m1_cat1[IDS_cat1] > 8.25), (m1_cat1[IDS_cat1] < 8.75))
id_SF_bin3_cat1 = np.logical_and((m1_cat1[~IDS_cat1] > 8.75), (m1_cat1[~IDS_cat1] < 9.25))
id_Q_bin3_cat1 = np.logical_and((m1_cat1[IDS_cat1] > 8.75), (m1_cat1[IDS_cat1] < 9.25))
id_SF_bin4_cat1 = np.logical_and((m1_cat1[~IDS_cat1] > 9.25), (m1_cat1[~IDS_cat1] < 9.75))
id_Q_bin4_cat1 = np.logical_and((m1_cat1[IDS_cat1] > 9.25), (m1_cat1[IDS_cat1] < 9.75))
id_SF_bin5_cat1 = np.logical_and((m1_cat1[~IDS_cat1] > 9.75), (m1_cat1[~IDS_cat1] < 10.25))
id_Q_bin5_cat1 = np.logical_and((m1_cat1[IDS_cat1] > 9.75), (m1_cat1[IDS_cat1] < 10.25))
id_SF_bin6_cat1 = np.logical_and((m1_cat1[~IDS_cat1] > 10.25), (m1_cat1[~IDS_cat1] < 10.75))
id_Q_bin6_cat1 = np.logical_and((m1_cat1[IDS_cat1] > 10.25), (m1_cat1[IDS_cat1] < 10.75))
id_SF_bin7_cat1 = np.logical_and((m1_cat1[~IDS_cat1] > 10.75), (m1_cat1[~IDS_cat1] < 11.25))
id_Q_bin7_cat1 = np.logical_and((m1_cat1[IDS_cat1] > 10.75), (m1_cat1[IDS_cat1] < 11.25))
id_SF_bin8_cat1 = np.logical_and((m1_cat1[~IDS_cat1] > 11.25), (m1_cat1[~IDS_cat1] < 11.75))
id_Q_bin8_cat1 = np.logical_and((m1_cat1[IDS_cat1] > 11.25), (m1_cat1[IDS_cat1] < 11.75))
id_SF_bin9_cat1 = np.logical_and((m1_cat1[~IDS_cat1] > 11.75), (m1_cat1[~IDS_cat1] < 12.25))
id_Q_bin9_cat1 = np.logical_and((m1_cat1[IDS_cat1] > 11.75), (m1_cat1[IDS_cat1] < 12.25))
id_SF_bin10_cat1 = np.logical_and((m1_cat1[~IDS_cat1] > 12.25), (m1_cat1[~IDS_cat1] < 12.75))
id_Q_bin10_cat1 = np.logical_and((m1_cat1[IDS_cat1] > 12.25), (m1_cat1[IDS_cat1] < 12.75))
id_SF_bin11_cat1 = np.logical_and((m1_cat1[~IDS_cat1] > 12.75), (m1_cat1[~IDS_cat1] < 13.25))
id_Q_bin11_cat1 = np.logical_and((m1_cat1[IDS_cat1] > 12.75), (m1_cat1[IDS_cat1] < 13.25))


nSF_cat1 = np.array([np.sum(id_SF_bin1_cat1),np.sum(id_SF_bin2_cat1), np.sum(id_SF_bin3_cat1), np.sum(id_SF_bin4_cat1),
        np.sum(id_SF_bin5_cat1), np.sum(id_SF_bin6_cat1), np.sum(id_SF_bin7_cat1), np.sum(id_SF_bin8_cat1),
        np.sum(id_SF_bin9_cat1), np.sum(id_SF_bin10_cat1), np.sum(id_SF_bin11_cat1)])

NQ_cat1 = np.array([np.sum(id_Q_bin1_cat1), np.sum(id_Q_bin2_cat1), np.sum(id_Q_bin3_cat1), np.sum(id_Q_bin4_cat1),
        np.sum(id_Q_bin5_cat1), np.sum(id_Q_bin6_cat1), np.sum(id_Q_bin7_cat1), np.sum(id_Q_bin8_cat1), 
        np.sum(id_Q_bin9_cat1), np.sum(id_Q_bin10_cat1),np.sum(id_Q_bin11_cat1)])

total_cat1 = nSF_cat1 + NQ_cat1

nSF_error1 = (nSF_cat1/total_cat1)*np.sqrt(1/nSF_cat1+1/total_cat1)
NQ_error1 = (NQ_cat1/total_cat1)*np.sqrt(1/NQ_cat1+1/total_cat1)

# If there is less than 50 galaxies the point is invalid because there is not enougth data
mask1 = nSF_cat1 > 30

#initialize layout of plot
fig1, ax = plt.subplots()

# G1 = ax.errorbar(Lambda, Gflux[:,0,0], yerr = Gfluxerr[:,0,0],fmt = "o", alpha = 0.7,
#             elinewidth = 0.7, capsize = 5, markeredgecolor = "r", markersize = 2)

#add scatterplot for cat1
nSF = ax.errorbar(mass_axe[mask1],nSF_cat1[mask1]/total_cat1[mask1], yerr = nSF_error1[mask1], fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "b", markersize = 2)
nQ = ax.errorbar(mass_axe[mask1],NQ_cat1[mask1]/total_cat1[mask1], yerr = NQ_error1[mask1], fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "r", markersize = 2)

ax.legend([nSF,nQ],['nSF','nQ'],loc="upper right")
ax.set_xlim(8,13)
ax.set_ylim(0,1)
ax.set_xlabel("$log_{10}(M_{\odot})$")
ax.set_ylabel("Fraction of galaxies")
ax.set_title("Fraction of galaxies with mass at 0 < z < 0.5")

# 0.5 < z < 1.0
m1_cat2 = np.log10(mass[idUVJ_cat3[0,:]])

id_SF_bin1_cat2 = np.logical_and((m1_cat2[~IDS_cat3] > 7.75), (m1_cat2[~IDS_cat3] < 8.25))
id_Q_bin1_cat2 = np.logical_and((m1_cat2[IDS_cat3] > 7.75), (m1_cat2[IDS_cat3] < 8.25))
id_SF_bin2_cat2 = np.logical_and((m1_cat2[~IDS_cat3] > 8.25), (m1_cat2[~IDS_cat3] < 8.75))
id_Q_bin2_cat2 = np.logical_and((m1_cat2[IDS_cat3] > 8.25), (m1_cat2[IDS_cat3] < 8.75))
id_SF_bin3_cat2 = np.logical_and((m1_cat2[~IDS_cat3] > 8.75), (m1_cat2[~IDS_cat3] < 9.25))
id_Q_bin3_cat2 = np.logical_and((m1_cat2[IDS_cat3] > 8.75), (m1_cat2[IDS_cat3] < 9.25))
id_SF_bin4_cat2 = np.logical_and((m1_cat2[~IDS_cat3] > 9.25), (m1_cat2[~IDS_cat3] < 9.75))
id_Q_bin4_cat2 = np.logical_and((m1_cat2[IDS_cat3] > 9.25), (m1_cat2[IDS_cat3] < 9.75))
id_SF_bin5_cat2 = np.logical_and((m1_cat2[~IDS_cat3] > 9.75), (m1_cat2[~IDS_cat3] < 10.25))
id_Q_bin5_cat2 = np.logical_and((m1_cat2[IDS_cat3] > 9.75), (m1_cat2[IDS_cat3] < 10.25))
id_SF_bin6_cat2 = np.logical_and((m1_cat2[~IDS_cat3] > 10.25), (m1_cat2[~IDS_cat3] < 10.75))
id_Q_bin6_cat2 = np.logical_and((m1_cat2[IDS_cat3] > 10.25), (m1_cat2[IDS_cat3] < 10.75))
id_SF_bin7_cat2 = np.logical_and((m1_cat2[~IDS_cat3] > 10.75), (m1_cat2[~IDS_cat3] < 11.25))
id_Q_bin7_cat2 = np.logical_and((m1_cat2[IDS_cat3] > 10.75), (m1_cat2[IDS_cat3] < 11.25))
id_SF_bin8_cat2 = np.logical_and((m1_cat2[~IDS_cat3] > 11.25), (m1_cat2[~IDS_cat3] < 11.75))
id_Q_bin8_cat2 = np.logical_and((m1_cat2[IDS_cat3] > 11.25), (m1_cat2[IDS_cat3] < 11.75))
id_SF_bin9_cat2 = np.logical_and((m1_cat2[~IDS_cat3] > 11.75), (m1_cat2[~IDS_cat3] < 12.25))
id_Q_bin9_cat2 = np.logical_and((m1_cat2[IDS_cat3] > 11.75), (m1_cat2[IDS_cat3] < 12.25))
id_SF_bin10_cat2 = np.logical_and((m1_cat2[~IDS_cat3] > 12.25), (m1_cat2[~IDS_cat3] < 12.75))
id_Q_bin10_cat2 = np.logical_and((m1_cat2[IDS_cat3] > 12.25), (m1_cat2[IDS_cat3] < 12.75))
id_SF_bin11_cat2 = np.logical_and((m1_cat2[~IDS_cat3] > 12.75), (m1_cat2[~IDS_cat3] < 13.25))
id_Q_bin11_cat2 = np.logical_and((m1_cat2[IDS_cat3] > 12.75), (m1_cat2[IDS_cat3] < 13.25))


nSF_cat2 = np.array([np.sum(id_SF_bin1_cat2),np.sum(id_SF_bin2_cat2), np.sum(id_SF_bin3_cat2), np.sum(id_SF_bin4_cat2),
        np.sum(id_SF_bin5_cat2), np.sum(id_SF_bin6_cat2), np.sum(id_SF_bin7_cat2), np.sum(id_SF_bin8_cat2),
        np.sum(id_SF_bin9_cat2), np.sum(id_SF_bin10_cat2), np.sum(id_SF_bin11_cat2)])

NQ_cat2 = np.array([np.sum(id_Q_bin1_cat2), np.sum(id_Q_bin2_cat2), np.sum(id_Q_bin3_cat2), np.sum(id_Q_bin4_cat2),
        np.sum(id_Q_bin5_cat2), np.sum(id_Q_bin6_cat2), np.sum(id_Q_bin7_cat2), np.sum(id_Q_bin8_cat2), 
        np.sum(id_Q_bin9_cat2), np.sum(id_Q_bin10_cat2),np.sum(id_Q_bin11_cat2)])

total_cat2 = nSF_cat2 + NQ_cat2

nSF_error2 = (nSF_cat2/total_cat2)*np.sqrt(1/nSF_cat2 + 1/total_cat2)
NQ_error2 = (NQ_cat2/total_cat2)*np.sqrt(1/NQ_cat2 + 1/total_cat2)

mask2 = nSF_cat2 > 30

#initialize layout of plot
fig1, ax = plt.subplots()

#add scatterplot for cat1
nSF = ax.errorbar(mass_axe[mask2],nSF_cat2[mask2]/total_cat2[mask2], yerr = nSF_error2[mask2], fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "b", markersize = 2, label="nSF")
nQ = ax.errorbar(mass_axe[mask2],NQ_cat2[mask2]/total_cat2[mask2], yerr = NQ_error2[mask2], fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "r", markersize = 2, label="nQ")

ax.legend([nSF,nQ],['nSF','nQ'],loc="upper right")
ax.set_xlim(8,13)
ax.set_ylim(0,1)
ax.set_xlabel("$log_{10}(M_{\odot})$")
ax.set_ylabel("Fraction of galaxies")
ax.set_title("Fraction of galaxies with mass at 0.5 < z < 1")

# 1 < z < 1.5
m1_cat3 = np.log10(mass[idUVJ_cat5[0,:]])

id_SF_bin1_cat3 = np.logical_and((m1_cat3[~IDS_cat5] > 7.75), (m1_cat3[~IDS_cat5] < 8.25))
id_Q_bin1_cat3 = np.logical_and((m1_cat3[IDS_cat5] > 7.75), (m1_cat3[IDS_cat5] < 8.25))
id_SF_bin2_cat3 = np.logical_and((m1_cat3[~IDS_cat5] > 8.25), (m1_cat3[~IDS_cat5] < 8.75))
id_Q_bin2_cat3 = np.logical_and((m1_cat3[IDS_cat5] > 8.25), (m1_cat3[IDS_cat5] < 8.75))
id_SF_bin3_cat3 = np.logical_and((m1_cat3[~IDS_cat5] > 8.75), (m1_cat3[~IDS_cat5] < 9.25))
id_Q_bin3_cat3 = np.logical_and((m1_cat3[IDS_cat5] > 8.75), (m1_cat3[IDS_cat5] < 9.25))
id_SF_bin4_cat3 = np.logical_and((m1_cat3[~IDS_cat5] > 9.25), (m1_cat3[~IDS_cat5] < 9.75))
id_Q_bin4_cat3 = np.logical_and((m1_cat3[IDS_cat5] > 9.25), (m1_cat3[IDS_cat5] < 9.75))
id_SF_bin5_cat3 = np.logical_and((m1_cat3[~IDS_cat5] > 9.75), (m1_cat3[~IDS_cat5] < 10.25))
id_Q_bin5_cat3 = np.logical_and((m1_cat3[IDS_cat5] > 9.75), (m1_cat3[IDS_cat5] < 10.25))
id_SF_bin6_cat3 = np.logical_and((m1_cat3[~IDS_cat5] > 10.25), (m1_cat3[~IDS_cat5] < 10.75))
id_Q_bin6_cat3 = np.logical_and((m1_cat3[IDS_cat5] > 10.25), (m1_cat3[IDS_cat5] < 10.75))
id_SF_bin7_cat3 = np.logical_and((m1_cat3[~IDS_cat5] > 10.75), (m1_cat3[~IDS_cat5] < 11.25))
id_Q_bin7_cat3 = np.logical_and((m1_cat3[IDS_cat5] > 10.75), (m1_cat3[IDS_cat5] < 11.25))
id_SF_bin8_cat3 = np.logical_and((m1_cat3[~IDS_cat5] > 11.25), (m1_cat3[~IDS_cat5] < 11.75))
id_Q_bin8_cat3 = np.logical_and((m1_cat3[IDS_cat5] > 11.25), (m1_cat3[IDS_cat5] < 11.75))
id_SF_bin9_cat3 = np.logical_and((m1_cat3[~IDS_cat5] > 11.75), (m1_cat3[~IDS_cat5] < 12.25))
id_Q_bin9_cat3 = np.logical_and((m1_cat3[IDS_cat5] > 11.75), (m1_cat3[IDS_cat5] < 12.25))
id_SF_bin10_cat3 = np.logical_and((m1_cat3[~IDS_cat5] > 12.25), (m1_cat3[~IDS_cat5] < 12.75))
id_Q_bin10_cat3 = np.logical_and((m1_cat3[IDS_cat5] > 12.25), (m1_cat3[IDS_cat5] < 12.75))
id_SF_bin11_cat3 = np.logical_and((m1_cat3[~IDS_cat5] > 12.75), (m1_cat3[~IDS_cat5] < 13.25))
id_Q_bin11_cat3 = np.logical_and((m1_cat3[IDS_cat5] > 12.75), (m1_cat3[IDS_cat5] < 13.25))


nSF_cat3 = np.array([np.sum(id_SF_bin1_cat3),np.sum(id_SF_bin2_cat3), np.sum(id_SF_bin3_cat3), np.sum(id_SF_bin4_cat3),
        np.sum(id_SF_bin5_cat3), np.sum(id_SF_bin6_cat3), np.sum(id_SF_bin7_cat3), np.sum(id_SF_bin8_cat3),
        np.sum(id_SF_bin9_cat3), np.sum(id_SF_bin10_cat3), np.sum(id_SF_bin11_cat3)])

NQ_cat3 = np.array([np.sum(id_Q_bin1_cat3), np.sum(id_Q_bin2_cat3), np.sum(id_Q_bin3_cat3), np.sum(id_Q_bin4_cat3),
        np.sum(id_Q_bin5_cat3), np.sum(id_Q_bin6_cat3), np.sum(id_Q_bin7_cat3), np.sum(id_Q_bin8_cat3), 
        np.sum(id_Q_bin9_cat3), np.sum(id_Q_bin10_cat3),np.sum(id_Q_bin11_cat3)])

total_cat3 = nSF_cat3 + NQ_cat3

nSF_error3 = (nSF_cat3/total_cat3)*np.sqrt(1/nSF_cat3 + 1/total_cat3)
NQ_error3 = (NQ_cat3/total_cat3)*np.sqrt(1/NQ_cat3 + 1/total_cat3)

mask3 = nSF_cat3 > 30

#initialize layout of plot
fig1, ax = plt.subplots()

#add scatterplot for cat1
nSF = ax.errorbar(mass_axe[mask3],nSF_cat3[mask3]/total_cat3[mask3], yerr = nSF_error3[mask3], fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "b", markersize = 2)
nQ = ax.errorbar(mass_axe[mask3],NQ_cat3[mask3]/total_cat3[mask3], yerr = NQ_error3[mask3], fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "r", markersize = 2)
ax.legend([nSF,nQ],['nSF','nQ'],loc="upper right")
ax.set_xlim(8,13)
ax.set_ylim(0,1)
ax.set_xlabel("$log_{10}(M_{\odot})$")
ax.set_ylabel("Fraction of galaxies")
ax.set_title("Fraction of galaxies with mass at 1 < z < 1.5")

# 1.5 < z < 2.5
# Lower 1.5 < z < 2 
m1_cat4_lower = np.log10(mass[idUVJ_cat7[0,:]])
m1_cat4_upper = np.log10(mass[idUVJ_cat9[0,:]])

id_SF_bin1_cat4 = np.logical_and((m1_cat4_lower[~IDS_cat7] > 7.75), (m1_cat4_lower[~IDS_cat7] < 8.25))
id_Q_bin1_cat4 = np.logical_and((m1_cat4_lower[IDS_cat7] > 7.75), ( m1_cat4_lower[IDS_cat7] < 8.25))
id_SF_bin2_cat4 = np.logical_and((m1_cat4_lower[~IDS_cat7] > 8.25 ), ( m1_cat4_lower[~IDS_cat7] < 8.75))
id_Q_bin2_cat4 = np.logical_and((m1_cat4_lower[IDS_cat7] > 8.25 ), ( m1_cat4_lower[IDS_cat7] < 8.75))
id_SF_bin3_cat4 = np.logical_and((m1_cat4_lower[~IDS_cat7] > 8.75 ), ( m1_cat4_lower[~IDS_cat7] < 9.25))
id_Q_bin3_cat4 = np.logical_and((m1_cat4_lower[IDS_cat7] > 8.75), ( m1_cat4_lower[IDS_cat7] < 9.25))
id_SF_bin4_cat4 = np.logical_and((m1_cat4_lower[~IDS_cat7] > 9.25), ( m1_cat4_lower[~IDS_cat7] < 9.75))
id_Q_bin4_cat4 = np.logical_and((m1_cat4_lower[IDS_cat7] > 9.25), ( m1_cat4_lower[IDS_cat7] < 9.75))
id_SF_bin5_cat4 = np.logical_and((m1_cat4_lower[~IDS_cat7] > 9.75), ( m1_cat4_lower[~IDS_cat7] < 10.25))
id_Q_bin5_cat4 = np.logical_and((m1_cat4_lower[IDS_cat7] > 9.75), ( m1_cat4_lower[IDS_cat7] < 10.25))
id_SF_bin6_cat4 = np.logical_and((m1_cat4_lower[~IDS_cat7] > 10.25), ( m1_cat4_lower[~IDS_cat7] < 10.75))
id_Q_bin6_cat4 = np.logical_and((m1_cat4_lower[IDS_cat7] > 10.25), ( m1_cat4_lower[IDS_cat7] < 10.75))
id_SF_bin7_cat4 = np.logical_and((m1_cat4_lower[~IDS_cat7] > 10.75), ( m1_cat4_lower[~IDS_cat7] < 11.25))
id_Q_bin7_cat4 = np.logical_and((m1_cat4_lower[IDS_cat7] > 10.75), ( m1_cat4_lower[IDS_cat7] < 11.25))
id_SF_bin8_cat4 = np.logical_and((m1_cat4_lower[~IDS_cat7] > 11.25), ( m1_cat4_lower[~IDS_cat7] < 11.75))
id_Q_bin8_cat4 = np.logical_and((m1_cat4_lower[IDS_cat7] > 11.25), ( m1_cat4_lower[IDS_cat7] < 11.75))
id_SF_bin9_cat4 = np.logical_and((m1_cat4_lower[~IDS_cat7] > 11.75), ( m1_cat4_lower[~IDS_cat7] < 12.25))
id_Q_bin9_cat4 = np.logical_and((m1_cat4_lower[IDS_cat7] > 11.75), ( m1_cat4_lower[IDS_cat7] < 12.25))
id_SF_bin10_cat4 = np.logical_and((m1_cat4_lower[~IDS_cat7] > 12.25), ( m1_cat4_lower[~IDS_cat7] < 12.75))
id_Q_bin10_cat4 = np.logical_and((m1_cat4_lower[IDS_cat7] > 12.25), ( m1_cat4_lower[IDS_cat7] < 12.75))
id_SF_bin11_cat4 = np.logical_and((m1_cat4_lower[~IDS_cat7] > 12.75), ( m1_cat4_lower[~IDS_cat7] < 13.25))
id_Q_bin11_cat4 = np.logical_and((m1_cat4_lower[IDS_cat7] > 12.75), ( m1_cat4_lower[IDS_cat7] < 13.25))

# upper
# 2 < z < 2.5
id_SF_bin1_cat4_upper = np.logical_and((m1_cat4_upper[~IDS_cat9] > 7.75), (m1_cat4_upper[~IDS_cat9] < 8.25))
id_Q_bin1_cat4_upper = np.logical_and((m1_cat4_upper[IDS_cat9] > 7.75), ( m1_cat4_upper[IDS_cat9] < 8.25))
id_SF_bin2_cat4_upper = np.logical_and((m1_cat4_upper[~IDS_cat9] > 8.25), ( m1_cat4_upper[~IDS_cat9] < 8.75))
id_Q_bin2_cat4_upper = np.logical_and((m1_cat4_upper[IDS_cat9] > 8.25), ( m1_cat4_upper[IDS_cat9] < 8.75))
id_SF_bin3_cat4_upper = np.logical_and((m1_cat4_upper[~IDS_cat9] > 8.75), ( m1_cat4_upper[~IDS_cat9] < 9.25))
id_Q_bin3_cat4_upper = np.logical_and((m1_cat4_upper[IDS_cat9] > 8.75), ( m1_cat4_upper[IDS_cat9] < 9.25))
id_SF_bin4_cat4_upper = np.logical_and((m1_cat4_upper[~IDS_cat9] > 9.25), ( m1_cat4_upper[~IDS_cat9] < 9.75))
id_Q_bin4_cat4_upper = np.logical_and((m1_cat4_upper[IDS_cat9] > 9.25), ( m1_cat4_upper[IDS_cat9] < 9.75))
id_SF_bin5_cat4_upper = np.logical_and((m1_cat4_upper[~IDS_cat9] > 9.75), ( m1_cat4_upper[~IDS_cat9] < 10.25))
id_Q_bin5_cat4_upper = np.logical_and((m1_cat4_upper[IDS_cat9] > 9.75), ( m1_cat4_upper[IDS_cat9] < 10.25))
id_SF_bin6_cat4_upper = np.logical_and((m1_cat4_upper[~IDS_cat9] > 10.25), ( m1_cat4_upper[~IDS_cat9] < 10.75))
id_Q_bin6_cat4_upper = np.logical_and((m1_cat4_upper[IDS_cat9] > 10.25), ( m1_cat4_upper[IDS_cat9] < 10.75))
id_SF_bin7_cat4_upper = np.logical_and((m1_cat4_upper[~IDS_cat9] > 10.75), ( m1_cat4_upper[~IDS_cat9] < 11.25))
id_Q_bin7_cat4_upper = np.logical_and((m1_cat4_upper[IDS_cat9] > 10.75), ( m1_cat4_upper[IDS_cat9] < 11.25))
id_SF_bin8_cat4_upper = np.logical_and((m1_cat4_upper[~IDS_cat9] > 11.25), ( m1_cat4_upper[~IDS_cat9] < 11.75))
id_Q_bin8_cat4_upper = np.logical_and((m1_cat4_upper[IDS_cat9] > 11.25), ( m1_cat4_upper[IDS_cat9] < 11.75))
id_SF_bin9_cat4_upper = np.logical_and((m1_cat4_upper[~IDS_cat9] > 11.75), ( m1_cat4_upper[~IDS_cat9] < 12.25))
id_Q_bin9_cat4_upper = np.logical_and((m1_cat4_upper[IDS_cat9] > 11.75), ( m1_cat4_upper[IDS_cat9] < 12.25))
id_SF_bin10_cat4_upper = np.logical_and((m1_cat4_upper[~IDS_cat9] > 12.25), ( m1_cat4_upper[~IDS_cat9] < 12.75))
id_Q_bin10_cat4_upper = np.logical_and((m1_cat4_upper[IDS_cat9] > 12.25), ( m1_cat4_upper[IDS_cat9] < 12.75))
id_SF_bin11_cat4_upper = np.logical_and((m1_cat4_upper[~IDS_cat9] > 12.75), ( m1_cat4_upper[~IDS_cat9] < 13.25))
id_Q_bin11_cat4_upper = np.logical_and((m1_cat4_upper[IDS_cat9] > 12.75), ( m1_cat4_upper[IDS_cat9] < 13.25))

nSF_cat4 = np.array([np.sum(id_SF_bin1_cat4),np.sum(id_SF_bin2_cat4), np.sum(id_SF_bin3_cat4), np.sum(id_SF_bin4_cat4),
        np.sum(id_SF_bin5_cat4), np.sum(id_SF_bin6_cat4), np.sum(id_SF_bin7_cat4), np.sum(id_SF_bin8_cat4),
        np.sum(id_SF_bin9_cat4), np.sum(id_SF_bin10_cat4), np.sum(id_SF_bin11_cat4)]) + np.array([np.sum(id_SF_bin1_cat4_upper),np.sum(id_SF_bin2_cat4_upper), np.sum(id_SF_bin3_cat4_upper), np.sum(id_SF_bin4_cat4_upper),
                np.sum(id_SF_bin5_cat4_upper), np.sum(id_SF_bin6_cat4_upper), np.sum(id_SF_bin7_cat4_upper), np.sum(id_SF_bin8_cat4_upper),
                np.sum(id_SF_bin9_cat4_upper), np.sum(id_SF_bin10_cat4_upper), np.sum(id_SF_bin11_cat4_upper)])

NQ_cat4 = np.array([np.sum(id_Q_bin1_cat4), np.sum(id_Q_bin2_cat4), np.sum(id_Q_bin3_cat4), np.sum(id_Q_bin4_cat4),
        np.sum(id_Q_bin5_cat4), np.sum(id_Q_bin6_cat4), np.sum(id_Q_bin7_cat4), np.sum(id_Q_bin8_cat4), 
        np.sum(id_Q_bin9_cat4), np.sum(id_Q_bin10_cat4),np.sum(id_Q_bin11_cat4)]) + np.array([np.sum(id_Q_bin1_cat4_upper), np.sum(id_Q_bin2_cat4_upper), np.sum(id_Q_bin3_cat4_upper), np.sum(id_Q_bin4_cat4_upper),
                np.sum(id_Q_bin5_cat4_upper), np.sum(id_Q_bin6_cat4_upper), np.sum(id_Q_bin7_cat4_upper), np.sum(id_Q_bin8_cat4_upper), 
                np.sum(id_Q_bin9_cat4_upper), np.sum(id_Q_bin10_cat4_upper),np.sum(id_Q_bin11_cat4_upper)])

total_cat4 = nSF_cat4 + NQ_cat4

nSF_error4 = (nSF_cat4/total_cat4)*np.sqrt(1/nSF_cat4 + 1/total_cat4)
NQ_error4 = (NQ_cat4/total_cat4)*np.sqrt(1/NQ_cat4 + 1/total_cat4)

mask4 = nSF_cat4 > 30

#initialize layout of plot
fig1, ax = plt.subplots()

#add scatterplot for category
SF = ax.errorbar(mass_axe[mask4],nSF_cat4[mask4]/total_cat4[mask4], yerr = nSF_error4[mask4], fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "b", markersize = 2, label='nSF')
nQ = ax.errorbar(mass_axe[mask4],NQ_cat4[mask4]/total_cat4[mask4], yerr = NQ_error4[mask4], fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "r", markersize = 2, label='nQ')
ax.legend([nSF,nQ],['nSF','nQ'],loc="upper right")
ax.set_xlim(8,13)
ax.set_ylim(0,1)
ax.set_xlabel("$log_{10}(M_{\odot})$")
ax.set_ylabel("Fraction of galaxies")
ax.set_title("Fraction of galaxies with mass at 1.5 < z < 2.5")

# 2.5 < z < 4
m1_cat5_lower = np.log10(mass[idUVJ_cat11[0,:]])
m1_cat5_mid = np.log10(mass[idUVJ_cat13[0,:]])
m1_cat5_upper = np.log10(mass[idUVJ_cat15[0,:]])

# 2.5 < z < 3
# lower
id_SF_bin1_cat5_lower = np.logical_and((m1_cat5_lower[~IDS_cat11] > 7.75), (m1_cat5_lower[~IDS_cat11] < 8.25))
id_Q_bin1_cat5_lower = np.logical_and((m1_cat5_lower[IDS_cat11] > 7.75), ( m1_cat5_lower[IDS_cat11] < 8.25))
id_SF_bin2_cat5_lower = np.logical_and((m1_cat5_lower[~IDS_cat11] > 8.25 ), ( m1_cat5_lower[~IDS_cat11] < 8.75))
id_Q_bin2_cat5_lower = np.logical_and((m1_cat5_lower[IDS_cat11] > 8.25 ), ( m1_cat5_lower[IDS_cat11] < 8.75))
id_SF_bin3_cat5_lower = np.logical_and((m1_cat5_lower[~IDS_cat11] > 8.75 ), ( m1_cat5_lower[~IDS_cat11] < 9.25))
id_Q_bin3_cat5_lower = np.logical_and((m1_cat5_lower[IDS_cat11] > 8.75), ( m1_cat5_lower[IDS_cat11] < 9.25))
id_SF_bin4_cat5_lower = np.logical_and((m1_cat5_lower[~IDS_cat11] > 9.25), ( m1_cat5_lower[~IDS_cat11] < 9.75))
id_Q_bin4_cat5_lower = np.logical_and((m1_cat5_lower[IDS_cat11] > 9.25), ( m1_cat5_lower[IDS_cat11] < 9.75))
id_SF_bin5_cat5_lower = np.logical_and((m1_cat5_lower[~IDS_cat11] > 9.75), ( m1_cat5_lower[~IDS_cat11] < 10.25))
id_Q_bin5_cat5_lower = np.logical_and((m1_cat5_lower[IDS_cat11] > 9.75), ( m1_cat5_lower[IDS_cat11] < 10.25))
id_SF_bin6_cat5_lower = np.logical_and((m1_cat5_lower[~IDS_cat11] > 10.25), ( m1_cat5_lower[~IDS_cat11] < 10.75))
id_Q_bin6_cat5_lower = np.logical_and((m1_cat5_lower[IDS_cat11] > 10.25), ( m1_cat5_lower[IDS_cat11] < 10.75))
id_SF_bin7_cat5_lower = np.logical_and((m1_cat5_lower[~IDS_cat11] > 10.75), ( m1_cat5_lower[~IDS_cat11] < 11.25))
id_Q_bin7_cat5_lower = np.logical_and((m1_cat5_lower[IDS_cat11] > 10.75), ( m1_cat5_lower[IDS_cat11] < 11.25))
id_SF_bin8_cat5_lower = np.logical_and((m1_cat5_lower[~IDS_cat11] > 11.25), ( m1_cat5_lower[~IDS_cat11] < 11.75))
id_Q_bin8_cat5_lower = np.logical_and((m1_cat5_lower[IDS_cat11] > 11.25), ( m1_cat5_lower[IDS_cat11] < 11.75))
id_SF_bin9_cat5_lower = np.logical_and((m1_cat5_lower[~IDS_cat11] > 11.75), ( m1_cat5_lower[~IDS_cat11] < 12.25))
id_Q_bin9_cat5_lower = np.logical_and((m1_cat5_lower[IDS_cat11] > 11.75), ( m1_cat5_lower[IDS_cat11] < 12.25))
id_SF_bin10_cat5_lower = np.logical_and((m1_cat5_lower[~IDS_cat11] > 12.25), ( m1_cat5_lower[~IDS_cat11] < 12.75))
id_Q_bin10_cat5_lower = np.logical_and((m1_cat5_lower[IDS_cat11] > 12.25), ( m1_cat5_lower[IDS_cat11] < 12.75))
id_SF_bin11_cat5_lower = np.logical_and((m1_cat5_lower[~IDS_cat11] > 12.75), ( m1_cat5_lower[~IDS_cat11] < 13.25))
id_Q_bin11_cat5_lower = np.logical_and((m1_cat5_lower[IDS_cat11] > 12.75), ( m1_cat5_lower[IDS_cat11] < 13.25))

# mid
# 3 < z < 3.5
id_SF_bin1_cat5_mid = np.logical_and((m1_cat5_mid[~IDS_cat13] > 7.75), (m1_cat5_mid[~IDS_cat13] < 8.25))
id_Q_bin1_cat5_mid = np.logical_and((m1_cat5_mid[IDS_cat13] > 7.75), ( m1_cat5_mid[IDS_cat13] < 8.25))
id_SF_bin2_cat5_mid = np.logical_and((m1_cat5_mid[~IDS_cat13] > 8.25), ( m1_cat5_mid[~IDS_cat13] < 8.75))
id_Q_bin2_cat5_mid = np.logical_and((m1_cat5_mid[IDS_cat13] > 8.25), ( m1_cat5_mid[IDS_cat13] < 8.75))
id_SF_bin3_cat5_mid = np.logical_and((m1_cat5_mid[~IDS_cat13] > 8.75), ( m1_cat5_mid[~IDS_cat13] < 9.25))
id_Q_bin3_cat5_mid = np.logical_and((m1_cat5_mid[IDS_cat13] > 8.75), ( m1_cat5_mid[IDS_cat13] < 9.25))
id_SF_bin4_cat5_mid = np.logical_and((m1_cat5_mid[~IDS_cat13] > 9.25), ( m1_cat5_mid[~IDS_cat13] < 9.75))
id_Q_bin4_cat5_mid = np.logical_and((m1_cat5_mid[IDS_cat13] > 9.25), ( m1_cat5_mid[IDS_cat13] < 9.75))
id_SF_bin5_cat5_mid = np.logical_and((m1_cat5_mid[~IDS_cat13] > 9.75), ( m1_cat5_mid[~IDS_cat13] < 10.25))
id_Q_bin5_cat5_mid = np.logical_and((m1_cat5_mid[IDS_cat13] > 9.75), ( m1_cat5_mid[IDS_cat13] < 10.25))
id_SF_bin6_cat5_mid = np.logical_and((m1_cat5_mid[~IDS_cat13] > 10.25), ( m1_cat5_mid[~IDS_cat13] < 10.75))
id_Q_bin6_cat5_mid = np.logical_and((m1_cat5_mid[IDS_cat13] > 10.25), ( m1_cat5_mid[IDS_cat13] < 10.75))
id_SF_bin7_cat5_mid = np.logical_and((m1_cat5_mid[~IDS_cat13] > 10.75), ( m1_cat5_mid[~IDS_cat13] < 11.25))
id_Q_bin7_cat5_mid = np.logical_and((m1_cat5_mid[IDS_cat13] > 10.75), ( m1_cat5_mid[IDS_cat13] < 11.25))
id_SF_bin8_cat5_mid = np.logical_and((m1_cat5_mid[~IDS_cat13] > 11.25), ( m1_cat5_mid[~IDS_cat13] < 11.75))
id_Q_bin8_cat5_mid = np.logical_and((m1_cat5_mid[IDS_cat13] > 11.25), ( m1_cat5_mid[IDS_cat13] < 11.75))
id_SF_bin9_cat5_mid = np.logical_and((m1_cat5_mid[~IDS_cat13] > 11.75), ( m1_cat5_mid[~IDS_cat13] < 12.25))
id_Q_bin9_cat5_mid = np.logical_and((m1_cat5_mid[IDS_cat13] > 11.75), ( m1_cat5_mid[IDS_cat13] < 12.25))
id_SF_bin10_cat5_mid = np.logical_and((m1_cat5_mid[~IDS_cat13] > 12.25), ( m1_cat5_mid[~IDS_cat13] < 12.75))
id_Q_bin10_cat5_mid = np.logical_and((m1_cat5_mid[IDS_cat13] > 12.25), ( m1_cat5_mid[IDS_cat13] < 12.75))
id_SF_bin11_cat5_mid = np.logical_and((m1_cat5_mid[~IDS_cat13] > 12.75), ( m1_cat5_mid[~IDS_cat13] < 13.25))
id_Q_bin11_cat5_mid = np.logical_and((m1_cat5_mid[IDS_cat13] > 12.75), ( m1_cat5_mid[IDS_cat13] < 13.25))

# upper
# 3.5 < z < 4
id_SF_bin1_cat5_upper = np.logical_and((m1_cat5_upper[~IDS_cat15] > 7.75), (m1_cat5_upper[~IDS_cat15] < 8.25))
id_Q_bin1_cat5_upper = np.logical_and((m1_cat5_upper[IDS_cat15] > 7.75), ( m1_cat5_upper[IDS_cat15] < 8.25))
id_SF_bin2_cat5_upper = np.logical_and((m1_cat5_upper[~IDS_cat15] > 8.25), ( m1_cat5_upper[~IDS_cat15] < 8.75))
id_Q_bin2_cat5_upper = np.logical_and((m1_cat5_upper[IDS_cat15] > 8.25), ( m1_cat5_upper[IDS_cat15] < 8.75))
id_SF_bin3_cat5_upper = np.logical_and((m1_cat5_upper[~IDS_cat15] > 8.75), ( m1_cat5_upper[~IDS_cat15] < 9.25))
id_Q_bin3_cat5_upper = np.logical_and((m1_cat5_upper[IDS_cat15] > 8.75), ( m1_cat5_upper[IDS_cat15] < 9.25))
id_SF_bin4_cat5_upper = np.logical_and((m1_cat5_upper[~IDS_cat15] > 9.25), ( m1_cat5_upper[~IDS_cat15] < 9.75))
id_Q_bin4_cat5_upper = np.logical_and((m1_cat5_upper[IDS_cat15] > 9.25), ( m1_cat5_upper[IDS_cat15] < 9.75))
id_SF_bin5_cat5_upper = np.logical_and((m1_cat5_upper[~IDS_cat15] > 9.75), ( m1_cat5_upper[~IDS_cat15] < 10.25))
id_Q_bin5_cat5_upper = np.logical_and((m1_cat5_upper[IDS_cat15] > 9.75), ( m1_cat5_upper[IDS_cat15] < 10.25))
id_SF_bin6_cat5_upper = np.logical_and((m1_cat5_upper[~IDS_cat15] > 10.25), ( m1_cat5_upper[~IDS_cat15] < 10.75))
id_Q_bin6_cat5_upper = np.logical_and((m1_cat5_upper[IDS_cat15] > 10.25), ( m1_cat5_upper[IDS_cat15] < 10.75))
id_SF_bin7_cat5_upper = np.logical_and((m1_cat5_upper[~IDS_cat15] > 10.75), ( m1_cat5_upper[~IDS_cat15] < 11.25))
id_Q_bin7_cat5_upper = np.logical_and((m1_cat5_upper[IDS_cat15] > 10.75), ( m1_cat5_upper[IDS_cat15] < 11.25))
id_SF_bin8_cat5_upper = np.logical_and((m1_cat5_upper[~IDS_cat15] > 11.25), ( m1_cat5_upper[~IDS_cat15] < 11.75))
id_Q_bin8_cat5_upper = np.logical_and((m1_cat5_upper[IDS_cat15] > 11.25), ( m1_cat5_upper[IDS_cat15] < 11.75))
id_SF_bin9_cat5_upper = np.logical_and((m1_cat5_upper[~IDS_cat15] > 11.75), ( m1_cat5_upper[~IDS_cat15] < 12.25))
id_Q_bin9_cat5_upper = np.logical_and((m1_cat5_upper[IDS_cat15] > 11.75), ( m1_cat5_upper[IDS_cat15] < 12.25))
id_SF_bin10_cat5_upper = np.logical_and((m1_cat5_upper[~IDS_cat15] > 12.25), ( m1_cat5_upper[~IDS_cat15] < 12.75))
id_Q_bin10_cat5_upper = np.logical_and((m1_cat5_upper[IDS_cat15] > 12.25), ( m1_cat5_upper[IDS_cat15] < 12.75))
id_SF_bin11_cat5_upper = np.logical_and((m1_cat5_upper[~IDS_cat15] > 12.75), ( m1_cat5_upper[~IDS_cat15] < 13.25))
id_Q_bin11_cat5_upper = np.logical_and((m1_cat5_upper[IDS_cat15] > 12.75), ( m1_cat5_upper[IDS_cat15] < 13.25))


nSF_cat5 = np.array([np.sum(id_SF_bin1_cat5_lower),np.sum(id_SF_bin2_cat5_lower), np.sum(id_SF_bin3_cat5_lower), np.sum(id_SF_bin4_cat5_lower),
        np.sum(id_SF_bin5_cat5_lower), np.sum(id_SF_bin6_cat5_lower), np.sum(id_SF_bin7_cat5_lower), np.sum(id_SF_bin8_cat5_lower),
        np.sum(id_SF_bin9_cat5_lower), np.sum(id_SF_bin10_cat5_lower), np.sum(id_SF_bin11_cat5_lower)]) + np.array([np.sum(id_SF_bin1_cat5_mid),np.sum(id_SF_bin2_cat5_mid), np.sum(id_SF_bin3_cat5_mid), np.sum(id_SF_bin4_cat5_mid),
                np.sum(id_SF_bin5_cat5_mid), np.sum(id_SF_bin6_cat5_mid), np.sum(id_SF_bin7_cat5_mid), np.sum(id_SF_bin8_cat5_mid),
                np.sum(id_SF_bin9_cat5_mid), np.sum(id_SF_bin10_cat5_mid), np.sum(id_SF_bin11_cat5_mid)]) + np.array([np.sum(id_SF_bin1_cat5_upper),np.sum(id_SF_bin2_cat5_upper), np.sum(id_SF_bin3_cat5_upper), np.sum(id_SF_bin4_cat5_upper),
                np.sum(id_SF_bin5_cat5_upper), np.sum(id_SF_bin6_cat5_upper), np.sum(id_SF_bin7_cat5_upper), np.sum(id_SF_bin8_cat5_upper),
                np.sum(id_SF_bin9_cat5_upper), np.sum(id_SF_bin10_cat5_upper), np.sum(id_SF_bin11_cat5_upper)])


                                                                                                                    
NQ_cat5 = np.array([np.sum(id_Q_bin1_cat5_lower), np.sum(id_Q_bin2_cat5_lower), np.sum(id_Q_bin3_cat5_lower), np.sum(id_Q_bin4_cat5_lower),
        np.sum(id_Q_bin5_cat5_lower), np.sum(id_Q_bin6_cat5_lower), np.sum(id_Q_bin7_cat5_lower), np.sum(id_Q_bin8_cat5_lower), 
        np.sum(id_Q_bin9_cat5_lower), np.sum(id_Q_bin10_cat5_lower),np.sum(id_Q_bin11_cat5_lower)]) + np.array([np.sum(id_Q_bin1_cat5_mid), np.sum(id_Q_bin2_cat5_mid), np.sum(id_Q_bin3_cat5_mid), np.sum(id_Q_bin4_cat5_mid),
                np.sum(id_Q_bin5_cat5_mid), np.sum(id_Q_bin6_cat5_mid), np.sum(id_Q_bin7_cat5_mid), np.sum(id_Q_bin8_cat5_mid), 
                np.sum(id_Q_bin9_cat5_mid), np.sum(id_Q_bin10_cat5_mid),np.sum(id_Q_bin11_cat5_mid)]) + np.array([np.sum(id_Q_bin1_cat5_upper), np.sum(id_Q_bin2_cat5_upper), np.sum(id_Q_bin3_cat5_upper), np.sum(id_Q_bin4_cat5_upper),
                        np.sum(id_Q_bin5_cat5_upper), np.sum(id_Q_bin6_cat5_upper), np.sum(id_Q_bin7_cat5_upper), np.sum(id_Q_bin8_cat5_upper), 
                        np.sum(id_Q_bin9_cat5_upper), np.sum(id_Q_bin10_cat5_upper),np.sum(id_Q_bin11_cat5_upper)])


total_cat5 = nSF_cat5 + NQ_cat5

nSF_error5 = (nSF_cat5/total_cat5)*np.sqrt(1/nSF_cat5 + 1/total_cat5)
NQ_error5 = (NQ_cat5/total_cat5)*np.sqrt(1/NQ_cat5 + 1/total_cat5)


mask5 = nSF_cat5 > 30

#initialize layout of plot
fig1, ax = plt.subplots()

#add scatterplot for cat1
SF = ax.errorbar(mass_axe[mask5],nSF_cat5[mask5]/total_cat5[mask5], yerr = nSF_error5[mask5], fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "b", markersize = 2, label='nSF')
nQ = ax.errorbar(mass_axe[mask5],NQ_cat5[mask5]/total_cat5[mask5], yerr = NQ_error5[mask5], fmt = "o", alpha = 0.7,
            elinewidth = 0.7, capsize = 5, markeredgecolor = "r", markersize = 2, label='nQ')
ax.legend([nSF,nQ],['nSF','nQ'],loc="upper right")
ax.set_xlim(8,13)
ax.set_ylim(0,1)
ax.set_xlabel("$log_{10}(M_{\odot})$")
ax.set_ylabel("Fraction of galaxies")
ax.set_title("Fraction of galaxies with mass at 2.5 < z < 4")
plt.show()

#####################################################
# plot(z, log10(mass))
mass = data['mass'][0,:]

     
########
# 0 < z < 0.5
mean_mass_cat1 = mass[idUVJ_cat1[0,:]]
mean_mass_nSF_cat1 = np.average(mean_mass_cat1[~IDS_cat1])
mean_mass_nQ_cat1 = np.average(mean_mass_cat1[IDS_cat1])
nSF_error1 = (mean_mass_nSF_cat1/cat1_total)*np.sqrt(1/mean_mass_nSF_cat1 + 1/cat1_total)
NQ_error1 = (mean_mass_nQ_cat1/cat1_total)*np.sqrt(1/mean_mass_nQ_cat1 + 1/cat1_total)

# 0.25 < z < 0.75
mean_mass_cat2_lower = mass[idUVJ_cat2_lower[0,:]]
mean_mass_cat2_upper = mass[idUVJ_cat2_upper[0,:]]

mean_mass_nSF_cat2 = (np.average(mean_mass_cat2_lower[~IDS_cat2_lower])+np.average(mean_mass_cat2_upper[~IDS_cat2_upper]))/2
mean_mass_nQ_cat2 = (np.average(mean_mass_cat2_lower[IDS_cat2_lower])+np.average(mean_mass_cat2_upper[IDS_cat2_upper]))/2
nSF_error2 = (mean_mass_nSF_cat2/cat2_total)*np.sqrt(1/mean_mass_nSF_cat2 + 1/cat2_total)
NQ_error2 = (mean_mass_nQ_cat2/cat2_total)*np.sqrt(1/mean_mass_nQ_cat2 + 1/cat2_total)

# 0.5 < z < 1
mean_mass_cat3 = mass[idUVJ_cat3[0,:]]
mean_mass_nSF_cat3 = np.average(mean_mass_cat3[~IDS_cat3])
mean_mass_nQ_cat3 = np.average(mean_mass_cat3[IDS_cat3])
nSF_error3 = (mean_mass_nSF_cat3/cat3_total)*np.sqrt(1/mean_mass_nSF_cat3 + 1/cat3_total)
NQ_error3 = (mean_mass_nQ_cat3/cat3_total)*np.sqrt(1/mean_mass_nQ_cat3 + 1/cat3_total)

# 0.75 < z < 1.25
mean_mass_cat4_lower = mass[idUVJ_cat4_lower[0,:]]
mean_mass_cat4_upper = mass[idUVJ_cat4_upper[0,:]]

mean_mass_nSF_cat4 = (np.average(mean_mass_cat4_lower[~IDS_cat4_lower])+np.average(mean_mass_cat4_upper[~IDS_cat4_upper]))/2
mean_mass_nQ_cat4 = (np.average(mean_mass_cat4_lower[IDS_cat4_lower])+np.average(mean_mass_cat4_upper[IDS_cat4_upper]))/2
nSF_error4 = (mean_mass_nSF_cat4/cat4_total)*np.sqrt(1/mean_mass_nSF_cat4 + 1/cat4_total)
NQ_error4 = (mean_mass_nQ_cat4/cat4_total)*np.sqrt(1/mean_mass_nQ_cat4 + 1/cat4_total)

# 1 < z < 1.5
mean_mass_cat5 = mass[idUVJ_cat5[0,:]]
mean_mass_nSF_cat5 = np.average(mean_mass_cat5[~IDS_cat5])
mean_mass_nQ_cat5 = np.average(mean_mass_cat5[IDS_cat5])
nSF_error5 = (mean_mass_nSF_cat5/cat5_total)*np.sqrt(1/mean_mass_nSF_cat5 + 1/cat5_total)
NQ_error5 = (mean_mass_nQ_cat5/cat5_total)*np.sqrt(1/mean_mass_nQ_cat5 + 1/cat5_total)

# 1.25 < z < 1.75
mean_mass_cat6 = mass[idUVJ_cat6[0,:]]
mean_mass_nSF_cat6 = np.average(mean_mass_cat6[~IDS_cat6])
mean_mass_nQ_cat6 = np.average(mean_mass_cat6[IDS_cat6])
nSF_error6 = (mean_mass_nSF_cat6/cat6_total)*np.sqrt(1/mean_mass_nSF_cat6 + 1/cat6_total)
NQ_error6 = (mean_mass_nQ_cat6/cat6_total)*np.sqrt(1/mean_mass_nQ_cat6 + 1/cat6_total)

# 1.5 < z < 2
mean_mass_cat7 = mass[idUVJ_cat7[0,:]]
mean_mass_nSF_cat7 = np.average(mean_mass_cat7[~IDS_cat7])
mean_mass_nQ_cat7 = np.average(mean_mass_cat7[IDS_cat7])
nSF_error7 = (mean_mass_nSF_cat7/cat7_total)*np.sqrt(1/mean_mass_nSF_cat7 + 1/cat7_total)
NQ_error7 = (mean_mass_nQ_cat7/cat7_total)*np.sqrt(1/mean_mass_nQ_cat7 + 1/cat7_total)

# 1.75 < z < 2.25
mean_mass_cat8_lower = mass[idUVJ_cat8_lower[0,:]]
mean_mass_cat8_upper = mass[idUVJ_cat8_upper[0,:]]

mean_mass_nSF_cat8 = (np.average(mean_mass_cat8_lower[~IDS_cat8_lower])+np.average(mean_mass_cat8_upper[~IDS_cat8_upper]))/2
mean_mass_nQ_cat8 = (np.average(mean_mass_cat8_lower[IDS_cat8_lower])+np.average(mean_mass_cat8_upper[IDS_cat8_upper]))/2
nSF_error8 = (mean_mass_nSF_cat8/cat8_total)*np.sqrt(1/mean_mass_nSF_cat8 + 1/cat8_total)
NQ_error8 = (mean_mass_nQ_cat8/cat8_total)*np.sqrt(1/mean_mass_nQ_cat8 + 1/cat8_total)

# 2 < z < 2.5
mean_mass_cat9 = mass[idUVJ_cat9[0,:]]
mean_mass_nSF_cat9 = np.average(mean_mass_cat9[~IDS_cat9])
mean_mass_nQ_cat9 = np.average(mean_mass_cat9[IDS_cat9])
nSF_error9 = (mean_mass_nSF_cat9/cat9_total)*np.sqrt(1/mean_mass_nSF_cat9 + 1/cat9_total)
NQ_error9 = (mean_mass_nQ_cat9/cat9_total)*np.sqrt(1/mean_mass_nQ_cat9 + 1/cat9_total)

# 2.25 < z < 2.75
mean_mass_cat10 = mass[idUVJ_cat10[0,:]]
mean_mass_nSF_cat10 = np.average(mean_mass_cat10[~IDS_cat10])
mean_mass_nQ_cat10 = np.average(mean_mass_cat10[IDS_cat10])
nSF_error10 = (mean_mass_nSF_cat10/cat10_total)*np.sqrt(1/mean_mass_nSF_cat10 + 1/cat10_total)
NQ_error10 = (mean_mass_nQ_cat10/cat10_total)*np.sqrt(1/mean_mass_nQ_cat10 + 1/cat10_total)

# 2.5 < z < 3
mean_mass_cat11 = mass[idUVJ_cat11[0,:]]
mean_mass_nSF_cat11 = np.average(mean_mass_cat11[~IDS_cat11])
mean_mass_nQ_cat11 = np.average(mean_mass_cat11[IDS_cat11])
nSF_error11 = (mean_mass_nSF_cat11/cat11_total)*np.sqrt(1/mean_mass_nSF_cat11 + 1/cat11_total)
NQ_error11 = (mean_mass_nQ_cat11/cat11_total)*np.sqrt(1/mean_mass_nQ_cat11 + 1/cat11_total)

# 2.75 < z < 3.25
mean_mass_cat12 = mass[idUVJ_cat12[0,:]]
mean_mass_nSF_cat12 = np.average(mean_mass_cat12[~IDS_cat12])
mean_mass_nQ_cat12 = np.average(mean_mass_cat12[IDS_cat12])
nSF_error12 = (mean_mass_nSF_cat12/cat12_total)*np.sqrt(1/mean_mass_nSF_cat12 + 1/cat12_total)
NQ_error12 = (mean_mass_nQ_cat12/cat12_total)*np.sqrt(1/mean_mass_nQ_cat12 + 1/cat12_total)

# 3 < z < 3.5
mean_mass_cat13 = mass[idUVJ_cat13[0,:]]
mean_mass_nSF_cat13 = np.average(mean_mass_cat13[~IDS_cat13])
mean_mass_nQ_cat13 = np.average(mean_mass_cat13[IDS_cat13])
nSF_error13 = (mean_mass_nSF_cat13/cat13_total)*np.sqrt(1/mean_mass_nSF_cat13 + 1/cat13_total)
NQ_error13 = (mean_mass_nQ_cat13/cat13_total)*np.sqrt(1/mean_mass_nQ_cat13 + 1/cat13_total)

# 3.25 < z < 3.75
mean_mass_cat14 = mass[idUVJ_cat14[0,:]]
mean_mass_nSF_cat14 = np.average(mean_mass_cat14[~IDS_cat14])
mean_mass_nQ_cat14 = np.average(mean_mass_cat14[IDS_cat14])
nSF_error14 = (mean_mass_nSF_cat14/cat14_total)*np.sqrt(1/mean_mass_nSF_cat14 + 1/cat14_total)
NQ_error14 = (mean_mass_nQ_cat14/cat14_total)*np.sqrt(1/mean_mass_nQ_cat14 + 1/cat14_total)

# 3.5 < z < 4
mean_mass_cat15 = mass[idUVJ_cat15[0,:]]
mean_mass_nSF_cat15 = np.average(mean_mass_cat15[~IDS_cat15])
mean_mass_nQ_cat15 = np.average(mean_mass_cat15[IDS_cat15])
nSF_error15 = (mean_mass_nSF_cat15/cat15_total)*np.sqrt(1/mean_mass_nSF_cat15 + 1/cat15_total)
NQ_error15 = (mean_mass_nQ_cat15/cat15_total)*np.sqrt(1/mean_mass_nQ_cat15 + 1/cat15_total)

mass_QSF = np.array([[0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75], 
         [mean_mass_nQ_cat1,mean_mass_nQ_cat2,mean_mass_nQ_cat3,mean_mass_nQ_cat4,mean_mass_nQ_cat5,mean_mass_nQ_cat6,mean_mass_nQ_cat7,mean_mass_nQ_cat8,mean_mass_nQ_cat9,
          mean_mass_nQ_cat10,mean_mass_nQ_cat11,mean_mass_nQ_cat12,mean_mass_nQ_cat13,mean_mass_nQ_cat14,mean_mass_nQ_cat15], 
         [mean_mass_nSF_cat1,mean_mass_nSF_cat2,mean_mass_nSF_cat3,mean_mass_nSF_cat4,mean_mass_nSF_cat5,mean_mass_nSF_cat6,mean_mass_nSF_cat7,mean_mass_nSF_cat8,
          mean_mass_nSF_cat9,mean_mass_nSF_cat10,mean_mass_nSF_cat11,mean_mass_nSF_cat12,mean_mass_nSF_cat13,mean_mass_nSF_cat14,mean_mass_nSF_cat15],
         [nSF_error1,nSF_error2,nSF_error3,nSF_error4,nSF_error5,nSF_error6,nSF_error7,nSF_error8,
          nSF_error9,nSF_error10,nSF_error11,nSF_error12,nSF_error13,nSF_error14,nSF_error15],
         [NQ_error1,NQ_error2,NQ_error3,NQ_error4,NQ_error5,NQ_error6,NQ_error7,NQ_error8,
          NQ_error9,NQ_error10,NQ_error11,NQ_error12,NQ_error13,NQ_error14,NQ_error15]])

"""
Following segment can be deleted when trendlines code is updated for all bins
"""
# plot(z, log10(mass))
#Define redshift intervals for the SFRS
SFR = data['sfr'][0,:]
# 0 < z < 0.5
SFR_cat1 = SFR[idUVJ_cat1[0,:]]
SFR_nSF_cat1 = SFR_cat1[~IDS_cat1]
SFR_nQ_cat1 = SFR_cat1[IDS_cat1]
MASS_nSF_cat1 = m1[~IDS_cat1]
MASS_nQ_cat1 = m1[IDS_cat1]

# 0.5 < z < 1
SFR_cat2 = SFR[idUVJ_cat3[0,:]]
SFR_nSF_cat2 = SFR_cat2[~IDS_cat3]
SFR_nQ_cat2 = SFR_cat2[IDS_cat3]
MASS_nSF_cat2 = m3[~IDS_cat3]
MASS_nQ_cat2 = m3[IDS_cat3]

# 1 < z < 1.5
SFR_cat3 = SFR[idUVJ_cat5[0,:]]
SFR_nSF_cat3 = SFR_cat3[~IDS_cat5]
SFR_nQ_cat3 = SFR_cat3[IDS_cat5]
MASS_nSF_cat3 = m5[~IDS_cat5]
MASS_nQ_cat3 = m5[IDS_cat5]

# 1.5 < z < 2
SFR_cat4 = SFR[idUVJ_cat7[0,:]]
SFR_nSF_cat4 = SFR_cat4[~IDS_cat7]
SFR_nQ_cat4 = SFR_cat4[IDS_cat7]
MASS_nSF_cat4 = m7[~IDS_cat7]
MASS_nQ_cat4 = m7[IDS_cat7]

# 2 < z < 2.5
SFR_cat5 = SFR[idUVJ_cat9[0,:]]
SFR_nSF_cat5 = SFR_cat5[~IDS_cat9]
SFR_nQ_cat5 = SFR_cat5[IDS_cat9]
MASS_nSF_cat5 = m9[~IDS_cat9]
MASS_nQ_cat5 = m9[IDS_cat9]

# 2.5 < z < 3
SFR_cat6 = SFR[idUVJ_cat11[0,:]]
SFR_nSF_cat6 = SFR_cat6[~IDS_cat11]
SFR_nQ_cat6 = SFR_cat6[IDS_cat11]
MASS_nSF_cat6 = m11[~IDS_cat11]
MASS_nQ_cat6 = m11[IDS_cat11]

# 3 < z < 3.5
SFR_cat7 = SFR[idUVJ_cat13[0,:]]
SFR_nSF_cat7 = SFR_cat7[~IDS_cat13]
SFR_nQ_cat7 = SFR_cat7[IDS_cat13]
MASS_nSF_cat7 = m13[~IDS_cat13]
MASS_nQ_cat7 = m13[IDS_cat13]

# 3.5 < z < 4
SFR_cat8 = SFR[idUVJ_cat15[0,:]]
SFR_nSF_cat8 = SFR_cat8[~IDS_cat15]
SFR_nQ_cat8 = SFR_cat8[IDS_cat15]
MASS_nSF_cat8 = m15[~IDS_cat15]
MASS_nQ_cat8 = m15[IDS_cat15]

SFR_nSF_all = [*SFR_nSF_cat1, *SFR_nSF_cat2, *SFR_nSF_cat3, *SFR_nSF_cat4, *SFR_nSF_cat5, *SFR_nSF_cat6, *SFR_nSF_cat7, *SFR_nSF_cat8]
SFR_nQ_all = [*SFR_nQ_cat1, *SFR_nQ_cat2, *SFR_nQ_cat3, *SFR_nQ_cat4, *SFR_nQ_cat5, *SFR_nQ_cat6, *SFR_nQ_cat7, *SFR_nQ_cat8]
MASS_nSF_all = [*MASS_nSF_cat1, *MASS_nSF_cat2, *MASS_nSF_cat3, *MASS_nSF_cat4, *MASS_nSF_cat5, *MASS_nSF_cat6, *MASS_nSF_cat7, *MASS_nSF_cat8]
MASS_nQ_all = [*MASS_nQ_cat1, *MASS_nQ_cat2, *MASS_nQ_cat3, *MASS_nQ_cat4, *MASS_nQ_cat5, *MASS_nQ_cat6, *MASS_nQ_cat7, *MASS_nQ_cat8]
"""
end of segment to be deleted
"""


#and the redshift vs stellar mass plot as a hexbin
#### 
mass_all = np.hstack((mean_mass_cat1,mean_mass_cat3,mean_mass_cat5,mean_mass_cat7,mean_mass_cat9,mean_mass_cat11,mean_mass_cat13,mean_mass_cat15))

z = data['z_phot'][0,:][selection_all[0,:]]
plt.clf()
fig, ax = plt.subplots()
#alpha = 1
SF_and_Q_MASS = ax.hexbin(x = z, y = mass_all, 
                          vmax = 500, gridsize=100, 
                          cmap='Greys',mincnt=2,yscale = 'log', bins = 'log')

SF = ax.errorbar(mass_QSF[0,:],mass_QSF[2,:], 
                 yerr = mass_QSF[3,:], fmt = "o", elinewidth = 1.7, 
                 capsize = 6, markeredgecolor = "b", 
                 markersize = 4, label='nSF',)

nQ = ax.errorbar(mass_QSF[0,:],mass_QSF[1,:], 
                 yerr = mass_QSF[4,:], fmt = "x", elinewidth = 1.7, 
                 capsize = 6, markeredgecolor = "r", 
                 markersize = 4, label='nQ')

cb_SF = fig.colorbar(SF_and_Q_MASS)
cb_SF.set_label("Number of galaxies", fontsize=12)
ax.set_xlim(0,4)
ax.legend([SF,nQ],['av. SF','av. Q'],loc="upper right")
ax.set_ylim(10**6,10**13)
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_xlabel("z", fontsize=14, font = "CAMBRIA", style = 'italic')
ax.set_ylabel("$M_{\odot}$", fontsize=12)

ax.set_title("Average mass at redshift z", fontsize=18)
plt.show()
#####



"""
09-05-2023
Next, we wish to output plots of the Star Formation Rate as a function of the galaxy's solar mass
It is of interest to distinguish between the star forming, denoted SF, and the Quiscent, denoted Q, which are categorized based on the boundaries in the Williams paper.
"""
#PLOTTE SFR vs SOLAR MASS

#Define redshift intervals for the SFRS
SFR = data['sfr'][0,:]
#define when to discount galaxies due to outlier SFR-values
idSFR = np.logical_and(np.log10(SFR)>-4,np.log10(SFR)<4)

# 0 < z < 0.5
SFR_cat1 = SFR[idUVJ_cat1[0,:]] #all the star formation rates values that are in 0<z<0.5


# 0.5 < z < 1
SFR_cat2 = SFR[idUVJ_cat3[0,:]]

# 1 < z < 1.5
SFR_cat3 = SFR[idUVJ_cat5[0,:]]

# 1.5 < z < 2.5
SFR_cat4_lower = SFR[idUVJ_cat7[0,:]]
SFR_cat4_upper = SFR[idUVJ_cat9[0,:]]
SFR_cat4 = np.hstack((SFR_cat4_lower,SFR_cat4_upper))


# 2.5 < z < 4
SFR_cat5_lower = SFR[idUVJ_cat11[0,:]]
SFR_cat5_mid = SFR[idUVJ_cat13[0,:]]
SFR_cat5_upper = SFR[idUVJ_cat15[0,:]]
SFR_cat5 = np.hstack((SFR_cat5_lower, SFR_cat5_mid, SFR_cat5_upper))



###
#We need intervals of approx same number of galaxies so the linear regression will be "weighted" correctly
#thus we define a function to split into segments
def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

########################################################
# #SFR BIN 1 - redshift 0-0.5 
# SFR_nSF_cat1

#make SFR averages for broad intervals of masses for SF for 0<z<0.5
z_bin1_SFR_Mass = np.vstack((np.log10(mean_mass_cat1[~IDS_cat1]),np.log10(SFR_nSF_cat1)))
#z_bin1_ave = np.average(xx,axis=1)



## for these plots we need to seperate out the outliers - ie the ones for which log10(SFR)>4 and <-4 
#otherwise our regression is gonna be skewed

#mass 8-8.5

ID_BIN1_cat1 = np.logical_and(z_bin1_SFR_Mass[0,:]>=8,z_bin1_SFR_Mass[0,:]<8.5)
INT1_sfr_cat1 = np.vstack((z_bin1_SFR_Mass[0,:][ID_BIN1_cat1],z_bin1_SFR_Mass[1,:][ID_BIN1_cat1]))

#mass 8.5-9
ID_BIN2_cat1 = np.logical_and(z_bin1_SFR_Mass[0,:]>=8.5,z_bin1_SFR_Mass[0,:]<9)
INT2_sfr_cat1 = np.vstack((z_bin1_SFR_Mass[0,:][ID_BIN2_cat1],z_bin1_SFR_Mass[1,:][ID_BIN2_cat1]))

#mass 9-9.5
ID_BIN3_cat1 = np.logical_and(z_bin1_SFR_Mass[0,:]>=9,z_bin1_SFR_Mass[0,:]<9.5)
INT3_sfr_cat1 = np.vstack((z_bin1_SFR_Mass[0,:][ID_BIN3_cat1],z_bin1_SFR_Mass[1,:][ID_BIN3_cat1]))

#mass 9.5-10.5
ID_BIN4_cat1 = np.logical_and(z_bin1_SFR_Mass[0,:]>=9.5,z_bin1_SFR_Mass[0,:]<10.5)
INT4_sfr_cat1 = np.vstack((z_bin1_SFR_Mass[0,:][ID_BIN4_cat1],z_bin1_SFR_Mass[1,:][ID_BIN4_cat1]))

#mass 10.5-12
ID_BIN5_cat1 = np.logical_and(z_bin1_SFR_Mass[0,:]>=10.5,z_bin1_SFR_Mass[0,:]<12)
INT5_sfr_cat1 = np.vstack((z_bin1_SFR_Mass[0,:][ID_BIN5_cat1],z_bin1_SFR_Mass[1,:][ID_BIN5_cat1]))

INT_mass_cat1 = np.array([np.mean(INT1_sfr_cat1[0,:]),np.mean(INT2_sfr_cat1[0,:]),np.mean(INT3_sfr_cat1[0,:]),np.mean(INT4_sfr_cat1[0,:]),np.mean(INT5_sfr_cat1[0,:])])
INT_SFR_cat1 = np.array([np.mean(INT1_sfr_cat1[1,:]),np.mean(INT2_sfr_cat1[1,:]),np.mean(INT3_sfr_cat1[1,:]),np.mean(INT4_sfr_cat1[1,:]),np.mean(INT5_sfr_cat1[1,:])])

#calculate equation for trendline
zINT1 = np.polyfit(INT_mass_cat1, INT_SFR_cat1, 1)
pINT1 = np.poly1d(zINT1)






################################################################
#SFR BIN 2 - redshift 0.5-1 

SFR_SF_cat2 = SFR_cat2[~IDS_cat3]
SFR_Q_cat2 = SFR_cat2[IDS_cat3]
mass_SF_cat2 = np.log10(mean_mass_cat3[~IDS_cat3])
mass_Q_cat2 = np.log10(mean_mass_cat3[IDS_cat3])

#combine needed SF data into one variable
SFR_mass_stacked_cat2 = np.transpose(np.vstack((np.log10(SFR_SF_cat2), mass_SF_cat2))) #mass_sf_cat2 is already logged

#sort rows by mass
ind = np.lexsort((SFR_mass_stacked_cat2[:,1], SFR_mass_stacked_cat2[:,0]))  #so it will go from least to largest mass
SFR_mass_stacked_cat2 = SFR_mass_stacked_cat2[ind] #sorts the entire 2 column matrix

#use the split command to make 10 approx equally long arrays with galaxies
SFR_segments_cat2 = list(split(SFR_mass_stacked_cat2[:,0],10))
Mass_segments_cat2 = list(split(SFR_mass_stacked_cat2[:,1],10))


#compute average values for each segment
SFR_ave_cat2 = np.array([
    np.average(SFR_segments_cat2[0]),
    np.average(SFR_segments_cat2[1]),
    np.average(SFR_segments_cat2[2]),
    np.average(SFR_segments_cat2[3]),
    np.average(SFR_segments_cat2[4]),
    np.average(SFR_segments_cat2[5]),
    np.average(SFR_segments_cat2[6]),
    np.average(SFR_segments_cat2[7]),
    np.average(SFR_segments_cat2[8]),
    np.average(SFR_segments_cat2[9])
    ])

Mass_ave_cat2 = np.array([
    np.average(Mass_segments_cat2[0]),
    np.average(Mass_segments_cat2[1]),
    np.average(Mass_segments_cat2[2]),
    np.average(Mass_segments_cat2[3]),
    np.average(Mass_segments_cat2[4]),
    np.average(Mass_segments_cat2[5]),
    np.average(Mass_segments_cat2[6]),
    np.average(Mass_segments_cat2[7]),
    np.average(Mass_segments_cat2[8]),
    np.average(Mass_segments_cat2[9])
    ])

#calculate equation for trendline
zINT2 = np.polyfit(Mass_ave_cat2, SFR_ave_cat2, 1) #1 means poly of 1st order ie linear
pINT2 = np.poly1d(zINT2)

################################################################
#SFR BIN 3 - redshift 1-1.5 
# SFR_nSF_cat3
z_bin3_SFR_Mass = np.vstack((np.log10(mean_mass_cat5[~IDS_cat5]),np.log10(SFR_nSF_cat3)))

#mass 8-8.5
ID_BIN1_cat3 = np.logical_and(z_bin3_SFR_Mass[0,:]>=8,z_bin3_SFR_Mass[0,:]<8.5)
INT1_sfr_cat3 = np.vstack((z_bin3_SFR_Mass[0,:][ID_BIN1_cat3],z_bin3_SFR_Mass[1,:][ID_BIN1_cat3]))

#mass 8.5-9
ID_BIN2_cat3 = np.logical_and(z_bin3_SFR_Mass[0,:]>=8.5,z_bin3_SFR_Mass[0,:]<9)
INT2_sfr_cat3 = np.vstack((z_bin3_SFR_Mass[0,:][ID_BIN2_cat3],z_bin3_SFR_Mass[1,:][ID_BIN2_cat3]))

#mass 9-9.5
ID_BIN3_cat3 = np.logical_and(z_bin3_SFR_Mass[0,:]>=9,z_bin3_SFR_Mass[0,:]<9.5)
INT3_sfr_cat3 = np.vstack((z_bin3_SFR_Mass[0,:][ID_BIN3_cat3],z_bin3_SFR_Mass[1,:][ID_BIN3_cat3]))

#mass 9.5-10.5
ID_BIN4_cat3 = np.logical_and(z_bin3_SFR_Mass[0,:]>=9.5,z_bin3_SFR_Mass[0,:]<10.5)
INT4_sfr_cat3 = np.vstack((z_bin3_SFR_Mass[0,:][ID_BIN4_cat3],z_bin3_SFR_Mass[1,:][ID_BIN4_cat3]))

#mass 10.5-12
ID_BIN5_cat3 = np.logical_and(z_bin3_SFR_Mass[0,:]>=10.5,z_bin3_SFR_Mass[0,:]<12)
INT5_sfr_cat3 = np.vstack((z_bin3_SFR_Mass[0,:][ID_BIN5_cat3],z_bin3_SFR_Mass[1,:][ID_BIN5_cat3]))

INT_mass_cat3 = np.array([np.mean(INT1_sfr_cat3[0,:]),np.mean(INT2_sfr_cat3[0,:]),np.mean(INT3_sfr_cat3[0,:]),np.mean(INT4_sfr_cat3[0,:]),np.mean(INT5_sfr_cat3[0,:])])
INT_SFR_cat3 = np.array([np.mean(INT1_sfr_cat3[1,:]),np.mean(INT2_sfr_cat3[1,:]),np.mean(INT3_sfr_cat3[1,:]),np.mean(INT4_sfr_cat3[1,:]),np.mean(INT5_sfr_cat3[1,:])])


#calculate equation for trendline
zINT3 = np.polyfit(INT_mass_cat3, INT_SFR_cat3, 1)
pINT3 = np.poly1d(zINT3)

################################################################
#SFR BIN 4 - redshift 1.5-2.5 
# SFR_nSF_cat4
mass_of_z1525 =  [*np.log10(mean_mass_cat7[~IDS_cat7]),*np.log10(mean_mass_cat9[~IDS_cat9])]
sfr_of_z1525 = [*np.log10(SFR_nSF_cat4),*np.log10(SFR_nSF_cat5)]
z_bin4_SFR_Mass = np.vstack((mass_of_z1525,sfr_of_z1525))

#mass 8-8.5
ID_BIN1_cat4 = np.logical_and(z_bin4_SFR_Mass[0,:]>=8,z_bin4_SFR_Mass[0,:]<8.5)
INT1_sfr_cat4 = np.vstack((z_bin4_SFR_Mass[0,:][ID_BIN1_cat4],z_bin4_SFR_Mass[1,:][ID_BIN1_cat4]))

#mass 8.5-9
ID_BIN2_cat4 = np.logical_and(z_bin4_SFR_Mass[0,:]>=8.5,z_bin4_SFR_Mass[0,:]<9)
INT2_sfr_cat4 = np.vstack((z_bin4_SFR_Mass[0,:][ID_BIN2_cat4],z_bin4_SFR_Mass[1,:][ID_BIN2_cat4]))

#mass 9-9.5
ID_BIN3_cat4 = np.logical_and(z_bin4_SFR_Mass[0,:]>=9,z_bin4_SFR_Mass[0,:]<9.5)
INT3_sfr_cat4 = np.vstack((z_bin4_SFR_Mass[0,:][ID_BIN3_cat4],z_bin4_SFR_Mass[1,:][ID_BIN3_cat4]))

#mass 9.5-10.5
ID_BIN4_cat4 = np.logical_and(z_bin4_SFR_Mass[0,:]>=9.5,z_bin4_SFR_Mass[0,:]<10.5)
INT4_sfr_cat4 = np.vstack((z_bin4_SFR_Mass[0,:][ID_BIN4_cat4],z_bin4_SFR_Mass[1,:][ID_BIN4_cat4]))

#mass 10.5-12
ID_BIN5_cat4 = np.logical_and(z_bin4_SFR_Mass[0,:]>=10.5,z_bin4_SFR_Mass[0,:]<12)
INT5_sfr_cat4 = np.vstack((z_bin4_SFR_Mass[0,:][ID_BIN5_cat4],z_bin4_SFR_Mass[1,:][ID_BIN5_cat4]))

INT_mass_cat4 = np.array([np.mean(INT1_sfr_cat4[0,:]),np.mean(INT2_sfr_cat4[0,:]),np.mean(INT3_sfr_cat4[0,:]),np.mean(INT4_sfr_cat4[0,:]),np.mean(INT5_sfr_cat4[0,:])])
INT_SFR_cat4 = np.array([np.mean(INT1_sfr_cat4[1,:]),np.mean(INT2_sfr_cat4[1,:]),np.mean(INT3_sfr_cat4[1,:]),np.mean(INT4_sfr_cat4[1,:]),np.mean(INT5_sfr_cat4[1,:])])
#calculate equation for trendline
zINT4 = np.polyfit(INT_mass_cat4, INT_SFR_cat4, 1)
pINT4 = np.poly1d(zINT4)

################################################################
#SFR BIN 5 - redshift 2.5 - 4 

mass_of_z2540 =  [*np.log10(mean_mass_cat11[~IDS_cat11]),*np.log10(mean_mass_cat13[~IDS_cat13]),*np.log10(mean_mass_cat15[~IDS_cat15])]
sfr_of_z2540 = [*np.log10(SFR_nSF_cat6),*np.log10(SFR_nSF_cat7),*np.log10(SFR_nSF_cat8)]
z_bin5_SFR_Mass = np.vstack((mass_of_z2540,sfr_of_z2540))

#mass 8-8.5
ID_BIN1_cat5 = np.logical_and(z_bin5_SFR_Mass[0,:]>=8,z_bin5_SFR_Mass[0,:]<8.5)
INT1_sfr_cat5 = np.vstack((z_bin5_SFR_Mass[0,:][ID_BIN1_cat5],z_bin5_SFR_Mass[1,:][ID_BIN1_cat5]))

#mass 8.5-9
ID_BIN2_cat5 = np.logical_and(z_bin5_SFR_Mass[0,:]>=8.5,z_bin5_SFR_Mass[0,:]<9)
INT2_sfr_cat5 = np.vstack((z_bin5_SFR_Mass[0,:][ID_BIN2_cat5],z_bin5_SFR_Mass[1,:][ID_BIN2_cat5]))

#mass 9-9.5
ID_BIN3_cat5 = np.logical_and(z_bin5_SFR_Mass[0,:]>=9,z_bin5_SFR_Mass[0,:]<9.5)
INT3_sfr_cat5 = np.vstack((z_bin5_SFR_Mass[0,:][ID_BIN3_cat5],z_bin5_SFR_Mass[1,:][ID_BIN3_cat5]))

#mass 9.5-10.5
ID_BIN4_cat5 = np.logical_and(z_bin5_SFR_Mass[0,:]>=9.5,z_bin5_SFR_Mass[0,:]<10.5)
INT4_sfr_cat5 = np.vstack((z_bin5_SFR_Mass[0,:][ID_BIN4_cat5],z_bin5_SFR_Mass[1,:][ID_BIN4_cat5]))

#mass 10.5-12
ID_BIN5_cat5 = np.logical_and(z_bin5_SFR_Mass[0,:]>=10.5,z_bin5_SFR_Mass[0,:]<12)
INT5_sfr_cat5 = np.vstack((z_bin5_SFR_Mass[0,:][ID_BIN5_cat5],z_bin5_SFR_Mass[1,:][ID_BIN5_cat5]))

INT_mass_cat5 = np.array([np.mean(INT1_sfr_cat5[0,:]),np.mean(INT2_sfr_cat5[0,:]),np.mean(INT3_sfr_cat5[0,:]),np.mean(INT4_sfr_cat5[0,:]),np.mean(INT5_sfr_cat5[0,:])])
INT_SFR_cat5 = np.array([np.mean(INT1_sfr_cat5[1,:]),np.mean(INT2_sfr_cat5[1,:]),np.mean(INT3_sfr_cat5[1,:]),np.mean(INT4_sfr_cat5[1,:]),np.mean(INT5_sfr_cat5[1,:])])
#calculate equation for trendline
zINT5 = np.polyfit(INT_mass_cat5, INT_SFR_cat5, 1)
pINT5 = np.poly1d(zINT5)



#######

#commence plot
plt.clf()
fig, ax1 = plt.subplots()

mass_array = range(10**7, 10**12, 10000000)
logmass = np.log10(mass_array)
lookback1 = 10.754
logSF1 = (0.84 - 0.026*lookback1)*np.log10(mass_array) - (6.51 - 0.11*lookback1)

#all the data points for the redshift bin
SFR_SF_cat1_p = ax1.hexbin(np.log10(mean_mass_cat1[~IDS_cat1]), np.log10(SFR_nSF_cat1),vmax = 1, cmap = "Blues", mincnt = 1, gridsize=(346,200))
SFR_Q_cat1_p = ax1.hexbin(np.log10(mean_mass_cat1[IDS_cat1]), np.log10(SFR_nQ_cat1),vmax = 1, cmap = "Reds", mincnt = 1, gridsize=(346,200))

trend_our_cat1 = ax1.plot(INT_mass_cat1,pINT1(INT_mass_cat1),color="limegreen", linewidth=2, linestyle="-", label = 'Simple regression')
trend_paper_cat1 = ax1.plot(logmass,logSF1,c="magenta",ls="-.",lw=2,label= 'Speagle et al')

trend_points_cat1 = ax1.scatter(INT_mass_cat1,INT_SFR_cat1, c='orange', s = 100, marker = 'x')

cb1 = fig.colorbar(SFR_SF_cat1_p)
cb2 = fig.colorbar(SFR_Q_cat1_p)
cb1.set_label("Number of SF galaxies", fontsize=10)
cb2.set_label("Number of Q galaxies", fontsize=10)
ax1.set_title('SFR as related to mass for 0 < z < 0.5 \n', fontsize = 14, loc = 'left')
ax1.set_xlabel('$log_{10}$ ($M_{\odot}$)')
ax1.set_ylabel('$log_{10}$(SFR)')
ax1.set_ylim(-2.75,2) # cut v -4 fordi der kommer linjer fra SED- fitting modellerne v lavere SFR
ax1.set_xlim(8,12)

plt.legend()
plt.show()

####
#SFR BIN 2 - corresponding to redshift interval earlier called cat3 
SFR_SF_cat2 = SFR_cat2[~IDS_cat3]
SFR_Q_cat2 = SFR_cat2[IDS_cat3]
mass_SF_cat2 = np.log10(mean_mass_cat3[~IDS_cat3])
mass_Q_cat2 = np.log10(mean_mass_cat3[IDS_cat3])

#commence plot
plt.clf()
fig, ax2 = plt.subplots()

mass_array = range(10**7, 10**12, 10000000)
logmass = np.log10(mass_array)
lookback2 = 7.069
logSF2 = (0.84 - 0.026*lookback2)*np.log10(mass_array) - (6.51 - 0.11*lookback2)

SFR_SF_cat2_p = ax2.hexbin(mass_SF_cat2, np.log10(SFR_SF_cat2),vmax = 1, cmap = "Blues", mincnt = 1, gridsize=(346,200))
SFR_Q_cat2_p = ax2.hexbin(mass_Q_cat2, np.log10(SFR_Q_cat2),vmax = 1, cmap = "Reds", mincnt = 1, gridsize=(346,200))
trendpoints_our_cat2 = ax2.scatter(Mass_ave_cat2,SFR_ave_cat2, c='gold', s = 60, marker = 'x')

trend_our_cat2 = ax2.plot(Mass_ave_cat2,pINT2(Mass_ave_cat2),color="limegreen", linewidth=2, linestyle="-",label='Simple regression',zorder=100)
trend_paper_cat2 = ax2.plot(logmass,logSF2,c="magenta",ls="-.",lw=2, label= 'Speagle et al')

cb1 = fig.colorbar(SFR_SF_cat2_p)
cb2 = fig.colorbar(SFR_Q_cat2_p)

cb1.set_label("Number of SF galaxies", fontsize=10)
cb2.set_label("Number of Q galaxies", fontsize=10)
ax2.set_title('SFR as related to mass for 0.5 < z < 1 \n',fontsize = 14,loc = 'left')
ax2.set_xlabel('$log_{10}$ ($M_{\odot}$)')
ax2.set_ylabel('$log_{10}$(SFR)')
ax2.set_ylim(-2.5,4)  
ax2.set_xlim(8,12)
plt.legend()
plt.show() #Once again it is seen how theres a significant increase in the average mass observed for the galaxies when redshift increases - observational bias?
#timescale is between approx 5 and 7.8 Gy back in time
#disregard log SFR <-4
#linear regression
#expect higher cross section value  as z increase
#bin logmass 0.3 ish - show trend. fit ave
#compare w relation from paper

####
#SFR BIN 3 - redshift 1-1.5
SFR_SF_cat3 = SFR_cat3[~IDS_cat5]
SFR_Q_cat3 = SFR_cat3[IDS_cat5]
mass_SF_cat3 = np.log10(mean_mass_cat5[~IDS_cat5])
mass_Q_cat3 = np.log10(mean_mass_cat5[IDS_cat5])

#commence plot
plt.clf() #we clear the plot so matplotlib doesnt accidentlly plot on top of previous
fig, ax3 = plt.subplots()

mass_array = range(10**7, 10**12, 10000000)
logmass = np.log10(mass_array)
lookback3 = 5.011
logSF3 = (0.84 - 0.026*lookback3)*np.log10(mass_array) - (6.51 - 0.11*lookback3)

SFR_SF_cat3_p = ax3.hexbin(mass_SF_cat3, np.log10(SFR_SF_cat3), vmax = 1, cmap = "Blues", mincnt = 1, gridsize=(346,200))
SFR_Q_cat3_p = ax3.hexbin(mass_Q_cat3, np.log10(SFR_Q_cat3), vmax = 1, cmap = "Reds", mincnt = 1, gridsize=(346,200))

trend_our_cat3 = ax3.plot(INT_mass_cat3,pINT3(INT_mass_cat3),color="limegreen", linewidth=2, linestyle="-", label = 'Simple regression')
trend_paper_cat3 = ax3.plot(logmass,logSF3,c="magenta",ls="-.",lw=2, label = 'Speagle et al')
trend_points_cat3 = ax3.scatter(INT_mass_cat3,INT_SFR_cat3, c='orange', s = 100, marker = 'x')

cb1 = fig.colorbar(SFR_SF_cat3_p)
cb2 = fig.colorbar(SFR_Q_cat3_p)
cb1.set_label("Number of SF galaxies", fontsize=10)
cb2.set_label("Number of Q galaxies", fontsize=10)
ax3.set_title('SFR as related to mass for 1 < z < 1.5 \n',fontsize = 14, loc = 'left')
ax3.set_xlabel('$log_{10}$ ($M_{\odot}$)')
ax3.set_ylabel('$log_{10}$(SFR)')
ax3.set_ylim(-2.25,4)
ax3.set_xlim(8,12)

plt.legend()
plt.show()  #only every 10th point is plotted bc otherwise it is difficult to discern tendency

####
#SFR BIN 4 - redshift 1.5-2.5
newIDS_cat4 = np.hstack((IDS_cat7,IDS_cat9)) #combining 2 redshift bins into 1 
SFR_SF_cat4 = SFR_cat4[~newIDS_cat4]
SFR_Q_cat4 = SFR_cat4[newIDS_cat4]
mass_SF_cat4 = np.log10(np.hstack((mean_mass_cat7[~IDS_cat7],mean_mass_cat9[~IDS_cat9])))
mass_Q_cat4 = np.log10(np.hstack((mean_mass_cat7[IDS_cat7],mean_mass_cat9[IDS_cat9])))

#commence plot
plt.clf() #we clear the plot so matplotlib doesnt accidentlly plot on top of prev
fig, ax4 = plt.subplots() #there's no pythonic reason for calling the figure ax"number" - its solely for readability

mass_array = range(10**7, 10**12, 10000000)
logmass = np.log10(mass_array)
lookback4 = 3.316
logSF4 = (0.84 - 0.026*lookback4)*np.log10(mass_array) - (6.51 - 0.11*lookback4)

SFR_SF_cat4_p = ax4.hexbin(mass_SF_cat4, np.log10(SFR_SF_cat4), vmax = 1, cmap = "Blues", mincnt = 1, gridsize=(346,200)) #only every 10 point is plotted
SFR_Q_cat4_p = ax4.hexbin(mass_Q_cat4, np.log10(SFR_Q_cat4), vmax = 1, cmap = "Reds", mincnt = 1, gridsize=(346,200))

trend_our_cat4 = ax4.plot(INT_mass_cat4,pINT4(INT_mass_cat4),color="limegreen", linewidth=2, linestyle="-", label = 'Simple regression')
trend_paper_cat4 = ax4.plot(logmass,logSF4,c="magenta",ls="-.",lw=2, label = 'Speagle et al')
trend_points_cat4 = ax4.scatter(INT_mass_cat4,INT_SFR_cat4, c='orange', s = 100, marker = 'x')

cb1 = fig.colorbar(SFR_SF_cat4_p)
cb2 = fig.colorbar(SFR_Q_cat4_p)
cb1.set_label("Number of SF galaxies", fontsize=10)
cb2.set_label("Number of Q galaxies", fontsize=10)
ax4.set_title('SFR as related to mass for 1.5 < z < 2.5 \n', loc = 'left')
ax4.set_xlabel('$log_{10}$ ($M_{\odot}$)')
ax4.set_ylabel('$log_{10}$(SFR)')
ax4.set_ylim(-2,4)
ax4.set_xlim(8,12)

plt.legend()
plt.show()

####
#SFR B IN 5 - redshift 2.5-4
newIDS_cat5 = np.hstack((IDS_cat11,IDS_cat13,IDS_cat15)) #combining 2 redshift bins into 1 
SFR_SF_cat5 = SFR_cat5[~newIDS_cat5]
SFR_Q_cat5 = SFR_cat5[newIDS_cat5]
mass_SF_cat5 = np.log10(np.hstack((mean_mass_cat11[~IDS_cat11],mean_mass_cat13[~IDS_cat13],mean_mass_cat15[~IDS_cat15])))
mass_Q_cat5 = np.log10(np.hstack((mean_mass_cat11[IDS_cat11],mean_mass_cat13[IDS_cat13],mean_mass_cat15[IDS_cat15])))

#commence plot
plt.clf() #we clear the plot so matplotlib doesnt accidentlly plot on top of prev
fig, ax5 = plt.subplots() #there's no pythonic reason for calling the figure ax"number" - its solely for readability

mass_array = range(10**7, 10**12, 10000000)
logmass = np.log10(mass_array)
lookback5 = 1.985
logSF5 = (0.84 - 0.026*lookback5)*np.log10(mass_array) - (6.51 - 0.11*lookback5)

SFR_SF_cat5_p = ax5.hexbin(mass_SF_cat5, np.log10(SFR_SF_cat5), vmax = 1, cmap = "Blues", mincnt = 1, gridsize=(346,200)) #only every 10 point is plotted
SFR_Q_cat5_p = ax5.hexbin(mass_Q_cat5, np.log10(SFR_Q_cat5), vmax = 1, cmap = "Reds", mincnt = 1, gridsize=(346,200))

trend_our_cat5 = ax5.plot(INT_mass_cat5,pINT5(INT_mass_cat5),color="limegreen", linewidth=2, linestyle="-", label = 'Simple regression')
trend_paper_cat5 = ax5.plot(logmass,logSF5,c="magenta",ls="-.",lw=2, label = 'Speagle et al')
trend_points_cat5 = ax5.scatter(INT_mass_cat5,INT_SFR_cat5, c='orange', s = 100, marker = 'x')

cb1 = fig.colorbar(SFR_SF_cat5_p)
cb2 = fig.colorbar(SFR_Q_cat5_p)
cb1.set_label("Number of SF galaxies", fontsize=10)
cb2.set_label("Number of Q galaxies", fontsize=10)
ax5.set_title('SFR as related to mass for 2.5 < z < 4 \n', fontsize = 14, loc = 'left')
ax5.set_xlabel('$log_{10}$ ($M_{\odot}$)')
ax5.set_ylabel('$log_{10}$(SFR)')
ax5.set_ylim(-1.75,4)
ax5.set_xlim(8,12)

plt.legend()
plt.show()


#making the code of shaded areas

dy = 0.2

fig, ax = plt.subplots()

speagle1 = ax.plot(logmass,logSF1,c="red",ls="-",lw=2,label= 'z = 0.25')
shade1 = ax.fill_between(logmass, logSF1 - dy, logSF1 + dy, color = "red", alpha = 0.3)
speagle2 = ax.plot(logmass,logSF2,c="green",ls="-",lw=2, label= 'z = 0.75')
shade2 = ax.fill_between(logmass, logSF2 - dy, logSF2 + dy, color = "green", alpha = 0.3)
speagle3 = ax.plot(logmass,logSF3,c="blue",ls="-",lw=2,label= 'z = 1.25')
shade3 = ax.fill_between(logmass, logSF3 - dy, logSF3 + dy, color = "blue", alpha = 0.3)
speagle4 = ax.plot(logmass,logSF4,c="brown",ls="-",lw=2, label= 'z = 2')
shade4 = ax.fill_between(logmass, logSF4 - dy, logSF4 + dy, color = "brown", alpha = 0.3)
speagle5 = ax.plot(logmass,logSF5,c="black",ls="-",lw=2,label= 'z = 3.25')
shade2 = ax.fill_between(logmass, logSF5 - dy, logSF5 + dy, color = "black", alpha = 0.3)
ax.set_xlabel('$log_{10}$ ($M_{\odot}$)')
ax.set_ylabel('$log_{10}$(SFR)')
ax.set_ylim(-1,3)
ax.set_xlim(9,11)

plt.legend()


# Close FITS file so it won't use up excess memory
hdu.close()