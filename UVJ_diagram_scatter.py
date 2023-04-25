# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 10:03:50 2023

@author: Nikolaj Lange Dons
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 00:24:57 2023

@author: Nikolaj Lange Dons
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import gridspec

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

# cat1: 0 < z < 0.5
idUVJ_cat1 = np.logical_and((data['z_phot'] > 0), (data['z_phot'] < 0.5), (data['CFHT_u_MAG'] < 50))
iduv1_cat1 = UVz[idUVJ_cat1[0,:]] >= 1.3
iduv2_cat1 = UVz[idUVJ_cat1[0,:]] <= 3
idvj1_cat1 = VJz[idUVJ_cat1[0,:]] <= 1.6
idvj2_cat1 = VJz[idUVJ_cat1[0,:]] >= -1
iduvj_cat1 = UVz[idUVJ_cat1[0,:]] >= 0.88*VJz[idUVJ_cat1[0,:]]+0.69
IDS_cat1 = iduv1_cat1 & idvj1_cat1 & iduvj_cat1 & idvj2_cat1 & iduv2_cat1
NQ_cat1 = sum(IDS_cat1)
NSF_cat1 = sum(~IDS_cat1)
cat1_total = NQ_cat1 + NSF_cat1

# cat2: 0.25 < z < 0.75
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

NQ_cat2 = sum(IDS_cat2_lower)+sum(IDS_cat2_upper)
NSF_cat2 = sum(~IDS_cat2_lower)+sum(~IDS_cat2_upper)
cat2_total = NQ_cat2 + NSF_cat2

# cat3: 0.5 < z < 1.0
idUVJ_cat3 = np.logical_and((data['z_phot'] > 0.5), (data['z_phot'] < 1), (data['CFHT_u_MAG'] < 50))
iduv1_cat3 = UVz[idUVJ_cat3[0,:]] >= 1.3
iduv2_cat3 = UVz[idUVJ_cat3[0,:]] <= 3
idvj1_cat3 = VJz[idUVJ_cat3[0,:]] <= 1.6
idvj2_cat3 = VJz[idUVJ_cat3[0,:]] >= -1
iduvj_cat3 = UVz[idUVJ_cat3[0,:]] >= 0.88*VJz[idUVJ_cat3[0,:]]+0.59
IDS_cat3 = iduv1_cat3 & idvj1_cat3 & iduvj_cat3 & idvj2_cat3 & iduv2_cat3
NQ_cat3 = sum(IDS_cat3)
NSF_cat3 = sum(~IDS_cat3)
cat3_total = NQ_cat3 + NSF_cat3

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

NQ_cat4 = sum(IDS_cat4_lower)+sum(IDS_cat4_upper)
NSF_cat4 = sum(~IDS_cat4_lower)+sum(~IDS_cat4_upper)
cat4_total = NQ_cat4 + NSF_cat4

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
idUVJ_cat7 = np.logical_and((data['z_phot'] > 1.25), (data['z_phot'] < 1.75), (data['CFHT_u_MAG'] < 50))
iduv1_cat7 = UVz[idUVJ_cat7[0,:]] >= 1.3
iduv2_cat7 = UVz[idUVJ_cat7[0,:]] <= 3
idvj1_cat7 = VJz[idUVJ_cat7[0,:]] <= 1.6
idvj2_cat7 = VJz[idUVJ_cat7[0,:]] >= -1
iduvj_cat7 = UVz[idUVJ_cat7[0,:]] >= 0.88*VJz[idUVJ_cat7[0,:]]+0.49
IDS_cat7 = iduv1_cat7 & idvj1_cat7 & iduvj_cat7 & idvj2_cat7 & iduv2_cat7
NQ_cat7 = sum(IDS_cat7)
NSF_cat7 = sum(~IDS_cat7)
cat7_total = NQ_cat7 + NSF_cat7

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
idUVJ_cat9 = np.logical_and((data['z_phot'] > 2), (data['z_phot'] < 2.25), (data['CFHT_u_MAG'] < 50))
iduv1_cat9 = UVz[idUVJ_cat9[0,:]] >= 1.3
iduv2_cat9 = UVz[idUVJ_cat9[0,:]] <= 3
idvj1_cat9 = VJz[idUVJ_cat9[0,:]] <= 1.6
idvj2_cat9 = VJz[idUVJ_cat9[0,:]] >= -1
iduvj_cat9 = UVz[idUVJ_cat9[0,:]] >= 0.88*VJz[idUVJ_cat9[0,:]]+0.39
IDS_cat9 = iduv1_cat9 & idvj1_cat9 & iduvj_cat9 & idvj2_cat9 & iduv2_cat9
NQ_cat9 = sum(IDS_cat9)
NSF_cat9 = sum(~IDS_cat9)
cat9_total = NQ_cat9 + NSF_cat9


# cat10: 2.25 < z < 2.75
idUVJ_cat10 = np.logical_and((data['z_phot'] > 2.25), (data['z_phot'] < 2.75), (data['CFHT_u_MAG'] < 50))
iduv1_cat10 = UVz[idUVJ_cat10[0,:]] >= 1.3
iduv2_cat10 = UVz[idUVJ_cat10[0,:]] <= 3
idvj1_cat10 = VJz[idUVJ_cat10[0,:]] <= 1.6
idvj2_cat10 = VJz[idUVJ_cat10[0,:]] >= -1
iduvj_cat10 = UVz[idUVJ_cat10[0,:]] >= 0.88*VJz[idUVJ_cat10[0,:]]+0.39
IDS_cat10 = iduv1_cat10 & idvj1_cat10 & iduvj_cat10 & idvj2_cat10 & iduv2_cat10
NQ_cat10 = sum(IDS_cat10)
NSF_cat10 = sum(~IDS_cat10)
cat10_total = NQ_cat10 + NSF_cat10

# cat11: 2.5 < z < 3
idUVJ_cat11 = np.logical_and((data['z_phot'] > 2.5), (data['z_phot'] < 3), (data['CFHT_u_MAG'] < 50))
iduv1_cat11 = UVz[idUVJ_cat11[0,:]] >= 1.3
iduv2_cat11 = UVz[idUVJ_cat11[0,:]] <= 3
idvj1_cat11 = VJz[idUVJ_cat11[0,:]] <= 1.6
idvj2_cat11 = VJz[idUVJ_cat11[0,:]] >= -1
iduvj_cat11 = UVz[idUVJ_cat11[0,:]] >= 0.88 * VJz[idUVJ_cat11[0,:]] + 0.39
IDS_cat11 = iduv1_cat11 & idvj1_cat11 & iduvj_cat11 & idvj2_cat11 & iduv2_cat11
NQ_cat11 = sum(IDS_cat11)
NSF_cat11 = sum(~IDS_cat11)
cat11_total = NQ_cat11 + NSF_cat11

# cat12: 2.75 < z < 3.25
idUVJ_cat12 = np.logical_and((data['z_phot'] > 2.75), (data['z_phot'] < 3.25), (data['CFHT_u_MAG'] < 50))
iduv1_cat12 = UVz[idUVJ_cat12[0,:]] >= 1.3
iduv2_cat12 = UVz[idUVJ_cat12[0,:]] <= 3
idvj1_cat12 = VJz[idUVJ_cat12[0,:]] <= 1.6
idvj2_cat12 = VJz[idUVJ_cat12[0,:]] >= -1
iduvj_cat12 = UVz[idUVJ_cat12[0,:]] >= 0.88 * VJz[idUVJ_cat12[0,:]] + 0.39
IDS_cat12 = iduv1_cat12 & idvj1_cat12 & iduvj_cat12 & idvj2_cat12 & iduv2_cat12
NQ_cat12 = sum(IDS_cat12)
NSF_cat12 = sum(~IDS_cat12)
cat12_total = NQ_cat12 + NSF_cat12

# cat13: 3 < z < 3.5
idUVJ_cat13 = np.logical_and((data['z_phot'] > 3), (data['z_phot'] < 3.5), (data['CFHT_u_MAG'] < 50))
iduv1_cat13 = UVz[idUVJ_cat13[0,:]] >= 1.3
iduv2_cat13 = UVz[idUVJ_cat13[0,:]] <= 3
idvj1_cat13 = VJz[idUVJ_cat13[0,:]] <= 1.6
idvj2_cat13 = VJz[idUVJ_cat13[0,:]] >= -1
iduvj_cat13 = UVz[idUVJ_cat13[0,:]] >= 0.88 * VJz[idUVJ_cat13[0,:]] + 0.39
IDS_cat13 = iduv1_cat13 & idvj1_cat13 & iduvj_cat13 & idvj2_cat13 & iduv2_cat13
NQ_cat13 = sum(IDS_cat13)
NSF_cat13 = sum(~IDS_cat13)
cat13_total = NQ_cat13 + NSF_cat13

# cat14: 3.25 < z < 3.75
idUVJ_cat14 = np.logical_and((data['z_phot'] > 3.25), (data['z_phot'] < 3.75), (data['CFHT_u_MAG'] < 50))
iduv1_cat14 = UVz[idUVJ_cat14[0,:]] >= 1.3
iduv2_cat14 = UVz[idUVJ_cat14[0,:]] <= 3
idvj1_cat14 = VJz[idUVJ_cat14[0,:]] <= 1.6
idvj2_cat14 = VJz[idUVJ_cat14[0,:]] >= -1
iduvj_cat14 = UVz[idUVJ_cat14[0,:]] >= 0.88 * VJz[idUVJ_cat14[0,:]] + 0.39
IDS_cat14 = iduv1_cat14 & idvj1_cat14 & iduvj_cat14 & idvj2_cat14 & iduv2_cat14
NQ_cat14 = sum(IDS_cat14)
NSF_cat14 = sum(~IDS_cat14)
cat14_total = NQ_cat14 + NSF_cat14

# cat15: 3.5 < z < 4
idUVJ_cat15 = np.logical_and((data['z_phot'] > 3.5), (data['z_phot'] < 4), (data['CFHT_u_MAG'] < 50))
iduv1_cat15 = UVz[idUVJ_cat15[0,:]] >= 1.3
iduv2_cat15 = UVz[idUVJ_cat15[0,:]] <= 3
idvj1_cat15 = VJz[idUVJ_cat15[0,:]] <= 1.6
idvj2_cat15 = VJz[idUVJ_cat15[0,:]] >= -1
iduvj_cat15 = UVz[idUVJ_cat15[0,:]] >= 0.88 * VJz[idUVJ_cat15[0,:]] + 0.39
IDS_cat15 = iduv1_cat15 & idvj1_cat15 & iduvj_cat15 & idvj2_cat15 & iduv2_cat15
NQ_cat15 = sum(IDS_cat15)
NSF_cat15 = sum(~IDS_cat15)
cat15_total = NQ_cat15 + NSF_cat15

# Create a matrix
QSF = np.array([[0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75], 
         [NQ_cat1,NQ_cat2,NQ_cat3,NQ_cat4,NQ_cat5,NQ_cat6,NQ_cat7,NQ_cat8,NQ_cat9,
          NQ_cat10,NQ_cat11,NQ_cat12,NQ_cat13,NQ_cat14,NQ_cat15], 
         [NSF_cat1,NSF_cat2,NSF_cat3,NSF_cat4,NSF_cat5,NSF_cat6,NSF_cat7,NSF_cat8,
          NSF_cat9,NSF_cat10,NSF_cat11,NSF_cat12,NSF_cat13,NSF_cat14,NSF_cat15],
         [cat1_total,cat2_total,cat3_total,cat4_total,cat5_total,cat6_total,
          cat7_total,cat8_total,cat9_total,cat10_total,cat11_total,cat12_total,
          cat13_total,cat14_total,cat15_total]])


#initialize layout of plot for nQ/nSF
fig1, ax = plt.subplots()

#add scatterplot for cat1
SF = ax.plot(QSF[0,:],QSF[2,:]/QSF[3,:], 'b', label="nSF")
Q = ax.plot(QSF[0,:],QSF[1,:]/QSF[3,:], 'r', label="nQ")
plt.legend()

ax.set_xlim(0,4)
ax.set_xlabel("Redshift z")
ax.set_ylabel("number of galaxies")
ax.set_title("Number of galaxies at redshift z")

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
ax[0,0].set_xlim(0,2.5)
ax[0,0].set_ylim(0,2.5)
ax[0,0].set_ylabel("UV", fontsize=12)
ax[0,0].text(0.25,0.925,"0 < z < 0.5",ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)

# cat2
hb2 = ax[0,1].hexbin(VJcat2, UVcat2, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[0,1].set_xlim(0,2.5)
ax[0,1].set_ylim(0,2.5)
ax[0,1].text(0.35,0.925,"0.25 < z <0.75",ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)

# cat3
hb3 = ax[1,0].hexbin(VJcat3, UVcat3, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,0].set_xlim(0,2.5)
ax[1,0].set_ylim(0,2.5)
ax[1,0].set_xlabel("VJ", fontsize=12)
ax[1,0].set_ylabel("UV", fontsize=12)
ax[1,0].text(0.25,0.925,"0.5 < z < 1",ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)

# cat4
hb4 = ax[1,1].hexbin(VJcat4, UVcat4, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,1].set_xlim(0,2.5)
ax[1,1].set_ylim(0,2.5)
ax[1,1].set_xlabel("VJ", fontsize=12)
ax[1,1].text(0.35,0.925,"0.75 < z < 1.25",ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)

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
ax[0,0].set_xlim(0,2.5)
ax[0,0].set_ylim(0,2.5)
ax[0,0].set_ylabel("UV", fontsize=12)
ax[0,0].text(0.25,0.925,"1 < z < 1.5",ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)

# cat6
hb6 = ax[0,1].hexbin(VJcat6, UVcat6, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[0,1].set_xlim(0,2.5)
ax[0,1].set_ylim(0,2.5)
ax[0,1].text(0.35,0.925,"1.25 < z <1.75",ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)

# cat7
hb7 = ax[1,0].hexbin(VJcat7, UVcat7, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,0].set_xlim(0,2.5)
ax[1,0].set_ylim(0,2.5)
ax[1,0].set_xlabel("VJ", fontsize=12)
ax[1,0].set_ylabel("UV", fontsize=12)
ax[1,0].text(0.25,0.925,"1.5 < z < 2",ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)

# cat8
hb8 = ax[1,1].hexbin(VJcat8, UVcat8, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,1].set_xlim(0,2.5)
ax[1,1].set_ylim(0,2.5)
ax[1,1].set_xlabel("VJ", fontsize=12)
ax[1,1].text(0.35,0.925,"1.75 < z < 2.25",ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)

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
ax[0,0].set_xlim(0,2.5)
ax[0,0].set_ylim(0,2.5)
ax[0,0].set_ylabel("UV", fontsize=12)
ax[0,0].text(0.25,0.925,"2 < z < 2.5",ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)

# cat10
hb10 = ax[0,1].hexbin(VJcat10, UVcat10, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[0,1].set_xlim(0,2.5)
ax[0,1].set_ylim(0,2.5)
ax[0,1].text(0.35,0.925,"2.25 < z <2.75",ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)

# cat11
hb11 = ax[1,0].hexbin(VJcat11, UVcat11, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,0].set_xlim(0,2.5)
ax[1,0].set_ylim(0,2.5)
ax[1,0].set_xlabel("VJ", fontsize=12)
ax[1,0].set_ylabel("UV", fontsize=12)
ax[1,0].text(0.25,0.925,"2.5 < z < 3",ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)

# cat12
hb12 = ax[1,1].hexbin(VJcat12, UVcat12, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,1].set_xlim(0,2.5)
ax[1,1].set_ylim(0,2.5)
ax[1,1].set_xlabel("VJ", fontsize=12)
ax[1,1].text(0.35,0.925,"2.75 < z < 3.25",ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)

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
ax[0,0].set_xlim(0,2.5)
ax[0,0].set_ylim(0,2.5)
ax[0,0].set_ylabel("UV", fontsize=12)
ax[0,0].text(0.25,0.925,"3 < z < 3.5",ha='right',va='bottom',transform=ax[0,0].transAxes,fontsize=12)

# cat14
hb14 = ax[0,1].hexbin(VJcat14, UVcat14, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[0,1].set_xlim(0,2.5)
ax[0,1].set_ylim(0,2.5)
ax[0,1].text(0.35,0.925,"3.25 < z <3.75",ha='right',va='bottom',transform=ax[0,1].transAxes,fontsize=12)

# cat15
hb11 = ax[1,0].hexbin(VJcat15, UVcat15, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(173,100))
ax[1,0].set_xlim(0,2.5)
ax[1,0].set_ylim(0,2.5)
ax[1,0].set_xlabel("VJ", fontsize=12)
ax[1,0].set_ylabel("UV", fontsize=12)
ax[1,0].text(0.25,0.925,"3.5 < z < 4",ha='right',va='bottom',transform=ax[1,0].transAxes,fontsize=12)

# cat16
hb16 = ax[1,1].hexbin(VJ, UV, vmax = 100, cmap = "binary", mincnt = 0, gridsize=(866,500))
ax[1,1].set_xlim(0,2.5)
ax[1,1].set_ylim(0,2.5)
ax[1,1].set_xlabel("VJ", fontsize=12)
ax[1,1].text(0.35,0.925,"0 < z < 4",ha='right',va='bottom',transform=ax[1,1].transAxes,fontsize=12)



# ADJUSTED COLORBAR FOR
figure, ax1 = plt.subplots()
#with adjusted colorbar sp it is not saturated

hb17 = ax1.hexbin(VJ, UV, vmax = 1000, cmap = "binary", mincnt = 0, gridsize=(173,100))
cb16 = figure.colorbar(hb17, ax = ax1)
ax1.set_xlim(0,2.5)
ax1.set_ylim(0,2.5)
ax1.set_xlabel("VJ", fontsize=12)
ax1.title.set_text('Non-saturated UVJ (all redshifts, 0-4)')
ax1.text(0.1,2.45,"0 < z < 4",ha='left',va='top',fontsize=12)


plt.show()



# Close FITS file so it won't use up excess memory
hdu.close()