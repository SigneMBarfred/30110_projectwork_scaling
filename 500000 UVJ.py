# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 00:24:57 2023

@author: Nikolaj Lange Dons
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

test0 = "%d galaxies"
int0 = len(UV)
print(test0%int0)

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

# UVJ diagram for cat1
# 0.1 < z < 0.3

selection_cat1 = np.logical_and((data['z_phot'] > 0.1), (data['z_phot'] < 0.3), (data['CFHT_u_MAG'] < 50))

UV_cat1 = U[selection_cat1[0,:]] - V[selection_cat1[0,:]]
VJ_cat1 = V[selection_cat1[0,:]] - J[selection_cat1[0,:]]

fig1, ax = plt.subplots(figsize = (8, 6))
hb1 = ax.hexbin(VJ_cat1, UV_cat1, vmax = 50, cmap = "hot", gridsize = (173,100), mincnt = 0)
cb1 = fig1.colorbar(hb1, ax = ax)
cb1.set_label("Number of galaxies")
ax.plot([-4, 0.693], [1.3, 1.3], color = "g", lw = 5) #draws the horizontal line
ax.plot([1.6, 1.6], [2.098, 4], color = "g", lw = 5) #draws the vertical line
ax.plot([0.693, 1.6], [1.3, 2.098], color = "g", lw = 5) #draws the tilted line
ax.set_xlim(0,2.5)
ax.set_ylim(0,2.5)
ax.set_xlabel("V - J [AB mag]", fontsize=12)
ax.set_ylabel("U - V [AB mag]", fontsize=12)
ax.set_title("UVJ diagram for 0.1 < z < 0.3", fontsize=18)

######################

# UVJ diagram for cat2
# 0.9 < z < 1.1

selection_cat2 = np.logical_and((data['z_phot'] > 0.9), (data['z_phot'] < 1.1), (data['CFHT_u_MAG'] < 50))

UV_cat2 = U[selection_cat2[0,:]] - V[selection_cat2[0,:]]
VJ_cat2 = V[selection_cat2[0,:]] - J[selection_cat2[0,:]]

fig2, ax = plt.subplots(figsize = (8, 6))
hb2 = ax.hexbin(VJ_cat2, UV_cat2, vmax = 50, cmap = "hot", gridsize = (173,100), mincnt = 0)
cb2 = fig2.colorbar(hb2, ax = ax)
cb2.set_label("Number of galaxies")
ax.plot([-4.0, 0.81], [1.3, 1.3], color = "g", lw = 5) #draws the horizontal line
ax.plot([1.6, 1.6], [2.0, 4], color = "g", lw = 5) #draws the vertical line
ax.plot([0.81, 1.6], [1.3, 2.0], color = "g", lw = 5) #draws the tilted line
ax.set_xlim(0,2.5)
ax.set_ylim(0,2.5)
ax.set_xlabel("V - J [AB mag]", fontsize=12)
ax.set_ylabel("U - V [AB mag]", fontsize=12)
ax.set_title("UVJ diagram for 0.9 < z < 1.1", fontsize=18)

######################

# UVJ diagram for cat3
# 1.9 < z < 2.1

selection_cat3 = np.logical_and((data['z_phot'] > 1.9), (data['z_phot'] < 2.1), (data['CFHT_u_MAG'] < 50))

UV_cat3 = U[selection_cat3[0,:]] - V[selection_cat3[0,:]]
VJ_cat3 = V[selection_cat3[0,:]] - J[selection_cat3[0,:]]

fig3, ax = plt.subplots(figsize = (8, 6))
hb3 = ax.hexbin(VJ_cat3, UV_cat3, vmax = 50, cmap = "hot", gridsize = (173,100), mincnt = 0)
cb3 = fig3.colorbar(hb3, ax = ax)
cb3.set_label("Number of galaxies")
ax.plot([-4.0, 0.92], [1.3, 1.3], color = "g", lw = 5) #draws the horizontal line
ax.plot([1.6, 1.6], [1.9, 4], color = "g", lw = 5) #draws the vertical line
ax.plot([0.92, 1.6], [1.3, 1.9], color = "g", lw = 5) #draws the tilted line
ax.set_xlim(0,2.5)
ax.set_ylim(0,2.5)
ax.set_xlabel("V - J [AB mag]", fontsize=12)
ax.set_ylabel("U - V [AB mag]", fontsize=12)
ax.set_title("UVJ diagram for 1.9 < z < 2.1", fontsize=18)

# Close FITS file so it won't use up excess memory
hdu.close()