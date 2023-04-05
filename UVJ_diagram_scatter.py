# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 00:24:57 2023

@author: Nikolaj Lange Dons
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Column
import matplotlib.pyplot as plt
from astropy.table import Table
import os

# Set directory
os.chdir('C://Users//Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet\Dokumenter/4 semester/Fagprojekt/Data behandling')

#loading data as recarrays
hdu = fits.open('filtered_SED.fits',memmap=True)
data = hdu[1].data

selection_all = np.logical_and((data['z_phot'] > 0), (data['z_phot'] < 4), (data['CFHT_u_MAG'] < 50))

U = data['restU'][0,:]
V = data['restV'][0,:]
J = data['restJ'][0,:]
sfr = data['sfr'][0,:]
mass = data['mass'][0,:]

UV = U[selection_all[0,:]] - V[selection_all[0,:]]
VJ = V[selection_all[0,:]] - J[selection_all[0,:]]
sfr_sel = sfr[selection_all[0,:]]

ssfr = np.log10(sfr[selection_all[0,:]])- mass[selection_all[0,:]]+9.

UV.reshape(476363)
VJ.reshape(476363)

test0 = "%d galaxies"
int0 = len(UV)
print(test0%int0)

# UVJ diagram for galaxies at all redshift
fig0 , ax = plt.subplots(figsize=(8,6))
hb0 = ax.scatter(VJ, UV, c=sfr_sel, vmin=0,vmax=20,cmap='jet',alpha=0.4)
cbar0 = fig0.colorbar(hb0, ax=ax ) # Adds a colorbar
cbar0.set_label('counts',rotation=270,labelpad=20, fontsize=12)
plt.xlim(-4,3)
plt.ylim(-4,3)
plt.xlabel('V - J [AB mag]',fontsize=12)
plt.ylabel('U - V [AB mag]',fontsize=12)
plt.title('UVJ diagram', fontsize=18)

######################

# UVJ diagram for cat1
# 0.1 < z < 0.3

selection_cat1 = np.logical_and((data['z_phot'] > 0.1), (data['z_phot'] < 0.3), (data['CFHT_u_MAG'] < 50))

UV_cat1 = U[selection_cat1[0,:]] - V[selection_cat1[0,:]]
VJ_cat1 = V[selection_cat1[0,:]] - J[selection_cat1[0,:]]
sfr_cat1 = sfr[selection_cat1[0,:]]

fig1, ax = plt.subplots(figsize=(8,6))
hb = ax.scatter(VJ_cat1, UV_cat1, c=sfr_cat1, vmin=0,vmax=20,cmap='jet',alpha=0.4)
ax.text(0.05,0.95,'Quiescent',ha='left',va='top',transform=ax.transAxes,fontsize=14)
cbar1 = fig1.colorbar(hb, ax=ax) # Adds a colorbar
cbar1.set_label('counts',rotation=270,labelpad=20, fontsize=12)
ax.plot([-4, 0.693], [1.3, 1.3], color = "r") #draws the horizontal line
ax.plot([1.6, 1.6], [2.098, 4], color = "r") #draws the vertical line
ax.plot([0.693, 1.6], [1.3, 2.098], color = "r") #draws the tilted line
ax.set_xlim(-4,3)
ax.set_ylim(-4,3)
ax.set_xlabel("V - J [AB mag]", fontsize=12)
ax.set_ylabel("U - V [AB mag]", fontsize=12)
ax.set_title("UVJ diagram for 0.1 < z < 0.3", fontsize=18)

######################

# UVJ diagram for cat2
# 0.9 < z < 1.1

selection_cat2 = np.logical_and((data['z_phot'] > 0.9), (data['z_phot'] < 1.1), (data['CFHT_u_MAG'] < 50))

UV_cat2 = U[selection_cat2[0,:]] - V[selection_cat2[0,:]]
VJ_cat2 = V[selection_cat2[0,:]] - J[selection_cat2[0,:]]
sfr_cat2 = sfr[selection_cat2[0,:]]

# P3 = plt.figure(3)
# plt.scatter(VJ_cat2, UV_cat2, c=sfr_cat2, vmin=0,vmax=20,cmap='jet',alpha=0.4)
# cbar3 = plt.colorbar() # Adds a colorbar
# cbar3.set_label('counts',rotation=270,labelpad=20)
# plt.xlim(-4,3)
# plt.ylim(-4,3)
# plt.xlabel('V-J')
# plt.ylabel('U-V')
# plt.title('UVJ diagram for 0.9 < z < 1.1')
# P3.show()

fig2, ax = plt.subplots(figsize=(8,6))
hb2 = ax.scatter(VJ_cat2, UV_cat2, c=sfr_cat2, vmin=0,vmax=20,cmap='jet',alpha=0.4)
ax.text(0.05,0.95,'Quiescent',ha='left',va='top',transform=ax.transAxes,fontsize=14)
cbar2 = fig2.colorbar(hb2, ax=ax) # Adds a colorbar
cbar2.set_label('counts',rotation=270,labelpad=20, fontsize=12)
ax.plot([-4.0, 0.81], [1.3, 1.3], color = "r") #draws the horizontal line
ax.plot([1.6, 1.6], [2.0, 4], color = "r") #draws the vertical line
ax.plot([0.81, 1.6], [1.3, 2.0], color = "r") #draws the tilted line
ax.set_xlim(-4,3)
ax.set_ylim(-4,3)
ax.set_xlabel("V - J [AB mag]", fontsize=12)
ax.set_ylabel("U - V [AB mag]", fontsize=12)
ax.set_title("UVJ diagram for 0.9 < z < 1.1", fontsize=18)

######################

# UVJ diagram for cat3
# 1.9 < z < 2.1

selection_cat3 = np.logical_and((data['z_phot'] > 1.9), (data['z_phot'] < 2.1), (data['CFHT_u_MAG'] < 50))

UV_cat3 = U[selection_cat3[0,:]] - V[selection_cat3[0,:]]
VJ_cat3 = V[selection_cat3[0,:]] - J[selection_cat3[0,:]]
sfr_cat3 = sfr[selection_cat3[0,:]]

fig3, ax = plt.subplots(figsize=(8,6))
hb3 = ax.scatter(VJ_cat3, UV_cat3, c=sfr_cat3, vmin=0,vmax=20,cmap='jet',alpha=0.4)
ax.text(0.05,0.95,'Quiescent',ha='left',va='top',transform=ax.transAxes,fontsize=14)
cbar3 = fig3.colorbar(hb3, ax=ax) # Adds a colorbar
cbar3.set_label('counts',rotation=270,labelpad=20, fontsize=12)
ax.plot([-4.0, 0.92], [1.3, 1.3], color = "r") #draws the horizontal line
ax.plot([1.6, 1.6], [1.9, 4], color = "r") #draws the vertical line
ax.plot([0.92, 1.6], [1.3, 1.9], color = "r") #draws the tilted line
ax.set_xlim(-4,3)
ax.set_ylim(-4,3)
ax.set_xlabel("V - J [AB mag]", fontsize=12)
ax.set_ylabel("U - V [AB mag]", fontsize=12)
ax.set_title("UVJ diagram for 1.9 < z < 2.1", fontsize=18)

## En anden plot metode
# P4 = plt.figure(4)
# plt.scatter(VJ_cat3, UV_cat3, c=sfr_cat3, vmin=0,vmax=20,cmap='jet',alpha=0.4)
# cbar4 = plt.colorbar() # Adds a colorbar
# cbar4.set_label('counts',rotation=270,labelpad=20)
# plt([-4.0, 0.92], [1.3, 1.3], color = "r") #draws the horizontal line
# plt([1.6, 1.6], [1.9, 4], color = "r") #draws the vertical line
# plt([0.92, 1.6], [1.3, 1.9], color = "r") #draws the tilted line
# plt.xlim(-4,3)
# plt.ylim(-4,3)
# plt.xlabel('V-J')
# plt.ylabel('U-V')
# plt.title('UVJ diagram for 1.9 < z < 2.1')

# Close FITS file so it won't use up excess memory
hdu.close()
