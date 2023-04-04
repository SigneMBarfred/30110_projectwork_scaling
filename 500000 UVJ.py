# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 09:38:48 2023

@author: Nikolaj Lange Dons
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table

#loading data as recarrays
hdu = fits.open('filtered_SED.fits')
data = hdu[1].data
t = Table.read('filtered_SED.fits', format='fits',memmap=True)

#defines UV and VJ
UV = np.transpose(data['restU']) - np.transpose(data['restV'])
VJ = np.transpose(data['restV']) - np.transpose(data['restJ'])

UV = UV.reshape(511006)
VJ = VJ.reshape(511006)

UV = np.float64(UV)
VJ = np.float64(VJ)

# UVJ diagram for all galaxies
fig, ax = plt.subplots(figsize=(4, 4))
h = ax.hexbin(VJ, UV, gridsize = 100, cmap = "inferno", vmax = 1000)
cb = fig.colorbar(h, ax=ax)
cb = cb.set_label('counts')
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_xlabel("V - J (AB mag)")
ax.set_ylabel("U - V (AB mag)")
ax.set_title("UVJ diagram")

"""
# Divide in to catogories
cat1_mask = np.logical_and(t['z_phot'] > 0.1, t['z_phot'] < 0.3)
cat1_index = [i for i, x in enumerate(cat1_mask[0,:]) if x]
cat1_restU = np.transpose(t['restU'][cat1_mask])
cat1_restV = np.transpose(t['restV'][cat1_mask])
cat1_restJ = np.transpose(t['restJ'][cat1_mask])

cat2_mask = np.logical_and(t['z_phot'] > 0.9, t['z_phot'] < 1.1)
cat2_index = [i for i, x in enumerate(cat2_mask[0,:]) if x]
cat2_restU = np.transpose(t['restU'][cat2_mask])
cat2_restV = np.transpose(t['restV'][cat2_mask])
cat2_restJ = np.transpose(t['restJ'][cat2_mask])

cat3_mask = np.logical_and(t['z_phot'] > 1.9, t['z_phot'] < 2.1)
cat3_index = [i for i, x in enumerate(cat3_mask[0,:]) if x]
cat3_restU = np.transpose(t['restU'][cat3_mask])
cat3_restV = np.transpose(t['restV'][cat3_mask])
cat3_restJ = np.transpose(t['restJ'][cat3_mask])

# Create UV and VJ arrays
UVcat1 = cat1_restU - cat1_restV
VJcat1 = cat1_restV - cat1_restJ
UVcat2 = cat2_restU - cat2_restV
VJcat2 = cat2_restV - cat2_restJ
UVcat3 = cat3_restU - cat3_restV
VJcat3 = cat3_restV - cat3_restJ

# sf boundary

cat1UVVJ = 0.88 * VJcat1 + 0.69

cat1_col = np.array([[-2.5*np.log10(t['CFHT_u_MAG'][cat1_mask]/t['SC_IB574_MAG'][cat1_mask])],
         [-2.5*np.log10(t['SC_IB574_MAG'][cat1_mask]/t['UVISTA_J_MAG'][cat1_mask])]])



# UVJ diagram for cat1
fig, ax = plt.subplots(figsize=(4,4))
ax.set_box_aspect(1)
ax.set_title("UVJ diagram for 0.1 < z > 0.3") #defines the title
ax.set_xlabel("V - J (AB mag)") #defines x-axis label
ax.set_ylabel("U - V (AB mag)") #defines y-axis label
ax.text(0.05,0.95,'Quiescent',ha='left',va='top',transform=ax.transAxes,fontsize=10)
hb = ax.hexbin(VJcat1, UVcat1, gridsize=80, cmap='inferno') #adds the points to the plot
cb = fig.colorbar(hb, ax=ax)
cb.set_label('counts')
ax.plot([-0.5, 0.693], [1.3, 1.3], color = "r") #draws the horizontal line
ax.plot([1.6, 1.6], [2.1, 4], color = "r") #draws the vertical line
ax.plot([0.693, 1.6], [1.3, 2.098], color = "r") #draws the tilted line
ax.set_xlim(0,4) #defines range of x-axis
ax.set_ylim(-0.5,3.5) #defines range of y-axis

# UVJ diagram for cat2
fig, ax = plt.subplots(figsize=(4,4))
ax.set_box_aspect(1)
ax.set_title("UVJ diagram for 0.9 < z > 1.1") #defines the title
ax.set_xlabel("V - J (AB mag)") #defines x-axis label
ax.set_ylabel("U - V (AB mag)") #defines y-axis label
ax.text(0.05,0.95,'Quiescent',ha='left',va='top',transform=ax.transAxes,fontsize=10)
hb = ax.hexbin(VJcat2, UVcat2, gridsize=80, cmap='inferno') #adds the points to the plot
cb = fig.colorbar(hb, ax=ax)
cb.set_label('counts')
ax.plot([-0.5, 0.693], [1.3, 1.3], color = "r") #draws the horizontal line
ax.plot([1.6, 1.6], [2.1, 4], color = "r") #draws the vertical line
ax.plot([0.693, 1.6], [1.3, 2.098], color = "r") #draws the tilted line
ax.set_xlim(0,4) #defines range of x-axis
ax.set_ylim(-0.5,3.5) #defines range of y-axis


# UVJ diagram for cat3
fig, ax = plt.subplots(figsize=(4,4))
ax.set_box_aspect(1)
ax.set_title("UVJ diagram for 1.9 < z > 2.1") #defines the title
ax.set_xlabel("V - J (AB mag)") #defines x-axis label
ax.set_ylabel("U - V (AB mag)") #defines y-axis label
ax.text(0.05,0.95,'Quiescent',ha='left',va='top',transform=ax.transAxes,fontsize=10)
hb = ax.hexbin(VJcat2, UVcat2, gridsize=80, cmap='inferno') #adds the points to the plot
cb = fig.colorbar(hb, ax=ax)
cb.set_label('counts')
ax.plot([-0.5, 0.693], [1.3, 1.3], color = "r") #draws the horizontal line
ax.plot([1.6, 1.6], [2.1, 4], color = "r") #draws the vertical line
ax.plot([0.693, 1.6], [1.3, 2.098], color = "r") #draws the tilted line
ax.set_xlim(0,4) #defines range of x-axis
ax.set_ylim(-0.5,3.5) #defines range of y-axis

"""

# Close FITS file so it won't use up excess memory
hdu.close()