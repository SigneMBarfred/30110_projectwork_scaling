# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 09:38:48 2023

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
t = Table.read('filtered_SED.fits', format='fits')

# Add UV and VJ column
UV = np.transpose(data['restU']) - np.transpose(data['restV'])
UVerr  = np.transpose(data['restU_err']) - np.transpose(data['restU_err'])
VJ = np.transpose(data['restV']) - np.transpose(data['restJ'])
VJerr = np.transpose(data['restV_err']) - np.transpose(data['restJ_err'])

t.add_column(np.transpose(UV), name='UV')
t.add_column(np.transpose(VJ), name='VJ')
t.add_column(np.transpose(UVerr), name='UV_err')
t.add_column(np.transpose(VJerr), name='VJ_err')

# Remove columns not important to UVJ
nimport = data.names[3:116]
t.remove_columns(nimport)

t_array = np.array([t['ID_1'],t['ALPHA_J2000'],t['DELTA_J2000'],t['z_phot'],
        t['restU'],t['restU_err'],t['restB'],t['restB_err'],t['restV'],
        t['restV_err'],t['restJ'],t['restJ_err'],t['mass'],t['sfr'],
        t['Av'],t['UV'],t['VJ'],t['UV_err'],t['VJ_err']])

t_array = np.transpose(t_array.reshape(19, 511006))

df = pd.DataFrame(t_array, columns = t.columns)

##############

df.plot.hexbin(x = 'VJ', y = 'UV', gridsize=50, cmap='inferno') #adds the points to the plot


fig, ax = df.plot(figsize=(4,4))
ax.set_box_aspect(1)
ax.set_title("UVJ diagram for 0.1 < z > 0.3") #defines the title
ax.set_xlabel("V - J (AB mag)") #defines x-axis label
ax.set_ylabel("U - V (AB mag)") #defines y-axis label
ax.text(0.05,0.95,'Quiescent',ha='left',va='top',transform=ax.transAxes,fontsize=10)
hb = ax.hexbin(x = 'VJ', y = 'UV', gridsize=80, cmap='inferno') #adds the points to the plot
cb = fig.colorbar(hb, ax=ax)
cb.set_label('counts')
ax.plot([-0.5, 0.693], [1.3, 1.3], color = "r") #draws the horizontal line
ax.plot([1.6, 1.6], [2.1, 4], color = "r") #draws the vertical line
ax.plot([0.693, 1.6], [1.3, 2.098], color = "r") #draws the tilted line
ax.set_xlim(0,4) #defines range of x-axis
ax.set_ylim(-0.5,3.5) #defines range of y-axis




#############

# UVJ diagram for alle gallaxer
plt.figure(figsize=(4, 4))
plt.title("UVJ diagram") #defines the title
plt.xlabel("V - J (AB mag)") #defines x-axis label
plt.ylabel("U - V (AB mag)") #defines y-axis label
plt.hexbin(VJ, UV, C=None, gridsize=100, bins=None, xscale='linear', yscale='linear', 
           extent=None, cmap=None, norm=None, vmin=None, vmax=None, alpha=None, 
           linewidths=None, edgecolors='face', mincnt=None, 
           marginals=False)
plt.xlim(0,4) #defines range of x-axis
plt.ylim(-0.5,3.5) #defines range of y-axis
plt.show() #shows the plot

# t.index_column('z_phot') # Index for that column

##############
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

# t_UVJ = Table([[data[data.names[0]]],[data[data.names[1]]],[data[data.names[2]]],
#     [data[data.names[3]]],[data[data.names[8]]],[data[data.names[9]]],
#     [data[data.names[10]]],[data[data.names[11]]],[data[data.names[36]]],
#     [data[data.names[37]]],[data[data.names[38]]],[data[data.names[39]]],
#     [data[data.names[116]]],[data[data.names[117]]],[data[data.names[118]]],
#     [data[data.names[119]]],[data[data.names[120]]],[data[data.names[121]]],
#     [data[data.names[122]]],[data[data.names[123]]],[data[data.names[124]]]
#     [data[data.names[125]]],[data[data.names[126]]],[data[data.names[127]]]], names=(data.names))

# Close FITS file so it won't use up excess memory
hdu.close()







