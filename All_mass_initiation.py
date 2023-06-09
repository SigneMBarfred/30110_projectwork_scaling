# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 09:38:22 2023

@author: Nikolaj Lange Dons
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table


#loading data as recarrays
hdul = fits.open('COSMOS2020_All_Mass.fits',memmap=True)
data = hdul[1].data

# Remove noise and unreliable values
for colname in data.names:
    data = data[data['z_phot'] > 0]
    data = data[data['mass'] < 10**16]
    data = data[data['mass'] > 10**4]
    data = data[data['UVISTA_Ks_FLUX'] != -99]
    data = data[data['CFHT_u_FLUX'] != -99]
    data = data[data['HSC_r_FLUX'] != -99]
    data = data[data['UVISTA_J_FLUX'] != -99]

# Replace -99 with nan
Mask = np.array([[data[data.names[4]]==-99],[data[data.names[5]]==-99],[data[data.names[6]]==-99],
                [data[data.names[7]]==-99],[data[data.names[8]]==-99],[data[data.names[9]]==-99],
                [data[data.names[10]]==-99],[data[data.names[11]]==-99],[data[data.names[12]]==-99],
                [data[data.names[13]]==-99],[data[data.names[14]]==-99],[data[data.names[15]]==-99],
                [data[data.names[16]]==-99],[data[data.names[17]]==-99],[data[data.names[18]]==-99],
                [data[data.names[19]]==-99],[data[data.names[20]]==-99],[data[data.names[21]]==-99],
                [data[data.names[22]]==-99],[data[data.names[23]]==-99],[data[data.names[24]]==-99],
                [data[data.names[25]]==-99],[data[data.names[26]]==-99],[data[data.names[27]]==-99],
                [data[data.names[28]]==-99],[data[data.names[29]]==-99],[data[data.names[30]]==-99]])

# for i in range(4,128):
#     data[data.names[i]][Mask[(i-4),0,:]] = np.nan

for i in range(4,31):
    data[data.names[i]][Mask[(i-4),0,:]] = np.nan

# Change back to HDU list to save as fits
t = Table([[data[data.names[0]]],[data[data.names[1]]],[data[data.names[2]]],[data[data.names[3]]],[data[data.names[4]]],[data[data.names[5]]],[data[data.names[6]]],
                [data[data.names[7]]],[data[data.names[8]]],[data[data.names[9]]],
                [data[data.names[10]]],[data[data.names[11]]],[data[data.names[12]]],
                [data[data.names[13]]],[data[data.names[14]]],[data[data.names[15]]],
                [data[data.names[16]]],[data[data.names[17]]],[data[data.names[18]]],
                [data[data.names[19]]],[data[data.names[20]]],[data[data.names[21]]],
                [data[data.names[22]]],[data[data.names[23]]],[data[data.names[24]]],
                [data[data.names[25]]],[data[data.names[26]]],[data[data.names[27]]],
                [data[data.names[28]]],[data[data.names[29]]],[data[data.names[30]]],], names=(data.names))

# Saves changes to new fits
t.write('All_mass_filtered.fits', format='fits', overwrite=True)


# Close FITS file so it won't use up excess memory
hdul.close()