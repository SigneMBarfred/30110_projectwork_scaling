# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 19:53:48 2023

@author: Nikolaj Lange Dons
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table

#loading data as recarrays
hdul = fits.open('COSMOS2020_best.fits',memmap=True)
data = hdul[1].data

# Remove noise and unreliable values
for colname in data.names:
    data = data[data['mass'] > 10**8]
    data = data[data['mass'] < 10**13]
    data = data[data['z_phot'] > 0]
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
                [data[data.names[28]]==-99],[data[data.names[29]]==-99],[data[data.names[30]]==-99],
                [data[data.names[31]]==-99],[data[data.names[32]]==-99],[data[data.names[33]]==-99],
                [data[data.names[34]]==-99],[data[data.names[35]]==-99],[data[data.names[36]]==-99],
                [data[data.names[37]]==-99],[data[data.names[38]]==-99],[data[data.names[39]]==-99],
                [data[data.names[40]]==-99],[data[data.names[41]]==-99],[data[data.names[42]]==-99],
                [data[data.names[43]]==-99],[data[data.names[44]]==-99],[data[data.names[45]]==-99],
                [data[data.names[46]]==-99],[data[data.names[47]]==-99],[data[data.names[48]]==-99],
                [data[data.names[49]]==-99],[data[data.names[50]]==-99],[data[data.names[51]]==-99],
                [data[data.names[52]]==-99],[data[data.names[53]]==-99],[data[data.names[54]]==-99],
                [data[data.names[55]]==-99],[data[data.names[56]]==-99],[data[data.names[57]]==-99],
                [data[data.names[58]]==-99],[data[data.names[59]]==-99],[data[data.names[60]]==-99],
                [data[data.names[61]]==-99],[data[data.names[62]]==-99],[data[data.names[63]]==-99],
                [data[data.names[64]]==-99],[data[data.names[65]]==-99],[data[data.names[66]]==-99],
                [data[data.names[67]]==-99],[data[data.names[68]]==-99],[data[data.names[69]]==-99],
                [data[data.names[70]]==-99],[data[data.names[71]]==-99],[data[data.names[72]]==-99],
                [data[data.names[73]]==-99],[data[data.names[74]]==-99],[data[data.names[75]]==-99],
                [data[data.names[76]]==-99],[data[data.names[77]]==-99],[data[data.names[78]]==-99],
                [data[data.names[79]]==-99],[data[data.names[80]]==-99],[data[data.names[81]]==-99],
                [data[data.names[82]]==-99],[data[data.names[83]]==-99],[data[data.names[84]]==-99],
                [data[data.names[85]]==-99],[data[data.names[86]]==-99],[data[data.names[87]]==-99],
                [data[data.names[88]]==-99],[data[data.names[89]]==-99],[data[data.names[90]]==-99],
                [data[data.names[91]]==-99],[data[data.names[92]]==-99],[data[data.names[93]]==-99],
                [data[data.names[94]]==-99],[data[data.names[95]]==-99],[data[data.names[96]]==-99],
                [data[data.names[97]]==-99],[data[data.names[98]]==-99],[data[data.names[99]]==-99],
                [data[data.names[100]]==-99],[data[data.names[101]]==-99],[data[data.names[102]]==-99],
                [data[data.names[103]]==-99],[data[data.names[104]]==-99],[data[data.names[105]]==-99],
                [data[data.names[106]]==-99],[data[data.names[107]]==-99],[data[data.names[108]]==-99],
                [data[data.names[109]]==-99],[data[data.names[110]]==-99],[data[data.names[111]]==-99],
                [data[data.names[112]]==-99],[data[data.names[113]]==-99],[data[data.names[114]]==-99],
                [data[data.names[115]]==-99],[data[data.names[116]]==-99],[data[data.names[117]]==-99],
                [data[data.names[118]]==-99],[data[data.names[119]]==-99],[data[data.names[120]]==-99],
                [data[data.names[121]]==-99],[data[data.names[122]]==-99],[data[data.names[123]]==-99],
                [data[data.names[124]]==-99],[data[data.names[125]]==-99],[data[data.names[126]]==-99],
                [data[data.names[127]]==-99]])

for i in range(4,128):
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
                [data[data.names[28]]],[data[data.names[29]]],[data[data.names[30]]],
                [data[data.names[31]]],[data[data.names[32]]],[data[data.names[33]]],
                [data[data.names[34]]],[data[data.names[35]]],[data[data.names[36]]],
                [data[data.names[37]]],[data[data.names[38]]],[data[data.names[39]]],
                [data[data.names[40]]],[data[data.names[41]]],[data[data.names[42]]],
                [data[data.names[43]]],[data[data.names[44]]],[data[data.names[45]]],
                [data[data.names[46]]],[data[data.names[47]]],[data[data.names[48]]],
                [data[data.names[49]]],[data[data.names[50]]],[data[data.names[51]]],
                [data[data.names[52]]],[data[data.names[53]]],[data[data.names[54]]],
                [data[data.names[55]]],[data[data.names[56]]],[data[data.names[57]]],
                [data[data.names[58]]],[data[data.names[59]]],[data[data.names[60]]],
                [data[data.names[61]]],[data[data.names[62]]],[data[data.names[63]]],
                [data[data.names[64]]],[data[data.names[65]]],[data[data.names[66]]],
                [data[data.names[67]]],[data[data.names[68]]],[data[data.names[69]]],
                [data[data.names[70]]],[data[data.names[71]]],[data[data.names[72]]],
                [data[data.names[73]]],[data[data.names[74]]],[data[data.names[75]]],
                [data[data.names[76]]],[data[data.names[77]]],[data[data.names[78]]],
                [data[data.names[79]]],[data[data.names[80]]],[data[data.names[81]]],
                [data[data.names[82]]],[data[data.names[83]]],[data[data.names[84]]],
                [data[data.names[85]]],[data[data.names[86]]],[data[data.names[87]]],
                [data[data.names[88]]],[data[data.names[89]]],[data[data.names[90]]],
                [data[data.names[91]]],[data[data.names[92]]],[data[data.names[93]]],
                [data[data.names[94]]],[data[data.names[95]]],[data[data.names[96]]],
                [data[data.names[97]]],[data[data.names[98]]],[data[data.names[99]]],
                [data[data.names[100]]],[data[data.names[101]]],[data[data.names[102]]],
                [data[data.names[103]]],[data[data.names[104]]],[data[data.names[105]]],
                [data[data.names[106]]],[data[data.names[107]]],[data[data.names[108]]],
                [data[data.names[109]]],[data[data.names[110]]],[data[data.names[111]]],
                [data[data.names[112]]],[data[data.names[113]]],[data[data.names[114]]],
                [data[data.names[115]]],[data[data.names[116]]],[data[data.names[117]]],
                [data[data.names[118]]],[data[data.names[119]]],[data[data.names[120]]],
                [data[data.names[121]]],[data[data.names[122]]],[data[data.names[123]]],
                [data[data.names[124]]],[data[data.names[125]]],[data[data.names[126]]],
                [data[data.names[127]]]], names=(data.names))

t.write('table1.fits', format='fits')

# Saves changes to new fits     
data.writeto('Filtered_hdul')

# Close FITS file so it won't use up excess memory
hdul.close()
