# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 12:49:55 2023

@author: Nikolaj Lange Dons

Mass vs redshift plot, showing observational bias
the redshift vs stellar mass plot as a hexbin
#### 
the data file we have doesnt work to display the full 
range of masses bc the mass boundaries have already been set 
with topcat when reading the file
- therefore get new file 
load the new data as follows 
loading data as recarrays
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os

# Set directory
os.chdir('C://Users//Nikolaj Lange Dons/OneDrive - Danmarks Tekniske Universitet\Dokumenter/4 semester/Fagprojekt/Data behandling')

# Load data
hdu = fits.open('All_mass_filtered.fits',memmap=True)
data = hdu[1].data

# Define U, V and J
selection_all = np.logical_and((data['z_phot'] > 0), (data['z_phot'] < 4))

# Convert restframe micro jansky to AB magnitude
U = 23.9 - 2.5*np.log10(data['restU'][0,:])
V = 23.9 - 2.5*np.log10(data['restV'][0,:])
J = 23.9 - 2.5*np.log10(data['restJ'][0,:])

# Define UV and VJ for 0 < z < 4
UV = U[selection_all[0,:]] - V[selection_all[0,:]]
VJ = V[selection_all[0,:]] - J[selection_all[0,:]]

# Define UV and VJ for all z
UVz = U - V
VJz = V - J

# Define mass, redshift and SFR array
SFR = data['sfr'][0,:]
mass = data['mass'][0,:]
z = data['z_phot'][0,:]

##################
# SFGs and QGs boundaries 

# cat1: 0 < z < 0.5
idUVJ_cat1 = np.logical_and((z > 0), (z < 0.5))
iduv1_cat1 = UVz[idUVJ_cat1] >= 1.3
iduv2_cat1 = UVz[idUVJ_cat1] <= 3
idvj1_cat1 = VJz[idUVJ_cat1] <= 1.6
idvj2_cat1 = VJz[idUVJ_cat1] >= -1
iduvj_cat1 = UVz[idUVJ_cat1] >= 0.88*VJz[idUVJ_cat1]+0.69
IDS_cat1 = iduv1_cat1 & idvj1_cat1 & iduvj_cat1 & idvj2_cat1 & iduv2_cat1

nQ_cat1 = sum(IDS_cat1)
nSF_cat1 = sum(~IDS_cat1)
cat1_total = nQ_cat1 + nSF_cat1
m1 = mass[idUVJ_cat1] # Mass

# cat2: 0.25 < z < 0.75
idUVJ_cat2_lower = np.logical_and((z > 0.25), (z < 0.5)) 
idUVJ_cat2_upper = np.logical_and((z > 0.5), (z < 0.75))

iduv1_cat2_lower = UVz[idUVJ_cat2_lower] >= 1.3
iduv2_cat2_lower = UVz[idUVJ_cat2_lower] <= 3
idvj1_cat2_lower = VJz[idUVJ_cat2_lower] <= 1.6
idvj2_cat2_lower = VJz[idUVJ_cat2_lower] >= -1
iduvj_cat2_lower = UVz[idUVJ_cat2_lower] >= 0.88*VJz[idUVJ_cat2_lower]+0.69
IDS_cat2_lower = iduv1_cat2_lower & idvj1_cat2_lower & iduvj_cat2_lower & idvj2_cat2_lower & iduv2_cat2_lower

iduv1_cat2_upper = UVz[idUVJ_cat2_upper] >= 1.3
iduv2_cat2_upper = UVz[idUVJ_cat2_upper] <= 3
idvj1_cat2_upper = VJz[idUVJ_cat2_upper] <= 1.6
idvj2_cat2_upper = VJz[idUVJ_cat2_upper] >= -1
iduvj_cat2_upper = UVz[idUVJ_cat2_upper] >= 0.88*VJz[idUVJ_cat2_upper]+0.59
IDS_cat2_upper = iduv1_cat2_upper & idvj1_cat2_upper & iduvj_cat2_upper & idvj2_cat2_upper & iduv2_cat2_upper

nQ_cat2 = sum(IDS_cat2_lower)+sum(IDS_cat2_upper)
nSF_cat2 = sum(~IDS_cat2_lower)+sum(~IDS_cat2_upper)
cat2_total = nQ_cat2 + nSF_cat2

# cat3: 0.5 < z < 1.0
idUVJ_cat3 = np.logical_and((z > 0.5), (z < 1))
iduv1_cat3 = UVz[idUVJ_cat3] >= 1.3
iduv2_cat3 = UVz[idUVJ_cat3] <= 3
idvj1_cat3 = VJz[idUVJ_cat3] <= 1.6
idvj2_cat3 = VJz[idUVJ_cat3] >= -1
iduvj_cat3 = UVz[idUVJ_cat3] >= 0.88*VJz[idUVJ_cat3]+0.59
IDS_cat3 = iduv1_cat3 & idvj1_cat3 & iduvj_cat3 & idvj2_cat3 & iduv2_cat3

nQ_cat3 = sum(IDS_cat3)
nSF_cat3 = sum(~IDS_cat3)
cat3_total = nQ_cat3 + nSF_cat3
m3 = mass[idUVJ_cat3] # Mass

# cat4: 0.75 < z < 1.25
idUVJ_cat4_lower = np.logical_and((z > 0.75), (z < 1))
idUVJ_cat4_upper = np.logical_and((z > 0.5), (z < 1.25))

iduv1_cat4_lower = UVz[idUVJ_cat4_lower] >= 1.3
iduv2_cat4_lower = UVz[idUVJ_cat4_lower] <= 3
idvj1_cat4_lower = VJz[idUVJ_cat4_lower] <= 1.6
idvj2_cat4_lower = VJz[idUVJ_cat4_lower] >= -1
iduvj_cat4_lower = UVz[idUVJ_cat4_lower] >= 0.88*VJz[idUVJ_cat4_lower]+0.59
IDS_cat4_lower = iduv1_cat4_lower & idvj1_cat4_lower & iduvj_cat4_lower & idvj2_cat4_lower & iduv2_cat4_lower

iduv1_cat4_upper = UVz[idUVJ_cat4_upper] >= 1.3
iduv2_cat4_upper = UVz[idUVJ_cat4_upper] <= 3
idvj1_cat4_upper = VJz[idUVJ_cat4_upper] <= 1.6
idvj2_cat4_upper = VJz[idUVJ_cat4_upper] >= -1
iduvj_cat4_upper = UVz[idUVJ_cat4_upper] >= 0.88*VJz[idUVJ_cat4_upper]+0.49
IDS_cat4_upper = iduv1_cat4_upper & idvj1_cat4_upper & iduvj_cat4_upper & idvj2_cat4_upper & iduv2_cat4_upper

nQ_cat4 = sum(IDS_cat4_lower)+sum(IDS_cat4_upper)
nSF_cat4 = sum(~IDS_cat4_lower)+sum(~IDS_cat4_upper)
cat4_total = nQ_cat4 + nSF_cat4

# cat5: 1.0 < z < 1.5
idUVJ_cat5 = np.logical_and((z > 1), (z < 1.5))
iduv1_cat5 = UVz[idUVJ_cat5] >= 1.3
iduv2_cat5 = UVz[idUVJ_cat5] <= 3
idvj1_cat5 = VJz[idUVJ_cat5] <= 1.6
idvj2_cat5 = VJz[idUVJ_cat5] >= -1
iduvj_cat5 = UVz[idUVJ_cat5] >= 0.88*VJz[idUVJ_cat5]+0.49
IDS_cat5 = iduv1_cat5 & idvj1_cat5 & iduvj_cat5 & idvj2_cat5 & iduv2_cat5

nQ_cat5 = sum(IDS_cat5)
nSF_cat5 = sum(~IDS_cat5)
cat5_total = nQ_cat5 + nSF_cat5
m5 = mass[idUVJ_cat5] # Mass

# cat6: 1.25 < z < 1.75
idUVJ_cat6 = np.logical_and((z > 1.25), (z < 1.75))
iduv1_cat6 = UVz[idUVJ_cat6] >= 1.3
iduv2_cat6 = UVz[idUVJ_cat6] <= 3
idvj1_cat6 = VJz[idUVJ_cat6] <= 1.6
idvj2_cat6 = VJz[idUVJ_cat6] >= -1
iduvj_cat6 = UVz[idUVJ_cat6] >= 0.88*VJz[idUVJ_cat6]+0.49
IDS_cat6 = iduv1_cat6 & idvj1_cat6 & iduvj_cat6 & idvj2_cat6 & iduv2_cat6

nQ_cat6 = sum(IDS_cat6)
nSF_cat6 = sum(~IDS_cat6)
cat6_total = nQ_cat6 + nSF_cat6

# cat7: 1.5 < z < 2
idUVJ_cat7 = np.logical_and((z > 1.5), (z < 2))
iduv1_cat7 = UVz[idUVJ_cat7] >= 1.3
iduv2_cat7 = UVz[idUVJ_cat7] <= 3
idvj1_cat7 = VJz[idUVJ_cat7] <= 1.6
idvj2_cat7 = VJz[idUVJ_cat7] >= -1
iduvj_cat7 = UVz[idUVJ_cat7] >= 0.88*VJz[idUVJ_cat7]+0.49
IDS_cat7 = iduv1_cat7 & idvj1_cat7 & iduvj_cat7 & idvj2_cat7 & iduv2_cat7

nQ_cat7 = sum(IDS_cat7)
nSF_cat7 = sum(~IDS_cat7)
cat7_total = nQ_cat7 + nSF_cat7
m7 = mass[idUVJ_cat7] # Mass

# cat8: 1.75 < z < 2.25
idUVJ_cat8_lower = np.logical_and((z > 1.75), (z < 2))
idUVJ_cat8_upper = np.logical_and((z > 2), (z < 2.25))

iduv1_cat8_lower = UVz[idUVJ_cat8_lower] >= 1.3
iduv2_cat8_lower = UVz[idUVJ_cat8_lower] <= 3
idvj1_cat8_lower = VJz[idUVJ_cat8_lower] <= 1.6
idvj2_cat8_lower = VJz[idUVJ_cat8_lower] >= -1
iduvj_cat8_lower = UVz[idUVJ_cat8_lower] >= 0.88*VJz[idUVJ_cat8_lower]+0.49
IDS_cat8_lower = iduv1_cat8_lower & idvj1_cat8_lower & iduvj_cat8_lower & idvj2_cat8_lower & iduv2_cat8_lower

iduv1_cat8_upper = UVz[idUVJ_cat8_upper] >= 1.3
iduv2_cat8_upper = UVz[idUVJ_cat8_upper] <= 3
idvj1_cat8_upper = VJz[idUVJ_cat8_upper] <= 1.6
idvj2_cat8_upper = VJz[idUVJ_cat8_upper] >= -1
iduvj_cat8_upper = UVz[idUVJ_cat8_upper] >= 0.88*VJz[idUVJ_cat8_upper]+0.49
IDS_cat8_upper = iduv1_cat8_upper & idvj1_cat8_upper & iduvj_cat8_upper & idvj2_cat8_upper & iduv2_cat8_upper

nQ_cat8 = sum(IDS_cat8_lower)+sum(IDS_cat8_upper)
nSF_cat8 = sum(~IDS_cat8_lower)+sum(~IDS_cat8_upper)
cat8_total = nQ_cat8 + nSF_cat8

# cat9: 2 < z < 2.5
idUVJ_cat9 = np.logical_and((z > 2), (z < 2.5))
iduv1_cat9 = UVz[idUVJ_cat9] >= 1.3
iduv2_cat9 = UVz[idUVJ_cat9] <= 3
idvj1_cat9 = VJz[idUVJ_cat9] <= 1.6
idvj2_cat9 = VJz[idUVJ_cat9] >= -1
iduvj_cat9 = UVz[idUVJ_cat9] >= 0.88*VJz[idUVJ_cat9]+0.49
IDS_cat9 = iduv1_cat9 & idvj1_cat9 & iduvj_cat9 & idvj2_cat9 & iduv2_cat9

nQ_cat9 = sum(IDS_cat9)
nSF_cat9 = sum(~IDS_cat9)
cat9_total = nQ_cat9 + nSF_cat9
m9 = mass[idUVJ_cat9] # Mass

# cat10: 2.25 < z < 2.75
idUVJ_cat10 = np.logical_and((z > 2.25), (z < 2.75))
iduv1_cat10 = UVz[idUVJ_cat10] >= 1.3
iduv2_cat10 = UVz[idUVJ_cat10] <= 3
idvj1_cat10 = VJz[idUVJ_cat10] <= 1.6
idvj2_cat10 = VJz[idUVJ_cat10] >= -1
iduvj_cat10 = UVz[idUVJ_cat10] >= 0.88*VJz[idUVJ_cat10]+0.49
IDS_cat10 = iduv1_cat10 & idvj1_cat10 & iduvj_cat10 & idvj2_cat10 & iduv2_cat10

nQ_cat10 = sum(IDS_cat10)
nSF_cat10 = sum(~IDS_cat10)
cat10_total = nQ_cat10 + nSF_cat10

# cat11: 2.5 < z < 3
idUVJ_cat11 = np.logical_and((z > 2.5), (z < 3))
iduv1_cat11 = UVz[idUVJ_cat11] >= 1.3
iduv2_cat11 = UVz[idUVJ_cat11] <= 3
idvj1_cat11 = VJz[idUVJ_cat11] <= 1.6
idvj2_cat11 = VJz[idUVJ_cat11] >= -1
iduvj_cat11 = UVz[idUVJ_cat11] >= 0.88 * VJz[idUVJ_cat11] + 0.49
IDS_cat11 = iduv1_cat11 & idvj1_cat11 & iduvj_cat11 & idvj2_cat11 & iduv2_cat11

nQ_cat11 = sum(IDS_cat11)
nSF_cat11 = sum(~IDS_cat11)
cat11_total = nQ_cat11 + nSF_cat11
m11 = mass[idUVJ_cat11] # Mass

# cat12: 2.75 < z < 3.25
idUVJ_cat12 = np.logical_and((z > 2.75), (z < 3.25))
iduv1_cat12 = UVz[idUVJ_cat12] >= 1.3
iduv2_cat12 = UVz[idUVJ_cat12] <= 3
idvj1_cat12 = VJz[idUVJ_cat12] <= 1.6
idvj2_cat12 = VJz[idUVJ_cat12] >= -1
iduvj_cat12 = UVz[idUVJ_cat12] >= 0.88 * VJz[idUVJ_cat12] + 0.49
IDS_cat12 = iduv1_cat12 & idvj1_cat12 & iduvj_cat12 & idvj2_cat12 & iduv2_cat12

nQ_cat12 = sum(IDS_cat12)
nSF_cat12 = sum(~IDS_cat12)
cat12_total = nQ_cat12 + nSF_cat12

# cat13: 3 < z < 3.5
idUVJ_cat13 = np.logical_and((z > 3), (z < 3.5))
iduv1_cat13 = UVz[idUVJ_cat13] >= 1.3
iduv2_cat13 = UVz[idUVJ_cat13] <= 3
idvj1_cat13 = VJz[idUVJ_cat13] <= 1.6
idvj2_cat13 = VJz[idUVJ_cat13] >= -1
iduvj_cat13 = UVz[idUVJ_cat13] >= 0.88 * VJz[idUVJ_cat13] + 0.49
IDS_cat13 = iduv1_cat13 & idvj1_cat13 & iduvj_cat13 & idvj2_cat13 & iduv2_cat13

nQ_cat13 = sum(IDS_cat13)
nSF_cat13 = sum(~IDS_cat13)
cat13_total = nQ_cat13 + nSF_cat13
m13 = mass[idUVJ_cat13] # Mass

# cat14: 3.25 < z < 3.75
idUVJ_cat14 = np.logical_and((z > 3.25), (z < 3.75))
iduv1_cat14 = UVz[idUVJ_cat14] >= 1.3
iduv2_cat14 = UVz[idUVJ_cat14] <= 3
idvj1_cat14 = VJz[idUVJ_cat14] <= 1.6
idvj2_cat14 = VJz[idUVJ_cat14] >= -1
iduvj_cat14 = UVz[idUVJ_cat14] >= 0.88 * VJz[idUVJ_cat14] + 0.49
IDS_cat14 = iduv1_cat14 & idvj1_cat14 & iduvj_cat14 & idvj2_cat14 & iduv2_cat14

nQ_cat14 = sum(IDS_cat14)
nSF_cat14 = sum(~IDS_cat14)
cat14_total = nQ_cat14 + nSF_cat14

# cat15: 3.5 < z < 4
idUVJ_cat15 = np.logical_and((z > 3.5), (z < 4))
iduv1_cat15 = UVz[idUVJ_cat15] >= 1.3
iduv2_cat15 = UVz[idUVJ_cat15] <= 3
idvj1_cat15 = VJz[idUVJ_cat15] <= 1.6
idvj2_cat15 = VJz[idUVJ_cat15] >= -1
iduvj_cat15 = UVz[idUVJ_cat15] >= 0.88 * VJz[idUVJ_cat15] + 0.49
IDS_cat15 = iduv1_cat15 & idvj1_cat15 & iduvj_cat15 & idvj2_cat15 & iduv2_cat15

nQ_cat15 = sum(IDS_cat15)
nSF_cat15 = sum(~IDS_cat15)
cat15_total = nQ_cat15 + nSF_cat15
m15 = mass[idUVJ_cat15] #mass

# Create a matrix containing intervals and no. of galaxies in each interval
QSF = np.array([[0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75], 
         [nQ_cat1,nQ_cat2,nQ_cat3,nQ_cat4,nQ_cat5,nQ_cat6,nQ_cat7,nQ_cat8,nQ_cat9,
          nQ_cat10,nQ_cat11,nQ_cat12,nQ_cat13,nQ_cat14,nQ_cat15], 
         [nSF_cat1,nSF_cat2,nSF_cat3,nSF_cat4,nSF_cat5,nSF_cat6,nSF_cat7,nSF_cat8,
          nSF_cat9,nSF_cat10,nSF_cat11,nSF_cat12,nSF_cat13,nSF_cat14,nSF_cat15],
         [cat1_total,cat2_total,cat3_total,cat4_total,cat5_total,cat6_total,
          cat7_total,cat8_total,cat9_total,cat10_total,cat11_total,cat12_total,
          cat13_total,cat14_total,cat15_total]])
# QSF = [z, num quiesent, num starforming, total]
# QSF = [0,       1     ,       2        ,  3   ]

##################
# plot(z, log10(mass))

# 0 < z < 0.5
mean_mass_nSF_cat1 = np.average(m1[~IDS_cat1])
mean_mass_nQ_cat1 = np.average(m1[IDS_cat1])
nSF_error1 = (mean_mass_nSF_cat1/cat1_total)*np.sqrt(1/mean_mass_nSF_cat1 + 1/cat1_total)
NQ_error1 = (mean_mass_nQ_cat1/cat1_total)*np.sqrt(1/mean_mass_nQ_cat1 + 1/cat1_total)

# 0.25 < z < 0.75
mean_mass_cat2_lower = mass[idUVJ_cat2_lower]
mean_mass_cat2_upper = mass[idUVJ_cat2_upper]

mean_mass_nSF_cat2 = (np.average(mean_mass_cat2_lower[~IDS_cat2_lower])+np.average(mean_mass_cat2_upper[~IDS_cat2_upper]))/2
mean_mass_nQ_cat2 = (np.average(mean_mass_cat2_lower[IDS_cat2_lower])+np.average(mean_mass_cat2_upper[IDS_cat2_upper]))/2
nSF_error2 = (mean_mass_nSF_cat2/cat2_total)*np.sqrt(1/mean_mass_nSF_cat2 + 1/cat2_total)
NQ_error2 = (mean_mass_nQ_cat2/cat2_total)*np.sqrt(1/mean_mass_nQ_cat2 + 1/cat2_total)

# 0.5 < z < 1
mean_mass_cat3 = mass[idUVJ_cat3]
mean_mass_nSF_cat3 = np.average(mean_mass_cat3[~IDS_cat3])
mean_mass_nQ_cat3 = np.average(mean_mass_cat3[IDS_cat3])
nSF_error3 = (mean_mass_nSF_cat3/cat3_total)*np.sqrt(1/mean_mass_nSF_cat3 + 1/cat3_total)
NQ_error3 = (mean_mass_nQ_cat3/cat3_total)*np.sqrt(1/mean_mass_nQ_cat3 + 1/cat3_total)

# 0.75 < z < 1.25
mean_mass_cat4_lower = mass[idUVJ_cat4_lower]
mean_mass_cat4_upper = mass[idUVJ_cat4_upper]

mean_mass_nSF_cat4 = (np.average(mean_mass_cat4_lower[~IDS_cat4_lower])+np.average(mean_mass_cat4_upper[~IDS_cat4_upper]))/2
mean_mass_nQ_cat4 = (np.average(mean_mass_cat4_lower[IDS_cat4_lower])+np.average(mean_mass_cat4_upper[IDS_cat4_upper]))/2
nSF_error4 = (mean_mass_nSF_cat4/cat4_total)*np.sqrt(1/mean_mass_nSF_cat4 + 1/cat4_total)
NQ_error4 = (mean_mass_nQ_cat4/cat4_total)*np.sqrt(1/mean_mass_nQ_cat4 + 1/cat4_total)

# 1 < z < 1.5
mean_mass_cat5 = mass[idUVJ_cat5]
mean_mass_nSF_cat5 = np.average(mean_mass_cat5[~IDS_cat5])
mean_mass_nQ_cat5 = np.average(mean_mass_cat5[IDS_cat5])
nSF_error5 = (mean_mass_nSF_cat5/cat5_total)*np.sqrt(1/mean_mass_nSF_cat5 + 1/cat5_total)
NQ_error5 = (mean_mass_nQ_cat5/cat5_total)*np.sqrt(1/mean_mass_nQ_cat5 + 1/cat5_total)

# 1.25 < z < 1.75
mean_mass_cat6 = mass[idUVJ_cat6]
mean_mass_nSF_cat6 = np.average(mean_mass_cat6[~IDS_cat6])
mean_mass_nQ_cat6 = np.average(mean_mass_cat6[IDS_cat6])
nSF_error6 = (mean_mass_nSF_cat6/cat6_total)*np.sqrt(1/mean_mass_nSF_cat6 + 1/cat6_total)
NQ_error6 = (mean_mass_nQ_cat6/cat6_total)*np.sqrt(1/mean_mass_nQ_cat6 + 1/cat6_total)

# 1.5 < z < 2
mean_mass_cat7 = mass[idUVJ_cat7]
mean_mass_nSF_cat7 = np.average(mean_mass_cat7[~IDS_cat7])
mean_mass_nQ_cat7 = np.average(mean_mass_cat7[IDS_cat7])
nSF_error7 = (mean_mass_nSF_cat7/cat7_total)*np.sqrt(1/mean_mass_nSF_cat7 + 1/cat7_total)
NQ_error7 = (mean_mass_nQ_cat7/cat7_total)*np.sqrt(1/mean_mass_nQ_cat7 + 1/cat7_total)

# 1.75 < z < 2.25
mean_mass_cat8_lower = mass[idUVJ_cat8_lower]
mean_mass_cat8_upper = mass[idUVJ_cat8_upper]

mean_mass_nSF_cat8 = (np.average(mean_mass_cat8_lower[~IDS_cat8_lower])+np.average(mean_mass_cat8_upper[~IDS_cat8_upper]))/2
mean_mass_nQ_cat8 = (np.average(mean_mass_cat8_lower[IDS_cat8_lower])+np.average(mean_mass_cat8_upper[IDS_cat8_upper]))/2
nSF_error8 = (mean_mass_nSF_cat8/cat8_total)*np.sqrt(1/mean_mass_nSF_cat8 + 1/cat8_total)
NQ_error8 = (mean_mass_nQ_cat8/cat8_total)*np.sqrt(1/mean_mass_nQ_cat8 + 1/cat8_total)

# 2 < z < 2.5
mean_mass_cat9 = mass[idUVJ_cat9]
mean_mass_nSF_cat9 = np.average(mean_mass_cat9[~IDS_cat9])
mean_mass_nQ_cat9 = np.average(mean_mass_cat9[IDS_cat9])
nSF_error9 = (mean_mass_nSF_cat9/cat9_total)*np.sqrt(1/mean_mass_nSF_cat9 + 1/cat9_total)
NQ_error9 = (mean_mass_nQ_cat9/cat9_total)*np.sqrt(1/mean_mass_nQ_cat9 + 1/cat9_total)

# 2.25 < z < 2.75
mean_mass_cat10 = mass[idUVJ_cat10]
mean_mass_nSF_cat10 = np.average(mean_mass_cat10[~IDS_cat10])
mean_mass_nQ_cat10 = np.average(mean_mass_cat10[IDS_cat10])
nSF_error10 = (mean_mass_nSF_cat10/cat10_total)*np.sqrt(1/mean_mass_nSF_cat10 + 1/cat10_total)
NQ_error10 = (mean_mass_nQ_cat10/cat10_total)*np.sqrt(1/mean_mass_nQ_cat10 + 1/cat10_total)

# 2.5 < z < 3
mean_mass_cat11 = mass[idUVJ_cat11]
mean_mass_nSF_cat11 = np.average(mean_mass_cat11[~IDS_cat11])
mean_mass_nQ_cat11 = np.average(mean_mass_cat11[IDS_cat11])
nSF_error11 = (mean_mass_nSF_cat11/cat11_total)*np.sqrt(1/mean_mass_nSF_cat11 + 1/cat11_total)
NQ_error11 = (mean_mass_nQ_cat11/cat11_total)*np.sqrt(1/mean_mass_nQ_cat11 + 1/cat11_total)

# 2.75 < z < 3.25
mean_mass_cat12 = mass[idUVJ_cat12]
mean_mass_nSF_cat12 = np.average(mean_mass_cat12[~IDS_cat12])
mean_mass_nQ_cat12 = np.average(mean_mass_cat12[IDS_cat12])
nSF_error12 = (mean_mass_nSF_cat12/cat12_total)*np.sqrt(1/mean_mass_nSF_cat12 + 1/cat12_total)
NQ_error12 = (mean_mass_nQ_cat12/cat12_total)*np.sqrt(1/mean_mass_nQ_cat12 + 1/cat12_total)

# 3 < z < 3.5
mean_mass_cat13 = mass[idUVJ_cat13]
mean_mass_nSF_cat13 = np.average(mean_mass_cat13[~IDS_cat13])
mean_mass_nQ_cat13 = np.average(mean_mass_cat13[IDS_cat13])
nSF_error13 = (mean_mass_nSF_cat13/cat13_total)*np.sqrt(1/mean_mass_nSF_cat13 + 1/cat13_total)
NQ_error13 = (mean_mass_nQ_cat13/cat13_total)*np.sqrt(1/mean_mass_nQ_cat13 + 1/cat13_total)

# 3.25 < z < 3.75
mean_mass_cat14 = mass[idUVJ_cat14]
mean_mass_nSF_cat14 = np.average(mean_mass_cat14[~IDS_cat14])
mean_mass_nQ_cat14 = np.average(mean_mass_cat14[IDS_cat14])
nSF_error14 = (mean_mass_nSF_cat14/cat14_total)*np.sqrt(1/mean_mass_nSF_cat14 + 1/cat14_total)
NQ_error14 = (mean_mass_nQ_cat14/cat14_total)*np.sqrt(1/mean_mass_nQ_cat14 + 1/cat14_total)

# 3.5 < z < 4
mean_mass_cat15 = mass[idUVJ_cat15]
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

# Select all valid galaxies
id_positive_mass = mass > 0
id_all_and_mass = selection_all[0,:] & id_positive_mass

z_all = z[id_all_and_mass]
mass_all = mass[id_all_and_mass]

fig, ax = plt.subplots()

SF_and_Q_MASS = ax.hexbin(x = z_all, y = mass_all,
                          vmax = 100, gridsize=(346,200), 
                          cmap='Greys',mincnt=2,yscale = 'log', bins = 'log')

SF = ax.errorbar(mass_QSF[0,:],mass_QSF[2,:], 
                 yerr = mass_QSF[3,:], fmt = "o", elinewidth = 1.7, 
                 capsize = 6, markeredgecolor = "b", 
                 markersize = 4, label='nSF',)

nQ = ax.errorbar(mass_QSF[0,:],mass_QSF[1,:], 
                 yerr = mass_QSF[4,:], fmt = "x", elinewidth = 1.7, 
                 capsize = 6, markeredgecolor = "r", 
                 markersize = 4, label='nQ')

# cb_SF = fig.colorbar(SF_and_Q_MASS)
# cb_SF.set_label("Number of galaxies", fontsize=12)
ax.set_xlim(0,4)
ax.legend([SF,nQ],['av. SF','av. Q'],loc="upper right")
ax.set_ylim(10**6,10**13)
ax.set_yscale('log')
ax.set_xlabel("z", fontsize=14, font = "CAMBRIA", style = 'italic')
ax.set_ylabel("$M_{\odot}$", fontsize=14)
ax.set_title("Average Stellar Mass at Redshift z \n", fontsize=14)
plt.show()



# Close FITS file so it won't use up excess memory
hdu.close()
