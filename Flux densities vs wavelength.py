# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 17:43:21 2023

@author: Rasmus
"""

import pandas as pd #imports pandas
import matplotlib.pyplot as plt #imports matplotlib

df = pd.read_fwf("seds.txt", header = None) #creates a dataframe from the data provided

u_central = 3709 * 10 ** -4 #defines a constant in micrometer
u_width = 518 * 10 ** -4/2 #defines a constant in micrometer
v_central = 5487 * 10 ** -4 #defines a constant in micrometer
v_width = 954 * 10 ** -4/2 #defines a constant in micrometer
j_central = 12525 * 10 ** -4 #defines a constant in micrometer
j_width = 1718 * 10 ** -4/2 #defines a constant in micrometer

u_band = df.loc[(df[0] <= u_central + u_width) & (df[0] >= u_central - u_width)] #defines the u-band range
v_band = df.loc[(df[0] <= v_central + v_width) & (df[0] >= v_central - v_width)] #defines the v-band range
j_band = df.loc[(df[0] <= j_central + j_width) & (df[0] >= j_central - j_width)] #defines the j-band range

u_average = u_band.mean()[1:len(df)] #finds average flux densities and removes wavelengths
v_average = v_band.mean()[1:len(df)] #finds average flux densities and removes wavelengths
j_average = j_band.mean()[1:len(df)] #finds average flux densities and removes wavelengths

df.columns = ["Wavelength (µm)","G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12"] #column names

#plots
ax = df.plot(x = "Wavelength (µm)", xlim = (10**-1, 10**1), y = 
        ["G1", "G2", "G3", "G4","G5","G6","G7","G8","G9","G10","G11","G12"], 
        ylim = (10**-6, 10**2), ylabel = "Flux density (W/m2)", lw = 0.5)

ax.set_yscale('log')
ax.set_xscale('log')
plt.title("Flux Density against Wavelength")
plt.grid()

ux = pd.Series(u_central, index = range(12))
uy = u_average
plt.plot(ux, uy, marker = ".", markersize = 5, markeredgecolor = "black",
         markerfacecolor = "black", ls = "None")

vx = pd.Series(v_central, index = range(12))
vy = v_average
plt.plot(vx, vy, marker = ".", markersize = 5, markeredgecolor = "black",
         markerfacecolor = "black", ls = "None")

jx = pd.Series(j_central, index = range(12))
jy = j_average
plt.plot(jx, jy, marker = ".", markersize = 5, markeredgecolor = "black",
         markerfacecolor = "black", ls = "None")

plt.show()