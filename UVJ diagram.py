# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 17:18:28 2023

@author: Rasmus
"""

import pandas as pd #imports pandas
import numpy as np #imports numpy
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

vj = (23.9 - 2.5 * np.log10(v_average)) - (23.9 - 2.5 * np.log10(j_average)) #finds v-j in mag (x-axis)
uv = (23.9 - 2.5 * np.log10(u_average)) - (23.9 - 2.5 * np.log10(v_average)) #finds u-v in mag (y-axis)

plt.figure(figsize=(4, 4))
plt.title("UVJ diagram") #defines the title
plt.xlabel("V - J (AB mag)") #defines x-axis label
plt.ylabel("U - V (AB mag)") #defines y-axis label
plt.scatter(vj,uv) #adds the points to the plot
plt.plot([-0.5, 0.693], [1.3, 1.3], color = "r") #draws the horizontal line
plt.plot([1.6, 1.6], [2.1, 4], color = "r") #draws the vertical line
plt.plot([0.693, 1.6], [1.3, 2.098], color = "r") #draws the tilted line
plt.xlim(0,4) #defines range of x-axis
plt.ylim(-0.5,3.5) #defines range of y-axis
plt.show() #shows the plot