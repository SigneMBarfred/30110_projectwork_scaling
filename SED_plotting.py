# -*- coding: utf-8 -*-
"""
Rasmus & Signe - SED and UVJ for 12 galaxy templates
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plottools


df = pd.read_csv('seds_1.txt', delim_whitespace=True, header=None)
df.columns = ["lambda", "G1", "G2", "G3", "G4","G5","G6","G7","G8","G9","G10","G11","G12"]

#evaluate all wavelengths in ångstrøm not microns
df['lambda'] = df['lambda']*(10**4)

#replace all nonexistent fluxes with not a number
df = df.replace(0, np.nan)

#averages on u, v and j band#define bands by cosmos 2020 paper
u = [3709]*(len(df.columns)-1)
u_int = df[(3709-(518/2) <= df['lambda']) & (df['lambda'] <= 3709+(518/2))]
u_ave = u_int.mean()[1:]

v = [5487]*(len(df.columns)-1)
v_int = df[(5487-(954/2) <= df['lambda']) & (df['lambda'] <= 5487+(954/2))]
v_ave = v_int.mean()[1:]

j = [12525]*(len(df.columns)-1)
j_int = df[(12525-(1718/2) <= df['lambda']) & (df['lambda'] <= 12525+(1718/2))]
j_ave = j_int.mean()[1:]

#convert averages to magn
df_fluxes = df.iloc[:,1:13] #select flux columns
df_magn = 23.9 - 2.5*np.log10(df_fluxes) #convert flux to AB-magn
#df_magn.replace([np.inf,-np.inf],0,inplace=True) #wait we cant do this
df_magn.columns = ["G1", "G2", "G3", "G4","G5","G6","G7","G8","G9","G10","G11","G12"]


#plot w calc averages in magn as well - does it fit?
u_ave_ab = 23.9 - 2.5*np.log10(u_ave)
v_ave_ab = 23.9 - 2.5*np.log10(v_ave)
j_ave_ab = 23.9 - 2.5*np.log10(j_ave)
#kunne gøres smartere ved at putte averages i liste
c = np.arange(0,12,1) 


result = pd.concat([df.iloc[:,0], df_magn], axis = 1)

ax = result.plot(x = "lambda", xlim = (10**2.9, 10**5), y = 
        ["G1", "G2", "G3", "G4","G5","G6","G7","G8","G9","G10","G11","G12"], 
        ylim = (20, 38), ylabel = "AB magnitude", lw = 0.5,cmap='viridis') #plots


plt.scatter(u, u_ave_ab.values,c=c,cmap='viridis',zorder=100) #plot averages for the 3 bands with different colours
plt.scatter(v, v_ave_ab.values,c=c,cmap='viridis',zorder=101)
plt.scatter(j, j_ave_ab.values,c=c,cmap='viridis',zorder=102)
ax_zoom = plottools.zoom_axes(ax,[10**3.2,10**3.3],[22,32],[10**4,10**6],[30,34])
ax_zoom.plot(x = "lambda", xlim = (10**2.9, 10**5), y = 
        ["G1", "G2", "G3", "G4","G5","G6","G7","G8","G9","G10","G11","G12"], 
        ylim = (20, 38), ylabel = "AB magnitude", lw = 0.5,cmap='viridis')
ax.invert_yaxis()
ax.set_xscale('log')
plt.title("AB magnitude against Wavelength")
plt.grid()
plt.show()

#uvj diagram
vj = (23.9 - 2.5 * np.log10(v_ave)) - (23.9 - 2.5 * np.log10(j_ave)) #finds v-j in mag (x-axis)
uv = (23.9 - 2.5 * np.log10(u_ave)) - (23.9 - 2.5 * np.log10(v_ave)) #finds u-v in mag (y-axis)
df_uvj = pd.concat([vj,uv],axis=1) #gather uvj data to one dataframe
classes = ["G1", "G2", "G3", "G4","G5","G6","G7","G8","G9","G10","G11","G12"]
fig, ax3 = plt.subplots() #make uvj plot
p = plt.scatter(vj,uv,c=c,cmap='viridis') 
plt.title("UVJ diagram") #defines the title
plt.xlabel("V - J (AB mag)") #defines x-axis label
plt.ylabel("U - V (AB mag)") #defines y-axis label
plt.plot([-0.5, 0.6931818182], [1.3, 1.3], color = "maroon") #draws the horizontal line
plt.plot([1.6, 1.6], [2.098, 4], color = "maroon") #draws the vertical line
plt.plot([0.6931818182, 1.6], [1.3, 2.098], color = "maroon") #draws the tilted line
plt.xlim(-0.5,4) #defines range of x-axis
plt.ylim(-0.5,4) #defines range of y-axis
plt.legend(handles=p.legend_elements()[0], labels=classes)
plt.show() #shows the plot