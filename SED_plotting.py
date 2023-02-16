# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

df = pd.read_csv('seds_1.txt', delim_whitespace=True, header=None)
df.columns = ["lambda", "G1", "G2", "G3", "G4","G5","G6","G7","G8","G9","G10","G11","G12"]
fig, ax = plt.subplots()


#flux - wavelength plot in microns
# ax.plot(df['lambda'],df['G1'], label='Line1')
# ax.plot(df['lambda'],df['G2'], label='Line2')
# ax.plot(df['lambda'],df['G3'], label='Line3')
# ax.plot(df['lambda'],df['G4'], label='Line4')
# tick_spacing = 10000
# ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
# plt.xlabel('wavelength (microns)')
# plt.ylabel('Flux')
# plt.title('Plot of SEDS')
# plt.legend()
# plt.show()


#evaluate all wavelengths in ångstrøm not microns
df['lambda'] = df['lambda']*(10**4)


ax.plot(df['lambda'],df['G1'], label='Line1')
ax.plot(df['lambda'],df['G2'], label='Line2')
ax.plot(df['lambda'],df['G3'], label='Line3')
ax.plot(df['lambda'],df['G4'], label='Line4')
tick_spacing = 10000
ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('wavelength (ångstrøm)')
plt.ylabel('Flux')
plt.title('Plot of SEDS')
plt.legend()
plt.show()

# log axes
#averages on u, v and j band
#convert averages to magn
#plot data in magn_ab
#plot w calc averages in magn as well - does it fit?

#print(seds_data)