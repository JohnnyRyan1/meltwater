#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Report albedo variability caused by meltwater ponding at different length scales.

"""

# Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Define user
user = 'jr555'

# Define path
path1 = '/Users/' + user + '/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'
path2 = '/Users/' + user + '/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/figures/'

# Import data
water = pd.read_csv(path1 + 'scales/water-2000.csv')
non = pd.read_csv(path1 + 'scales/non-water-2000.csv')
area = pd.read_csv(path1 + 'scales/water-area.csv')

# Define windows
windows = np.array((1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90,
                    100))

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

# Get 600-1600 m [3:8]
water_mid = water.iloc[:, 3:8].mean(axis=1)
non_mid = non.iloc[:, 3:8].mean(axis=1)

#%%

"""
P1
"""
# Define colour map
colormap = plt.get_cmap('coolwarm_r', water.shape[1])
v = '#E05861'

# Plot the area-frequency plot
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10,4), 
                                    layout='constrained')

for i in range(non.shape[0]):
    ax1.plot(1 - (non.iloc[i, :] / water.iloc[i,:]), elevations[:-1],
             color=colormap(i), lw=2, alpha=0.5, zorder=0)
ax1.set_xlim(0, 0.6)
ax1.set_ylim(0, 3200)

ax1.set_xlabel("Albedo variability due to meltwater ponding", fontsize=13)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.set_ylabel("Elevation", fontsize=13)
ax1.set_yticks(elevations[:-1][::2])
ax1.set_yticklabels(elevations[:-1][::2])

# Plot legend
custom_lines = [Line2D([0], [0], color=colormap(0.), lw=2),
                Line2D([0], [0], color=colormap(.5), lw=2),
                Line2D([0], [0], color=colormap(1.), lw=2)]
ax1.legend(custom_lines, ['1 km', '25 km', '100 km'], fontsize=12)

ax2.barh(range(len(area)), area['water_area'], align='edge', alpha=0.5, color=v, edgecolor='k')
ax2.set_ylim(0,16)
ax2.grid(linestyle='dotted', lw=1, zorder=1)
ax2.tick_params(axis='both', which='major', labelsize=12)
ax2.set_yticks(np.arange(0, len(area), 2))
ax2.set_yticklabels([])
ax2.set_xlabel("Meltwater extent (km$^2$)", fontsize=13)
ax2.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)


ax1.text(0.03, 0.89, "a", fontsize=20, transform=ax1.transAxes)
ax2.text(0.03, 0.87, "b", fontsize=20, transform=ax2.transAxes)

#plt.savefig(path2 + 'fig1-variability.png', dpi=300)

#%%

"""
Note:
    Ryan et al. (2018) found that meltwater ponding accounted for 15% of albedo
    variability at 1000 m and length scale 25 km.
    
    Length scale of 25 km is index 10
    
"""

print('Meltwater can account for up to %.1f %% of albedo variability' %(np.max(1 - (non / water))*100))

print('At length scales of 25 km, we find %.1f %%' %(((1 - (non.iloc[10, :] / water.iloc[10,:])).iloc[5])*100))

print('At lower elevations, we find %.1f %%' %((1 - (non / water)).iloc[10, 0:4].mean()*100))
 
print('At higher elevations, we find %.1f %%' %((1 - (non / water)).iloc[10, 11:].mean()*100))

#%%











