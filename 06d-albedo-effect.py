#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate albedo lowering effect of water.

"""

# Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Define user
user = 'johnnyryan'

# Define path
path1 = '/Users/' + user + '/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

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
colormap = plt.get_cmap('viridis', water.shape[1])
v = '#440154FF'

# Plot the area-frequency plot
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12,4), 
                                    layout='constrained')


ax1.barh(range(len(area)), area['water_area'], align='edge', alpha=0.5, color=v, edgecolor='k')
ax1.set_ylim(0,17)
ax1.grid(linestyle='dotted', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=13)
ax1.set_yticks(np.arange(0, len(area), 2))
ax1.set_yticklabels(elevations[:-1][::2])

for i in range(non.shape[0]):
    ax2.plot(1 - (non.iloc[i, :] / water.iloc[i,:]), elevations[:-1],
             color=colormap(i), lw=2, alpha=0.5, zorder=0)
#ax2.scatter(0.15, 1000, s=100, zorder=1, color=v)
ax2.set_xlim(0, 0.6)
ax2.set_ylim(0, 3200)

ax1.set_xlabel("Water area (km$^2$)", fontsize=14)
ax1.set_ylabel("Elevation", fontsize=14)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax1.tick_params(axis='both', which='major', labelsize=12)

ax2.set_xlabel("Albedo variability due to meltwater ponding", fontsize=14)
ax2.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax2.tick_params(axis='both', which='major', labelsize=12)
ax2.set_yticklabels([])

from matplotlib.lines import Line2D

custom_lines = [Line2D([0], [0], color=colormap(0.), lw=2),
                Line2D([0], [0], color=colormap(.5), lw=2),
                Line2D([0], [0], color=colormap(1.), lw=2)]

ax2.legend(custom_lines, ['1 km', '25 km', '100 km'], fontsize=14)

ax1.text(0.03, 0.89, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.03, 0.87, "b", fontsize=24, transform=ax2.transAxes)

#%%

"""
Note:
    Ryan et al. (2018) found that meltwater ponding accounted for 15% of albedo
    variability at 1000 m and length scale 25 km.
    
"""
print('We find %.1f %%' %(((1 - (non.iloc[10, :] / water.iloc[10,:]))[5])*100))



#%%

# Import data
water_albedo_2018 = pd.read_csv(path1 + 'albedo-effects-2018.csv')
non_water_albedo_2018 = pd.read_csv(path1 + 'non-water-albedo-effect-2018.csv')
water_albedo_2019 = pd.read_csv(path1 + 'albedo-effects-2019.csv')
non_water_albedo_2019 = pd.read_csv(path1 + 'non-water-albedo-effect-2019.csv')

# Combine datasets
water_albedo = pd.concat([water_albedo_2018, water_albedo_2019])
non_water_albedo = pd.concat([non_water_albedo_2018, non_water_albedo_2019])

#%%

print("Albedo variability in 5x5 windows with water is +/- %.4f" % (water_albedo['std'].mean()))

print("Albedo variability in 5x5 windows without water is +/- %.4f" % (non_water_albedo['std'].mean()))

print("Local albedo variability is %.1f higher in 5x5 windows with water" % (water_albedo['std'].mean() / non_water_albedo['std'].mean()))

print('Meltwater ponding is highly correlated with albedo with a mean correlation coefficient (r) of %.2f' %(water_albedo['r'].mean()))

print('Mean albedo effect of meltwater ponding is  %.4f +/- %.4f' % (water_albedo['slope'].mean(), water_albedo['slope'].std()))

print('About %.1f%% of local albedo variability is likely not caused by meltwater ponding' % ((non_water_albedo['std'].mean() / water_albedo['std'].mean()) * 100))

print('It is not clear whether there this source of albedo variability increases or \
decreases our albedo effects but likely captured by the std. dev.')

#%%

# Mean albedo effect by elevation











