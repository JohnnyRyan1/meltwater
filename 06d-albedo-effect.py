#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate albedo lowering effect of water.

"""

# Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

# Import data
water = pd.read_csv(path1 + 'scales/water.csv')
non = pd.read_csv(path1 + 'scales/non.csv')

# Define windows
windows = np.array((1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90,
                    100))

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

# Get 600-1600 m [3:8]
water_mid = water.iloc[:, 3:8].mean(axis=1)
non_mid = non.iloc[:, 3:8].mean(axis=1)

# By elevation
1 - (non.iloc[10, :] / water.iloc[10,:])

#%%

# Report stats
1 - (non_mid[2:5].mean() / water_mid[2:5].mean())

1 - (non_mid[15:].mean() / water_mid[15:].mean())



#%%

"""
P1
"""

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

# Plot the area-frequency plot
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12,4), 
                                    layout='constrained')


ax1.plot(windows[:-1], 1- (non_mid / water_mid),
         color=c1, lw=2, label='Elevation = 600-1600 m')
ax1.set_xlim(0, 90)
#ax1.plot(windows[:-1], water['std'] / water['std'].max(), label='All grid cells',
#         color=c2, lw=2)

#ax1.set_ylim(0, 0.7)

ax2.plot(1 - (non.iloc[10, :] / water.iloc[10,:]), elevations,
         color=c2, lw=2, label='Length scale = 25 km')

ax1.set_xlabel("Length scale (km)", fontsize=12)
ax1.set_ylabel("Albedo variability due to \n meltwater ponding", fontsize=12)
ax1.legend(fontsize=12, loc=1)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax1.tick_params(axis='both', which='major', labelsize=12)

ax2.set_xlabel("Albedo variability due to meltwater ponding", fontsize=12)
ax2.set_ylabel("Elevation", fontsize=12)
ax2.legend(fontsize=12, loc=1)
ax2.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax2.tick_params(axis='both', which='major', labelsize=12)

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











