#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Report albedo lowering effect of water.

"""

# Import packages
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Define user 
user = 'johnnyryan'

# Define path
path1 = '/Users/' + user + '/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'
path2 = '/Users/' + user + '/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/figures/'

# List DataFrames 
df_2018_list = sorted(glob.glob(path1 + 'effects/albedo-effects-2018-*.csv'))
df_2019_list = sorted(glob.glob(path1 + 'effects/albedo-effects-2019-*.csv'))

r_values = []
slope_values = []
for f in range(len(df_2018_list)):
    # Read both DataFrames
    df1 = pd.read_csv(df_2018_list[f])
    df2 = pd.read_csv(df_2019_list[f])
    
    # Concatenate
    df = pd.concat([df1, df2])
    
    # Append
    r_values.append(df['r'].values)
    slope_values.append(df['slope'].values)

# Define window size from 1 to 100 km
windows = np.array((2, 3, 4, 5, 6, 7, 8, 9, 10))

# Find mean r and slope
mean_r = []
std_r = []
mean_slope = []
std_slope = []
for i in range(len(r_values)):
    mean_r.append(r_values[i].mean())
    std_r.append(r_values[i].std())
    mean_slope.append(slope_values[i].mean())
    std_slope.append(slope_values[i].std())
    
#%%
idx = np.array(mean_r).argmin()

print('Highest correlation coefficients are %.2f' % (mean_r[3]))

print('Mean correlation coefficients is %.2f' % (np.mean(mean_r[0:9])))

print('Standard deviation of correlation coefficients is %.2f' % (np.mean(std_r[0:9])))

print('Mean slope at length scales <10 km are %.2f' % (np.mean(mean_slope[0:9])))

print('Standard deviation of slopes at length scales <10 km is %.2f' % (np.mean(std_slope[0:9])))


print('Correlation coefficients at 25 km is %.2f' % (mean_r[9]))
print('Standard deviation of correlation coefficients is %.2f' % (std_r[9]))

print('Correlation coefficients at 50 km is %.2f' % (mean_r[10]))
print('Standard deviation of correlation coefficients is %.2f' % (std_r[10]))

print('Correlation coefficients at 100 km is %.2f' % (mean_r[11]))
print('Standard deviation of correlation coefficients is %.2f' % (std_r[11]))

#%%

# Plot showing correlation coefficient with length scale
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(7, 5),
                                             layout='constrained', sharex=True)

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

labels = windows
bp1 = ax1.boxplot(r_values[0:9], showfliers=False, patch_artist=True)
bp2 = ax2.boxplot(slope_values[0:9], showfliers=False, patch_artist=True)
ax2.set_ylim(-0.35, 0.01)
ax1.set_ylabel('Correlation coefficient (r)', fontsize=12)
ax2.set_ylabel('Slope', fontsize=12)
ax2.set_xlabel('Length scale (km)', fontsize=12)

# Apply single color to all boxes
for box in bp1['boxes']:
    box.set(facecolor=c2, edgecolor='black', linewidth=1.2, alpha=0.7)

for box in bp2['boxes']:
    box.set(facecolor=c1, edgecolor='black', linewidth=1.2, alpha=0.7)
    
# Customize median line appearance
for median in bp1['medians']:
    median.set_color('k')
    median.set_linewidth(1)

for median in bp2['medians']:
    median.set_color('k')
    median.set_linewidth(1)

for ax in [ax1, ax2]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=12)

ax2.set_xticks(range(1, len(labels) + 1))
ax2.set_xticklabels(labels)


ax1.text(0.02, 0.89, "a", fontsize=20, transform=ax1.transAxes)
ax2.text(0.02, 0.07, "b", fontsize=20, transform=ax2.transAxes)

plt.savefig(path2 + 'sfig-2-correlation-length.png', dpi=300)


#%%

# Combine into a single array
slopes = np.concatenate(slope_values[0:9])

c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

# Plot showing histogram to slopes
fig, ax1 = plt.subplots(figsize=(6,4), layout='constrained')

# Create histogram
n, bins, patches = ax1.hist(slopes, bins=65, edgecolor='black', zorder=1, 
                            color=c2, alpha=0.8)

# Add titles and labels
ax1.set_xlabel('Albedo lowering effect of meltwater ponding', fontsize=12)
ax1.set_ylabel('Frequency', fontsize=12)
ax1.set_xlim(-0.3, 0)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.axvline(-0.11, color='k', ls='dashed', lw=1.5, zorder=2)
ax1.axvline(-0.11+0.06, color='k', ls='dashed', lw=1, zorder=2)
ax1.axvline(-0.11-0.06, color='k', ls='dashed', lw=1, zorder=2)

# Show plot
plt.savefig(path2 + 'sfig-7-histogram-slopes.png', dpi=300)


#%%























