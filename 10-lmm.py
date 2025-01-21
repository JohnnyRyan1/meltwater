#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Linear mixing model to theorize relationship between surface water and albedo

"""

# Import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

#%%

# Define albedo values
water = 0.15
ice = 0.55
snow = 0.84

fraction_water = np.arange(0, 1, 0.01)

ice_albedo, snow_albedo = [], []
for i in fraction_water:
    fraction_ice = 1 - i
    ice_albedo.append((fraction_ice * ice) + (i * (ice-water)))
    snow_albedo.append((fraction_ice * snow) + (i * (snow-water)))

slope1, intercept, r, p, se = linregress(fraction_water, ice_albedo)
slope2, intercept, r, p, se = linregress(fraction_water, snow_albedo)

#%%

fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(6,4), 
                                    layout='constrained')

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.plot(fraction_water, ice_albedo, zorder=2, color=c2, lw=2,
            label='water on ice')
ax1.plot(fraction_water, snow_albedo, zorder=2, color=c3, lw=2,
            label='water on snow')
ax1.legend(fontsize=12)
ax1.set_xlabel('Fraction of surface water', fontsize=12)
ax1.set_ylabel('Albedo (unitless)', fontsize=12)
ax1.grid(linestyle='dotted', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=12)

#%%







