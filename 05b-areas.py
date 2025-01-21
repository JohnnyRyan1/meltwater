#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Areas of meltwater in WQ7.

"""

# Import packages
import rioxarray as rio
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt

# Define path to files
drone = gpd.read_file('/Users/jr555/Documents/research/hydrology/drone/20150721/wq7-20150721-polys-within.shp')
zhang_2018 = gpd.read_file('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/zhang/zhang_polys_2018_wq7.shp')
zhang_2019 = gpd.read_file('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/zhang/zhang_polys_2019_wq7.shp')
outline = gpd.read_file('/Users/jr555/Documents/research/hydrology/drone/20150721/Outlines/wq7-20150721-outline.shp')

# Compute areas
drone['area'] = drone.area
zhang_2018['area'] = zhang_2018.area
zhang_2019['area'] = zhang_2019.area

# Import raster for total area calculation
classified = rio.open_rasterio('/Users/jr555/Documents/research/hydrology/drone/20150721/wq7-20150721-classified-within.tif')

# Count non-NaN values across all dimensions
non_nan_count_total = classified.count()
study_area = non_nan_count_total * classified.rio.transform()[0] * classified.rio.transform()[0]
study_area_km2 = study_area / 1000000

#%%

print('Area of water in WQ7 in 2018 from Zhang et al. = %.2f km2' %np.sum(zhang_2018['area'] / 1000000))
print('Area of water in WQ7 in 2019 from Zhang et al. = %.2f km2' %np.sum(zhang_2019['area'] / 1000000))

print('Area of water in WQ7 in 2015 from Ryan et al. = %.2f km2' %np.sum(drone['area'] / 1000000))

#%%

# Prepare the histogram
bins = np.logspace(np.log10(drone['area'].min()), np.log10(drone['area'].max()), num=30)  # Log-spaced bins for power-law data
hist1, bin_edges = np.histogram(drone['area'], bins=bins)  # Use density=True for probability density
hist2, bin_edges = np.histogram(zhang_2019['area'], bins=bins)  # Use density=True for probability density
hist3, bin_edges = np.histogram(zhang_2018['area'], bins=bins)  # Use density=True for probability density

# Calculate bin centers
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

#%%
# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

# Plot the area-frequency plot
fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(7,5), 
                                    layout='constrained')
ax1.fill_between(bin_centers, hist1, step="post", alpha=0.5, color=c4, label="This study")
ax1.fill_between(bin_centers, hist2, step="post", alpha=0.5, color=c1, label="Zhang (2019)")
ax1.fill_between(bin_centers, hist3, step="post", alpha=0.5, color=c2, label="Zhang (2018)")

ax1.axvline(x=100, ls='dashed', lw=1.5, color='k')

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("Area (m$^2$)", fontsize=12)
ax1.set_ylabel("Frequency", fontsize=12)
ax1.legend(fontsize=12)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax1.tick_params(axis='both', which='major', labelsize=12)

#%%

print('Minimum feature in Zhang et al. is %.0f m2' % np.min(zhang_2019['area']))

fraction_less_100 = np.sum(drone[drone['area'] < 100]['area']) / np.sum(drone['area'])
print('%.1f %% of water area consists of features that are less than 100 m2 in area' % (fraction_less_100*100))

fraction_less_100 = np.sum(drone[(drone['area'] > 100) & (drone['area'] < 1000)]['area']) / np.sum(drone['area'])
print('%.1f %% of water area consists of features that are < 1000 m2 and > 100 m2 in area' % (fraction_less_100*100))


drone_1000 = np.sum(drone[drone['area'] > 100]['area'])
zhang_1000 = np.sum(zhang_2019[zhang_2019['area'] > 100]['area'])

print('Zhang classifies %.2f %% of water > 100 m2' % ((zhang_1000/drone_1000)*100))

drone_less_1000 = np.sum(drone[drone['area'] < 100]['area'])
zhang_less_1000 = np.sum(zhang_2019[zhang_2019['area'] < 100]['area'])

print('Zhang classifies %.2f %% of water < 100 m2' % ((zhang_less_1000/drone_less_1000)*100))

smallest = bins[np.argmax(hist1)]
print('The smallest reliably derived features are %.2f m2' % smallest)


#%%



#%%




















