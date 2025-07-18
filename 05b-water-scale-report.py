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

# Define user
user = 'johnnyryan'

# Define path
path1 = '/Users/' + user + '/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'
path2 = '/Users/' + user + '/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/03-final-revision/'

# Define path to files
drone = gpd.read_file(path1 + 'drone/wq7/wq7-20150721-polys-within.shp')
russell = gpd.read_file(path1 + 'drone/russell/rus-20150712-polys.shp')
zhang_2018 = gpd.read_file(path1 + 'zhang/zhang_polys_2018_wq7.shp')
zhang_2019 = gpd.read_file(path1 + 'zhang/zhang_polys_2019_wq7.shp')

# Compute areas
drone['area'] = drone.area
russell['area'] = russell.area
zhang_2018['area'] = zhang_2018.area
zhang_2019['area'] = zhang_2019.area

# Import raster for total area calculation
classified = rio.open_rasterio(path1 + 'drone/wq7/wq7-20150721-classified-within.tif')

# Count non-NaN values across all dimensions
non_nan_count_total = classified.count()
study_area = non_nan_count_total * classified.rio.transform()[0] * classified.rio.transform()[0]
study_area_km2 = study_area / 1000000

# Area of Russell/Isunguata
russell_total = gpd.read_file(path1 + 'drone/wq7/wq7-outline.shp')

1721 + 1460 + 1745 + 1683 + 1442 + 1681

#%%

print('Area of water in WQ7 in 2018 from Zhang et al. = %.2f km2' %np.sum(zhang_2018['area'] / 1000000))
print('Area of water in WQ7 in 2019 from Zhang et al. = %.2f km2' %np.sum(zhang_2019['area'] / 1000000))

print('Area of water in WQ7 in 2015 from Ryan et al. = %.2f km2' %np.sum(drone['area'] / 1000000))


print('Fraction of water in WQ7 in 2018 from Zhang et al. = %.2f %%' %((np.sum(zhang_2018['area'] / 1000000)/study_area_km2)*100))
print('Fraction of water in WQ7 in 2019 from Zhang et al. = %.2f %%' %((np.sum(zhang_2019['area'] / 1000000)/study_area_km2)*100))

print('Fraction of water in WQ7 in 2015 from Ryan et al. = %.2f %%' %((np.sum(drone['area'] / 1000000)/study_area_km2)*100))
print('Fraction of water at Russell in 2015 from Ryan et al. = %.3f %%' %((np.sum(russell['area'])/russell_total.area.iloc[0])))

#%%

# Prepare the histogram
bins = np.logspace(np.log10(drone['area'].min()), np.log10(drone['area'].max()), num=30)  # Log-spaced bins for power-law data
hist1, bin_edges = np.histogram(drone['area'], bins=bins)  # Use density=True for probability density
hist2, bin_edges = np.histogram(zhang_2019['area'], bins=bins)  # Use density=True for probability density
hist3, bin_edges = np.histogram(zhang_2018['area'], bins=bins)  # Use density=True for probability density
hist4, bin_edges = np.histogram(russell['area'], bins=bins)

# Calculate bin centers
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

#%%
# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

# Plot the area-frequency plot
fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(6,4), 
                                    layout='constrained')
ax1.fill_between(bin_centers, hist1, step="post", alpha=0.5, color=c4, label="Drone (2015)")
ax1.fill_between(bin_centers, hist2, step="post", alpha=0.5, color=c1, label="Satellite (2019)")
ax1.fill_between(bin_centers, hist3, step="post", alpha=0.5, color=c2, label="Satellite (2018)")

ax1.axvline(x=100, ls='dashed', lw=1.5, color='k')
ax1.axvline(x=1000, ls='dashed', lw=1.5, color='k')

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("Area (m$^2$)", fontsize=12)
ax1.set_ylabel("Frequency", fontsize=12)
ax1.legend(fontsize=11)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0, alpha=0.5)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.set_xlim(0.2, 300000)
plt.savefig(path2 + 'fig-5-meltwater-histograms.pdf')

#%%

print('Minimum feature in Zhang et al. is %.0f m2' % np.min(zhang_2019['area']))

fraction_less_100 = np.sum(drone[drone['area'] < 100]['area']) / np.sum(drone['area'])
print('%.1f %% of water area consists of features that are less than 100 m2 in area' % (fraction_less_100*100))

fraction_less_1000 = np.sum(drone[drone['area'] < 1000]['area']) / np.sum(drone['area'])
print('%.1f %% of water area consists of features that are less than 1000 m2 in area' % (fraction_less_1000*100))

fraction_100_1000 = np.sum(drone[(drone['area'] > 100) & (drone['area'] < 1000)]['area']) / np.sum(drone['area'])
print('%.1f %% of water area consists of features that are < 1000 m2 and > 100 m2 in area' % (fraction_100_1000*100))


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




















