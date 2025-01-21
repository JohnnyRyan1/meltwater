#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate radiative feedback of surface water by elevation.

"""

# Import packages
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

# Read raster data
sw_2018 = xr.open_dataset(path1 + 'zhang/surface_water_mask_2018_1km.tif')
sw_2019 = xr.open_dataset(path1 + 'zhang/surface_water_mask_2019_1km.tif')

# Set NaNs to zero
sw_2018 = sw_2018.fillna(0)
sw_2019 = sw_2019.fillna(0)

# Convert to array
sw_2018 = sw_2018['band_data'].values[0,:,:]
sw_2019 = sw_2019['band_data'].values[0,:,:]

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

# Define MERRA files
merra = xr.open_dataset(path1 + 'allwave-t2m-downscaled.nc')
merra_climatology = merra['swd_allsky'].mean(axis=2).values
merra_climatology[merra_climatology == 0] = np.nan

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

# Import data
water_albedo_2018 = pd.read_csv(path1 + 'albedo-effects-2018.csv')
water_albedo_2019 = pd.read_csv(path1 + 'albedo-effects-2019.csv')

# Combine datasets
water_albedo = pd.concat([water_albedo_2018, water_albedo_2019])

# Water effect
water_effect, effect_uncert = np.abs(water_albedo['slope'].mean()), water_albedo['slope'].std()

# Other radiative forcing processes
other = xr.open_dataset(path1 + 'final-surface-forcing-grids.nc')

#%%

# Compute albedo change
albedo_change_2018 = sw_2018 * water_effect
albedo_change_2019 = sw_2019 * water_effect

albedo_change_2018_min = sw_2018 * (water_effect - effect_uncert)
albedo_change_2019_min = sw_2019 * (water_effect - effect_uncert)

albedo_change_2018_max = sw_2018 * (water_effect + effect_uncert)
albedo_change_2019_max = sw_2019 * (water_effect + effect_uncert)

# Compute radiative forcing
swnet_mean_2018 = albedo_change_2018 * merra_climatology
swnet_mean_2019 = albedo_change_2019 * merra_climatology

swnet_mean_2018_min = albedo_change_2018_min * merra_climatology
swnet_mean_2019_min = albedo_change_2019_min * merra_climatology

swnet_mean_2018_max = albedo_change_2018_max * merra_climatology
swnet_mean_2019_max = albedo_change_2019_max * merra_climatology

#%%

area, area_water_2018, area_water_2019 = [], [], []
forcing_2018, min_2018, max_2018 = [], [], []
forcing_2019, min_2019, max_2019 = [], [], []
temp_2018, temp_2019 = [], []

for e in range(len(elevations) - 1):
    print(e)
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    forcing_2018.append(np.nanmean(swnet_mean_2018[elevation_mask]))
    min_2018.append(np.nanmean(swnet_mean_2018_min[elevation_mask]))
    max_2018.append(np.nanmean(swnet_mean_2018_max[elevation_mask]))

    forcing_2019.append(np.nanmean(swnet_mean_2019[elevation_mask]))
    min_2019.append(np.nanmean(swnet_mean_2019_min[elevation_mask]))
    max_2019.append(np.nanmean(swnet_mean_2019_max[elevation_mask]))

    area.append(elevation_mask.sum())
    area_water_2018.append(np.nansum(sw_2018[elevation_mask]))
    area_water_2019.append(np.nansum(sw_2019[elevation_mask]))
    temp_2018.append(np.nanmean(merra['t2m'].values[:,:,16][elevation_mask]))
    temp_2019.append(np.nanmean(merra['t2m'].values[:,:,17][elevation_mask]))

area = np.array(area)
area_water_2018 = np.array(area_water_2018)
area_water_2019 = np.array(area_water_2019)

forcing_2018 = np.array(forcing_2018)
min_2018 = np.array(min_2018)
max_2018 = np.array(max_2018)
temp_2018 = np.array(temp_2018)

forcing_2019 = np.array(forcing_2019)
min_2019 = np.array(min_2019)
max_2019 = np.array(max_2019)
temp_2019 = np.array(temp_2019)

# Save as DataFrame
df = pd.DataFrame(list(zip(forcing_2018, min_2018, max_2018,
                           forcing_2019, min_2019, max_2019,
                           area_water_2018, area_water_2019, 
                           temp_2018, temp_2019, area)))
df.columns = ['forcing_2018', 'min_2018', 'max_2018', 
              'forcing_2019', 'min_2019', 'max_2019',
              'water_2018', 'water_2019', 'temp_2018', 'temp_2019', 'area']

# Compute energy in PJ (power * area * time)
df['energy_pj_2018'] = (df['forcing_2018'] * (df['area'] * 1000000) * 86400) / 1e15
df['energy_pj_2019'] = (df['forcing_2019'] * (df['area'] * 1000000) * 86400) / 1e15

df['min_energy_pj_2018'] = (df['min_2018'] * (df['area'] * 1000000) * 86400) / 1e15
df['min_energy_pj_2019'] = (df['min_2019'] * (df['area'] * 1000000) * 86400) / 1e15

df['max_energy_pj_2018'] = (df['max_2018'] * (df['area'] * 1000000) * 86400) / 1e15
df['max_energy_pj_2019'] = (df['max_2019'] * (df['area'] * 1000000) * 86400) / 1e15

#%%


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(8, 8),
                                             layout='constrained')

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.barh(range(len(area_water_2019)), area_water_2019, align='edge',  alpha=0.4, color=c1, edgecolor='k')
ax1.barh(range(len(area_water_2018)), area_water_2018, align='edge',  alpha=0.4, color=c2, edgecolor='k')
ax1.set_ylim(0,17)
ax1.grid(linestyle='dotted', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=13)
ax1.yaxis.set_ticks(np.arange(0, len(area_water_2018), 2))
ax1.set_yticklabels(elevations[:-1][::2])

ax2.plot(temp_2019, elevations[:-1], color=c1, zorder=2, lw=2, 
         alpha=0.8, label='')
ax2.plot(temp_2018, elevations[:-1], color=c2, zorder=2, lw=2, 
         alpha=0.8, label='')
#ax2.set_xlim(0, 6)
ax2.yaxis.set_ticks(elevations[:-1][::2])
ax2.set_yticklabels([])

ax3.plot(forcing_2018, elevations[:-1], color=c2, zorder=2, lw=2, 
         alpha=0.8)
ax3.fill_betweenx(elevations[:-1],
                 min_2018,
                 max_2018,
                 zorder=1, color=c2, alpha=0.2)
ax3.plot(forcing_2019, elevations[:-1], color=c1, zorder=2, lw=2, 
         alpha=0.8)
ax3.fill_betweenx(elevations[:-1],
                 min_2019,
                 max_2019,
                 zorder=1, color=c1, alpha=0.2)
ax3.set_xlim(0, 2)
ax3.yaxis.set_ticks(elevations[:-1][::2])
ax3.set_yticklabels(elevations[:-1][::2])

ax4.plot(df['energy_pj_2018'], elevations[:-1], color=c2, zorder=2, lw=2, 
         alpha=0.8, label='')
ax4.plot(df['energy_pj_2019'], elevations[:-1], color=c1, zorder=2, lw=2, 
         alpha=0.8, label='')
ax4.fill_betweenx(elevations[:-1],
                 df['min_energy_pj_2018'],
                 df['max_energy_pj_2018'],
                 zorder=1, color=c2, alpha=0.2)
ax4.fill_betweenx(elevations[:-1],
                 df['min_energy_pj_2019'],
                 df['max_energy_pj_2019'],
                 zorder=1, color=c1, alpha=0.2)
ax4.set_xlim(0, 6)
ax4.yaxis.set_ticks(elevations[:-1][::2])
ax4.set_yticklabels([])

ax1.set_xlabel('Surface water area (km$^2$)', fontsize=14)
ax2.set_xlabel('Summer air temperature (K)', fontsize=14)
ax3.set_xlabel('Radiative forcing (W m$^{-2}$)', fontsize=14)
ax4.set_xlabel('Radiative forcing (PJ d$^{-1}$)', fontsize=14)

ax1.set_ylabel('Elevation (m a.s.l.)', fontsize=14)

for ax in [ax2, ax3, ax4]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.set_ylim(0, 3400)

ax1.text(0.03, 0.91, "a", fontsize=20, transform=ax1.transAxes, zorder=2)
ax2.text(0.03, 0.91, "b", fontsize=20, transform=ax2.transAxes)
ax3.text(0.03, 0.91, "c", fontsize=20, transform=ax3.transAxes, zorder=1)
ax4.text(0.03, 0.91, "d", fontsize=20, transform=ax4.transAxes, zorder=1)


#%%


"""
P1
"""

# Mean radiative forcing between 0-1600 m
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values <= 1600)
print(np.nanmean(swnet_mean_2018[elevation_mask]))
print(np.nanmean(swnet_mean_2019[elevation_mask]))

#%%

"""
P2
"""

# Energy absorbed by ice sheet due to meltwater ponding
print(df['energy_pj_2018'].sum())
print(df['energy_pj_2019'].sum())

#%%
"""
P3

"""

swnet_diff = np.nanmean(swnet_mean_2019) - np.nanmean(swnet_mean_2018)
temp_diff = np.nanmean(merra['t2m'].values[:,:,17]) - np.nanmean(merra['t2m'].values[:,:,16])
energy_2018 = np.nansum(((swnet_mean_2018 * 1000000) * 86400) / 1e15)
energy_2019 = np.nansum(((swnet_mean_2019 * 1000000) * 86400) / 1e15)

snowline_2018 = np.nansum(((other['snowline'][:,:,16].values.astype(float) * 1000000) * 86400) / 1e15)
snowline_2019 = np.nansum(((other['snowline'][:,:,17].values.astype(float) * 1000000) * 86400) / 1e15)

snow_2018 = np.nansum(((other['snow'][:,:,16].values.astype(float) * 1000000) * 86400) / 1e15)
snow_2019 = np.nansum(((other['snow'][:,:,17].values.astype(float) * 1000000) * 86400) / 1e15)

ice_2018 = np.nansum(((other['ice'][:,:,16].values.astype(float) * 1000000) * 86400) / 1e15)
ice_2019 = np.nansum(((other['ice'][:,:,17].values.astype(float) * 1000000) * 86400) / 1e15)

# Compute feedback
df['temp_change'] = df['temp_2019'] - df['temp_2018']
df['water_feedback'] = (df['water_2019'] - df['water_2018']) / df['temp_change']
df['forcing_feedback'] = (df['forcing_2019'] - df['forcing_2018']) / df['temp_change']
df['energy_feedback'] =  (df['energy_pj_2019'] - df['energy_pj_2018']) / df['temp_change']

print('Summer air temperatures were %0.2f K warmer in 2019 than 2018' % (np.nanmean(merra['t2m'].values[:,:,17]) - np.nanmean(merra['t2m'].values[:,:,16])))
print('Extent of water increases by %0.2f km2 K-1' % np.nansum(df['water_feedback']))

print('Radiative forcing increases by %0.2f W m-2 K-1' % (swnet_diff / temp_diff))
print('Radiative forcing increases by %0.2f PJ d-1 K-1' % ((energy_2019 - energy_2018) / temp_diff))

#%%

# What would quadruple the amount of water look like?
sw_2019[~mask] = np.nan
sw_2019_double = np.minimum(sw_2019 * 4, 1)
albedo_change_double_2019 = sw_2019_double * water_effect
swnet_double_2019 = albedo_change_double_2019 * merra_climatology
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values <= 1600)
print(np.nanmean(swnet_double_2019[elevation_mask]))
energy_double_2019 = np.nansum(((swnet_double_2019 * 1000000) * 86400) / 1e15)
print(energy_double_2019)

sw_2018[~mask] = np.nan
sw_2018_double = np.minimum(sw_2018 * 4, 1)
albedo_change_double_2018 = sw_2018_double * water_effect
swnet_double_2018 = albedo_change_double_2018 * merra_climatology
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values <= 1600)
print(np.nanmean(swnet_double_2018[elevation_mask]))
energy_double_2018 = np.nansum(((swnet_double_2018 * 1000000) * 86400) / 1e15)
print(energy_double_2018)

#%%

water_double = []

for e in range(len(elevations) - 1):
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    water_double.append(np.nansum(sw_2019_double[elevation_mask]))

water_double = np.array(water_double)

# Save as DataFrame
double_df = pd.DataFrame(list(zip(water_double)))
double_df.columns = ['water']


#%%












