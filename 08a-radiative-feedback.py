#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate radiative feedback of surface water by elevation.

"""

# Import packages
import xarray as xr
import rioxarray as rio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Define mean albedo lowering effect and uncertainty
water_effect = 0.11
effect_uncert = 0.06

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'
path2 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/figures/'

# Read raster data
sw_2018 = xr.open_dataset(path1 + 'zhang/surface_water_mask_2018_1km.tif')
sw_2019 = xr.open_dataset(path1 + 'zhang/surface_water_mask_2019_1km.tif')

# Regions
regions_file = xr.open_dataset(path1 + 'temp_albedo_summer_climatologies.nc')
regions = regions_file['regions']

# Set NaNs to zero
sw_2018 = sw_2018.fillna(0)
sw_2019 = sw_2019.fillna(0)

# Convert to array
sw_2018 = sw_2018['band_data'].values[0,:,:]
sw_2019 = sw_2019['band_data'].values[0,:,:]

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values
elev = ismip_1km['SRF'].rio.write_crs('EPSG:3413')

# Define MERRA files
merra = xr.open_dataset(path1 + 'allwave-t2m-downscaled.nc')
merra_climatology = merra['swd_allsky'].mean(axis=2).values
merra_climatology[merra_climatology == 0] = np.nan

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

# Other radiative forcing processes
other = xr.open_dataset(path1 + 'final-surface-forcing-grids.nc')

# Import MODIS albedo data
mcd_2018 = rio.open_rasterio(path1 + 'mcd-2018.tif')
mcd_2019 = rio.open_rasterio(path1 + 'mcd-2019.tif')
mcd_2018_match = mcd_2018.rio.reproject_match(elev)
mcd_2019_match = mcd_2019.rio.reproject_match(elev)

# Import data
water = pd.read_csv(path1 + 'scales/water-2000.csv')
non = pd.read_csv(path1 + 'scales/non-water-2000.csv')
total_water_area = pd.read_csv(path1 + 'scales/water-area.csv')

# Define windows
windows = np.array((1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90,
                    100))

# Get 600-1600 m [3:8]
water_mid = water.iloc[:, 3:8].mean(axis=1)
non_mid = non.iloc[:, 3:8].mean(axis=1)

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

# Compute feedback
df['temp_change'] = df['temp_2019'] - df['temp_2018']
df['water_feedback'] = (df['water_2019'] - df['water_2018']) / df['temp_change']
df['forcing_feedback'] = (df['forcing_2019'] - df['forcing_2018']) / df['temp_change']
df['energy_feedback'] =  (df['energy_pj_2019'] - df['energy_pj_2018']) / df['temp_change']
df['energy_feedback_min'] =  (df['min_energy_pj_2019'] - df['max_energy_pj_2018']) / df['temp_change']
df['energy_feedback_max'] =  (df['max_energy_pj_2019'] - df['min_energy_pj_2018']) / df['temp_change']

#%%

r_list = np.arange(1, 9)
forcing = []
forcing_min, forcing_max = [], []
albedo = []
for r in r_list:
    data_list = []
    data_list_min = []
    data_list_max = []
    albedo_list = []
    for e in range(len(elevations) - 1):
        elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) &\
            (ismip_1km['SRF'].values < elevations[e+1]) & (regions == r)
        data_list.append(np.nanmean(swnet_mean_2019[elevation_mask]))
        data_list_min.append(np.nanmean(swnet_mean_2019_min[elevation_mask]))
        data_list_max.append(np.nanmean(swnet_mean_2019_max[elevation_mask]))
        albedo_list.append(np.nanmean(mcd_2019_match.values[0,:,:][elevation_mask]))
    forcing.append(data_list)
    forcing_min.append(data_list_min)
    forcing_max.append(data_list_max)
    albedo.append(albedo_list)
    
df_regions = pd.DataFrame(forcing)
df_regions_min = pd.DataFrame(forcing_min)
df_regions_max = pd.DataFrame(forcing_max)
df_albedo = pd.DataFrame(albedo)

#%%


fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(12, 8),
                                             layout='constrained')

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

# Define colour map
colormap = plt.get_cmap('coolwarm_r', water.shape[1])
v = '#E05861'

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

ax2.barh(range(len(area_water_2019)), area_water_2019, align='edge', 
         alpha=0.4, color=c1, edgecolor='k', label='2019')
ax2.barh(range(len(area_water_2018)), area_water_2018, align='edge', 
         alpha=0.4, color=c2, edgecolor='k', label='2018')
ax2.set_ylim(0,17)
ax2.grid(linestyle='dotted', lw=1, zorder=1)
ax2.tick_params(axis='both', which='major', labelsize=13)
ax2.yaxis.set_ticks(np.arange(0, len(area_water_2018), 2))
ax2.set_yticklabels([])
ax2.legend(loc=1, fontsize=12)

ax3.plot(forcing_2019, elevations[:-1], color=c1, zorder=2, lw=2, 
         alpha=0.8, label='2019')
ax3.fill_betweenx(elevations[:-1],
                 min_2019,
                 max_2019,
                 zorder=1, color=c1, alpha=0.2)
ax3.plot(forcing_2018, elevations[:-1], color=c2, zorder=2, lw=2, 
         alpha=0.8, label='2018')
ax3.fill_betweenx(elevations[:-1],
                 min_2018,
                 max_2018,
                 zorder=1, color=c2, alpha=0.2)
ax3.set_xlim(0, 1.2)
ax3.yaxis.set_ticks(elevations[:-1][::2])
ax3.set_yticklabels(elevations[:-1][::2])
ax3.set_yticklabels([])
ax3.legend(loc=1, fontsize=12)

ax4.plot(df['energy_pj_2019'], elevations[:-1], color=c1, zorder=2, lw=2, 
         alpha=0.8, label='2019')
ax4.plot(df['energy_pj_2018'], elevations[:-1], color=c2, zorder=2, lw=2, 
         alpha=0.8, label='2018')
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
ax4.legend(loc=1, fontsize=12)

ax5.plot(df_regions.iloc[0, :], elevations[:-1],
             color=c4, lw=2, alpha=1, zorder=0, label='N')
ax5.fill_betweenx(elevations[:-1],
                 df_regions_min.iloc[0, :],
                 df_regions_max.iloc[0, :],
                 zorder=1, color=c4, alpha=0.2)
ax5.plot(df_regions.iloc[5, :], elevations[:-1],
             color=c3, lw=2, alpha=1, zorder=0, label='SW')
ax5.fill_betweenx(elevations[:-1],
                 df_regions_min.iloc[5, :],
                 df_regions_max.iloc[5, :],
                 zorder=1, color=c3, alpha=0.2)
ax5.axhline(y=1250, ls='dashed', color='k')
ax5.axhline(y=400, ls='dashed', color='k')
ax5.legend(loc=1, fontsize=12)
ax5.yaxis.set_ticks(elevations[:-1][::2])
ax5.set_yticklabels(elevations[:-1][::2])
ax5.set_xlim(0, 3)
ax5.set_yticklabels([])

ax6.plot(df['energy_feedback'], elevations[:-1], color=c4, zorder=2, lw=2, 
         alpha=0.8, label='')
ax6.set_xlim(0, 2.3)
ax6.fill_betweenx(elevations[:-1],
                 df['energy_feedback_min'],
                 df['energy_feedback_max'],
                 zorder=1, color=c4, alpha=0.2)
ax6.yaxis.set_ticks(elevations[:-1][::2])
ax6.set_yticklabels(elevations[:-1][::2])
ax6.set_yticklabels([])

ax1.set_xlabel('Albedo variability due to meltwater', fontsize=14)
ax2.set_xlabel('Meltwater extent (km$^{2}$)', fontsize=14)
ax3.set_xlabel('Radiative forcing (W m$^{-2}$)', fontsize=14)
ax4.set_xlabel('Radiative forcing (PJ d$^{-1}$)', fontsize=14)
ax5.set_xlabel('Radiative forcing (W m$^{-2}$)', fontsize=14)
ax6.set_xlabel('Radiative feedback (PJ d$^{-1}$ K$^{-1}$)', fontsize=14)

ax1.set_ylabel('Elevation (m a.s.l.)', fontsize=14)
ax4.set_ylabel('Elevation (m a.s.l.)', fontsize=14)


for ax in [ax3, ax4, ax5, ax6]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.set_ylim(0, 3400)

ax1.text(0.05, 0.91, "a", fontsize=20, transform=ax1.transAxes, zorder=2)
ax2.text(0.03, 0.91, "b", fontsize=20, transform=ax2.transAxes)
ax3.text(0.03, 0.91, "c", fontsize=20, transform=ax3.transAxes, zorder=1)
ax4.text(0.03, 0.91, "d", fontsize=20, transform=ax4.transAxes, zorder=1)
ax5.text(0.03, 0.91, "e", fontsize=20, transform=ax5.transAxes, zorder=1)
ax6.text(0.03, 0.91, "f", fontsize=20, transform=ax6.transAxes, zorder=1)

ax5.text(0.83, 0.14, "site 1", fontsize=12, transform=ax5.transAxes, zorder=1)
ax5.text(0.83, 0.39, "site 2", fontsize=12, transform=ax5.transAxes, zorder=1)

plt.savefig(path2 + 'fig-3-radiative-forcing.png', dpi=300)

#%%


"""
P1a
"""

# Mean radiative forcing between 0-1600 m
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 1000) & (ismip_1km['SRF'].values <= 2000)

print(np.nanmean(swnet_mean_2019[elevation_mask]))
print(np.nanmean(swnet_mean_2019_min[elevation_mask]))
print(np.nanmean(swnet_mean_2019_max[elevation_mask]))

print('\n')

elevation_mask = (mask == True) & (ismip_1km['SRF'].values <=600)
print(np.nanmean(swnet_mean_2019[elevation_mask]))
print(np.nanmean(swnet_mean_2019_min[elevation_mask]))
print(np.nanmean(swnet_mean_2019_max[elevation_mask]))

print('\n')

elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 2000)
print(np.nanmean(swnet_mean_2019[elevation_mask]))
print(np.nanmean(swnet_mean_2019_min[elevation_mask]))
print(np.nanmean(swnet_mean_2019_max[elevation_mask]))

print('\n')

print(np.nanmean(swnet_mean_2018[elevation_mask]))

#%%
"""
P1b - radiative forcing in the dark zone (Sector 6)

"""

elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 1000) &\
    (ismip_1km['SRF'].values < 1400) & (regions == 6)
print(np.nanmean(swnet_mean_2019[elevation_mask]))
print(np.nanmean(swnet_mean_2019_min[elevation_mask]))
print(np.nanmean(swnet_mean_2019_max[elevation_mask]))

print('\n')

elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 1600) &\
    (ismip_1km['SRF'].values <= 1800) & (regions == 6)
print(np.nanmean(swnet_mean_2019[elevation_mask]))
print(np.nanmean(swnet_mean_2019_min[elevation_mask]))
print(np.nanmean(swnet_mean_2019_max[elevation_mask]))

print('\n')

elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 200) &\
    (ismip_1km['SRF'].values < 400) & (regions == 1)
print(np.nanmean(swnet_mean_2019[elevation_mask]))
print(np.nanmean(swnet_mean_2019_min[elevation_mask]))
print(np.nanmean(swnet_mean_2019_max[elevation_mask]))

print('\n')

elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 1000) &\
    (ismip_1km['SRF'].values <= 1200) & (regions == 1)
print(np.nanmean(swnet_mean_2019[elevation_mask]))
print(np.nanmean(swnet_mean_2019_min[elevation_mask]))
print(np.nanmean(swnet_mean_2019_max[elevation_mask]))

#%%
print(df_regions.iloc[0,5])
print(df_regions_min.iloc[0,5])

print(df_regions.iloc[5,5:7].mean())
print(df_regions_min.iloc[5,5:7].mean())
print(df_regions.iloc[5,8])
print(df_regions_min.iloc[5,8])

print(df_albedo.iloc[5,5:7].mean())
print(df_albedo.iloc[5,8].mean())
print(df_albedo.iloc[6,7].mean())
print(df_albedo.iloc[1,5].mean())

print(df_regions.iloc[5,6])
print(df_regions.iloc[5,8])


#%%

"""
P2
"""

# Energy absorbed by ice sheet due to meltwater ponding

print(df['energy_pj_2019'].sum())
print(df['min_energy_pj_2019'].sum())
print(df['max_energy_pj_2019'].sum())

print(df['energy_pj_2018'].sum())
print(df['min_energy_pj_2018'].sum())
print(df['max_energy_pj_2018'].sum())

# Convert to melt potential
((df['energy_pj_2019'].sum() * 1e+12) / 333.55) * 1e-12
((df['energy_pj_2018'].sum() * 1e+12) / 333.55) * 1e-12

snowline_2018 = np.nansum(((other['snowline'][:,:,16].values.astype(float) * 1000000) * 86400) / 1e15)
snowline_2019 = np.nansum(((other['snowline'][:,:,17].values.astype(float) * 1000000) * 86400) / 1e15)

snow_2018 = np.nansum(((other['snow'][:,:,16].values.astype(float) * 1000000) * 86400) / 1e15)
snow_2019 = np.nansum(((other['snow'][:,:,17].values.astype(float) * 1000000) * 86400) / 1e15)

ice_2018 = np.nansum(((other['ice'][:,:,16].values.astype(float) * 1000000) * 86400) / 1e15)
ice_2019 = np.nansum(((other['ice'][:,:,17].values.astype(float) * 1000000) * 86400) / 1e15)

bulk = snowline_2019+snow_2019+ice_2019

meltwater_2019 = df['energy_pj_2019'].sum()

print('Bulk radiative forcing is %.1f PJ d-1' % (bulk))

print('Meltwater ponding accounted for %.1f %%' % ((meltwater_2019/bulk)*100))

print('Snowline contributed %.1f %%' % ((snowline_2019 / bulk)*100))
print('Snow albedo contributed %.1f %%' % ((snow_2019 / bulk)*100))
print('Glacier ice albedo contributed %.1f %%' % ((ice_2019 / bulk)*100))



#%%
"""
P3

"""

swnet_diff = np.nanmean(swnet_mean_2019) - np.nanmean(swnet_mean_2018)
temp_diff = np.nanmean(merra['t2m'].values[:,:,17]) - np.nanmean(merra['t2m'].values[:,:,16])
energy_2018 = np.nansum(((swnet_mean_2018 * 1000000) * 86400) / 1e15)
energy_2019 = np.nansum(((swnet_mean_2019 * 1000000) * 86400) / 1e15)

min_2018 = np.nansum(((swnet_mean_2018_min * 1000000) * 86400) / 1e15)
max_2018 = np.nansum(((swnet_mean_2018_max * 1000000) * 86400) / 1e15)

min_2019 = np.nansum(((swnet_mean_2019_min * 1000000) * 86400) / 1e15)
max_2019 = np.nansum(((swnet_mean_2019_max * 1000000) * 86400) / 1e15)

print('Summer air temperatures were %0.2f K warmer in 2019 than 2018' % (np.nanmean(merra['t2m'].values[:,:,17]) - np.nanmean(merra['t2m'].values[:,:,16])))
print('Extent of water increases by %0.2f km2 K-1' % np.nansum(df['water_feedback']))

elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 1000) & (ismip_1km['SRF'].values <= 1800)
feedback_1000 = np.nanmean((swnet_mean_2019-swnet_mean_2018)[elevation_mask])
uncert_1000 = np.nanmean((swnet_mean_2019_max-swnet_mean_2018_min)[elevation_mask])

print('Radiative forcing increases by %0.2f PJ d-1 K-1' % ((energy_2019 - energy_2018) / temp_diff))
print('Uncertainty in radiative forcing %0.2f PJ d-1 K-1' % ((np.sqrt(((max_2019 - energy_2019) / temp_diff)**2 + ((max_2018 - energy_2018) / temp_diff)**2))))

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












