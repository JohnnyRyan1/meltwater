#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute radiative forcing 

"""

# Import modules
import xarray as xr
import numpy as np
import pandas as pd
import netCDF4

#%%

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

# Define year
year = str(2019)

# Define time period
start = 30+29
end = 30+31+5

# Import MODIS
modis = xr.open_dataset(path1 + 'mod10a1-albedo-' + year + '.nc')

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')

# Define MERRA files
merra = xr.open_dataset(path1 + 'allwave-t2m-downscaled.nc')
merra_climatology = merra['swd_allsky'].mean(axis=2).values
merra_climatology[merra_climatology == 0] = np.nan

# Define mask
mask = ismip_1km['GIMP'].values

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)

# Define ice threshold
i = 55
 
# Some preprocessing
albedo = modis['albedo'].values.astype(np.float32)
albedo[mask3d == 0] = np.nan
albedo[albedo == 0] = np.nan
        
# Add max snow albedo
albedo[albedo > 84] = 84
albedo[albedo <= 30] = 30

# Classify
classified = np.copy(albedo)[:,:,start:end]
classified[classified <= i] = 1
classified[classified > i] = 2

#%%

#######################################################################
# Observed albedo
#######################################################################
observed = np.copy(albedo)[:,:,start:end]

# Compute SWnet
swnet_mean = (1 - (observed.mean(axis=2) / 100)) * merra_climatology

# Mask ice sheet
swnet_mean[mask == 0] = np.nan
       
#######################################################################
# Fix snow albedo
#######################################################################

# Fix snow albedo to fresh snow
fix_snow = np.copy(albedo)[:,:,start:end]   

fix_snow[classified == 2] = 84
fix_snow[classified == 1] = fix_snow[classified == 1]

# Compute SWnet
swnet_fixed_snow_mean = (1 - (fix_snow.mean(axis=2) / 100)) * merra_climatology
                       
# Mask ice sheet
swnet_fixed_snow_mean[mask == 0] = np.nan

#######################################################################
# Fix glacier ice albedo
#######################################################################

# Fix glacier ice albedo to the ice threshold
fix_ice = np.copy(albedo)[:,:,start:end]
fix_ice[classified == 1] = i
fix_ice[classified == 2] = fix_ice[classified == 2]

# Compute SWnet
swnet_fixed_ice_mean = (1 - (fix_ice.mean(axis=2) / 100)) * merra_climatology
                        
# Mask ice sheet
swnet_fixed_ice_mean[mask == 0] = np.nan

#######################################################################
# Fix snow albedo to that 84 and glacier ice albedo to 55
#######################################################################

# Fix both
fix_both = np.copy(albedo)[:,:,start:end]
fix_both[classified == 1] = i
fix_both[classified == 2] = 84

# Compute SWnet
swnet_fixed_both_mean = (1 - (fix_both.mean(axis=2) / 100)) * merra_climatology
                        
# Mask ice sheet
swnet_fixed_both_mean[mask == 0] = np.nan

#######################################################################
# Fix all the ice sheet to 84
#######################################################################

# Fix all with snow and glacier ice albedo observed on June 1
fix_all = np.copy(albedo)[:,:,start:end]
fix_all[classified > 0] = 84

# Compute SWnet
swnet_fixed_all_mean = (1 - (fix_all.mean(axis=2) / 100)) * merra_climatology

# Mask ice sheet
swnet_fixed_all_mean[mask == 0] = np.nan

# Compute radiative forcing
ice_diff = swnet_mean - swnet_fixed_ice_mean
snow_diff = swnet_mean - swnet_fixed_snow_mean
snowline_diff = swnet_fixed_both_mean - swnet_fixed_all_mean
bulk_diff = swnet_mean - swnet_fixed_all_mean

#%%
###############################################################################
# Make DataFrame
###############################################################################

# Make DataFrame
df = pd.DataFrame(list(zip(exp1, exp2, exp3, exp4, exp5)))

df.columns = ['fixed_snow_all', 'observed_albedo', 'fixed_snow_ice', 
              'fixed_ice', 'fixed_snow']

# Divide by area
df = df / np.sum(mask == 1)

# Melt-albedo feedbacks increase SWnet by...
total_albedo = (df['observed_albedo'] - df['fixed_snow_all']) / df['fixed_snow_all']
print(total_albedo.mean())

# Compute radiative forcing in W m-2

# SWnet due to reducing glacier ice albedo 
df['ice_forcing'] = df['observed_albedo'] - df['fixed_ice']

# SWnet due to reducing snow albedo
df['snow_forcing'] = df['observed_albedo'] - df['fixed_snow']

# SWnet due to snowline fluctuations
df['snowline_forcing'] = df['fixed_snow_ice'] - df['fixed_snow_all']

# Save as csv
df.to_csv(path + 'positive-forcing-results.csv')

# Remove first layer
snowline_diff = snowline_diff[:,:,1:]
ice_diff = ice_diff[:,:,1:]
snow_diff = snow_diff[:,:,1:]
bulk_diff = bulk_diff[:,:,1:]
bulk_swnet = bulk_swnet[:,:,1:]

# Save grids as NetCDF
lats, lons = modis['latitude'].values, modis['longitude'].values
    













