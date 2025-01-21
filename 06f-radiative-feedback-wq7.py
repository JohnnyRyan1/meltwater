#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute radiative feedback over WQ7

"""

# Import packages
import xarray as xr
import rioxarray as rio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

# Read raster data
drone = rio.open_rasterio(path1 + 'wq7-20150721-classified-1km.tif')
zhang = rio.open_rasterio(path1 + 'zhang/surface_water_mask_2018_1km.tif')

# Match projections
drone_match = drone.rio.reproject_match(zhang)

# Convert to array
drone_match = drone_match.values[0,:,:]
zhang = zhang.values[0,:,:]

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

# Define MERRA files
merra = xr.open_dataset(path1 + 'allwave-t2m-downscaled.nc')
merra_climatology = merra['swd_allsky'].mean(axis=2).values
merra_climatology[merra_climatology == 0] = np.nan

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
albedo_drone = drone_match * water_effect
albedo_zhang = zhang * water_effect

albedo_drone_min = drone_match * (water_effect - effect_uncert)
albedo_zhang_min = zhang * (water_effect - effect_uncert)

albedo_drone_max = drone_match * (water_effect + effect_uncert)
albedo_zhang_max = zhang * (water_effect + effect_uncert)

# Compute radiative forcing
swnet_drone = albedo_drone * merra_climatology
swnet_zhang = albedo_zhang * merra_climatology

swnet_drone_min = albedo_drone_min * merra_climatology
swnet_zhang_min = albedo_zhang_min * merra_climatology

swnet_drone_max = albedo_drone_max * merra_climatology
swnet_zhang_max = albedo_zhang_max * merra_climatology

#%%

print('Area of water from drone = %.2f in 2015' % (np.nansum(drone_match)))
print('Area of water from Zhang = %.2f in 2018' % (np.nansum(zhang[np.isfinite(drone_match)])))

print('We find %.1f more water coverage frome the drone imagery than Zhang' % (np.nansum(drone_match) / np.nansum(zhang[np.isfinite(drone_match)])))

forcing_drone = np.mean(swnet_drone[np.isfinite(swnet_drone)])
forcing_zhang = np.mean(swnet_zhang[np.isfinite(swnet_drone)])

print('Radiative forcing caused by melting ponding is %.2f W m-2 in drone' % (forcing_drone))
print('Radiative forcing caused by melting ponding is %.2f W m-2 in Zhang' % (forcing_zhang))

#%%


