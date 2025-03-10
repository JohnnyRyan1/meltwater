#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute radiative feedback over WQ7

"""

# Import packages
import xarray as xr
import rioxarray as rio
import numpy as np
import matplotlib.pyplot as plt

# Define mean albedo lowering effect and uncertainty
water_effect = 0.11
effect_uncert = 0.06

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

# Read raster data
drone = rio.open_rasterio(path1 + 'wq7-20150721-classified-1km.tif')
zhang_2018 = rio.open_rasterio(path1 + 'zhang/surface_water_mask_2018_1km.tif')
zhang_2019 = rio.open_rasterio(path1 + 'zhang/surface_water_mask_2019_1km.tif')

# Match projections
drone_match = drone.rio.reproject_match(zhang_2018)

# Convert to array
drone_match = drone_match.values[0,:,:]
zhang_2018 = zhang_2018.values[0,:,:]
zhang_2019 = zhang_2019.values[0,:,:]

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

# Define MERRA files
merra = xr.open_dataset(path1 + 'allwave-t2m-downscaled.nc')
merra_climatology = merra['swd_allsky'].mean(axis=2).values
merra_climatology[merra_climatology == 0] = np.nan

# Other radiative forcing processes
other = xr.open_dataset(path1 + 'final-surface-forcing-grids.nc')

# Import MODIS
mcd_2015 = rio.open_rasterio('/Users/jr555/Documents/research/hydrology/modis/2015/MCD43A3.A2015202.h16v02.061.2021329232024.tif')
mcd_2015_match = mcd_2015.rio.reproject_match(ismip_1km['SRF'].rio.write_crs('EPSG:3413'))
mcd_2015_match = mcd_2015_match.values[0,:,:]
mcd_2015_match = mcd_2015_match.astype(float) / 1000
mcd_2015_match[mcd_2015_match == 32.767] = np.nan

# Regions
regions_file = xr.open_dataset(path1 + 'temp_albedo_summer_climatologies.nc')
regions = regions_file['regions']

#%%

# Compute albedo change
albedo_drone = drone_match * water_effect
albedo_zhang_2018 = zhang_2018 * water_effect
albedo_zhang_2019 = zhang_2019 * water_effect

albedo_drone_min = drone_match * (water_effect - effect_uncert)
albedo_zhang_min_2018 = zhang_2018 * (water_effect - effect_uncert)
albedo_zhang_min_2019 = zhang_2019 * (water_effect - effect_uncert)

albedo_drone_max = drone_match * (water_effect + effect_uncert)
albedo_zhang_max_2018 = zhang_2018 * (water_effect + effect_uncert)
albedo_zhang_max_2019 = zhang_2019 * (water_effect + effect_uncert)

# Compute radiative forcing
swnet_drone = albedo_drone * merra_climatology
swnet_zhang_2018 = albedo_zhang_2018 * merra_climatology
swnet_zhang_2019 = albedo_zhang_2019 * merra_climatology

swnet_drone_min = albedo_drone_min * merra_climatology
swnet_zhang_min_2018 = albedo_zhang_min_2018 * merra_climatology
swnet_zhang_min_2019 = albedo_zhang_min_2019 * merra_climatology

swnet_drone_max = albedo_drone_max * merra_climatology
swnet_zhang_max_2018 = albedo_zhang_max_2018 * merra_climatology
swnet_zhang_max_2019 = albedo_zhang_max_2019 * merra_climatology

fixed_ice = np.copy(mcd_2015_match)
fixed_ice[np.isfinite(fixed_ice)] = 0.55
swnet_fixed_ice = (1 - fixed_ice) * merra_climatology
swnet_ice = (1 - mcd_2015_match) * merra_climatology
swnet_ice_albedo = swnet_ice - swnet_fixed_ice

fixed_snow =  np.copy(mcd_2015_match)
fixed_snow[np.isfinite(fixed_snow)] = 0.84
swnet_fixed_snow = (1 - fixed_snow) * merra_climatology
swnet_snow = (1 - mcd_2015_match) * merra_climatology
swnet_snow_albedo = swnet_snow - swnet_fixed_snow

#%%

forcing_drone = np.mean(swnet_drone[np.isfinite(swnet_drone)])
forcing_drone_min = np.mean(swnet_drone_min[np.isfinite(swnet_drone)])
forcing_drone_max = np.mean(swnet_drone_max[np.isfinite(swnet_drone)])

forcing_zhang_2018 = np.mean(swnet_zhang_2018[np.isfinite(swnet_drone)])
forcing_zhang_min_2018 = np.mean(swnet_zhang_min_2018[np.isfinite(swnet_drone)])
forcing_zhang_max_2018 = np.mean(swnet_zhang_max_2018[np.isfinite(swnet_drone)])

forcing_zhang_2019 = np.mean(swnet_zhang_2019[np.isfinite(swnet_drone)])
forcing_zhang_min_2019 = np.mean(swnet_zhang_min_2019[np.isfinite(swnet_drone)])
forcing_zhang_max_2019 = np.mean(swnet_zhang_max_2019[np.isfinite(swnet_drone)])

print('Radiative forcing caused by melting ponding is %.2f W m-2 in drone' % (forcing_drone))
print('Radiative forcing caused by melting ponding is %.2f W m-2 in drone' % (forcing_drone_min))
print('Radiative forcing caused by melting ponding is %.2f W m-2 in drone' % (forcing_drone_max))

print('\n')

print('Radiative forcing caused by melting ponding is %.2f W m-2 in Zhang' % (forcing_zhang_2018))
print('Radiative forcing caused by melting ponding is %.2f W m-2 in Zhang' % (forcing_zhang_min_2018))
print('Radiative forcing caused by melting ponding is %.2f W m-2 in Zhang' % (forcing_zhang_max_2018))

print('\n')

print('Radiative forcing caused by melting ponding is %.2f W m-2 in Zhang' % (forcing_zhang_2019))
print('Radiative forcing caused by melting ponding is %.2f W m-2 in Zhang' % (forcing_zhang_min_2019))
print('Radiative forcing caused by melting ponding is %.2f W m-2 in Zhang' % (forcing_zhang_max_2019))

#%%

# Mask
ice_2019_wq7 = swnet_ice_albedo[np.isfinite(drone_match)]

print('Radiative forcing due to glacier ice albedo is %.2f W m-2' % (np.nanmean(ice_2019_wq7)))

elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 1600) &\
    (ismip_1km['SRF'].values < 1800) & (regions == 6)
    

print('Radiative forcing due snow albedo in the percolation zone is is %.2f W m-2' %(np.nanmean(swnet_snow_albedo[elevation_mask])))

print('Therefore %.2f %%' % (((forcing_drone*4) / (np.nanmean(swnet_snow_albedo[elevation_mask])))*100))










