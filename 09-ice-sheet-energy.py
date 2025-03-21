#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
How much shortwave energy does the ice sheet absorb per day during the summer?

"""

# Import packages
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rioxarray as rio

#%%

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

# Define MERRA files
merra = xr.open_dataset(path1 + 'allwave-t2m-downscaled.nc')

# Compute SW climatology
merra_climatology = merra['swd_allsky'].mean(axis=2).values
merra_climatology[merra_climatology == 0] = np.nan

# Define MODIS files
modis_2018 = xr.open_dataset(path1 + 'mod10a1-albedo-2018.nc')
modis_2019 = xr.open_dataset(path1 + 'mod10a1-albedo-2019.nc')

# Define MCD43A3 data for Zhang time periods
mcd_2018 = rio.open_rasterio(path1 + 'mcd-2018.tif')
mcd_2019 = rio.open_rasterio(path1 + 'mcd-2019.tif')

mod_2018 = rio.open_rasterio(path1 + 'mod-2018.tif')
mod_2019 = rio.open_rasterio(path1 + 'mod-2019.tif')

#%%

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].rio.write_crs('EPSG:3413').astype(int)
mask_match_2018 = mask.rio.reproject_match(mcd_2018)
mask_match_2019 = mask.rio.reproject_match(mcd_2019)

# Mask
mcd_2018 = xr.where(mask_match_2018, mcd_2018, np.nan)
mcd_2019 = xr.where(mask_match_2019, mcd_2019, np.nan)

# Some preprocessing
mcd_2018 = np.where(mcd_2018 == 0, np.nan, mcd_2018)
mcd_2018 = np.where(mcd_2018 <= 0.3, 0.3, mcd_2018)
mcd_2018 = np.where(mcd_2018 > 0.84, 0.84, mcd_2018)

mcd_2019 = np.where(mcd_2019 == 0, np.nan, mcd_2019)
mcd_2019 = np.where(mcd_2019 <= 0.3, 0.3, mcd_2019)
mcd_2019 = np.where(mcd_2019 > 0.84, 0.84, mcd_2019)

#%%
# Some preprocessing
modis_2018 = modis_2018['albedo'].values.astype(np.float32)
modis_2018 = np.nanmean(modis_2018, axis=2)
modis_2018[modis_2018 == 0] = np.nan
modis_2018[modis_2018 > 84] = 84
modis_2018[modis_2018 <= 30] = 30

modis_2019 = modis_2019['albedo'].values.astype(np.float32)
modis_2019 = np.nanmean(modis_2019, axis=2)
modis_2019[modis_2019 == 0] = np.nan
modis_2019[modis_2019 > 84] = 84
modis_2019[modis_2019 <= 30] = 30

#%%

# Import data
water_albedo_2018 = pd.read_csv(path1 + 'albedo-effects-2018.csv')
water_albedo_2019 = pd.read_csv(path1 + 'albedo-effects-2019.csv')

# Combine datasets
water_albedo = pd.concat([water_albedo_2018, water_albedo_2019])

# Water effect
water_effect, effect_uncert = np.abs(water_albedo['slope'].mean()), water_albedo['slope'].std()

#%%

# Compute albedo lowering due to meltwater ponding
print('Meltwater covered %.2f %% of the ice in the summer of 2018' %((4903/mask.sum())*100))
print('Meltwater covered %.2f %% of the ice in the summer of 2019' %((9988/mask.sum())*100))

print('MCD43A3 albedo averages %.2f during late July in 2018' %(np.nanmean(mcd_2018)))
print('MCD43A3 albedo averages %.2f during late July in 2019' %(np.nanmean(mcd_2019)))

print('Meltwater reduced ice sheet albedo by %.2f %% in late July 2018' % (((0.3*0.001)/0.78)*100))
print('Meltwater reduced ice sheet albedo by %.2f %% in late July 2019' % (((0.62*0.001)/0.72)*100))


# Compute SWnet
swnet_2018 = (1 - (modis_2018 / 100)) * merra_climatology
swnet_2019 = (1 - (modis_2019 / 100)) * merra_climatology

# Compute energy in PJ (power * area * time)
swnet_energy_2018 = ((swnet_2018 * 1000000) * 86400) / 1e15
swnet_energy_2019 = ((swnet_2019 * 1000000) * 86400) / 1e15

print('SW energy absorbed by ice sheet in 2018 was %.2f PJ d-1' % (np.nansum(swnet_energy_2018)))
print('SW energy absorbed by ice sheet in 2019 was %.2f PJ d-1' % (np.nansum(swnet_energy_2019)))



#%%











