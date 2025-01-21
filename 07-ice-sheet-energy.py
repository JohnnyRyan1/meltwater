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


#%%

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'


# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

# Define MERRA files
merra = xr.open_dataset(path1 + 'allwave-t2m-downscaled.nc')

# Compute SW climatology
merra_climatology = merra['swd_allsky'].mean(axis=2).values
merra_climatology[merra_climatology == 0] = np.nan

# Define MODIS files
modis_2018 = xr.open_dataset(path1 + 'mod10a1-albedo-2018.nc')
modis_2019 = xr.open_dataset(path1 + 'mod10a1-albedo-2019.nc')

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

# Compute SWnet
swnet_2018 = (1 - (modis_2018 / 100)) * merra_climatology
swnet_2019 = (1 - (modis_2019 / 100)) * merra_climatology

# Compute energy in PJ (power * area * time)
swnet_energy_2018 = ((swnet_2018 * 1000000) * 86400) / 1e15
swnet_energy_2019 = ((swnet_2019 * 1000000) * 86400) / 1e15

print('Energy absorbed by ice sheet in 2018 was %.2f PJ d-1' % (np.nansum(swnet_energy_2018)))
print('Energy absorbed by ice sheet in 2019 was %.2f PJ d-1' % (np.nansum(swnet_energy_2019)))



#%%











