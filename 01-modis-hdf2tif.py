#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

DESCRIPTION

1. Read and stack MOD10A1 HDFs

"""

# Import modules
import rioxarray as rio
import os
import glob

#%%

# Define path
path = '/Users/jr555/Documents/research/hydrology/modis/'

# Define location of MODIS data
modis_files = sorted(glob.glob(path + 'MOD*.hdf'))


modis_data = rio.open_rasterio(modis_files[0])

data_variable = modis_data["Snow_Albedo_Daily_Tile"]

# Reproject to EPSG:3413
data_array_reprojected = data_variable.rio.reproject("EPSG:3413")

data_array_reprojected.rio.to_raster(path + 'MOD10A1.A2015201.h16v02.061.2021328053529.tif')


#%%


# Define path
path = '/Users/jr555/Documents/research/hydrology/modis/'

# Define location of MODIS data
modis_files = sorted(glob.glob(path + 'MCD*.hdf'))


modis_data = rio.open_rasterio(modis_files[0])

data_variable = modis_data["Albedo_BSA_shortwave"]

# Reproject to EPSG:3413
data_array_reprojected = data_variable.rio.reproject("EPSG:3413")

data_array_reprojected.rio.to_raster(path + 'MCD43A3.A2015202.h16v02.061.2021329232024.tif')





