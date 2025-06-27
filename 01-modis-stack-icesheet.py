#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

DESCRIPTION

1. Read and stack MOD10A1 or MCD43A3 HDFs

"""

# Import modules
import rioxarray as rio
import xarray as xr
import numpy as np
import os
import glob

#%%

# Define year
year = str(2019)

# Define path
path = '/Users/jr555/Documents/research/hydrology/modis/' + year + '/'

# Define location of MODIS data
modis_files = sorted(glob.glob(path + 'MOD*.hdf'))

# Define tiles
tiles = ['15v01', '15v02', '16v00', '16v01', '16v02', '17v00', '17v01', '17v02']

# Import surface water raster for matching
sw = rio.open_rasterio('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/zhang/surface_water_mask_' + year + '_500m.tif')

#%%
mean_albedo_tiles = []
for tile in tiles:
    
    modis_files_tile = []
    for file in modis_files:
        
        # Get the path and filename separately
        infilepath, infilename = os.path.split(file)
        
        # Get the short name (filename without extension)
        infileshortname, extension = os.path.splitext(infilename)
        
        if infileshortname[18:23] == tile:
            modis_files_tile.append(file)
    
    # Produce mean albedo
    data_arrays = []
    for file in modis_files_tile:
        # Read
        modis_data = rio.open_rasterio(file)
        
        # Get snow albedo
        data_arrays.append(modis_data["Snow_Albedo_Daily_Tile"])

    # Stack the DataArrays along a new dimension (e.g., 'time')
    stacked_da = xr.concat(data_arrays, dim="time")    
    
    # Set values greater than 100 to NaN
    stacked_da = stacked_da.where(stacked_da <= 100, np.nan)
    
    # Set values less than 20 to NaN
    stacked_da = stacked_da.where(stacked_da >= 20, np.nan)
    
    # Scale
    stacked_da = stacked_da*0.01
    
    # Compute nanmean along the time dimension
    mean_albedo = stacked_da.mean(dim="time", skipna=True)
    
    # Append
    mean_albedo_tiles.append(mean_albedo)
    
#%%

# Merge the DataArrays into one
merged_da = xr.combine_by_coords(mean_albedo_tiles)

# Reproject to EPSG:3413
merged_da = merged_da.rio.reproject("EPSG:3413")

# Match to surface water extent 500m
merged_da_match = merged_da.rio.reproject_match(sw)

# Export as GeoTiff
merged_da_match['Snow_Albedo_Daily_Tile'].rio.to_raster('/Users/jr555/Documents/research/hydrology/modis/mod-' + year +'.tif')

#%%

# Define location of MODIS data
modis_files = sorted(glob.glob(path + 'MCD*.hdf'))


#%%
mean_albedo_tiles = []
for tile in tiles:
    
    modis_files_tile = []
    for file in modis_files:
        
        # Get the path and filename separately
        infilepath, infilename = os.path.split(file)
        
        # Get the short name (filename without extension)
        infileshortname, extension = os.path.splitext(infilename)
        
        if infileshortname[18:23] == tile:
            modis_files_tile.append(file)
    
    # Produce mean albedo
    data_arrays = []
    for file in modis_files_tile:
        # Read
        modis_data = rio.open_rasterio(file)
        
        # Get snow albedo
        data_arrays.append(modis_data["Albedo_BSA_shortwave"])

    # Stack the DataArrays along a new dimension (e.g., 'time')
    stacked_da = xr.concat(data_arrays, dim="time")    
    
    # Set values equal to 32767 to NaN
    stacked_da = stacked_da.where(stacked_da != 32767, np.nan)
    
    # Scale
    stacked_da = stacked_da*0.001
    
    # Set values greater than 1 to NaN
    stacked_da = stacked_da.where(stacked_da <= 1, np.nan)
    
    # Set values less than 0.2 to NaN
    stacked_da = stacked_da.where(stacked_da >= 0.2, np.nan)
    
    # Compute nanmean along the time dimension
    mean_albedo = stacked_da.mean(dim="time", skipna=True)
    
    # Append
    mean_albedo_tiles.append(mean_albedo)

#%%
# Merge the DataArrays into one
merged_da = xr.combine_by_coords(mean_albedo_tiles)

# Reproject to EPSG:3413
merged_da = merged_da.rio.reproject("EPSG:3413")

# Match to surface water extent 500m
merged_da_match = merged_da.rio.reproject_match(sw)
    
# Export as GeoTiff
merged_da_match['Albedo_BSA_shortwave'].rio.to_raster('/Users/jr555/Documents/research/hydrology/modis/mcd-' + year + '.tif')

#%%








