#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Bunch of scripts used to reproject the 200 m surface waters from 
Chudley et al. (2022) to 1 km.

"""

# Import packages
import numpy as np
import rioxarray
from rasterio.enums import Resampling

#%%

# Define user
user = 'johnnyryan'

# Define path
path1 = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
path2 = '/Users/'  + user + '/Dropbox (University of Oregon)/research/hydrology/data/'

year = '2019'

in_path = '/Users/johnnyryan/Documents/' 
input_filename = year + '_water_fraction_200m.tif'

out_path = '/Users/johnnyryan/Documents/tiles-v2/' 
output_filename = 'tile_'

tile_size_x = 10000
tile_size_y = 10000

ds = rioxarray.open_rasterio(in_path + input_filename, masked=True)
xsize = ds.shape[2]
ysize = ds.shape[1]

# Import ISMIP 1 km grid
ismip_1km = rioxarray.open_rasterio(path1 + '1km-ISMIP6-GIMP.nc', variable='GIMP')
ismip_1km = ismip_1km.rio.write_crs('EPSG:3413')

# Get ISMIP6 lat lons
x_1km, y_1km = np.meshgrid(ismip_1km['x'].values, ismip_1km['y'].values)

#%%

# Downscale from 200 m to 1 km
scale_factor = 0.2
    
# Define new height and width
new_width = int(ds.rio.width * scale_factor)
new_height = int(ds.rio.height * scale_factor)

# Downsample
downsample = ds.rio.reproject(
    ds.rio.crs,
    shape=(new_height, new_width),
    resampling=Resampling.average,
)

# Save
downsample.astype('float').rio.to_raster(in_path + 'reproject-v2/' + input_filename + '.tif', tiled=True)


#%%

# Resample to match extent of the ISMIP 1 km grid
ds = rioxarray.open_rasterio(in_path + 'reproject-v2/' + input_filename + '.tif', masked=True)

#%%

ds_match = ds.rio.reproject_match(ismip_1km)
ds_match = ds_match.where(ds_match != 1.7976931348623157e+308)
ds_match.rio.to_raster(in_path + year + '_water_fraction_1km.nc', tiled=True)
ds_match.rio.to_raster(in_path + year + '_water_fraction_1km.tif', tiled=True)

#%%










