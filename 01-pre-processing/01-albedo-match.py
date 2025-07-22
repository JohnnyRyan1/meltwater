#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Local albedo variability and water area.

"""

# Import packages
import rioxarray as rio
import geopandas as gpd
import numpy as np
from rasterio.enums import Resampling

#%%

# Define path
path1 = '/Users/jr555/Documents/research/hydrology/drone/20150721/'
path2 = '/Users/jr555/Documents/research/hydrology/modis/'

# Read raster data
wq7 = rio.open_rasterio(path1 + 'wq7-20150721-classified.tif')
mcd = rio.open_rasterio(path2 + 'MCD43A3.A2015202.h16v02.061.2021329232024.tif')
mod = rio.open_rasterio(path2 + 'MOD10A1.A2015201.h16v02.061.2021328053529.tif')

# Read outline
outline = gpd.read_file(path1 + 'Outlines/wq7-20150721-outline.shp')

#%%

# Do some value manipulation
wq7 = wq7.where(wq7 != 1, 3)
wq7 = wq7.where(wq7 != 0, 1)
wq7 = wq7.where(wq7 != 2, 0)

# Clip the raster using the shapefile
wq7_clip = wq7.rio.clip(outline.geometry, outline.crs)


# Set zeros to NaNs
wq7_clip = wq7_clip.where(wq7_clip != 0, np.nan)
wq7_clip = wq7_clip.where(wq7_clip != 1, 0)
wq7_clip = wq7_clip.where(wq7_clip != 3, 1)

# Check
wq7_clip.rio.to_raster(path1 + 'wq7-20150721-classified-within.tif')

#%%

# Downscale from 0.3 m to 500 m
scale_factor = 0.0006
    
# Define new height and width
new_width = int(wq7_clip.rio.width * scale_factor)
new_height = int(wq7_clip.rio.height * scale_factor)

# Downsample
downsample = wq7_clip.rio.reproject(
    wq7_clip.rio.crs,
    shape=(new_height, new_width),
    resampling=Resampling.average,
)

# Check
downsample.rio.to_raster(path1 + 'wq7-20150721-classified-500m.tif')

#%% 

# Clip MODIS
mcd_clip = mcd.rio.clip(outline.geometry, outline.crs)
mod_clip = mod.rio.clip(outline.geometry, outline.crs)

#%%

# Match
mcd_clip_match = mcd_clip.rio.reproject_match(downsample)
mod_clip_match = mod_clip.rio.reproject_match(downsample)

# Save
mcd_clip_match.rio.to_raster(path2 + 'MCD43A3.A2015202.h16v02.061.2021329232024_matched.tif')
mod_clip_match.rio.to_raster(path2 + 'MOD10A1.A2015201.h16v02.061.2021328053529_matched.tif')




#%%













