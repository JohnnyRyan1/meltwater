#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Make footprints from WQ7 ortho.

"""

# Import packages
import rioxarray as rio
from geocube.vector import vectorize

#%%

# Define path
path = '/Users/jr555/Documents/research/hydrology/drone/20150721/'

# Open the raster file
image = rio.open_rasterio(path + 'wq7-20150721-ortho.tif')

# Get one band
target_band = image.isel(band=0)
   
# Create a mask of non-zero values
mask = target_band.where(target_band != 0).astype('uint8')

# Vectorize
gdf = vectorize(mask)

# Filter 
gdf = gdf[gdf['_data'] == 1]

# Define projection
gdf.crs = "EPSG:3413"

# Save to file
gdf.to_file(path + 'Outlines/wq7-20150721-ortho.shp')


#%%




























