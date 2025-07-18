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
path = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/zhang/'

# Open the raster file
image = rio.open_rasterio(path + 'surface_water_mask_2019_wq7.tif')

# Get one band
target_band = image.isel(band=0)
   
# Vectorize
gdf = vectorize(target_band)

# Filter 
gdf = gdf[gdf['_data'] == 1]

# Define projection
gdf.crs = "EPSG:3413"

# Save to file
gdf.to_file(path + 'zhang_polys_2019_wq7.shp')


#%%




























