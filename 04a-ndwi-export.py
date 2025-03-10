#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export NDWI for Russell Glacier

"""

# Import packages
import glob
import os
import rioxarray as rio

#%%

# Define path
path = '/Users/jr555/Documents/research/hydrology/drone/isunguata/'

# Define training shapefiles
ortho_files = sorted(glob.glob(path +  'Orthomosaics/*.tif*'))


#%%

for file in ortho_files:
    
    # Get the path and filename separately
    infilepath, infilename = os.path.split(file)
    
    # Get the short name (filename without extension)
    infileshortname, extension = os.path.splitext(infilename)
    
    # Import raster
    image = rio.open_rasterio(file)

    # Convert to three-band image
    image = image.isel(band=[0,1,2]).astype(float)
    
    # Make NDWI image
    ndwi = (image[2, :, :] - image[0,:, :]) / \
            (image[2, :, :] + image[0,:, :])
    ndwi.rio.to_raster(path + 'ndwi/' + infileshortname + '.tif')
    
#%%