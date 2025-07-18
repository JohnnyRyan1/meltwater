#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Make footprints from orthos.

"""

# Import packages
import rioxarray as rio
from geocube.vector import vectorize
import glob
import os

#%%

# Define path
path = '/Volumes/Extreme SSD/hydrology/data/'

# Define orthomosaics
orthos = sorted(glob.glob(path + 'drone/*/Orthomosaics/*.tif'))

def raster_to_polygon(raster_file, infileshortname, date):
    
    # Open the raster file
    image = rio.open_rasterio(raster_file)
    
    # Get one band
    target_band = image.isel(band=0)
       
    # Create the mask (e.g., mask for non-zero values)
    mask = (target_band < 255).astype("uint8")
    
    # Vectorize
    gdf = vectorize(mask)
    
    # Filter 
    gdf = gdf[gdf['_data'] == 1]
    
    # Define projection
    gdf.crs = "EPSG:3413"
    
    # Save to file
    gdf.to_file(path + 'outlines/' + date + '-' + infileshortname[8:] + '.shp')


for ortho in orthos:
    
    # Get the path and filename separately
    infilepath, infilename = os.path.split(ortho)
    
    # Get the short name (filename without extension)
    infileshortname, extension = os.path.splitext(infilename)
    
    # Get grandparent folder
    date = os.path.basename(os.path.dirname(os.path.dirname(ortho)))
    
    if os.path.exists(path + 'outlines/' + date + '-' + infileshortname[8:] + '.shp'):
        pass
    else:
        raster_to_polygon(ortho, infileshortname, date)


#%%




























