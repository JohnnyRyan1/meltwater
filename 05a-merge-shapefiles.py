#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Size distribution of surface meltwater features.

"""

# Import packages
import geopandas as gpd
import glob

# Define path to files
path = '/Users/jr555/Documents/research/hydrology/drone/20150721/Polys/'

# Define files
files = sorted(glob.glob(path + '*.shp'))

#%%

# Initialize an empty list to hold the GeoDataFrames
geometries = []

# Loop through each shapefile and read it into a GeoDataFrame
for shp in files:
    gdf = gpd.read_file(shp)
    geometries.extend(list(gdf.geometry))

    
#%%

gdf = gpd.GeoDataFrame(geometries, columns=['geometry'], geometry='geometry',
                       crs='EPSG:3413')

gdf.to_file('/Users/jr555/Documents/research/hydrology/drone/20150721/wq7-20150721-polys.shp')

#%%