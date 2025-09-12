#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Size distribution of surface meltwater features.

"""

# Import packages
import geopandas as gpd
import glob

# Define path to files
path = '/Users/jr555/Documents/research/hydrology/drone/'

# Define files
files1 = sorted(glob.glob(path + 'isunguata/Polys/*.shp'))
files2 = sorted(glob.glob(path + 'russell/Polys/*.shp'))

files = files1+files2

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

gdf.to_file(path + 'russell/merged-polys.shp')

#%%

# Read
shp = gpd.read_file(path + 'russell/merged-polys.shp')

# Union
shp_union = shp.union_all()

# Put back into GeoDataFrame
gdf = gpd.GeoDataFrame(geometry=[shp_union])

# Explode
gdf_explode = gdf.explode()

# Set CRS
gdf_explode = gdf_explode.set_crs('EPSG:3413')

# Save
gdf_explode.to_file(path + 'russell/merged-polys-union.shp')
