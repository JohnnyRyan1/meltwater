#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Bunch of scripts used to reproject the 10 m surface water map to 500 m or 1 km.

Note:
    Have to match the extents of the surface water and water-filled crevasse
    rasters before this works. In QGIS --> Right-click --> Export --> Save As
    --> Extent (Calculate from layer) --> surface_water_mask_201X.tif


Note: 
    Can convert the ISMIP6_1km.nc file to GeoTiff with the following code.
    
    import rioxarray as rio
    import xarray as xr
    ismip_1km = xr.open_dataset("1km-ISMIP6-GIMP.nc", decode_coords="all")  
    gimp_1km = ismip_1km["GIMP"]
    gimp_1km.rio.write_crs('EPSG:3413')
    gimp_1km.rio.to_raster("ISMIP6_1km.tif")

"""

# Import packages
import os
import numpy as np
import rioxarray as rio
import xarray as xr
import glob
from rasterio.enums import Resampling

#%%

# Define year
year = str(2018)

# Define resolution in m
resolution = 1000

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

in_path = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/zhang/' 
water_filename = 'surface_water_mask_' + year + '.tif'
crevasse_filename = 'water-filled_crevasses_' + year + '_expanded.tif'

path = '/Users/jr555/Documents/'
out_path1 = path + 'water_tiles/' 
out_path2 = path + 'crevasse_tiles/' 
output_filename = 'tile_'

tile_size_x = 10000
tile_size_y = 10000

water = rio.open_rasterio(in_path + water_filename, masked=True)
crevasse = rio.open_rasterio(in_path + crevasse_filename, masked=True)

# Import ISMIP 1 km grid
ismip_1km = rio.open_rasterio(path1 + '1km-ISMIP6-GIMP.tif')
ismip_1km = ismip_1km.rio.write_crs('EPSG:3413')

xsize = water.shape[2]
ysize = water.shape[1]

#%%

# Tile the grid to make resampling possible
for i in range(0, xsize, tile_size_x):
    for j in range(0, ysize, tile_size_y):
        if os.path.exists(str(out_path2) + str(output_filename) + str(i) + "_" + str(j) + ".tif"):
            pass
        else:
            com_string = "gdal_translate -of GTIFF -srcwin " + str(i)+ ", " + str(j) + ", " + str(tile_size_x) + ", " + str(tile_size_y) + " " + str(in_path) + str(crevasse_filename) + " " + str(out_path2) + str(output_filename) + str(i) + "_" + str(j) + ".tif"
            os.system(com_string)

# Tile the grid to make resampling possible
for i in range(0, xsize, tile_size_x):
    for j in range(0, ysize, tile_size_y):
        if os.path.exists(str(out_path1) + str(output_filename) + str(i) + "_" + str(j) + ".tif"):
            pass
        else:
            com_string = "gdal_translate -of GTIFF -srcwin " + str(i)+ ", " + str(j) + ", " + str(tile_size_x) + ", " + str(tile_size_y) + " " + str(in_path) + str(water_filename) + " " + str(out_path1) + str(output_filename) + str(i) + "_" + str(j) + ".tif"
            os.system(com_string)

#%%

# Define tiles
tiles1 = sorted(glob.glob(out_path1 + '*.tif'))
tiles2 = sorted(glob.glob(out_path2 + '*.tif'))

# Downscale
scale_factor = 10 / resolution

# Reproject
for t in range(len(tiles1)):
    
    # Get path and filename seperately 
    infilepath, infilename = os.path.split(tiles1[t]) 
    # Get file name without extension            
    infileshortname, extension = os.path.splitext(infilename)
    
    if os.path.exists(path + 'reproject/' + infileshortname + '.tif'):
        pass
    else:
    
        # Read
        w = rio.open_rasterio(tiles1[t], masked=True)
        c = rio.open_rasterio(tiles2[t], masked=True)
        
        if w.isnull().all().values:
            
            # Define new height and width
            new_width = int(w.rio.width * scale_factor)
            new_height = int(w.rio.height * scale_factor)
            
            # Downsample
            downsample = w.rio.reproject(
                w.rio.crs,
                shape=(new_height, new_width),
                resampling=Resampling.average,
            )
            
            # Save
            downsample.astype('float').rio.to_raster(path + 'reproject/' + infileshortname + '.tif', tiled=True)
    
        else:
            
            # Subtract
            data = xr.where(np.isnan(c), w.values, w.values - c.values)
            
            # Retain the x and y dimensions and coordinates
            x_coords = w["x"]
            y_coords = w["y"]
    
            ds = xr.DataArray(data=data, dims=["band", "y", "x"], 
                                  coords={"x": x_coords, "y": y_coords})
            
            # Assign CRS to the new DataArray
            ds = ds.rio.write_crs(w.rio.crs)
            
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
            downsample.astype('float').rio.to_raster(path + 'reproject/' + infileshortname + '.tif', tiled=True)


#%%

# Mosaic tiles

"""
Generate a text file containing files:
    
    ls -1 *.tif > /Users/jr555/Documents/tiff_list.txt

Then:
    
    gdal_merge.py -o /Users/jr555/Documents/surface_water_mask_2018_1km.tif --optfile /Users/jr555/Documents/tiff_list.txt

"""



#%%

# Resample to match extent of the ISMIP 1 km grid
ds = rio.open_rasterio('/Users/jr555/Documents/surface_water_mask_' + year + '_1km.tif', masked=True)
ds = ds.rio.write_crs('EPSG:3413')

#%%

ds_match = ds.rio.reproject_match(ismip_1km)
#ds_match = ds_match.where(ds_match != 1.7976931348623157e+308)
ds_match.rio.to_raster(path + 'surface_water_mask_' + year + '_1km.tif', tiled=True)
ds_match.to_netcdf(path + 'surface_water_mask_' + year + '_1km.nc')


#%%











