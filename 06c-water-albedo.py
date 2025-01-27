#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Local albedo variability explained by surface water 

Note:
    This is a random sample 2000 pixels in each elevation bin.

"""

# Import packages
import rioxarray as rio
import xarray as xr
import numpy as np
import pandas as pd

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

# Define year
year = str(2019)

# Define sample size
sample_size = 2000

# Read raster data
sw = rio.open_rasterio(path1 + 'zhang/surface_water_mask_' + year +'_500m.tif')

# Replace NaN values with zeros
sw = sw.fillna(0)

mcd = rio.open_rasterio(path1 + 'mcd-' + year + '.tif')
mod = rio.open_rasterio(path1 + 'mod-' + year + '.tif')
x_grid, y_grid = np.meshgrid(mcd['x'], mcd['y'])

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path1 + '1km-ISMIP6-GIMP.nc')
elev = ismip_1km['SRF'].rio.write_crs('EPSG:3413')
elev_match = elev.rio.reproject_match(mcd)
mask = ismip_1km['GIMP'].rio.write_crs('EPSG:3413').astype(int)
mask_match = mask.rio.reproject_match(mcd)

# Define which MODIS data to use
modis = mcd

# Mask
modis = xr.where(mask_match, modis, np.nan)
sw = xr.where(mask_match, sw, np.nan)

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

# Compute water area by elevation
water_area, ice_area = [], []
for e in range(len(elevations) - 1):
    elevation_mask = (elev_match > elevations[e]) &\
        (elev_match < elevations[e+1])
    water_area.append(np.nansum(sw.values[elevation_mask.values])*0.25)
    ice_area.append(elevation_mask.values.sum()*0.25)

# Make a DataFrame
df = pd.DataFrame(list(zip(elevations, water_area, ice_area)), 
                  columns=['elevation', 'water_area', 'ice_area'])

# Export
df.to_csv(path1 + 'scales/water-area.csv', index=False)

#%%

# Define window size from 1 to 50 km
windows = np.array((1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90,
                    100))
kernel_size = (windows[:-1] + windows[1:])**2

# Define lists
final_data_list = []

for w in range(len(windows[:-1])):
    
    print('Window size is %.0f pixels' % ((windows[w] +  windows[w+1])**2))

    data_list = []
    
    for e in range(len(elevations) - 1):
        #print('Processing %.2f to %.2f m' %(bin_edges[i], bin_edges[i+1]))
        
        # Get first elevation band
        band = (elev_match > elevations[e]) & (elev_match < elevations[e+1]) &\
            (mask_match.values == 1)
        
        # Find cells
        x, y = np.where(band == 1)
        
        # Select number of samples
        np.random.seed(42)
        samples = np.random.choice(x.shape[0], size=sample_size, replace=False)
        #print('There are %.0f samples' % binned_counts[i])
        
        data_list1 = []

        for n in samples:
            
            # Get row/col
            coords = x[n], y[n]
            
            # Define window
            x0 = np.maximum(coords[0]-windows[w], 0)
            x1 = coords[0]+windows[w+1]
            y0 = np.maximum(coords[1]-windows[w], 0)
            y1 = coords[1]+windows[w+1]
                    
            # Identify nearest neighbors
            neighbors = np.array(sw[:,:,0][x0:x1,y0:y1]).flatten()
            albedo = np.array(modis[:,:,0][x0:x1,y0:y1]).flatten()
            
            # Only compute linear regression if number of albedo values is greater than 9
            if np.isfinite(albedo).sum() < 9:
                pass
            else:
                # Find grid cells where surface water is NaN
                indices = np.where(np.isnan(neighbors))
        
                # Remove elements with NaN surface water
                albedo_filtered = np.where(np.isnan(neighbors), np.nan, albedo)
                
                # Remove elements with NaN surface water AND albedo
                albedo_non_nan = albedo_filtered[np.isfinite(albedo_filtered)]
                water_non_nan = neighbors[np.isfinite(albedo_filtered)]
                
                # Append to list
                data_list.append((y_grid[coords],x_grid[coords], np.nanstd(albedo_non_nan),
                                  np.nanmean(albedo_non_nan), elev_match[coords].values, 
                                  albedo_non_nan.shape[0]))
            
    final_data_list.append(data_list)
    
stats_list = []

for i in range(len(final_data_list)):
    # Make a DataFrame
    df = pd.DataFrame(final_data_list[i])
    df.columns = ['y', 'x', 'std', 'mean_albedo', 'elevation', 'count']
    df['bin'] = pd.cut(df['elevation'], bins=elevations, right=False)
    group = df.groupby('bin', observed=True)['std'].mean()
    stats_list.append(group.values)

# Make a DataFrame
df = pd.DataFrame(stats_list)

# Export
df.to_csv(path1 + 'scales/water-2000.csv', index=False)

    
#%%

# Make a DataFrame
df = pd.DataFrame(data_list)
df.columns = ['y', 'x', 'std', 'mean_albedo', 'intercept', 'slope', 'r', 'p', 'elevation']
df['elevation'] = df['elevation'].astype(float)

# Filter
df = df[df['slope'] < 0]
df = df[df['p'] < 0.001]

df.to_csv(path1 + 'albedo-effects-' + year + '.csv')


#%%






