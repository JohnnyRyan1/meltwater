#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Downscale MERRA-2 shortwave downward radiation.

"""

# Import modules
import xarray as xr
from scipy import stats
import numpy as np
from pyresample import kd_tree, geometry
import netCDF4

#%%

# Define path
path = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

# Import MERRA-2
merra = xr.open_dataset(path + 'merra/swd_2018.nc')

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Define mask
mask = ismip_1km['GIMP'].values

# Define target projection
target_def = geometry.SwathDefinition(lons=ismip_1km['lon'].values, lats=ismip_1km['lat'].values)

# Define MERRA-2 elevations
elev_file = xr.open_dataset(path + 'MERRA2_101.const_2d_asm_Nx.00000000.nc4.nc4')

# Convert geopotential height to elevation
elev = elev_file['PHIS'].values[0,:,:] / 9.81

# Extract lat/lon
merra_lon, merra_lat = np.meshgrid(merra['longitude'].values, merra['latitude'].values)

# Define source projection
source_def = geometry.SwathDefinition(lons=merra_lon, lats=merra_lat)

#%%

def produce_stats(climatology):
        
    p, a, b, r = np.zeros((55,105)), np.zeros((55,105)), np.zeros((55,105)), np.zeros((55,105))
    
    for i in range(b.shape[0]):
        for j in range(b.shape[1]):
            
            imin = i-2
            imax = i+3
            jmin = j-2
            jmax = j+3
            
            if imin < 0:
                imin=0
            if jmin < 0:
                jmin=0
                   
            # Elevations
            elevations = elev[imin:imax,jmin:jmax]
            
            if np.sum(elevations) == 0:
                p[i, j], a[i, j],b[i, j], r[i, j] = np.nan, np.nan, np.nan, np.nan
            else:
                           
                # Find 24 neighboring pixels - about ~100 km given that MERRA-2 grid cell size is ~50 x 50 km
                # (40075*np.cos(np.deg2rad(70)))/360
                neighbors = climatology[imin:imax,jmin:jmax]
                
                # Linear regression
                coeffs = stats.linregress(elevations.flatten(), neighbors.flatten())
                
                # Append to grid
                p[i, j], a[i, j],b[i, j], r[i, j] = coeffs[3], coeffs[1], coeffs[0], coeffs[2]
    
    return p, a, b, r

def resample_grids(p, a, b, r, orig):
    
    # Reproject
    p_resample = kd_tree.resample_nearest(source_def, p, target_def, radius_of_influence=50000)
    a_resample = kd_tree.resample_nearest(source_def, a, target_def, radius_of_influence=50000)
    b_resample = kd_tree.resample_nearest(source_def, b, target_def, radius_of_influence=50000)
    r_resample = kd_tree.resample_nearest(source_def, r, target_def, radius_of_influence=50000)
    orig_resample = kd_tree.resample_nearest(source_def, orig, target_def, radius_of_influence=50000)
    
    # Mask
    p_resample[mask == 0] = np.nan
    a_resample[mask == 0] = np.nan
    b_resample[mask == 0] = np.nan
    r_resample[mask == 0] = np.nan
    orig_resample[mask == 0] = np.nan
    
    return p_resample, a_resample, b_resample, r_resample, orig_resample

def downscale(climatology):
    
    # Downscale stats
    p, a, b, r = produce_stats(climatology)
    
    # Resample
    p_resample, a_resample, b_resample, r_resample, orig_resample = resample_grids(p, a, b, r, climatology)

    # Downscale
    downscale = a_resample + (ismip_1km['SRF'].values * b_resample)
    
    # Replace non-significant values to the original grid value
    downscale[p_resample > 0.1] = orig_resample[p_resample > 0.1]

    return downscale, p_resample, r_resample

#%%

# Produce climatology for allsky shortwave between Jul 25 and Jul 30
swd_allsky = np.nanmean(merra['swd_allsky'].values[54:59,:,:], axis=0)

# Downscale daily climatology for allsky shortwave
swd_allsky_inter, swd_allsky_inter_p, swd_allsky_inter_r = downscale(swd_allsky)

# Produce climatology for allsky shortwave between Jul 25 and Jul 30
swd_clrsky = np.nanmean(merra['swd_clrsky'].values[54:59, :, :], axis=0)

# Downscale daily climatology for allsky shortwave
swd_clrsky_inter, swd_clrsky_inter_p, swd_clrsky_inter_r = downscale(swd_clrsky)

#%%
###############################################################################
# Save 1 km dataset to NetCDF
###############################################################################
dataset = netCDF4.Dataset(path + 'swd-downscaled-2018.nc', 'w', format='NETCDF4_CLASSIC')
print('Creating %s' % path + 'swd-downscaled-2018.nc')
dataset.Title = "Downscaled allsky/clearsky downward shortwave radiation from MERRA-2"
import time
dataset.History = "Created " + time.ctime(time.time())
dataset.Projection = "WGS 84"
dataset.Reference = "Ryan, J. C. et al. (unpublished)"
dataset.Contact = "jonathan.ryan@duke.edu"
    
# Create new dimensions
lat_dim = dataset.createDimension('y', swd_clrsky_inter.shape[0])
lon_dim = dataset.createDimension('x', swd_clrsky_inter.shape[1])
    
# Define variable types
Y = dataset.createVariable('latitude', np.float32, ('y','x'))
X = dataset.createVariable('longitude', np.float32, ('y','x'))

y = dataset.createVariable('y', np.float32, ('y'))
x = dataset.createVariable('x', np.float32, ('x'))
    
# Define units
Y.units = "degrees"
X.units = "degrees"
   
# Create the actual 3D variable
allsky_swd_nc = dataset.createVariable('swd_allsky', np.float32, ('y','x'))
clrsky_swd_nc = dataset.createVariable('swd_clrsky', np.float32, ('y','x'))

# Write data to layers
Y[:] = ismip_1km['lat'].values
X[:] = ismip_1km['lon'].values
y[:] = ismip_1km['lat'].values[:,0]
x[:] = ismip_1km['lon'].values[0,:]
allsky_swd_nc[:] = swd_allsky_inter.astype(np.float32)
clrsky_swd_nc[:] = swd_clrsky_inter.astype(np.float32)

print('Writing data to %s' % path + 'swd-downscaled-2018.nc')
    
# Close dataset
dataset.close()

#%%