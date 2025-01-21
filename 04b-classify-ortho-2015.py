#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Classify WQ7 orthomosaic

"""

# Import packages
import rioxarray as rio
import glob
import os
import numpy as np
import pandas as pd
from _snic.lib import SNIC_main
from PIL import Image
from timeit import default_timer as timer
from cffi import FFI
import xarray as xr
from geocube.vector import vectorize
from tqdm import tqdm
from joblib import Parallel, delayed

#%%

# Define path
path = '/Users/jr555/Documents/research/hydrology/'

# Define files
files = sorted(glob.glob(path + 'drone/20150721/Orthomosaics/*.tif'))

# Read
df = pd.read_csv(path + 'ndwi-thresholds-2015.csv')

# Threhold defined as the highest combination of F1-score and precision
ndwi_threshold = df['threshold'][np.argmax(df['accuracy'] + df['precision'])]

# Define desired superpixel area size
super_pixel_area = 10

# Set parameters and call the C function
compactness = 1
doRGBtoLAB = True # only works if it is a three channel image

#%%

# Define functions
def segment(imgname,numsuperpixels,compactness,doRGBtoLAB):
	#--------------------------------------------------------------
	# read image and change image shape from (h,w,c) to (c,h,w)
	#--------------------------------------------------------------
	img = Image.open(imgname)
	img = np.asarray(img)
    	
    	# Drop alpha channel
	img = img[:,:,0:3]
	print(img.shape)

	dims = img.shape
	h,w,c = dims[0],dims[1],1
	if len(dims) > 1:
		c = dims[2]
		img = img.transpose(2,0,1)
		print(c, "channels")
	
	#--------------------------------------------------------------
	# Reshape image to a single dimensional vector
	#--------------------------------------------------------------
	img = img.reshape(-1).astype(np.double)
	labels = np.zeros((h,w), dtype = np.int32)
	numlabels = np.zeros(1,dtype = np.int32)
	#--------------------------------------------------------------
	# Prepare the pointers to pass to the C function
	#--------------------------------------------------------------
	ffibuilder = FFI()
	pinp = ffibuilder.cast("double*", ffibuilder.from_buffer(img))
	plabels = ffibuilder.cast("int*", ffibuilder.from_buffer(labels.reshape(-1)))
	pnumlabels = ffibuilder.cast("int*", ffibuilder.from_buffer(numlabels))

	
	start = timer()
	SNIC_main(pinp,w,h,c,numsuperpixels,compactness,doRGBtoLAB,plabels,pnumlabels)
	end = timer()

	#--------------------------------------------------------------
	# Collect labels
	#--------------------------------------------------------------
	print("number of superpixels: ", numlabels[0])
	print("time taken in seconds: ", end-start)

	return labels.reshape(h,w),numlabels[0]

# Function to compute the median for a given label
def compute_median_for_label(label):
    # Get the pixel coordinates for the current region
    coords = np.where(labels == label)

    # Extract pixel values for each channel
    red_values = red_channel[coords]
    green_values = green_channel[coords]
    blue_values = blue_channel[coords]

    # Compute the median values for each channel
    median_red = np.median(red_values)
    median_green = np.median(green_values)
    median_blue = np.median(blue_values)

    return label, coords, median_red, median_green, median_blue


#%%
for file in files:
    
    # Get the path and filename separately
    infilepath, infilename = os.path.split(file)
    
    # Get the short name (filename without extension)
    infileshortname, extension = os.path.splitext(infilename)
    
    if os.path.exists(path + '/drone/20150721/Classified/' + infileshortname + '.tif'):
        pass
    else:
        
        # Import raster
        image = rio.open_rasterio(file)
    
        # Convert to three-band image
        image = image.isel(band=[0,1,2])
        
        if np.max(image[0,:,:]).values == 0:
            pass
        else:
        
            # Define pixel size
            transform = image.rio.transform()
            pixel_width = transform[0]
            pixel_area = pixel_width * pixel_width
            
            # Define number of segments
            numsuperpixels = int(image[0,:,:].size * pixel_area / super_pixel_area)
            
            # Segment
            labels, numlabels = segment(file, numsuperpixels, compactness, doRGBtoLAB)
            
            # Create a new image to hold the results
            median_image = np.zeros((3, image.values.shape[1], image.values.shape[2]))
            
            # Extract RGB channels beforehand to avoid redundant access
            red_channel, green_channel, blue_channel = image.values
            
            
            # Parallel processing for each label
            results = Parallel(n_jobs=-1)(delayed(compute_median_for_label)(label) 
                                          for label in tqdm(np.unique(labels), desc="Processing"))
        
            # Assign the computed median values to the new image
            for label, coords, median_red, median_green, median_blue in results:
                median_image[0, coords[0], coords[1]] = median_red
                median_image[1, coords[0], coords[1]] = median_green
                median_image[2, coords[0], coords[1]] = median_blue
                 
            # Make an segmented image
            median_image = median_image.astype(float)
            ndwi = (median_image[2, :, :] - median_image[0,:, :]) / \
                    (median_image[2, :, :] + median_image[0,:, :])
            ndwi_da = xr.DataArray(np.repeat(ndwi[np.newaxis, :, :], 3, axis=0), 
                                   dims=image.dims, coords=image.coords)
            ndwi_da = ndwi_da[2:,:,:]
        
            # Classify image
            classified = xr.where(ndwi_da > ndwi_threshold, 1, xr.where(ndwi_da <= ndwi_threshold, 0, 2))
            classified = classified.astype('uint8')
            classified = classified.rio.write_crs("EPSG:3413")
            classified.rio.to_raster(path + '/drone/20150721/Classified/' + infileshortname + '.tif')
            
            # Vectorize
            classified = xr.where(classified == 2, 0, classified)
            gdf = vectorize(classified)
            
            # Filter 
            gdf = gdf[gdf['_data'] == 1]
            gdf = gdf.set_crs("EPSG:3413")
            
            # Save to file
            gdf.to_file(path + 'drone/20150721/Polys/' + infileshortname + '.shp')
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        