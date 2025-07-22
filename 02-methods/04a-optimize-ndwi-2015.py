#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Optimize NDWI threshold.

"""

# Import packages
import glob
import os
import numpy as np
import rioxarray as rio
import geopandas as gpd
import pandas as pd
from _snic.lib import SNIC_main
from PIL import Image
from timeit import default_timer as timer
from cffi import FFI
import xarray as xr
from tqdm import tqdm
from joblib import Parallel, delayed
from sklearn.metrics import precision_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
import matplotlib.pyplot as plt

#%%

# Define training shapefiles
water_files = sorted(glob.glob('/Users/jr555/Documents/research/hydrology/evaluation-water/2015/*.shp'))
ice_files = sorted(glob.glob('/Users/jr555/Documents/research/hydrology/evaluation-ice/2015/*.shp'))

# Define some paths
raster_path = '/Users/jr555/Documents/research/hydrology/drone/20150721/Orthomosaics/'

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

water_labels = []
water_ndwi = []
ice_labels = []
ice_ndwi = []

for i in range(len(water_files)):
    
    # Get the path and filename separately
    infilepath, infilename = os.path.split(water_files[i])
    
    # Get the short name (filename without extension)
    infileshortname, extension = os.path.splitext(infilename)
    
    # Import raster
    file = raster_path + infileshortname + '.tif'
    image = rio.open_rasterio(file)

    # Import Segment Anything water
    water = gpd.read_file(water_files[i])
    
    # Import manually digittized ice
    ice = gpd.read_file(ice_files[i])

    # Convert to three-band image
    image = image.isel(band=[0,1,2])
    
    # Define pixel size
    transform = image.rio.transform()
    pixel_width = transform[0]
    pixel_area = pixel_width * pixel_width
    
    # Define number of segments
    numsuperpixels = int(image[0,:,:].size * pixel_area / super_pixel_area)
    
    # Segment
    labels, numlabels = segment(file,numsuperpixels,compactness,doRGBtoLAB)
    
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
    #ndwi_da.rio.to_raster(path1 + 'ndwi_image.tif')
    
    w_labels = []
    w_ndwi = []
    # Make training data for water
    for i in range(water.shape[0]):
        ndwi_values = ndwi_da.rio.clip([water.geometry.iloc[i]], water.crs).values.flatten()
        ndwi_values = ndwi_values[~np.isnan(ndwi_values)]
        w_labels.extend([1] * ndwi_values.shape[0])
        w_ndwi.extend(list(ndwi_values))
        
    i_labels = []
    i_ndwi = []
    # Make training data for ice
    for i in range(ice.shape[0]):
        ndwi_values = ndwi_da.rio.clip([ice.geometry.iloc[i]], ice.crs).values.flatten()
        ndwi_values = ndwi_values[~np.isnan(ndwi_values)]
        i_labels.extend([0] * ndwi_values.shape[0])
        i_ndwi.extend(list(ndwi_values))
        
    water_labels.extend(w_labels)
    water_ndwi.extend(w_ndwi)
    ice_labels.extend(i_labels)
    ice_ndwi.extend(i_ndwi)
    
#%%
# Make DataFrame
labels = water_labels + ice_labels
ndwi = water_ndwi + ice_ndwi
df = pd.DataFrame(list(zip(labels, ndwi)), columns=['y_true', 'X'])

threshold = np.arange(0, 0.31, 0.005)

f1 = []
precision = []
accuracy = []
recall = []
for t in threshold:
    
    # Make predictions
    df['y_pred'] = (df['X'] > t).astype(int)
    
    # Calculate F1-score 
    f1.append(f1_score(df['y_true'], df['y_pred']))
    precision.append(precision_score(df['y_true'], df['y_pred']))
    accuracy.append(accuracy_score(df['y_true'], df['y_pred']))
    recall.append(recall_score(df['y_true'], df['y_pred']))
    
# New DataFrame
new_df = pd.DataFrame(list(zip(threshold, f1, precision, accuracy, recall)), 
                      columns=['threshold', 'f1', 'precision', 'accuracy',
                               'recall'])

# Save
new_df.to_csv('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/ndwi-thresholds-2015.csv')


#%%

# Import data
new_df = pd.read_csv('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/ndwi-thresholds-2015.csv')

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

# Plot
fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(7,4), 
                                    layout='constrained')

ax1.plot(new_df['threshold'], new_df['recall'], lw=2, zorder=0,
         label='Recall', color=c1)
ax1.plot(new_df['threshold'], new_df['accuracy'], lw=2, zorder=0, 
         label='Accuracy', color=c2)
ax1.plot(new_df['threshold'], new_df['precision'], lw=2, zorder=0,
         label='Precision', color=c3)

ax1.set_xlabel("NDWI threshold", fontsize=12)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5, zorder=0)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.set_ylabel("Score", fontsize=12)
ax1.set_xlim(0, 0.3)
ax1.set_ylim(0.4, 1.03)
ax1.legend(loc=3, fontsize=12)


#%%










