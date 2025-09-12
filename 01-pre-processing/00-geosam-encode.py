#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Encode images using Geo-SAM

"""

# Import packages
from geosam import ImageEncoder
import glob
import pandas as pd

# Choose glacier
glacier = 'isunguata'
#glacier = 'russell'

# Define paths
path = '/Users/jr555/Documents/research/hydrology/drone/' + glacier + '/'
checkpoint_path = '/Users/jr555/Documents/research/hydrology/sam_vit_l_0b3195.pth'

# Define savepath
feature_dir = path + 'GeoSAM'

# Define ortho files
ortho_files = sorted(glob.glob(path + 'Orthomosaics/' + '*.tif'))

# Define list of orthos with water in them
water = pd.read_csv(path + 'water-orthos.csv')

# Init ImageEncoder
img_encoder = ImageEncoder(checkpoint_path)

#%%

for i in range(len(water)):
    
    process_file = []
    for file in ortho_files:
        if water['file'].iloc[i] in file:
            process_file.append(file)
            
    # Encode image
    img_encoder.encode_image(process_file[0], feature_dir)



