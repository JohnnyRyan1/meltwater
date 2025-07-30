#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 10:58:43 2025

@author: jr555
"""

import geopandas as gpd
from shapely.geometry import LineString
import pandas as pd

gdf = gpd.read_file('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/drone/wq7/wq7_merged.shp')


lineStringObj = LineString( [[a.x, a.y] for a in gdf.geometry.values] )

line_df = pd.DataFrame()
line_df['Attrib'] = [1,]
line_gdf = gpd.GeoDataFrame(line_df, geometry=[lineStringObj,])



line_gdf.to_file('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/drone/wq7/wq7_merged-line.shp')


#%%

gdf = gpd.read_file('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/drone/russell/merged-images.shp')


lineStringObj = LineString( [[a.x, a.y] for a in gdf.geometry.values] )

line_df = pd.DataFrame()
line_df['Attrib'] = [1,]
line_gdf = gpd.GeoDataFrame(line_df, geometry=[lineStringObj,])



line_gdf.to_file('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/drone/russell/merged-images-line.shp')


#%%