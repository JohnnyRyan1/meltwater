#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Intersect albedo effects with suprglacial lake shapefile to compute depths.

"""

import geopandas as gpd
from scipy import stats

points = gpd.read_file('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/albedo-effects-2019.shp')
lakes = gpd.read_file('/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/zhang/supraglacial_lakes_2019.shp')

# Perform the spatial join
combined = gpd.sjoin(points, lakes, how="inner", predicate="intersects")

combined['area'] = (combined['Area'] / 1e+6)


slope, intercept, r, p, se = stats.linregress(combined['area'], combined['slope'])


c = combined[['volume', 'area', 'slope']].groupby('volume').mean()

slope, intercept, r, p, se = stats.linregress(c['area'], c['slope'])


slope, intercept, r, p, se = stats.linregress(points['elevation'], points['slope'])
slope, intercept, r, p, se = stats.linregress(points['y'], points['slope'])








