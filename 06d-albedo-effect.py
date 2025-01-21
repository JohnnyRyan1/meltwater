#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate albedo lowering effect of water.

"""

# Import packages
import pandas as pd

# Define path
path1 = '/Users/jr555/Library/CloudStorage/OneDrive-DukeUniversity/research/hydrology/data/'

# Import data
water_albedo_2018 = pd.read_csv(path1 + 'albedo-effects-2018.csv')
non_water_albedo_2018 = pd.read_csv(path1 + 'non-water-albedo-effect-2018.csv')
water_albedo_2019 = pd.read_csv(path1 + 'albedo-effects-2019.csv')
non_water_albedo_2019 = pd.read_csv(path1 + 'non-water-albedo-effect-2019.csv')

# Combine datasets
water_albedo = pd.concat([water_albedo_2018, water_albedo_2019])
non_water_albedo = pd.concat([non_water_albedo_2018, non_water_albedo_2019])

#%%

print("Albedo variability in 5x5 windows with water is +/- %.4f" % (water_albedo['std'].mean()))

print("Albedo variability in 5x5 windows without water is +/- %.4f" % (non_water_albedo['std'].mean()))

print("Local albedo variability is %.1f higher in 5x5 windows with water" % (water_albedo['std'].mean() / non_water_albedo['std'].mean()))

print('Meltwater ponding is highly correlated with albedo with a mean correlation coefficient (r) of %.2f' %(water_albedo['r'].mean()))

print('Mean albedo effect of meltwater ponding is  %.4f +/- %.4f' % (water_albedo['slope'].mean(), water_albedo['slope'].std()))

print('About %.1f%% of local albedo variability is likely not caused by meltwater ponding' % ((non_water_albedo['std'].mean() / water_albedo['std'].mean()) * 100))

print('It is not clear whether there this source of albedo variability increases or \
decreases our albedo effects but likely captured by the std. dev.')

#%%

# Mean albedo effect by elevation

















