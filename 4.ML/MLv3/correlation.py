#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np

md_feats   = pd.read_csv('./dssr_ml_features_extended_statistics7.csv')
shape_file = glob.glob('*.shape')[0]                # Gets first (only) .shape file
shape_data = pd.read_csv(shape_file, names=['1m7']) #shape_data = pd.read_csv('./1GID.shape', names=['1m7'])
print("Original shape data:")
print(shape_data.head())

# validation: numeric AND between 0-10
def is_valid_shape(value):
    try:
        num = float(value)
        return 0 <= num <= 10
    except:
        return False

valid_mask = shape_data['1m7'].apply(is_valid_shape)
print(f"Valid nucleotides: {valid_mask.sum()} out of {len(shape_data)}")

# Filter both datasets
md_feats_valid = md_feats[valid_mask]
shape_valid    = shape_data[valid_mask]['1m7']


numeric_feats = md_feats_valid.drop(['nucleotide'], axis=1).select_dtypes(include=[np.number])

pearson_corr = numeric_feats.corrwith(shape_valid, method='pearson')
spearman_corr = numeric_feats.corrwith(shape_valid, method='spearman')

corr_df = pd.DataFrame({
    'pearson_r': pearson_corr,
    'spearman_rho': spearman_corr,
    'n_nucleotides': len(shape_valid)
}).sort_values(by='pearson_r', key=lambda s: s.abs(), ascending=False)

print("\nTop correlations (clean data):")
print(corr_df.head(10))

corr_df.to_csv('correlations7_clean.csv')
print("\nSaved to 'correlations7_clean.csv'")


'''
numeric_feats = md_feats.drop(['nucleotide'], axis=1)
pearson_corr = numeric_feats.corrwith(shape_data['1m7'], method='pearson')
spearman_corr = numeric_feats.corrwith(shape_data['1m7'], method='spearman')

corr_df = pd.DataFrame({
    'pearson_r': pearson_corr,
    'spearman_rho': spearman_corr
}).sort_values(by='pearson_r', key=lambda s: s.abs(), ascending=False)

print(corr_df.head(10))

corr_df.to_csv('correlations7.csv')
'''

