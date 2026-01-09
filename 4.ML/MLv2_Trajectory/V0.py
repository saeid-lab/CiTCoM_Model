# Importing a few libraries:
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import re
from pathlib import Path
import mdtraj as md
import nglview as nv
import argparse
from collections import defaultdict
pd.set_option('display.max_columns', None)
pd.set_option('display.max_row', None)

# Converting trajectory to pdb
# gmx trjconv -s ../md-2.tpr -f ../md-1EHZ_whole_fit_Rep123.xtc -n ../index.ndx -o pdbs.pdb -skip 1000

# analyzing the trajectory in pdb format using dssr:
# x3dna-dssr --input=pdbs.pdb --output=pdbs.json --json --prefix=traj --md

json_url = './pdbs.json'
with open(json_url) as f:
    dssr_json = json.load(f)
    
    
    
models = dssr_json['models']

nts_list = []
for nt in dssr_json['models'][0]['parameters']['nts']:
    index = nt['index']
    name = nt['nt_id'].split(':')[1]
    nts_list.append([index, name])

nt_ids = [nt for _, nt in nts_list]
nt_index = {nt: i for i, nt in enumerate(nt_ids)}

n_nt = len(nts_list)
n_frames = len(models)


pair_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
for f, m in enumerate(models):
    for pair in m['parameters'].get('pairs', []):
        i = nt_index[pair['nt1'].split(':')[1]]
        j = nt_index[pair['nt2'].split(':')[1]]
        pair_matrix[i, f] = 1
        pair_matrix[j, f] = 1

pair_df = pd.DataFrame(pair_matrix, index=nt_ids)


wc_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
non_wc_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
for f, m in enumerate(models):
    for pair in m['parameters'].get('pairs', []):
        i = nt_index[pair['nt1'].split(':')[1]]
        j = nt_index[pair['nt2'].split(':')[1]]
        
        # Check if WC pair (cWW classification in LW or DSSR field)
        lw_class = pair.get('LW', '')
        dssr_class = pair.get('DSSR', '')
        
        if 'cWW' in lw_class or 'cW-W' in dssr_class:
            wc_matrix[i, f] = 1
            wc_matrix[j, f] = 1
        else:
            non_wc_matrix[i, f] = 1
            non_wc_matrix[j, f] = 1

wc_df = pd.DataFrame(wc_matrix, index=nt_ids)
non_wc_df = pd.DataFrame(non_wc_matrix, index=nt_ids)


multi_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
for f, m in enumerate(models):
    for mp in m['parameters'].get('multiplets', []):
        nts = mp['nts_long'].split(',')
        for nt in nts:
            nt_name = nt.split(':')[1]
            multi_matrix[nt_index[nt_name], f] = 1

multi_df = pd.DataFrame(multi_matrix, index=nt_ids)



helices_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
for f, m in enumerate(models):
    for mp in m['parameters'].get('helices', []):
        index = mp['pairs']
        for nt in index:
            i = nt_index[nt['nt1'].split(':')[1]]
            j = nt_index[nt['nt2'].split(':')[1]]
            helices_matrix[i, f] = 1
            helices_matrix[j, f] = 1

helices_df = pd.DataFrame(helices_matrix, index=nt_ids)



stems_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
for f, m in enumerate(models):
    for mp in m['parameters'].get('stems', []):
        index = mp['pairs']
        for nt in index:
            i = nt_index[nt['nt1'].split(':')[1]]
            j = nt_index[nt['nt2'].split(':')[1]]
            stems_matrix[i, f] = 1
            stems_matrix[j, f] = 1

stems_df = pd.DataFrame(stems_matrix, index=nt_ids)


hairpin_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
for f, m in enumerate(models):
    for hairpin in m['parameters'].get('hairpins', []):
        nts_in_loop = hairpin.get('nts_long', '').split(',')
        for nt_id in nts_in_loop:
            nt_name = nt_id.split(':')[1] if ':' in nt_id else nt_id.strip()
            if nt_name in nt_index:
                i = nt_index[nt_name]
                hairpin_matrix[i, f] = 1

hairpin_df = pd.DataFrame(hairpin_matrix, index=nt_ids)


internal_loop_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
for f, m in enumerate(models):
    # Process internal loops
    for iloop in m['parameters'].get('iloops', []):
        nts_in_loop = iloop.get('nts_long', '').split(',')
        for nt_id in nts_in_loop:
            nt_name = nt_id.split(':')[1] if ':' in nt_id else nt_id.strip()
            if nt_name in nt_index:
                i = nt_index[nt_name]
                internal_loop_matrix[i, f] = 1
    
    # Process bulges (a special type of internal loop)
    for bulge in m['parameters'].get('bulges', []):
        nts_in_loop = bulge.get('nts_long', '').split(',')
        for nt_id in nts_in_loop:
            nt_name = nt_id.split(':')[1] if ':' in nt_id else nt_id.strip()
            if nt_name in nt_index:
                i = nt_index[nt_name]
                internal_loop_matrix[i, f] = 1

internal_loop_df = pd.DataFrame(internal_loop_matrix, index=nt_ids)



junction_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
for f, m in enumerate(models):
    for junction in m['parameters'].get('junctions', []):
        nts_in_loop = junction.get('nts_long', '').split(',')
        for nt_id in nts_in_loop:
            nt_name = nt_id.split(':')[1] if ':' in nt_id else nt_id.strip()
            if nt_name in nt_index:
                i = nt_index[nt_name]
                junction_matrix[i, f] = 1

junction_df = pd.DataFrame(junction_matrix, index=nt_ids)



fraction_paired = pair_df.mean(axis=1)
fraction_wc_paired = wc_df.mean(axis=1)
fraction_non_wc_paired = non_wc_df.mean(axis=1)
fraction_unpaired = 1 - fraction_paired
fraction_in_multiplet = multi_df.mean(axis=1)
fraction_in_helix = helices_df.mean(axis=1)
fraction_in_stem = stems_df.mean(axis=1)
fraction_in_hairpin = hairpin_df.mean(axis=1)
fraction_in_internal_loop = internal_loop_df.mean(axis=1)
fraction_in_junction = junction_df.mean(axis=1)


features_df = pd.DataFrame({
    'nucleotide': nt_ids,
    'fraction_paired': fraction_paired,
    'fraction_wc_paired': fraction_wc_paired,
    'fraction_non_wc_paired': fraction_non_wc_paired,
    'fraction_unpaired': fraction_unpaired,
    'fraction_in_helix': fraction_in_helix,
    'fraction_in_stem': fraction_in_stem,
    'fraction_in_hairpin': fraction_in_hairpin,
    'fraction_in_internal_loop': fraction_in_internal_loop,
    'fraction_in_junction': fraction_in_junction,
})


# Define categories and parameters
MINIMAL_CATEGORIES = {
    'puckering': ["C3'-endo", "C2'-endo", "C4'-exo", "C1'-exo", "O4'-endo", 
                  "C4'-endo", "C1'-endo", "C2'-exo", "C3'-exo", "O4'-exo", "other"],
    'sugar_class': ["~C3'-endo", "~C2'-endo", "other"],
    'form': ['A', 'B', 'other'],
    'glyco_bond': ['anti', 'syn', 'other'],
    'bb_type': ['BI', 'BII', 'other']
}

TORSION_PARAMS = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta',
                  'epsilon_zeta', 'chi', 'eta', 'theta', 'phase_angle', 'amplitude']

# Initialize storage
conformational_features = {field: {cat: np.zeros((n_nt, n_frames), dtype=np.uint8) 
                                    for cat in categories} 
                           for field, categories in MINIMAL_CATEGORIES.items()}
torsion_matrices = {param: np.full((n_nt, n_frames), np.nan) for param in TORSION_PARAMS}

# Single loop to extract both conformational states and torsion angles
for f, m in enumerate(models):
    for nt in m['parameters'].get('nts', []):
        nt_name = nt['nt_id'].split(':')[1]
        if nt_name not in nt_index:
            continue
        
        i = nt_index[nt_name]
        
        # Extract conformational categories
        for field, categories in MINIMAL_CATEGORIES.items():
            value = nt.get(field, '')
            category = value if value in categories else 'other'
            conformational_features[field][category][i, f] = 1
        
        # Extract torsion angles
        for param in TORSION_PARAMS:
            value = nt.get(param)
            if value is not None and value != '---':
                try:
                    torsion_matrices[param][i, f] = float(value)
                except (ValueError, TypeError):
                    pass

# Build conformational DataFrame
conformational_probs = {}
for field, matrices in conformational_features.items():
    for category, matrix in matrices.items():
        col_name = f"{field}_{category}".replace("'", "").replace("-", "_").replace("~", "")
        conformational_probs[col_name] = matrix.mean(axis=1)

conformational_df = pd.DataFrame(conformational_probs, index=nt_ids)

# Build torsion statistics DataFrame
torsion_stats = {}
for param, matrix in torsion_matrices.items():
    torsion_stats[f'{param}_mean'] = np.nanmean(matrix, axis=1)
    torsion_stats[f'{param}_std'] = np.nanstd(matrix, axis=1)
    torsion_stats[f'{param}_min'] = np.nanmin(matrix, axis=1)
    torsion_stats[f'{param}_max'] = np.nanmax(matrix, axis=1)

torsion_df = pd.DataFrame(torsion_stats, index=nt_ids)

# Combine all features
conformational_df.insert(0, 'nucleotide', nt_ids)
torsion_df.insert(0, 'nucleotide', nt_ids)

ml_features_complete = features_df.merge(conformational_df, on='nucleotide').merge(torsion_df, on='nucleotide')
ml_features_complete.to_csv('dssr_ml_features_complete.csv', index=False)

print(f"Complete ML dataset: {ml_features_complete.shape[0]} nucleotides × {ml_features_complete.shape[1]-1} features")
print(f"  - Structural: 11")
print(f"  - Conformational: {len(conformational_df.columns)-1}")
print(f"  - Torsion statistics: {len(torsion_df.columns)-1}")


