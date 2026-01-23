#!/usr/bin/env python3
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

# gmx trjconv -s ../md-2.tpr -f ../md-1GID_whole_fit_Rep123.xtc -n ../index.ndx -o pdbs.pdb -skip 1000 (every 2ns)
# x3dna-dssr --input=pdbs.pdb --output=pdbs.json --json --prefix=traj --md

# === Load DSSR JSON ===
json_url = './pdbs.json'
with open(json_url) as f:
    dssr_json = json.load(f)

models = dssr_json['models']

# === Extract nucleotide list ===
nts_list = []
for nt in dssr_json['models'][0]['parameters']['nts']:
    index = nt['index']
    name = nt['nt_id'].split(':')[1]
    nts_list.append([index, name])

nt_ids = [nt for _, nt in nts_list]
nt_index = {nt: i for i, nt in enumerate(nt_ids)}
n_nt = len(nts_list)
n_frames = len(models)

print(f"DSSR data: {n_nt} nucleotides, {n_frames} frames from JSON")

# ============================================================
# ORIGINAL DSSR-BASED FEATURES (unchanged)
# ============================================================

# === Base pairing features ===
pair_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
for f, m in enumerate(models):
    for pair in m['parameters'].get('pairs', []):
        i = nt_index[pair['nt1'].split(':')[1]]
        j = nt_index[pair['nt2'].split(':')[1]]
        pair_matrix[i, f] = 1
        pair_matrix[j, f] = 1

pair_df = pd.DataFrame(pair_matrix, index=nt_ids)

# === Watson-Crick vs non-WC pairing ===
wc_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
non_wc_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)

for f, m in enumerate(models):
    for pair in m['parameters'].get('pairs', []):
        i = nt_index[pair['nt1'].split(':')[1]]
        j = nt_index[pair['nt2'].split(':')[1]]

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

# === Multiplets ===
multi_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
for f, m in enumerate(models):
    for mp in m['parameters'].get('multiplets', []):
        nts = mp['nts_long'].split(',')
        for nt in nts:
            nt_name = nt.split(':')[1]
            multi_matrix[nt_index[nt_name], f] = 1

multi_df = pd.DataFrame(multi_matrix, index=nt_ids)

# === Helices ===
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

# === Stems ===
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

# === Hairpins ===
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

# === Internal loops and bulges ===
internal_loop_matrix = np.zeros((n_nt, n_frames), dtype=np.uint8)
for f, m in enumerate(models):
    for iloop in m['parameters'].get('iloops', []):
        nts_in_loop = iloop.get('nts_long', '').split(',')
        for nt_id in nts_in_loop:
            nt_name = nt_id.split(':')[1] if ':' in nt_id else nt_id.strip()
            if nt_name in nt_index:
                i = nt_index[nt_name]
                internal_loop_matrix[i, f] = 1

    for bulge in m['parameters'].get('bulges', []):
        nts_in_loop = bulge.get('nts_long', '').split(',')
        for nt_id in nts_in_loop:
            nt_name = nt_id.split(':')[1] if ':' in nt_id else nt_id.strip()
            if nt_name in nt_index:
                i = nt_index[nt_name]
                internal_loop_matrix[i, f] = 1

internal_loop_df = pd.DataFrame(internal_loop_matrix, index=nt_ids)

# === Junctions ===
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

# === Compute fractions ===
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

# === Conformational and torsion angle features ===
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

conformational_features = {field: {cat: np.zeros((n_nt, n_frames), dtype=np.uint8)
                                   for cat in categories}
                           for field, categories in MINIMAL_CATEGORIES.items()}

torsion_matrices = {param: np.full((n_nt, n_frames), np.nan) for param in TORSION_PARAMS}

for f, m in enumerate(models):
    for nt in m['parameters'].get('nts', []):
        nt_name = nt['nt_id'].split(':')[1]
        if nt_name not in nt_index:
            continue
        i = nt_index[nt_name]

        for field, categories in MINIMAL_CATEGORIES.items():
            value = nt.get(field, '')
            category = value if value in categories else 'other'
            conformational_features[field][category][i, f] = 1

        for param in TORSION_PARAMS:
            value = nt.get(param)
            if value is not None and value != '---':
                try:
                    torsion_matrices[param][i, f] = float(value)
                except (ValueError, TypeError):
                    pass

conformational_probs = {}
for field, matrices in conformational_features.items():
    for category, matrix in matrices.items():
        col_name = f"{field}_{category}".replace("'", "").replace("-", "_").replace("~", "")
        conformational_probs[col_name] = matrix.mean(axis=1)

conformational_df = pd.DataFrame(conformational_probs, index=nt_ids)

torsion_stats = {}
for param, matrix in torsion_matrices.items():
    torsion_stats[f'{param}_mean'] = np.nanmean(matrix, axis=1)
    torsion_stats[f'{param}_std'] = np.nanstd(matrix, axis=1)
    torsion_stats[f'{param}_min'] = np.nanmin(matrix, axis=1)
    torsion_stats[f'{param}_max'] = np.nanmax(matrix, axis=1)

torsion_df = pd.DataFrame(torsion_stats, index=nt_ids)

# ============================================================
# GEOMETRICAL FEATURES WITH EXTENDED FLUCTUATION STATISTICS
# ============================================================

print("\n=== Loading MD trajectory for geometrical features ===")
traj_url = './pdbs.pdb'  # ← USER: Update this path
traj = md.load(filename_or_filenames=traj_url)

print(f'Trajectory info:')
print(f'  - Chains: {traj.n_chains}')
print(f'  - Residues: {traj.n_residues}')
print(f'  - Frames: {traj.n_frames}')

# Define atom groups for centroids
nucleobase_atoms = {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"}
sugar_atoms = {"C1'", "C2'", "C3'", "C4'", "O4'"}
phosphate_atoms = {"P"}
O2p_atom = {"O2'"}

residues = list(traj.top.residues)
num_nucleotides = len(residues)

def compute_centroid(residue, traj, atom_names):
    """Compute geometric center for a given set of atoms in a residue.
    Returns: (n_frames, 3) - one centroid per frame"""
    atom_indices = [atom.index for atom in residue.atoms if atom.name in atom_names]
    if len(atom_indices) == 0:
        return np.full((traj.n_frames, 3), np.nan)
    positions = traj.xyz[:, atom_indices, :] * 10  # nm → Å
    return np.mean(positions, axis=1)

# Compute centroids for each residue across ALL frames
print(f"\nComputing centroids for {num_nucleotides} residues across {traj.n_frames} frames...")
nucleobase_centroids = np.array([compute_centroid(res, traj, nucleobase_atoms) for res in residues])
sugar_centroids = np.array([compute_centroid(res, traj, sugar_atoms) for res in residues])
phosphate_centroids = np.array([compute_centroid(res, traj, phosphate_atoms) for res in residues])
O2p_coords = np.array([
    traj.xyz[:, [atom.index for atom in res.atoms if atom.name in O2p_atom], :] * 10
    if any(atom.name in O2p_atom for atom in res.atoms)
    else np.full((traj.n_frames, 1, 3), np.nan)
    for res in residues
]).squeeze(2)

print(f"Centroid array shapes:")
print(f"  nucleobase_centroids: {nucleobase_centroids.shape}")
print(f"  sugar_centroids: {sugar_centroids.shape}")
print(f"  phosphate_centroids: {phosphate_centroids.shape}")
print(f"  O2p_coords: {O2p_coords.shape}")

def compute_distances_per_frame(arr1, arr2, shift):
    """Compute distance between residue i and i+shift for each frame.
    Returns: (n_residues-shift, n_frames)"""
    return np.linalg.norm(arr1[:-shift] - arr2[shift:], axis=2)

def compute_all_distance_stats(distances_per_frame):
    """Compute comprehensive statistics for distance array.
    Args: distances_per_frame shape (n_residues, n_frames)
    Returns: dict with mean, std, ptp, p90, cv, var"""
    stats = {}
    stats['mean'] = np.nanmean(distances_per_frame, axis=1)
    stats['std'] = np.nanstd(distances_per_frame, axis=1)
    stats['var'] = np.nanvar(distances_per_frame, axis=1)
    stats['ptp'] = np.ptp(distances_per_frame, axis=1)  # peak-to-peak (range)
    stats['p90'] = np.nanpercentile(distances_per_frame, 90, axis=1)
    # Coefficient of variation: std/mean (normalized fluctuation)
    stats['cv'] = stats['std'] / (stats['mean'] + 1e-8)
    return stats

# === Distance features with 6 statistics each ===
print(f"\nComputing distance features with extended statistics...")
geo_features = {}

for k in [1, 2, 3]:
    # RiboseCentroid - RiboseCentroid
    dists = compute_distances_per_frame(sugar_centroids, sugar_centroids, k)
    stats = compute_all_distance_stats(dists)
    for stat_name, values in stats.items():
        geo_features[f"RiboseCentroid_RiboseCentroid(+{k})_{stat_name}"] = values

    # baseCentroid - baseCentroid
    dists = compute_distances_per_frame(nucleobase_centroids, nucleobase_centroids, k)
    stats = compute_all_distance_stats(dists)
    for stat_name, values in stats.items():
        geo_features[f"baseCentroid_baseCentroid(+{k})_{stat_name}"] = values

    # baseCentroid - P
    dists = compute_distances_per_frame(nucleobase_centroids, phosphate_centroids, k)
    stats = compute_all_distance_stats(dists)
    for stat_name, values in stats.items():
        geo_features[f"baseCentroid_P(+{k})_{stat_name}"] = values

    # O2' - P
    dists = compute_distances_per_frame(O2p_coords, phosphate_centroids, k)
    stats = compute_all_distance_stats(dists)
    for stat_name, values in stats.items():
        geo_features[f"O2'_P(+{k})_{stat_name}"] = values

# === Additional atomic distance features (C1'-C1'+1 and C2-C2+1) ===
print(f"Computing C1'-C1' and C2-C2 distances...")

# Extract C1' and C2 atoms for each residue
def get_atom_coords(residues, traj, atom_name):
    """Get coordinates for specific atom across all residues and frames.
    Returns: (n_residues, n_frames, 3) or NaN if atom not found"""
    coords = []
    for res in residues:
        atom_indices = [atom.index for atom in res.atoms if atom.name == atom_name]
        if len(atom_indices) == 1:
            positions = traj.xyz[:, atom_indices[0], :] * 10  # nm → Å
            coords.append(positions)
        else:
            coords.append(np.full((traj.n_frames, 3), np.nan))
    return np.array(coords)

# Get C1' coordinates (sugar atom)
C1p_coords = get_atom_coords(residues, traj, "C1'")

# Get C2 coordinates (base carbonyl for pyrimidines, or C2 for purines)
# Note: C2 exists in all nucleotides (purines and pyrimidines)
C2_coords = get_atom_coords(residues, traj, "C2")

print(f"  C1' coordinates shape: {C1p_coords.shape}")
print(f"  C2 coordinates shape: {C2_coords.shape}")

# Compute C1'-C1'+1 distances
def compute_sequential_atom_distance(coords, shift=1):
    """Compute distance between atom i and atom i+shift.
    Returns: (n_residues-shift, n_frames)"""
    if coords.shape[0] <= shift:
        return np.array([])
    distances = np.linalg.norm(coords[:-shift] - coords[shift:], axis=2)
    return distances

# C1'-C1'+1 distance with all statistics
C1p_C1p_dists = compute_sequential_atom_distance(C1p_coords, shift=1)
if len(C1p_C1p_dists) > 0:
    stats = compute_all_distance_stats(C1p_C1p_dists)
    for stat_name, values in stats.items():
        geo_features[f"C1'_C1'(+1)_{stat_name}"] = values

# C2-C2+1 distance with all statistics
C2_C2_dists = compute_sequential_atom_distance(C2_coords, shift=1)
if len(C2_C2_dists) > 0:
    stats = compute_all_distance_stats(C2_C2_dists)
    for stat_name, values in stats.items():
        geo_features[f"C2_C2(+1)_{stat_name}"] = values

print(f"  Added C1'-C1'(+1) and C2-C2(+1) features with 6 statistics each")

# === Local density counts ===
print(f"Computing local density features...")
heavy_atom_coords = traj.xyz[:, [a.index for a in traj.top.atoms if a.element.symbol != 'H'], :] * 10
heavy_atom_resid = np.array([a.residue.index for a in traj.top.atoms if a.element.symbol != 'H'])

def count_within_cutoff(res_idx, cutoff):
    """Count heavy atoms within cutoff from this residue's sugar centroid."""
    center = sugar_centroids[res_idx]
    dists = np.linalg.norm(heavy_atom_coords - center[:, None, :], axis=2)
    mask = heavy_atom_resid != res_idx
    counts = np.sum((dists < cutoff) & mask[None, :], axis=1)
    return np.mean(counts)

for cutoff in [4, 6, 8]:
    geo_features[f"count_within_{cutoff}A"] = [count_within_cutoff(i, cutoff) for i in range(num_nucleotides)]

# === Angle features with 6 statistics each ===
print(f"Computing angle features with extended statistics...")

def compute_angle_all_frames(a, b, c):
    """Compute angle ABC (vertex at B) for all frames.
    Args: a, b, c are (n_frames, 3)
    Returns: (n_frames,) array of angles in degrees"""
    ab = b - a
    cb = b - c

    ab_norm = np.linalg.norm(ab, axis=-1, keepdims=True)
    cb_norm = np.linalg.norm(cb, axis=-1, keepdims=True)

    ab = ab / (ab_norm + 1e-8)
    cb = cb / (cb_norm + 1e-8)

    cosang = np.clip(np.sum(ab * cb, axis=-1), -1.0, 1.0)
    angles = np.degrees(np.arccos(cosang))

    return angles

def compute_all_angle_stats(angles_per_frame):
    """Compute comprehensive statistics for angle array.
    Args: angles_per_frame shape (n_frames,)
    Returns: dict with mean, std, ptp, p90, cv, var"""
    stats = {}
    stats['mean'] = np.nanmean(angles_per_frame)
    stats['std'] = np.nanstd(angles_per_frame)
    stats['var'] = np.nanvar(angles_per_frame)
    stats['ptp'] = np.ptp(angles_per_frame)
    stats['p90'] = np.nanpercentile(angles_per_frame, 90)
    stats['cv'] = stats['std'] / (stats['mean'] + 1e-8)
    return stats

angle_features_dict = {stat: [] for stat in ['mean', 'std', 'var', 'ptp', 'p90', 'cv']}

for i in range(num_nucleotides):
    # S-B-P angle
    angles = compute_angle_all_frames(
        sugar_centroids[i],
        nucleobase_centroids[i],
        phosphate_centroids[i]
    )
    stats = compute_all_angle_stats(angles)
    for stat_name in angle_features_dict.keys():
        angle_features_dict[stat_name].append({f"S-B-P_{stat_name}": stats[stat_name]})

    # O2'-P-P(+1) angle
    if i < num_nucleotides - 1:
        angles = compute_angle_all_frames(
            O2p_coords[i],
            phosphate_centroids[i],
            phosphate_centroids[i + 1]
        )
        stats = compute_all_angle_stats(angles)
        for stat_name in angle_features_dict.keys():
            angle_features_dict[stat_name][i][f"O2'-P-P(+1)_{stat_name}"] = stats[stat_name]
    else:
        for stat_name in angle_features_dict.keys():
            angle_features_dict[stat_name][i][f"O2'-P-P(+1)_{stat_name}"] = np.nan

    # O2'-P-S angle
    angles = compute_angle_all_frames(
        O2p_coords[i],
        phosphate_centroids[i],
        sugar_centroids[i]
    )
    stats = compute_all_angle_stats(angles)
    for stat_name in angle_features_dict.keys():
        angle_features_dict[stat_name][i][f"O2'-P-S_{stat_name}"] = stats[stat_name]

    # O2'-S-B angle
    angles = compute_angle_all_frames(
        O2p_coords[i],
        sugar_centroids[i],
        nucleobase_centroids[i]
    )
    stats = compute_all_angle_stats(angles)
    for stat_name in angle_features_dict.keys():
        angle_features_dict[stat_name][i][f"O2'-S-B_{stat_name}"] = stats[stat_name]

    # B-P-B+1 angle
    if i < num_nucleotides - 1:
        angles = compute_angle_all_frames(
            nucleobase_centroids[i],
            phosphate_centroids[i],
            nucleobase_centroids[i+1]
        )
        stats = compute_all_angle_stats(angles)
        for stat_name in angle_features_dict.keys():
            angle_features_dict[stat_name][i][f"B-P-B+1_{stat_name}"] = stats[stat_name]
    else:
        for stat_name in angle_features_dict.keys():
            angle_features_dict[stat_name][i][f"B-P-B+1_{stat_name}"] = np.nan

    # B-P-B+2 angle
    if i < num_nucleotides - 2:
        angles = compute_angle_all_frames(
            nucleobase_centroids[i],
            phosphate_centroids[i],
            nucleobase_centroids[i+2]
        )
        stats = compute_all_angle_stats(angles)
        for stat_name in angle_features_dict.keys():
            angle_features_dict[stat_name][i][f"B-P-B+2_{stat_name}"] = stats[stat_name]
    else:
        for stat_name in angle_features_dict.keys():
            angle_features_dict[stat_name][i][f"B-P-B+2_{stat_name}"] = np.nan

    # S-P-S+1 angle
    if i < num_nucleotides - 1:
        angles = compute_angle_all_frames(
            sugar_centroids[i],
            phosphate_centroids[i],
            sugar_centroids[i+1]
        )
        stats = compute_all_angle_stats(angles)
        for stat_name in angle_features_dict.keys():
            angle_features_dict[stat_name][i][f"S-P-S+1_{stat_name}"] = stats[stat_name]
    else:
        for stat_name in angle_features_dict.keys():
            angle_features_dict[stat_name][i][f"S-P-S+1_{stat_name}"] = np.nan

    # S-P-S+2 angle
    if i < num_nucleotides - 2:
        angles = compute_angle_all_frames(
            sugar_centroids[i],
            phosphate_centroids[i],
            sugar_centroids[i+2]
        )
        stats = compute_all_angle_stats(angles)
        for stat_name in angle_features_dict.keys():
            angle_features_dict[stat_name][i][f"S-P-S+2_{stat_name}"] = stats[stat_name]
    else:
        for stat_name in angle_features_dict.keys():
            angle_features_dict[stat_name][i][f"S-P-S+2_{stat_name}"] = np.nan

# Combine all angle statistics into DataFrames
angles_combined_df = pd.DataFrame()
for stat_name in ['mean', 'std', 'var', 'ptp', 'p90', 'cv']:
    df_stat = pd.DataFrame(angle_features_dict[stat_name])
    angles_combined_df = pd.concat([angles_combined_df, df_stat], axis=1)

# Align geometrical features with nucleotide list
geo_features_aligned = {}
for key, values in geo_features.items():
    if len(values) < num_nucleotides:
        values = list(values) + [np.nan] * (num_nucleotides - len(values))
    geo_features_aligned[key] = values

geo_distances_df = pd.DataFrame(geo_features_aligned)


# === Additional atomic-level angles from paper ===
print(f"Computing 8 additional atomic angles from literature...")

# Get atomic coordinates for new angles
C1p_coords = get_atom_coords(residues, traj, "C1'")  # Already have this
C4p_coords = get_atom_coords(residues, traj, "C4'")
C2_coords = get_atom_coords(residues, traj, "C2")    # Already have this
O2p_coords = O2p_coords  # Already have this from earlier
O3p_coords = get_atom_coords(residues, traj, "O3'")
O5p_coords = get_atom_coords(residues, traj, "O5'")
P_coords = phosphate_centroids  # Already have this

print(f"  C4' coordinates shape: {C4p_coords.shape}")
print(f"  O3' coordinates shape: {O3p_coords.shape}")
print(f"  O5' coordinates shape: {O5p_coords.shape}")

# Helper function for angles across residues
def compute_angle_sequential(atom1_i, atom2_vertex, atom3_i_plus_k, shift=1):
    """Compute angle A-B-C where C is at position i+shift.
    Args:
        atom1_i: (n_res, n_frames, 3) - atom A at position i
        atom2_vertex: (n_res, n_frames, 3) - atom B (vertex)
        atom3_i_plus_k: (n_res, n_frames, 3) - atom C at position i+k
        shift: offset for atom3
    Returns: (n_res-shift, n_frames) array of angles
    """
    if atom1_i.shape[0] <= shift:
        return np.array([])

    # For angle at residue i, need: A(i), B(i or i+shift), C(i+shift)
    angles_all = []
    for i in range(atom1_i.shape[0] - shift):
        # Determine which residue the vertex is at
        if atom2_vertex.shape[0] > i + shift:
            # Vertex at i+shift (most common case)
            a = atom1_i[i]  # A at i
            b = atom2_vertex[i + shift]  # B at i+shift
            c = atom3_i_plus_k[i + shift]  # C at i+shift
        else:
            # Vertex at i
            a = atom1_i[i]
            b = atom2_vertex[i]
            c = atom3_i_plus_k[i + shift]

        angles = compute_angle_all_frames(a, b, c)
        angles_all.append(angles)

    return np.array(angles_all)

# Compute all 8 angles
atomic_angle_features = {}

# 1. C1' - C4' - P(+1)
for i in range(num_nucleotides - 1):
    angles = compute_angle_all_frames(
        C1p_coords[i],
        C4p_coords[i],
        P_coords[i + 1]
    )
    stats = compute_all_angle_stats(angles)
    if i == 0:
        for stat_name in stats.keys():
            atomic_angle_features[f"C1'-C4'-P(+1)_{stat_name}"] = []
    for stat_name, val in stats.items():
        atomic_angle_features[f"C1'-C4'-P(+1)_{stat_name}"].append(val)

# 2. C1' - P(+1) - C1'(+1)
for i in range(num_nucleotides - 1):
    angles = compute_angle_all_frames(
        C1p_coords[i],
        P_coords[i + 1],
        C1p_coords[i + 1]
    )
    stats = compute_all_angle_stats(angles)
    if i == 0:
        for stat_name in stats.keys():
            atomic_angle_features[f"C1'-P(+1)-C1'(+1)_{stat_name}"] = []
    for stat_name, val in stats.items():
        atomic_angle_features[f"C1'-P(+1)-C1'(+1)_{stat_name}"].append(val)

# 3. C2 - C1' - P(+1)  [CORRECTED: vertex at i, not i+1]
for i in range(num_nucleotides - 1):
    angles = compute_angle_all_frames(
        C2_coords[i],
        C1p_coords[i],
        P_coords[i + 1]
    )
    stats = compute_all_angle_stats(angles)
    if i == 0:
        for stat_name in stats.keys():
            atomic_angle_features[f"C2-C1'-P(+1)_{stat_name}"] = []
    for stat_name, val in stats.items():
        atomic_angle_features[f"C2-C1'-P(+1)_{stat_name}"].append(val)

# 4. C2 - C4' - P(+1)
for i in range(num_nucleotides - 1):
    angles = compute_angle_all_frames(
        C2_coords[i],
        C4p_coords[i],
        P_coords[i + 1]
    )
    stats = compute_all_angle_stats(angles)
    if i == 0:
        for stat_name in stats.keys():
            atomic_angle_features[f"C2-C4'-P(+1)_{stat_name}"] = []
    for stat_name, val in stats.items():
        atomic_angle_features[f"C2-C4'-P(+1)_{stat_name}"].append(val)

# 5. C2 - P(+1) - C2(+1)  [Complements C2_C2(+1) distance!]
for i in range(num_nucleotides - 1):
    angles = compute_angle_all_frames(
        C2_coords[i],
        P_coords[i + 1],
        C2_coords[i + 1]
    )
    stats = compute_all_angle_stats(angles)
    if i == 0:
        for stat_name in stats.keys():
            atomic_angle_features[f"C2-P(+1)-C2(+1)_{stat_name}"] = []
    for stat_name, val in stats.items():
        atomic_angle_features[f"C2-P(+1)-C2(+1)_{stat_name}"].append(val)

# 6. O2' - P(+1) - O5'(+1)  [Direct 2'-OH environment!]
for i in range(num_nucleotides - 1):
    angles = compute_angle_all_frames(
        O2p_coords[i],
        P_coords[i + 1],
        O5p_coords[i + 1]
    )
    stats = compute_all_angle_stats(angles)
    if i == 0:
        for stat_name in stats.keys():
            atomic_angle_features[f"O2'-P(+1)-O5'(+1)_{stat_name}"] = []
    for stat_name, val in stats.items():
        atomic_angle_features[f"O2'-P(+1)-O5'(+1)_{stat_name}"].append(val)

# 7. O3' - P(+1) - O2'(+1)  [Backbone-2'-OH angle]
for i in range(num_nucleotides - 1):
    angles = compute_angle_all_frames(
        O3p_coords[i],
        P_coords[i + 1],
        O2p_coords[i + 1]
    )
    stats = compute_all_angle_stats(angles)
    if i == 0:
        for stat_name in stats.keys():
            atomic_angle_features[f"O3'-P(+1)-O2'(+1)_{stat_name}"] = []
    for stat_name, val in stats.items():
        atomic_angle_features[f"O3'-P(+1)-O2'(+1)_{stat_name}"].append(val)

# 8. O5'(+1) - P(+1) - O2'(+1)  [CORRECTED: O5' at i+1]
for i in range(num_nucleotides - 1):
    angles = compute_angle_all_frames(
        O5p_coords[i + 1],
        P_coords[i + 1],
        O2p_coords[i + 1]
    )
    stats = compute_all_angle_stats(angles)
    if i == 0:
        for stat_name in stats.keys():
            atomic_angle_features[f"O5'(+1)-P(+1)-O2'(+1)_{stat_name}"] = []
    for stat_name, val in stats.items():
        atomic_angle_features[f"O5'(+1)-P(+1)-O2'(+1)_{stat_name}"].append(val)


# 9. C2 - C1'(+1) - P(+1)  [EXACT AS PUBLISHED - Cross-residue base-sugar-phosphate]
for i in range(num_nucleotides - 1):
    angles = compute_angle_all_frames(
        C2_coords[i],        # C2 base at residue i
        C1p_coords[i + 1],   # C1' sugar at residue i+1 (VERTEX)
        P_coords[i + 1]      # P phosphate at residue i+1
    )
    stats = compute_all_angle_stats(angles)
    if i == 0:
        for stat_name in stats.keys():
            atomic_angle_features[f"C2-C1'(+1)-P(+1)_{stat_name}"] = []
    for stat_name, val in stats.items():
        atomic_angle_features[f"C2-C1'(+1)-P(+1)_{stat_name}"].append(val)

# 10. O5' - P(+1) - O2'(+1)  [EXACT AS PUBLISHED - Cross-residue O5'-phosphate-2'OH]
for i in range(num_nucleotides - 1):
    angles = compute_angle_all_frames(
        O5p_coords[i],       # O5' at residue i
        P_coords[i + 1],     # P phosphate at residue i+1 (VERTEX)
        O2p_coords[i + 1]    # O2' at residue i+1
    )
    stats = compute_all_angle_stats(angles)
    if i == 0:
        for stat_name in stats.keys():
            atomic_angle_features[f"O5'-P(+1)-O2'(+1)_{stat_name}"] = []
    for stat_name, val in stats.items():
        atomic_angle_features[f"O5'-P(+1)-O2'(+1)_{stat_name}"].append(val)

# Pad all atomic angle features to num_nucleotides
for key in atomic_angle_features:
    values = atomic_angle_features[key]
    if len(values) < num_nucleotides:
        values = values + [np.nan] * (num_nucleotides - len(values))
    atomic_angle_features[key] = values

atomic_angles_df = pd.DataFrame(atomic_angle_features)

print(f"  Added 10 atomic angles × 6 statistics = 60 features")
print(f"  Total atomic angle features: {len(atomic_angles_df.columns)}")


# Combine distance and angle geometrical features
geo_complete_df = pd.concat([geo_distances_df, angles_combined_df, atomic_angles_df], axis=1)
geo_complete_df.insert(0, 'nucleotide', nt_ids)

n_geo_features = len(geo_complete_df.columns) - 1
print(f"\nGeometrical features computed: {n_geo_features} features")
print(f"  - Distance features: 14 types × 6 stats = 84")  # Added C1'-C1' and C2-C2
print(f"  - Angle features (centroids): 8 types × 6 stats = 48")
print(f"  - Angle features (atomic): 10 types × 6 stats = 60")
print(f"  - Density counts: 3")
print(f"  Total geometrical: {n_geo_features}")
print(f"\n  NEW: C1'-C1'(+1) and C2-C2(+1) features")
print(f"    - Expected to achieve r ~ 0.69-0.70 based on literature!")
print(f"    - C1': Sugar-sugar distance (backbone flexibility)")
print(f"    - C2: Base-base distance (base stacking/opening)")

# ============================================================
# MERGE ALL FEATURES
# ============================================================

conformational_df.insert(0, 'nucleotide', nt_ids)
torsion_df.insert(0, 'nucleotide', nt_ids)

ml_features_complete = (features_df
                        .merge(conformational_df, on='nucleotide')
                        .merge(torsion_df, on='nucleotide')
                        .merge(geo_complete_df, on='nucleotide'))

# Save to CSV
ml_features_complete.to_csv('dssr_ml_features_extended_statistics7.csv', index=False)

print(f"\n{'='*70}")
print(f"=== COMPLETE ML DATASET WITH EXTENDED STATISTICS ===")
print(f"{'='*70}")
print(f"Total: {ml_features_complete.shape[0]} nucleotides × {ml_features_complete.shape[1]-1} features")
print(f" - Structural: 9")
print(f" - Conformational: {len(conformational_df.columns)-1}")
print(f" - Torsion statistics: {len(torsion_df.columns)-1}")
print(f" - Geometrical (extended): {n_geo_features}")
print(f"\nOutput file: dssr_ml_features_extended_statistics.csv")
print(f"\nPer distance/angle, now computing 6 statistics:")
print(f"  - mean:  average value")
print(f"  - std:   standard deviation (fluctuation)")
print(f"  - var:   variance (fluctuation squared)")
print(f"  - ptp:   peak-to-peak (max - min range)")
print(f"  - p90:   90th percentile (high excursion events)")
print(f"  - cv:    coefficient of variation (std/mean)")
print(f"\n** All computed over {traj.n_frames} trajectory frames **")
print(f"{'='*70}")
