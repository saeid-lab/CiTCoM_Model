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
# pd.set_option('display.max_columns', None)

# === Managing the files ===

parser = argparse.ArgumentParser(description="Geometric feature extraction")

parser.add_argument("--pdb_id", required=True)
parser.add_argument("--pdb_file", required=True)
parser.add_argument("--json_file", required=True)
parser.add_argument("--rmsf_file", default=None)
parser.add_argument("--exp_1m7_file", default=None)
parser.add_argument("--exp_dms_file", default=None)

args = parser.parse_args()



pdb_id   = args.pdb_id
pdb_url  = args.pdb_file
json_url = args.json_file   # json file from DSSR analysis "x3dna-dssr --input={pdb_path} --output={json_path} --json --prefix={prefix_path}

RMSF_URL     = args.rmsf_file
Data_1M7_URL = args.exp_1m7_file
Data_DMS_URL = args.exp_dms_file
# === End of managing the files ===


traj = md.load(filename_or_filenames=pdb_url)

RMSF = None
try:
    RMSF = np.loadtxt(RMSF_URL, comments=('#', '@'))
except Exception as e:
    print("RMSF error:", type(e), e)

exp_1m7_nomg_nosam = None
try:
    exp_1m7_nomg_nosam = np.loadtxt(fname=Data_1M7_URL, usecols=1)
except Exception as e:
    print("1M7 error:", type(e), e)

exp_dms_nomg_nosam = None
try:
    exp_dms_nomg_nosam = np.loadtxt(fname=Data_DMS_URL, usecols=1)
except Exception as e:
    print("DMS error:", type(e), e)

# === Print some basics ===
print('Number of chains:', traj.n_chains,
      '\nNumber of residues:', traj.n_residues,
      '\nNumber of frames:', traj.n_frames,)

# # Get residues from the topology
# residues = list(traj.top.residues)
# num_nucleotides = len(residues)
# print(f'number of nucleotides:{num_nucleotides} \nnucleotides: {residues}')




# Define atom groups for calculating centroids. Other atoms won't be used.
nucleobase_atoms = {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"}  # Common base atoms
sugar_atoms = {"C1'", "C2'", "C3'", "C4'", "O4'"}  # Ribose atoms
phosphate_atoms = {"P"}  # Just the phosphate atom
O2p_atom = {"O2'"}
residues = list(traj.top.residues)
num_nucleotides = len(residues)


def compute_centroid(residue, traj, atom_names):
    """Compute geometric center for a given set of atoms in a residue."""
    atom_indices = [atom.index for atom in residue.atoms if atom.name in atom_names]
    positions = traj.xyz[:, atom_indices, :]*10   # nm → A
    return np.mean(positions, axis=1)  # Shape: (n_frames, 3)

# Compute nucleobase, sugar, phosphate,and O2' centroids for each residue
nucleobase_centroids = np.array([compute_centroid(res, traj, nucleobase_atoms) for res in residues])
sugar_centroids = np.array([compute_centroid(res, traj, sugar_atoms) for res in residues])
phosphate_centroids = np.array([compute_centroid(res, traj, phosphate_atoms) for res in residues])

O2p_coords = np.array([
    traj.xyz[:, [atom.index for atom in res.atoms if atom.name in O2p_atom], :] * 10
    if any(atom.name in O2p_atom for atom in res.atoms)
    else np.full((traj.n_frames, 1, 3), np.nan)
    for res in residues
]).squeeze(2)


# === Utility for pairwise distances ===
def distance(a, b):
    """Euclidean distance across all frames"""
    return np.linalg.norm(a - b, axis=2)  # (num_res-Δ, n_frames)

def shift_distance(arr1, arr2, shift):
    """Compute distance between residue i and i+shift."""
    return np.linalg.norm(arr1[:-shift] - arr2[shift:], axis=2)

# === Compute all distance features ===
features = {}

for k in [1, 2, 3]:
    features[f"RiboseCentroid_RiboseCentroid(+{k})"] = np.nanmean(shift_distance(sugar_centroids, sugar_centroids, k), axis=1)
    features[f"baseCentroid_baseCentroid(+{k})"] = np.nanmean(shift_distance(nucleobase_centroids, nucleobase_centroids, k), axis=1)
    features[f"baseCentroid_P(+{k})"] = np.nanmean(shift_distance(nucleobase_centroids, phosphate_centroids, k), axis=1)
    features[f"O2'_P(+{k})"] = np.nanmean(shift_distance(O2p_coords, phosphate_centroids, k), axis=1)
    
# === Compute local density counts (within 4, 6, 8 Å) ===
# Use all heavy atoms for simplicity
heavy_atom_coords = traj.xyz[:, [a.index for a in traj.top.atoms if a.element.symbol != 'H'], :] * 10
heavy_atom_resid = np.array([a.residue.index for a in traj.top.atoms if a.element.symbol != 'H'])

def count_within_cutoff(res_idx, cutoff):
    """Count heavy atoms within cutoff from this residue's sugar centroid (mean over frames)."""
    center = sugar_centroids[res_idx]  # (n_frames, 3)
    # Compute all distances for each frame
    dists = np.linalg.norm(heavy_atom_coords - center[:, None, :], axis=2)  # (n_frames, n_atoms)
    mask = heavy_atom_resid != res_idx  # exclude atoms from same residue
    counts = np.sum((dists < cutoff) & mask[None, :], axis=1)
    return np.mean(counts)

for cutoff in [4, 6, 8]:
    features[f"count_within_{cutoff}A"] = [count_within_cutoff(i, cutoff) for i in range(num_nucleotides)]
    
    
features_aligned = {}
for key, values in features.items():
    # if shorter than num_res, pad with NaN at the end
    if len(values) < num_nucleotides:
        values = list(values) + [np.nan]*(num_nucleotides - len(values))
    features_aligned[key] = values

# === Build final DataFrame ===
df_features = pd.DataFrame(features_aligned)
df_features.insert(0, "resid", [res.index for res in residues])
df_features.insert(1, "resname", [res.name for res in residues])



# ---------- Helper functions ----------

def vector_angle(v1, v2):
    """Return angle between vectors v1 and v2 in radians."""
    dot = np.einsum('ij,ij->i', v1, v2)
    norm = np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)
    cosang = np.clip(dot / np.maximum(norm, 1e-8), -1.0, 1.0)
    return np.degrees(np.arccos(cosang))

def get_angle(A, B, C):
    """Angle ABC (at vertex B) over frames."""
    return vector_angle(A - B, C - B)

def mean_angle_over_frames(A, B, C):
    """Return mean angle over all frames, ignoring NaNs."""
    return np.nanmean(get_angle(A, B, C))

def base_plane_normal(base_coords, base_atom_masks=None):
    """
    Compute geometric base normals for each residue and frame.
    If multiple frames exist, compute per-frame normals and average.
    base_coords: (n_res, n_frames, n_atoms, 3) or (n_res, n_frames, 3)
    """
    n_res = base_coords.shape[0]

    # If only centroids were given (no atom-level detail)
    if base_coords.ndim == 3 and base_coords.shape[1] == 1:
        return np.full((n_res, 3), np.nan)  # can't derive plane from one point
    elif base_coords.ndim == 3:  # centroids only, multi-frame trajectory
        # fallback to motion-based gradient
        diffs = np.gradient(base_coords, axis=1)
        v1 = diffs[:, :-1, :]
        v2 = np.roll(v1, -1, axis=1)
        normals = np.cross(v1.mean(axis=1), v2.mean(axis=1))
        return normals / (np.linalg.norm(normals, axis=1, keepdims=True) + 1e-8)

    elif base_coords.ndim == 4:
        # True geometric normals from base atoms
        normals = []
        for i in range(n_res):
            frame_normals = []
            for f in range(base_coords.shape[1]):
                pts = base_coords[i, f]
                pts -= pts.mean(axis=0)
                u, s, vh = np.linalg.svd(pts)
                normal = vh[-1]
                frame_normals.append(normal)
            normals.append(np.mean(frame_normals, axis=0))
        normals = np.array(normals)
        return normals / (np.linalg.norm(normals, axis=1, keepdims=True) + 1e-8)


def compute_rna_angles(
    nucleobase_centroids,
    sugar_centroids,
    phosphate_centroids,
    O2p_coords,
    resnames=None
):
    """
    Compute key geometric angles for each RNA residue.
    Works for both single-frame PDBs and multi-frame trajectories.
    """

    # --- Ensure 3D shape (n_res, n_frames, 3) ---
    def ensure_3d(x):
        if x.ndim == 2:
            return x[:, np.newaxis, :]
        return x

    nucleobase_centroids = ensure_3d(nucleobase_centroids)
    sugar_centroids = ensure_3d(sugar_centroids)
    phosphate_centroids = ensure_3d(phosphate_centroids)
    O2p_coords = ensure_3d(O2p_coords)

    n_res, n_frames, _ = sugar_centroids.shape
    results = []

    # --- Compute base normals (geometric or fallback motion-based) ---
    base_normals = base_plane_normal(nucleobase_centroids)

    # --- Helper for angle averaging over frames ---
    def mean_angle_over_frames(a, b, c):
        # a,b,c each (n_frames, 3)
        ab = b - a
        cb = b - c
        ab /= np.linalg.norm(ab, axis=-1, keepdims=True) + 1e-8
        cb /= np.linalg.norm(cb, axis=-1, keepdims=True) + 1e-8
        cosang = np.clip(np.sum(ab * cb, axis=-1), -1.0, 1.0)
        ang = np.degrees(np.arccos(cosang))
        return np.nanmean(ang)

    # --- Main loop per residue ---
    for i in range(n_res):
        row = {
            "resid": i,
            "resname": resnames[i] if resnames is not None else None,
        }

        # 1. Sugar–Base–Phosphate (S–B–P)
        row["S-B-P"] = mean_angle_over_frames(
            sugar_centroids[i], nucleobase_centroids[i], phosphate_centroids[i]
        )

        # 2. O2'–P–(next)P (+1)
        if i < n_res - 1:
            row["O2'-P-P(+1)"] = mean_angle_over_frames(
                O2p_coords[i], phosphate_centroids[i], phosphate_centroids[i + 1]
            )
        else:
            row["O2'-P-P(+1)"] = np.nan

        # 3. O2'–P–Sugar
        row["O2'-P-S"] = mean_angle_over_frames(
            O2p_coords[i], phosphate_centroids[i], sugar_centroids[i]
        )

        # 4. O2'–S–B
        row["O2'-S-B"] = mean_angle_over_frames(
            O2p_coords[i], sugar_centroids[i], nucleobase_centroids[i]
        )
        
        # 5. Base-Phosphate-Base+1 (B-P-B+1)
        if i < n_res - 1:
            row["B-P-B+1"] = mean_angle_over_frames(
                nucleobase_centroids[i], phosphate_centroids[i], nucleobase_centroids[i+1]
            )
        else:
            row["B-P-B+1"] = np.nan

        # 6. Base-Phosphate-Base+2 (B-P-B+2)
        if i < n_res - 2:
            row["B-P-B+2"] = mean_angle_over_frames(
                nucleobase_centroids[i], phosphate_centroids[i], nucleobase_centroids[i+2]
            )
        else:
            row["B-P-B+2"] = np.nan
            
        # 7. Sugar-Phosphate-Sugar+1 (S-P-S+1)
        if i < n_res - 1:
            row["S-P-S+1"] = mean_angle_over_frames(
                sugar_centroids[i], phosphate_centroids[i], sugar_centroids[i+1]
            )
        else:
            row["S-P-S+1"] = np.nan
            
        # 8. Sugar-Phosphate-Sugar+1 (S-P-S+1)
        if i < n_res - 2:
            row["S-P-S+2"] = mean_angle_over_frames(
                sugar_centroids[i], phosphate_centroids[i], sugar_centroids[i+2]
            )
        else:
            row["S-P-S+2"] = np.nan
            
        results.append(row)

    return pd.DataFrame(results)

angles_df = compute_rna_angles(
    nucleobase_centroids,
    sugar_centroids,
    phosphate_centroids,
    O2p_coords,
    resnames=residues
)


# ---------- Merging distances and angles together ---------- 
angles_df_c = angles_df.iloc[:, 2:] 
df_merged = pd.concat([df_features, angles_df_c], axis=1)


# ---------- DSSR stuff ---------- 
with open(json_url) as f:
    dssr_json = json.load(f)

nts_list = []
for nt in dssr_json['nts']:
    index = nt['index']
    name = nt['nt_id'].split('.')[1]
    nts_list.append([index, name])

df = pd.DataFrame(nts_list, columns=['nuc_index', 'nt_name'])
df.set_index('nt_name', inplace=True)
df.insert(0, "pdb_id", pdb_id)

for base in ['A', 'U', 'C', 'G']:
    df[f'base_{base}'] = 0

# One-hot for bases
for idx, row in df.iterrows():
    base = idx[0]  # first letter of nt_name
    df.at[idx, f'base_{base}'] = 1
    
df['paired'] = 0
df['paired_to'] = -1
df['LW'] = None

# Fill pairing information
if 'pairs' in dssr_json and dssr_json['num_pairs'] > 0:
    for i in range(dssr_json['num_pairs']):
        nt = dssr_json['pairs'][i]
        nt1, nt2, LW = nt['nt1'], nt['nt2'], nt['LW']
        chain1, name1 = nt1.split('.')
        chain2, name2 = nt2.split('.')
        
        idx1 = df.at[name1, 'nuc_index']
        idx2 = df.at[name2, 'nuc_index']
        
        df.at[name1, 'paired'] = 1
        df.at[name1, 'paired_to'] = idx2
        df.at[name1, 'LW'] = LW
        
        df.at[name2, 'paired'] = 1
        df.at[name2, 'paired_to'] = idx1
        df.at[name2, 'LW'] = LW

# Ensure paired_to is numeric
df['paired_to'] = pd.to_numeric(df['paired_to'], errors='coerce').fillna(-1).astype(int)

# Fill LW NaN with 'other'
df['LW'] = df['LW'].fillna('other')

# --- One-hot LW features: cis/trans + edges ---
edges = ['WW', 'WH', 'WS', 'HW', 'HH', 'HS', 'SW', 'SH', 'SS', 'other']

# cis/trans
df['LW_cis'] = df['LW'].str.startswith('c').astype(int)
df['LW_trans'] = df['LW'].str.startswith('t').astype(int)

# edges
for edge in edges:
    df[f'LW_edge_{edge}'] = df['LW'].str[-2:].apply(lambda x: int(x == edge))

# Drop original LW column
df = df.drop(columns=['LW'])


df['is_in_multiplet'] = 0
df['multiplet_ids'] = 0
df['multiplet_size'] = 0

if 'multiplets' in dssr_json and dssr_json['num_multiplets'] > 0:
    for multiplet in dssr_json['multiplets']:
        m_id = multiplet['index']
        nts_long = multiplet['nts_long'].split(',')
        m_size = multiplet['num_nts']
    
        for nt_full in nts_long:
            nt_name = nt_full.split('.')[1]
    
            df.at[nt_name, 'is_in_multiplet'] = 1
            df.at[nt_name, 'multiplet_size'] = m_size
            # df.at[nt_name, 'multiplet_ids'] = m_id
            prev = df.at[nt_name, 'multiplet_ids']
            df.at[nt_name, 'multiplet_ids'] = f"{prev},{m_id}" if prev else str(m_id)
            
df['is_in_helix'] = 0
df['helix_id'] = 0
df['helix_size'] = 0

if 'helices' in dssr_json and dssr_json['num_helices'] > 0:
    # Loop over helices and fill info
    for helix in dssr_json['helices']:
        h_id = helix['index']
        h_size = helix['num_pairs']
        
        for pair in helix['pairs']:
            nt1_name = pair['nt1'].split('.')[1]
            nt2_name = pair['nt2'].split('.')[1]
            
            # Update nt1
            df.at[nt1_name, 'is_in_helix'] = 1
            df.at[nt1_name, 'helix_id'] = h_id
            df.at[nt1_name, 'helix_size'] = h_size
            
            # Update nt2
            df.at[nt2_name, 'is_in_helix'] = 1
            df.at[nt2_name, 'helix_id'] = h_id
            df.at[nt2_name, 'helix_size'] = h_size
            
            
df['is_in_stem'] = 0
df['stem_id'] = 0
df['stem_size'] = 0

if 'stems' in dssr_json and dssr_json['num_stems'] > 0:
    for stem in dssr_json['stems']:
        stem_id = stem['index']
        stem_size = stem['num_pairs']

        for pair in stem['pairs']:
            nt1 = pair['nt1'].split('.')[1]
            nt2 = pair['nt2'].split('.')[1]

            for nt_name in [nt1, nt2]:
                df.at[nt_name, 'is_in_stem'] = 1
                df.at[nt_name, 'stem_size'] = stem_size

                prev = df.at[nt_name, 'stem_id']
                df.at[nt_name, 'stem_id'] = f"{prev},{stem_id}" if prev else str(stem_id)
                
df['is_in_stack'] = 0
df['stack_id'] = 0
df['stack_size'] = 0

if 'stacks' in dssr_json and dssr_json['num_stacks'] > 0:
    for stack in dssr_json['stacks']:
        stack_id = stack['index']
        stack_size = stack['num_nts']
        
        for nt_full in stack['nts_long'].split(','):
            nt_name = nt_full.split('.')[1]

            df.at[nt_name, 'is_in_stack'] = 1
            df.at[nt_name, 'stack_size'] = stack_size

            prev = df.at[nt_name, 'stack_id']
            df.at[nt_name, 'stack_id'] = f"{prev},{stack_id}" if prev else str(stack_id)
            
            
df['is_in_hairpin'] = 0
df['hairpin_id'] = 0
df['hairpin_size'] = 0
df['hairpin_bridging_nt'] = 0
df['hairpin_bridging_id'] = 0


if 'hairpins' in dssr_json and dssr_json['num_hairpins'] > 0:
    for hpin in dssr_json['hairpins']:
        hairpin_id = hpin['index']
        hairpin_size = hpin['num_nts']

        # Loop over all nucleotides in hairpin
        for nt_full in hpin['nts_long'].split(','):
            nt_name = nt_full.split('.')[1]

            df.at[nt_name, 'is_in_hairpin'] = 1
            df.at[nt_name, 'hairpin_id'] = hairpin_id
            df.at[nt_name, 'hairpin_size'] = hairpin_size

        # Loop over bridging nucleotides
        for bridge in hpin.get('bridges', []):
            bridge_id = bridge['index']
            for nt_full in bridge['nts_long'].split(','):
                nt_name = nt_full.split('.')[1]
                if nt_name not in df.index:
                    continue
                df.at[nt_name, 'hairpin_bridging_nt'] = 1
                df.at[nt_name, 'hairpin_bridging_id'] = bridge_id
                
                
df['is_in_bulge'] = 0
df['bulge_id'] = 0
df['bulge_size'] = 0
df['bulge_bridging_nt'] = 0
df['bulge_bridging_id'] = 0

if 'bulges' in dssr_json and dssr_json['num_bulges'] > 0:
    for bulge in dssr_json['bulges']:
        bulge_id = bulge['index']
        bulge_size = bulge['num_nts']

        # Loop over all nucleotides in the bulge
        for nt_full in bulge['nts_long'].split(','):
            nt_name = nt_full.split('.')[1]
            
            df.at[nt_name, 'is_in_bulge'] = 1
            df.at[nt_name, 'bulge_id'] = bulge_id
            df.at[nt_name, 'bulge_size'] = bulge_size

        # Loop over bridging nucleotides inside the bulge
        for bridge in bulge.get('bridges', []):
            bridge_id = bridge['index']
            for nt_full in bridge.get('nts_long', '').split(','):
                if not nt_full:  # skip empty strings
                    continue
                nt_name = nt_full.split('.')[1]
                if nt_name not in df.index:
                    continue
                df.at[nt_name, 'bulge_bridging_nt'] = 1
                df.at[nt_name, 'bulge_bridging_id'] = bridge_id
                
                
df['is_in_junction'] = 0
df['junction_id'] = 0
df['junction_size'] = 0
df['junction_bridging_nt'] = 0
df['junction_bridging_id'] = 0

if 'junctions' in dssr_json and dssr_json['num_junctions'] > 0:
    for junction in dssr_json['junctions']:
        junction_id = junction['index']
        junction_size = junction['num_nts']

        # All nucleotides in the junction
        for nt_full in junction['nts_long'].split(','):
            nt_name = nt_full.split('.')[1]

            df.at[nt_name, 'is_in_junction'] = 1
            df.at[nt_name, 'junction_id'] = junction_id
            df.at[nt_name, 'junction_size'] = junction_size

        # Bridging nucleotides inside the junction
        for bridge in junction.get('bridges', []):
            bridge_id = bridge['index']
            for nt_full in bridge.get('nts_long', '').split(','):
                if not nt_full:  # skip empty strings
                    continue
                nt_name = nt_full.split('.')[1]
                if nt_name not in df.index:
                    continue
                df.at[nt_name, 'junction_bridging_nt'] = 1
                df.at[nt_name, 'junction_bridging_id'] = bridge_id
                

df['is_in_ssSegment'] = 0
df['ssSegment_id'] = 0
df['ssSegment_size'] = 0

if 'ssSegments' in dssr_json and dssr_json['num_ssSegments'] > 0:
    for ssSegment in dssr_json['ssSegments']:
        ssSegment_id = ssSegment['index']
        ssSegment_size = ssSegment['num_nts']

        for nt_full in ssSegment['nts_long'].split(','):
            nt_name = nt_full.split('.')[1]
            if nt_name not in df.index:
                continue

            df.at[nt_name, 'is_in_ssSegment'] = 1
            df.at[nt_name, 'ssSegment_id'] = ssSegment_id
            df.at[nt_name, 'ssSegment_size'] = ssSegment_size
            
            
torsions = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'epsilon_zeta', 'bb_type', 'chi', 'glyco_bond']
puckers = ['v0', 'v1', 'v2', 'v3', 'v4', 'amplitude', 'phase_angle', 'puckering']
angles = ['eta', 'theta', 'eta_prime', 'theta_prime', 'eta_base', 'theta_base']
base_stacks = ['ssZp', 'Dp', 'splay_angle', 'splay_distance', 'splay_ratio']

for col in torsions + puckers + angles + base_stacks:
    if col in ['bb_type', 'glyco_bond', 'puckering']:
        df[col] = None  # or pd.NA for nullable string
    else:
        df[col] = np.nan


for nt in dssr_json['nts']:
    nt_name = nt['nt_id'].split('.')[1]
    for col in torsions + puckers + angles + base_stacks:
        df.at[nt_name, col] = nt.get(col)
        
df['bb_type'] = df['bb_type'].replace('--', 'unknown')
df = pd.get_dummies(df, columns=['bb_type'], prefix='bb', dtype=int)

df = pd.get_dummies(df, columns=['glyco_bond'], prefix='glyco', dtype=int)

df['puckering_c3-endo'] = (df['puckering'] == "C3'-endo").astype(int)
df['puckering_c2-endo'] = (df['puckering'] == "C2'-endo").astype(int)
df['puckering_others'] = (~df['puckering'].isin(["C3'-endo", "C2'-endo"])).astype(int)
df = df.drop('puckering', axis=1)


df['num_standard_hbonds'] = 0
df['num_acceptable_hbonds'] = 0
df['num_questionable_hbonds'] = 0

if 'hbonds' in dssr_json and dssr_json['num_hbonds'] > 0:
    for hbond in dssr_json['hbonds']:
        hbond_type = hbond.get('donAcc_type', 'unknown')
        nt1 = hbond['atom1_id'].split('@')[1].split('.')[1]
        nt2 = hbond['atom2_id'].split('@')[1].split('.')[1]

        for nt_name in [nt1, nt2]:
            if hbond_type == 'standard':
                df.at[nt_name, 'num_standard_hbonds'] += 1
            elif hbond_type == 'acceptable':
                df.at[nt_name, 'num_acceptable_hbonds'] += 1
            elif hbond_type == 'questionable':
                df.at[nt_name, 'num_questionable_hbonds'] += 1
                
                
if exp_1m7_nomg_nosam is not None:
    exp_1m7_nomg_nosam = np.where(exp_1m7_nomg_nosam == -999, np.nan, exp_1m7_nomg_nosam)
    df['shape_1m7'] = exp_1m7_nomg_nosam
    
if exp_dms_nomg_nosam is not None:
    exp_dms_nomg_nosam = np.where(exp_dms_nomg_nosam == -999, np.nan, exp_dms_nomg_nosam)
    df['dms'] = exp_dms_nomg_nosam
    
if RMSF is not None:
    df['RMSF'] = RMSF[:,1:2]
    

# ---------- Merging with previous df ----------
df_merged_df_c = df_merged.iloc[:, 2:] 

df_merged_all = pd.concat([
    df.reset_index(drop=True),
    df_merged_df_c.reset_index(drop=True)
], axis=1)

    
# Save the CSV file
output_filename = f'{pdb_id}_allfeatures.csv'
df_merged_all.to_csv(output_filename)
print(f"  Successfully saved: {output_filename}")
