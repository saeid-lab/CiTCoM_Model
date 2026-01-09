import json
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', None)

pdb_id = '1C2X_C1'
json_url = '../Graph/C1-f1850-s3104/heavy_C1-f1850-s3104.json'
RMSF  = np.loadtxt('../1C2X_ttclust/clusters_xtc/rmsf_cluster_1.xvg', comments=('#', '@'))

exp_1m7_nomg_nosam = np.loadtxt(fname='../Graph/exp_data/PDB_1M7_noSAM_noMg_T30C_all_RNAPDB022.shape', usecols=1)
exp_dms_nomg_nosam = np.loadtxt(fname='../Graph/exp_data/PDB_DMS_noSAM_noMg_T30C_all_RNAPDB022.shape', usecols=1)


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
                
                
exp_1m7_nomg_nosam = np.where(exp_1m7_nomg_nosam == -999, np.nan, exp_1m7_nomg_nosam)
exp_dms_nomg_nosam = np.where(exp_dms_nomg_nosam == -999, np.nan, exp_dms_nomg_nosam)

df['shape_1m7'] = exp_1m7_nomg_nosam
df['dms'] = exp_dms_nomg_nosam
df['RMSF'] = RMSF[:,1:2]

df.to_csv(f'{pdb_id}_features.csv')
