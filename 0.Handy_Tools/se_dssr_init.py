#!/home/saeid/.conda/envs/ml2/bin/python3

import subprocess
import argparse
import os
import tempfile
import shutil
import json
import pandas as pd
from Bio.Align import PairwiseAligner

def split_nt_id(nt_id):
    """Splits nt_id into (chain, name) e.g., 'B.G19' -> ('B', 'G19')."""
    return nt_id.split('.') if '.' in nt_id else ['', nt_id]

def extract_dssr_info(json_path):
    """Extracts sequence and structure information from DSSR JSON."""
    try:
        with open(json_path, 'r') as f:
            dssr_json = json.load(f)
        dbn = dssr_json.get('dbn', {})
        all_chains = dbn.get('all_chains', {})
        bseq = all_chains.get('bseq', 'N/A')
        sstr = all_chains.get('sstr', 'N/A')
        return bseq, sstr
    except Exception as e:
        print(f"Error extracting DSSR info: {e}")
        return 'N/A', 'N/A'

def generate_csv(json_path, csv_path):
    """Extracts nucleotide info from DSSR JSON and saves to CSV using pandas."""
    try:
        with open(json_path, 'r') as f:
            dssr_json = json.load(f)

        if 'nts' not in dssr_json:
            print("Warning: 'nts' key not found in JSON. CSV not generated.")
            return

        nts_list = []
        for nt in dssr_json['nts']:
            index = nt['index']
            # nt_id is usually like 'B.G1' or 'A.C25'
            full_nt_id = nt['nt_id']
            name = full_nt_id.split('.')[1] if '.' in full_nt_id else full_nt_id
            nts_list.append([index, name])

        df = pd.DataFrame(nts_list, columns=['nuc_index', 'nt_name'])
        df.set_index('nt_name', inplace=True)

        # Initialize new features
        df['is_in_helix'] = 0
        df['helix_size'] = 0
        df['is_in_stem'] = 0
        df['stem_size'] = 0

        # Extract Helix Information
        if 'helices' in dssr_json and dssr_json.get('num_helices', 0) > 0:
            for helix in dssr_json['helices']:
                h_size = helix['num_pairs']
                for pair in helix['pairs']:
                    nt1_name = split_nt_id(pair['nt1'])[1]
                    nt2_name = split_nt_id(pair['nt2'])[1]
                    for nt_name in [nt1_name, nt2_name]:
                        if nt_name in df.index:
                            df.at[nt_name, 'is_in_helix'] = 1
                            df.at[nt_name, 'helix_size'] = h_size

        # Extract Stem Information
        if 'stems' in dssr_json and dssr_json.get('num_stems', 0) > 0:
            for stem in dssr_json['stems']:
                s_size = stem['num_pairs']
                for pair in stem['pairs']:
                    nt1_name = split_nt_id(pair['nt1'])[1]
                    nt2_name = split_nt_id(pair['nt2'])[1]
                    for nt_name in [nt1_name, nt2_name]:
                        if nt_name in df.index:
                            df.at[nt_name, 'is_in_stem'] = 1
                            df.at[nt_name, 'stem_size'] = s_size

        df.to_csv(csv_path)
        print(f"Successfully generated CSV: {csv_path}")

    except Exception as e:
        print(f"Error generating CSV: {e}")

def main():
    parser = argparse.ArgumentParser(description="Run DSSR on a PDB file using a temporary sandbox and export CSV with pandas.")
    parser.add_argument("--pdb_id", required=True, help="PDB ID to use as prefix and output filename")
    parser.add_argument("--pdb_file", required=True, help="Path to the input PDB file")
    parser.add_argument("--no_csv", action="store_true", help="Disable CSV generation")
    parser.add_argument("--fasta", help="RNA sequence string to print at end")

    args = parser.parse_args()

    pdb_id = args.pdb_id
    pdb_file = os.path.abspath(args.pdb_file)
    
    if not os.path.exists(pdb_file):
        print(f"Error: File {pdb_file} not found.")
        return

    output_dir = os.path.dirname(pdb_file)
    final_json_path = os.path.join(output_dir, f"{pdb_id}.json")
    final_csv_path = os.path.join(output_dir, f"{pdb_id}.csv")

    # Create a temporary directory for DSSR intermediate files
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"Created temporary sandbox: {tmpdir}")
        
        # DSSR output file in the temp directory
        tmp_json_path = os.path.join(tmpdir, f"{pdb_id}.json")
        
        # Run DSSR inside the temporary directory
        cmd = [
            'x3dna-dssr', 
            f'--input={pdb_file}', 
            f'--output={tmp_json_path}', 
            '--json', 
            f'--prefix={pdb_id}'
        ]
        
        try:
            print(f"Running DSSR on {pdb_file}...")
            subprocess.run(cmd, check=True, cwd=tmpdir)
            
            # Move the final JSON back to the original directory
            if os.path.exists(tmp_json_path):
                shutil.move(tmp_json_path, final_json_path)
                print(f"Successfully generated and moved: {final_json_path}")
                
                # Generate CSV if requested
                if not args.no_csv:
                    generate_csv(final_json_path, final_csv_path)
            else:
                print("Error: DSSR did not generate the expected JSON file.")
                
        except subprocess.CalledProcessError as e:
            print(f"Error running DSSR: {e}")
        except FileNotFoundError:
            print("Error: x3dna-dssr command not found. Please ensure it's in your PATH.")

    if args.fasta:
        bseq, sstr = extract_dssr_info(final_json_path)
        print('-' * 20)
        print("DSSR-inferred sequence (bseq):", bseq)
        print("DSSR structure (dbn/sstr):", sstr)
        print('-' * 20)
        print("User-provided RNA sequence:", args.fasta)
        print('-' * 20)
        
        try:
            aligner = PairwiseAligner(match_score=1.0)
            alignments = aligner.align(args.fasta, bseq)
            print("Alignment Result:")
            print(alignments[0])
        except Exception as e:
            print(f"Error performing alignment: {e}")

if __name__ == "__main__":
    main()
