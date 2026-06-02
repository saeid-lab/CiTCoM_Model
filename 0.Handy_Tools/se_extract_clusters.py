#!/home/saeid/.conda/envs/gnn/bin/python3

import re
import argparse
import shutil
import numpy as np
import MDAnalysis as mda
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Extract cluster frames and sub-trajectories from ttclust output.")
    parser.add_argument("--ttclust_info", type=str, required=True, help="Path to the ttclust.txt information file.")
    parser.add_argument("--topology", type=str, required=True, help="Path to the topology file (e.g., .gro, .pdb).")
    parser.add_argument("--trajectory", type=str, required=True, help="Path to the full trajectory file (e.g., .xtc).")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Define and create output directory
    output_dir = Path("sub-trajectory")
    output_dir.mkdir(exist_ok=True)
    
    ttclust_file = Path(args.ttclust_info)
    if not ttclust_file.exists():
        print(f"Error: ttclust_info file '{args.ttclust_info}' not found.")
        return

    # Phase 1: Parse ttclust info and save frame lists
    print(f"Reading {ttclust_file}...")
    text = ttclust_file.read_text()
    
    # Regex to capture cluster number and all frames
    pattern = re.compile(r"cluster\s+(\d+).*?Frames\s*:\s*\[(.*?)\]", re.S)
    matches = pattern.findall(text)
    
    if not matches:
        print("No clusters found in the provided ttclust_info file.")
        return

    print(f"Found {len(matches)} clusters. Saving frame lists to {output_dir}...")
    
    cluster_frames = {}
    for cluster_id, frames_str in matches:
        frames = [int(f.strip()) for f in frames_str.split(",") if f.strip()]
        cluster_frames[cluster_id] = np.array(frames)
        
        # Save frame list to file
        frame_file = output_dir / f"cluster_{cluster_id}_frames.txt"
        with frame_file.open("w") as f:
            f.write("\n".join(map(str, frames)))
        print(f"  ✅ cluster_{cluster_id}_frames.txt saved ({len(frames)} frames)")

    # Phase 2: Extract sub-trajectories
    print(f"\nLoading trajectory: {args.trajectory} with topology: {args.topology}...")
    try:
        u = mda.Universe(args.topology, args.trajectory)
    except Exception as e:
        print(f"Error loading trajectory with MDAnalysis: {e}")
        return

    for cluster_id, frames in cluster_frames.items():
        # ttclust uses 1-based indexing; MDAnalysis uses 0-based indexing.
        # Shift frames by -1
        mda_frames = frames - 1
        
        # Filter out negative indices just in case (though unlikely with 1-based input)
        mda_frames = mda_frames[mda_frames >= 0]
        
        print(f"Extracting cluster {cluster_id} with {len(mda_frames)} frames...")
        out_xtc = output_dir / f"cluster_{cluster_id}.xtc"
        
        with mda.Writer(str(out_xtc), n_atoms=u.atoms.n_atoms) as W:
            for ts in u.trajectory[mda_frames]:
                W.write(u.atoms)
        
        print(f"  ✅ Saved {out_xtc}")

    # Phase 3: Copy topology file to sub-trajectory folder
    topology_src = Path(args.topology)
    if topology_src.exists():
        topology_dest = output_dir / topology_src.name
        shutil.copy2(topology_src, topology_dest)
        print(f"\n✅ Topology file '{topology_src.name}' copied to {output_dir}")
    else:
        print(f"\n⚠️ Warning: Topology file '{args.topology}' not found for copying.")

    print("\n🎉 All clusters processed successfully!")

if __name__ == "__main__":
    main()
