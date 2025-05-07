## Requied files for preparing the simulation:
1. amber14sb_OL15.ff (Force Field folder)
2. min.mdp (Minimzation input file for step 9)

**top** directory: we prepare the gromacs input files (gro,top, etc.) (Remmber to Always Log). 
**mdp_single** : mdp files are gromacs settings for running MD.
**md** : first replica of MD for the system.
**rep_2** and **rep_3**: second and third replicas.


### Example of executed commands for 1DUL:
1. Downloading 1DUL.pdb
2. Removing "standard amino acids" and "CCC" residues in Chimera > 1dul_cleaned.pdb
3. awk 'NR >= 1643 && NR <= 1941 { $0 = substr($0, 1, 13) "OW" substr($0, 16, 2) "SOL" substr($0, 21, 1) "C" substr($0, 23) } { print }' 1dul_cleaned.pdb > temp.pdb && mv temp.pdb 1dul_cleaned_OW.pdb
4. gmx pdb2gmx -f 1dul_cleaned_OW.pdb -ignh  # 1: OL15 nucleic and 3: TIP4P-Ew
5. gmx editconf -f conf.gro -o 1dul_box.gro -bt dodecahedron -d 1.5 -c
6. gmx solvate -cp 1dul_box.gro -cs tip4p.gro -o 1dul_sol.gro -p topol.top
7. gmx grompp -f min.mdp -c 1dul_sol.gro -r 1dul_sol.gro -p topol.top -o genion_input.tpr -maxwarn 1
8. gmx genion -s genion_input.tpr -p topol.top -o 1dul_sol_kcl.gro -pname K -pq 1 -np 100 -nname CL -nq -1  -nn 64 # Group 6 - SOL **calc** 23852 * 0.15 / 55.5 = 64
9. gmx grompp -f min.mdp -c 1dul_sol_kcl.gro -p topol.top -o input_min.tpr -maxwarn 1 -r 1dul_sol_kcl.gro
10. gmx mdrun -s input_min.tpr -deffnm min -v

11. gmx make_ndx -f min.gro  # Enter 1 & ! a H* (1 for RNA then heavy atom)
```
12. gmx genrestr -f conf.gro  -fc  1000 1000 1000 -o posres1_RNA.itp  -n index.ndx # Select 10 for RNA_&_!H*
12. gmx genrestr -f conf.gro  -fc  500 500 500 -o posres500_RNA.itp  -n index.ndx # Select 10 for RNA_&_!H*
12. gmx genrestr -f conf.gro  -fc  100 100 100 -o posres100_RNA.itp  -n index.ndx # Select 10 for RNA_&_!H*
12. gmx genrestr -f conf.gro  -fc  10 10 10 -o posres10_RNA.itp  -n index.ndx # Select 10 for RNA_&_!H*
```

### Then add all the posres*_RNA.itp to topol.top (OR topol_RNA_chain_B.itp) file, like the example below:
```
; Include Position restraint file
#ifdef POSRES1
#include "posres1_RNA.itp"
#endif

; Include Position restraint file
#ifdef POSRES500
#include "posres500_RNA.itp"
#endif

; Include Position restraint file
#ifdef POSRES100
#include "posres100_RNA.itp"
#endif

; Include Position restraint file
#ifdef POSRES10
#include "posres10_RNA.itp"
#endif

etc.
```

### Copy min.gro from top dir to md.
12. Upload all the folders to Adastra.
13. Go to md folder and run this script: auto_md_eq.sh
14. Once it is finished in adastra: 
- module load  .archive/CCE-GPU-3.0.0 gromacs/2023_amd-mpi-omp-plumed-python3
- gmx_mpi make_ndx -f min.gro # just q
- gmx_mpi grompp -f ./../mdp_single/md-2.mdp -c md_eq6.gro -p ../top/topol.top -o md-2 -n index -maxwarn 1 
- sbatch job_restart.sh

### for Replica 2 and 3:
- cd rep_2 && sbatch eq_rep.sh --> once finished sbatch job_restart.sh 
- cd rep_3 && sbatch eq_rep.sh --> once finished sbatch job_restart.sh 
