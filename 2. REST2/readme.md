# Requied files for running a simulation:
1. amber14sb_OL15.ff (Force Field folder)	(for generating processed.top)
2. topol.top and its dependencies (.itp files)	(for generating processed.top)
3. Last equilibrium gro file: e.g. md_eq6.gro	(for generating processed.top)
4. HotAtomIndicator.sh	(adding _ for atoms)
5. rest2_eq.mdp		(mdp file for a short equilibrium)
6. Replicas.sh		(creating 24 replica folders with necessary data)
7. eq_rest2.sh		(slurm job for the short equilibrium)
8. REST2.sh		(replica exchange process)
9. md-2-rest2.mdp	(mdp file for replica exchange)
10. empty plumed.dat	(is required)

# produce a processed topology
gmx grompp -f rest2_eq.mdp -c md_eq6.gro -p topol.top -o topol-unscaled.tpr -pp processed.top -maxwarn 1


# choose the "hot" atoms by appending "_" in the atom names. You can choose whatever.
- vim processed.top
- OR using the HotAtomIndicator.sh script. Just choose the first line/last lines and file names... (bash HotAtomIndicator.sh)

# generate the actual topology scaled with factor 1.0 untill...
bash Replicas.sh
``` 
It does create N replicas with different lambda -> different topol.top
Then copies the eq_rest2.sh and plumed.dat to each folder.
```

# Upload REST2 folder to adastra
# Launching a short equilibrium
bash auto_eq_rest2.sh


# Launching production
sbatch REST2.sh (for single run 24h)
sbatch REST2_restart.sh (for restart job)
