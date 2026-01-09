# CiTCoM Modeling Team
---

## Introduction
 In this repository, you'll find necessary input files, examples, tutorials, etc. to prepare your input files and analyze the outputs.
 
Before starting, make sure you have the following softwares installed on your computer:
- [Gromacs 2023.5](https://manual.gromacs.org/)
- [VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)
- [Chimera](https://www.cgl.ucsf.edu/chimera/download.html)
- [PLUMED](https://www.plumed.org/)

** Please note that PLUMED must be installed first, then Gromacs. 

---
## List of directories
1. Prepare and run plain MD in Gromacs.
2. Prepare and run Replica Exchange with Solute Scaling (REST2)
3. Prepare and run Umbrella sampling



## Useful commands:
gmx_mpi trjcat -f md-2.1.part0001.xtc etc.xtc -o total.xtc
