
# CiTCoM Modeling Team
---

## Introduction
 This repository includes the technical aspects of my postdoc work on 2025-2026. In here, you'll find input files, examples, tutorials, etc. to prepare your input files and analyze the outputs. Most of the scripts are optimized to run with "Adastra - CINES" HPC.
 
Before starting, make sure you have the following softwares installed on your computer:
- [Gromacs 2023.5](https://manual.gromacs.org/)
- [VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)
- [Chimera](https://www.cgl.ucsf.edu/chimera/download.html)
- [PLUMED](https://www.plumed.org/)

** Please note that PLUMED must be installed first, then Gromacs. 

---
## List of directories
0. Just a bunch of handy tools
1. Prepare and run plain MD in Gromacs.
2. Prepare and run Replica Exchange with Solute Scaling (REST2)
3. Prepare and run Umbrella sampling
4. A few tests regarding the ML approach

## Installing requirements
To prepare the MD setup, first the PLUMED must be installed. As of now, the latest version ([v2.10.x](http://plumed.github.io/doc-v2.10/user-doc/html/index.html)) is available on the official website.

### Plumed:
Download and extract the zip file. Run:

     ./configure --enable-modules=all --prefix=/path/to/install/dir --enable-shared
     make
     make install

Update your .bashrc with: `source /path/to/install/sourceme.sh`

If you get hwloc warnings, you can use [this approach](https://github.com/openwall/john/issues/5088) to hide them.

### Gromacs:
PLUMED works with handful number and versions of softwares. Make sure you choose the right version. In case of gromacs, in the latest version (v2.10.x), the followings are supported:
 -  gromacs-2022-5
-   gromacs-2023-5
-   gromacs-2024-3
-   gromacs-2025-0

After extracting the gromacs source:

    cd /path/to/gromacs/folder
    plumed patch -p

and you can continue with compiling it using the gormacs [installation guide](https://manual.gromacs.org/documentation/2023.5/install-guide/index.html).

## Useful commands:
gmx_mpi trjcat -f md-2.1.part0001.xtc etc.xtc -o total.xtc

