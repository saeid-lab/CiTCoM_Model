### Working with files and directories:

- Get dir(s) size:
```
du -h --max-depth=1 | sort -hr
du -sh dir_path
```
Give you a list of directories and their corresponding size.

- Make a zip file: (with and without compression)
```
tar -czvf $STOREDIR/PATH/To/File.tar.gz Dir_Name
tar -cvf $STOREDIR/PATH/To/File.tar Dir_Name
```
Create File.tar.gz in the the mentioned path from Dir_name

- Extract just file names (for heavy and big zip files):
```
tar ztvf File.tar.gz > list_prot_RNA.txt
```
If the archive is valid, this will silently complete. If there is corruption, you'll get an error like gzip: stdin: unexpected end of file.

- Sort the xtc files in order:
```
ls -v md-2*.xtc | tr '\n' ' '
```

- Get size of the .xtc files in the current dir:
```
du -ch *.xtc | tail -n 1
```

- Attaching partial xtc files:
```
module load  .archive/CCE-GPU-3.0.0 gromacs/2023_amd-mpi-omp-plumed-python3
gmx_mpi trjcat -f md-2.1.part0001.xtc md-2.2.part0002.xtc md-2.3.part0003.xtc -o mm.xtc
```
