PDB_ID='1NBS'

echo -e 'RNA\n' | gmx trjconv -f md-"$PDB_ID"_Rep1.xtc -o md-"$PDB_ID"_Rep1_whole.xtc -s md-2.tpr -n index.ndx -pbc whole
echo -e 'RNA\n' | gmx trjconv -f md-"$PDB_ID"_Rep2.xtc -o md-"$PDB_ID"_Rep2_whole.xtc -s md-2.tpr -n index.ndx -pbc whole
echo -e 'RNA\n' | gmx trjconv -f md-"$PDB_ID"_Rep3.xtc -o md-"$PDB_ID"_Rep3_whole.xtc -s md-2.tpr -n index.ndx -pbc whole

echo -e 'RNA\nRNA\n' | gmx trjconv -f md-"$PDB_ID"_Rep1_whole.xtc -o md-"$PDB_ID"_Rep1_fit.xtc -s md-2.tpr -fit rot+trans -n index.ndx
echo -e 'RNA\nRNA\n' | gmx trjconv -f md-"$PDB_ID"_Rep2_whole.xtc -o md-"$PDB_ID"_Rep2_fit.xtc -s md-2.tpr -fit rot+trans -n index.ndx
echo -e 'RNA\nRNA\n' | gmx trjconv -f md-"$PDB_ID"_Rep3_whole.xtc -o md-"$PDB_ID"_Rep3_fit.xtc -s md-2.tpr -fit rot+trans -n index.ndx

echo -e 'c\nc\nc\n' | gmx trjcat -f md-"$PDB_ID"_Rep1_fit.xtc md-"$PDB_ID"_Rep2_fit.xtc md-"$PDB_ID"_Rep3_fit.xtc -o md-"$PDB_ID"_whole_fit_Rep123.xtc -cat -settime

echo -e 'RNA\n' | gmx trjconv -f min.gro -s md-2.tpr -n index.ndx -pbc whole -o min_RNA.gro
echo -e 'RNA\n' | gmx convert-tpr -s md-2.tpr -n index.ndx -o md-2_RNA.tpr
echo -e 'RNA_&_!H*\n' | gmx convert-tpr -s md-2.tpr -n index.ndx -o md-2_RNA_heavy.tpr


 #Choosing only RNA(1)
ttclust -f md-"$PDB_ID"_whole_fit_Rep123.xtc -t min_RNA.gro  -sr "all" -sa "all" -s 100 -aa n -l "$PDB_ID"_ttclust.txt



# x3dna-dssr --input= --json --output=4ue5.json --prefix=4ue5

