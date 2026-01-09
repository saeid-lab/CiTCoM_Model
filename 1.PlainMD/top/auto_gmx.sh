#!/bin/bash
echo "_______       _____            ______________  ______  __"
echo "___    |___  ___  /______      __  ____/__   |/  /_  |/ /"
echo "__  /| |  / / /  __/  __ \\     _  / __ __  /|_/ /__    / "
echo "_  ___ / /_/ // /_ / /_/ /     / /_/ / _  /  / / _    |  "
echo "/_/  |_\\__,_/ \\__/ \\____/      \\____/  /_/  /_/  /_/|_|  "
echo "                                                         "
echo "                                                         "
echo "    welcome to GMX auto script"
echo "         By Saeid E"

read -p "Are you ready to proceed? (y/n): " -n 1 -r REPLY
echo

if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    echo "Script canceled by user."
    exit 1 
fi

read -p "What is this RNA PDB ID? " pdbid

echo "Proceeding with the script..."

echo "##### Defining a simulation box #####"
sleep 3
gmx editconf -f conf.gro -o "${pdbid}_box.gro" -bt dodecahedron -d 1.5 -c > log_autogmx.log 2>&1

echo "##### Adding Water molecules to the simulation box #####"
sleep 3
gmx solvate -cp "${pdbid}_box.gro" -cs tip4p.gro -o "${pdbid}_sol.gro" -p topol.top >> log_autogmx.log 2>&1

echo "##### Adding Ions to the simulation box #####"
sleep 3
gmx grompp -f min.mdp -c "${pdbid}_sol.gro" -r "${pdbid}_sol.gro" -p topol.top -o genion_input.tpr -maxwarn 1 >> log_autogmx.log 2>&1

echo "##### Capturing number of water molecules and ions #####"
sleep 3


total_charge=$(grep "System has non-zero total charge" log_autogmx.log | awk '{printf "%.0f\n", $NF}')
n_waters=$(grep "Water residues" log_autogmx.log | awk '{print $3}')

# Compute c_ions = n_waters * 0.15 / 55.5
c_ions=$(awk -v n="$n_waters" 'BEGIN { printf "%.0f", (n * 0.15) / 55.5 }')

# Compute absolute charge
abs_charge=$(awk -v c="$total_charge" 'BEGIN { printf "%.0f", (c < 0 ? -c : c) }')


echo "System Total Charges: $total_charge"
echo "Number of Water Molecules: $n_waters"
echo "Number of Calculated Ions: $c_ions"
echo "Absolute charge: $abs_charge"

echo "##### Adding ions to the simulation box #####"
sleep 5
# Run gmx genion with appropriate parameters
if (( total_charge < 0 )); then
    total_pos=$((c_ions + abs_charge))
    echo 'SOL' | gmx genion -s genion_input.tpr -p topol.top -o "${pdbid}_sol_kcl.gro" \
        -pname K -pq 1 -np "$total_pos" \
        -nname CL -nq -1 -nn "$c_ions" \
        >> log_autogmx.log 2>&1
        
elif (( total_charge > 0 )); then
    total_neg=$((c_ions + abs_charge))
     echo 'SOL' | gmx genion -s genion_input.tpr -p topol.top -o "${pdbid}_sol_kcl.gro" \
        -pname K -pq 1 -np "$c_ions" \
        -nname CL -nq -1 -nn "$total_neg" \
        >> log_autogmx.log 2>&1
else
    # Neutral charge case
     echo 'SOL' | gmx genion -s genion_input.tpr -p topol.top -o "${pdbid}_sol_kcl.gro" \
        -pname K -pq 1 -np "$c_ions" \
        -nname CL -nq -1 -nn "$c_ions"\
        
        >> log_autogmx.log 2>&1
fi

echo "##### System is now neutralized #####"
sleep 3

echo "##### Creating minimzation input file #####"
sleep 3
gmx grompp -f min.mdp -c "${pdbid}_sol_kcl.gro" -p topol.top -o input_min.tpr -maxwarn 1 -r "${pdbid}_sol_kcl.gro" >> log_autogmx.log 2>&1




read -p "Are you ready to launch the minimization? (y/n): " -n 1 -r REPLYY
echo

if [[ ! $REPLYY =~ ^[Yy]$ ]]
then
    echo "Script canceled by user."
    exit 1 
fi
echo "Proceeding with the script... minimization"
sleep 3
gmx mdrun -s input_min.tpr -deffnm min -v


echo "##### Minimization is finished #####"
echo "##### Creating index file for RNA's heavy atoms (1 & ! a H* ) #####"
sleep 3

echo -e '1 & ! a H*\nq' | gmx make_ndx -f min.gro  >> log_autogmx.log 2>&1


echo "##### Creating position restraint files #####"
sleep 3
echo 'RNA_&_!H*' | gmx genrestr -f conf.gro  -fc  1000 1000 1000 -o posres1_RNA.itp  -n index.ndx >> log_autogmx.log 2>&1
echo 'RNA_&_!H*' | gmx genrestr -f conf.gro  -fc  500 500 500 -o posres500_RNA.itp  -n index.ndx >> log_autogmx.log 2>&1
echo 'RNA_&_!H*' | gmx genrestr -f conf.gro  -fc  100 100 100 -o posres100_RNA.itp  -n index.ndx >> log_autogmx.log 2>&1
echo 'RNA_&_!H*' | gmx genrestr -f conf.gro  -fc  10 10 10 -o posres10_RNA.itp  -n index.ndx >> log_autogmx.log 2>&1


read -p "Enter the name of the RNA topology file (.top or .itp): " topo_file
sleep 3

cat << 'EOF' >> "$topo_file"

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

EOF

echo "Position restraint includes appended to $topo_file"
sleep 3

echo "done"

cp min.gro ../md
# You can then use these variables for further processing in your script
# For example:
# if [[ -n "$total_charge" ]]; then
#   echo "The system charge is $total_charge."
# fi
# https://gromacs.org-gmx-users.maillist.sys.kth.narkive.com/gFWFJd1K/gmx-users-how-to-automate-genion-completely
