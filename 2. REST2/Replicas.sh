# number of replicas
nrep=24

# "effective" temperature range
tmin=300
tmax=1000

# build geometric progression
list=$(
awk -v n=$nrep \
    -v tmin=$tmin \
    -v tmax=$tmax \
  'BEGIN{for(i=0;i<n;i++){
    t=tmin*exp(i*log(tmax/tmin)/(n-1));
    printf(t); if(i<n-1)printf(",");
  }
}'
)

# clean directory
rm -fr \#*
rm lambdas.txt
#rm -fr topol*

for((i=0;i<nrep;i++))
do

# choose lambda as T[0]/T[i]
# remember that high temperature is equivalent to low lambda
  lambda=$(echo $list | awk 'BEGIN{FS=",";}{print $1/$'$((i+1))';}')
  echo 'writing topology for lambda' $lambda 'in folder:'$i 
  echo $lambda >> lambdas.txt
# process topology
# (if you are curious, try "diff topol0.top topol1.top" to see the changes)
  mkdir -p $i
  cp eq_rest2.sh $i
  cp plumed.dat $i
  cd $i
  plumed partial_tempering $lambda < ../processed.top > topol.top
# prepare tpr file
# -maxwarn is often needed because box could be charged
# gmx grompp -c ../md_eq6.gro -o topol$i.tpr -f ../rest2_eq.mdp -p topol$i.top -maxwarn 2
  
  cd ..
done

#  mpirun -np $nrep gmx_mpi mdrun_d -v -plumed plumed.dat -multi $nrep -replex 100 -nsteps 15000000 -hrex -dlb no
