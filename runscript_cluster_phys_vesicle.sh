#!/bin/sh

# # Path for saving simulation data
#pathname=/home/khatami.physics.sharif/membrane_project/
pathname=/home/maryam/Membrane_Project/phys/
cd $pathname
## k = spring constant(membrane), R = initial radius of membrane, r = radius of chiral motion of active beads, phi = initial packing fraction of active beads, cluster_max*ppn = maximum number of nodes allowed = 24
k=500
R=10
r=5
phi=0.4
seed_begin_value=1
seed_end_value=10
cluster_max=4

## ############################################################################################
# # This part creates run scripts
for i in $(seq 1 $cluster_max)
do
  touch $pathname/$i.sh
done

## ############################################################################################
# # This part creates run scripts
for i in $(seq 1 $cluster_max)
do
  echo "#!/bin/bash" > $pathname/$i.sh
  printf  "\n">> $pathname/$i.sh
## REMEMBER!!! to correct the following line properly, "phi0.1" ->...??
  echo "#PBS -N R${R}_r${r}_phi0.1_sh$i" >> $pathname/$i.sh
  echo "#PBS -l nodes=1:ppn=2" >> $pathname/$i.sh
  echo "#PBS -m abe" >> $pathname/$i.sh
  echo "#PBS -M khatami@physics.sharif.edu" >> $pathname/$i.sh
  echo "#PBS -l walltime=192:00:00" >> $pathname/$i.sh
  echo "#PBS -l cput=1536:00:00" >> $pathname/$i.sh
#  echo "#PBS -q default" >> $pathname/$i.sh
  printf  "\n">> $pathname/$i.sh
  echo "cd \"\$PBS_O_WORKDIR\"" >> $pathname/$i.sh
  printf  "\n">> $pathname/$i.sh
  echo "seed_i=$[$seed_begin_value+$i-1]">> $pathname/$i.sh
  echo "seed_f=$seed_end_value">> $pathname/$i.sh
  printf  "\n">> $pathname/$i.sh
  echo "seed=$[$seed_begin_value+$i-1]">> $pathname/$i.sh
  echo "while [ \"\$seed\" -le \"\$seed_f\" ]" >> $pathname/$i.sh
  echo "do" >> $pathname/$i.sh
  echo "  mpirun -np 2 ./a.out ${k} ${R} ${r} ${phi} \$seed" >> $pathname/$i.sh
  echo "  seed=\$[\$seed+$cluster_max]" >> $pathname/$i.sh
  echo "done" >> $pathname/$i.sh
done

## ############################################################################################
# # This part submits jobs to the queue
for i in $(seq 1 $cluster_max)
do
  qsub $i.sh
done

## ############################################################################################
