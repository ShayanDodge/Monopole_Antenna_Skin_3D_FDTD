#!/bin/sh
#PBS -l nodes=1:ppn=32
##PBS -j oe
#PBS -o $PBS_JOBID.out
#PBS -e $PBS_JOBID.err
#PBS -N jet1_10
#PBS -q batch
#PBS -l walltime=00:10:00
#Please config the Sanple_parfor.m for parallel compute 
export MATLAB14=/share/matlabInput/
export MATLAB18=/share/Application/MATLAB/R2018a/
export PATH=$PATH:/usr/local/bin
cd $PBS_O_WORKDIR

# modify this section
export MFILE=three_layer_dispersive_skin_fdtd_cpml_3D_sinesoidal.m

time $MATLAB18/bin/matlab < $MFILE &> ${PBS_JOBID}.Matout
