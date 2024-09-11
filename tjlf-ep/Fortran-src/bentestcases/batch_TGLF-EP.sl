#!/bin/bash

#SBATCH -p medium
#SBATCH -n 3
#SBATCH -t 24:00:00
#SBATCH -o ./run.out
#SBATCH -e ./run.err
#SBATCH -J TGLF-EP_OMFIT
#SBATCH --export=all

module purge
module load atom
module unload gcc/8.x
module unload env/gcc8.x
module load env/gcc8.x
module swap env/gcc8.x env/intel2020

env > env.txt

mpirun -n 3 /home/towlej/.julia/dev/TJLF.jl/tjlf-ep/Fortran-src/TGLF-EP/TGLFEP_driver

