#!/bin/sh

#PBS -N gem_gpu
#PBS -A fus123
#PBS -l nodes=8
#PBS -l walltime=01:00:00
#PBS -q debug

cd $PBS_O_WORKDIR

mkdir -p matrix
mkdir -p out
mkdir -p dump

export PGI_ACC_TIME=1
export OMP_NUM_THREADS=8
#aprun -n 4 pgprof -o gem.%p.prof ./gem_main >& run.out
aprun -n 8 -d 8 -j 1 ./gem_main >& run.out
