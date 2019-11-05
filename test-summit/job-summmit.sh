#!/bin/sh

#BSUB -P phy122
#BSUB -W 0:30
#BSUB -nnodes 64
#BSUB -alloc_flags gpumps
#BSUB -J gem_test
#BSUB -o gem_test.%J
#BSUB -e gem_test.%J

cd $LS_SUBCWD

mkdir -p matrix
mkdir -p out
mkdir -p dump

#cp ~walkup/bin/numactl .
#OpenACC option
#export PGI_ACC_TIME=1
#export RANKS_PER_NODE=42
#aprun -n 4 pgprof -o gem.%p.prof ./gem_main >& run.out
##nvprof option: p means number of mpi, openmp-profiling off means only gpu timeline will be reported
#jsrun -n 96 -r 6 -a 7 -g 1 -c 7 nvprof -s --openmp-profiling off --output-profile gem.%p.prof ./gem_main > run.out 2> run.err
#pure option
jsrun -n 384 -r 6 -a 7 -g 1 -c 7 ./gem_main > run.out 2>run.err
#jsrun --smpiargs="-gpu" -n 24 -r 6 -a 7 -g 1 -c 7 ~walkup/bin/h6smt1.sh ./gem_main > run.out 2>run.err
#jsrun --smpiargs="-gpu" -n 8 -r 1 -a 42 -g 6 -c 42 ~walkup/bin/h6smt1.sh ./gem_main > run.out 2>run.err
