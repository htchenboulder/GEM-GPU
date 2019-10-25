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

#OpenACC option
#export PGI_ACC_TIME=1
#aprun -n 4 pgprof -o gem.%p.prof ./gem_main >& run.out
#nvprof option
#jsrun -n 24 -a 7 -g 1 -c 7 nvprof -s --openmp-profiling off --output-profile gem.%p.prof ./gem_main > run.out 2> run.err
#pure option
jsrun -n 384 -a 7 -g 1 -c 7 ./gem_main  > run.out 2>run.err
