#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 2
#SBATCH -n 128

mpiexec -n 128 ~/projects/yinglong/yinglong-opt -i test1_gbIsotropy.i