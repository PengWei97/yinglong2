#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 64

mpiexec -n 64 ~/projects/yinglong/yinglong-opt -i old_tripleJunction.i

# /public3/home/scg6211/projects/yinglong/grainGrowth/new_old_moppf/old/old_tripleJunction.i