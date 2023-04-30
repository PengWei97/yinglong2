#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 2
#SBATCH -n 128

# mpiexec -n 128 ~/projects/yinglong/yinglong-opt -i test5_stored_gbAniso.i --recover
mpiexec -n 128 ~/projects/yinglong/yinglong-opt -i test5_stored_gbAniso.i

# tar -cvf - ex_case5_1*/*.e.* ex_case5_1*/*.e-s0??[258].* | pigz -9 -p 20 > ex_case5_1.tgz
