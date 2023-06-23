# mpiexec -np 30 ~/projects/yinglong/yinglong-opt -i test0.i 

# tar -cvf - ex_test2/*.e.* ex_test2/*.e-s???2.* | pigz -9 -p 20 > ex_test2.tgz

# METHOD=dbg make -j30

# mpiexec -n 30 ~/projects/yinglong/yinglong-dbg -i polycrystal_2.i --stop-for-debugger 
# --start-in-debugger='gdb'

# gdb ~/projects/yinglong/yinglong-opt -i inputfile.i

# mamba activate moose
# mamba update --all
# make clobberall
# make -j30