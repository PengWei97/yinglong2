mpiexec -np 38 ~/projects/yinglong/yinglong-opt -i test4_ebsd_stored_GBAniso.i > 01.log

# ~/projects/yinglong/yinglong-dbg -i test2_GBAnisotropyMisori_vt.i
# tar -cvf - ex_gbAnisotropyMisori_vt8_1/*.e-s00?[2].* | pigz -9 -p 20 > ex_gbAnisotropyMisori_vt.tgz