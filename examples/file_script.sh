# tar -cvf - case1_gbAnisotropy_v1/*.e-s??0?.* | pigz -9 -p 20 > case1_gbAnisotropy_v1_ex.tgz

tar -cvf - case1_gbAnisotropy/*.e.* | pigz -9 -p 20 > case1_gbAnisotropy_ex.tgz

# ACSEDGPolyEBSD PolycrystalStoredEnergyEBSDAction DeformedGrainEBSDMaterial

in_dir_path="/home/pw-moose/projects/moose/modules/phase_field"
out_dir_path="/home/pw-moose/projects/yinglong"

# mkdir ${out_dir_path}/src/kernels/
# mkdir ${out_dir_path}/include/kernels/

/home/pw-moose/projects/moose/modules/phase_field/src/materials/DeformedGrainMaterial.C

cp ${in_dir_path}/src/materials/DeformedGrainMaterial.C ${out_dir_path}/src/materials/DeformedGrainEBSDMaterial.C
cp ${in_dir_path}/include/materials/DeformedGrainMaterial.h ${out_dir_path}/include/materials/DeformedGrainEBSDMaterial.h

code ${out_dir_path}/src/materials/DeformedGrainEBSDMaterial.C
code ${out_dir_path}/include/materials/DeformedGrainEBSDMaterial.h