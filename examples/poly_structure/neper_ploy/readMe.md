# beginning
1. to build neper operating environment
   1. [搬运】NEPER安装教程（包括Ubuntu虚拟机，前置环境的安装）](https://www.bilibili.com/video/BV1cq4y1G7vt/?spm_id_from=333.337.search-card.all.click&vd_source=0c01f528c72373cb8a1352ad0923483f)
   2. 
2. Reference URL 1：[Generating twin and multiple twinning in a grain #697](https://github.com/neperfepx/neper/discussions/697)
3. Reference URL 2：[Converting Neper Mesh Format to Dream3D Inl Format #698](https://github.com/neperfepx/neper/discussions/698)

# Operating procedures
1. First, use the following command to create neper microstructure in neper:
```bash
neper -T -dim 2 -n "10::msfile(n_msfile)" -morpho "gg::msfile(morpho_msfile)" -domain "square(10,10)" -o d -ori "uniform" -oridescriptor 'euler-bunge:passive'
neper -V d.tess -datacellcol "scaleid(1)" -datacelltrs 0.1 -print img1

neper -M d.tess -statelset id,nodes -statnode id,x,y,z -elttype "quad" -cl 0.1
```
ps: to create the necessary files, including .msh, .stelset, and .stdnode.

2. Second, run `msh_to_inl.m` to getting the `case_msh2inl.inl` file.


