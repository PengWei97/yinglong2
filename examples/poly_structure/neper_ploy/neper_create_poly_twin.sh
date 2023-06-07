#!/bin/bash

# 清除旧的文件
rm -rf d.*

# 循环
for i in {3..3}; do
    # 根据循环索引设置参数值
    id_value="$i"
    img_value="img16_$i"

    # 运行neper命令
    neper -T -dim 2 -n "10::msfile(n_msfile)" -morpho "gg::msfile(morpho_msfile)" -domain "square(100,100)" -o d -ori "uniform" -oridescriptor 'euler-bunge:passive' -id "16"

    neper -V d.tess -datacellcol "scaleid(1)" -datacelltrs 0.1 -print "$img_value"

    neper -M d.tess -statelset id,nodes -statnode id,x,y,z -elttype "quad" -cl 1.0
done


# neper -M d.tess -statelset id,nodes -statnode id,x,y,z -elttype "quad" -cl 0.1 -id 1

# lamellar(w=0.1:0.02:0.08:0.02:1.0,v=(1,-0.5,0),pos=start)
# neper -T -n 1 -ori 'fiber(1,1,1,0,0,1)' -oridescriptor 'euler-bunge:passive' -format 'ori'
# -ori "fiber(0,0,1,0,0,1)"
#  -ori "random"
# lamellar(w=34.0:5.0:50.0,v=(0,1,0),pos=start)
# 3 - lamellar(w=20.0:5.0:7.5:5.0:7.5:5.0:20.0,v=(0,1,0),pos=start)
# 4 - lamellar(w=10.0:5.0:7.5:5.0:7.5:5.0:7.5:5.0:20.0,v=(0,1,0),pos=start)
# 5 - lamellar(w=5.0:5.0:6.0:5.0:6.0:5.0:6.0:5.0:6.0:5.0:20.0,v=(0,1,0),pos=start)