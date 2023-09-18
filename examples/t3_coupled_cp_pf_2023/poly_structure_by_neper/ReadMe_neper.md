> 用于记录使用neper创建多晶微观结构的学习笔记以及命令行

# 学习笔记


# 命令行

```bash

neper -T -n 100 -domain "cube(100,100,100)"
neper -T -n 20 -morpho graingrowth -domain "cube(20,20,20)"
neper -T -n 200 -dim 2 -morpho graingrowth -oridescriptor euler-bunge -domain "square(100,100)"
#  -periodicity all

neper -M n20-id1.tess -elttype hex
neper -M n200-id1.tess -elttype tri -rcl "body>5?0.35:1"
neper -M n20-id1.tess -elttype tri -rcl 1.0
# tri quad

neper -V n200-id1.msh -datacellcol id -print n200_msh

cp n200-id1.msh n200-id2.msh
```