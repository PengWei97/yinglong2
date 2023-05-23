neper -T -dim 2 -n "10::msfile(n_msfile)" -morpho "gg::msfile(morpho_msfile)" -domain "square(10,10)" -o d -ori "uniform" -oridescriptor 'euler-bunge:passive'
neper -V d.tess -datacellcol "scaleid(1)" -datacelltrs 0.1 -print img1

neper -M d.tess -statelset id,nodes -statnode id,x,y,z -elttype "quad" -cl 0.1

# lamellar(w=0.1:0.02:0.08:0.02:1.0,v=(1,-0.5,0),pos=start)
# neper -T -n 1 -ori 'fiber(1,1,1,0,0,1)' -oridescriptor 'euler-bunge:passive' -format 'ori'
# -ori "fiber(0,0,1,0,0,1)"
#  -ori "random"