wipe

file mkdir data

model basic -ndm 2 -ndf 3

node 1 0. 0.
node 2 3. 0.


fix 1 1 1 1



geomTransf Linear 1


set E 2.e11


uniaxialMaterial Steel01 1 2.354e8 $E 0.02

section Fiber 1 {

patch quad 1 2 8 -0.15 0.125 -0.15 -0.125 0.15 -0.125 0.15 0.125

}

element nonlinearBeamColumn 1 1 2 5 1 1



recorder Node -file data/node2disp.out -time -node 2 -dof 1 2 3 disp
recorder Node -file data/node1reac.out -time -node 1 -dof 1 2 3 reaction


puts "Model Done!"