wipe
file mkdir opensees_results
model basic -ndm 2 -ndf 3
set e 2.e11
set i 6.572e-5
set a 4.265e-3
node 1 0. 0.
node 2 0. 3.
node 3 0. 3.
node 4 4. 3.
node 5 4. 0.
node 6 2. 3.

fix 1 1 1 1
fix 5 1 1 1

equalDOF 2 3 1 2

set transftag 1
geomTransf Linear $transftag
element elasticBeamColumn 1 1 2 $a $e $i $transftag



element elasticBeamColumn 21 3 6 $a $e $i $transftag
element elasticBeamColumn 22 6 4 $a $e $i $transftag

element elasticBeamColumn 3 4 5 $a $e $i $transftag

#recorder Node -file data/node1.out -time -node 1 -dof 1 2 reaction
#recorder Node -file data/node5.out -time -node 5 -dof 1 2 reaction
recorder Node -file opensees_results/opensees_disp_6.txt -time -node 6 -dof 1 2 3 disp
#recorder Node -file data/node2disp.out -time -node 2 -dof 1 2 3 disp

puts "model build!"