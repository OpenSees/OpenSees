# Material model: FiniteDeformationElastic3D
# Element model: TotalLagrangianFD20NodeBrick
# Pure Unixial Test
# Zhao Cheng and Boris Jeremic [zcheng,jeremic]@ucdavis.edu

# wipe previous analysis
wipeAnalysis

# create the modelbuilder
model BasicBuilder -ndm 3 -ndf 3
# build the model
node  1 +1.0 +1.0 +1.0
node  2 -1.0 +1.0 +1.0
node  3 -1.0 -1.0 +1.0
node  4 +1.0 -1.0 +1.0
node  5 +1.0 +1.0 -1.0
node  6 -1.0 +1.0 -1.0
node  7 -1.0 -1.0 -1.0
node  8 +1.0 -1.0 -1.0
node  9  0.0 +1.0 +1.0
node 10 -1.0  0.0 +1.0
node 11  0.0 -1.0 +1.0
node 12 +1.0  0.0 +1.0
node 13  0.0 +1.0 -1.0
node 14 -1.0  0.0 -1.0
node 15  0.0 -1.0 -1.0
node 16 +1.0  0.0 -1.0
node 17 +1.0 +1.0  0.0
node 18 -1.0 +1.0  0.0
node 19 -1.0 -1.0  0.0
node 20 +1.0 -1.0  0.0

# boundary condition
fix   7  1 1 1
fix   6  1 0 1
fix   5  0 0 1
fix   8  0 1 1
fix  15  0 1 1
fix  14  1 0 1
fix  13  0 0 1
fix  16  0 0 1

# specific large deformation material model
nDMaterial FiniteDeformationElastic3D 1 NeoHookean3D 1971.67 422.5 0.0
#nDMaterial FiniteDeformationElastic3D 1 DecoupledNH3D 1971.67 422.5 0.0
#nDMaterial FiniteDeformationElastic3D 1 DecoupledLog3D 1971.67 422.5 0.0
#nDMaterial FiniteDeformationElastic3D 1 DecoupledOS3D 3 6.3e5 0.012e5 -0.1e5 1.3 5.0 -2.0 19.7167e5 0.0
#nDMaterial FiniteDeformationElastic3D 1 DecoupledMRS3D 184.84375 26.40625 1971.67 0.0


# total lagrangian finite deformation node brick
element TLFD20nBrick 1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 1 0.0 0.0 0.0

set load1 -200
set load2 800

pattern Plain 1 "Linear"     {
           load   1   0 0 $load1
           load   2   0 0 $load1
           load   3   0 0 $load1
           load   4   0 0 $load1
           load   9   0 0 $load2
           load  10  0 0 $load2
           load  11  0 0 $load2
           load  12  0 0 $load2  }

#print node 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
#print info.out
#recorder Element 1 -file element.out stress

set ndz 0.02

integrator LoadControl 0.01
#integrator DisplacementControl 1 3 $ndz 2 $ndz $ndz
test NormDispIncr 1e-6 30  1
system SparseGeneral
constraints Plain
algorithm KrylovNewton
numberer Plain
analysis Static

set NN2 20

analyze $NN2

recorder Element 1 -file -time element.out stress
recorder Node  -file node1.out  -time -node 1 -node 2 -dof 1 2 3 disp

print node 1 2   17 18

wipe
