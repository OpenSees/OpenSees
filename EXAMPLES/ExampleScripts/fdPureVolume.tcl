# Material model: FiniteDeformationElastic3D
# Element model: TotalLagrangianFD20NodeBrick
# Pure Volume Change Test
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
fix  7   1 1 1
fix  8   0 1 1
fix 15   0 1 1
fix  6   1 0 1
fix 14   1 0 1
fix  5   0 0 1
fix 13   0 0 1
fix 16   0 0 1

fix  3   1 1 0
fix  4   0 1 0
fix 11   0 1 0
fix  2   1 0 0
fix 10   1 0 0

fix 19   1 1 0
fix 20   0 1 0
fix 18   1 0 0

# specific large deformation material model
# finite deformation elastic 3D material model
nDMaterial FiniteDeformationElastic3D 1 NeoHookean3D 1971.67 422.5 0.0
#nDMaterial FiniteDeformationElastic3D 1 DecoupledNH3D 1971.67 422.5 0.0
#nDMaterial FiniteDeformationElastic3D 1 DecoupledLog3D 1971.67 422.5 0.0

# total lagrangian finite deformation node brick
element TLFD20nBrick 1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 1 0.0 0.0 0.0

set load 2.0e6
set load1p [expr $load]
set load2p [expr -4.0 * $load]

pattern Plain 1 "Linear"      {
	load   1  $load1p  $load1p  $load1p
	load   2        0  $load1p  $load1p
	load   3        0        0  $load1p
	load   4  $load1p        0  $load1p
	load   9        0  $load2p  $load2p
	load  10        0        0  $load2p
	load  11        0        0  $load2p
	load  12  $load2p        0  $load2p

	load   5  $load1p  $load1p        0
	load   6        0  $load1p        0
	load  13        0  $load2p        0
	load   8  $load1p        0        0
	load  16  $load2p        0        0

	load  17  $load2p  $load2p	  0
	load  18        0  $load2p	  0
	load  20  $load2p        0	  0 }

set ndz 0.01

integrator DisplacementControl 1 1 $ndz 2 $ndz $ndz
test NormDispIncr 1.0e-5 20  1
system SparseGeneral
constraints Plain
numberer Plain
analysis Static

set NN2 20

recorder Element 1 -file PVelement.out stress
recorder Element 1 -file PVforce.out force

analyze $NN2

print node 1 2  19 20

wipe
