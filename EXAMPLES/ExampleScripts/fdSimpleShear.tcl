# Material model: FiniteDeformationElastic3D
# Element model: TotalLagrangianFD20NodeBrick
# Simple Shear Test
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
fix  7 1 1 1
fix  6 1 1 1
fix 14 1 1 1
fix  5 1 1 1
fix  8 1 1 1
fix 15 1 1 1
fix 13 1 1 1
fix 16 1 1 1

fix  1 0 1 1
fix  2 0 1 1
fix  3 0 1 1
fix  4 0 1 1
fix  9 0 1 1
fix 10 0 1 1
fix 11 0 1 1
fix 12 0 1 1
fix 17 0 1 1
fix 18 0 1 1
fix 19 0 1 1
fix 20 0 1 1

# specific large deformation material model

# Neo-Hookean
#set WE "-NH 11.83e5 0.4"

## Logarithmic
set WE "-Log 11.83e5 0.4"

## Ogden
#set WE "-Ogden 11.83e5 0.4 3 6.3e5 0.012e5 -0.1e5 1.3 5.0 -2.0"

#Mooney-Rivlin
#set WE "-MR 11.83e5 0.4 1.8484375e5 0.2640625e5"

# finite deformation elastic 3D material model
nDMaterial FiniteDeformationElastic3D 1 -WEnergy $WE 0.0

# total lagrangian finite deformation node brick
element TLFD20nBrick 1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 1 0.0 0.0 0.0

set load1 -2.0e5
set load2 8.0e5

pattern Plain 1 "Linear"      {
           load   1  $load1 0 0
           load   2  $load1 0 0
           load   3  $load1 0 0
           load   4  $load1 0 0
           load   9  $load2 0 0
           load  10  $load2 0 0
           load  11  $load2 0 0
           load  12  $load2 0 0  }

## epualize displacements
## make it to be simple shear
equalDOF  1   2  1  1
equalDOF  1   3  1  1
equalDOF  1   4  1  1
equalDOF  9  10  1  1
equalDOF  9  11  1  1
equalDOF  9  12  1  1
equalDOF 17  18  1  1
equalDOF 17  19  1  1
equalDOF 17  20  1  1
#


set ndz 0.02

#integrator LoadControl 0.01
integrator DisplacementControl 1 1 $ndz 2 $ndz $ndz
test NormDispIncr 1.0e-6 50  1
system SparseGeneral
constraints Plain
algorithm Newton
numberer Plain
analysis Static

set NN2 100

#recorder Element 1 -file SSelement.out stress
#recorder Node  -file SSnode1.out  -node 1 -node 2 -dof 1 disp
#recorder Element 1 -file SSforce.out force

analyze $NN2

print node 1 2   17 18

wipe
