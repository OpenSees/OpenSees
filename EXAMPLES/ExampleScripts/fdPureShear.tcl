# Material model: FiniteDeformationElastic3D
# Element model: TotalLagrangianFD20NodeBrick
# Pure Shear Test
# Zhao Cheng and Boris Jeremic [zcheng,jeremic]@ucdavis.edu

# wipe previous analysis
wipeAnalysis

set S2 1.41421356237309504880168872

# create the modelbuilder
model BasicBuilder -ndm 3 -ndf 3

# build the model
node  1   [expr  $S2/2]              0.0          1.0
node  2             0.0    [expr  $S2/2]          1.0
node  3   [expr -$S2/2]              0.0          1.0
node  4             0.0    [expr -$S2/2]          1.0
node  5   [expr  $S2/2]              0.0          0.0
node  6             0.0    [expr  $S2/2]          0.0
node  7   [expr -$S2/2]              0.0          0.0
node  8             0.0    [expr -$S2/2]          0.0
node  9   [expr  $S2/4]	   [expr  $S2/4]	  1.0
node 10   [expr -$S2/4]	   [expr  $S2/4]	  1.0
node 11   [expr -$S2/4]	   [expr -$S2/4]	  1.0
node 12   [expr  $S2/4]	   [expr -$S2/4]	  1.0
node 13   [expr  $S2/4]	   [expr  $S2/4]	  0.0
node 14   [expr -$S2/4]	   [expr  $S2/4]	  0.0
node 15   [expr -$S2/4]	   [expr -$S2/4]	  0.0
node 16   [expr  $S2/4]	   [expr -$S2/4]	  0.0
node 17   [expr  $S2/2]              0.0          0.5
node 18             0.0    [expr  $S2/2]          0.5
node 19   [expr -$S2/2]              0.0          0.5
node 20             0.0    [expr -$S2/2]          0.5


# boundary condition
fix  7  1  1  1
fix  8  0  0  1
fix  6  0  0  1
fix  5  0  1  1

fix  3  1  1  1
fix  4  0  0  1
fix  2  0  0  1
fix  1  0  1  1

fix 19  1  1  1
fix 20  0  0  1
fix 18  0  0  1
fix 17  0  1  1

fix  9  0  0  1
fix 10  0  0  1
fix 11  0  0  1
fix 12  0  0  1
fix 13  0  0  1
fix 14  0  0  1
fix 15  0  0  1
fix 16  0  0  1

# specific large deformation material model

# Neo-Hookean
set WE "-NH 11.83e5 0.4"

## Logarithmic
#set WE "-Log 11.83e5 0.4"

## Ogden
#set WE "-Ogden 11.83e5 0.4 3 6.3e5 0.012e5 -0.1e5 1.3 5.0 -2.0"

#Mooney-Rivlin
#set WE "-MR 11.83e5 0.4 1.8484375e5 0.2640625e5"


# finite deformation elastic 3D material model
nDMaterial FiniteDeformationElastic3D 1 -WEnergy $WE 0.0

# total lagrangian finite deformation node brick
element TLFD20nBrick 1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 1 0.0 0.0 0.0

set ld0 1.0e5
set ld1p  [expr     $ld0*$S2/3]
set ld1n  [expr         -$ld1p]
set ld2p  [expr   4*$ld0*$S2/3]
set ld2n  [expr         -$ld2p]
set ld3p  [expr        $ld2p/2]
set ld3n  [expr        $ld2n/2]

pattern Plain 1 "Linear" {
        load   1  $ld1n      0      0
	load   5  $ld1n      0      0
	load  17  $ld2p      0      0

	load   2      0  $ld1p 	    0
	load   6      0  $ld1p      0
	load  18      0  $ld2n      0

	load   4      0  $ld1n      0
	load   8      0  $ld1n      0
	load  20      0  $ld2p      0

	load   9  $ld3p  $ld3n      0  
	load  13  $ld3p  $ld3n      0 	 
	load  12  $ld3p  $ld3p      0 	
	load  16  $ld3p  $ld3p      0

	load  10  $ld3n  $ld3n      0  
	load  14  $ld3n  $ld3n      0 	 
	load  11  $ld3n  $ld3p      0 	
	load  15  $ld3n  $ld3p      0 }
	 
set ndz 0.01

#integrator LoadControl 0.02
integrator DisplacementControl 1 1 $ndz 2 $ndz $ndz
test NormDispIncr 1.0e-5 50  1
system SparseGeneral
constraints Plain
algorithm Newton
numberer Plain
analysis Static

set NN2 100

#recorder Element 1 -file element.out stress
#recorder Element 1 -file force.out force
#recorder Node  -file node1.out  -node 1 -node 2 -dof 1 disp

analyze $NN2

print node 1 2 3 4  17 18 19 20


wipe
