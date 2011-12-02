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

#equalDOF	1	3	2	3
#equalDOF	1	3	3	3
#equalDOF	1	3	4	3
#equalDOF	1	3	9	3
#equalDOF	1	3	10	3
#equalDOF	1	3	11	3
#equalDOF	1	3	12	3

# specific large deformation material model
nDMaterial FiniteDeformationElastic3D 1 NeoHookean3D 1971.67 422.5 0.0

## Von Mises *****************************************
## finite deformation yield function
set fdy "-VM 60"
## finite deformation flow rule
set fdf "-VM 60"
## finite deformation isotropic (scalar) evolution law
set fdes "-LS 500.0"
## finite deformation kinematic (tensor) evolution law
set fdet "-Linear 100.0"

## Druke-Prager **************************************
# fd yield surface
#set fdy "-DP 0 30 "
# fd flow rule
#set fdf "-DP 0 "
# fd isotropic (scalar) evolution law
#set fdes "-LS 100.0"

# finite deformation material
nDMaterial FiniteDeformationEP3D 2  1 -fdY $fdy -fdF $fdf -fdES $fdes -fdET $fdet 
#nDMaterial FiniteDeformationEP3D 2  1 -fdY $fdy -fdF $fdf -fdES $fdes 
#nDMaterial FiniteDeformationEP3D 2  1 -fdY $fdy -fdF $fdf -fdET $fdet 
#nDMaterial FiniteDeformationEP3D 2  1 -fdY $fdy -fdF $fdf 

# total lagrangian finite deformation node brick
element TLFD20nBrick 1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 2 0.0 0.0 0.0

set load1 -200.0
set load2 800.0

pattern Plain 1 "Linear"     {
           load   1   0 0 $load1
           load   2   0 0 $load1
           load   3   0 0 $load1
           load   4   0 0 $load1
           load   9   0 0 $load2
           load  10  0 0 $load2
           load  11  0 0 $load2
           load  12  0 0 $load2  }

set ndz 0.01

set NN2 20
	
#integrator LoadControl 0.01
integrator DisplacementControl 1 3 $ndz 2 $ndz $ndz
test NormDispIncr 1.0e-7 50  1
system SparseGeneral
constraints Plain
algorithm KrylovNewton
numberer Plain
analysis Static
	
recorder Element 1 -file PUstress.out stress
recorder Element 1 -file PUforce.out force
recorder Node  -file PUnode1.out  -time -node 1  -dof 3

for {set i 1} {$i <=$NN2} {incr i} {
    puts $i
    analyze 1    
}

print node 1 2   17 18

wipe
