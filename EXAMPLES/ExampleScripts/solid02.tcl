#jeremic@ucdavis.edu zhyang@ucdavis.edu
#Example with 8 brick elements and one pile
# 20Jan2001


# tcl version of main_T3Dep_1_dy.cpp


# ################################
# create the modelbuilder
# #################################

model BasicBuilder -ndm 3 -ndf 3

set g 9.81
set rho_s 1.8

# Elastic modulus of stiff soil kPa
set Es 100000

# ################################
# build the model
# #################################

node   1  -2.000000E+00 -2.000000E+00 -2.000000E+00
node   2   0.000000E+00 -2.000000E+00 -2.000000E+00
node   3   2.000000E+00 -2.000000E+00 -2.000000E+00
node   4  -2.000000E+00  0.000000E+00 -2.000000E+00
node   5   0.000000E+00  0.000000E+00 -2.000000E+00
node   6   2.000000E+00  0.000000E+00 -2.000000E+00
node   7  -2.000000E+00  2.000000E+00 -2.000000E+00
node   8   0.000000E+00  2.000000E+00 -2.000000E+00
node   9   2.000000E+00  2.000000E+00 -2.000000E+00
node  10  -2.000000E+00 -2.000000E+00  0.000000E+00
node  11   0.000000E+00 -2.000000E+00  0.000000E+00
node  12   2.000000E+00 -2.000000E+00  0.000000E+00
node  13  -2.000000E+00  0.000000E+00  0.000000E+00
node  14   0.000000E+00  0.000000E+00  0.000000E+00
node  15   2.000000E+00  0.000000E+00  0.000000E+00
node  16  -2.000000E+00  2.000000E+00  0.000000E+00
node  17   0.000000E+00  2.000000E+00  0.000000E+00
node  18   2.000000E+00  2.000000E+00  0.000000E+00
node  19  -2.000000E+00 -2.000000E+00  2.000000E+00
node  20   0.000000E+00 -2.000000E+00  2.000000E+00
node  21   2.000000E+00 -2.000000E+00  2.000000E+00
node  22  -2.000000E+00  0.000000E+00  2.000000E+00
node  23   0.000000E+00  0.000000E+00  2.000000E+00
node  24   2.000000E+00  0.000000E+00  2.000000E+00
node  25  -2.000000E+00  2.000000E+00  2.000000E+00
node  26   0.000000E+00  2.000000E+00  2.000000E+00
node  27   2.000000E+00  2.000000E+00  2.000000E+00

fix  1 1 1 1
fix  2 1 1 1  
fix  3 1 1 1
fix  4 1 1 1
fix  5 1 1 1
fix  6 1 1 1
fix  7 1 1 1
fix  8 1 1 1
fix  9 1 1 1
fix 10 1 0 1
fix 11 1 0 1
fix 12 1 0 1
fix 13 1 0 1
fix 15 1 0 1
fix 16 1 0 1
fix 17 1 0 1
fix 18 1 0 1
fix 19 1 0 1
fix 20 1 0 1
fix 21 1 0 1
fix 22 1 0 1
fix 24 1 0 1
fix 25 1 0 1
fix 26 1 0 1
fix 27 1 0 1
    
#add beam nodes
model BasicBuilder -ndm 3 -ndf 6

node  28   0.000000E+00  0.000000E+00  0.000000E+00
node  29   0.000000E+00  0.000000E+00  2.000000E+00
node  30   0.000000E+00  0.000000E+00  10.000000E+00

#add the mass of superstructure
mass 30 0.00334 0.00334 0.00334	 0.00334 0.00334 0.00334

equalDOF 10 16 2
equalDOF 11 17 2
equalDOF 12 18 2
equalDOF 19 25 2
equalDOF 20 26 2
equalDOF 21 27 2

equalDOF 14 28 1 2 3
equalDOF 23 29 1 2 3

# elastic material
nDMaterial ElasticIsotropic3D 1 $Es 0.3 $rho_s
# the template material of yours
#sset YS {DruckerPrager }
#set PS {DruckerPrager 0.05}
#
#set startstress {50 0.001 0.0}
#set otherstress {0 0 0}
#set scalars {0.05 0 0.85}
#set tensors {0.0 0.0 0.0}
#set NOS 3
#set NOT 3
#set EPS {3000.0 3000.0 0.3 1.8 $startstress $otherstrain $otherstrain 
#         $otherstrain $otherstrain $NOS $scalars $NOT $tensors}
#
#nDMaterial Template 1 -YS $YS -PS $PS -EPS $EPS ......

#            tag            8 nodes                  matID    bforce1,2&3    	 massdensity
#element brick  1   5    6    7    8    1    2    3   4 1  0.0 0.0 -9.81     
element Brick8N  1    2    5    4    1   11   14   13  10 1  0.0 0.0 [expr -$g]	 1.8
element Brick8N  2    3    6    5    2   12   15   14  11 1  0.0 0.0 [expr -$g]	 1.8
element Brick8N  3    5    8    7    4   14   17   16  13 1  0.0 0.0 [expr -$g]	 1.8
element Brick8N  4    6    9    8    5   15   18   17  14 1  0.0 0.0 [expr -$g]	 1.8
element Brick8N  5   11   14   13   10   20   23   22  19 1  0.0 0.0 [expr -$g]	 1.8
element Brick8N  6   12   15   14   11   21   24   23  20 1  0.0 0.0 [expr -$g]	 1.8
element Brick8N  7   14   17   16   13   23   26   25  22 1  0.0 0.0 [expr -$g]	 1.8
element Brick8N  8   15   18   17   14   24   27   26  23 1  0.0 0.0 [expr -$g]	 1.8

# Elastic modulus of aluminium kPa
set Ea 7.0e7

# Shear modulus of aluminium kPa
set G 2.63e7

# Beam cross-sectional area
set Abeam 1.493e-4 

# Beam polar moment of inertia
set Jbeam 1.463e-8

# Beam moment of inertia
set Ibeam 0.732e-8

# Geometric transformation
#                tag  vecxz
geomTransf Linear 1   0 1 0

# Define beam elements
#                         tag ndI ndJ  A    E  G     Jx    Iy    Iz   transf
element elasticBeamColumn  9  28  29 $Abeam $Ea $G $Jbeam $Ibeam $Ibeam  1
element elasticBeamColumn 10  29  30 $Abeam $Ea $G $Jbeam $Ibeam $Ibeam  1


set Series "Path -filePath NR228.txt -dt 0.005 -factor 1"

pattern UniformExcitation  1   2  -accel $Series 

# create the recorder
#recorder Node Node.out disp -time -node  1 2 3 4 5 6 7 8 -dof 1 2 3
recorder Node -file node1.out -time -node 14 23 30 -dof 2 disp

recorder plot node.out Middel_of_soil 10 10 300 300 -columns 2 1
recorder plot node1.out "PEER workshop, solid02.tcl: Top_of_soil" 0 0 500 100 -columns 1 3 
recorder plot node1.out "PEER workshop, solid02.tcl: Superstructure" 0 130 500 100 -columns 1 4 

# #################################
# create the transient analysis
# #################################

integrator Newmark  0.55  0.2756
numberer RCM
#constraints Plain
constraints Penalty 1e12 1e12
#constraints Transformation    
test NormDispIncr 2.0e-5 20 0

#constraints Lagrange 1.0 1.0
#test NormDispIncr 1.0e-10 10 1

algorithm Newton
numberer RCM
system UmfPack

analysis Transient

# ################################
# perform the analysis
# #################################

analyze 375 0.04


wipe
