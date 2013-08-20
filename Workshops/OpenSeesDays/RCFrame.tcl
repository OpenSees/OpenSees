model Basic -ndm 2 -ndf 3

set ft 12.0
set in 1.0
set cm 0.3937

# Set parameters for overall model geometry
set width    [expr 42.0*$ft]
set height   [expr 36.0*$ft]

# Create nodes
node  1       0.0     0.0 
node  2    $width     0.0 
node  3       0.0 $height
node  4    $width $height

# set boundary conditions
fix   1     1    1    1
fix   2     1    1    1

# create materials for sections
uniaxialMaterial Concrete01  1  -6.0 -0.004 -5.0 -0.014; # core
uniaxialMaterial Concrete01  2  -5.0 -0.002  0.0 -0.006; # cover
uniaxialMaterial Steel01  3 60.0 30000.0 0.01

# create sections
set bWidth  [expr 5.0*$ft]
set bDepth  [expr 5.0*$ft]
set cover  1.5
set As     0.60;     # area of no. 7 bars

source RCsection2D.tcl
RCsection2D 1 $bWidth $bDepth $cover 1 2 3 6 3 \
    $As 10 10 2

# create geometric transformations
geomTransf PDelta 1  
#geomTransf Linear 1  
geomTransf Linear 2  

# Create the columns using distributed plasticisty
set np 5
set eleType forceBeamColumn; # forceBeamColumn od dispBeamColumn
element $eleType  1   1   3   $np    1       1 
element $eleType  2   2   4   $np    1       1 

# Create the beam element using elastic element
set d [expr 8.0*$ft];		# Beam Depth
set b [expr 5.0*$ft];		# Beam Width
set Eb 3600.0; set Ab [expr $b*$d];; set Izb [expr ($b*pow($d,3))/12.0]; 
section Elastic 2   $Eb $Ab $Izb;	# elastic beam section
#element elasticBeamColumn   3   3   4   $Ab   $Eb  $Izb    2
element $eleType  3   3   4   $np    2       2 
