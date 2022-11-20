model basic -ndm 2 -ndf 3

# Define materials 1) confined, 2) unconfined, 3 steel
uniaxialMaterial Concrete01  1  -6.0  -0.004   -5.0     -0.014
uniaxialMaterial Concrete01  2  -5.0   -0.002   0.0     -0.006
set fy 60.0;      # Yield stress
set E 30000.0;    # Young's modulus
uniaxialMaterial Steel01  3  $fy $E 0.01

# set some parameters
set width 15
set depth 24 
set cover  1.5
set As    0.60;     # area of no. 7 bars

source RCsection2D.tcl
RCsection2D 1 $depth $width $cover 1 2 3 3 1 $As 10 10 2

# Estimate yield curvature
# (Assuming no axial load and only top and bottom steel)
set d [expr $depth-$cover]	;# d -- from cover to rebar
set epsy [expr $fy/$E]	;# steel yield strain
set Ky [expr $epsy/(0.7*$d)]

# Print estimate to standard output
puts "Estimated yield curvature: $Ky"

# Set axial load 
set P -180

set mu 15;		# Target ductility for analysis
set numIncr 100;	# Number of analysis increments

# Call the section analysis procedure
source MomentCurvature.tcl
MomentCurvature 1 $P [expr $Ky*$mu] $numIncr