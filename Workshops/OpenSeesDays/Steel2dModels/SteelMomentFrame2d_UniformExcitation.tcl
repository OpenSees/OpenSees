# written: fmk
# units: kip & in


# following are required variables that need to be set prior to running script:
# floorOffsets : list of floor offsets
# colOffsets   : list of col offstes
# colSizes     : list of col W Sections for interior columns
# colExtSizes  : list of col W sections for exterior columns
# beamSizes    : list of beam W sections for each floor. 
# floorWeight  : total weight of typ floor for mass & beam load calculation
# roofWeight   : total weight of typ floor for mass & beam load calculation
# percentLoadFrame: percentage of total weight taken by frame for loading calc
# Fy : steel yield strength
# E  : steel elastic modulus
# b  : steel hardening ratio
# gMotion  : name of earthquake acceleration record, PEER AT2 format required
# scale   : scale factor to be applied to earthquake gMotion

wipe;

# set some constants
set PI 3.14159
set in 1.0;
set g 386.4;	

# load procedures in other files
source Steel2d.tcl
source ReadRecord.tcl;

# determine some properties based on inputs
set numFloor [expr [llength $floorOffsets]+1]
set numCline [expr [llength $colOffsets]+1]
set width 0.
for {set i 0; set width 0;} {$i < [expr $numCline-1]} {incr i 1} {
    set width [expr $width + [lindex $colOffsets $i]]
}
set massAtFloorNode [expr $floorWeight*$percentMassFrame/($g*$numCline*1.0)]
set massAtRoofNode  [expr $roofWeight*$percentMassFrame/($g*$numCline*1.0)]
set uniformRoofLoad  [expr $roofWeight*$percentLoadFrame/$width]
set uniformFloorLoad [expr $floorWeight*$percentLoadFrame/$width]

# check list dimensions for errors
if {[llength $colSizes] != [expr $numFloor-1]} {puts "ERROR: colSizes"; quit}
if {$numCline >= 10} {puts "ERROR: too many column lines, reprogram"; quit}
if {$numFloor >= 10} {puts "ERROR: too many floors, reprogram"; quit}

#
# Model Generation
#

# Define the model builder, ndm = #dimension, ndf = #dofs at a node
model BasicBuilder -ndm 2 -ndf 3;  

# Build the Nodes
for {set floor 1; set floorLoc 0} {$floor <= $numFloor} {incr floor 1} {
    if {$floor == $numFloor} {
	set massX $massAtRoofNode
	set massY [expr $massX*$percentLoadFrame]; # gravity cols take vertical
    } else {
	set massX $massAtFloorNode
	set massY [expr $massX*$percentLoadFrame]; 
    }
    for {set colLine 1; set colLoc 0;} {$colLine <= $numCline} {incr colLine 1} {
	node $colLine$floor $colLoc $floorLoc -mass $massX $massY 0.
	if {$floor == 1} {
	    fix $colLine$floor 1 1 1
	}
	if {$colLine < $numCline} {
	    set colLoc [expr $colLoc + [lindex $colOffsets [expr $colLine-1]]]
	}
    }
    if {$floor < $numFloor} {
	set floorLoc [expr $floorLoc + [lindex $floorOffsets [expr $floor-1]]]
    }
}
	    
# define material 
uniaxialMaterial Steel02 1 $Fy $E $b 20 0.925 0.15

# build the columns
geomTransf PDelta 1
for {set colLine 1} {$colLine <= $numCline} {incr colLine 1} {
    for {set floor1 1; set floor2 2} {$floor1 < $numFloor} {incr floor1 1; incr floor2 1} {
	if {$colLine == 1 || $colLine == $numCline} {
	    set theSection [lindex $colExtSizes [expr $floor1 -1]]
	} else {
	    set theSection [lindex $colSizes [expr $floor1 -1]]
	}
	ForceBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection 1 1 -nip 5
	#	ElasticBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection 29000 1 
    }
}

# build the beams
geomTransf Linear 2
for {set colLine1  1; set colLine2 2} {$colLine1 < $numCline} {incr colLine1 1; incr colLine2 1} {
    for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
	set theSection [lindex $beamSizes [expr $floor -2]]
	ForceBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor $colLine2$floor $theSection 1 2
	#	ElasticBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor $colLine2$floor $theSection 29000 2
    }
}

#
# add the gravity loads
#

pattern Plain 101 Linear {
    for {set colLine1  1} {$colLine1 < $numCline} {incr colLine1 1} {
	set colLine2 [expr $colLine1 + 1]
	for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
	    if {$floor == $numFloor} {
		eleLoad -ele $colLine1$floor$colLine2$floor -type beamUniform -$uniformRoofLoad
	    } else {
		eleLoad -ele $colLine1$floor$colLine2$floor -type beamUniform -$uniformFloorLoad
	    }
	}
    }
}

#
#  Perform Gravity-analysis: load-controlled static analysis
#
set Tol 1.0e-6;    
constraints Plain;	
numberer RCM;		
system BandGeneral;	
test NormDispIncr $Tol 10;
algorithm Newton;		
integrator LoadControl 0.1;
analysis Static;		
analyze 10

#
# set gravity loads constant and reset time to zero
#

loadConst -time 0.0

# add some damping 
#  NOTE damping mass and initial stiffness as opposed to currnet stiffness

set lambda [eigen $mode2]
set omegaI [expr pow([lindex $lambda [expr $mode1-1]],0.5)];
set omegaJ [expr pow([lindex $lambda [expr $mode2-1]],0.5)];
set alphaM [expr $dampRatio*(2*$omegaI*$omegaJ)/($omegaI+$omegaJ)];	
set betaKcomm [expr 2.*$dampRatio/($omegaI+$omegaJ)]; 
rayleigh $alphaM 0. 0. $betaKcomm 

# puts " Period T1: [expr 2*$PI/sqrt([lindex $lambda 0])] sec"

# create a load patrern for uniform excitation
ReadRecord $gMotion.AT2 $gMotion.g3 dt nPt;
timeSeries Path 10 -filePath $gMotion.g3 -dt $dt -factor [expr $g*$scale]
pattern UniformExcitation 1 1 -accel 10;

#
# create some recorders to monitor nodal disp and accelerations
# for nodes column line 1
#

set nodeList []
for {set floor 1} {$floor <= $numFloor} {incr floor 1} {
    lappend nodeList 1$floor
}

set cmd "recorder EnvelopeNode -file floorAccEnv.out  -timeSeries 10 -node $nodeList -dof 1 accel"; eval $cmd;
set cmd "recorder EnvelopeNode -file floorAccEnv.out  -timeSeries 10 -node $nodeList -dof 1 accel"; eval $cmd;

set cmd "recorder EnvelopeNode -file floorDispEnv.out  -node $nodeList -dof 1 disp"; eval $cmd;
set cmd "recorder Node -file floorDisp.out  -node $nodeList -dof 1 disp"; eval $cmd;

set tFinal [expr $dt*$nPt];  # maximum duration of ground-motion analysis
constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-6 10 
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient

set ok 0
set currentTime 0.0
while {$ok == 0 && $currentTime < $tFinal} {
    set ok [analyze 1 $dt]
    if {$ok != 0} {
	test NormDispIncr 1.0e-6 2000 1
	algorithm ModifiedNewton -initial
	set ok [analyze 1 $dt]
	test NormDispIncr 1.0e-6 10 
	algorithm Newton
    } 
    set currentTime [getTime]
}

remove recorders

# print out some results from recorded data
set a [open floorDispEnv.out r]
set line [gets $a]; set line [gets $a]; set line [gets $a]
puts "gMotion: $gMotion scale: $scale colLine 1 Max Node Disp: $line"
close $a

set a [open floorAccEnv.out r]
set line [gets $a]; set line [gets $a]; set line [gets $a]
#puts "MAX ACCEL: $line"
close $a

