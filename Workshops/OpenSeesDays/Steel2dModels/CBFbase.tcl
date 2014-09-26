# EXAMPLE inverted V BRACING

# load procedures in other files
source Steel2d.tcl
source ReadRecord.tcl;

# set up my structure
set in 1
set floorOffsets {180. 180. 180.}
set colOffsets   {180. 180.} 
set masses      {0. 0.419 0.419 0.400}

set colSizes     {W14X176 W14X176 W14X176};
set beamSizes    {W36X210  W30X116 W27X84};
set braceSizes   {HSS12X12X1/2 HSS12X12X1/2 HSS10X10X3/8}

set floorLoad -0.11238
set roofLoad -0.1026

# build colLocations and floorLocations & set some variables
set numFloor [expr [llength $floorOffsets]+1]
set numStory [expr $numFloor-1]
set numCline [expr [llength $colOffsets]+1]

# check of list dimensions for errors
if {[llength $masses] != $numFloor} {puts "ERROR: massX"; quit}
if {[llength $colSizes] != [expr $numFloor-1]} {puts "ERROR: colSizes"; quit}
if {[llength $beamSizes] != [expr $numFloor-1]} {puts "ERROR: beamSizes"; quit}
if {[llength $braceSizes] != [expr $numFloor-1]} {puts "ERROR: beamSizes"; quit}

wipe;
model BasicBuilder -ndm 2 -ndf 3;  # Define the model builder, ndm = #dimension, ndf = #dofs

# Build the Nodes
for {set floor 1; set floorLoc 0.} {$floor <= $numFloor} {incr floor 1} {
    set mass [lindex $masses [expr $floor-1]]
    for {set colLine 1; set colLoc 0.} {$colLine <= $numCline} {incr colLine 1} {
	node $colLine$floor $colLoc $floorLoc -mass $mass $mass 0.
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
	    
# define material for col and beams
set Es 29000.0;  # modulus of elasticity for steel
set Fy 50.0; 	 # yield stress of steel
set b 0.003;	 # strain hardening ratio
uniaxialMaterial Steel02 1 $Fy $Es $b 20 0.925 0.15 0.0005 0.01 0.0005 0.01

# add the columns
geomTransf PDelta 1
for {set colLine 1} {$colLine <= $numCline} {incr colLine 2} {
    for {set floor1 1; set floor2 2} {$floor1 < $numFloor} {incr floor1 1; incr floor2 1} {
	set theSection [lindex $colSizes [expr $floor1 -1]]
	ForceBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection 1 1 -nip 5
    }
}

# add the beams, pinned connection at column end
geomTransf Linear 2
for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
    set colLine1 1; set colLine2 2; set colLine3 3;
    set theSection [lindex $beamSizes [expr $floor -2]]
    DispBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor $colLine2$floor $theSection 1 2 -release1
    DispBeamWSection2d $colLine2$floor$colLine3$floor $colLine2$floor $colLine3$floor $theSection 1 2 -release2
}


# add uniform loads to beams
pattern Plain 101 Linear {
    for {set colLine1  1} {$colLine1 < $numCline} {incr colLine1 1} {
	set colLine2 [expr $colLine1 + 1]
	for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
	    if {$floor == $numFloor} {
		eleLoad -ele $colLine1$floor$colLine2$floor -type beamUniform $roofLoad
	    } else {
		eleLoad -ele $colLine1$floor$colLine2$floor -type beamUniform $floorLoad
	    }
	}
    }
}

# define material for braces
set Fy_b 46.0; 	 
set E0 0.095
set m -0.3
uniaxialMaterial Steel02 2 $Fy_b $Es $b 20 0.925 0.15 0.0005 0.01 0.0005 0.01
uniaxialMaterial Fatigue 3 2 -E0 $E0 -m $m -min -1.0 -max 0.04

set imperfection 0.001

# proc HSSbrace {eleTag iNode jNode secType matTag numSeg Im transfTag args}





