# written: fmk
# units: kip & in

# load procedures in other files
source Steel2d.tcl
source ReadRecord.tcl;

# set some variables
set motion el_centro
#set motion Oak_2_50_5_FN
set PI 3.14159
set in 1.0;
set g 386.4;	

set roofWeight  [expr 80*120.*72./1000.]; #kips; 
set floorWeight [expr 95*120.*72./1000.];
 
set numFrameResisting 2.0; # lateral load resisting frames
set percentLoadFrame [expr 15./120.]

# set up my lists
set floorOffsets {216. 150. 150. 150. 150. 150.}
set colOffsets   {288. 288. 288.} 
set colSizes     {W30X173 W30X173 W27X146 W27X146 W24X104 W24X104};
set colExtSizes  {W14X193 W14X193 W14X159 W14X159 W14X109 W14X109};
set beamSizes    {W30X99 W30X99 W27X94 W27X94 W24X76 W24X76};

#calculated properties
set numFloor [expr [llength $floorOffsets]+1]
set numCline [expr [llength $colOffsets]+1]

for {set i 0; set width 0;} {$i < [expr $numCline-1]} {incr i 1} {set width [expr $width + [lindex $colOffsets $i]]}
set massAtFloorNode [expr $floorWeight/($g*$numFrameResisting*$numCline*1.0)]
set massAtRoofNode  [expr $roofWeight/($g*$numFrameResisting*$numCline*1.0)]
set uniformRoofLoad  [expr $roofWeight*$percentLoadFrame/$width]
set uniformFloorLoad [expr $floorWeight*$percentLoadFrame/$width]

# check of list dimensions for errors
if {[llength $colSizes] != [expr $numFloor-1]} {puts "ERROR: colSizes"; quit}
if {$numCline >= 10} {puts "ERROR: too many column lines, reprogram"; quit}
if {$numFloor >= 10} {puts "ERROR: too many floors, reprogram"; quit}

wipe;
model BasicBuilder -ndm 2 -ndf 3;  # Define the model builder, ndm = #dimension, ndf = #dofs

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
#uniaxialMaterial Steel02 1 50.0 29000 0.003 20 0.925 0.15
uniaxialMaterial Steel02 1 36.0 29000 0.003 20 0.925 0.15

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

# Gravity-analysis: load-controlled static analysis
set Tol 1.0e-6;    
constraints Plain;	
numberer RCM;		
system BandGeneral;	
test NormDispIncr $Tol 10;
algorithm Newton;		
integrator LoadControl 0.1;
analysis Static;		
analyze 10

# maintain constant gravity loads and reset time to zero
loadConst -time 0.0

constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-8 10 
algorithm Newton
integrator Newmark 0.5 0.25
#integrator HHT 0.9
analysis Transient


# add some damping
set pDamp 0.03
set lambda [eigen 3]
set omegaI [expr pow([lindex $lambda 0],0.5)];
set omegaJ [expr pow([lindex $lambda 2],0.5)];
set alphaM [expr $pDamp*(2*$omegaI*$omegaJ)/($omegaI+$omegaJ)];	
set betaKcomm [expr 2.*$pDamp/($omegaI+$omegaJ)]; 

#rayleigh $alphaM 0. 0. $betaKcomm
#rayleigh $alphaM $betaKcomm 0. 0.

set lambda [eigen 10]
set omegaI [expr pow([lindex $lambda 0],0.5)];
set alpha0 [expr 2.0 * $omegaI * $pDamp]
#rayleigh $alpha0 0. 0. 0.
#modalDamping  $pDamp

puts " Period T1: [expr 2*$PI/sqrt([lindex $lambda 0])] sec"

# create a load patrern for uniform excitation
ReadRecord $motion.AT2 $motion.g3 dt nPt;
timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor [expr 1. * $g]
pattern UniformExcitation 1 1 -accel 10;

set nodeList []
for {set floor 1} {$floor <= $numFloor} {incr floor 1} {
    lappend nodeList 1$floor
}
set cmd "recorder EnvelopeNode -file floorAccEnv.out  -timeSeries 10 -node $nodeList -dof 1 accel"; eval $cmd;

set cmd "recorder EnvelopeNode -file floorDispEnv.out  -node $nodeList -dof 1 disp"; eval $cmd;
set cmd "recorder Node -file floorDisp.out  -time -node $nodeList -dof 1 disp"; eval $cmd;
set cmd "recorder Node -file floorReaction.out  -time -node $nodeList -dof 3 reaction"; eval $cmd;


puts $nodeList

set tFinal	[expr $dt*$nPt];	# maximum duration of ground-motion analysis

set ok 0
set currentTime 0.0
while {$ok == 0 && $currentTime < $tFinal} {
    set ok [analyze 1 $dt]
    if {$ok != 0} {
	test NormDispIncr 1.0e-5 2000 1
	algorithm ModifiedNewton -initial
	set ok [analyze 1 $dt]
	test NormDispIncr 1.0e-6 10 
	algorithm Newton
    } 
    if {$ok != 0} {
	algorithm NewtonLineSearch 0.8
	set ok [analyze 1 $dt]
	algorithm Newton
    } 
    set currentTime [getTime]
}

remove recorders

puts "\n Results"
set a [open floorDispEnv.out r]
set line [gets $a]; set line [gets $a]; set line [gets $a]
puts "MAX DISP:  $line"
close $a

set a [open floorAccEnv.out r]
set line [gets $a]; set line [gets $a]; set line [gets $a]
#puts "MAX ACCEL: $line"
close $a

