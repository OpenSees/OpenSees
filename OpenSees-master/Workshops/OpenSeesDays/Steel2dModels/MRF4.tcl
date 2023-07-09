# load procedures in other files
source Steel2d.tcl
source ReadRecord.tcl;

# set some variables
#set motion el_centro
set motion Oak_2_50_5_FN
set in 1.0;
set g 386.4;				# acceleration due to gravity


# set up my structure
set floorOffsets {204. 180. 180.}
set colOffsets   {360. 360. 360. 360. 360.} 
set massesX      {0. 0.419 0.419 0.400}
set massesY      {0. 0.105 0.105 0.096}

set colSizes     {W14X370 W14X370 W14X211};
set beamSizes    {W33X141 W33X130 W27X102};

set floorLoad -0.11238
set roofLoad -0.1026

# build colLocations and floorLocations & set some variables
set numFloor [expr [llength $floorOffsets]+1]
set numCline [expr [llength $colOffsets]+1]

set floorLocations 0; set floorLoc 0; 
set colLocations 0; set colLoc 0; 

for {set i 1} {$i < $numFloor} {incr i 1} {
    set floorLoc [expr $floorLoc + [lindex $floorOffsets [expr $i-1]]]
    lappend floorLocations $floorLoc;
}
for {set i 1} {$i < $numCline} {incr i 1} {
    set colLoc [expr $colLoc + [lindex $colOffsets [expr $i-1]]]
    lappend colLocations $colLoc;
}

# following in case num floors or cols exceed 10
if {$numFloor < 10} {
    set floorStart 1;
} elseif {$numFloor < 100} {
    set floorStart 10;
} else {
    set floorStart 100;
}

if {$numCline < 10} {
    set clineStart 1;
} elseif {$numCline < 100} {
    set clineStart 10;
} else {
    set clineStart 100;
}

# check of list dimensions for errors
if {[llength $massesX] != $numFloor} {puts "ERROR: massX"; quit}
if {[llength $massesY] != $numFloor} {puts "ERROR: massY"; quit}
if {[llength $colSizes] != [expr $numFloor-1]} {puts "ERROR: colSizes"; quit}
if {[llength $beamSizes] != [expr $numFloor-1]} {puts "ERROR: beamSizes"; quit}

wipe;
model BasicBuilder -ndm 2 -ndf 3;  # Define the model builder, ndm = #dimension, ndf = #dofs

# Build the Nodes
for {set floor 1} {$floor <= $numFloor} {incr floor 1} {
    set floorN [expr $floorStart + $floor-1]
    set floorLoc [lindex $floorLocations [expr $floor-1]]
    set massX [lindex $massesX [expr $floor-1]]
    set massY [lindex $massesY [expr $floor-1]]
    for {set colLine 1} {$colLine <= $numCline} {incr colLine 1} {
	set colLineN [expr $clineStart + $colLine -1]
	set colLoc [lindex $colLocations [expr $colLine-1]]
	node $colLineN$floorN $colLoc $floorLoc -mass $massX $massY 0.
	if {$floor == 1} {
	    fix $colLineN$floorN 1 1 1
	}
    }
}
	    
# define material 
set Es 29000.0;  # modulus of elasticity for steel
set Fy 50.0; 	 # yield stress of steel
set b 0.003;	 # strain hardening ratio
uniaxialMaterial Steel02 1 $Fy $Es $b 20 0.925 0.15

# build the columns
geomTransf PDelta 1
for {set colLine 1} {$colLine <= $numCline} {incr colLine 1} {
    set colLineN [expr $clineStart + $colLine -1]
    for {set floor1 1} {$floor1 < $numFloor} {incr floor1 1} {
	set floor2 [expr $floor1+1]
	set floor1N [expr $floorStart + $floor1-1]
	set floor2N [expr $floorStart + $floor1]
	set theSection [lindex $colSizes [expr $floor1 -1]]
	ForceBeamWSection2d $colLineN$floor1N$colLineN$floor2N $colLineN$floor1N $colLineN$floor2N $theSection 1 1 -nip 5
    }
}

# build the beams
geomTransf Linear 2
for {set colLine1  1} {$colLine1 < $numCline} {incr colLine1 1} {
    set colLine2 [expr $colLine1 + 1]
    set colLine1N [expr $clineStart + $colLine1 -1]
    set colLine2N [expr $clineStart + $colLine1]
    for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
	set floorN [expr $floorStart + $floor-1]
	set theSection [lindex $beamSizes [expr $floor -2]]
	ForceBeamWSection2d $colLine1N$floorN$colLine2N$floorN $colLine1N$floorN $colLine2N$floorN $theSection 1 2
    }
}

# add uniform loads to beams
pattern Plain 101 Linear {
    for {set colLine1  1} {$colLine1 < $numCline} {incr colLine1 1} {
	set colLine2 [expr $colLine1 + 1]
	set colLine1N [expr $clineStart + $colLine1 -1]
	set colLine2N [expr $clineStart + $colLine1]
	for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
	    set floorN [expr $floorStart + $floor-1]
	    if {$floor == $numFloor} {
		eleLoad -ele $colLine1N$floorN$colLine2N$floorN -type beamUniform $roofLoad
	    } else {
		eleLoad -ele $colLine1N$floorN$colLine2N$floorN -type beamUniform $floorLoad
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

# add some damping
set pDamp 0.03
set lambda [eigen 3]
set omegaI [expr pow([lindex $lambda 0],0.5)];
set omegaJ [expr pow([lindex $lambda 2],0.5)];
set alphaM [expr $pDamp*(2*$omegaI*$omegaJ)/($omegaI+$omegaJ)];	
set betaKcomm [expr 2.*$pDamp/($omegaI+$omegaJ)]; 
rayleigh $alphaM 0. 0. $betaKcomm

# create a load patrern for uniform excitation
ReadRecord $motion.AT2 $motion.g3 dt nPt;
timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor $g
pattern UniformExcitation 1 1 -accel 10;


set nodeList []
for {set floor 1} {$floor <= $numFloor} {incr floor 1} {
    lappend nodeList 1$floor
}
set cmd "recorder EnvelopeNode -file floorAccEnv.out  -timeSeries 10 -node $nodeList -dof 1 accel"; eval $cmd;
set cmd "recorder EnvelopeNode -file floorAccEnv.out  -timeSeries 10 -node $nodeList -dof 1 accel"; eval $cmd;

set cmd "recorder EnvelopeNode -file floorDispEnv.out  -node $nodeList -dof 1 disp"; eval $cmd;
set cmd "recorder Node -file floorDisp.out  -node $nodeList -dof 1 disp"; eval $cmd;

set tFinal	[expr $dt*$nPt];	# maximum duration of ground-motion analysis
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

wipe

set a [open floorDispEnv.out r]
set line [gets $a]; set line [gets $a]; set line [gets $a]
puts "MAX DISP:  $line"
close $a

set a [open floorAccEnv.out r]
set line [gets $a]; set line [gets $a]; set line [gets $a]
puts "MAX ACCEL: $line"
close $a

