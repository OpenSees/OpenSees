# load procedures in other files
source Steel2d.tcl
source ReadRecord.tcl;

# set some variables
#set motion el_centro
set motion Oak_2_50_5_FN
set in 1.0;
set g 38.4;				# acceleration due to gravity

# set up my lists
set floorOffsets {268. 160. 160. 160. 160. 160. 160.}
set colOffsets   {300. 300. 300. 300. 300.} 

set roofWeight 1537.0; set roofMass [expr $roofWeight/($g*2.*6)]; # kips 2 frames per dirn 6 col line
set floorWeight 1920.0; set floorMass [expr $floorWeight/($g*2.*6)]

set massesCMD  "set masses {0. $floorMass $floorMass $floorMass $floorMass $floorMass $floorMass $roofMass}"
eval $massesCMD
puts $massesCMD

set colSizes       {W24X146 W24X146 W24X76 W24X76 W14X61 W14X61 W24X55};
set beamSizes      {W24X207 W21X62 W21X50 W21X44 W21X44 W21X44 W21X44};

set numFloor [expr [llength $floorOffsets]+1]
set numCline [expr [llength $colOffsets]+1]


# build colLocations and floorLocations
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


# check of list dimensions for errors
if {[llength $masses] != $numFloor} {puts "ERROR: massX"; quit}
if {[llength $colSizes] != [expr $numFloor-1]} {puts "ERROR: colSizes"; quit}
if {[llength $beamSizes] != [expr $numFloor-1]} {puts "ERROR: beamSizes"; quit}

wipe;
model BasicBuilder -ndm 2 -ndf 3;  # Define the model builder, ndm = #dimension, ndf = #dofs

# Build the Nodes
for {set floor 1} {$floor <= $numFloor} {incr floor 1} {
    set floorLoc [lindex $floorLocations [expr $floor-1]]
    set mass [lindex $masses [expr $floor-1]]
    for {set colLine 1} {$colLine <= $numCline} {incr colLine 1} {
	set colLoc [lindex $colLocations [expr $colLine-1]]
	node $colLine$floor $colLoc $floorLoc -mass $mass $mass 0.
	if {$floor == 1} {
	    fix $colLine$floor 1 1 1
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
    for {set floor1 1} {$floor1 < $numFloor} {incr floor1 1} {
	set floor2 [expr $floor1+1]
	set theSection [lindex $colSizes [expr $floor1 -1]]
	DispBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection 1 1 -nip 5
    }
}

# build the beams
geomTransf Linear 2
for {set colLine1  1} {$colLine1 < $numCline} {incr colLine1 1} {
    set colLine2 [expr $colLine1 + 1]
    for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
	set theSection [lindex $beamSizes [expr $floor -2]]
	DispBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor $colLine2$floor $theSection 1 2
    }
}

# add uniform loads to beams
set floorLoad -0.11238
set roofLoad -0.1026
pattern Plain 101 Linear {
    for {set colLine1  1} {$colLine1 < $numCline} {incr colLine1 1} {
	set colLine2 [expr $colLine1 + 1]
	for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
	    if {$floor == 4} {
		eleLoad -ele $colLine1$floor$colLine2$floor -type beamUniform $roofLoad
	    } else {
		eleLoad -ele $colLine1$floor$colLine2$floor -type beamUniform $floorLoad
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
    draw
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

set a [open floorDispEnv.out r]
set line [gets $a]; set line [gets $a]; set line [gets $a]
puts "MAX DISP:  $line"
close $a

set a [open floorAccEnv.out r]
set line [gets $a]; set line [gets $a]; set line [gets $a]
puts "MAX ACCEL: $line"
close $a

