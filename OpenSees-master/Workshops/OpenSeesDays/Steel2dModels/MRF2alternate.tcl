# load procedures in other files
source Steel2d.tcl
source ReadRecord.tcl;

# set some variables
set motion el_centro
set in 1.0;
set g 600.4;				# acceleration due to gravity

foreach beamType {ElasticBeam ElasticForceBeam ForceBeam DispBeam BeamWithConcentratedHinge} {
    foreach colType {ElasticBeam ForceBeam DispBeam BeamWithFiberHinge} {

	puts "\ncolType: $colType beamType: $beamType"


	wipe;
	model BasicBuilder -ndm 2 -ndf 3;  # Define the model builder, ndm = #dimension, ndf = #dofs
	
	
	set floorLocations {0 204. 384. 564.}
	set colLocations   {0. 360. 720. 1080. 1440. 1800.} 
	
	set massesX        {0. 0.419 0.419 0.400}
	set massesY        {0. 0.105 0.105 0.096}
	set colSizes       {W14X370 W14X370 W14X211};
	set beamSizes      {W33X141 W33X130 W27X102};
	
	set numFloor [llength $floorLocations]
	set numCline [llength $colLocations]
	
	# check of list dimensions
	if {[llength $massesX] != $numFloor} {puts "ERROR: massX"; quit}
	if {[llength $massesY] != $numFloor} {puts "ERROR: massY"; quit}
	if {[llength $colSizes] != [expr $numFloor-1]} {puts "ERROR: colSizes"; quit}
	if {[llength $beamSizes] != [expr $numFloor-1]} {puts "ERROR: beamSizes"; quit}
	
	# Build the Nodes
	for {set floor 1} {$floor <= $numFloor} {incr floor 1} {
	    set floorLoc [lindex $floorLocations [expr $floor-1]]
	    set massX [lindex $massesX [expr $floor-1]]
	    set massY [lindex $massesY [expr $floor-1]]
	    for {set colLine 1} {$colLine <= $numCline} {incr colLine 1} {
		set colLoc [lindex $colLocations [expr $colLine-1]]
		node $colLine$floor $colLoc $floorLoc -mass $massX $massY 0.
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
		if {$colType == "ForceBeam"} {
		    ForceBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection 1 1 -nip 5
		} elseif {$colType == "DispBeam"} {
		    DispBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection 1 1 -nip 5 -nFlange 10 -nWeb 10
		} elseif {$colType == "ElasticBeam"} {
		    ElasticBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection $Es 1 
		} elseif {$colType == "ElasticForceBeam"} {		
		    ForceBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection 1 1 -nip 5 -nFlange 10 -nWeb 10 -elasticSection $Es 
		} elseif {$colType == "BeamWithHinges"} {		
		    BeamWithHingesWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection 1 1 -nip 5 -nFlange 10 -nWeb 10 -elasticSection $Es 
		}
	    }
	}
	
	# build the beams
	geomTransf Linear 2
	for {set colLine1  1} {$colLine1 < $numCline} {incr colLine1 1} {
	    set colLine2 [expr $colLine1 + 1]
	    for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
		set theSection [lindex $beamSizes [expr $floor -2]]
		
		if {$beamType == "ForceBeam"} {
		    ForceBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor $colLine2$floor $theSection 1 2 -nFlange 10 -nWeb 10
		} elseif {$beamType == "DispBeam"} {
		    DispBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor $colLine2$floor $theSection 1 2
		} elseif {$beamType == "ElasticBeam"} {
		    ElasticBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor $colLine2$floor $theSection $Es  2
		} elseif {$beamType == "ElasticForceBeam"} {		
		    ForceBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor $colLine2$floor $theSection 1 2 -nFlange 10 -nWeb 10 -elasticSection $Es
		}
		} elseif {$beamType == "BeamWithConcentratedHinge"} {		
		    BeamWithPlasticHingesWSection2d $colLine1$floor$colLine2$floor $colLine1$floor $colLine2$floor $theSection $Es 2 
		}
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

	puts "omegaI: $omegaI"
	
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
    }
}
