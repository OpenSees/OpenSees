# load procedures in other files
source Steel2d.tcl
source ReadRecord.tcl;

# set some variables
#set motion el_centro
set motion Oak_2_50_5_FN
set in 1.0;
set g 386.4;				# acceleration due to gravity


foreach nFlange {1 2 3 10 10 20 35} nWeb {4 4 4 5 10 20 30} {
    
    puts "nFlange: $nFlange nWeb: $nWeb"

    wipe;
    model BasicBuilder -ndm 2 -ndf 3;  # Define the model builder, ndm = #dimension, ndf = #dofs
    
    # set some lists containing floor and col line locations and nodal masses
    set floorLocations {0 204. 384. 564.}
    set colLocations   {0. 360. 720. 1080. 1440. 1800.} 
    set massesX        {0. 0.419 0.419 0.400}
    set massesY        {0. 0.105 0.105 0.096}
    
    # add nodes at each floor at each column line location & fix nodes if floor 1
    foreach floor {1 2 3 4} floorLoc $floorLocations massX $massesX massY $massesY {
	foreach colLine {1 2 3 4 5 6} colLoc $colLocations {
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
    
    # set some lists for col and beam sizes
    
    set colSizes  {W14X370 W14X370 W14X211};  #col story 1, 2 and 3
    set beamSizes {W33X141 W33X130 W27X102};  #beams floor 1, 2, and 3
    
    # add column at each column line between all floors
    geomTransf PDelta 1
    foreach colLine {1 2 3 4 5 6} {
	foreach floor1 {1 2 3} floor2 { 2 3 4} {
	    set theSection [lindex $colSizes [expr $floor1 -1]]
	    ForceBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection 1 1 -nip 5 -nFlange $nFlange -nWeb $nWeb
	}
    }
    
    # add beams between column lines at each floor
    geomTransf Linear 2
    foreach colLine1 {1 2 3 4 5} colLine2 {2 3 4 5 6} {
	foreach floor {2 3 4} {
	    set theSection [lindex $beamSizes [expr $floor -2]]
	    ForceBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor $colLine2$floor $theSection 1 2 -nFlange $nFlange -nWeb $nWeb
	}
    }
    
    # add uniform loads to beams
    set floorLoad -0.11238
    set roofLoad -0.1026
    pattern Plain 101 Linear {
	foreach colLine1 {1 2 3 4 5} colLine2 {2 3 4 5 6} {
	    foreach floor {2 3 4} {
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
    
    recorder EnvelopeNode -file floorAccEnv.out  -timeSeries 10 -node 11 12 13 14 -dof 1 accel
    recorder Node -file floorAcc.out -timeSeries 10 -node 11 12 13 14 -dof 1 accel	
    
    recorder EnvelopeNode -file floorDispEnv.out  -node 11 12 13 14 -dof 1 disp
    recorder Node -file floorDisp.out  -node 11 12 13 14 -dof 1 disp	
    
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

