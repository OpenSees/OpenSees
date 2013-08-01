
set numP [getNP]
set pid  [getPID]
set count 0;


set numMat [lindex $argv 0]
set numMotion [lindex $argv 1]

for {set matID 1} {$matID<=$numMat} {incr matID} {
    for {set motionID 1} {$motionID<=$numMotion} {incr motionID} {
	if {[expr $count % $numP] == $pid}  {

	    puts "$pid DOING mat: $matID motion: $motionID"
	    model BasicBuilder -ndm 2 -ndf 3

	    source matProperties.tcl.$matID
	    set nodes [source model.tcl]
	    source motion.tcl.$motionID
	    set node1 [lindex $nodes 0]
	    set a "recorder Node -file NodeDisp.out.$matID.$motionID -time -node $node1 -dof 1 disp"
	    eval $a
	    set a "recorder Node -file NodeAccel.out.$matID.$motionID -time -timeSeries 2 -node $node1 -dof 1 accel"
	    eval $a
	    
	    system BandGeneral
	    constraints Plain
	    test NormDispIncr 1.0e-12  10 
	    algorithm Newton
	    numberer RCM
	    integrator Newmark  0.5  0.25 
	    analysis Transient
	    
	    
	    set t 0.0
	    set maxT [expr (1+$nPts)*$dT];
	    set ok 0.0
	    set maxD 0.0
	    
	    while {$ok == 0 && $t < $maxT} {
		set ok [analyze 1 $dT]
		if {$ok != 0} {
		    test NormDispIncr 1.0e-8  10 1
		    algorithm ModifiedNewton -initial
		    set ok [analyze 1 $dT]
		    test NormDispIncr 1.0e-8  10 
		    algorithm Newton
		}
		
		set t [getTime]
	    }

	    wipe
	}
	incr count 1
    }
}
