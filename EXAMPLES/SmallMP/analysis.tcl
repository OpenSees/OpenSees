
proc doGravity {} {
    system BandGeneral
    constraints Transformation
    numberer RCM
    test NormDispIncr 1.0e-12  10 3
    algorithm Newton
    integrator LoadControl 0.1
    
    analysis Static
    
    set ok [analyze 10]
    return $ok
}


proc doDynamic {dT nPts} {

    system BandGeneral
    constraints Plain
    test NormDispIncr 1.0e-12  10 
    algorithm Newton
    numberer RCM
    integrator Newmark  0.5  0.25 

    analysis Transient

    # set some variables
    set tFinal [expr $nPts * $dT]
    set tCurrent [getTime]
    set ok 0

    # Perform the transient analysis
    while {$ok == 0 && $tCurrent < $tFinal} {
	set ok [analyze 1 $dT]
	
	# if the analysis fails try initial tangent iteration
	if {$ok != 0} {
	    puts "regular newton failed .. lets try an initial stiffness for this step"
	    test NormDispIncr 1.0e-12  100 0
	    algorithm ModifiedNewton -initial
	    set ok [analyze 1 $dT]
	    if {$ok == 0} {puts "that worked .. back to regular newton"}
	    test NormDispIncr 1.0e-12  10 
	    algorithm Newton
	}
	
	set tCurrent [getTime]
    }
    
    return $ok
}

