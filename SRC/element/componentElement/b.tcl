logFile a.out

set anaType Pushover
set stiff 1e3;
set stiff2 1e3;
set b 0.1;
set Fy 2500.
foreach eleType {Hinge Component} {

    model basic -ndm 2 -ndf 3
    set height   100.
    node  1       0.0     0.0 
    node  3       0.0 $height

    fix   1     1    1    1

    geomTransf Linear 2
    uniaxialMaterial Steel01 2 $Fy $stiff $b
    if {$eleType == "Elastic"} {
	element elasticBeamColumn   1   1  3  100. 100. 1000.    2
    } elseif {$eleType == "Hinge"} {
	node 31 0. 0
	element elasticBeamColumn   1   31  3  100.  100. 1000.    2
	element zeroLength 33 1 31 -mat 2 -dir 6
	fix 31 1 1 0
    } else {
	uniaxialMaterial Elastic 1 $stiff2
	element componentElement2d   1   1   3 100. 100. 1000.  2 2 1
    }

    if {$anaType == "Pushover"} {

	timeSeries Linear 1
	pattern Plain 1 1 {
	    load 3 10. 0. 0.
#	    eleLoad 3 125 0 
#	    eleLoad -ele 1 -type beamUniform 10.
#	    eleLoad -ele 2 -type beamUniform 10.
       }

	system ProfileSPD
	constraints Plain
	test NormDispIncr 1.0e-6  3 4
	algorithm Newton 
	numberer RCM
	integrator LoadControl 1
	analysis Static
	analyze 5
	
    } else {
	set P 280.
	set g 386.4
	set m [expr $P/$g];       # expr command to evaluate an expression
	mass  3    $m   $m    0
	mass  4    $m   $m    0
	set outFile ARL360.g3
	
	# Source in TCL proc to read PEER SMD record
	source ReadSMDFile.tcl
	
	# Permform the conversion from SMD record to OpenSees record
	#              inFile     outFile dt
	ReadSMDFile ARL360.at2 $outFile dt
	timeSeries Path 1 -filePath $outFile -dt $dt -factor $g
	
	pattern UniformExcitation  2   1  -accel 1
	
	rayleigh 0.0 0.0 0.0 0.0
	
	wipeAnalysis
	
	system FullGeneral
	constraints Plain
	test NormDispIncr 1.0e-6  10 0
	algorithm Newton
	numberer RCM
	integrator Newmark  0.5  0.25 
	analysis Transient
	
	recorder Node -time -file disp.out -node 3 4 -dof 1 2 3 disp
	recorder Node -time -file a2.out -timeSeries 1 -node 3 4 -dof 1 accel
	recorder Node -time -file a1.out -node 3 4 -dof 1 2 3 accel
	set lambda [eigen 2]
	set PI 3.14159
	set T1 [expr 2*$PI/sqrt([lindex $lambda 0])] 
	if {$eleType == "Elastic" } {
	    set eE $T1
	} elseif {$eleType == "Hinge"} {
	    set hE $T1
	} else {
	    set cE $T1
	}
	
	
	set tFinal [expr 2000 * 0.01]
	#    set tFinal [expr 27 * 0.01]
	set tCurrent [getTime]
	set ok 0
	
	while {$ok == 0 && $tCurrent < $tFinal} {
	    
	    set ok [analyze 1 .01]

	    
	    # if the analysis fails try initial tangent iteration
	    if {$ok == 1000} {
		puts "regular newton failed .. lets try an initail stiffness for this step"
		test NormDispIncr 1.0e-6  1000 1
		algorithm Newton -initial
		set ok [analyze 1 .01]
		if {$ok == 0} {puts "that worked .. back to regular newton"}
		test NormDispIncr 1.0e-8  10 0
		algorithm Newton
	    }
	    
	    set tCurrent [getTime]
	}
	
	if {$ok == 0} {
	    puts "Transient analysis completed SUCCESSFULLY";
	} else {
	    puts "Transient analysis completed FAILED";    
	}
    }

    if {$eleType == "Elastic" } {
	set eDisp [nodeDisp 3]
	set eForce [eleResponse 1 forces]
	set eTime [getTime]
	print ele 1
    } elseif {$eleType == "Hinge"} {
	set hDisp [nodeDisp 3]
	set hForce {}
	lappend hForce [eleResponse 1 forces] 
	lappend hForce [eleResponse 33 forces]
	set hTime [getTime]
    } else {
	set cDisp [nodeDisp 3]
	set cForce [eleResponse 1 forces]
	set cTime [getTime]
    }
    wipe
}
#puts "eTime: $eTime"
puts "hTime: $hTime"
puts "cTime: $cTime"

puts "hDisp: $hDisp"
puts "cDisp: $cDisp"

#puts "eForce: $eForce"
puts "hForce: $hForce"
puts "cForce $cForce"


