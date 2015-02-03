
#  PDelta Effects 

# REFERENCES:
# R.C.Kaehler, D.W.White, Y.D.Kim, "Frame Design Using Web-Tapered Members", AISC 2011

puts "AISC - Design Guide 25 - Frame Design Using Web-Tapered Members"

set H  10.0
set L 196.0
set PI [expr 2.0*asin(1.0)]

set ok 0
set results []
set counter 0

puts "Prismatic Beam Benchmark Problems\n"
puts "    - Case 1 (Single Curvature)     - elasticBeamColumn"
puts "------+--------+-------------------------+-------------------------"
puts "      |        |     Tip Displacement    |      Base Moment        "
puts "------+--------+--------+---------+------+---------+--------+------"
set formatString {%5s|%8s|%8s|%9s|%6s|%9s|%8s|%6s}
puts [format $formatString numEle alpha Exact OpenSees %Error Exact OpenSees %Error]
puts "------+--------+--------+---------+------+---------+--------+------"

foreach numEle {1 2 4 10} {
    foreach alpha {0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.67} {
	wipe
	model Basic -ndm 2
	
	set E 29500.0
	set A 51.7
	set I 2150.0
	
	set Pel [expr $PI*$PI*$E*$I/($L*$L)]
	set Pcr [expr $Pel/4.0]
	set Pr  $Pcr
	
	set dY [expr $L/$numEle]
	for {set i 0} {$i <= $numEle} {incr i 1} {
	    node [expr $i +1] 0. [expr $i * $dY]
	}
	
	geomTransf PDelta 1
	section Elastic 1 $E $A $I
	set eleTag 1; set iNode 1; set jNode 2;
	for {set i 0} {$i < $numEle} {incr i 1} {
	    element elasticBeamColumn $eleTag $iNode $jNode $A $E $I 1
	    incr eleTag; incr iNode; incr jNode;
	}
	
	fix 1 1 1 1
	
	set u    [expr $PI/2.0*sqrt($alpha * $Pr/$Pel)]
	
	if {$u != 0} {
	    set resU [expr ($H*$L*$L*$L/(3.0*$E*$I)) * (3.0*(tan(2.0*$u)-2.0*$u)/(8.0*$u*$u*$u))]
	    set resM [expr $H*$L*(tan(2.*$u)/(2.*$u))]
	} else {
	    set resU [expr ($H*$L*$L*$L/(3.0*$E*$I))]
	    set resM [expr $H*$L]
	}
	
	timeSeries Linear 1
	pattern Plain 1 1 {
	    load  [expr $numEle+1] $H [expr -$alpha*$Pr] 0.
	}
	
	constraints Plain
	system ProfileSPD
	numberer Plain
	integrator LoadControl 1
	test NormDispIncr 1.0e-12 6 0
	algorithm Newton
	analysis Static
	analyze 1
	
	set delta [nodeDisp [expr $numEle + 1] 1]
	set moment [lindex [eleResponse 1 forces] 2]
	set formatString {%6.0f|%8.2f|%8.4f|%9.4f|%6.1f|%9.2f|%8.2f|%6.1f}
	puts [format $formatString $numEle $alpha $resU $delta [expr 100*($resU-$delta)/$delta] $resM $moment [expr 100*($resM-$moment)/$moment] ]
    }
    puts "------+--------+--------+---------+------+---------+--------+------"
}

# test on last one
if {[expr abs(100*($resU-$delta)/$delta)] > 0.5 || [expr abs(100*($resM-$moment)/$moment)] > 0.5} {
    set ok 1
    puts "[expr abs(100*($resU-$delta)/$delta)] > 0.5 || [expr abs(100*($resM-$moment)/$moment)] > 0.5"
}

puts "Prismatic Beam Benchmark Problems\n"
puts "    - Case 1 (Single Curvature) - forceBeamColumnCBDI element"
puts "------+--------+-------------------------+-------------------------"
puts "      |        |     Tip Displacement    |      Base Moment        "
puts "------+--------+--------+---------+------+---------+--------+------"
set formatString {%5s|%8s|%8s|%9s|%6s|%9s|%8s|%6s}
puts [format $formatString numEle alpha Exact OpenSees %Error Exact OpenSees %Error]
puts "------+--------+--------+---------+------+---------+--------+------"
foreach numEle {1 } {
    foreach alpha {0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.67} {
	wipe
	model Basic -ndm 2
	
	set E 29500.0
	set A 51.7
	set I 2150.0
	
	set Pel [expr $PI*$PI*$E*$I/($L*$L)]
	set Pcr [expr $Pel/4.0]
	set Pr  $Pcr
	
	set dY [expr $L/$numEle]
	for {set i 0} {$i <= $numEle} {incr i 1} {
	    node [expr $i +1] 0. [expr $i * $dY]
	}
	
	geomTransf PDelta 1
	section Elastic 1 $E $A $I
	set eleTag 1; set iNode 1; set jNode 2;
	for {set i 0} {$i < $numEle} {incr i 1} {
	    element forceBeamColumnCBDI $eleTag $iNode $jNode 1 "Legendre 1 4"
	    incr eleTag; incr iNode; incr jNode;
	}
	
	fix 1 1 1 1
	
	set u    [expr $PI/2.0*sqrt($alpha * $Pr/$Pel)]
	
	if {$u != 0} {
	    set resU [expr ($H*$L*$L*$L/(3.0*$E*$I)) * (3.0*(tan(2.0*$u)-2.0*$u)/(8.0*$u*$u*$u))]
	    set resM [expr $H*$L*(tan(2.*$u)/(2.*$u))]
	} else {
	    set resU [expr ($H*$L*$L*$L/(3.0*$E*$I))]
	    set resM [expr $H*$L]
	}
	
	timeSeries Linear 1
	pattern Plain 1 1 {
	    load  [expr $numEle+1] $H [expr -$alpha*$Pr] 0.
	}
	
	constraints Plain
	system ProfileSPD
	numberer Plain
	integrator LoadControl 1
	test NormDispIncr 1.0e-12 6 0
	algorithm Newton
	analysis Static
	analyze 1
	
	set delta [nodeDisp [expr $numEle + 1] 1]
	set moment [lindex [eleResponse 1 forces] 2]
	set formatString {%6.0f|%8.2f|%8.4f|%9.4f|%6.1f|%9.2f|%8.2f|%6.1f}
	puts [format $formatString $numEle $alpha $resU $delta [expr 100*($resU-$delta)/$delta] $resM $moment [expr 100*($resM-$moment)/$moment] ]
    }
    puts "------+--------+--------+---------+------+---------+--------+------"
}


# test on last one
if {[expr abs(100*($resU-$delta)/$delta)] > 0.5 || [expr abs(100*($resM-$moment)/$moment)] > 0.5} {
    set ok 1
    puts "[expr abs(100*($resU-$delta)/$delta)] > 0.5 || [expr abs(100*($resM-$moment)/$moment)] > 0.5"
}


puts "\n\n    - Case 2 (Double Curvature)  - elasticBeamColumn"
puts "------+--------+-------------------------+-------------------------"
puts "      |        |     Tip Displacement    |      Base Moment        "
puts "------+--------+--------+---------+------+---------+--------+------"
set formatString {%5s|%8s|%8s|%9s|%6s|%9s|%8s|%6s}
puts [format $formatString numEle alpha Exact OpenSees %Error Exact OpenSees %Error]
puts "------+--------+--------+---------+------+---------+--------+------"
foreach numEle {1 2 10} {
    foreach alpha {0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.67} {
	wipe
	model Basic -ndm 2
	
	set E 29500.0
	set A 51.7
	set I 2150.0
	
	set Pel [expr $PI*$PI*$E*$I/($L*$L)]
	set Pcr [expr $Pel/4.0]
	set Pr  $Pcr
	
	set dY [expr $L/$numEle]
	for {set i 0} {$i <= $numEle} {incr i 1} {
	    node [expr $i +1] 0. [expr $i * $dY]
	}
	
	geomTransf PDelta 1
	set eleTag 1; set iNode 1; set jNode 2;
	for {set i 0} {$i < $numEle} {incr i 1} {
	    element elasticBeamColumn $eleTag $iNode $jNode $A $E $I 1
	    incr eleTag; incr iNode; incr jNode;
	}
	
	fix 1 1 1 1
	fix [expr $numEle+1] 0 0 1
	
	set u    [expr $PI/2.0*sqrt($alpha * $Pr/$Pel)]
	
	if {$u != 0} {
	    set resU [expr ($H*$L*$L*$L/(12.0*$E*$I)) * (3.0*(tan($u)-$u)/($u*$u*$u))]
	    set resM [expr $H*$L/2.0*(tan($u)/($u))]
	} else {
	    set resU [expr ($H*$L*$L*$L/(12.0*$E*$I))]
	    set resM [expr $H*$L/2.0]
	}
	
	timeSeries Linear 1
	pattern Plain 1 1 {
	    load  [expr $numEle+1] $H [expr -$alpha*$Pr] 0.
	}
	
	constraints Plain
	system ProfileSPD
	numberer Plain
	integrator LoadControl 1
	test NormDispIncr 1.0e-12 6 0
	algorithm Newton
	analysis Static
	analyze 1
	
	set delta [nodeDisp [expr $numEle + 1] 1]
	set moment [lindex [eleResponse 1 forces] 2]
	set formatString {%6.0f|%8.2f|%8.4f|%9.4f|%6.1f|%9.2f|%8.2f|%6.1f}
	puts [format $formatString $numEle $alpha $resU $delta [expr 100*($resU-$delta)/$delta] $resM $moment [expr 100*($resM-$moment)/$moment] ]
    }
    puts "------+--------+--------+---------+------+---------+--------+------"
}

# test on last one
if {[expr abs(100*($resU-$delta)/$delta)] > 0.5 || [expr abs(100*($resM-$moment)/$moment)] > 0.5} {
    set ok 1
    puts "[expr abs(100*($resU-$delta)/$delta)] > 0.5 || [expr abs(100*($resM-$moment)/$moment)] > 0.5"
}


set results [open results.out a+]
if {$ok == 0} {
    puts "PASSED Verification Test AISC25.tcl \n\n"
    puts $results "PASSED : AISC25.tcl"
} else {
    puts "FAILED Verification Test AISC25.tcl \n\n"
    puts $results "FAILED : AISC25.tcl"
}

close $results
