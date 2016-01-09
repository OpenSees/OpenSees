# This Extends EigenFrame.tcl verification test to:
#   1) run different element options to test: ForceBeamColumn, DspBeamColumn, ElasticSection and FiberSection2d.
#   2) run different solver options to test:

# REFERENCES
#   as per EigenFrame.tcl

puts "EigenFrame.Extra.tcl: Verification 2d Bathe & Wilson original Elastic Frame - with other options"
puts "  - eigenvalue "

set eleTypes {elasticBeam forceBeamElasticSection dispBeamElasticSection forceBeamFiberSectionElasticMaterial dispBeamFiberSectionElasticMaterial}

foreach eleType $eleTypes {

    wipe

    model Basic -ndm 2
    
    #    units kip, ft                                                                                                                              
    
    # properties  
    set bayWidth 20.0;
    set storyHeight 10.0;
    
    set numBay 10
    set numFloor 9

    set A 3.0;         #area = 3ft^2    
    set E 432000.0;    #youngs mod = 432000 k/ft^2  
    set I 1.0;         #second moment of area I=1ft^4       
    set M 3.0;         #mas/length = 4 kip sec^2/ft^2           
    set coordTransf "Linear";  # Linear, PDelta, Corotational
    set massType "-lMass";  # -lMass, -cMass


    set nPts 3;      # numGauss Points

    # an elastic material
    uniaxialMaterial Elastic 1 $E

    # an elastic section
    section Elastic 1 $E $A $I
        
    # a fiber section with A=3 and I = 1 (b=1.5, d=2) 2d bending about y-y axis
 #   set b 1.5; set d 2.0;
    set y 2.0; set z 1.5;
    set numFiberY 2000;  # note we only need so many to get the required accuracy on eigenvalue 1e-7!
    set numFiberZ 1;
    section Fiber 2 {
#	patch rect 1 $numFiberY $numFiberZ 0.0 0.0 $z $y
	patch quad 1 $numFiberY $numFiberZ [expr -$y/2.0] [expr -$z/2.0] [expr $y/2.0] [expr -$z/2.0] [expr $y/2.0] [expr $z/2.0] [expr -$y/2.0] [expr $z/2.0]
    }
    

    # add the nodes         
    #  - floor at a time    
    set nodeTag 1
    set yLoc 0.
    for {set j 0} {$j <= $numFloor} {incr j 1} {
	set xLoc 0.
	for {set i 0} {$i <=$numBay} {incr i 1} {
	    node $nodeTag $xLoc $yLoc
	    set xLoc [expr $xLoc + $bayWidth]
	    incr nodeTag 1
	}
	set yLoc [expr $yLoc + $storyHeight]
    }
    
    # fix base nodes        
    for {set i 1} {$i <= [expr $numBay+1]} {incr i 1} {
	fix $i 1 1 1
    }
    
    # add column element    
    geomTransf $coordTransf 1
    set eleTag 1
    for {set i 0} {$i <=$numBay} {incr i 1} {
	set end1 [expr $i+1]
	set end2 [expr $end1 + $numBay +1]
	for {set j 0} {$j<$numFloor} {incr j 1} {

	    if {$eleType == "elasticBeam"} {
		element elasticBeamColumn $eleTag $end1 $end2 $A $E $I 1 -mass $M $massType
	    } elseif {$eleType == "forceBeamElasticSection"} {
		element forceBeamColumn $eleTag $end1 $end2 $nPts 1 1 -mass $M
	    } elseif {$eleType == "dispBeamElasticSection"} {
		element dispBeamColumn $eleTag $end1 $end2 $nPts 1 1 -mass $M $massType
	    } elseif {$eleType == "forceBeamFiberSectionElasticMaterial"} {
		element forceBeamColumn $eleTag $end1 $end2 $nPts 2 1 -mass $M
	    } elseif {$eleType == "dispBeamFiberSectionElasticMaterial"} {
		element dispBeamColumn $eleTag $end1 $end2 $nPts 2 1 -mass $M $massType
	    } else {
		puts "BARF"
	    }
	    set end1 $end2
	    set end2 [expr $end1 + $numBay +1]
	    incr eleTag 1
	}
    }
    

    # add beam elements     
    for {set j 1} {$j<=$numFloor} {incr j 1} {
	set end1 [expr ($numBay+1)*$j+1]
	set end2 [expr $end1 + 1]
	for {set i 0} {$i <$numBay} {incr i 1} {
	    if {$eleType == "elasticBeam"} {
		element elasticBeamColumn $eleTag $end1 $end2 $A $E $I 1 -mass $M $massType
	    } elseif {$eleType == "forceBeamElasticSection"} {
		element forceBeamColumn $eleTag $end1 $end2 $nPts 1 1 -mass $M
	    } elseif {$eleType == "dispBeamElasticSection"} {
		element dispBeamColumn $eleTag $end1 $end2 $nPts 1 1 -mass $M $massType
	    } elseif {$eleType == "forceBeamFiberSectionElasticMaterial"} {
		element forceBeamColumn $eleTag $end1 $end2 $nPts 2 1 -mass $M
	    } elseif {$eleType == "dispBeamFiberSectionElasticMaterial"} {
		element dispBeamColumn $eleTag $end1 $end2 $nPts 2 1 -mass $M $massType
	    } else {
		puts "BARF"
	    }
#	    element elasticBeamColumn $eleTag $end1 $end2 $A $E $I 1 -mass $M
	    set end1 $end2
	    set end2 [expr $end1 + 1]
	    incr eleTag 1
	}
    }
    
    # calculate eigenvalues
    set numEigen 3
    set eigenValues [eigen $numEigen]
    set PI [expr 2*asin(1.0)]
    
    # determine PASS/FAILURE of test
    set testOK 0

    # print table of camparsion
    #                         Bathe & Wilson               Peterson                    SAP2000                  SeismoStruct

    set comparisonResults {{0.589541 5.52695 16.5878} {0.589541 5.52696 16.5879} {0.589541 5.52696 16.5879} {0.58955 5.527 16.588}}
    puts "\n\nEigenvalue Comparisons for eleType: $eleType"
    set tolerances {9.99e-7 9.99e-6 9.99e-5}
    set formatString {%15s%15s%15s%15s%15s}
    puts [format $formatString OpenSees Bathe&Wilson Peterson SAP2000 SeismoStruct]
    set formatString {%15.5f%15.4f%15.4f%15.4f%15.3f}
    for {set i 0} {$i<$numEigen} {incr i 1} {
	set lambda [lindex $eigenValues $i]
	puts [format $formatString $lambda  [lindex [lindex $comparisonResults 0] $i] [lindex [lindex $comparisonResults 1] $i] [lindex [lindex $comparisonResults 2] $i] [lindex [lindex $comparisonResults 3] $i]]
	set resultOther [lindex [lindex $comparisonResults 2] $i]
	set tol [lindex $tolerances $i]
	if {[expr abs($lambda-$resultOther)] > $tol} {
	    set testOK -1;
	    puts "failed-> [expr abs($lambda-$resultOther)] $tol"
	}
    }
    
    set results [open results.out a+]
    if {$testOK == 0} {
	puts "PASSED Verification Test EigenFrame.Extra.tcl $eleType  \n\n"
	puts $results "PASSED : EigenFrame.Extra.tcl eleType: $eleType"
    } else {
	puts "FAILED Verification Test EigenFrame.Extra.tcl $eleType  \n\n"
	puts $results "FAILED : EigenFrame.Extra.tcl eleType: $eleType"
    }
    close $results
}



set solverTypes {-genBandArpack -fullGenLapack -UmfPack -SuperLU -ProfileSPD}

foreach solverType $solverTypes {

    set eleType elasticBeam

    wipe

    model Basic -ndm 2
    
    #    units kip, ft                                           
    
    # properties  
    set bayWidth 20.0;
    set storyHeight 10.0;
    
    set numBay 10
    set numFloor 9

    set A 3.0;         #area = 3ft^2    
    set E 432000.0;    #youngs mod = 432000 k/ft^2  
    set I 1.0;         #second moment of area I=1ft^4       
    set M 3.0;         #mas/length = 4 kip sec^2/ft^2           
    set coordTransf "Linear";  # Linear, PDelta, Corotational
    set massType "-lMass";  # -lMass, -cMass


    set nPts 3;      # numGauss Points

    # an elastic material
    uniaxialMaterial Elastic 1 $E

    # an elastic section
    section Elastic 1 $E $A $I
        
    # a fiber section with A=3 and I = 1 (b=1.5, d=2) 2d bending about y-y axis
 #   set b 1.5; set d 2.0;
    set y 2.0; set z 1.5;
    set numFiberY 2000;  # note we only need so many to get the required accuracy on eigenvalue 1e-7!
    set numFiberZ 1;
    section Fiber 2 {
#	patch rect 1 $numFiberY $numFiberZ 0.0 0.0 $z $y
	patch quad 1 $numFiberY $numFiberZ [expr -$y/2.0] [expr -$z/2.0] [expr $y/2.0] [expr -$z/2.0] [expr $y/2.0] [expr $z/2.0] [expr -$y/2.0] [expr $z/2.0]
    }
    

    # add the nodes         
    #  - floor at a time    
    set nodeTag 1
    set yLoc 0.
    for {set j 0} {$j <= $numFloor} {incr j 1} {
	set xLoc 0.
	for {set i 0} {$i <=$numBay} {incr i 1} {
	    node $nodeTag $xLoc $yLoc
	    set xLoc [expr $xLoc + $bayWidth]
	    incr nodeTag 1
	}
	set yLoc [expr $yLoc + $storyHeight]
    }
    
    # fix base nodes        
    for {set i 1} {$i <= [expr $numBay+1]} {incr i 1} {
	fix $i 1 1 1
    }
    
    # add column element    
    geomTransf $coordTransf 1
    set eleTag 1
    for {set i 0} {$i <=$numBay} {incr i 1} {
	set end1 [expr $i+1]
	set end2 [expr $end1 + $numBay +1]
	for {set j 0} {$j<$numFloor} {incr j 1} {

	    if {$eleType == "elasticBeam"} {
		element elasticBeamColumn $eleTag $end1 $end2 $A $E $I 1 -mass $M $massType
	    } elseif {$eleType == "forceBeamElasticSection"} {
		element forceBeamColumn $eleTag $end1 $end2 $nPts 1 1 -mass $M
	    } elseif {$eleType == "dispBeamElasticSection"} {
		element dispBeamColumn $eleTag $end1 $end2 $nPts 1 1 -mass $M $massType
	    } elseif {$eleType == "forceBeamFiberSectionElasticMaterial"} {
		element forceBeamColumn $eleTag $end1 $end2 $nPts 2 1 -mass $M
	    } elseif {$eleType == "dispBeamFiberSectionElasticMaterial"} {
		element dispBeamColumn $eleTag $end1 $end2 $nPts 2 1 -mass $M $massType
	    } else {
		puts "BARF"
	    }
	    set end1 $end2
	    set end2 [expr $end1 + $numBay +1]
	    incr eleTag 1
	}
    }
    

    # add beam elements     
    for {set j 1} {$j<=$numFloor} {incr j 1} {
	set end1 [expr ($numBay+1)*$j+1]
	set end2 [expr $end1 + 1]
	for {set i 0} {$i <$numBay} {incr i 1} {
	    if {$eleType == "elasticBeam"} {
		element elasticBeamColumn $eleTag $end1 $end2 $A $E $I 1 -mass $M $massType
	    } elseif {$eleType == "forceBeamElasticSection"} {
		element forceBeamColumn $eleTag $end1 $end2 $nPts 1 1 -mass $M
	    } elseif {$eleType == "dispBeamElasticSection"} {
		element dispBeamColumn $eleTag $end1 $end2 $nPts 1 1 -mass $M $massType
	    } elseif {$eleType == "forceBeamFiberSectionElasticMaterial"} {
		element forceBeamColumn $eleTag $end1 $end2 $nPts 2 1 -mass $M
	    } elseif {$eleType == "dispBeamFiberSectionElasticMaterial"} {
		element dispBeamColumn $eleTag $end1 $end2 $nPts 2 1 -mass $M $massType
	    } else {
		puts "BARF"
	    }
#	    element elasticBeamColumn $eleTag $end1 $end2 $A $E $I 1 -mass $M
	    set end1 $end2
	    set end2 [expr $end1 + 1]
	    incr eleTag 1
	}
    }
    
    # calculate eigenvalues
    set numEigen 3
    set eigenValues [eigen $solverType $numEigen]
    set PI [expr 2*asin(1.0)]
    
    # determine PASS/FAILURE of test
    set testOK 0

    # print table of camparsion
    #                         Bathe & Wilson               Peterson                    SAP2000                  SeismoStruct

    set comparisonResults {{0.589541 5.52695 16.5878} {0.589541 5.52696 16.5879} {0.589541 5.52696 16.5879} {0.58955 5.527 16.588}}
    puts "\n\nEigenvalue Comparisons for solverType: $solverType"
    set tolerances {9.99e-7 9.99e-6 9.99e-5}
    set formatString {%15s%15s%15s%15s%15s}
    puts [format $formatString OpenSees Bathe&Wilson Peterson SAP2000 SeismoStruct]
    set formatString {%15.5f%15.4f%15.4f%15.4f%15.3f}
    for {set i 0} {$i<$numEigen} {incr i 1} {
	set lambda [lindex $eigenValues $i]
	puts [format $formatString $lambda  [lindex [lindex $comparisonResults 0] $i] [lindex [lindex $comparisonResults 1] $i] [lindex [lindex $comparisonResults 2] $i] [lindex [lindex $comparisonResults 3] $i]]
	set resultOther [lindex [lindex $comparisonResults 2] $i]
	set tol [lindex $tolerances $i]
	if {[expr abs($lambda-$resultOther)] > $tol} {
	    set testOK -1;
	    puts "failed-> [expr abs($lambda-$resultOther)] $tol"
	}
    }
    
    set results [open results.out a+]
    if {$testOK == 0} {
	puts "PASSED Verification Test EigenFrame.Extra.tcl $eleType  \n\n"
	puts $results "PASSED : EigenFrame.Extra.tcl solverType: $solverType"
    } else {
	puts "FAILED Verification Test EigenFrame.Extra.tcl $eleType  \n\n"
	puts $results "FAILED : EigenFrame.Extra.tcl solverType: $solverType"
    }
    close $results
}

