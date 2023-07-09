
# Two dimensional Frame: Eigenvalue & Static Loads


# REFERENCES:
# used in verification by SAP2000:
# SAP2000 Integrated Finite Element Analysis and Design of Structures, Verification Manual, 
# Computers and Structures, 1997. Example 1.
# and seismo-struct (Example 10)
# SeismoStruct, Verification Report For Version 6, 2012. Example 11.



# set some properties

puts "PortalFrame2d.tcl: Verification 2d Elastic Frame"
puts "  - eigenvalue and static pushover analysis"

wipe

model Basic -ndm 2

# properties  

#    units kip, ft

set numBay 2
set numFloor 7

set bayWidth 360.0;
set storyHeights {162.0 162.0 156.0 156.0 156.0 156.0 156.0}

set E 29500.0
set massX 0.49
set M 0.
set coordTransf "Linear";  # Linear, PDelta, Corotational
set massType "-lMass";  # -lMass, -cMass


set beams   {W24X160 W24X160 W24X130 W24X130 W24X110 W24X110 W24X110}
set eColumn {{W14X246 W14X246 W14X246 W14X211 W14X211 W14X176 W14X176}}
set iColumn {{W14X287 W14X287 W14X287 W14X246 W14X246 W14X211 W14X211}}
set columns [concat $eColumn $iColumn $eColumn]

array set WSection {
    W14X176  	"51.7 2150."
    W14X211  	"62.1 2670."
    W14X246  	"72.3 3230."
    W14X287  	"84.4 3910."
    W24X110  	"32.5 3330."
    W24X130  	"38.3 4020."
    W24X160  	"47.1 5120."
}

# procedure to read 

proc ElasticBeamColumn {eleTag iNode jNode sectType E transfTag M massType} {
    global WSection
    global in
    set found 0

    set sectType [string map {x X} $sectType]
    set sectType [string map {w W} $sectType]

    foreach {section prop} [array get WSection $sectType] {
	set propList [split $prop]
	set A [lindex $propList 0]
	set I [lindex $propList 1]
	element elasticBeamColumn $eleTag $iNode $jNode $A $E $I $transfTag -mass $M $massType
	set found 1
    }
    
    if {$found == 0} {
	puts "FiberSteelWSection2d sectType: $sectType not found for sectTag: $sectTag"
    }
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
    if {$j < $numFloor} {
	set storyHeight [lindex $storyHeights $j]
    }
    set yLoc [expr $yLoc + $storyHeight]
}

# fix first floor
fix 1 1 1 1
fix 2 1 1 1
fix 3 1 1 1

#rigid floor constraint & masses
set nodeTagR 5
set nodeTag  4
for {set j 1} {$j <= $numFloor} {incr j 1} {
    for {set i 0} {$i <=$numBay} {incr i 1} {
	if {$nodeTag != $nodeTagR} {
	    equalDOF $nodeTagR $nodeTag 1
	} else {
	    mass $nodeTagR $massX 1.0e-10 1.0e-10
	}
	incr nodeTag 1
    }
    incr nodeTagR [expr $numBay+1]
}

# add the columns
# add column element    
geomTransf $coordTransf 1
set eleTag 1
for {set j 0} {$j <= $numBay} {incr j 1} {
    set end1 [expr $j+1]
    set end2 [expr $end1 + $numBay +1]
    set thisColumn [lindex $columns $j]
    for {set i 0} {$i<$numFloor} {incr i 1} {
	set secType [lindex $thisColumn $i]
	ElasticBeamColumn $eleTag $end1 $end2 $secType $E 1 $M &massType
	set end1 $end2
	set end2 [expr $end1 + $numBay +1]
	incr eleTag 1
    }
}

# add beam elements     
for {set j 1} {$j<=$numFloor} {incr j 1} {
    set end1 [expr ($numBay+1)*$j+1]
    set end2 [expr $end1 + 1]
    set secType [lindex $beams [expr $j-1]]
    for {set i 0} {$i <$numBay} {incr i 1} {
	ElasticBeamColumn $eleTag $end1 $end2 $secType $E 1 $M &massType
        set end1 $end2
	set end2 [expr $end1 + 1]
        incr eleTag 1
    }
}

# calculate eigenvalues & print results     
set numEigen 7
set eigenValues [eigen $numEigen]
set PI [expr 2*asin(1.0)]

#
# apply loads for static analysis & perform analysis
#

timeSeries Linear 1
pattern Plain 1 1 {
   load  22 20.0 0. 0.
   load  19 15.0 0. 0.
   load  16 12.5 0. 0.
   load  13 10.0 0. 0.
   load  10  7.5 0. 0.
   load   7  5.0 0. 0.
   load   4  2.5 0. 0. 
}
integrator LoadControl 1.0
algorithm Linear
analysis Static
analyze 1

# determine PASS/FAILURE of test
set ok 0

#
# print pretty output of comparisons
#

#               SAP2000   SeismoStruct
set comparisonResults {{1.2732 0.4313 0.2420 0.1602 0.1190 0.0951 0.0795} {1.2732 0.4313 0.2420 0.1602 0.1190 0.0951 0.0795}}
puts "\n\nPeriod Comparisons:"
set formatString {%10s%15s%15s%15s}
puts [format $formatString Period OpenSees SAP2000 SeismoStruct]
set formatString {%10s%15.5f%15.4f%15.4f}
for {set i 0} {$i<$numEigen} {incr i 1} {
    set lambda [lindex $eigenValues $i]
    set period [expr 2*$PI/sqrt($lambda)];
    puts [format $formatString [expr $i+1] $period [lindex [lindex $comparisonResults 0] $i] [lindex [lindex $comparisonResults 1] $i]]
    set resultOther [lindex [lindex $comparisonResults 0] $i]
    if {[expr abs($period-$resultOther)] > 9.99e-5} {
	set ok -1;
    }
}


# print table of camparsion
#       Parameter          SAP2000   SeismoStruct
set comparisonResults {{"Disp Top" "Axial Force Bottom Left" "Moment Bottom Left"} {1.45076 69.99 2324.68} {1.451 70.01 2324.71}}
set tolerances {9.99e-6 9.99e-3 9.99e-3}

puts "\n\nSatic Analysis Result Comparisons:"
set formatString {%30s%15s%15s%15s}
puts [format $formatString Parameter OpenSees SAP2000 SeismoStruct]
set formatString {%30s%15.3f%15.2f%15.2f}
for {set i 0} {$i<3} {incr i 1} {
    if {$i == 0} {
	set result [nodeDisp 22 1]
    } elseif {$i == 1} {
	set result [expr abs([lindex [eleResponse 1 forces] 1])]
    } else {
	set result [lindex [eleResponse 1 forces] 2]	
    }
    puts [format $formatString [lindex [lindex $comparisonResults 0] $i] $result [lindex [lindex $comparisonResults 1] $i] [lindex [lindex $comparisonResults 2] $i]]
    set resultOther [lindex [lindex $comparisonResults 1] $i]
    set tol [lindex $tolerances $i]
    if {[expr abs($result-$resultOther)] > $tol} {
	set ok -1;
	puts "failed-> $i [expr abs($result-$resultOther)] $tol"
    }
}

set results [open results.out a+]
if {$ok == 0} {
    puts "PASSED Verification Test PortalFrame2d.tcl \n\n"
    puts $results "PASSED : PortalFrame2d.tcl"
} else {
    puts "FAILED Verification Test PortalFrame2d.tcl \n\n"
    puts $results "FAILED : PortalFrame2d.tcl"
}
close $results

