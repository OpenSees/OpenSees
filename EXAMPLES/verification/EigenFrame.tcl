# Bathe & Wilson eigenvalue problem
#   Results presented by Bathe & wilson in 1972
#   and again by Peterson in 1981

#REFERENCES: 
# 1) Bathe, K.J and Wilson, E.L.Large Eigenvalue Problems in Synamic Analysis, ASCE,
# Journal of Eng. Mech.. Division, 98(6), 1471-1485, 
# 2) Peterson, F.E. EASE2, Elastic Analysis for Structural Engineering, Example Problem Manual, 
# Engineering Analysis Corporation, Berkeley, CA 1981.

# used in verification by SAP2000:
# SAP2000 Integrated Finite Element Analysis and Design of Structures, Verification Manual, 
# Computers and Structures, 1997. Example 2.
# and seismo-struct (Example 10)
# SeismoStruct, Verification Report For Version 6, 2012. Example 10.


puts "EigenFrame.tcl: Verification 2d Bathe & Wilson original Elastic Frame"
puts "  - eigenvalue "

wipe

model Basic -ndm 2

#    units kip, ft                                                                                                                              

# properties  
set bayWidth 20.0;
set storyHeight 10.0;

set numBay 10
set numFloor 9
set A 3.0;         #area = 3ft^2    
set E 432000.0;   #youngs mod = 432000 k/ft^2  
set I 1.0;         #second moment of area I=1ft^4       
set M 3.0;      #mas/length = 4 kip sec^2/ft^2       
set coordTransf "Linear";  # Linear, PDelta, Corotational
set massType "-lMass";  # -lMass, -cMass


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
	element elasticBeamColumn $eleTag $end1 $end2 $A $E $I 1 -mass $M $massType
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
        element elasticBeamColumn $eleTag $end1 $end2 $A $E $I 1 -mass $M $massType
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
puts "\n\nEigenvalue Comparisons:"
set tolerances {9.99e-6 9.99e-6 9.99e-5}; # tolerances prescribed by documented precision
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
    puts "PASSED Verification Test EigenFrame.tcl \n\n"
    puts $results "PASSED : EigenFrame.tcl"
} else {
    puts "FAILED Verification Test EigenFrame.tcl \n\n"
    puts $results "FAILED : EigenFrame.tcl"
}
close $results
