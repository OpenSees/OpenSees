
# Pinched Shell Cylindrical Problem
#Lindberg, G. M. M., D. Olson, and G. R. Cowper, “New Developments in the Finite Element
#AnalysisofShells,”QuarterlyBulletinoftheDivisionofMechanicalEngineeringandtheNational
#AeronauticalEstablishment,NationalResearchCouncilofCanada,vol. 4,1969.

# R = Radius, t= thickness, v = 0.3
# for R/t = 100, v = 0.3, L/R = 2 displacement under load = 164.24 * P / (E * t)

puts "PinchedCantiliver: Validation of Shell Elements with Elastic Sections"

set P 1
set R 300
set L 600
set E 3e6
set thickness [expr $R/100.]

set uExact [expr -164.24*$P/($E*$thickness)]

set formatString {%20s%15.5e}
puts "\n  Displacement Under Applied Load:\n"

set formatString {%20s%10s%15s%15s%15s}
puts [format $formatString "Element Type" "   mesh  " "OpenSees" "Exact" "%Error"]

foreach shellType {ShellMITC4 ShellDKGQ ShellNLDKGQ} {
    foreach numEle {4 16 32} {

	wipe

	# ----------------------------
	# Start of model generation
	# ----------------------------
	model basic -ndm 3 -ndf 6
	
	set radius $R
	set length [expr $L/2.]

	
	set E 3.0e6
	set v 0.3
	set PI 3.14159
	
	nDMaterial ElasticIsotropic 1 $E $v
	nDMaterial PlateFiber 2 1
	section PlateFiber 1 2 $thickness
	#section ElasticMembranePlateSection  1   $E $v $thickness 0.
	set eleArgs "1"
	
	set nR $numEle
	set nY $numEle

	set tipNode [expr ($nR+1)*($nY+1)]
	
	#create nodes
	for {set i 0; set nodeTag 1} {$i<= $nR} {incr i 1} {
	    set theta [expr $i*$PI/(2.0*$nR)]
	    set xLoc [expr 300*cos($theta)]
	    set zLoc [expr 300*sin($theta)]
	    for {set j 0;} {$j <= $nY} {incr j 1; incr nodeTag; } {
		set yLoc [expr $j*$length/(1.0*$nY)]
		node $nodeTag $xLoc $yLoc $zLoc
		
	    }
	}
	
	#create elements
	
	for {set i 0; set eleTag 1} {$i< $nR} {incr i 1} {
	    set iNode [expr $i*($nY+1)+1]; 
	    set jNode [expr $iNode +1];
	    set lNode [expr $iNode+($nY+1)]
	    set kNode [expr $lNode+1]
	    for {set j 0;} {$j < $nY} {incr j 1; incr eleTag; incr iNode; incr jNode; incr kNode; incr lNode} {
		element $shellType $eleTag $iNode $jNode $kNode $lNode 1
	    }
	}
	# define the boundary conditions
	fixX $radius  0 0 1 1 1 0 -tol 1.0e-2
	fixZ $radius  1 0 0 0 1 1 -tol 1.0e-2
	fixY      0   1 0 1 0 0 0 -tol 1.0e-2
	fixY $length  0 1 0 1 0 1 -tol 1.0e-2
	
	#define loads
	timeSeries Linear 1
	pattern Plain 1 1 {
	    load $tipNode 0 0 [expr -1./4.0] 0 0 0; 
	}
	
	integrator LoadControl  1.0 
	test EnergyIncr     1.0e-10    20       0
	algorithm Newton
	numberer RCM
	constraints Plain 
	system Umfpack
	analysis Static 
	
	analyze 1
	set res [nodeDisp $tipNode 3]
	set err [expr abs(100*($uExact-$res)/$uExact)]
	set formatString {%20s%5d%3s%2d%15.5e%15.5e%15.2f}
	puts [format $formatString $shellType $numEle " x " $numEle $res $uExact $err]
    }
}

set tol 5.0; 
if {[expr abs(100*($uExact-$res)/$uExact)] > $tol} {
    set testOK 1
} else {
    set testOK 0
}

set results [open results.out a+]
if {$testOK == 0} {
    puts "\nPASSED Validation Test PinchedCylinder.tcl \n\n"
    puts $results "PASSED : PinchedCylinder.tcl"
} else {
    puts "\nFAILED Validation Test PinchedCylinder.tcl \n\n"
    puts $results "FAILED : PinchedCylinder.tcl"
}

close $results