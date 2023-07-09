
# Newmark Integrators - Linear & Nonlinear Examples

#REFERENCES: 
# 1) Chopra, A.K. "Dynamics of Structures: Theory and Applications"
# Prentice Hall, 4th Edition, 2012.
#   - Sections 5:
#        Linear: Examples 5.3 and 5.4
#        Nonlinear: Examples 5.5 and 5.6


puts "Newmark.tcl: Verification of Newmark Integrators (Chopra)"

# global variables
set PI [expr 2.0*asin(1.0)]
set testOK 0;    # variable used to keep track of SUCCESS or FAILURE
set tol 1.0e-3

# model properties
set m 0.2533
set k 10.0
set dampRatio 0.05

set integratorCmds {"Newmark Average Acceleration" "integrator Newmark 0.5 0.25" "Newmark Linear Acceleration" "integrator Newmark 0.5 [expr 1.0/6.0]"}

# procedure to build a linear model
#   input args: K - desired stiffness
#               periodStruct - desired structure period (used to compute mass)
#               dampRatio (zeta) - desired damping ratio
proc buildModel {k m dampRatio yieldDisp} {

    global PI

    set wn [expr sqrt($k/$m)]
    
    wipe
    model basic -ndm 1 -ndf 1
    
    node  1  0.
    node  2  0. -mass $m

    if {$yieldDisp == 0.} {
	uniaxialMaterial Elastic 1 $k
    } else {
	uniaxialMaterial ElasticPP 1 $k $yieldDisp
    }
    element zeroLength 1 1 2 -mat 1 -dir 1
    
    fix 1 1

    # add damping using rayleigh damping on the mass term
    set a0 [expr 2.0*$wn*$dampRatio]
    rayleigh $a0 0. 0. 0.

}

#
# procedure to build a transient analysis
#    input args: integrator command
#                algoType: Linear do Linear analysis, else do Nonlinear Newton
proc buildAnalysis {integratorCommand algoCmd} {

    constraints Plain
    numberer Plain
    eval $integratorCommand; 
    test NormDispIncr 1.0e-4 6 0
    eval $algoCmd
    system ProfileSPD
    analysis Transient
}

# Section 5.4 - Newmark Linear System
puts "Linear System"

set resultsD {
    {0.0437 0.2326 0.6121 1.0825 1.4309 1.4231 0.9622 0.1908 -0.6044 -1.1442}
    {0.0300 0.2193 0.6166 1.1130 1.4782 1.4625 0.9514 0.1273 -0.6954 -1.2208}}
set resultsV {
    {0.8733 2.9057 4.6833 4.4260 2.2421 -2.3996 -6.8192 -8.6092 -7.2932 -3.5026}
    {0.8995 2.9819 4.7716 4.7419 2.1802 -2.6911 -7.1468 -8.7758 -7.1539 -3.0508}}
set resultsA {
    {17.4666 23.1801 12.3719 -11.5175 -38.1611 -54.6722 -33.6997 -2.1211 28.4423 47.3701}
    {17.9904 23.6566 12.1372 -12.7305 -39.9425 -56.0447 -33.0689  0.4892 31.9491 50.1114}} 


set count 0
foreach {integratorName integratorCmd} $integratorCmds {

    puts "\n - $integratorName"
    set formatString {%20s%20s%20s}    
    puts [format $formatString Displacement Velocity Acceleration]
    set formatString {%10s%10s%10s%10s%10s%10s}    
    puts [format $formatString OpenSees Hand OpenSees Hand OpenSees Hand]
    set formatString {%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f}
    set resultD [lindex $resultsD $count]
    set resultV [lindex $resultsV $count]
    set resultA [lindex $resultsA $count]

    # build the model
    buildModel $k $m $dampRatio 0.

    # add load pattern
    timeSeries Trig 1 0.0 0.6 1.2 -factor 10.0
    #    timeSeries Trig 1 -dt 0.1 -path {0.0 5.0}
    pattern Plain 1 1 {
	load 2 1.0 
    }

    # build analysis
    buildAnalysis $integratorCmd "algorithm Linear"

    # perform analysis, checking at every step
    for {set i 0} {$i< 10} {incr i 1} {
	analyze 1 0.1
	set tCurrent [getTime]
	set uOpenSees [nodeDisp 2 0]
	set uComputed [lindex $resultD $i]
	set vOpenSees [nodeVel 2 0]
	set vComputed [lindex $resultV $i]
	set aOpenSees [nodeAccel 2 0]
	set aComputed [lindex $resultA $i]
	if {[expr abs($uComputed-$uOpenSees)] > $tol} {
	    set testOK -1;
	    puts [format $formatString $uOpenSees $uComputed]
	    puts "failed  abs($uOpenSees - $uComputed) = [expr abs($uComputed-$uOpenSees)]> $tol"
	    set i 11
	} else {
#
	    puts [format $formatString $uOpenSees $uComputed $vOpenSees $vComputed $aOpenSees $aComputed]
	}
    }
    incr count
}

# Section 5.7 - Newmark Nonlinear System
puts "\n\nNonlinear System - Newton Average Acceleration With Differing Nonlinear Algorithms"

set integratorCmd "integrator Newmark 0.5 0.25" 

set algorithmCmds {
    "Newton" "algorithm Newton" 
    "Modified Newton"     "algorithm ModifiedNewton"}

set resultsD {
    {0.0437 0.2326 0.6121 1.1143 1.6214 1.9891 2.0951 1.9240 1.5602 1.415}
    {0.0437 0.2326 0.6121 1.1143 1.6214 1.9891 2.0951 1.9240 1.5602 1.414}}
set resultsV {
    {0.8733 2.9057 4.6833 5.3624 4.7792 2.5742 -0.4534 -2.960 -4.3075 -4.0668}
    {0.8733 2.9057 4.6833 5.3623 4.7791 2.5741 -0.4534 -2.960 -4.3076 -4.0668}}
set resultsA {
    {17.4666 23.1801 12.3719 1.2103 -12.8735 -31.2270 -29.3242 -20.9876 -5.7830 10.5962}
    {17.4666 23.1801 12.3719 1.2095 -12.8734 -31.2270 -29.3242 -20.9879 -5.7824 10.5969}}


set count 0
foreach {algoName algoCmd} $algorithmCmds {

    puts "\n - $algoCmd"
    set formatString {%20s%20s%20s}    
    puts [format $formatString Displacement Velocity Acceleration]
    set formatString {%10s%10s%10s%10s%10s%10s}    
    puts [format $formatString OpenSees Hand OpenSees Hand OpenSees Hand]
    set formatString {%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f}

    set resultD [lindex $resultsD $count]
    set resultV [lindex $resultsV $count]
    set resultA [lindex $resultsA $count]

    # build the model
    buildModel $k $m $dampRatio 0.75

    # add load pattern
    timeSeries Trig 1 0.0 0.6 1.2 -factor 10.0
    #    timeSeries Trig 1 -dt 0.1 -path {0.0 5.0}
    pattern Plain 1 1 {
	load 2 1.0 
    }

    # build analysis
    buildAnalysis $integratorCmd $algoCmd

    # perform analysis, checking at every step
    for {set i 0} {$i< 9} {incr i 1} {
	analyze 1 0.1
	set tCurrent [getTime]
	set uOpenSees [nodeDisp 2 0]
	set uComputed [lindex $resultD $i]
	set vOpenSees [nodeVel 2 0]
	set vComputed [lindex $resultV $i]
	set aOpenSees [nodeAccel 2 0]
	set aComputed [lindex $resultA $i]
	if {[expr abs($uComputed-$uOpenSees)] > $tol} {
	    set testOK -1;
	    puts [format $formatString $uOpenSees $uComputed]
	    puts "failed  abs($uOpenSees - $uComputed) = [expr abs($uComputed-$uOpenSees)]> $tol"
	    set i 11
	} else {
	    puts [ format $formatString $uOpenSees $uComputed $vOpenSees $vComputed $aOpenSees $aComputed]
	}
    }
    incr count
}

set results [open results.out a+]
if {$testOK == 0} {
    puts "\nPASSED Verification Test NewmarkIntegratir.tcl \n\n"
    puts $results "PASSED : NewmarkIntegrator.tcl"
} else {
    puts "\nFAILED Verification Test NewmarkIntegrator.tcl \n\n"
    puts $results "FAILED : NewmarkIntegrator.tcl"
}
close $results



