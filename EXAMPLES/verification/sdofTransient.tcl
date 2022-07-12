# Linear Elastic SINGLE DOF Model Transient Analysis

#REFERENCES: 
# 1) Chopra, A.K. "Dynamics of Structures: Theory and Applications"
# Prentice Hall, 1995.
#   - Sections 3.1, Section 3.2 and Section 6.4

puts "sdofTransient.tcl: Verification of Elastic SDOF systems (Chopra)"

#
# global variables
#

set PI [expr 2.0*asin(1.0)]
set g 386.4
set testOK 0;    # variable used to keep track of SUCCESS or FAILURE
set tol 1.0e-2

# procedure to build a linear model
#   input args: K - desired stiffness
#               periodStruct - desired structure period (used to compute mass)
#               dampRatio (zeta) - desired damping ratio

proc buildModel {K periodStruct dampRatio} {

    global PI

    set wn [expr 2.0 * $PI / $periodStruct]
    set m [expr $K/($wn * $wn)]
    
    wipe
    model basic -ndm 1 -ndf 1
    
    node  1  0.
    node  2  0. -mass $m
    
    uniaxialMaterial Elastic 1 $K
    element zeroLength 1 1 2 -mat 1 -dir 1
    
    fix 1 1

    # add damping using rayleigh damping on the mass term
    set a0 [expr 2.0*$wn*$dampRatio]
    rayleigh $a0 0. 0. 0.
}

#
# procedure to build a linear transient analysis
#    input args: integrator command

proc buildLinearAnalysis {integratorCommand} {
    # do analysis
    constraints Plain
    numberer Plain
    algorithm Linear
    eval $integratorCommand; # integrator Newmark 0.5 [expr 1.0/6.0]
    system ProfileSPD
    analysis Transient
}

# Section 3.1 - Harmonic Vibrartion of Undamped Elastic SDOF System
puts "   - Undamped System Harmonic Exciatation (Section 3.1)"

# harmonic force propertires
set P 2.0
set periodForce 5.0
set tFinal [expr 2.251*$periodForce]

# model properties
set periodStruct 0.8
set K 2.0

# derived quantities
set w [expr 2.0 * $PI / $periodForce]
set wn [expr 2.0 * $PI / $periodStruct]

# build the model
buildModel $K $periodStruct 0.0

# add load pattern
timeSeries Trig 1 0.0 [expr 100.0*$periodForce] $periodForce -factor $P
pattern Plain 1 1 {
    load 2 1.0 
}

# build analysis
buildLinearAnalysis "integrator Newmark 0.5 [expr 1.0/6.0]"

# perform analysis, checking at every step
set dt [expr $periodStruct/1.0e4]; # something small for accuracy
set tCurrent 0.
puts "\n  for 1000 time steps computing solution and checking against exact solution"
set count 0

while {$tCurrent < $tFinal} {
    analyze 1 $dt
    set tCurrent [getTime]
    set uOpenSees [nodeDisp 2 0]
    set uExact [expr $P/$K * 1.0/(1 - ($w*$w)/($wn*$wn)) * (sin($w*$tCurrent) - ($w/$wn)*sin($wn*$tCurrent))]

    if {[expr abs($uExact-$uOpenSees)] > $tol} {
	set testOK -1;
	puts "failed  undamped harmonic> [expr abs($uExact-$uOpenSees)]> $tol at time $tCurrent"
	set tCurrent $tFinal
    }
}

set formatString {%20s%15.5f%10s%15.5f}
puts "\n  example results for last step at $tCurrent (sec):"
puts [format $formatString OpenSees: $uOpenSees Exact: $uExact]

if {[expr abs($uExact-$uOpenSees)] > $tol} {
    set testOK -1;
    puts "failed  undamped harmonic> [expr abs($uExact-$uOpenSees)] $tol"
}

# Section 3.2 - Harmonic Vibrartion of Damped Elastic SDOF System
puts "\n\n   - Damped System Harmonic Excitation (Section 3.2)"

set dampRatio 0.05

# build the model
buildModel $K $periodStruct $dampRatio

# add load pattern
timeSeries Trig 1 0.0 [expr 100.0*$periodForce] $periodForce -factor $P
pattern Plain 1 1 {
    load 2 1.0 
}

# build analysis
buildLinearAnalysis "integrator Newmark 0.5 [expr 1.0/6.0]"

# some variables needed in exact computation
set wd [expr $wn*sqrt(1-$dampRatio*$dampRatio)]
set wwn2 [expr ($w*$w)/($wn*$wn)]
set det [expr (1.0-$wwn2)*(1-$wwn2) + 4.0 * $dampRatio*$dampRatio*$wwn2]
set ust [expr $P/$K]

set C [expr $ust/$det * (1.-$wwn2)]
set D [expr $ust/$det * (-2.*$dampRatio*$w/$wn)]
set A -$D;
set B [expr $ust/$det * (1.0/$wd) * ((-2. * $dampRatio * $w/$wn) - $w * (1.0 - $wwn2))]


set t 0.
while {$t < $tFinal} {
    analyze 1 $dt
    set t [getTime]
    set uOpenSees [nodeDisp 2 0]
    set uExact [expr exp(-$dampRatio*$wn*$t)*($A * cos($wd * $t) + $B*sin($wd*$t)) + $C*sin($w*$t) + $D*cos($w*$t)]
    if {[expr abs($uExact-$uOpenSees)] > $tol} {
	set testOK -1;
	puts "failed  damped harmonic> [expr abs($uExact-$uOpenSees)]> $tol at time $t"
	set t $tFinal
    }
}

set formatString {%20s%15.5f%10s%15.5f}

puts "\nDisplacement Comparison at $tCurrent (sec):"
puts [format $formatString OpenSees: $uOpenSees Exact: $uExact]

if {[expr abs($uExact-$uOpenSees)] > $tol} {
    set testOK -1;
    puts "failed  undamped harmonic> [expr abs($uExact-$uOpenSees)] $tol"
}

# Section 6.4 - Earthquake Response of Linear System

puts "\n\n   - Earthquake Response (Section 6.4)\n"


set tol 3.0e-2; 
set results {2.67 5.97 7.47 9.91 7.47 5.37}

# read earthquake record, setting dt and nPts variables with data in the file elCentro.at2
source ReadRecord.tcl
ReadRecord elCentro.at2 elCentro.dat dt nPts

# print table header
set formatString {%15s%15s%15s%15s}
puts [format $formatString Period dampRatio OpenSees Exact]

# perform analysis for bunch of periods and damping ratio's
set counter 0
foreach {period dampRatio} {0.5 0.02 1.0 0.02 2.0 0.02 2.0 0.0 2.0 0.02 2.0 0.05} {

    # build the model
    buildModel $K $period $dampRatio

    # add load pattern
    timeSeries Path 1 -filePath elCentro.dat -dt $dt -factor $g
    pattern UniformExcitation 1 1 -accel 1

    # build analysis
    buildLinearAnalysis "integrator Newmark 0.5 [expr 1.0/6.0]"

    set maxU 0.0; 
    for {set i 0} {$i < $nPts} {incr i 1} {
	analyze 1 $dt
	set u [expr abs([nodeDisp 2 1])]
	if {$u > $maxU} {
	    set maxU $u
	}
    }
    set formatString {%15.2f%15.2f%15.2f%15.2f}
    set uExact [lindex $results $counter]
    puts [format $formatString $period $dampRatio $maxU $uExact]
    if {[expr abs($maxU-$uExact)] > $tol} {
	set testOK -1;
	puts "failed  earthquake record period: $period dampRatio: $dampRatio: $maxU $uExact [expr abs($uExact-$maxU)] $tol"
    }
    incr counter
}

set results [open results.out a+]
if {$testOK == 0} {
    puts "\nPASSED Verification Test sdofTransient.tcl \n\n"
    puts $results "PASSED : sdofTransient.tcl"
} else {
    puts "\nFAILED Verification Test sdofTransient.tcl \n\n"
    puts $results "FAILED : sdofTransient.tcl"
}
close $results



