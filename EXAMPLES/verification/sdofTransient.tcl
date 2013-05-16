# SINGLE  DOF

#REFERENCES: 
# 1) Chopra, A.K. "Dynamics of Structures: Theory and Applications"
# Prentice Hall, 1995.
#   - Sections 3.1, Section 3.2 and Section 

puts "sdofTransient.tcl: Verification of Elastic SDOF systems (Chopra)"

set testOK 0;

# Section 3.1 - Harmonic Vibrartion of Undamped Elastic SDOF System
puts "   - Undamped System Harmonic Exciatation (Section 3.1)"

# harmonic force propertires
set P 2.0
set periodForce 5.0

# sdof system
set periodStruct 0.8
set m 2.0;

set tFinal [expr 2.251*$periodForce]

# derived quantaties
set PI [expr 2.0*asin(1.0)]
set w [expr 2.0 * $PI / $periodForce]
set wn [expr 2.0 * $PI / $periodStruct]

set K [expr $wn * $wn * $m]

wipe
model basic -ndm 1 -ndf 1

node  1  0.
node  2  0. -mass $m

uniaxialMaterial Elastic 1 $K
element zeroLength 1 1 2 -mat 1 -dir 1

fix 1 1

timeSeries Trig 1 0.0 [expr 100.0*$periodForce] $periodForce -factor $P
pattern Plain 1 1 {
    load 2 1.0 
}

constraints Plain
numberer Plain
algorithm Linear
integrator Newmark 0.5 [expr 1.0/6.0]
system ProfileSPD
analysis Transient

set dt [expr $periodStruct/1.0e4]; # something small for accuracy

set tCurrent 0.
while {$tCurrent < $tFinal} {
    analyze 1 $dt
    set tCurrent [getTime]
}

set uOpenSees [nodeDisp 2 0]

# exact
set tFinal [getTime]

set uExact [expr $P/$K * 1.0/(1 - ($w*$w)/($wn*$wn)) * (sin($w*$tFinal) - ($w/$wn)*sin($wn*$tFinal))]

set formatString {%20s%15.5f%10s%15.5f}

puts "\nDisplacement Comparison at $tFinal (sec):"
puts [format $formatString OpenSees: $uOpenSees Exact: $uExact]

set tol 1.0e-4;
if {[expr abs($uExact-$uOpenSees)] > $tol} {
    set testOK -1;
    puts "failed  undamped harmonic> [expr abs($uExact-$uOpenSees)] $tol"
}

# Section 3.2 - Harmonic Vibrartion of Damped Elastic SDOF System
puts "\n\n   - Damped System Harmonic Excitation (Section 3.2)"

set dampRatio 0.01

# derived quantaties
set c [expr 2.0*$m*$wn*$dampRatio]
set a0 [expr 2.0*$wn*$dampRatio]

reset

rayleigh $a0 0. 0. 0.

set tCurrent 0.
while {$tCurrent < $tFinal} {
    analyze 1 $dt
    set tCurrent [getTime]
}

set tFinal [getTime]

set uOpenSees [nodeDisp 2 0]

# exact
set wwn2 [expr ($w*$w)/($wn*$wn)]
set C [expr $P/$K * (1.-$wwn2)/((1.-$wwn2)*(1.-$wwn2) + 4.*$dampRatio*$dampRatio*$wwn2)]

set D [expr $P/$K * (-2.*$dampRatio*$w/$wn)/((1.-$wwn2)*(1.-$wwn2) + 4.*$dampRatio*$dampRatio*$wwn2)]
set B [expr (-$P/$K) * (($w/$wn) / (1.-$wwn2))]

set wd [expr $wn*sqrt(1-$dampRatio*$dampRatio)]

set uExact [expr exp(-$dampRatio*$wn*$tFinal)*$B*sin($wd*$tFinal) + $C*sin($w*$tFinal) + $D*cos($w*$tFinal)]

set formatString {%20s%15.5f%10s%15.5f}

puts "\nDisplacement Comparison at $tFinal (sec):"
puts [format $formatString OpenSees: $uOpenSees Exact: $uExact]

set tol 1.0e-4;
if {[expr abs($uExact-$uOpenSees)] > $tol} {
    set testOK -1;
    puts "failed  undamped harmonic> [expr abs($uExact-$uOpenSees)] $tol"
}

# Section 6.4 - Earthquake Response of Linear System

puts "\n\n   - Earthquake Response (Section 6.4)\n"

    # sdof system
set L 2.0
set A 10.0;
set m 3.0;
set g 386.4;

# derived quantaties
set PI [expr 2.0*asin(1.0)]

source ReadRecord.tcl
ReadRecord elCentro.at2 elCentro.dat dt nPts

set results {2.67 5.97 7.47 9.91 7.47 5.37}

# print table header
set formatString {%15s%15s%15s%15s}
puts [format $formatString Period dampRatio OpenSees Exact]

# perform analysis for bunch of periods and damping ratio's
set counter 0
foreach {period dampRatio} {0.5 0.02 1.0 0.02 2.0 0.02 2.0 0.0 2.0 0.02 2.0 0.05} {

    # model derived properties
    set wn [expr 2.0 * $PI / $period]
    set K [expr $wn * $wn * $m]
    set E [expr $L * $K / $A]
    set a0 [expr 2.0*$wn*$dampRatio]
    set a1 [expr 2.0*$dampRatio/$wn]

    wipe
    model basic -ndm 1 -ndf 1
    
    node  1       0.0
    node  2       $L -mass $m

    uniaxialMaterial Elastic 1 $E 
    element truss 1 1 2 $A 1 

    rayleigh $a0 0.0 0.0 0.0 
    
    fix 1 1

    timeSeries Path 1 -filePath elCentro.dat -dt $dt -factor $g
    pattern UniformExcitation 1 1 -accel 1
    
    constraints Plain
    numberer Plain
    algorithm Linear

    integrator Newmark 0.5 [expr 1.0/6.0]
    system ProfileSPD
    analysis Transient

    set maxU 0.0; 
    for {set i 0} {$i < $nPts} {incr i 1} {
	analyze 1 $dt
	set u [expr abs([nodeDisp 2 1])]
	if {$u > $maxU} {
	    set maxU $u
	}
    }
    set formatString {%15.2f%15.2f%15.2f%15.2f}
    puts [format $formatString $period $dampRatio $maxU [lindex $results $counter]]
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



