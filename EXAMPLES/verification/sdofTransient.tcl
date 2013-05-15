# Single DOF

#REFERENCES: 
# 1) Chopra, A.K. "Dynamics of Structures: Theory and Applications"
# Prentice Hall, 1995.
#   - Sections 3.1 and Scetin 3.2


set testOK 0;

# Section 3.1 - Harmonic Vibrartion of Undamped Elastic SDOF System

puts "sdofTransient.tcl: Verification of Harmonic Vibrartion of Elastic SDOF systems (Chopra)"
puts "   - Undamped System (Section 3.1)"

# harmonic force
set P 1.0
set periodForce 10.

# sdof system
set periodStruct 2.
set L 1.0
set A 1.0;
set m 1.0;

# derived quantaties
set PI [expr 2.0*asin(1.0)]
set w [expr 2.0 * $PI / $periodForce]
set wn [expr 2.0 * $PI / $periodStruct]

set K [expr $wn * $wn * $m]
set E [expr $L * $K / $A]

wipe
model basic -ndm 1 -ndf 1

node  1       0.0
node  2       $L -mass $m
uniaxialMaterial Elastic 1 $E
element truss 1 1 2 $A 1

fix 1 1

timeSeries Trig 1 0.0 [expr 100.0*$periodForce] $periodForce -factor $P
pattern Plain 1 1 {
    load 2 1.0 
}

constraints Plain
numberer Plain
algorithm Linear
integrator Newmark 0.5 0.25
system ProfileSPD
analysis Transient

set tFinal [expr 4*$periodStruct]
set dt [expr $periodStruct/1.0e4]; # something small for accuracy

set tCurrent 0.
while {$tCurrent < $tFinal} {
    analyze 1 $dt
    set tCurrent [getTime]
}

set uOpenSees [nodeDisp 2 0]

# exact
set uExact [expr $P/$K * 1.0/(1 - ($w*$w)/($wn*$wn)) * (sin($w*$tFinal - ($w/$wn)*sin($wn*$tFinal)))]


set formatString {%20s%15.6f%10s%15.6f}

puts "\nDisplacement Comparison at $tFinal:"
puts [format $formatString OpenSees: $uOpenSees Exact: $uExact]

set tol 1.0e-4;
if {[expr abs($uExact-$uOpenSees)] > $tol} {
    set testOK -1;
    puts "failed  undamped harmonic> [expr abs($uExact-$uOpenSees)] $tol"
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



