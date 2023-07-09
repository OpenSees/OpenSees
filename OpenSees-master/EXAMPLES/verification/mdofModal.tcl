
#
# file to compare modal & rayleigh for 2 dof problem
#  NOTE: to get eigenvalue to work for 2 dof problem 
#        we had to add a dummy dof (so problem really 3 dof, 
#        though that last mode should have small contribution)
#

# idea to test model with 2 different loadings
#  case 1: harmonic
#  case 2: earthquake (elCentro)
#

#
# some parameters
#

set zeta 0.05
set periodForce 1.
set P 100.
set g 386.1
set myTol 1.0e-8
#
# function to build a model
# 

proc buildModel {} {

    global PI
    set m1 2.
    set m2 4.
    set k1 5.
    set k2 3.
    
    wipe
    model Basic -ndm 1 -ndf 1
    node 1 0.
    node 2 0. -mass $m1
    node 3 0. -mass $m2
    node 4 0. -mass [expr $m2/1e8]

    uniaxialMaterial Elastic 1 $k1
    uniaxialMaterial Elastic 2 $k2
    uniaxialMaterial Elastic 3 [expr $k2*1e8]

    element zeroLength 1 1 2 -mat 1 -dir 1 -doRayleigh 1
    element zeroLength 2 2 3 -mat 1 -dir 1 -doRayleigh 1
    element zeroLength 3 3 4 -mat 1 -dir 1 
    
    fix 1 1
}


#
# trig motio
#

buildModel 

# add load pattern
timeSeries Trig 1 0.0 [expr 100.0*$periodForce] $periodForce -factor $P
pattern Plain 1 1 {
    load 2 1.0 
}

constraints Plain
system FullGeneral
numberer Plain
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient

set lambda [eigen -fullGenLapack 2]
set w1 [expr sqrt([lindex $lambda 0])]
set w2 [expr sqrt([lindex $lambda 1])]

set a0 [expr $zeta*2.0*$w1*$w2/($w1 + $w2)];
set a1 [expr $zeta*2.0/($w1 + $w2)]
rayleigh $a0 $a1 0. 0.

recorder Node -file nodeR1.out -time -node 2 3 -dof 1 disp
analyze 1000 .01

set rayleighRES1 [nodeDisp 3 1]

buildModel 

# add load pattern
timeSeries Trig 1 0.0 [expr 100.0*$periodForce] $periodForce -factor $P
pattern Plain 1 1 {
    load 2 1.0 
}

constraints Plain
system FullGeneral
numberer Plain
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient

set lambda [eigen 2]
modalDamping $zeta

recorder Node -file m1.out -time -node 2 3 -dof 1 disp
analyze 1000 .01

set modalRES1 [nodeDisp 3 1]

#
# elCentro motion
#

source ReadRecord.tcl
ReadRecord elCentro.at2 elCentro.dat dt nPts

buildModel 

# add load pattern
timeSeries Path 1 -filePath elCentro.dat -dt $dt -factor $g
pattern UniformExcitation 1 1 -accel 1

set lambda [eigen -fullGenLapack 2]
set w1 [expr sqrt([lindex $lambda 0])]
set w2 [expr sqrt([lindex $lambda 1])]

set a0 [expr $zeta*2.0*$w1*$w2/($w1 + $w2)];
set a1 [expr $zeta*2.0/($w1 + $w2)]
rayleigh $a0 $a1 0. 0.

recorder Node -file r2.out -time -node 2 3 -dof 1 disp
analyze 1000 .01

set rayleighRES2 [nodeDisp 3 1]


buildModel 

timeSeries Path 1 -filePath elCentro.dat -dt $dt -factor $g
pattern UniformExcitation 1 1 -accel 1

constraints Plain
system FullGeneral
numberer Plain
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient

set lambda [eigen 2]
modalDamping $zeta


recorder Node -file m2.out -time -node 2 3 -dof 1 disp
analyze 1000 .01
set modalRES2 [nodeDisp 3 1]

#
# print results
#



#
# check & output results
#

set error1 [expr $modalRES1-$rayleighRES1]
set error2 [expr $modalRES2-$rayleighRES2]
puts "RESULTS Comparing Rayleigh & Modal"
set formatString {%20s%15.5f%15.5f%10s%15.5f}
puts [format $formatString Harmonic: $rayleighRES1 $modalRES1 Diff: $error1]
puts [format $formatString Earthquake: $rayleighRES2 $modalRES2 Diff: $error2]

set results [open results.out a+]
if {[expr abs($rayleighRES1-$modalRES1)] > $myTol || [expr abs($rayleighRES2-$modalRES2)] > $myTol} {
    puts $results "FAILED : mdofModal.tcl"
    puts "FAILED : mdofModal.tcl"
} else {
    puts $results "SUCCESS : mdofModal.tcl"
    puts "SUCCESS : mdofModal.tcl"
}
close $results
    
