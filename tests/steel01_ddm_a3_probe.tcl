# Minimal Steel01 a3 sensitivity probe

wipe
model basic -ndm 1 -ndf 1

node 1 0.0
node 2 1.0
fix 1 1

set fy 50.0
set E0 20000.0
set b 0.01
set a1 0.0
set a2 1.0
set a3 0.02
if {$argc > 0} {
    set a3 [lindex $argv 0]
}
set a4 1.0

uniaxialMaterial Steel01 1 $fy $E0 $b $a1 $a2 $a3 $a4
element truss 1 1 2 1.0 1

# Cyclic history: decrease to -80, then reverse to +120
set loadValues {}

for {set p 0} {$p >= -80} {incr p -1} {
    lappend loadValues $p
}

for {set p -79} {$p <= 120} {incr p} {
    lappend loadValues $p
}

timeSeries Path 1 -dt 1.0 -useLast -values {*}$loadValues
pattern Plain 1 1 {
    load 2 1.0
}

constraints Plain
numberer Plain
system BandGeneral
test NormUnbalance 1.0e-10 20 0
algorithm NewtonLineSearch -type Bisection -tol 0.8 -maxIter 20 -minEta 0.001 -maxEta 10.0
integrator LoadControl 1.0
analysis Static

reliability
parameter 1 element 1 a3
sensitivityIntegrator -static
sensitivityAlgorithm -computeAtEachStep

set nSteps [expr {[llength $loadValues] - 1}]
set ok [analyze $nSteps]

set u [nodeDisp 2 1]
set finalTime [getTime]
set finalLoad [getLoadFactor 1]
set finalStress [eleResponse 1 material stress]
set finalStrain [eleResponse 1 material strain]
set duDDM [sensNodeDisp 2 1 1]

puts "a3=$a3"
puts "finalTime=$finalTime"
puts "finalLoad=$finalLoad"
puts "finalStress=$finalStress"
puts "finalStrain=$finalStrain"
puts "analysisCode=$ok"
puts "finalDisp=$u"
puts "DDM_du_da3=$duDDM"
