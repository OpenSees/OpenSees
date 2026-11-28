# Rollup3D.tcl
#
# Cantilever rolled into a circle by an end moment (NAFEMS NL5 type).
# L = 10, EI = 100, 10 elements, end moment M = frac * 2*pi*EI/L about the
# local z axis. Exact tip displacement for curvature k = 2*pi*frac/L:
#   ux = sin(kL)/k - L,  uy = (1 - cos(kL))/k
# frac = 1 rolls the beam into a full circle: tip (-10, 0, 0).
#
# Source: adapted from Section 4.2.1 of Perez, C.M. and Filippou, F.C. (2024)
# "On nonlinear geometric transformations of finite elements", IJNME 125(17),
# e7506 (examples/perez2024nonlinear.py in xara), which uses the same
# cantilever (L = 10, EI = 100, 10 elements) with frac = 1/8 and
# shear-deformable elements; here Euler-Bernoulli elements and, by default,
# the full circle.
#
# usage: OpenSees Rollup3D.tcl ?elastic|force? ?nsteps? ?frac?

set eleType elastic
set nsteps  3
set frac    1.0
if {[llength $argv] > 0} {set eleType [lindex $argv 0]}
if {[llength $argv] > 1} {set nsteps  [lindex $argv 1]}
if {[llength $argv] > 2} {set frac    [lindex $argv 2]}

proc runRollup {transf reverse eleType nsteps frac} {
    wipe
    model basic -ndm 3 -ndf 6

    set ne 10; set L 10.0
    set E 1.0e4; set G 1.0e4; set A 1.0; set I 1.0e-2; set J 1.0e-2
    set pi [expr acos(-1.0)]

    for {set n 1} {$n <= [expr $ne+1]} {incr n} {
        node $n [expr ($n-1)*$L/$ne] 0.0 0.0
    }
    fix 1 1 1 1 1 1 1

    geomTransf $transf 1 0 0 1
    if {$eleType == "force"} {
        section Elastic 1 $E $A $I $I $G $J
        beamIntegration Lobatto 1 1 5
    }
    for {set e 1} {$e <= $ne} {incr e} {
        set i $e; set j [expr $e+1]
        if {$reverse} {set t $i; set i $j; set j $t}
        if {$eleType == "force"} {
            element forceBeamColumn $e $i $j 1 1
        } else {
            element elasticBeamColumn $e $i $j $A $E $G $J $I $I 1
        }
    }

    set tip [expr $ne+1]
    set M [expr $frac*2*$pi*$E*$I/$L]
    pattern Plain 1 Linear {
        load $tip 0. 0. 0. 0. 0. $M
    }

    system FullGeneral
    numberer Plain
    constraints Plain
    test NormDispIncr 1.0e-10 60 0
    algorithm Newton
    integrator LoadControl [expr 1.0/$nsteps]
    analysis Static

    set its {}
    for {set s 0} {$s < $nsteps} {incr s} {
        if {[analyze 1] != 0} {
            return [format "%-20s reverse=%d  diverged in step %d" $transf $reverse [expr $s+1]]
        }
        lappend its [testIter]
    }
    return [format "%-20s reverse=%d  tip = %10.6f %10.6f   iterations: %s" \
                $transf $reverse [nodeDisp $tip 1] [nodeDisp $tip 2] $its]
}

set pi [expr acos(-1.0)]
set k  [expr 2*$pi*$frac/10.0]
puts "Cantilever roll-up, frac = $frac, $eleType elements, $nsteps load steps"
puts [format "exact                           tip = %10.6f %10.6f" \
          [expr sin($k*10.0)/$k - 10.0] [expr (1-cos($k*10.0))/$k]]
foreach transf {Corotational CorotationalBattini} {
    foreach reverse {0 1} {
        puts [runRollup $transf $reverse $eleType $nsteps $frac]
    }
}
