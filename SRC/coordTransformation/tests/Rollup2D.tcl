# Rollup2D.tcl
#
# 2d cantilever rolled into a circle by an end moment (NAFEMS NL5 type).
# L = 10, EI = 100, end moment M = frac * 2*pi*EI/L.
# Exact tip displacement for curvature k = 2*pi*frac/L:
#   ux = sin(kL)/k - L,  uy = (1 - cos(kL))/k,  rotation = kL
# frac = 1 rolls the beam into a full circle: tip (-10, 0), rotation 2*pi.
#
# Source: 2d version of the cantilever of Section 4.2.1 of Perez, C.M. and
# Filippou, F.C. (2024) "On nonlinear geometric transformations of finite
# elements", IJNME 125(17), e7506 (examples/perez2024nonlinear.py in xara).
#
# A corotational element with a constant moment reproduces the circle
# exactly at the nodes, so with frac = 1 the exact tip is obtained for any
# number of elements (ne >= 2). Small ne gives large rotations per element.
#
# Compares the 2d transformations Corotational and PDelta, with the
# element node order as given and reversed.
#
# usage: OpenSees Rollup2D.tcl ?elastic|force? ?nsteps? ?frac? ?ne?

set eleType elastic
set nsteps  3
set frac    1.0
set ne      10
if {[llength $argv] > 0} {set eleType [lindex $argv 0]}
if {[llength $argv] > 1} {set nsteps  [lindex $argv 1]}
if {[llength $argv] > 2} {set frac    [lindex $argv 2]}
if {[llength $argv] > 3} {set ne      [lindex $argv 3]}

proc runRollup2D {transf reverse eleType nsteps frac ne} {
    wipe
    model basic -ndm 2 -ndf 3

    set L 10.0; set E 1.0e4; set A 1.0; set I 1.0e-2
    set pi [expr acos(-1.0)]

    for {set n 1} {$n <= [expr $ne+1]} {incr n} {
        node $n [expr ($n-1)*$L/$ne] 0.0
    }
    fix 1 1 1 1

    geomTransf $transf 1
    if {$eleType == "force"} {
        section Elastic 1 $E $A $I
        beamIntegration Lobatto 1 1 5
    }
    for {set e 1} {$e <= $ne} {incr e} {
        set i $e; set j [expr $e+1]
        if {$reverse} {set t $i; set i $j; set j $t}
        if {$eleType == "force"} {
            element forceBeamColumn $e $i $j 1 1
        } else {
            element elasticBeamColumn $e $i $j $A $E $I 1
        }
    }

    set tip [expr $ne+1]
    set M [expr $frac*2*$pi*$E*$I/$L]
    pattern Plain 1 Linear {
        load $tip 0. 0. $M
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
            return [format "%-14s reverse=%d  diverged in step %d" $transf $reverse [expr $s+1]]
        }
        lappend its [testIter]
    }
    return [format "%-14s reverse=%d  tip = %10.6f %10.6f %10.6f   iterations: %s" \
                $transf $reverse [nodeDisp $tip 1] [nodeDisp $tip 2] [nodeDisp $tip 3] $its]
}

set pi [expr acos(-1.0)]
set k  [expr 2*$pi*$frac/10.0]
puts "2d cantilever roll-up, frac = $frac, $ne $eleType elements, $nsteps load steps"
puts [format "exact                     tip = %10.6f %10.6f %10.6f" \
          [expr sin($k*10.0)/$k - 10.0] [expr (1-cos($k*10.0))/$k] [expr $k*10.0]]
foreach transf {Corotational PDelta} {
    foreach reverse {0 1} {
        puts [runRollup2D $transf $reverse $eleType $nsteps $frac $ne]
    }
}
