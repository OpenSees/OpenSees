# Bathe3D.tcl
#
# Bathe & Bolourchi 45-degree bend: a cantilever curved in the x-z plane
# (radius 100, 45 degrees, 8 elements) with an out-of-plane tip load of 600.
# Compares the 3d corotational transformations, with the element node order
# as given and reversed.
#
# Source: Bathe, K.-J. and Bolourchi, S. (1979), IJNME 14(7):961-986; used
# (Table 8) by Perez, C.M. and Filippou, F.C. (2024) "On nonlinear geometric
# transformations of finite elements", IJNME 125(17), e7506, to compare
# geometric transformations (examples/perez2024nonlinear.py in xara).
#
# Reference (Corot3DBattini / xara Corotational03, Euler-Bernoulli elements):
#   tip (ux, uy, uz) = (-23.7823, 53.6463, -13.6962)
#
# usage: OpenSees Bathe3D.tcl ?elastic|force? ?nsteps?

set eleType elastic
set nsteps  3
if {[llength $argv] > 0} {set eleType [lindex $argv 0]}
if {[llength $argv] > 1} {set nsteps  [lindex $argv 1]}

proc runBathe {transf reverse eleType nsteps} {
    wipe
    model basic -ndm 3 -ndf 6

    set ne 8
    set R  100.0
    set pi [expr acos(-1.0)]
    for {set n 1} {$n <= [expr $ne+1]} {incr n} {
        set a [expr ($n-1)*$pi/4.0/$ne]
        node $n [expr $R*sin($a)] 0.0 [expr $R*(1-cos($a))]
    }
    fix 1 1 1 1 1 1 1

    # unit square section
    set E 1.0e7; set G 0.5e7; set A 1.0
    set I [expr 1.0/12.0]; set J [expr 2*$I*0.8436]

    geomTransf $transf 1 0 1 0
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
    pattern Plain 1 Linear {
        load $tip 0. 600. 0. 0. 0. 0.
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
    return [format "%-20s reverse=%d  tip = %10.4f %10.4f %10.4f   iterations: %s" \
                $transf $reverse [nodeDisp $tip 1] [nodeDisp $tip 2] [nodeDisp $tip 3] $its]
}

puts "Bathe 45-degree bend, $eleType elements, $nsteps load steps"
foreach transf {Corotational CorotationalBattini} {
    foreach reverse {0 1} {
        puts [runBathe $transf $reverse $eleType $nsteps]
    }
}
