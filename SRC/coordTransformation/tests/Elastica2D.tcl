# Elastica2D.tcl
#
# 2d cantilever with a large transverse tip load P (fixed direction),
# the classic elastica benchmark (Mattiasson 1981). L = 10, EI = 100.
# The load is increased to P L^2/EI = 10 and the tip displacement is
# compared with the exact solution at P L^2/EI = 1, 2, 3, 5, 10.
#
# Exact values (computed by shooting on EI theta'' = -P cos(theta);
# they agree with Mattiasson's published table):
#   PL^2/EI   u/L (shortening)   w/L (deflection)   tip rotation
#      1         0.05643            0.30172           0.46135
#      2         0.16064            0.49346           0.78175
#      3         0.25442            0.60325           0.98602
#      5         0.38763            0.71379           1.21537
#     10         0.55500            0.81061           1.43029
#
# Compares the 2d transformations Corotational, PDelta and Linear, with the
# element node order as given and reversed.
#
# usage: OpenSees Elastica2D.tcl ?elastic|force? ?ne? ?stepsPerUnit?

set eleType elastic
set ne      10
set spu     1
if {[llength $argv] > 0} {set eleType [lindex $argv 0]}
if {[llength $argv] > 1} {set ne      [lindex $argv 1]}
if {[llength $argv] > 2} {set spu     [lindex $argv 2]}

set exact {
    1  0.05643 0.30172 0.46135
    2  0.16064 0.49346 0.78175
    3  0.25442 0.60325 0.98602
    5  0.38763 0.71379 1.21537
    10 0.55500 0.81061 1.43029
}

proc runElastica2D {transf reverse eleType ne spu exact} {
    wipe
    model basic -ndm 2 -ndf 3

    set L 10.0; set E 1.0e4; set A 1.0e3; set I 1.0e-2
    set EI [expr $E*$I]

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

    # reference load: P L^2/EI = 1 at load factor 1
    set tip [expr $ne+1]
    set P [expr $EI/($L*$L)]
    pattern Plain 1 Linear {
        load $tip 0. $P 0.
    }

    system FullGeneral
    numberer Plain
    constraints Plain
    test NormDispIncr 1.0e-10 60 0
    algorithm Newton
    integrator LoadControl [expr 1.0/$spu]
    analysis Static

    puts [format "%s  (reverse=%d)" $transf $reverse]
    puts "   PL2/EI      u/L     exact      w/L     exact    rotation   exact    iterations"
    set lambda 0.0
    set its {}
    foreach {k u0 w0 t0} $exact {
        while {$lambda < $k - 1.0e-9} {
            if {[analyze 1] != 0} {
                puts [format "   diverged at PL2/EI = %.3f" [expr $lambda + 1.0/$spu]]
                return
            }
            set lambda [getTime]
            lappend its [testIter]
        }
        puts [format "   %6.1f   %8.5f %8.5f   %8.5f %8.5f   %8.5f %8.5f   %s" $k \
                  [expr -[nodeDisp $tip 1]/$L] $u0 [expr [nodeDisp $tip 2]/$L] $w0 \
                  [nodeDisp $tip 3] $t0 $its]
        set its {}
    }
}

puts "2d cantilever elastica, $ne $eleType elements, $spu load step(s) per unit of PL^2/EI"
foreach transf {Corotational PDelta Linear} {
    foreach reverse {0 1} {
        runElastica2D $transf $reverse $eleType $ne $spu $exact
    }
}
