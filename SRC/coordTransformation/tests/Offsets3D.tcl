# Offsets3D.tcl
#
# Rigid joint offsets (-jntOffset) with CorotationalBattini.
#
# A cantilever of 8 elements along X; every element has rigid offsets
#   voI = (0.10, 0.20, 0.00), voJ = (-0.10, 0.20, 0.05)
# (global, from node to element end), so the element axis does not pass
# through the nodes. The tip node is loaded by forces and moments; load
# factor 1 gives tip rotations of about 0.6 rad, 3 about 1.6 rad.
#
# Checks:
#  - small loads (load factor <= 0.001): Linear with the same -jntOffset, and an
#    explicit model in which each offset is a very stiff element between the
#    node and an extra element-end node, must give the same result (checks the
#    offset convention: global vector from node to element end).
#  - large loads: reference values from the Corotational03 (Battini)
#    transformation of xara with the same -jntOffset, which agree with
#    CorotationalBattini to all printed digits:
#      lf    tip (ux, uy, uz)                   rotation dofs
#      1    (-0.023872, 0.743612, 0.533402)     (0.090314, -0.076588, 0.130908)
#      3    (-0.333480, 2.119159, 1.459035)     (0.316327, -0.163282, 0.405430)
#      8    (-1.844844, 4.730908, 2.572505)     (0.834396,  0.049590, 1.169442)
#      15   (-4.485872, 6.635610, 2.607671)     (1.065884,  0.606682, 2.220649)
#    (1, 10 steps; 3, 20; 8, 40; 15, 60). The stiff-link model is not used for
#    large loads: very stiff corotational elements make Newton fail (the first
#    linear predictor produces huge spurious link forces), with this or any
#    corotational transformation for the links.
#
# usage: OpenSees Offsets3D.tcl ?loadFactor? ?nsteps?   (default 1.0 10)

set lf     1.0
set nsteps 10
if {[llength $argv] > 0} {set lf     [lindex $argv 0]}
if {[llength $argv] > 1} {set nsteps [lindex $argv 1]}

set ne 8; set L 10.0
set E 1.0e4; set G 4.0e3; set A 1.0; set Iy 1.0e-2; set Iz 2.0e-2; set J 1.5e-2
set oI {0.10 0.20 0.00}
set oJ {-0.10 0.20 0.05}

proc analyzeIt {tip nsteps} {
    system UmfPack
    numberer RCM
    constraints Plain
    test NormDispIncr 1.0e-8 60 0
    algorithm Newton
    integrator LoadControl [expr 1.0/$nsteps]
    analysis Static
    set its {}
    for {set s 0} {$s < $nsteps} {incr s} {
        if {[analyze 1] != 0} {return "diverged at step [expr $s+1]"}
        lappend its [testIter]
    }
    return [format "tip = %10.6f %10.6f %10.6f  rot = %9.6f %9.6f %9.6f  its: %s" \
                [nodeDisp $tip 1] [nodeDisp $tip 2] [nodeDisp $tip 3] \
                [nodeDisp $tip 4] [nodeDisp $tip 5] [nodeDisp $tip 6] $its]
}

proc loads {tip lf} {
    pattern Plain 1 Linear {
        load $tip 0.0 [expr 0.3*$lf] [expr 0.2*$lf] [expr 0.5*$lf] 0.0 [expr 1.5*$lf]
    }
}

# model 1: offsets given to the transformation
proc modelOffsets {transf ne L E G A Iy Iz J oI oJ lf nsteps} {
    wipe
    model basic -ndm 3 -ndf 6
    for {set n 1} {$n <= [expr $ne+1]} {incr n} { node $n [expr ($n-1)*$L/$ne] 0.0 0.0 }
    fix 1 1 1 1 1 1 1
    geomTransf $transf 1 0 0 1 -jntOffset {*}$oI {*}$oJ
    for {set e 1} {$e <= $ne} {incr e} {
        element elasticBeamColumn $e $e [expr $e+1] $A $E $G $J $Iy $Iz 1
    }
    set tip [expr $ne+1]
    loads $tip $lf
    return [analyzeIt $tip $nsteps]
}

# model 2: offsets as explicit stiff elements
proc modelLinks {transf ne L E G A Iy Iz J oI oJ lf nsteps} {
    wipe
    model basic -ndm 3 -ndf 6
    for {set n 1} {$n <= [expr $ne+1]} {incr n} { node $n [expr ($n-1)*$L/$ne] 0.0 0.0 }
    fix 1 1 1 1 1 1 1
    geomTransf $transf 1 0 0 1
    geomTransf $transf 2 0 0 1
    geomTransf $transf 3 0 1 0
    set big 1.0e4
    set tag 100
    for {set e 1} {$e <= $ne} {incr e} {
        set xI [expr ($e-1)*$L/$ne]; set xJ [expr $e*$L/$ne]
        set a [expr 1000+$e]; set b [expr 2000+$e]
        node $a [expr $xI+[lindex $oI 0]] [lindex $oI 1] [lindex $oI 2]
        node $b [expr $xJ+[lindex $oJ 0]] [lindex $oJ 1] [lindex $oJ 2]
        element elasticBeamColumn $e $a $b $A $E $G $J $Iy $Iz 1
        # stiff links node -> element end (choose a transf not parallel to link)
        foreach {n m o} [list $e $a $oI [expr $e+1] $b $oJ] {
            set t 2
            if {abs([lindex $o 2]) > 0.5*sqrt([lindex $o 0]**2+[lindex $o 1]**2+[lindex $o 2]**2)} {set t 3}
            element elasticBeamColumn [incr tag] $n $m 1.0 [expr $E*$big] [expr $G*$big] 1.0 1.0 1.0 $t
        }
    }
    set tip [expr $ne+1]
    loads $tip $lf
    return [analyzeIt $tip $nsteps]
}

puts "Rigid offsets: load factor $lf, $nsteps steps"
if {$lf <= 0.001} {
    puts [format "%-34s %s" "Linear, -jntOffset"             [modelOffsets Linear $ne $L $E $G $A $Iy $Iz $J $oI $oJ $lf $nsteps]]
}
puts [format "%-34s %s" "CorotationalBattini, -jntOffset" [modelOffsets CorotationalBattini $ne $L $E $G $A $Iy $Iz $J $oI $oJ $lf $nsteps]]
if {$lf <= 0.001} {
    puts [format "%-34s %s" "CorotationalBattini, stiff links"  [modelLinks CorotationalBattini $ne $L $E $G $A $Iy $Iz $J $oI $oJ $lf $nsteps]]
}
