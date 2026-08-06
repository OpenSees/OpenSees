# Retained-node SP + equalDOF under Transformation (multi-step).
#
# Usage:
#   OpenSees test_retained_sp.tcl partial   ;# default: Uy fixed, Ux free on retained
#   OpenSees test_retained_sp.tcl full      ;# both DOFs fixed on retained
#
# One mesh for both modes: inclined trusses 1--2 and 3--4 meeting at a
# coincident apex (equalDOF 2->3), plus a horizontal tip truss 2--5.
# Horizontal load is applied at the free tip (node 5). Five static
# LoadControl steps under constraints Transformation.
#
# partial: retained node 2 has UY fixed and UX free. Node 3 must track
# node 2 every step; final ux2 = P/(EA/L)_inclined = 100*sqrt(2)/10000
# and ux5 = ux2 + P*L_tip/(EA) = ux2 + 0.02.
#
# full: retained node 2 fully fixed. Apex stays at 0; tip still moves
# ux5 = P*L_tip/(EA) = 0.02 (nontrivial free DOF under load).
wipe

set mode partial
if {$argc >= 1} {
    set mode [lindex $argv 0]
}
if {$mode ni {partial full}} {
    puts stderr "ERROR: mode must be partial or full (got '$mode')"
    exit 1
}

model basic -ndm 2 -ndf 2

node 1 0.0 0.0
node 2 1.0 1.0
node 3 1.0 1.0
node 4 2.0 0.0
node 5 3.0 1.0

fix 1 1 1
fix 4 1 1
fix 5 0 1
if {$mode eq "full"} {
    fix 2 1 1
} else {
    # partial fix on the RETAINED node: Uy fixed, Ux free
    fix 2 0 1
}

equalDOF 2 3   1 2

uniaxialMaterial Elastic 1 1000.0
element truss 1 1 2 10.0 1
element truss 2 3 4 10.0 1
element truss 3 2 5 10.0 1

timeSeries Linear 1
pattern Plain 1 1 {
    load 5 100.0 0.0
}

constraints Transformation
numberer Plain
system UmfPack
test NormUnbalance 1.0e-8 2 0
algorithm Newton
integrator LoadControl 0.2
analysis Static

set tol 1.0e-10
set fail 0
set nSteps 5
for {set step 1} {$step <= $nSteps} {incr step} {
    if {[analyze 1] != 0} {
        puts "FAIL: analysis did not converge at step $step (mode=$mode)"
        exit 1
    }
    set ux2 [nodeDisp 2 1]
    set uy2 [nodeDisp 2 2]
    set ux3 [nodeDisp 3 1]
    set uy3 [nodeDisp 3 2]
    set ux5 [nodeDisp 5 1]
    set diffX [expr {abs($ux2 - $ux3)}]
    set diffY [expr {abs($uy2 - $uy3)}]
    puts [format "mode=%s step %d: retained=(% .8e,% .8e) constrained=(% .8e,% .8e) tip_ux=% .8e |dX|=%.3e |dY|=%.3e" \
        $mode $step $ux2 $uy2 $ux3 $uy3 $ux5 $diffX $diffY]
    if {$diffX > $tol || $diffY > $tol} { set fail 1 }
    if {[expr {abs($uy2)}] > $tol || [expr {abs($uy3)}] > $tol} { set fail 1 }
    if {$mode eq "full"} {
        if {[expr {abs($ux2)}] > $tol || [expr {abs($ux3)}] > $tol} { set fail 1 }
    }
}

if {$fail} {
    puts "FAIL: constrained node does not track retained node (mode=$mode)"
    exit 1
}

set ux2 [nodeDisp 2 1]
set ux5 [nodeDisp 5 1]
set ux_tip_rel [expr {100.0 * 2.0 / 10000.0}]
if {$mode eq "partial"} {
    # Inclined pair: K_xx = EA/L = 10000/sqrt(2); tip spring carries P with ux5-ux2 = P*L/EA
    set ux_apex [expr {100.0 * sqrt(2.0) / 10000.0}]
    set ux_tip [expr {$ux_apex + $ux_tip_rel}]
    if {[expr {abs($ux2 - $ux_apex)}] > 1.0e-9} {
        puts [format "FAIL: wrong apex ux2=%.6e expected %.6e" $ux2 $ux_apex]
        exit 1
    }
    if {[expr {abs($ux5 - $ux_tip)}] > 1.0e-9} {
        puts [format "FAIL: wrong tip ux5=%.6e expected %.6e" $ux5 $ux_tip]
        exit 1
    }
} else {
    if {[expr {abs($ux5 - $ux_tip_rel)}] > 1.0e-9} {
        puts [format "FAIL: wrong tip ux5=%.6e expected %.6e (apex should be fixed)" $ux5 $ux_tip_rel]
        exit 1
    }
    if {[expr {abs($ux5)}] < 1.0e-12} {
        puts "FAIL: expected nonzero tip motion in full mode"
        exit 1
    }
}

puts "PASS mode=$mode"
exit 0
