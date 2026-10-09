# Steel01 DDM sensitivity check: isotropic hardening (a1-a4) plus fy, E0, b
#
# Compares DDM (sensNodeDisp) against centered finite differences at EVERY
# load step of force-controlled histories. Force control matters: after the
# first cycle the peak strain depends on the hardening parameters, so the
# extrema/shift sensitivities become nonzero and the history propagation
# (SHVs rows 2-5) is actually exercised.
#
# Scenarios
#   fullCycles : every reversal sets a NEW extreme (+72,-92,+100,-112).
#                Exercises the "new extremum" branch of the extrema update.
#   partialRev : reversals that happen BEFORE the previous extreme is reached
#                (unload to +20, reload only to +88; unload to -32, reload
#                only to -100). The extrema must be left unchanged and the
#                shift re-evaluated from the old extrema. Exercises the
#                "no new extremum" branch on both the tension and the
#                compression side.
#   noHardening: a1 = a3 = 0, regression for fy, E0, b.
#
# Run:  OpenSees steel01_ddm_isotropic_check.tcl
# Exit code 0 = all checks passed, 1 = at least one failed.
#
# NOTE: the names passed to "parameter ... element 1 <name>" must match what
# Steel01::setParameter accepts (fy / E / b / a1..a4 assumed here).

set convTol 1.0e-12  ;# NormUnbalance tolerance (tight, to keep FD noise low)
set uNoise  1.0e-13  ;# assumed displacement noise per run (about 20x the
                      # strain error implied by convTol and the smallest tangent)
set rtol    1.0e-4   ;# local relative tolerance
set hRel    1.0e-3   ;# FD step, relative to max(|p|, 0.01)
set dP      2.0      ;# load increment per step (peaks must be multiples of it)

# ---- load history from a list of peak loads (yield load = fy = 50.5) ----
proc loadHistory {peaks dP} {
    set vals {0.0}
    set cur 0.0
    foreach pk $peaks {
        set n   [expr {int(round(abs($pk - $cur) / $dP))}]
        set dir [expr {$pk > $cur ? 1.0 : -1.0}]
        for {set k 1} {$k <= $n} {incr k} {
            lappend vals [expr {$cur + $dir * $k * $dP}]
        }
        set cur $pk
    }
    return $vals
}

# ---- one analysis; returns {dispList sensList} (sensList empty if no DDM) ----
proc runModel {p sensName peaks} {
    global dP convTol
    set vals [loadHistory $peaks $dP]
    wipeReliability
    wipe
    model basic -ndm 1 -ndf 1
    node 1 0.0
    node 2 1.0
    fix 1 1
    uniaxialMaterial Steel01 1 [dict get $p fy] [dict get $p E] [dict get $p b] \
        [dict get $p a1] [dict get $p a2] [dict get $p a3] [dict get $p a4]
    element truss 1 1 2 1.0 1
    timeSeries Path 1 -dt 1.0 -useLast -values $vals
    pattern Plain 1 1 {
        load 2 1.0
    }
    constraints Plain
    numberer Plain
    system BandGeneral
    test NormUnbalance $convTol 20 0
    algorithm NewtonLineSearch -type Bisection -tol 0.8 -maxIter 20 -minEta 0.001 -maxEta 10.0
    integrator LoadControl 1.0
    analysis Static

    if {$sensName ne ""} {
        reliability
        parameter 1 element 1 $sensName
        sensitivityIntegrator -static
        sensitivityAlgorithm -computeAtEachStep
    }

    set nSteps [expr {[llength $vals] - 1}]
    set u {}
    set du {}
    for {set i 1} {$i <= $nSteps} {incr i} {
        if {[analyze 1] != 0} {
            error "analysis failed at step $i (param=$sensName)"
        }
        lappend u [nodeDisp 2 1]
        if {$sensName ne ""} {
            lappend du [sensNodeDisp 2 1 1]
        }
    }
    return [list $u $du]
}

proc absmax {lst} {
    set m 0.0
    foreach x $lst { if {abs($x) > $m} {set m [expr {abs($x)}]} }
    return $m
}

# ---- scenarios: label, base parameters, parameters to test, peak loads ----
set hard   {fy 50.5 E 20000.0 b 0.01 a1 0.02 a2 1.0 a3 0.03 a4 1.5}
set noHard {fy 50.5 E 20000.0 b 0.01 a1 0.0  a2 1.0 a3 0.0  a4 1.0}
set allP   {fy E b a1 a2 a3 a4}

set scenarios [list \
    [list fullCycles  $hard   $allP        {72.0 -92.0 100.0 -112.0}] \
    [list partialRev  $hard   $allP        {72.0 -92.0 100.0 20.0 88.0 -112.0 -32.0 -100.0 120.0}] \
    [list noHardening $noHard {fy E b}     {72.0 -92.0 100.0 -112.0}] \
]

set nFail 0
set totalSteps 0
foreach sc $scenarios {
    lassign $sc label base names peaks
    set nSteps [expr {[llength [loadHistory $peaks $dP]] - 1}]
    puts "--- $label: $nSteps steps per run"

    foreach name $names {
        set p0 [dict get $base $name]
        set h  [expr {$hRel * max(abs($p0), 0.01)}]

        # fixed absolute tolerance = FD noise floor for THIS parameter
        # (displacement noise / step), independent of the response history.
        # Always > 0, so err/lim below can never divide by zero.
        set atol [expr {$uNoise / $h}]

        lassign [runModel $base $name $peaks] uBase duDDM
        lassign [runModel [dict replace $base $name [expr {$p0 + $h}]] "" $peaks] uPlus  _
        lassign [runModel [dict replace $base $name [expr {$p0 - $h}]] "" $peaks] uMinus _
        incr totalSteps [expr {3 * $nSteps}]

        set fd {}
        foreach a $uPlus b $uMinus { lappend fd [expr {($a - $b) / (2.0 * $h)}] }

        # guard FIRST: a response at or below the noise floor proves nothing
        set scale [expr {max([absmax $fd], [absmax $duDDM])}]
        if {$scale <= 10.0 * $atol} {
            puts "FAIL  $label  $name : sensitivity is at the noise floor (max [format %.3e $scale]), test does not exercise the parameter"
            incr nFail
            continue
        }

        # local tolerance: fixed floor + relative part based on BOTH values
        set worst 0.0
        set badStep -1
        set i 0
        foreach d $duDDM f $fd {
            incr i
            set err   [expr {abs($d - $f)}]
            set lim   [expr {$atol + $rtol * max(abs($d), abs($f))}]
            set ratio [expr {$err / $lim}]
            if {$ratio > $worst} {set worst $ratio}
            if {$ratio > 1.0 && $badStep < 0} {set badStep $i}
        }

        if {$badStep >= 0} {
            puts "FAIL  $label  $name : first mismatch at step $badStep (worst err/limit = [format %.3g $worst])"
            incr nFail
        } else {
            puts "PASS  $label  $name : max |du/dp| = [format %.4e $scale], worst err/limit = [format %.3g $worst]"
        }
    }
}

puts "total analysis steps: $totalSteps"
if {$nFail > 0} {
    puts "$nFail check(s) FAILED"
    exit 1
}
puts "All Steel01 DDM checks passed"
exit 0
