# PDeltaColumn2D.tcl
#
# Second-order (P-Delta / P-delta) effects in a 2d cantilever column.
# The column (L = 10, EI = 1e4) carries an axial compression P = frac*Pcr,
# Pcr = pi^2 EI / (4 L^2) (Euler load), applied first and held constant.
# A small lateral tip load H is then added. The amplification of the
# lateral tip displacement over the first-order value H L^3/(3 EI) is
# compared with the exact value
#
#   amp = 3 (tan(u) - u) / u^3,   u = (pi/2) sqrt(frac)
#
# and the implied buckling load P/(1 - 1/amp) is compared with Pcr.
#
# PDelta captures only the drift between element ends (P-Delta); the effect
# of member curvature (P-delta) needs several elements per member. The same
# holds for Corotational. Linear gives no amplification.
#
# usage: OpenSees PDeltaColumn2D.tcl ?frac?

set frac 0.5
if {[llength $argv] > 0} {set frac [lindex $argv 0]}

proc runColumn {transf ne frac} {
    wipe
    model basic -ndm 2 -ndf 3

    set L 10.0; set E 1.0e4; set A 100.0; set I 1.0
    set pi  [expr acos(-1.0)]
    set Pcr [expr $pi*$pi*$E*$I/(4*$L*$L)]

    for {set n 1} {$n <= [expr $ne+1]} {incr n} {
        node $n 0.0 [expr ($n-1)*$L/$ne]
    }
    fix 1 1 1 1

    geomTransf $transf 1
    for {set e 1} {$e <= $ne} {incr e} {
        element elasticBeamColumn $e $e [expr $e+1] $A $E $I 1
    }
    set tip [expr $ne+1]

    # 1) axial load, held constant
    pattern Plain 1 Linear {
        load $tip 0. [expr -$frac*$Pcr] 0.
    }
    system FullGeneral
    numberer Plain
    constraints Plain
    test NormDispIncr 1.0e-12 50 0
    algorithm Newton
    integrator LoadControl 0.1
    analysis Static
    analyze 10
    loadConst -time 0.0

    # 2) small lateral load
    set H 1.0e-3
    pattern Plain 2 Linear {
        load $tip $H 0. 0.
    }
    integrator LoadControl 1.0
    analyze 1

    set d0  [expr $H*$L*$L*$L/(3*$E*$I)]
    set amp [expr [nodeDisp $tip 1]/$d0]
    if {abs($amp - 1.0) < 1.0e-12} {
        set ratio "-"
    } else {
        set ratio [format "%.3f" [expr $frac/(1.0 - 1.0/$amp)]]
    }
    return [format "%-13s ne=%d   amplification = %.4f   implied Pcr/Pcr_exact = %s" \
                $transf $ne $amp $ratio]
}

set pi [expr acos(-1.0)]
set u  [expr $pi/2.0*sqrt($frac)]
puts "2d cantilever column, P = $frac Pcr, small lateral tip load"
puts [format "exact amplification = %.4f" [expr 3*(tan($u)-$u)/pow($u,3)]]
foreach transf {Linear PDelta Corotational} {
    foreach ne {1 2 4 8} {
        puts [runColumn $transf $ne $frac]
    }
}
