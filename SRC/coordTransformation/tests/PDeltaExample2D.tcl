# PDeltaExample2D.tcl
#
# Example of the PDelta documentation page (geomTransf/PDelta.rst):
# moment amplification and implied buckling load of a cantilever column,
# 1, 2, 4 and 8 elements. Expected output:
#   exact amplification = 1.986
#   elements 1 / 2 / 4 / 8: amplification 1.698 / 1.893 / 1.961 / 1.980,
#   implied buckling load / Pcr 1.216 / 1.060 / 1.020 / 1.010
#
# cantilever column under axial load P = 0.5 Pcr plus a small lateral tip load H,
# divided into 1, 2, 4 and 8 elements

set L 10.0;  set E 1.0e4;  set A 100.0;  set I 1.0
set pi  [expr acos(-1.0)]
set Pcr [expr $pi*$pi*$E*$I/(4*$L*$L)]    ;# Euler load of a cantilever
set P   [expr 0.5*$Pcr]
set H   1.0e-3

# exact amplification of the lateral tip displacement
set u [expr $pi/2*sqrt($P/$Pcr)]
puts [format "exact amplification = %.3f\n" [expr 3*(tan($u)-$u)/pow($u,3)]]
puts "elements   amplification   implied buckling load / Pcr"

foreach ne {1 2 4 8} {
    wipe
    model basic -ndm 2 -ndf 3

    for {set i 0} {$i <= $ne} {incr i} {
        node [expr $i+1] 0.0 [expr $i*$L/$ne]
    }
    fix 1 1 1 1

    geomTransf PDelta 1
    for {set i 1} {$i <= $ne} {incr i} {
        element elasticBeamColumn $i $i [expr $i+1] $A $E $I 1
    }
    set tip [expr $ne+1]

    # 1) axial load, applied in 10 steps and then held constant
    pattern Plain 1 Linear {
        load $tip 0.0 [expr -$P] 0.0
    }
    system BandGeneral
    numberer Plain
    constraints Plain
    test NormDispIncr 1.0e-12 20
    algorithm Newton
    integrator LoadControl 0.1
    analysis Static
    analyze 10
    loadConst -time 0.0

    # 2) small lateral load at the tip
    pattern Plain 2 Linear {
        load $tip $H 0.0 0.0
    }
    integrator LoadControl 1.0
    analyze 1

    # amplification over the first-order displacement H L^3 / 3EI
    set amp [expr [nodeDisp $tip 1]/($H*pow($L,3)/(3*$E*$I))]
    puts [format "%8d   %13.3f   %27.3f" $ne $amp [expr ($P/$Pcr)/(1.0-1.0/$amp)]]
}
