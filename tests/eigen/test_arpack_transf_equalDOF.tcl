# ArpackSOE::addM with Transformation-condensed equalDOF.
#
# 2D plane-strain quad column (H=10, N=100), base fixed, left/right nodes at
# every level joined by equalDOF so the mesh is a pure 1D shear beam. Closed
# form: T1 = 4H/Vs with Vs=sqrt(G/rho). G=1e5, rho=2 => Texact = 0.178885 s.
#
# Default `eigen` uses ArpackSOE's diagonal-mass cache. When Transformation
# puts both ends of an element on the same equation, that cache used to
# double-count m(i,i) and inflate T1 by exactly sqrt(2). `eigen -fullGenLapack`
# assembles M correctly and is the control. Assert:
#   1) ARPACK T1 matches Texact
#   2) ARPACK T1 matches fullGenLapack T1
# under constraints Transformation.
#
# eigen needs no analysis of its own: constraints/numberer/system are enough,
# and the eigen command fills in the remaining defaults.

set H   10.0
set N   100
set dy  [expr {$H/$N}]
set w   1.0
set nu  0.0
set G   100000.0
set E   [expr {2.0*$G*(1.0+$nu)}]
set rho 2.0
set Vs  [expr {sqrt($G/$rho)}]
set Texact [expr {4.0*$H/$Vs}]
set pi  [expr {acos(-1.0)}]

# --- ARPACK (the buggy path) ---
wipe
model BasicBuilder -ndm 2 -ndf 2
nDMaterial ElasticIsotropic 1 $E $nu $rho
for {set j 0} {$j <= $N} {incr j} {
    set y [expr {-$H + $j*$dy}]
    node [expr {2*$j+1}] 0.0 $y
    node [expr {2*$j+2}] $w  $y
}
fix 1 1 1
fix 2 1 1
for {set j 1} {$j <= $N} {incr j} {
    element quad $j [expr {2*($j-1)+1}] [expr {2*($j-1)+2}] \
                    [expr {2*$j+2}]     [expr {2*$j+1}] 1.0 PlaneStrain 1
    equalDOF [expr {2*$j+1}] [expr {2*$j+2}] 1 2
}
constraints Transformation
numberer Plain
system UmfPack
set T_arpack [expr {2.0*$pi/sqrt([lindex [eigen 3] 0])}]

# --- fullGenLapack control ---
wipe
model BasicBuilder -ndm 2 -ndf 2
nDMaterial ElasticIsotropic 1 $E $nu $rho
for {set j 0} {$j <= $N} {incr j} {
    set y [expr {-$H + $j*$dy}]
    node [expr {2*$j+1}] 0.0 $y
    node [expr {2*$j+2}] $w  $y
}
fix 1 1 1
fix 2 1 1
for {set j 1} {$j <= $N} {incr j} {
    element quad $j [expr {2*($j-1)+1}] [expr {2*($j-1)+2}] \
                    [expr {2*$j+2}]     [expr {2*$j+1}] 1.0 PlaneStrain 1
    equalDOF [expr {2*$j+1}] [expr {2*$j+2}] 1 2
}
constraints Transformation
numberer Plain
system UmfPack
set T_lapack [expr {2.0*$pi/sqrt([lindex [eigen -fullGenLapack 3] 0])}]

set rExact  [expr {$T_arpack/$Texact}]
set rLapack [expr {$T_arpack/$T_lapack}]
set sqrt2   [expr {sqrt(2.0)}]
puts [format "T_arpack=%.6f  T_lapack=%.6f  Texact=%.6f  T_arpack/Texact=%.6f  T_arpack/T_lapack=%.6f" \
      $T_arpack $T_lapack $Texact $rExact $rLapack]

# Discretization vs continuum exact is ~0.1%; the bug produced exactly sqrt(2).
# ARPACK must match the assembled-matrix path and must not sit near sqrt(2).
if {abs($rLapack - 1.0) > 1.0e-6} {
    puts [format "FAIL: ARPACK/fullGenLapack = %.6f (want ~1)" $rLapack]
    exit 1
}
if {abs($rExact - $sqrt2) < 1.0e-2} {
    puts [format "FAIL: ARPACK T1/Texact = %.6f looks like the sqrt(2) mass-doubling bug" $rExact]
    exit 1
}
if {abs($rExact - 1.0) > 1.0e-3} {
    puts [format "FAIL: ARPACK T1/Texact = %.6f (want ~1 within 0.1%%)" $rExact]
    exit 1
}
puts "PASS"
exit 0
