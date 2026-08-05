# 9_4_QuadUP element with equalDOF on top corners and on opposing midsides.
#
# One 9-node QuadUP element with mixed ndf (corners 3: UX,UY,P; midsides and
# centre 2: UX,UY). The base is fully restrained, the top corners are tied to
# each other (UX,UY) and the two side midsides are tied to each other (UX,UY),
# with a pressure datum on both top corners. Loading is purely vertical.
#
# Analysis: one pseudo-static Newmark step under constraints Transformation,
# Newton, NormUnbalance. Checks that both equalDOF pairs stay exactly matched,
# that the column settles vertically, and that no lateral drift appears under
# a symmetric vertical load.
#
# Exercises TransformationFE's Ti^T K(i,j) Tj products. K(i,j) is the element
# tangent block between nodes i and j; with mixed ndf that block is rectangular
# (e.g. 3x2 between a corner and a midside). The old
# Matrix::addMatrixTripleProduct assumed a square middle factor and malformed
# the condensed tangent, so this analysis diverged on the first step. (equalDOF
# T matrices here are still square; the non-square piece is K(i,j).)
#
# The constraint layout deliberately keeps every retained node free of SPs on
# its tied DOFs, so the lateral-drift check is pure roundoff and this test
# depends on no other fix.

wipe
model basic -ndm 2 -ndf 3
node 1 0.0 0.0
node 2 1.0 0.0
node 3 1.0 1.0
node 4 0.0 1.0

model basic -ndm 2 -ndf 2
node 5 0.5 0.0
node 6 1.0 0.5
node 7 0.5 1.0
node 8 0.0 0.5
node 9 0.5 0.5

model basic -ndm 2 -ndf 3
fix 1 1 1 0
fix 2 1 1 0
fix 3 0 0 1
fix 4 0 0 1
model basic -ndm 2 -ndf 2
fix 5 1 1

model basic -ndm 2 -ndf 3
equalDOF 4 3 1 2
model basic -ndm 2 -ndf 2
equalDOF 8 6 1 2

model basic -ndm 2 -ndf 3
nDMaterial ElasticIsotropic 1 1.0e5 0.3 2.0
element 9_4_QuadUP 1 1 2 3 4 5 6 7 8 9 1.0 1 2.2e6 1.0 1.0e-3 1.0e-3 0.0 -9.81

pattern Plain 1 Constant {
    load 4 0.0 -50.0 0
    load 3 0.0 -50.0 0
}
model basic -ndm 2 -ndf 2
pattern Plain 2 Constant {
    load 7 0.0 -50.0
}

numberer Plain
system UmfPack
test NormUnbalance 1.0e-6 2 0
constraints Transformation
# QuadUP has no static formulation (the fluid coupling sits in the damping
# matrix), so gravity is applied as a pseudo-static transient: gamma=1.5 with
# beta=1.0 is the most dissipative unconditionally stable Newmark pair
# (stability needs beta >= (gamma+0.5)^2/4), and one large step lands on the
# steady drained solution.
integrator Newmark 1.5 1.0
algorithm Newton
analysis Transient

set ok [analyze 1 1000]
if {$ok != 0} {
    puts "FAIL test_quadup_shear_column: analyze returned $ok"
    exit 1
}

set u3 [nodeDisp 3]
set u4 [nodeDisp 4]
set u6 [nodeDisp 6]
set u8 [nodeDisp 8]

# tied pairs must match to roundoff
set dux43 [expr {abs([lindex $u3 0]-[lindex $u4 0])}]
set duy43 [expr {abs([lindex $u3 1]-[lindex $u4 1])}]
set dux86 [expr {abs([lindex $u6 0]-[lindex $u8 0])}]
set duy86 [expr {abs([lindex $u6 1]-[lindex $u8 1])}]
if {$dux43 > 1.0e-12 || $duy43 > 1.0e-12 || $dux86 > 1.0e-12 || $duy86 > 1.0e-12} {
    puts [format "FAIL test_quadup_shear_column: equalDOF violated (%.3e %.3e %.3e %.3e)" \
          $dux43 $duy43 $dux86 $duy86]
    exit 1
}

# the column must actually settle under the vertical load
set uy [lindex $u3 1]
if {$uy > -1.0e-4} {
    puts "FAIL test_quadup_shear_column: no vertical settlement, uy=$uy"
    exit 1
}

# symmetric vertical load on a symmetric mesh: no lateral drift anywhere
set uxmax 0.0
foreach n {1 2 3 4} {
    set v [expr {abs([lindex [nodeDisp $n] 0])}]
    if {$v > $uxmax} {set uxmax $v}
}
if {$uxmax > 1.0e-12} {
    puts "FAIL test_quadup_shear_column: spurious lateral drift ux=$uxmax"
    exit 1
}

puts [format "PASS test_quadup_shear_column (uy=%.6e max|ux|=%.3e)" $uy $uxmax]
exit 0
