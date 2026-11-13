# =============================================================================
# Verification tests for the EmbeddedNodeContact element
#   run:  OpenSees test_EmbeddedNodeContact.tcl
# Each check prints PASS/FAIL; the script ends with a summary and exits with
# code 0 only if all checks pass.
# =============================================================================

set ::nPass 0
set ::nFail 0
proc check {name val ref tol} {
    set err [expr {abs($val - $ref)}]
    if {$err <= $tol} {
        incr ::nPass
        puts [format "  PASS  %-55s value=% .6e  ref=% .6e" $name $val $ref]
    } else {
        incr ::nFail
        puts [format "  FAIL  %-55s value=% .6e  ref=% .6e  err=%.3e" $name $val $ref $err]
    }
}

# a free dummy dof (spring to ground) so that fully prescribed models
# never produce an empty system of equations
proc addDummy {ndm} {
    if {$ndm == 3} {
        node 900 10 10 10; node 901 10 10 10
        fix 901 1 1 1
    } else {
        node 900 10 10; node 901 10 10
        fix 901 1 1
    }
    uniaxialMaterial Elastic 900 1.0
    if {$ndm == 3} {
        element zeroLength 900 901 900 -mat 900 900 900 -dir 1 2 3
    } else {
        element zeroLength 900 901 900 -mat 900 900 -dir 1 2
    }
}

proc staticAnalysis {dt} {
    constraints Transformation
    numberer Plain
    system FullGeneral
    test NormDispIncr 1.0e-10 50 0
    algorithm Newton
    integrator LoadControl $dt
    analysis Static
}

# advance the analysis up to pseudo-time t (steps of size dt)
proc runTo {t dt} {
    while {[getTime] < $t - 0.5*$dt} {
        if {[analyze 1] != 0} {
            puts "  analysis failed at time [getTime]"
            return -1
        }
    }
    return 0
}

# ---------------------------------------------------------------------------
# Local frame for -orient 0 0 1 (same as ZeroLengthContactASDimplex): n=+z, t1=+y, t2=-x
# 3D host tetrahedron (fixed) + embedded node 5, prescribed displacement path
#   Kn = Kt = 1e4, mu = 0.5
#   z: 0 -> -0.01 (t=1) ... hold ... -> +0.01 (t=7, separation)
#   x: 0 (t=1) -> 0.02 (t=3) -> 0.0 (t=5)
# ---------------------------------------------------------------------------
proc buildTet3D {intType} {
    wipe
    model BasicBuilder -ndm 3 -ndf 3
    node 1 0 0 0; node 2 1 0 0; node 3 0 1 0; node 4 0 0 1
    node 5 0.2 0.2 0.2
    fix 1 1 1 1; fix 2 1 1 1; fix 3 1 1 1; fix 4 1 1 1
    element EmbeddedNodeContact 1 5 1 2 3 4 1.0e4 1.0e4 0.5 -orient 0 0 1 -intType $intType
    addDummy 3
    timeSeries Path 1 -time {0 1 3 5 6 7 100} -values {0 0 0.02 0.0 0.0 0.0 0.0}
    timeSeries Path 2 -time {0 1 6 7 100} -values {0 -0.01 -0.01 0.01 0.01}
    timeSeries Constant 3
    pattern Plain 1 1 { sp 5 1 1.0 }
    pattern Plain 2 2 { sp 5 3 1.0 }
    pattern Plain 3 3 { sp 5 2 0.0 }
}

puts "\n=== T1: 3D tetra host, implicit: compression, stick, slip, reversal, separation ==="
buildTet3D 0
set dt 0.05
staticAnalysis $dt
runTo 1.0 $dt
check "N at full compression (Kn*dz = -100)" [eleResponse 1 normalContactForce] -100.0 1e-6
check "contact state = stick (1)" [eleResponse 1 contactState] 1 0
reactions
set N {0.4 0.2 0.2 0.2}
set sumZ [nodeReaction 5 3]
for {set i 1} {$i <= 4} {incr i} {
    set Ni [lindex $N [expr {$i-1}]]
    check "host node $i reaction z = N_i*100" [nodeReaction $i 3] [expr {$Ni*100.0}] 1e-6
    set sumZ [expr {$sumZ + [nodeReaction $i 3]}]
}
check "self-equilibrium (sum of z reactions)" $sumZ 0.0 1e-8
runTo 1.25 $dt
check "stick: T = Kt*ux = 25" [eleResponse 1 tangentialContactForce] 25.0 1e-6
runTo 3.0 $dt
check "slip: T = mu*|N| = 50" [eleResponse 1 tangentialContactForce] 50.0 1e-6
check "contact state = slip (2)" [eleResponse 1 contactState] 2 0
check "local jump along t2 (= -x) = -0.02" [lindex [eleResponse 1 localDisplacement] 2] -0.02 1e-10
runTo 3.5 $dt
check "elastic unloading: T2 = 0 at ux = 0.015" [lindex [eleResponse 1 localForce] 2] 0.0 1e-6
runTo 5.0 $dt
check "reverse slip: T2 = +50 (t2 = -x)" [lindex [eleResponse 1 localForce] 2] 50.0 1e-6
runTo 7.0 $dt
check "separation: N = 0" [eleResponse 1 normalContactForce] 0.0 1e-8
check "separation: |T| = 0" [eleResponse 1 tangentialContactForce] 0.0 1e-6
check "contact state = open (0)" [eleResponse 1 contactState] 0 0
check "separation: global force on Cnode z = 0" [lindex [eleResponse 1 force] 2] 0.0 1e-8

puts "\n=== T2: IMPL-EX (-intType 1), same path ==="
buildTet3D 1
staticAnalysis $dt
runTo 1.0 $dt
check "N at full compression (IMPL-EX)" [eleResponse 1 normalContactForce] -100.0 1e-6
runTo 3.0 $dt
check "slip T ~ 50 (IMPL-EX, 1% tol)" [eleResponse 1 tangentialContactForce] 50.0 0.5
runTo 5.0 $dt
check "reverse slip T2 ~ +50 (IMPL-EX, 1% tol)" [lindex [eleResponse 1 localForce] 2] 50.0 0.5
runTo 7.0 $dt
runTo 7.5 $dt
check "separation N = 0 (IMPL-EX)" [eleResponse 1 normalContactForce] 0.0 1e-8

puts "\n=== T3: linear displacement field of host + embedded node -> zero force ==="
# u = A*X with A = [[0.01, 0.02, 0.0], [-0.03, 0.0, 0.01], [0.02, 0.01, -0.02]] + rigid (0.1,-0.2,0.3)
wipe
model BasicBuilder -ndm 3 -ndf 3
set XYZ {{0 0 0} {1 0 0} {0 1 0} {0 0 1} {0.2 0.2 0.2}}
for {set i 1} {$i <= 5} {incr i} {
    eval node $i [lindex $XYZ [expr {$i-1}]]
}
element EmbeddedNodeContact 1 5 1 2 3 4 1.0e4 1.0e4 0.5 -orient 0.3 0.4 0.866
addDummy 3
timeSeries Linear 1
pattern Plain 1 1 {
    for {set i 1} {$i <= 5} {incr i} {
        lassign [lindex $XYZ [expr {$i-1}]] x y z
        sp $i 1 [expr { 0.1 + 0.01*$x + 0.02*$y}]
        sp $i 2 [expr {-0.2 - 0.03*$x + 0.01*$z}]
        sp $i 3 [expr { 0.3 + 0.02*$x + 0.01*$y - 0.02*$z}]
    }
}
staticAnalysis 1.0
analyze 1
set f [eleResponse 1 force]
set fmax 0.0
foreach v $f { set fmax [expr {max($fmax, abs($v))}] }
check "max |element force| under linear field" $fmax 0.0 1e-8
set jmax 0.0
foreach v [eleResponse 1 localDisplacement] { set jmax [expr {max($jmax, abs($v))}] }
check "max |local jump| under linear field" $jmax 0.0 1e-12

puts "\n=== T4: element added after initial displacement (initial state removal) ==="
wipe
model BasicBuilder -ndm 3 -ndf 3
node 1 0 0 0; node 2 1 0 0; node 3 0 1 0; node 4 0 0 1
node 5 0.2 0.2 0.2
fix 1 1 1 1; fix 2 1 1 1; fix 3 1 1 1; fix 4 1 1 1
addDummy 3
timeSeries Linear 1
pattern Plain 1 1 { sp 5 1 0.0; sp 5 2 0.0; sp 5 3 -0.05 }
staticAnalysis 1.0
analyze 1
loadConst -time 0.0
remove loadPattern 1
wipeAnalysis
element EmbeddedNodeContact 1 5 1 2 3 4 1.0e4 1.0e4 0.5 -orient 0 0 1
timeSeries Path 2 -time {0 1 100} -values {-0.05 -0.06 -0.06}
pattern Plain 2 2 { sp 5 3 1.0 }
timeSeries Constant 3
pattern Plain 3 3 { sp 5 1 0.0; sp 5 2 0.0 }
staticAnalysis 0.5
analyze 2
check "only the new 0.01 penetration counts: N = -100" [eleResponse 1 normalContactForce] -100.0 1e-6

puts "\n=== T5: mixed DOFs (embedded node ndf 6, host ndf 3) ==="
wipe
model BasicBuilder -ndm 3 -ndf 3
node 1 0 0 0; node 2 1 0 0; node 3 0 1 0; node 4 0 0 1
node 5 0.2 0.2 0.2 -ndf 6
fix 1 1 1 1; fix 2 1 1 1; fix 3 1 1 1; fix 4 1 1 1
fix 5 0 0 0 1 1 1
element EmbeddedNodeContact 1 5 1 2 3 4 1.0e4 1.0e4 0.5 -orient 0 0 1
addDummy 3
timeSeries Linear 1
pattern Plain 1 1 { sp 5 1 0.001; sp 5 2 0.0; sp 5 3 -0.01 }
staticAnalysis 1.0
analyze 1
check "ndf6: N = -100" [eleResponse 1 normalContactForce] -100.0 1e-6
check "ndf6: T2 = -10 (stick, t2 = -x)" [lindex [eleResponse 1 localForce] 2] -10.0 1e-6
check "ndf6: force vector size = 6+4*3" [llength [eleResponse 1 force]] 18 0

# ---------------------------------------------------------------------------
# T6: 2D, three continuum elements: two tri31 soil elements + one quad block.
# The two bottom nodes of the block are embedded in the top soil triangle.
#   W = 100 (gravity at top nodes), mu = 0.4 -> sliding base shear = 40
# ---------------------------------------------------------------------------
proc buildBlock2D {intType} {
    wipe
    model BasicBuilder -ndm 2 -ndf 2
    nDMaterial ElasticIsotropic 1 1.0e5 0.3
    nDMaterial ElasticIsotropic 2 1.0e7 0.2
    # soil
    node 1 0.0 0.0; node 2 2.0 0.0; node 3 2.0 1.0; node 4 0.0 1.0
    fix 1 1 1; fix 2 1 1
    element tri31 1 1 2 3 1.0 PlaneStrain 1
    element tri31 2 1 3 4 1.0 PlaneStrain 1
    # block (separate nodes, coincident with the soil top edge)
    node 5 0.5 1.0; node 6 1.5 1.0; node 7 1.5 2.0; node 8 0.5 2.0
    element quad 3 5 6 7 8 1.0 PlaneStrain 2
    # contact (block bottom nodes embedded in soil triangle 1-3-4)
    element EmbeddedNodeContact 11 5 1 3 4 1.0e5 1.0e5 0.4 -orient 0 1 -intType $intType
    element EmbeddedNodeContact 12 6 1 3 4 1.0e5 1.0e5 0.4 -orient 0 1 -intType $intType
    # gravity
    timeSeries Linear 1
    pattern Plain 1 1 { load 7 0.0 -50.0; load 8 0.0 -50.0 }
    constraints Plain
    numberer Plain
    system FullGeneral
    test NormDispIncr 1.0e-10 50 0
    algorithm Newton
    integrator LoadControl 0.1
    analysis Static
    analyze 10
    loadConst -time 0.0
}

foreach intType {0 1} {
    puts "\n=== T6: 2D tri31 + tri31 + quad, gravity then lateral push (intType $intType) ==="
    buildBlock2D $intType
    set Nsum [expr {[eleResponse 11 normalContactForce] + [eleResponse 12 normalContactForce]}]
    check "gravity: sum of normal contact forces = -W" $Nsum -100.0 1e-6
    # lateral push: prescribed ux = 0.05 at the block top nodes
    # (a perfectly sliding block has zero lateral stiffness, so force or
    #  displacement control on a load factor would be singular)
    timeSeries Linear 2
    pattern Plain 2 2 { sp 7 1 0.05; sp 8 1 0.05 }
    wipeAnalysis
    constraints Transformation; numberer Plain; system FullGeneral
    test NormDispIncr 1.0e-10 50 0
    algorithm Newton
    # note: the implicit ASDimplex law (also in zeroLengthContactASDimplex) needs
    # small steps at the stick/slip transition; IMPL-EX converges with large steps
    set dtp [expr {$intType == 0 ? 0.002 : 0.02}]
    set np [expr {int(round(1.0/$dtp))}]
    integrator LoadControl $dtp
    analysis Static
    set ok 0
    for {set i 0} {$i < $np && $ok == 0} {incr i} { set ok [analyze 1] }
    check "push converged ($np steps)" $ok 0 0
    reactions
    set H [expr {[nodeReaction 7 1] + [nodeReaction 8 1]}]
    set tol [expr {$intType == 0 ? 1e-6 : 0.4}]
    check "sliding base shear = mu*W = 40" $H 40.0 $tol
    # IMPL-EX: equilibrium holds for the extrapolated (localForceImplex) forces
    set resp [expr {$intType == 0 ? "localForce" : "localForceImplex"}]
    set Tsum [expr {[lindex [eleResponse 11 $resp] 1] + [lindex [eleResponse 12 $resp] 1]}]
    # localForce is the element resisting force in the contact frame (t1 = +x here),
    # so it balances the applied push H
    check "sum of contact tangential forces = H" $Tsum $H [expr {1e-6*abs($H)}]
    check "both contacts slipping" [expr {[eleResponse 11 contactState] + [eleResponse 12 contactState]}] 4 0
    check "block slid (base node ux > 0.03)" [expr {[nodeDisp 5 1] > 0.03}] 1 0
}

# ---------------------------------------------------------------------------`r`n# T7: EmbeddedNodeContact with the embedded node on a host vertex (N = [1,0,0])
# must be identical to zeroLengthContactASDimplex between the two nodes.
# (2D: the parent element needs 3 values after -orient)
# ---------------------------------------------------------------------------
puts "\n=== T7: equivalence with zeroLengthContactASDimplex (embedded node on a host vertex) ==="
proc buildEquiv {kind intType} {
    wipe
    model BasicBuilder -ndm 2 -ndf 2
    nDMaterial ElasticIsotropic 1 1.0e5 0.3
    nDMaterial ElasticIsotropic 2 1.0e7 0.2
    # soil: top nodes at x = 0, 0.5, 1.5, 2
    node 1 0.0 0.0; node 2 2.0 0.0
    node 3 0.0 1.0; node 4 0.5 1.0; node 9 1.5 1.0; node 10 2.0 1.0
    fix 1 1 1; fix 2 1 1
    element tri31 1 1 2 10 1.0 PlaneStrain 1
    element tri31 2 1 10 9 1.0 PlaneStrain 1
    element tri31 3 1 9 4 1.0 PlaneStrain 1
    element tri31 4 1 4 3 1.0 PlaneStrain 1
    node 5 0.5 1.0; node 6 1.5 1.0; node 7 1.5 2.0; node 8 0.5 2.0
    element quad 5 5 6 7 8 1.0 PlaneStrain 2
    if {$kind == "ENC"} {
        element EmbeddedNodeContact 11 5 4 3 1 1.0e5 1.0e5 0.4 -orient 0 1 -intType $intType
        element EmbeddedNodeContact 12 6 9 4 1 1.0e5 1.0e5 0.4 -orient 0 1 -intType $intType
    } else {
        element zeroLengthContactASDimplex 11 4 5 1.0e5 1.0e5 0.4 -orient 0 1 0 -intType $intType
        element zeroLengthContactASDimplex 12 9 6 1.0e5 1.0e5 0.4 -orient 0 1 0 -intType $intType
    }
    timeSeries Linear 1
    pattern Plain 1 1 { load 7 0.0 -50.0; load 8 0.0 -50.0 }
    constraints Plain; numberer Plain; system FullGeneral
    test NormDispIncr 1.0e-10 50 0; algorithm Newton
    integrator LoadControl 0.1; analysis Static
    analyze 10
    loadConst -time 0.0
    timeSeries Linear 2
    pattern Plain 2 2 { sp 7 1 0.05; sp 8 1 0.05 }
    wipeAnalysis
    constraints Transformation; numberer Plain; system FullGeneral
    test NormDispIncr 1.0e-10 50 0; algorithm Newton
}

foreach intType {0 1} {
  foreach dt [expr {$intType == 0 ? {0.002} : {0.02}}] {
    set hist {}
    foreach kind {ENC ZL} {
        buildEquiv $kind $intType
        integrator LoadControl $dt
        analysis Static
        set n [expr {int(round(1.0/$dt))}]
        set ok 0; set i 0; set H {}
        for {} {$i < $n && $ok == 0} {incr i} {
            set ok [analyze 1]
            reactions
            lappend H [expr {[nodeReaction 7 1] + [nodeReaction 8 1]}]
        }
        check "$kind push converged (intType $intType)" $ok 0 0
        lappend hist $H
    }
    set a [lindex $hist 0]; set b [lindex $hist 1]
    set m [expr {min([llength $a],[llength $b])}]
    set dmax 0.0
    for {set k 0} {$k < $m} {incr k} { set dmax [expr {max($dmax, abs([lindex $a $k]-[lindex $b $k]))}] }
    check "same base shear history as zeroLengthContactASDimplex (intType $intType)" $dmax 0.0 1e-8
  }
}

puts "\n=== SUMMARY: $::nPass passed, $::nFail failed ==="
wipe
if {$::nFail > 0} { exit 1 }
exit 0
