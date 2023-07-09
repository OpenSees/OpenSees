set E  30000.0
set A  20.0
set Iz 1400.0
set Iy 1200.0
set G  18000.0
set J  [expr $Iz+$Iy]

set L  100.0
set beta 0.2
set lp [expr $beta*$L]

# Point loads
set Py 200.0
set Pz 250.0
set N  150.0
set aL 0.6
set a  [expr $aL*$L]

# Distributed loads
set wy 3.0
set wz 4.0
set wx 2.0

set nIP 3

set elements {1 2 3 4}

foreach element $elements {
	
    model basic -ndm 3 -ndf 6
	
    node 1 0.0 0.0 0.0
    node 2  $L 0.0 0.0

    fix 1 1 1 1 1 0 0
    fix 2 0 1 1 0 0 0  

    section Elastic 1 $E $A $Iz $Iy $G $J
    
    geomTransf Linear 1 0 0 1
    
    switch $element {
	1 {
	    puts "Element: ElasticBeamColumn"
	    element elasticBeamColumn 1 1 2 $A $E $G $J $Iy $Iz 1
	}
	2 {
	    puts "Element: DispBeamColumn"
	    element dispBeamColumn 1 1 2 $nIP 1 1
	}
	3 {
	    puts "Element: NonlinearBeamColumn"
	    element nonlinearBeamColumn 1 1 2 $nIP 1 1
	}
	4 {
	    puts "Element: BeamWithHinges"
	    element beamWithHinges 1 1 2 1 $lp 1 $lp $E $A $Iz $Iy $G $J 1
	}
    }
    
    pattern Plain 1 "Constant" {
	set Py2 [expr $Py/2]
	set Pz2 [expr $Pz/2]
	set N2  [expr $N/2]
	
	eleLoad -ele 1 -type -beamPoint $Py2 $Pz2 $aL $N2
	eleLoad -ele 1 -type -beamPoint $Py2 $Pz2 $aL $N2

	set wy2 [expr $wy/2]
	set wz2 [expr $wz/2]
	set wx2 [expr $wx/2]
	
	eleLoad -ele 1 -type -beamUniform $wy2 $wz2 $wx2
	eleLoad -ele 1 -type -beamUniform $wy2 $wz2 $wx2
    }
    
    test NormUnbalance 1.0e-10 10 1
    algorithm Newton
    integrator LoadControl 1.0
    constraints Plain
    system ProfileSPD
    numberer Plain
    analysis Static
    
    analyze 2
    
    print node 1 2

    set d1 [expr $aL*$N*$L/($E*$A)]
    set d2 [expr $wx*pow($L,2)/(2*$E*$A)]
    puts "Exact axial displacement: [expr $d1+$d2]"

    set Vy [expr $Py*(1-$aL)]
    set Vz [expr $Pz*(1-$aL)]

    set d1 [expr -$Vy/(6*$E*$Iz*$L)*($a*$a*$L-2*$a*$L*$L)]
    set d2 [expr $wy*pow($L,3)/(24*$E*$Iz)]
    puts "Exact rotation (node 1): [expr $d1+$d2]"

    set d1 [expr -$Vy/(6*$E*$Iz*$L)*($a*$a*$L+$a*$L*$L)]
    puts "Exact rotation (node 2): [expr $d1-$d2]"

    set d1 [expr -$Vz/(6*$E*$Iy*$L)*($a*$a*$L-2*$a*$L*$L)]
    set d2 [expr $wz*pow($L,3)/(24*$E*$Iy)]
    puts "Exact rotation (node 1): [expr $d1+$d2]"

    set d1 [expr -$Vz/(6*$E*$Iy*$L)*($a*$a*$L+$a*$L*$L)]
    puts "Exact rotation (node 2): [expr $d1-$d2]"
    puts ""

    wipe
}
