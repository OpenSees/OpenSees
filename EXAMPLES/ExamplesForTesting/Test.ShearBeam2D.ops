set E 30000.0
set A 20.0
set I 1400.0
set G 18000.0
set alpha [expr 5.0/6.0]

set L  100.0
set beta 0.2
set lp [expr $beta*$L]

set P 200.0

set nIP 3

set elements {1 2}

foreach element $elements {
	
    model basic -ndm 2 -ndf 3
	
    node 1 0.0 0.0
    node 2  $L 0.0
    
    fix 1 1 1 1

    section Elastic 1 $E $A $I

    uniaxialMaterial Elastic 1 [expr $alpha*$G*$A]
    section Uniaxial 2 1 Vy

    section Aggregator 3 1 Vy -section 1

    geomTransf Linear 1
    
    switch $element {
	1 {
	    puts "Element: NonlinearBeamColumn"
	    element nonlinearBeamColumn 1 1 2 $nIP 3 1
	}
	2 {
	    puts "Element: BeamWithHinges"
	    element beamWithHinges 1 1 2 3 $lp 1 $lp $E $A $I 1
	}
	3 {
	    puts "Element: BeamWithHinges2"
	    element beamWithHinges2 1 1 2 1 $lp 1 $lp $E $A $I 1 -constHinge 2
	}
    }
    
    pattern Plain 1 "Constant" {
	load 2 0.0 $P 0.0
    }
    
    test NormUnbalance 1.0e-10 10 1
    algorithm Newton
    integrator LoadControl 1.0
    constraints Plain
    system ProfileSPD
    numberer Plain
    analysis Static
    
    analyze 1
    
    print node 2

    puts "Exact displacement:   [expr $P*pow($L,3)/(3*$E*$I) + $P*$L/($alpha*$G*$A)]"
    puts "Bending contribution: [expr $P*pow($L,3)/(3*$E*$I)]"
    puts "Shear contribution:   [expr $P*$L/($alpha*$G*$A)]"
    puts ""

    wipe
}
