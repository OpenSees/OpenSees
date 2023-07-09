proc printEigenvalues {E A I L} {

    # Compute element eigenvalues
    set eigenValues [eigen standard 6]

    set e1 [expr 2*$E*$A/$L]
    set e2 [expr 2*$E*$I/$L]
    set e3 [expr 6*$E*$I/$L + 2*12*$E*$I/($L*$L*$L)]

    set exact "$e1 $e2 $e3"
    set exact [lsort -real $exact]

    puts "Eigenvalues"
    puts "Computed      Exact"
    for {set i 0} {$i < 3} {incr i} {
	puts "[lindex $eigenValues [expr $i+3]]  [lindex $exact $i]"
    }
}
proc printDisplacements {E A I L P H M} {

    puts "Nodal Displacements"
    puts "Computed      Exact"
    puts "[nodeDisp 2 1]  [expr $P*$L/($E*$A)]"
    puts "[nodeDisp 2 2]  [expr $H*pow($L,3)/(3*$E*$I) + $M*pow($L,2)/(2*$E*$I)]"
    puts "[nodeDisp 2 3]  [expr $H*pow($L,2)/(2*$E*$I) + $M*$L/($E*$I)]"
    puts ""
}

set E 1.0
set b 2.0
set d 12.0

set A [expr $b*$d]
set I [expr 1.0/12*$b*$d*$d*$d]

set L  100.0
set beta 0.1
set lp [expr $beta*$L]

set nIP 3

set elements {1 2 3 4}
set sections {1 2 3 4 5 6 7}

foreach element $elements {
    foreach section $sections {
	
	model basic -ndm 2 -ndf 3
	
	node 1 0.0 0.0
	node 2  $L 0.0
	
	switch $section {
	    1 {
		puts "Section: Elastic"
		section Elastic 1 $E $A $I
	    }
	    2 {
		puts "Section: Fiber"
		uniaxialMaterial Elastic 1 $E
		section Fiber 1 {
		    set b2 [expr $b/2]
		    set d2 [expr $d/2]
		    patch rect 1 10 1 $d2 $b2 -$d2 -$b2
		}
	    }
	    3 {
		puts "Section: Aggregator (two uniaxial)"
		uniaxialMaterial Elastic 1 [expr $E*$A]
		uniaxialMaterial Elastic 2 [expr $E*$I]
		section Aggregator 1 1 P 2 Mz
	    }
	    4 {
		puts "Section: Aggregator (uniaxial plus section)"
		uniaxialMaterial Elastic 1 [expr $E*$A]
		uniaxialMaterial Elastic 2 [expr $E*$I]
		section Uniaxial 2 1 P
		section Aggregator 1 2 Mz -section 2
	    }
            5 {
                puts "Section: YieldSurface01 (Orbison)"
                ysEvolutionModel null 100 1.0 1.0
                yieldSurface_BC Orbison2D 100 1.0e12 1.0e12 100
                section YS_Section2D01 1 $E $A $I 100
            }
            6 {
                puts "Section: YieldSurface01 (Attalla)"
                ysEvolutionModel null 100 1.0 1.0
                yieldSurface_BC Attalla2D 100 1.0e12 1.0e12 100
                section YS_Section2D01 1 $E $A $I 100
            }
            7 {
                puts "Section: YieldSurface01 (ElTawil)"
                ysEvolutionModel null 100 1.0 1.0
                yieldSurface_BC ElTawil2D 100 1.0e12 0.0 1.0e12 -1.0e12 100
                section YS_Section2D01 1 $E $A $I 100
            }
	}

	geomTransf Linear 1

	switch $element {
	    1 {
		puts "Element: ElasticBeamColumn"
		element elasticBeamColumn 1 1 2 $A $E $I 1
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
		element beamWithHinges 1 1 2 1 $lp 1 $lp $E $A $I 1
	    }
	}
	
	printEigenvalues $E $A $I $L

	fix 1 1 1 1

	set M 10.0      ;# End moment
	set H 1.0       ;# Transverse load
	set P -100.0    ;# Axial load
	
	pattern Plain 1 Linear {
	    load 2 $P $H $M
	}
	
	test NormUnbalance 1.0e-12 10
	algorithm Newton
	integrator LoadControl 1
	constraints Plain
	system ProfileSPD
	numberer Plain
	analysis Static
	
	analyze 1

	printDisplacements $E $A $I $L $P $H $M

	wipe
    }
}
