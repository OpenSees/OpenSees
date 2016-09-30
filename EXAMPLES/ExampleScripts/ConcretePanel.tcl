## Panel tests by Adebar/Collins (1989)
#
# nxm mesh with uniformly distributed tractions at ends
#
# Department of Civil and Environmental Engineering
# University of California, Berkeley


foreach sectType {"LayeredShell"} {

    foreach numEle {2} {

	wipe
	model basic -ndm 3 -ndf 6
	
	#panel dimensions
	set L 1524.0; # mm
	set d 1524.0; # mm
	set t  310.0; # mm
	
	for {set i 0; set nodeTag 1} {$i <= $numEle} {incr i 1} {
	    for {set j 0} {$j <= $numEle} {incr j 1} {
		node $nodeTag [expr $j*$d/$numEle] [expr $i*$L/$numEle] 0.
		incr nodeTag 1
	    }
	}
	
	fixY 0. 1 1 1 1 1 1
        
	# nDMaterial ElasticIsotropic   1   100   0.25  1.27
	nDMaterial PlasticDamageConcrete3d   1   31000 0.2 4. 26. 0.5 1.0 2.0 0.8
	uniaxialMaterial Steel01 1 21000. 485. [expr 2500./210000.]
	# nDMaterial Damage2p 1 15 -fct 3 -E 32000 -ni 0.2
	nDMaterial PlateRebar 101 1 45
	nDMaterial PlateRebar 102 1 -45
	nDMaterial PlateFromPlaneStress 104 1 0.3
	

#	nDMaterial ElasticIsotropic 200 20000 0.3
	
	if {$sectType == "PlateFiber"} {
	    # section PlateFiber $secTag $matTag $h
	    nDMaterial PlateFiber 4 1
	    section PlateFiber    7 4 $t
	} 

	if {$sectType == "LayeredShell"} {
	    # section LayeredShell $sectionTag $nLayers $matTag1 $thickness1...$matTagn $thicknessn
	    set numCLayers 4;
	    set tCover 50;
	    set tSteel [expr (.0358*$t)/4.0]; # 3.58% steel
	    set tLayer [expr ($t-4*$tSteel-2*$tCover)/$numCLayers];
	    set numLayers [expr 2 + 4 + $numCLayers];

	    set tLayer [expr $tLayer/($numLayers*1.)]
puts "section  LayeredShell $numLayers 104 $tCover 101 $tSteel 102 $tSteel 1 $tLayer 1 $tLayer 1 $tLayer 1 $tLayer 101 $tSteel 102 $tSteel 104 $tCover "

section  LayeredShell $numLayers $tCover 1 $tSteel 101 $tSteel 102 $tLayer 104 $tLayer 104 $tLayer 104 $tLayer 104 $tSteel 102 $tSteel 101 $tCover 1
	
	}	
	for {set i 0; set eleTag 1} {$i < $numEle} {incr i 1} {
	    set iNode [expr $i*($numEle+1)+1]
	    set lNode [expr $iNode+$numEle+1]
	    for {set j 0} {$j < $numEle} {incr j 1} {
		element ShellMITC4 $eleTag $iNode [expr $iNode+1] [expr $lNode+1] $lNode 7
		incr iNode 1
		incr lNode 1
		incr eleTag 1
	    }
	}
	
	set P [expr 10.0/($numEle+1.0)]
	
	timeSeries Linear 1
	pattern Plain 1 1 {
	    set iNode [expr ($numEle)*($numEle+1)+1]
	    for {set j 0} {$j <= $numEle} {incr j 1} {
		load $iNode $P 0 0 0. 0. 0.
		incr iNode 1
	    }
	}
	
	numberer Plain
	constraints Plain
	system ProfileSPD
	integrator LoadControl 1.
	test NormDispIncr 1.0e-6 4 0
	algorithm Newton
	analysis Static
	analyze 1
	
	set disp [nodeDisp [expr ($numEle+1)*($numEle+1)] 1]
	puts "$sectType numEle: $numEle Disp: $disp"
    }
}




