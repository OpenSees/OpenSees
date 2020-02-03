# Linear Elastic Planar Shear Walls

#REFERENCES: 
# 1) ETABS Software Verification Examples, Computers and Structures, Inc, 2003 (Example 15A)

puts "PlanarShearWall.tcl: Verification of Linear Elastic Planar Shear Wall"
puts "   NOTE: using SAP2000 results for verification"

# Multiple Shear Wall building models with lengths of 120",240" and 360" lengths. Buildings of 1 story, 3 Story and 
# 6 story are modelled. Each buildings story height is 120" All walls 12". All materials elastic with modulus of elasticity 
# of 3000 ksi and poisson's ration of 0.2. At each floor in the buildings al nodes at the floor level are constrained to
# move horizontally together.
#
# Loading: For each building a 100k load is applied at top left node.
# Results: compare nodal displacement at node were load is applied to Etabs and SAP, results verified using SAP.

# NOTE: The discretization of the SAP and ETABS models are not known at this time

set testOK 0
set tol 5.0e-3

set resultsSAP   {2.4287 0.1031 0.0186 0.3205 0.0187 0.0052 0.0185 0.0029 0.0013}
set resultsETABS {2.3926 0.0985 0.0172 0.3068 0.0169 0.0046 0.0144 0.0024 0.0011}


# wall properties
set E 3000.
set v 0.2
set t 12.0

set floorHeight 120.
set bayWidth 120.

# modeling parameters # of elements in floorHeight * bayWidth block

#set nxBay   16
#set nyFloor 16
set nxBay   16
set nyFloor 16

# ----------------------------
# Start of model generation
# ----------------------------

# QUAD

#foreach numFloor {6 3 1} {


foreach eleType {quad SSPquad} {

    set counter 0
    puts "\n - using $eleType elements"

    set formatString {%12s%12s%12s%12s%12s%12s%12s}
    puts [format $formatString "# Stories" "Wall Height" "Wall Length" "ETABS" "SAP2000" "OpenSees" Difference]
    set formatString {%12.0f%12.0f%12.0f%12.4f%12.4f%12.4f%12.4f}

    foreach numFloor {6 3 1} {
	foreach numBay {1 3 6} {
	    wipe
	    
	    model basic -ndm 2 -ndf 2
	    
	    nDMaterial ElasticIsotropic 1  $E $v
	    
	    # set some parameters for node and element generation
	    set eleType quad
	    
	    set eleArgs "{$t PlaneStress 1}"
	    
	    set nx [expr $numBay * $nxBay];     # number of elements along building length
	    set ny [expr $numFloor * $nyFloor]; # number of elements along building height
	    
	    set nodeTop [expr ($nx+1)*$ny + 1]
	    
	    # generate the nodes and elements
	    set blockCmd "block2D $nx $ny 1 1 $eleType $eleArgs {
                          1 0. 0.
                          2 [expr $bayWidth * $numBay] 0.
                          3 [expr $bayWidth * $numBay] [expr $floorHeight * $numFloor]
                          4 0. [expr $floorHeight * $numFloor]
            }"
	    eval $blockCmd
	    
	    # add some loads
	    pattern Plain 1 Linear {
		load $nodeTop  100.0  0.0  
	    }
	    
	    fixY 0.0   1 1 
	    
	    #floor constraints
	    for {set i 1} {$i <= $numFloor} {incr i 1} {
		set mNode [expr ($nx+1)*$nyFloor*$i + 1]
		for {set j 1} {$j <= $nx} {incr j 1} {
		    equalDOF $mNode [expr $mNode+$j] 1
		}
	    }
	    
	    integrator LoadControl  1.0  
	    algorithm Linear
	    numberer RCM
	    constraints Plain 
	    system Umfpack
	    analysis Static 
	    
	    analyze 1
	    
	    set disp [nodeDisp $nodeTop 1]
	    set dispETABS [lindex $resultsETABS $counter]
	    set dispSAP   [lindex $resultsSAP $counter]
	    set diffR  [expr abs($dispSAP-$disp)]
	    puts [format $formatString $numFloor [expr $floorHeight * $numFloor] [expr $bayWidth * $numBay] $dispETABS $dispSAP $disp $diffR]
	    
	    
	    # verify result
	    if {[expr abs($disp-$dispSAP)] > $tol} {
		set testOK -1;
		puts "failed  quad: $disp - $dispSAP [expr abs($disp-$dispSAP)] > $tol"
	    }
	    
	    incr counter
	}
    }
}

# Shell
puts "\n - using Shell elements"

set formatString {%12s%12s%12s%12s%12s%12s%12s}
puts [format $formatString "# Stories" "Wall Height" "Wall Length" "ETABS" "SAP2000" "OpenSees" Difference]
set formatString {%12.0f%12.0f%12.0f%12.4f%12.4f%12.4f%12.4f}

set counter 0
foreach numFloor {6 3 1} {
    foreach numBay { 1 3 6} {
	wipe

	model basic -ndm 3 -ndf 6

	# create the material
#	set t 12.0
#	nDMaterial ElasticIsotropic 1 $E $v
#	nDMaterial PlateFiber 4 1
#	section PlateFiber 1 4 $t

	section ElasticMembranePlateSection  1 $E $v $t 0.0
	
	# set some parameters for node and element generation
	set Plate shell
	
	set eleArgs "1"

	set nx [expr $numBay * $nxBay];     # number of elements along building length
	set ny [expr $numFloor * $nyFloor]; # number of elements along building height	

	set nodeTop [expr ($nx+1)*$ny + 1]
	
	# generate the nodes and elements
	set blockCmd "block2D $nx $ny 1 1 $Plate $eleArgs {
                          1             0.                                0.              0.
                          2 [expr $bayWidth * $numBay]                    0.              0.
                          3 [expr $bayWidth * $numBay]   [expr $floorHeight * $numFloor]  0.
                          4             0.               [expr $floorHeight * $numFloor]  0.   
         }"

	eval $blockCmd
	
	# add some loads
	pattern Plain 1 Linear {
	    load $nodeTop  100.0  0.0  0.0   0.0   0.0  0.0
	}
	
	fixY 0.0   1 1 1 0 0 0


	#floor constraints
	for {set i 1} {$i <= $numFloor} {incr i 1} {
	    set mNode [expr ($nx+1)*$nyFloor*$i + 1]
	    for {set j 1} {$j <= $nx} {incr j 1} {
		equalDOF $mNode [expr $mNode+$j] 1
	    }
	}
	
	integrator LoadControl  1.0  
	algorithm Linear
	numberer RCM
	constraints Plain
	system UmfPack
	analysis Static 
	
	analyze 1
	
	set disp [nodeDisp $nodeTop 1]

	set dispETABS [lindex $resultsETABS $counter]
	set dispSAP   [lindex $resultsSAP $counter]
	set diffR  [expr abs($dispSAP-$disp)]
	puts [format $formatString $numFloor [expr $floorHeight * $numFloor] [expr $bayWidth * $numBay] $dispETABS $dispSAP $disp $diffR]

	# verify result
	if {[expr abs($disp-$dispSAP)] > $tol} {
	    set testOK -1;
	    puts "failed  shell: $disp - $dispSAP [expr abs($disp-$dispSAP)] > $tol"
	}

	incr counter	
    }
}


# Brick


foreach eleType {stdBrick SSPbrick} {
    set counter 0

    puts "\n - using $eleType elements"

    set formatString {%12s%12s%12s%12s%12s%12s%12s}
    puts [format $formatString "#Stories" "Wall Height" "Wall Length" "ETABS" "SAP2000" "OpenSees" "Difference"]
    set formatString {%12.0f%12.0f%12.0f%12.4f%12.4f%12.4f%12.4f}


    foreach numFloor {6 3 1} {
	foreach numBay { 1 3 6} {
	    wipe
	    
	    model basic -ndm 3 -ndf 3
	    
	    # create the material
	    
	    nDMaterial ElasticIsotropic 1 $E $v
	    
	    # set some parameters for node and element generation
	    set eleArgs "1"
	    
	    set nx [expr $numBay * $nxBay];     # number of elements along building length
	    set ny [expr $numFloor * $nyFloor]; # number of elements along building height	
	    set nz 1
	    
	    set nodeTop [expr ($nx+1)*$ny + 1]
	    
	    # generate the nodes and elements
	    set blockCmd "block3D $nx $ny $nz 1 1 SSPbrick $eleArgs {
                          1             0.                                0.              0.
                          2 [expr $bayWidth * $numBay]                    0.              0.
                          3 [expr $bayWidth * $numBay]   [expr $floorHeight * $numFloor]  0.
                          4             0.               [expr $floorHeight * $numFloor]  0.   
                          5             0.                                0.              $t
                          6 [expr $bayWidth * $numBay]                    0.              $t
                          7 [expr $bayWidth * $numBay]   [expr $floorHeight * $numFloor]  $t
                          8             0.               [expr $floorHeight * $numFloor]  $t   
             }"
	    
	    eval $blockCmd
	    
	    # add some loads
	    pattern Plain 1 Linear {
		load $nodeTop  50.0  0.0  0.0 
		load [expr $nodeTop + ($nx+1)*($ny+1)]  50.0  0.0  0.0 
	    }
	    
	    fixY 0.0   1 1 1 
	    
	    
	    #floor constraints
	    for {set i 1} {$i <= $numFloor} {incr i 1} {
		set mNode1 [expr ($nx+1)*$nyFloor*$i + 1]
		set mNode2 [expr $mNode1 + ($nx+1)*($ny+1)]
		for {set j 1} {$j <= $nx} {incr j 1} {
		    equalDOF $mNode1 [expr $mNode1+$j] 1
		    equalDOF $mNode2 [expr $mNode2+$j] 1
		}
	    }
	    
	    integrator LoadControl  1.0  
	    algorithm Linear
	    numberer RCM
	    constraints Plain 
	    system Umfpack
	    analysis Static 
	    
	    analyze 1
	    
	    set disp [nodeDisp $nodeTop 1]
	    
	    set dispETABS [lindex $resultsETABS $counter]
	    set dispSAP   [lindex $resultsSAP $counter]
	    set diffR  [expr abs($dispSAP-$disp)]
	    puts [format $formatString $numFloor [expr $floorHeight * $numFloor] [expr $bayWidth * $numBay] $dispETABS $dispSAP $disp $diffR]
	    
	    
	    # verify result
	    if {[expr abs($disp-$dispSAP)] > $tol} {
		set testOK -1;
		puts "failed  brick: $disp - $dispSAP [expr abs($disp-$dispSAP)] > $tol"
	    }
	    
	    incr counter	
	}
    }
}


set results [open results.out a+]
if {$testOK == 0} {
    puts "\nPASSED Verification Test PlanarShearWall.tcl \n\n"
    puts $results "PASSED : PlanarShearWall.tcl"
} else {
    puts "\nFAILED Verification Test PlanarShearWall.tcl \n\n"
    puts $results "FAILED : PlanarShearWall.tcl"
}
close $results