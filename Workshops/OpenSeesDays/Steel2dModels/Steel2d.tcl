# SteelWSection.tcl

# Written: fmk

# Major Contributions: 
#   WSection: Remo DeSouza, UC Berkeley.
#   SteelWSectionMR: Dimitrios Lignos, McGill University.
#   SteelWSectionMR02: Filipe Ribeiro and Andre Barbosa, Oregon State University.
#   SteelWSectionMChi02: Filipe Ribeiro and Andre Barbosa, Oregon State University.
#   HSSbrace: Dimitrios Lignos, McGill University.
#   

#
# 1. PROCEDURAL PROTOTYPES FOR CREATING ELEMENTS:
#

# ElasticBeamWSection2d $eleTag $iNode $jNode $sectType $E $transfTag $args
#   args: <YY> <-release1> <-release2>
#
# ForceBeamWSection2d $eleTag $iNode $jNode $sectType $matTag $transfTag $args 
#    args: <YY> <-nFlange $x> <-nWeb $x> <-nip $x> <-elasticSection $x> <-release1> <-release2>
#
# DispBeamWSection2d $eleTag $iNode $jNode $sectType $matTag $transfTag $args 
#    args: <YY> <-nFlange $x> <-nWeb $x> <-nip $x> <-elasticSection $E> <-release1> <-release2>
#
# BeamWithHingesWSection2d $eleTag $iNode $jNode $sectType $matTag $transfTag $args 
#    args: <YY> <-release1> <-release2> <-nFlange $x> <-nWeb $x> <-hingeLength $x>
#									
# BeamWithPlasticHingesWSection2d $eleTag $iNode $jNode $sectType $E $Fy $H $Lb $Com_Type $Comp_Action $Lp $transfTag $args
#    args: <YY> <-metric> <-release1> <-release2>
#
# ElasticBeamHSSection2d $eleTag $iNode $jNode $sectType $E $transfTag $args
#    args: <YY> <-release1> <-release2>
#
# HSSbrace $eleTag $iNode $jNode $sectType $matTag $numSeg $Im $transfTag $args
#    args: <YY> <-nip $x> <-elasticSection $x> 
#
# BeamWithConcentratedHingesWSection2d $eleTag $iNode $jNode $sectType $E $Fy $H $Lb $Com_Type $Comp_Action $nFactor $transfTag $args
#    args: <YY> <-metric> <-release1> <-release2> 

#
# 2. PROCEDURAL PROTOTYPES FOR CREATING SECTIONS:
#

# ElasticSteelWSection2d $sectTag $sectType $E $args
#    args: <YY>
#
# FiberSteelWSection2d $sectTag $sectType $matTag $nFlange $nWeb $args
#    args: <YY>
#
# Wsection $secID $matID $d $bf $tf $tw $nfdw $nftw $nfbf $nftf {$Orient XX} 
#
# ElasticBeamWSection2d $eleTag $iNode $jNode $sectType $E $transfTag $args 
#    args: <YY>
#
# Wsection $secID $matID $d $bf $tf $tw $nfdw $nftw $nfbf $nftf 
#    args: <YY>
#
# ElasticHSSection2d {sectTag sectType E args} {
#    args: <YY>
#
# FiberHSSection2d $sectTag $sectType $matTag $nFlange $nWeb $args
#    args: <YY>
#
# HSSectionD $secID $matID $d $b $t $nfdw $nftw
#
# SteelWSectionMR $matTag $E $Fy $Ix $Sx $H $L $d $tw $bf $tf $Lb $ry $Com_Type $Comp_Action $args
#    args: <-hLength $x> <-metric>
#
# SteelWSectionMR02 $matTag $E $Fy $H $L $Lb $sectType $Com_Type $Comp_Action $args
#    args: <-nFactor $x> <-metric>
#
# SteelWSectionMChi02 $matTag $E $Fy $H $L $Lb $sectType $Com_Type $Comp_Action $Lp $args
#    args: <-metric>

#
# 3. PROCEDURES FOR CREATING ELEMENTS:
#

proc ElasticBeamWSection2d {eleTag iNode jNode sectType E transfTag args} {
    global WSection
    global in
    set found 0

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
        set Orient "YY"
    }

    if {[lsearch $args "-release1"] != -1} {
	set hingeEnd1 1
	node $eleTag$hingeEnd1 [nodeCoord $iNode 1] [nodeCoord $iNode 2]
	equalDOF $iNode $eleTag$hingeEnd1 1 2
	set iNode $eleTag$hingeEnd1
    }

    if {[lsearch $args "-release2"] != -1} {
	set hingeEnd2 2
	node $eleTag$hingeEnd2 [nodeCoord $jNode 1] [nodeCoord $jNode 2]
	equalDOF $jNode $eleTag$hingeEnd2 1 2
	set jNode $eleTag$hingeEnd2
    }

    foreach {section prop} [array get WSection $sectType] {
	set propList [split $prop]
	#AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
	set A [expr [lindex $propList 0]*$in*$in]
	set Ixx [expr [lindex $propList 5]*$in*$in*$in*$in]
	set Iyy [expr [lindex $propList 6]*$in*$in*$in*$in]
	if {$Orient == "YY" } {
	    element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iyy $transfTag
	} else {
	    element elasticBeamColumn $eleTag $iNode $jNode $A $E $Ixx $transfTag
	}
	set found 1
    }

    if {$found == 0} {
	puts "ElasticBeamWSection2d sectType: $sectType not found for ee: $eleTag"
    }
}

proc ForceBeamWSection2d {eleTag iNode jNode sectType matTag transfTag args} {

    global FiberSteelWSection2d
    global ElasticSteelWSection2d

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
        set Orient "YY"
    }
   
    set nFlange 10
    if {[lsearch $args "-nFlange"] != -1} {
	set loc [lsearch $args "-nFlange"]
        set nFlange [lindex $args [expr $loc+1]]
    }

    set nWeb  5
    if {[lsearch $args "-nWeb"] != -1} {
	set loc [lsearch $args "-nWeb"]
        set nWeb [lindex $args [expr $loc+1]]
    }

    set nip 4
    if {[lsearch $args "-nip"] != -1} {
	set loc [lsearch $args "-nip"]
        set nip [lindex $args [expr $loc+1]]
    }

    if {[lsearch $args "-release1"] != -1} {
	set hingeEnd1 1
	node $eleTag$hingeEnd1 [nodeCoord $iNode 1] [nodeCoord $iNode 2]
	equalDOF $iNode $eleTag$hingeEnd1 1 2
	set iNode $eleTag$hingeEnd1
    }

    if {[lsearch $args "-release2"] != -1} {
	set hingeEnd2 2
	node $eleTag$hingeEnd2 [nodeCoord $jNode 1] [nodeCoord $jNode 2]
	equalDOF $jNode $eleTag$hingeEnd2 1 2
	set jNode $eleTag$hingeEnd2
    }

    if {[lsearch $args "-elasticSection"] != -1} {
	set loc [lsearch $args "-elasticSection"] 
        set E [lindex $args [expr $loc+1]]
	ElasticSteelWSection2d $eleTag $sectType $E  $Orient
    } else {
	FiberSteelWSection2d $eleTag $sectType $matTag $nFlange $nWeb $Orient
    }

    element forceBeamColumn $eleTag $iNode $jNode $nip $eleTag $transfTag
}

proc DispBeamWSection2d {eleTag iNode jNode sectType matTag transfTag args} {

    global FiberSteelWSection2d

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
        set Orient "YY"
    }

    set nFlange 10
    if {[lsearch $args "-nFlange"] != -1} {
	set loc [lsearch $args "-nFlange"]
        set nFlange [lindex $args [expr $loc+1]]
    }

    set nWeb  5
    if {[lsearch $args "-nWeb"] != -1} {
	set loc [lsearch $args "-nWeb"]
        set nWeb [lindex $args [expr $loc+1]]
    }

    set nip 4
    if {[lsearch $args "-nip"] != -1} {
	set loc [lsearch $args "-nip"]
        set nip [lindex $args [expr $loc+1]]
    }

    set intType "Lobatto"
    if {[lsearch $args "-int"] != -1} {
	set loc [lsearch $args "-int"]
        set intType [lindex $args [expr $loc+1]]
    }    

    if {[lsearch $args "-release1"] != -1} {
	set hingeEnd1 1
	node $eleTag$hingeEnd1 [nodeCoord $iNode 1] [nodeCoord $iNode 2]
	equalDOF $iNode $eleTag$hingeEnd1 1 2
	set iNode $eleTag$hingeEnd1
    }

    if {[lsearch $args "-release2"] != -1} {
	set hingeEnd2 2
	node $eleTag$hingeEnd2 [nodeCoord $jNode 1] [nodeCoord $jNode 2]
	equalDOF $jNode $eleTag$hingeEnd2 1 2
	set jNode $eleTag$hingeEnd2
    }

    set eleType dispBeamColumn

    if {[lsearch $args "-elasticSection"] != -1} {
	set loc [lsearch $args "-elasticSection"] 
        set E [lindex $args [expr $loc+1]]
	ElasticSteelWSection2d $eleTag $sectType $E  $Orient
    } else {
	FiberSteelWSection2d $eleTag $sectType $matTag $nFlange $nWeb $Orient
    }

    element $eleType $eleTag $iNode $jNode $nip $eleTag $transfTag -integration $intType
}

proc BeamWithHingesWSection2d {eleTag iNode jNode sectType matTag transfTag args} {

    global FiberSteelWSection2d
    global WSection
    global in

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
        set Orient "YY"
    }

    set nFlange 10
    if {[lsearch $args "-nFlange"] != -1} {
	set loc [lsearch $args "-nFlange"]
        set nFlange [lindex $args [expr $loc+1]]
    }

    set nWeb    10
    if {[lsearch $args "-nWeb"] != -1} {
	set loc [lsearch $args "-nWeb"]
        set nWeb [lindex $args [expr $loc+1]]
    }

    if {[lsearch $args "-release1"] != -1} {
	set hingeEnd1 1
	node $eleTag$hingeEnd1 [nodeCoord $iNode 1] [nodeCoord $iNode 2]
	equalDOF $iNode $eleTag$hingeEnd1 1 2
	set iNode $eleTag$hingeEnd1
    }

    if {[lsearch $args "-release2"] != -1} {
	set hingeEnd2 2
	node $eleTag$hingeEnd2 [nodeCoord $jNode 1] [nodeCoord $jNode 2]
	equalDOF $jNode $eleTag$hingeEnd2 1 2
	set jNode $eleTag$hingeEnd2
    }

    set found 0
    set d 0
    foreach {section prop} [array get WSection $sectType] {
	set propList [split $prop]
	set  d [expr [lindex $propList 1]*$in]
	set found 1
    }
    set Lp $d

    if {[lsearch $args "-hingeLength"] != -1} {
	set loc [lsearch $args "-hingeLength"]
        set hingeLength [lindex $args [expr $loc+1]]
    } 

    FiberSteelWSection2d $eleTag $sectType $matTag $nFlange $nWeb
    element forceBeamColumn $eleTag $iNode $jNode $transfTag "HingeRadau $eleTag $Lp $eleTag $Lp $eleTag"
}

proc BeamWithSteel01HingesWSection2d {eleTag iNode jNode sectType E Fy b transfTag args} {

    global FiberSteelWSection2d
    global WSection
    global in

    set n 10
    if {[lsearch $args "n"] != -1} {
	set loc [lsearch $args "n"]
        set n [lindex $args [expr $loc+1]]
    }

    set doRayleigh 0
    if {[lsearch $args "-doRayleigh"] != -1} {
	set loc [lsearch $args "-doRayleigh"]
        set doRayleigh [lindex $args [expr $loc+1]]
    }

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
        set Orient "YY"
    }

    set nFlange 10
    if {[lsearch $args "-nFlange"] != -1} {
	set loc [lsearch $args "-nFlange"]
        set nFlange [lindex $args [expr $loc+1]]
    }

    set nWeb    10
    if {[lsearch $args "-nWeb"] != -1} {
	set loc [lsearch $args "-nWeb"]
        set nWeb [lindex $args [expr $loc+1]]
    }

    set found 0
    set d 0

    foreach {section prop} [array get WSection $sectType] {
	set propList [split $prop]

	#AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
	set A [expr [lindex $propList 0]*$in*$in]
	set  d [expr [lindex $propList 1]*$in]
	set  bf [expr [lindex $propList 2]*$in]
	set  tw [expr [lindex $propList 3]*$in]
	set  tf [expr [lindex $propList 4]*$in]
	set Ixx [expr [lindex $propList 5]*$in*$in*$in*$in]
	set Iyy [expr [lindex $propList 6]*$in*$in*$in*$in]
	set Zx [expr [lindex $propList 7]*$in*$in*$in]
	set found 1
    }

    set dX [expr [nodeCoord $jNode 1] - [nodeCoord $iNode 1]]
    set dY [expr [nodeCoord $jNode 2] - [nodeCoord $iNode 2]]
    set L [expr sqrt($dX*$dX+$dY*$dY)]
    set Lp $d


    # create 2 additional nodes at either end, constrain to move in 1 and 2 with nodes at end
    set hingeEnd1 1
    node $eleTag$hingeEnd1 [nodeCoord $iNode 1] [nodeCoord $iNode 2]
    equalDOF $iNode $eleTag$hingeEnd1 1 2

    set hingeEnd2 2
    node $eleTag$hingeEnd2 [nodeCoord $jNode 1] [nodeCoord $jNode 2]
    equalDOF $jNode $eleTag$hingeEnd2 1 2

    #    set iNodeCrds [nodeCoord $iNode]; set x1 [lindex $iNodeCrds 0]; set y1 [lindex $iNodeCrds 1]
    #    set jNodeCrds [nodeCoord $jNode]; set x2 [lindex $jNodeCrds 0]; set y2 [lindex $jNodeCrds 1]
    #    set L [expr sqrt(($x2-$x1)*($x2-$x1)+($y2-$y1)*($y2-$y1))]
    
    # create a material for the hinge
    set Ehinge [expr $n*6.0*$Ixx*$E/$L]
    set bhinge [expr $b/$n]
    set Mb [expr $Fy*$Zx]
    uniaxialMaterial Steel01 $eleTag $Mb $Ehinge $bhinge

    # add two zero length rotational hinges
    set one 1
    set two 2
    element zeroLength $eleTag$one $iNode $eleTag$hingeEnd1 -mat $eleTag -dir 6 -doRayleigh $doRayleigh
    element zeroLength $eleTag$two $eleTag$hingeEnd2 $jNode -mat $eleTag -dir 6 -doRayleigh $doRayleigh

    # add an elastic element in between
    ElasticBeamWSection2d $eleTag $eleTag$hingeEnd1 $eleTag$hingeEnd2 $sectType $E $transfTag
}


proc BeamWithPlasticHingesWSection2d {eleTag iNode jNode sectType E Fy H Lb Com_Type Comp_Action Lp transfTag args} {   
 #########################################################################################################
 #                                                                                                    
 # Creates a Finite-Length Plastic-Hinge (FLPH) element with two discrete hinges at both ends and a   
 #   linear elastic segment in between                                                                 
 #                                                                                                    
 # Flexural behavior of plastic hinge sections is defined using the Bilin02 model using procedure     
 #   SteelWSectionMChi02 (see below)                                                                  
 #                                                                                                    
 # Based on paper "DETERIORATION MODELING OF STEEL MOMENT RESISTING FRAMES USING FINITE-LENGTH             			 
 #                 PLASTIC HINGE FORCE-BASED BEAM-COLUMN ELEMENTS" 			              
 #  by:  F.L.A. Ribeiro, A.R. Barbosa, M.H. Scott, L.C. Neves    			             
 #  URL: http://web.engr.oregonstate.edu/~barbosa/products/ribeiro_barbosa_scott_neves.pdf 
 #
 # Procedure written by: F.L.A. Ribeiro and Andre Barbosa, FEB-12-2015                                
 # Contact: f.ribeiro@fct.unl.pt; andre.barbosa@oregonstate.edu           
 #
 # eleTag      - Element ID
 # iNode       - First node
 # jNode       - Second node
 # sectType    - Wsection used (see list below)
 # E           - Young's modulus (in MPa or ksi)               
 # Fy          - Yield stress (in MPa or ksi)                  
 # H           - Member Length without considering the panel zones (in mm or in)               
 # Lb          - Unbraced length from point of plastic hinge location to point of zero moment (in mm or in)   
 # Com_Type    - Type of component (Use: "other-than-RBS" for this procedure)                               
 # Comp_Action - Composite Action flag (Use: 1 (yes), 0 (No) )                                              
 # Lp          - Plastic Hinge length (assumed to be the same for both ends)   				
 # transfTag   - geometric transformation                                                                
 # args        - <YY> <-metric> <-release1> <-release2>                                                  
 #             -metric - activate this option if  in and ksi are used                                    
 #	           -release1 and/or release2 - activate this option if beam ends are hinged (no bending moment) 
 #                                                                                                          
 #########################################################################################################
    global SteelWSectionMChi02
    global in
	global WSection

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
	puts "YY orientation not handled - uses XX!"
    }
   
    if {[lsearch $args "-release1"] != -1} {
	set hingeEnd1 1
	node $eleTag$hingeEnd1 [nodeCoord $iNode 1] [nodeCoord $iNode 2]
	equalDOF $iNode $eleTag$hingeEnd1 1 2
	set iNode $eleTag$hingeEnd1
    }

    if {[lsearch $args "-release2"] != -1} {
	set hingeEnd2 2
	node $eleTag$hingeEnd2 [nodeCoord $jNode 1] [nodeCoord $iNode 2]
	equalDOF $iNode $eleTag$hingeEnd2 1 2
	set jNode $eleTag$hingeEnd2
    }

    set node1Crds [nodeCoord $iNode]; set x1 [lindex $node1Crds 0]; set y1 [lindex $node1Crds 1]
    set node2Crds [nodeCoord $jNode]; set x2 [lindex $node2Crds 0]; set y2 [lindex $node2Crds 1]

    set L [expr sqrt(($x2-$x1)*($x2-$x1)+($y2-$y1)*($y2-$y1))]
    set beta1 [expr -6.0*(3*$L*$L*$Lp - 24*$L*$Lp*$Lp +32*$Lp*$Lp*$Lp)/($L*($L-8.0*$Lp)*($L-8.0*$Lp))]
    set beta2 [expr 3.0*(3*$L*$L*$L - 48*$L*$L*$Lp + 224*$L*$Lp*$Lp - 256*$Lp*$Lp*$Lp)/($L*(3*$L-16.0*$Lp)*(3*$L-16.0*$Lp))]
    
    set secTag  $eleTag
    set countExtra 2
    set secTag2 $eleTag$countExtra;
    incr countExtra 3
    set secTag3 $eleTag$countExtra;
    
    set found 0
    foreach {section prop} [array get WSection $sectType] {
	#AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
	set propList [split $prop]
	set A [expr [lindex $propList 0]*$in*$in]
	set Ix [expr [lindex $propList 5]*$in*$in*$in*$in]
	set Iy [expr [lindex $propList 6]*$in*$in*$in*$in]
	set found 1
    }
    
    section Elastic $secTag2 $E $A [expr $beta1*$Ix]
    section Elastic $secTag3 $E $A [expr $beta2*$Ix]
    
    if {[lsearch $args "-metric"] != -1} {
	SteelWSectionMChi02 $secTag $E $Fy $H $L $Lb $sectType $Com_Type $Comp_Action $Lp  -metric
    } else {
	SteelWSectionMChi02 $secTag $E $Fy $H $L $Lb $sectType $Com_Type $Comp_Action $Lp
    }
    
    uniaxialMaterial Elastic [expr $secTag*100] [expr $E*$A];
    section Aggregator $secTag [expr $secTag*100] P $secTag Mz;
    
    set Locations "0 [expr (8.0/3.0*$Lp)/$L] [expr (4.0*$Lp+($L-8*$Lp)/2*(1-1/sqrt(3)))/$L] [expr (4.0*$Lp+($L-8*$Lp)/2*(1+1/sqrt(3)))/$L] [expr ($L-8.0/3.0*$Lp)/$L] 1.0";
    set weights "[expr $Lp/$L] [expr 3.0*$Lp/$L] [expr (($L-8.0*$Lp)/2)/$L] [expr (($L-8.0*$Lp)/2)/$L] [expr 3.0*$Lp/$L] [expr $Lp/$L]"; 
    set secTags "$secTag $secTag2 $secTag3 $secTag3 $secTag2 $secTag";
    set integration "LowOrder 6 $secTags $Locations $weights";
    element forceBeamColumn $eleTag $iNode $jNode $transfTag $integration 
}


proc ElasticBeamHSSection2d {eleTag iNode jNode sectType E transfTag args} {
    global HSSection
    global in
    set found 0
    
    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
        set Orient "YY"
    }
    
    if {[lsearch $args "-release1"] != -1} {
	set hingeEnd1 1
	node $eleTag$hingeEnd1 [nodeCoord $iNode 1] [nodeCoord $iNode 2]
	equalDOF $iNode $eleTag$hingeEnd1 1 2
	set iNode $eleTag$hingeEnd1
    }
    
    if {[lsearch $args "-release2"] != -1} {
	set hingeEnd2 2
	node $eleTag$hingeEnd2 [nodeCoord $jNode 1] [nodeCoord $jNode 2]
	equalDOF $jNode $eleTag$hingeEnd2 1 2
	set jNode $eleTag$hingeEnd2
    }
    
    foreach {section prop} [array get HSSection $sectType] {
	set propList [split $prop]

	set A [expr [lindex $propList 1]*$in*$in]
	set Ixx [expr [lindex $propList 5]*$in*$in*$in*$in]
	set Iyy [expr [lindex $propList 9]*$in*$in*$in*$in]
	
	if {$Orient == "YY" } {
	    element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iyy $transfTag
	} else {
	    element elasticBeamColumn $eleTag $iNode $jNode $A $E $Ixx $transfTag
	}
	set found 1
    }
    
    if {$found == 0} {
	puts "ElasticBeamWSection2d sectType: $sectType not found for ee: $eleTag"
    }
}

proc HSSbrace {eleTag iNode jNode sectType matTag numSeg Im transfTag args} {
    
    # This procedure develops a 2D brace element in a 2D/3D system (Z coordinates were set to 0).
    # Corotational Transformation is used by default
    #
    # Developed by Dimitrios G. Lignos, PhD
    # Contact: dimitrios.lignos@mcgill.ca
    #
    # Uses a displacement based element
    # Last Modified: 06/08/2014
    
    # args: 
    #  eleTag - element number (needed to provide node and ele tags for each beam segment and nodes
    #  iNode
    #  jNode
    #  sectType - HSS section design
    #  numSeg - num ele divisions
    #  Im - offset
    #  numInt - numIntegration points in beams
    #  transfTag - transformation tag
    
    global FiberHSSection2d
    global ElasticSteelHSSection2d
    
    set nip 3
    if {[lsearch $args "-nip"] != -1} {
	set loc [lsearch $args "-nip"]
        set nip [lindex $args [expr $loc+1]]
    }

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
        set Orient "YY"
    }

    set nFlange 4
    set nWeb 5
    if {[lsearch $args "-elasticSection"] != -1} {
	set loc [lsearch $args "-elasticSection"] 
        set E [lindex $args [expr $loc+1]]
	ElasticHSSSection2d $eleTag $sectType $E  $Orient
    } else {
	FiberHSSection2d $eleTag $sectType $matTag $nFlange $nWeb $Orient
    }

    set PI [expr 2*asin(1.0)];	# define constant pi

    # To get the coordinates of these 2 points
    set X1 [nodeCoord $iNode 1]
    set Y1 [nodeCoord $iNode 2]
    set X2 [nodeCoord $jNode 1]
    set Y2 [nodeCoord $jNode 2]

    # add nodes & boundary conditions for end releases
    set hingeEnd1 0
    node $eleTag$hingeEnd1 $X1 $Y1
    equalDOF $iNode $eleTag$hingeEnd1 1 2
    set iNode $eleTag$hingeEnd1

    set hingeEnd2 $numSeg
    node $eleTag$hingeEnd2 $X2 $Y2
    equalDOF $jNode $eleTag$hingeEnd2 1 2
    set jNode $eleTag$hingeEnd2

    # Set nodeID for the first intermediate node .. simply iNode with a 1 at the end
    set nodeID $iNode

    # Get the distance between the given points
    set L [expr sqrt(pow(($X2-$X1),2)+ pow(($Y2-$Y1),2))]

    # Get the sin and cos of the inclined angle
    set Cos [expr ($X2-$X1)/$L]
    set Sin [expr ($Y2-$Y1)/$L]
    
    for {set i 1} {$i <= [expr $numSeg-1]} {incr i 1} {
	# get the coordinates of each intermediate node in local system
	set nodeid [expr $eleTag$i]
	set xLocal [expr $L/$numSeg*$i]
	set yLocal [expr sin($PI*$i/$numSeg)*$Im*$L]
	
	set xRotX [expr $xLocal]
	set yRotX [expr $yLocal*0.707]
	
	set xRotZ [expr $xRotX*$Cos-$yRotX*$Sin]
	set yRotZ [expr $xRotX*$Sin+$yRotX*$Cos]
	
	# Use transformation matrix to convert the coordinate from local system to global system 
	set xGlobal [expr $X1+$xRotZ]
	set yGlobal [expr $Y1+$yRotZ]

	# add node
	node $nodeid $xGlobal $yGlobal
    }

    # Define segments
    set ElementID $eleTag

    # Define first element
    element dispBeamColumn $ElementID $iNode [expr $nodeID+1] $nip $eleTag $transfTag
    
    # Define internal elements #
    for {set i 1} {$i <[expr $numSeg-1]} {incr i 1} {
	set ElementID $eleTag$i
	# set some parameters #
	set iIntNode [expr $i +$nodeID]
	set jIntNode [expr $i +$nodeID+1]
	# add the Brace Element #
	element dispBeamColumn $ElementID $iIntNode $jIntNode $nip $eleTag $transfTag
    }

    # Define last element
    element dispBeamColumn $eleTag$numSeg $jIntNode $jNode $nip $eleTag $transfTag
}

proc BeamWithConcentratedHingesWSection2d {eleTag iNode jNode sectType E Fy H Lb Com_Type Comp_Action nFactor transfTag args} { 
 #########################################################################################################
 #                                                                                                    
 # Creates a Concentrated Plastic Hinge (CPH) element with two zero-length springs at both ends and a
 #                                    linear elastic element in between                                    
 #                                                                                                    
 # Zero-length springs moment-rotation behavior is defined using the Bilin02 modelusing procedure     
 #   SteelWSectionMR02 (see below)                                                                  
 #                                                                                                    
 # Based on paper under revision "Implementation and calibration of              			 
 #                 finite-length plastic hinge elements for use in seismic structural collapse analysis."  			              
 #  by:  F.L.A. Ribeiro, L.C. Neves, A.R. Barbosa
 #  URL: TBD
 #
 # Procedure written by: F.L.A. Ribeiro and Andre Barbosa, FEB-12-2015                                
 # Contact: f.ribeiro@fct.unl.pt; andre.barbosa@oregonstate.edu    
 #
 #  eleTag      - Element ID
 #  iNode       - First node
 #  jNode       - Second node
 #  sectType    - Wsection used (see list below)    
 #  E           - Young's modulus (in MPa or ksi)   
 #  Fy          - Yield stress (in MPa or ksi)      
 #  H           - Member Length without considering the panel zones (in mm or in)  
 #  Lb          - Unbraced length from point of plastic hinge location to point of zero moment (in mm or in) 
 #  Com_Type    - Type of component (Use: "other-than-RBS" for this procedure)
 #  Comp_Action - Composite Action flag (Use: 1 (yes), 0 (No) )               
 #  nFactor     - Elastic stiffness amplification factor (to make the springs rigid - plastic) - suggested value: 1000 
 #  transfTag   - geometric transformation     
 #  args        - <YY> <-metric> <-release1> <-release2>  
 #              -metric - activate this option if  in and ksi are used	  
 #              -release1 and/or release2 - activate this option if beam ends are hinged (no bending moment)   
 #                                                                                                          
 #########################################################################################################	
 
    global SteelWSectionMR02
    global in
	global WSection

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
	puts "YY orientation not handled - uses XX!"
    }
	
    set hingeEnd1 1
    node $eleTag$hingeEnd1 [nodeCoord $iNode 1] [nodeCoord $iNode 2]
    equalDOF $iNode $eleTag$hingeEnd1 1 2
    set hingeEnd2 2
    node $eleTag$hingeEnd2 [nodeCoord $jNode 1] [nodeCoord $jNode 2]
    equalDOF $jNode $eleTag$hingeEnd2 1 2
    
    set node1Crds [nodeCoord $iNode]; set x1 [lindex $node1Crds 0]; set y1 [lindex $node1Crds 1]
    set node2Crds [nodeCoord $jNode]; set x2 [lindex $node2Crds 0]; set y2 [lindex $node2Crds 1]
    
    set L [expr sqrt(($x2-$x1)*($x2-$x1)+($y2-$y1)*($y2-$y1))]
    
    if {[lsearch $args "-metric"] != -1} {
	SteelWSectionMR02 $eleTag $E $Fy $H $L $Lb $sectType $Com_Type $Comp_Action -metric -nFactor $nFactor
    } else {
	SteelWSectionMR02 $eleTag $E $Fy $H $L $Lb $sectType $Com_Type $Comp_Action -nFactor $nFactor
    }
    
    if {[lsearch $args "-release1"] == -1} {
	element zeroLength $eleTag$hingeEnd1 $iNode  $eleTag$hingeEnd1 -mat $eleTag -dir 6
    }
    
    if {[lsearch $args "-release2"] == -1} {
	element zeroLength $eleTag$hingeEnd2 $jNode  $eleTag$hingeEnd2 -mat $eleTag -dir 6
    }
    
    set found 0
    foreach {section prop} [array get WSection $sectType] {
	set propList [split $prop]
	#AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
	set A [expr [lindex $propList 0]*$in*$in]
	set Ix [expr [lindex $propList 5]*$in*$in*$in*$in]
	set found 1
    }

    element elasticBeamColumn $eleTag $eleTag$hingeEnd1 $eleTag$hingeEnd2 $A $E $Ix $transfTag
}


#
# 4. PROCEDURES TO CREATE SECTIONS
#


proc ElasticSteelWSection2d {sectTag sectType E args} {
    global WSection
    global in

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
        set Orient "YY"
    }

    set found 0
    foreach {section prop} [array get WSection $sectType] {
	set propList [split $prop]
	#AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
	set A [expr [lindex $propList 0]*$in*$in]
	set Ixx [expr [lindex $propList 5]*$in*$in*$in*$in]
	set Iyy [expr [lindex $propList 6]*$in*$in*$in*$in]

	if {$Orient == "YY" } {
	    section Elastic $sectTag $E $A $Iyy
	} else {
	    section Elastic $sectTag $E $A $Ixx
	}
	set found 1
    }

    if {$found == 0} {
	puts "FiberSteelWSection2d sectType: $sectType not found for sectTag: $sectTag"
    }
}

proc FiberSteelWSection2d {sectTag sectType matTag nFlange nWeb args} {
    global WSection
    global in

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
        set Orient "YY"
    }

    set found 0
    foreach {section prop} [array get WSection $sectType] {
	set propList [split $prop]

	set  d [expr [lindex $propList 1]*$in]
	set bf [expr [lindex $propList 2]*$in]
	set tw [expr [lindex $propList 3]*$in]
	set tf [expr [lindex $propList 4]*$in]

	Wsection $sectTag $matTag $d $bf $tf $tw $nFlange 1 1 $nWeb $Orient
	set found 1
    }
    if {$found == 0} {
	puts "FiberSteelWSection2d sectType: $sectType not found for sectTag: $sectTag"
    }
}


proc Wsection {secID matID d bf tf tw nfdw nftw nfbf nftf {Orient XX}} {
    # ###################################################################
    # Wsection  $secID $matID $d $bf $tf $tw $nfdw $nftw $nfbf $nftf
    # ###################################################################
    # create a standard W section given the nominal section properties
    # written: Remo M. de Souza
    # date: 06/99
    # modified: 08/99  (according to the new general modelbuilder)
    # input parameters
    # secID - section ID number
    # matID - material ID number 
    # d  = nominal depth
    # tw = web thickness
    # bf = flange width
    # tf = flange thickness
    # nfdw = number of fibers along web depth 
    # nftw = number of fibers along web thickness
    # nfbf = number of fibers along flange width
    # nftf = number of fibers along flange thickness

    set dw [expr $d - 2 * $tf]
    set y1 [expr -$d/2.0]
    set y2 [expr -$dw/2.0]
    set y3 [expr  $dw/2.0]
    set y4 [expr  $d/2.0]
    
    set z1 [expr -$bf/2.0]
    set z2 [expr -$tw/2.0]
    set z3 [expr  $tw/2.0]
    set z4 [expr  $bf/2.0]

    if {$Orient == "Weak" || $Orient == "YY" } {
	set dw [expr $d - 2 * $tf]
	set z1 [expr -$d/2.0]
	set z2 [expr -$dw/2.0]
	set z3 [expr  $dw/2.0]
	set z4 [expr  $d/2.0]

	set y1 [expr  $bf/2.0]
	set y2 [expr  $tw/2.0]
	set y3 [expr -$tw/2.0]
	set y4 [expr -$bf/2.0]
	
	section fiberSec  $secID  {
	    patch quadr  $matID  $nfbf $nftf   $y1 $z3   $y1 $z4   $y4 $z4   $y4 $z3
	    patch quadr  $matID  $nftw $nfdw   $y2 $z3   $y3 $z3   $y3 $z2   $y2 $z2
	    patch quadr  $matID  $nfbf $nftf   $y1 $z1   $y1 $z2   $y4 $z2   $y4 $z1
	}
	
    } else {
	set dw [expr $d - 2 * $tf]
	set y1 [expr -$d/2.0]
	set y2 [expr -$dw/2.0]
	set y3 [expr  $dw/2.0]
	set y4 [expr  $d/2.0]

	set z1 [expr -$bf/2.0]
	set z2 [expr -$tw/2.0]
	set z3 [expr  $tw/2.0]
	set z4 [expr  $bf/2.0]
	
	section fiberSec  $secID  {
	    #                     nfIJ  nfJK    yI  zI    yJ  zJ    yK  zK    yL  zL
	    patch quadr  $matID  $nfbf $nftf   $y1 $z4   $y1 $z1   $y2 $z1   $y2 $z4
	    patch quadr  $matID  $nftw $nfdw   $y2 $z3   $y2 $z2   $y3 $z2   $y3 $z3
	    patch quadr  $matID  $nfbf $nftf   $y3 $z4   $y3 $z1   $y4 $z1   $y4 $z4
	}
    }
}

proc ElasticHSSection2d {sectTag sectType E args} {
    global HSSection
    global in

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
        set Orient "YY"
    }

    set found 0
    foreach {section prop} [array get HSSection $sectType] {
	set propList [split $prop]
	# AISC_Manual_Label "W A h b tdes Ix Zx Sx rx Iy Zy Sy ry J C
	set A [expr [lindex $propList 1]*$in*$in]
	set Ixx [expr [lindex $propList 5]*$in*$in*$in*$in]
	set Iyy [expr [lindex $propList 9]*$in*$in*$in*$in]

	if {$Orient == "YY" } {
	    section Elastic $sectTag $E $A $Iyy
	} else {
	    section Elastic $sectTag $E $A $Ixx
	}
	set found 1
    }

    if {$found == 0} {
	puts "FiberSteelWSection2d sectType: $sectType not found for sectTag: $sectTag"
    }
}

proc FiberHSSection2d {sectTag sectType matTag nFlange nWeb args} {
    global HSSection
    global HSSectionD
    global in

    set Orient "XX"
    if {[lsearch $args "YY"] != -1} {
        set Orient "YY"
    }

    set found 0
    foreach {section prop} [array get HSSection $sectType] {
	set propList [split $prop]
	# AISC_Manual_Label "W A h b tdes Ix Zx Sx rx Iy Zy Sy ry J C
	set  d [expr [lindex $propList 2]*$in]
	set  b [expr [lindex $propList 3]*$in]
	set  t [expr [lindex $propList 4]*$in]

	if {$Orient == "XX"} {
	    HSSectionD $sectTag $matTag $d $b $t $nFlange $nWeb 
	} else {
	    HSSectionD $sectTag $matTag $b $d $t $nFlange $nWeb 
	}
	set found 1
    }
    if {$found == 0} {
	puts "FiberSteelWSection2d sectType: $sectType not found for sectTag: $sectTag"
    }
}


proc HSSectionD {secID matID d b t nfdw nftw} {
    # create a standard HSS section given the nominal section properties
    # written: Dimitrios G. Lignos
    # date: 01/21/2010
    # modified: 01/22/2010  (according to the new general modelbuilder)

    # FMK: NOTE TO SELF THIS ONLY WORKS FOR SQUARE HSS's SEND LIGNOS AN EMAIL!
    # ANOTHER NOTE TO SELF - A LOT OF WASTED FIBER FOR 2D CASE - REWRITE
    
    set D $d; 
    set tf $t;
    
    section  fiberSec    $secID  { 
	# PatchAISC :    matTag    NSIJ  NSJK                  Iy                   Iz                   Jy             Jz             Ky             Kz             Ly             Lz 
	patch  quadr     $matID   $nfdw $nftw  +[expr $D/2. - $tf]        +[expr $D/2.]  +[expr $D/2. - $tf]        -[expr $D/2.]        +[expr $D/2.]        -[expr $D/2.]        +[expr $D/2.]        +[expr $D/2.] 
	patch  quadr     $matID   $nftw $nfdw  -[expr $D/2. - $tf]        +[expr $D/2.]  -[expr $D/2. - $tf]  +[expr $D/2. - $tf]  +[expr $D/2. - $tf]  +[expr $D/2. - $tf]  +[expr $D/2. - $tf]        +[expr $D/2.] 
	patch  quadr     $matID   $nftw $nfdw  -[expr $D/2. - $tf]  -[expr $D/2. - $tf]  -[expr $D/2. - $tf]        -[expr $D/2.]  +[expr $D/2. - $tf]        -[expr $D/2.]  +[expr $D/2. - $tf]  -[expr $D/2. - $tf] 
	patch  quadr     $matID   $nfdw $nftw        -[expr $D/2.]        +[expr $D/2.]        -[expr $D/2.]        -[expr $D/2.]  -[expr $D/2. - $tf]        -[expr $D/2.]  -[expr $D/2. - $tf]        +[expr $D/2.]
    } 
}


proc SteelWSectionMR {matTag E Fy Ix Sx H L d tw bf tf Lb ry Com_Type Comp_Action args} {

    ############################################################################################
    # Procedure to Construct the Modified IMK Material with Moment-Rotation curve for steel                              # 
    #                                                                                                                                                            # 
    # The input parameters for bare steel components (beams, columns) are based on the following papers:         # 
    #                                                                                                                                                            # 
    # 1. Lignos, D.G., Krawinkler, H. (2011). “Deterioration Modeling of Steel Components in Support of               # 
    #    Collapse Prediction of Steel Moment Frames under Earthquake Loading",                                              # 
    #    ASCE, Journal of Structural Engineering, Vol. 137 (11), 1291-1302.                                                        # 
    #                                                                                                                                                            # 
    # 2. Lignos, D.G., Krawinkler, H. (2013). “Development and Utilization of Structural                                       # 
    #    Component Databases for Performance-Based Earthquake Engineering”,                                                # 
    #    ASCE, Journal of Structural Engineering, Vol. 139 (NEES 2), 1382-1394.                                                # 
    #                                                                                                                                                            # 
    # The input parameters for composite steel beams are based on the following paper:                                     # 
    #                                                                                                                                                            # 
    # 1. Elkady, A., Lignos, D.G. (2014). “Modeling of the Composite Action in Fully Restrained                          # 
    #    Beam-to-Column Connections: Implications in the Seismic Design and Collapse Capacity of Steel           # 
    #    Special Moment Frames", Earthquake Engineering and Structural Dynamics, doi: 10.1002/eqe.2430.       # 
    #                                                                                                                                                            # 
    # Input Variables for Procedure                                                                                                                 # 
    #                                                                                                                                                            # 
    # SpringID    - Spring zerolength ID                                                                                                           # 
    # Node_i      - First node                                                                                                                           # 
    # Node_j      - Second node                                                                                                                      # 
    # E           - Young's modulus                                                                                                                   # 
    # Fy          - Yield stress                                                                                                                          # 
    # Ix          - Moment of inertia of section                                                                                                     # 
    # Sx          - Plastic modulus of section                                                                                                     # 
    # H           - Member Length without considering the panel zones                                                                 # 
    # L           - Shear Span                                                                                                                           # 
    # d           - Section depth                                                                                                                        # 
    # tw          - Web thickness                                                                                                                      # 
    # bf          - Flange width                                                                                                                          # 
    # tf          - Flange thickness                                                                                                                     # 
    # Lb          - Unbraced length from point of plastic hinge location to point of zero moment                               # 
    # ry          - radius of gyration with respect to the weak axis of the cross section                                           # 
    # Com_Type    - Type of component (Use: 'RBS' or 'other-than-RBS')                                                           # 
    # Comp_Action - Composite Action flag (Use: 1 (yes), 0 (No) )                                                                    # 
    ###########################################################################################

    # parameters c1, c2 for unit conversion if Imperial Units are used else these variables should be set equal to 1.0
    set c1 1.0
    set c2 1.0

    if {[lsearch $args "-metric"] != -1} {
	set c1 25.4; 
	set c2 6.895;
    }

    set hLength 1.0;
    if {[lsearch $args "-hLength"] != -1} {
	set loc [lsearch $args "-hLength"]
        set hLength [lindex $args [expr $loc+1]]
    }
    
    # Element flexural stiffness assuming that the element is in double curvature 
    set K [expr 6.*$E* $Ix / $H]; 
    
    if {$Com_Type == 'other-than-RBS'} {
	# Pre-capping plastic rotation
	set theta_p   [expr 0.0865 * pow(($d/$tw),-0.365)  * pow(($bf/2./$tf),-0.140) *  pow(($L/$d),0.340) * pow(($c1 * $d/533.),-0.721) * pow(($c2 * $Fy/355.),-0.230)];
	
	# Post-capping plastic rotation
	set theta_pc  [expr 5.63 * pow(($d/$tw),-0.565)  * pow(($bf/2./$tf),-0.800) *  pow(($c1 * $d/533.),-0.280)  * pow(($c2 * $Fy/355.),-0.430)];
	
	# Reference Cumulative Energy
	set Lmda      [expr 495.0 * pow(($d/$tw),-1.340)  * pow(($bf/2./$tf),-0.595) *  pow(($c2 * $Fy/355.),-0.360)];
    }
    
    if {$Com_Type == 'RBS'} {
	# Pre-capping plastic rotation 
	set theta_p   [expr 0.19 * pow(($d/$tw),-0.314) * pow(($bf/2./$tf),-0.100) *  pow(($Lb/$ry),-0.185) * pow(($L/$d),0.113) * pow(($c1 * $d/533.),-0.760) * pow(($c2 * $Fy/355.),-0.070)];
	
	# Post-capping plastic rotation
	set theta_pc  [expr 9.52 * pow(($d/$tw),-0.513) * pow(($bf/2./$tf),-0.863) *  pow(($Lb/$ry),-0.108) * pow(($c2 * $Fy/355.),-0.360)];
	
	# Reference Cumulative Energy
	set Lmda      [expr 585. * pow(($d/$tw),-1.140) * pow(($bf/2./$tf),-0.632) *  pow(($Lb/$ry),-0.205) * pow(($c2 * $Fy/355),-0.391)]
    }
    
    # Ultimate Chord Rotation
    set theta_u 0.2;
    
    # Residual Strength Factor
    set Res 0.4;
    
    # Effective Yield Moment for Positive and Negative Loading Direction
    set My_P [expr  1.1 * $Sx * $Fy]; 
    set My_N [expr -1.1 * $Sx * $Fy];
    
    set c_S 1.0; set c_C 1.0; set c_A 1.0; set c_K 1.0;
    
    # Check for Composite Action for Beam Springs: If yes adjust the spring parameters based on the ones proposed by Elkady and Lignos (2014)
    if {$Com_Type == 1} { 
	set theta_p_P   [expr 1.80*$theta_p];
	set theta_p_N   [expr 0.95*$theta_p];
	set theta_pc_P  [expr 1.35*$theta_pc];
	set theta_pc_N  [expr 0.95*$theta_pc];
	set D_P 1.15; set D_N 1.00;
	
	# Capping-to-Yield Flexural Strength for Positive and Negative Loading Directions
	set McMyP 1.30; set McMyN 1.05;
	
	# Effective Yield Moment for Positive and Negative Loading Direction after Slab Adjustment
	set My_P      [expr  1.35 * $My]; 
	set My_N      [expr -1.25 * $My];
	
	# If No composite Action is considered (Columns and Bare Beam cases)
    } else {
	set theta_p_P   $theta_p;
	set theta_p_N   $theta_p;
	set theta_pc_P  $theta_pc;
	set theta_pc_N  $theta_pc;
	set D_P 1.0; set D_N 1.0;
	
	# Capping-to-Yield Flexural Strength for Positive and Negative Loading Directions
	set McMyP 1.00; set McMyN 1.00;
    }
    
    # Strain hardening ratios for Positive and Negative Loading Directions
    set as_mem_p [expr  ($McMyP-1.)*$My_P/($theta_p_P * 6.*$E * $Ix/$H)];
    set as_mem_n [expr -($McMyN-1.)*$My_N/($theta_p_N * 6.*$E * $Ix/$H)];

    
    # Define Uniaxial Material Modified Ibarra-Medina-Krawinkler (IMK) Model with Bilinear Hysteretic response
    uniaxialMaterial Bilin $matTag $K $as_mem_p $as_mem_n $My_P $My_N $Lmda $Lmda $Lmda $Lmda $c_S $c_C $c_A $c_K $theta_p_P $theta_p_N $theta_pc_P $theta_pc_N $Res $Res $theta_u $theta_u $D_P $D_N
}


proc SteelWSectionMR02 {matTag E Fy H L Lb sectType Com_Type Comp_Action args} {
 ##############################################################################################                                         
 # Procedure to Construct a Moment-Rotation Curve Using the Bilin02 Model (for steel)
 #	
 # Written by: D. Lignos, Ph.D.                                                     
 # Adapted by: F.L.A. Ribeiro and Andre Barbosa, FEB-12-2015                        
 #	
 # The implementation of this model follows papers:   								
 #  
 # 1. Ribeiro, F., Neves, L., and Barbosa, A.(2015). "Implementation and calibration of finite-length 
 #    plastic hinge elements for use in seismic structural collapse analysis." Submitted to Journal of Earthquake	              
 #    Engineering         
 #	
 # The input parameters for bare steel components (beams, columns) are based on the following papers:            
 #	
 # 1. Lignos, D.G., Krawinkler, H. (2011). “Deterioration Modeling of Steel Components in Support of 
 #    Collapse Prediction of Steel Moment Frames under Earthquake Loading",                          
 #    ASCE, Journal of Structural Engineering, Vol. 137 (11), 1291-1302.                             
 # 
 # 2. Lignos, D.G., Krawinkler, H. (2013). “Development and Utilization of Structural
 #    Component Databases for Performance-Based Earthquake Engineering”,             
 #    ASCE, Journal of Structural Engineering, Vol. 139 (NEES 2), 1382-1394.         
 # 
 # The input parameters for composite steel beams are based on the following paper: 
 # 
 # 1. Elkady, A., Lignos, D.G. (2014). “Modeling of the Composite Action in Fully Restrained
 #    Beam-to-Column Connections: Implications in the Seismic Design and Collapse Capacity of Steel
 #    Special Moment Frames", Earthquake Engineering and Structural Dynamics, doi: 10.1002/eqe.2430.
 #
 # Input Variables for Procedure 
 # 
 # matTag      - Deterioration model ID     
 # E           - Young's modulus (in MPa or ksi)  
 # Fy          - Yield stress (in MPa or ksi)    
 # H           - Member Length without considering the panel zones (in mm or in) 
 # L           - Shear Span (in mm or in)                            
 # Lb          - Unbraced length from point of plastic hinge location to point of zero moment (in mm or in)   
 # Com_Type    - Type of component (Use: 'RBS' or 'other-than-RBS')   
 # Comp_Action - Composite Action flag (Use: 1 (yes), 0 (No) )      
 # args        - <-nFactor $x> <-metric>							
 #				-nFactor - elastic stiffness amplification factor   
 #				-metric - activate this option if  in and ksi are used
 ###############################################################################################
	
    global in
    global WSection
    
    # parameters c1, c2 for unit conversion if Imperial Units are used else these variables should be set equal to 1.0
    set c1 1.0
    set c2 1.0
    
    if {[lsearch $args "-metric"] != -1} {
	set c1 25.4; 
	set c2 6.895;
    }
    
    set found 0.
    foreach {section prop} [array get WSection $sectType] {
	set propList [split $prop]                              
	#AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
	set A [expr [lindex $propList 0]*$in*$in]
	set d [expr [lindex $propList 1]*$in]
	set bf [expr [lindex $propList 2]*$in]
	set tw [expr [lindex $propList 3]*$in]
	set tf [expr [lindex $propList 4]*$in]
	set Ix [expr [lindex $propList 5]*$in*$in*$in*$in]
	set Iy [expr [lindex $propList 6]*$in*$in*$in*$in]
	set Sx [expr [lindex $propList 8]*$in*$in*$in]
	set rx [expr [lindex $propList 9]*$in]
	set Sy [expr [lindex $propList 11]*$in*$in*$in]
	set ry [expr [lindex $propList 12]*$in]             
	
	set found 1.	
    }
    
    if {[lsearch $args "-nFactor"] != -1} {
	set loc [lsearch $args "-nFactor"]
	set nFactor [lindex $args [expr $loc+1]]
	set Ixx [expr ($nFactor + 1.)/$nFactor * $Ix]
    } else {
	set nFactor 0.
	set Ixx $Ix
    }		
    
    # Element flexural stiffness assuming that the element is in double curvature 
    set K [expr 6.*$E* $Ixx / $H]; 
    
    if {$Com_Type == "other-than-RBS"} {
	# Pre-capping plastic rotation
	set theta_p   [expr 0.0865 * pow(($d/$tw),-0.365)  * pow(($bf/2./$tf),-0.140) *  pow(($L/$d),0.340) * pow(($c1 * $d/533.),-0.721) * pow(($c2 * $Fy/355.),-0.230)];
	
	# Post-capping plastic rotation
	set theta_pc  [expr 5.63 * pow(($d/$tw),-0.565)  * pow(($bf/2./$tf),-0.800) *  pow(($c1 * $d/533.),-0.280)  * pow(($c2 * $Fy/355.),-0.430)];
	
	# Reference Cumulative Energy
	set Lmda      [expr 495.0 * pow(($d/$tw),-1.340)  * pow(($bf/2./$tf),-0.595) *  pow(($c2 * $Fy/355.),-0.360)];
    }
    
    if {$Com_Type == "RBS"} {
	# Pre-capping plastic rotation 
	set theta_p   [expr 0.19 * pow(($d/$tw),-0.314) * pow(($bf/2./$tf),-0.100) *  pow(($Lb/$ry),-0.185) * pow(($L/$d),0.113) * pow(($c1 * $d/533.),-0.760) * pow(($c2 * $Fy/355.),-0.070)];
	
	# Post-capping plastic rotation
	set theta_pc  [expr 9.52 * pow(($d/$tw),-0.513) * pow(($bf/2./$tf),-0.863) *  pow(($Lb/$ry),-0.108) * pow(($c2 * $Fy/355.),-0.360)];
	
	# Reference Cumulative Energy
	set Lmda      [expr 585. * pow(($d/$tw),-1.140) * pow(($bf/2./$tf),-0.632) *  pow(($Lb/$ry),-0.205) * pow(($c2 * $Fy/355),-0.391)];
    }
    
    # Ultimate Chord Rotation
    set theta_u 0.4;
    
    # Residual Strength Factor
    set Res 0.4;
    
    # Effective Yield Moment for Positive and Negative Loading Direction
    set My_P [expr  1.1 * $Sx * $Fy]; 
    set My_N [expr  -1.1 * $Sx * $Fy];
    
    set c_S 1.0; set c_C 1.0; set c_A 1.0; set c_K 1.0;
    
    # Check for Composite Action for Beam Springs: If yes adjust the spring parameters based on the ones proposed by Elkady and Lignos (2014)
    if {$Com_Type == 1} { 
	set theta_p_P   [expr 1.80*$theta_p];
	set theta_p_N   [expr 0.95*$theta_p];
	set theta_pc_P  [expr 1.35*$theta_pc];
	set theta_pc_N  [expr 0.95*$theta_pc];
	set D_P 1.15; set D_N 1.00;
	
	# Capping-to-Yield Flexural Strength for Positive and Negative Loading Directions
	set McMyP 1.30; set McMyN 1.05;
	
	# Effective Yield Moment for Positive and Negative Loading Direction after Slab Adjustment
	set My_P      [expr  1.35 * $My_P]; 
	set My_N      [expr  1.25 * $My_N];
	
	# If No composite Action is considered (Columns and Bare Beam cases)
    } else {
	set theta_p_P   $theta_p;
	set theta_p_N   $theta_p;
	set theta_pc_P  $theta_pc;
	set theta_pc_N  $theta_pc;
	set D_P 1.0; set D_N 1.0;
	
	# Capping-to-Yield Flexural Strength for Positive and Negative Loading Directions
	set McMyP 1.05; set McMyN 1.05;
    }
    
    # Strain hardening ratios for Positive and Negative Loading Directions
    set as_mem_p [expr   ($McMyP-1.)*$My_P/($theta_p_P * $K)];
    set as_mem_n [expr  -($McMyN-1.)*$My_N/($theta_p_N * $K)];
    
    # Define Uniaxial Material Modified Ibarra-Medina-Krawinkler (IMK) Model with Bilinear Hysteretic response
    uniaxialMaterial Bilin02 $matTag $K $as_mem_p $as_mem_n $My_P $My_N $Lmda $Lmda $Lmda $Lmda $c_S $c_C $c_A $c_K $theta_p_P $theta_p_N $theta_pc_P $theta_pc_N $Res $Res $theta_u $theta_u $D_P $D_N $nFactor
}


proc SteelWSectionMChi02 {matTag E Fy H L Lb sectType Com_Type Comp_Action Lp args} {
    ################################################################################################
    # Procedure to Construct a Moment-Curvature Curve Using the Bilin02 Model (for steel)  
    #	
    # Written by: F.L.A. Ribeiro and Andre Barbosa, FEB-12-2015                                
    # 	
    # The implementation of this model follows papers:   
    # 
    # 1. Ribeiro, F., Neves, L., and Barbosa, A.(2015). “Implementation and calibration of finite-length  
    #    plastic hinge elements for use in seismic structural collapse analysis.” Submitted to Journal of Earthquake
    #    Engineering
    #	
    # The input parameters for bare steel components (beams, columns) are based on the following papers: 
    #	
    # 1. Lignos, D.G., Krawinkler, H. (2011). “Deterioration Modeling of Steel Components in Support of  
    #    Collapse Prediction of Steel Moment Frames under Earthquake Loading",                           
    #    ASCE, Journal of Structural Engineering, Vol. 137 (11), 1291-1302.                              
    # 
    # 2. Lignos, D.G., Krawinkler, H. (2013). “Development and Utilization of Structural    
    #    Component Databases for Performance-Based Earthquake Engineering”,     
    #    ASCE, Journal of Structural Engineering, Vol. 139 (NEES 2), 1382-1394. 
    # 
    # The input parameters for composite steel beams are based on the following paper: 
    # 
    # 1. Elkady, A., Lignos, D.G. (2014). “Modeling of the Composite Action in Fully Restrained 
    #    Beam-to-Column Connections: Implications in the Seismic Design and Collapse Capacity of Steel  
    #    Special Moment Frames", Earthquake Engineering and Structural Dynamics, doi: 10.1002/eqe.2430. 
    # 
    # Input Variables for Procedure        
    # 
    # matTag      - Deterioration model ID               
    # E           - Young's modulus (in MPa or ksi)      
    # Fy          - Yield stress (in MPa or ksi)         
    # H           - Member Length without considering the panel zones (in mm or in) 
    # L           - Shear Span (in mm or in)                                        
    # Lb          - Unbraced length from point of plastic hinge location to point of zero moment (in mm or in) 
    # Com_Type    - Type of component (Use: 'RBS' or 'other-than-RBS')    
    # Comp_Action - Composite Action flag (Use: 1 (yes), 0 (No) )         
    # Lp          - Plastic hinge length									
    # args        - <-metric>										
    #				metric - activate this option if  in and ksi are used 
    ###############################################################################################
    
    global in
    global WSection
    
    # parameters c1, c2 for unit conversion if Imperial Units are used else these variables should be set equal to 1.0
    set c1 1.0
    set c2 1.0
    
    if {[lsearch $args "-metric"] != -1} {
	set c1 25.4; 
	set c2 6.895;
    }
    
    set found 0.                                                                      
    foreach {section prop} [array get WSection $sectType] {       

	set propList [split $prop]                                                    
	
	#AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
	set A [expr [lindex $propList 0]*$in*$in]
	set d [expr [lindex $propList 1]*$in]
	set bf [expr [lindex $propList 2]*$in]
	set tw [expr [lindex $propList 3]*$in]
	set tf [expr [lindex $propList 4]*$in]
	set Ix [expr [lindex $propList 5]*$in*$in*$in*$in]
	set Iy [expr [lindex $propList 6]*$in*$in*$in*$in]
	set Sx [expr [lindex $propList 8]*$in*$in*$in]
	set rx [expr [lindex $propList 9]*$in]
	set Sy [expr [lindex $propList 11]*$in*$in*$in]
	set ry [expr [lindex $propList 12]*$in]             

	set found 1.                                                                      
    }
	
    # Element flexural stiffness assuming that the element is in double curvature 
    set K [expr 6.*$E* $Ix / $H * $Lp]; 
    
    if {$Com_Type == "other-than-RBS"} {
	# Pre-capping plastic rotation
	set theta_p   [expr 1./$Lp * 0.0865 * pow(($d/$tw),-0.365)  * pow(($bf/2./$tf),-0.140) *  pow(($L/$d),0.340) * pow(($c1 * $d/533.),-0.721) * pow(($c2 * $Fy/355.),-0.230)];
	
	# Post-capping plastic rotation
	set theta_pc  [expr 1./$Lp * 5.63 * pow(($d/$tw),-0.565)  * pow(($bf/2./$tf),-0.800) *  pow(($c1 * $d/533.),-0.280)  * pow(($c2 * $Fy/355.),-0.430)];
	
	# Reference Cumulative Energy
	set Lmda      [expr 1./$Lp * 495.0 * pow(($d/$tw),-1.340)  * pow(($bf/2./$tf),-0.595) *  pow(($c2 * $Fy/355.),-0.360)];
    }
    
    if {$Com_Type == "RBS"} {
	# Pre-capping plastic rotation 
	set theta_p   [expr 1./$Lp * 0.19 * pow(($d/$tw),-0.314) * pow(($bf/2./$tf),-0.100) *  pow(($Lb/$ry),-0.185) * pow(($L/$d),0.113) * pow(($c1 * $d/533.),-0.760) * pow(($c2 * $Fy/355.),-0.070)];
	
	# Post-capping plastic rotation
	set theta_pc  [expr 1./$Lp * 9.52 * pow(($d/$tw),-0.513) * pow(($bf/2./$tf),-0.863) *  pow(($Lb/$ry),-0.108) * pow(($c2 * $Fy/355.),-0.360)];
	
	# Reference Cumulative Energy
	set Lmda      [expr 1./$Lp * 585. * pow(($d/$tw),-1.140) * pow(($bf/2./$tf),-0.632) *  pow(($Lb/$ry),-0.205) * pow(($c2 * $Fy/355),-0.391)]
    }
    
    # Ultimate Chord Rotation
    set theta_u [expr 1./$Lp * 0.4];
    
    # Residual Strength Factor
    set Res 0.4;
    
    # Effective Yield Moment for Positive and Negative Loading Direction
    set My_P [expr  1.1 * $Sx * $Fy]; 
    set My_N [expr  -1.1 * $Sx * $Fy];
    
    set c_S 1.0; set c_C 1.0; set c_A 1.0; set c_K 1.0;
    
    # Check for Composite Action for Beam Springs: If yes adjust the spring parameters based on the ones proposed by Elkady and Lignos (2014)
    if {$Com_Type == 1} { 
	set theta_p_P   [expr 1.80*$theta_p];
	set theta_p_N   [expr 0.95*$theta_p];
	set theta_pc_P  [expr 1.35*$theta_pc];
	set theta_pc_N  [expr 0.95*$theta_pc];
	set D_P 1.15; set D_N 1.00;
	
	# Capping-to-Yield Flexural Strength for Positive and Negative Loading Directions
	set McMyP 1.30; set McMyN 1.05;
	
	# Effective Yield Moment for Positive and Negative Loading Direction after Slab Adjustment
	set My_P      [expr  1.35 * $My_P]; 
	set My_N      [expr  1.25 * $My_N];
	
	# If No composite Action is considered (Columns and Bare Beam cases)
    } else {
	set theta_p_P   $theta_p;
	set theta_p_N   $theta_p;
	set theta_pc_P  $theta_pc;
	set theta_pc_N  $theta_pc;
	set D_P 1.0; set D_N 1.0;
	
	# Capping-to-Yield Flexural Strength for Positive and Negative Loading Directions
	set McMyP 1.05; set McMyN 1.05;
    }
    
    # Strain hardening ratios for Positive and Negative Loading Directions
    set as_mem_p [expr  ($McMyP-1.)*$My_P/($theta_p_P * 6.*$E * $Ix/$H * $Lp)];
    set as_mem_n [expr -($McMyN-1.)*$My_N/($theta_p_N * 6.*$E * $Ix/$H * $Lp)];
    
    # Define Uniaxial Material Modified Ibarra-Medina-Krawinkler (IMK) Model with Bilinear Hysteretic response
    uniaxialMaterial Bilin02 $matTag $K $as_mem_p $as_mem_n $My_P $My_N $Lmda $Lmda $Lmda $Lmda $c_S $c_C $c_A $c_K $theta_p_P $theta_p_N $theta_pc_P $theta_pc_N $Res $Res $theta_u $theta_u $D_P $D_N
}


#
# 5. AISC W Section Table
#


#AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
array set WSection {
W44X335 "98.5 44.0 15.9 1.03 1.77 31100 1200 1620 1410 17.8 236 150 3.49 74.7"
W44X290 "85.4 43.6 15.8 0.865 1.58 27000 1040 1410 1240 17.8 205 132 3.49 50.9"
W44X262 "77.2 43.3 15.8 0.785 1.42 24100 923 1270 1110 17.7 182 117 3.47 37.3"
W44X230 "67.8 42.9 15.8 0.710 1.22 20800 796 1100 971 17.5 157 101 3.43 24.9"
W40X593 "174 43.0 16.7 1.79 3.23 50400 2520 2760 2340 17.0 481 302 3.80 445"
W40X503 "148 42.1 16.4 1.54 2.76 41600 2040 2320 1980 16.8 394 249 3.72 277"
W40X431 "127 41.3 16.2 1.34 2.36 34800 1690 1960 1690 16.6 328 208 3.65 177"
W40X397 "117 41.0 16.1 1.22 2.20 32000 1540 1800 1560 16.6 300 191 3.64 142"
W40X372 "110 40.6 16.1 1.16 2.05 29600 1420 1680 1460 16.5 277 177 3.60 116"
W40X362 "106 40.6 16.0 1.12 2.01 28900 1380 1640 1420 16.5 270 173 3.60 109"
W40X324 "95.3 40.2 15.9 1.00 1.81 25600 1220 1460 1280 16.4 239 153 3.58 79.4"
W40X297 "87.3 39.8 15.8 0.930 1.65 23200 1090 1330 1170 16.3 215 138 3.54 61.2"
W40X277 "81.5 39.7 15.8 0.830 1.58 21900 1040 1250 1100 16.4 204 132 3.58 51.5"
W40X249 "73.5 39.4 15.8 0.750 1.42 19600 926 1120 993 16.3 182 118 3.55 38.1"
W40X215 "63.5 39.0 15.8 0.650 1.22 16700 803 964 859 16.2 156 101 3.54 24.8"
W40X199 "58.8 38.7 15.8 0.650 1.07 14900 695 869 770 16.0 137 88.2 3.45 18.3"
W40X392 "116 41.6 12.4 1.42 2.52 29900 803 1710 1440 16.1 212 130 2.64 172"
W40X331 "97.7 40.8 12.2 1.22 2.13 24700 644 1430 1210 15.9 172 106 2.57 105"
W40X327 "95.9 40.8 12.1 1.18 2.13 24500 640 1410 1200 16.0 170 105 2.58 103"
W40X294 "86.2 40.4 12.0 1.06 1.93 21900 562 1270 1080 15.9 150 93.5 2.55 76.6"
W40X278 "82.3 40.2 12.0 1.03 1.81 20500 521 1190 1020 15.8 140 87.1 2.52 65.0"
W40X264 "77.4 40.0 11.9 0.960 1.73 19400 493 1130 971 15.8 132 82.6 2.52 56.1"
W40X235 "69.1 39.7 11.9 0.830 1.58 17400 444 1010 875 15.9 118 74.6 2.54 41.3"
W40X211 "62.1 39.4 11.8 0.750 1.42 15500 390 906 786 15.8 105 66.1 2.51 30.4"
W40X183 "53.3 39.0 11.8 0.650 1.20 13200 331 774 675 15.7 88.3 56.0 2.49 19.3"
W40X167 "49.3 38.6 11.8 0.650 1.03 11600 283 693 600 15.3 76.0 47.9 2.40 14.0"
W40X149 "43.8 38.2 11.8 0.630 0.830 9800 229 598 513 15.0 62.2 38.8 2.29 9.36"
W36X652 "192 41.1 17.6 1.97 3.54 50600 3230 2910 2460 16.2 581 367 4.10 593"
W36X529 "156 39.8 17.2 1.61 2.91 39600 2490 2330 1990 16.0 454 289 4.00 327"
W36X487 "143 39.3 17.1 1.50 2.68 36000 2250 2130 1830 15.8 412 263 3.96 258"
W36X441 "130 38.9 17.0 1.36 2.44 32100 1990 1910 1650 15.7 368 235 3.92 194"
W36X395 "116 38.4 16.8 1.22 2.20 28500 1750 1710 1490 15.7 325 208 3.88 142"
W36X361 "106 38.0 16.7 1.12 2.01 25700 1570 1550 1350 15.6 293 188 3.85 109"
W36X330 "96.9 37.7 16.6 1.02 1.85 23300 1420 1410 1240 15.5 265 171 3.83 84.3"
W36X302 "89.0 37.3 16.7 0.945 1.68 21100 1300 1280 1130 15.4 241 156 3.82 64.3"
W36X282 "82.9 37.1 16.6 0.885 1.57 19600 1200 1190 1050 15.4 223 144 3.80 52.7"
W36X262 "77.2 36.9 16.6 0.840 1.44 17900 1090 1100 972 15.3 204 132 3.76 41.6"
W36X247 "72.5 36.7 16.5 0.800 1.35 16700 1010 1030 913 15.2 190 123 3.74 34.7"
W36X231 "68.2 36.5 16.5 0.760 1.26 15600 940 963 854 15.1 176 114 3.71 28.7"
W36X256 "75.3 37.4 12.2 0.960 1.73 16800 528 1040 895 14.9 137 86.5 2.65 52.9"
W36X232 "68.0 37.1 12.1 0.870 1.57 15000 468 936 809 14.8 122 77.2 2.62 39.6"
W36X210 "61.9 36.7 12.2 0.830 1.36 13200 411 833 719 14.6 107 67.5 2.58 28.0"
W36X194 "57.0 36.5 12.1 0.765 1.26 12100 375 767 664 14.6 97.7 61.9 2.56 22.2"
W36X182 "53.6 36.3 12.1 0.725 1.18 11300 347 718 623 14.5 90.7 57.6 2.55 18.5"
W36X170 "50.0 36.2 12.0 0.680 1.10 10500 320 668 581 14.5 83.8 53.2 2.53 15.1"
W36X160 "47.0 36.0 12.0 0.650 1.02 9760 295 624 542 14.4 77.3 49.1 2.50 12.4"
W36X150 "44.3 35.9 12.0 0.625 0.940 9040 270 581 504 14.3 70.9 45.1 2.47 10.1"
W36X135 "39.9 35.6 12.0 0.600 0.790 7800 225 509 439 14.0 59.7 37.7 2.38 7.00"
W33X387 "114 36.0 16.2 1.26 2.28 24300 1620 1560 1350 14.6 312 200 3.77 148"
W33X354 "104 35.6 16.1 1.16 2.09 22000 1460 1420 1240 14.5 282 181 3.74 115"
W33X318 "93.7 35.2 16.0 1.04 1.89 19500 1290 1270 1110 14.5 250 161 3.71 84.4"
W33X291 "85.6 34.8 15.9 0.960 1.73 17700 1160 1160 1020 14.4 226 146 3.68 65.1"
W33X263 "77.4 34.5 15.8 0.870 1.57 15900 1040 1040 919 14.3 202 131 3.66 48.7"
W33X241 "71.1 34.2 15.9 0.830 1.40 14200 933 940 831 14.1 182 118 3.62 36.2"
W33X221 "65.3 33.9 15.8 0.775 1.28 12900 840 857 759 14.1 164 106 3.59 27.8"
W33X201 "59.1 33.7 15.7 0.715 1.15 11600 749 773 686 14.0 147 95.2 3.56 20.8"
W33X169 "49.5 33.8 11.5 0.670 1.22 9290 310 629 549 13.7 84.4 53.9 2.50 17.7"
W33X152 "44.9 33.5 11.6 0.635 1.06 8160 273 559 487 13.5 73.9 47.2 2.47 12.4"
W33X141 "41.5 33.3 11.5 0.605 0.960 7450 246 514 448 13.4 66.9 42.7 2.43 9.70"
W33X130 "38.3 33.1 11.5 0.580 0.855 6710 218 467 406 13.2 59.5 37.9 2.39 7.37"
W33X118 "34.7 32.9 11.5 0.550 0.740 5900 187 415 359 13.0 51.3 32.6 2.32 5.30"
W30X391 "115 33.2 15.6 1.36 2.44 20700 1550 1450 1250 13.4 310 198 3.67 173"
W30X357 "105 32.8 15.5 1.24 2.24 18700 1390 1320 1140 13.3 279 179 3.64 134"
W30X326 "95.9 32.4 15.4 1.14 2.05 16800 1240 1190 1040 13.2 252 162 3.60 103"
W30X292 "86.0 32.0 15.3 1.02 1.85 14900 1100 1060 930 13.2 223 144 3.58 75.2"
W30X261 "77.0 31.6 15.2 0.930 1.65 13100 959 943 829 13.1 196 127 3.53 54.1"
W30X235 "69.3 31.3 15.1 0.830 1.50 11700 855 847 748 13.0 175 114 3.51 40.3"
W30X211 "62.3 30.9 15.1 0.775 1.32 10300 757 751 665 12.9 155 100 3.49 28.4"
W30X191 "56.1 30.7 15.0 0.710 1.19 9200 673 675 600 12.8 138 89.5 3.46 21.0"
W30X173 "50.9 30.4 15.0 0.655 1.07 8230 598 607 541 12.7 123 79.8 3.42 15.6"
W30X148 "43.6 30.7 10.5 0.650 1.18 6680 227 500 436 12.4 68.0 43.3 2.28 14.5"
W30X132 "38.8 30.3 10.5 0.615 1.00 5770 196 437 380 12.2 58.4 37.2 2.25 9.72"
W30X124 "36.5 30.2 10.5 0.585 0.930 5360 181 408 355 12.1 54.0 34.4 2.23 7.99"
W30X116 "34.2 30.0 10.5 0.565 0.850 4930 164 378 329 12.0 49.2 31.3 2.19 6.43"
W30X108 "31.7 29.8 10.5 0.545 0.760 4470 146 346 299 11.9 43.9 27.9 2.15 4.99"
W30X99  "29.0 29.7 10.5 0.520 0.670 3990 128 312 269 11.7 38.6 24.5 2.10 3.77"
W30X90  "26.3 29.5 10.4 0.470 0.610 3610 115 283 245 11.7 34.7 22.1 2.09 2.84"
W27X539 "159 32.5 15.3 1.97 3.54 25600 2110 1890 1570 12.7 437 277 3.65 496"
W27X368 "109 30.4 14.7 1.38 2.48 16200 1310 1240 1060 12.2 279 179 3.48 170"
W27X336 "99.2 30.0 14.6 1.26 2.28 14600 1180 1130 972 12.1 252 162 3.45 131"
W27X307 "90.2 29.6 14.4 1.16 2.09 13100 1050 1030 887 12.0 227 146 3.41 101"
W27X281 "83.1 29.3 14.4 1.06 1.93 11900 953 936 814 12.0 206 133 3.39 79.5"
W27X258 "76.1 29.0 14.3 0.980 1.77 10800 859 852 745 11.9 187 120 3.36 61.6"
W27X235 "69.4 28.7 14.2 0.910 1.61 9700 769 772 677 11.8 168 108 3.33 47.0"
W27X217 "63.9 28.4 14.1 0.830 1.50 8910 704 711 627 11.8 154 100 3.32 37.6"
W27X194 "57.1 28.1 14.0 0.750 1.34 7860 619 631 559 11.7 136 88.1 3.29 27.1"
W27X178 "52.5 27.8 14.1 0.725 1.19 7020 555 570 505 11.6 122 78.8 3.25 20.1"
W27X161 "47.6 27.6 14.0 0.660 1.08 6310 497 515 458 11.5 109 70.9 3.23 15.1"
W27X146 "43.2 27.4 14.0 0.605 0.975 5660 443 464 414 11.5 97.7 63.5 3.20 11.3"
W27X129 "37.8 27.6 10.0 0.610 1.10 4760 184 395 345 11.2 57.6 36.8 2.21 11.1"
W27X114 "33.6 27.3 10.1 0.570 0.930 4080 159 343 299 11.0 49.3 31.5 2.18 7.33"
W27X102 "30.0 27.1 10.0 0.515 0.830 3620 139 305 267 11.0 43.4 27.8 2.15 5.28"
W27X94  "27.6 26.9 10.0 0.490 0.745 3270 124 278 243 10.9 38.8 24.8 2.12 4.03"
W27X84  "24.7 26.7 10.0 0.460 0.640 2850 106 244 213 10.7 33.2 21.2 2.07 2.81"
W24X370 "109 28.0 13.7 1.52 2.72 13400 1160 1130 957 11.1 267 170 3.27 201"
W24X335 "98.3 27.5 13.5 1.38 2.48 11900 1030 1020 864 11.0 238 152 3.23 152"
W24X306 "89.7 27.1 13.4 1.26 2.28 10700 919 922 789 10.9 214 137 3.20 117"
W24X279 "81.9 26.7 13.3 1.16 2.09 9600 823 835 718 10.8 193 124 3.17 90.5"
W24X250 "73.5 26.3 13.2 1.04 1.89 8490 724 744 644 10.7 171 110 3.14 66.6"
W24X229 "67.2 26.0 13.1 0.960 1.73 7650 651 675 588 10.7 154 99.4 3.11 51.3"
W24X207 "60.7 25.7 13.0 0.870 1.57 6820 578 606 531 10.6 137 88.8 3.08 38.3"
W24X192 "56.5 25.5 13.0 0.810 1.46 6260 530 559 491 10.5 126 81.8 3.07 30.8"
W24X176 "51.7 25.2 12.9 0.750 1.34 5680 479 511 450 10.5 115 74.3 3.04 23.9"
W24X162 "47.8 25.0 13.0 0.705 1.22 5170 443 468 414 10.4 105 68.4 3.05 18.5"
W24X146 "43.0 24.7 12.9 0.650 1.09 4580 391 418 371 10.3 93.2 60.5 3.01 13.4"
W24X131 "38.6 24.5 12.9 0.605 0.960 4020 340 370 329 10.2 81.5 53.0 2.97 9.50"
W24X117 "34.4 24.3 12.8 0.550 0.850 3540 297 327 291 10.1 71.4 46.5 2.94 6.72"
W24X104 "30.7 24.1 12.8 0.500 0.750 3100 259 289 258 10.1 62.4 40.7 2.91 4.72"
W24X103 "30.3 24.5 9.00 0.550 0.980 3000 119 280 245 10.0 41.5 26.5 1.99 7.07"
W24X94 "27.7 24.3 9.07 0.515 0.875 2700 109 254 222 9.87 37.5 24.0 1.98 5.26"
W24X84 "24.7 24.1 9.02 0.470 0.770 2370 94.4 224 196 9.79 32.6 20.9 1.95 3.70"
W24X76 "22.4 23.9 8.99 0.440 0.680 2100 82.5 200 176 9.69 28.6 18.4 1.92 2.68"
W24X68 "20.1 23.7 8.97 0.415 0.585 1830 70.4 177 154 9.55 24.5 15.7 1.87 1.87"
W24X62 "18.2 23.7 7.04 0.430 0.590 1550 34.5 153 131 9.23 15.7 9.80 1.38 1.71"
W24X55 "16.2 23.6 7.01 0.395 0.505 1350 29.1 134 114 9.11 13.3 8.30 1.34 1.18"
W21X201 "59.3 23.0 12.6 0.910 1.63 5310 542 530 461 9.47 133 86.1 3.02 40.9"
W21X182 "53.6 22.7 12.5 0.830 1.48 4730 483 476 417 9.40 119 77.2 3.00 30.7"
W21X166 "48.8 22.5 12.4 0.750 1.36 4280 435 432 380 9.36 108 70.0 2.99 23.6"
W21X147 "43.2 22.1 12.5 0.720 1.15 3630 376 373 329 9.17 92.6 60.1 2.95 15.4"
W21X132 "38.8 21.8 12.4 0.650 1.04 3220 333 333 295 9.12 82.3 53.5 2.93 11.3"
W21X122 "35.9 21.7 12.4 0.600 0.960 2960 305 307 273 9.09 75.6 49.2 2.92 8.98"
W21X111 "32.6 21.5 12.3 0.550 0.875 2670 274 279 249 9.05 68.2 44.5 2.90 6.83"
W21X101 "29.8 21.4 12.3 0.500 0.800 2420 248 253 227 9.02 61.7 40.3 2.89 5.21"
W21X93 "27.3 21.6 8.42 0.580 0.930 2070 92.9 221 192 8.70 34.7 22.1 1.84 6.03"
W21X83 "24.4 21.4 8.36 0.515 0.835 1830 81.4 196 171 8.67 30.5 19.5 1.83 4.34"
W21X73 "21.5 21.2 8.30 0.455 0.740 1600 70.6 172 151 8.64 26.6 17.0 1.81 3.02"
W21X68 "20.0 21.1 8.27 0.430 0.685 1480 64.7 160 140 8.60 24.4 15.7 1.80 2.45"
W21X62 "18.3 21.0 8.24 0.400 0.615 1330 57.5 144 127 8.54 21.7 14.0 1.77 1.83"
W21X55 "16.2 20.8 8.22 0.375 0.522 1140 48.4 126 110 8.40 18.4 11.8 1.73 1.24"
W21X48 "14.1 20.6 8.14 0.350 0.430 959 38.7 107 93.0 8.24 14.9 9.52 1.66 0.803"
W21X57 "16.7 21.1 6.56 0.405 0.650 1170 30.6 129 111 8.36 14.8 9.35 1.35 1.77"
W21X50 "14.7 20.8 6.53 0.380 0.535 984 24.9 110 94.5 8.18 12.2 7.64 1.30 1.14"
W21X44 "13.0 20.7 6.50 0.350 0.450 843 20.7 95.4 81.6 8.06 10.2 6.37 1.26 0.770"
W18X311 "91.6 22.3 12.0 1.52 2.74 6970 795 754 624 8.72 207 132 2.95 176"
W18X283 "83.3 21.9 11.9 1.40 2.50 6170 704 676 565 8.61 185 118 2.91 134"
W18X258 "76.0 21.5 11.8 1.28 2.30 5510 628 611 514 8.53 166 107 2.88 103"
W18X234 "68.6 21.1 11.7 1.16 2.11 4900 558 549 466 8.44 149 95.8 2.85 78.7"
W18X211 "62.3 20.7 11.6 1.06 1.91 4330 493 490 419 8.35 132 85.3 2.82 58.6"
W18X192 "56.2 20.4 11.5 0.960 1.75 3870 440 442 380 8.28 119 76.8 2.79 44.7"
W18X175 "51.4 20.0 11.4 0.890 1.59 3450 391 398 344 8.20 106 68.8 2.76 33.8"
W18X158 "46.3 19.7 11.3 0.810 1.44 3060 347 356 310 8.12 94.8 61.4 2.74 25.2"
W18X143 "42.0 19.5 11.2 0.730 1.32 2750 311 322 282 8.09 85.4 55.5 2.72 19.2"
W18X130 "38.3 19.3 11.2 0.670 1.20 2460 278 290 256 8.03 76.7 49.9 2.70 14.5"
W18X119 "35.1 19.0 11.3 0.655 1.06 2190 253 262 231 7.90 69.1 44.9 2.69 10.6"
W18X106 "31.1 18.7 11.2 0.590 0.940 1910 220 230 204 7.84 60.5 39.4 2.66 7.48"
W18X97 "28.5 18.6 11.1 0.535 0.870 1750 201 211 188 7.82 55.3 36.1 2.65 5.86"
W18X86 "25.3 18.4 11.1 0.480 0.770 1530 175 186 166 7.77 48.4 31.6 2.63 4.10"
W18X76 "22.3 18.2 11.0 0.425 0.680 1330 152 163 146 7.73 42.2 27.6 2.61 2.83"
W18X71 "20.9 18.5 7.64 0.495 0.810 1170 60.3 146 127 7.50 24.7 15.8 1.70 3.49"
W18X65 "19.1 18.4 7.59 0.450 0.750 1070 54.8 133 117 7.49 22.5 14.4 1.69 2.73"
W18X60 "17.6 18.2 7.56 0.415 0.695 984 50.1 123 108 7.47 20.6 13.3 1.68 2.17"
W18X55 "16.2 18.1 7.53 0.390 0.630 890 44.9 112 98.3 7.41 18.5 11.9 1.67 1.66"
W18X50 "14.7 18.0 7.50 0.355 0.570 800 40.1 101 88.9 7.38 16.6 10.7 1.65 1.24"
W18X46 "13.5 18.1 6.06 0.360 0.605 712 22.5 90.7 78.8 7.25 11.7 7.43 1.29 1.22"
W18X40 "11.8 17.9 6.02 0.315 0.525 612 19.1 78.4 68.4 7.21 10.0 6.35 1.27 0.810"
W18X35 "10.3 17.7 6.00 0.300 0.425 510 15.3 66.5 57.6 7.04 8.06 5.12 1.22 0.506"
W16X100 "29.4 17.0 10.4 0.585 0.985 1490 186 198 175 7.10 54.9 35.7 2.51 7.73"
W16X89 "26.2 16.8 10.4 0.525 0.875 1300 163 175 155 7.05 48.1 31.4 2.49 5.45"
W16X77 "22.6 16.5 10.3 0.455 0.760 1110 138 150 134 7.00 41.1 26.9 2.47 3.57"
W16X67 "19.6 16.3 10.2 0.395 0.665 954 119 130 117 6.96 35.5 23.2 2.46 2.39"
W16X57 "16.8 16.4 7.12 0.430 0.715 758 43.1 105 92.2 6.72 18.9 12.1 1.60 2.22"
W16X50 "14.7 16.3 7.07 0.380 0.630 659 37.2 92.0 81.0 6.68 16.3 10.5 1.59 1.52"
W16X45 "13.3 16.1 7.04 0.345 0.565 586 32.8 82.3 72.7 6.65 14.5 9.34 1.57 1.11"
W16X40 "11.8 16.0 7.00 0.305 0.505 518 28.9 73.0 64.7 6.63 12.7 8.25 1.57 0.794"
W16X36 "10.6 15.9 6.99 0.295 0.430 448 24.5 64.0 56.5 6.51 10.8 7.00 1.52 0.545"
W16X31 "9.13 15.9 5.53 0.275 0.440 375 12.4 54.0 47.2 6.41 7.03 4.49 1.17 0.461"
W16X26 "7.68 15.7 5.50 0.250 0.345 301 9.59 44.2 38.4 6.26 5.48 3.49 1.12 0.262"
W14X730 "215 22.4 17.9 3.07 4.91 14300 4720 1660 1280 8.17 816 527 4.69 1450"
W14X665 "196 21.6 17.7 2.83 4.52 12400 4170 1480 1150 7.98 730 472 4.62 1120"
W14X605 "178 20.9 17.4 2.60 4.16 10800 3680 1320 1040 7.80 652 423 4.55 869"
W14X550 "162 20.2 17.2 2.38 3.82 9430 3250 1180 931 7.63 583 378 4.49 669"
W14X500 "147 19.6 17.0 2.19 3.50 8210 2880 1050 838 7.48 522 339 4.43 514"
W14X455 "134 19.0 16.8 2.02 3.21 7190 2560 936 756 7.33 468 304 4.38 395"
W14X426 "125 18.7 16.7 1.88 3.04 6600 2360 869 706 7.26 434 283 4.34 331"
W14X398 "117 18.3 16.6 1.77 2.85 6000 2170 801 656 7.16 402 262 4.31 273"
W14X370 "109 17.9 16.5 1.66 2.66 5440 1990 736 607 7.07 370 241 4.27 222"
W14X342 "101 17.5 16.4 1.54 2.47 4900 1810 672 558 6.98 338 221 4.24 178"
W14X311 "91.4 17.1 16.2 1.41 2.26 4330 1610 603 506 6.88 304 199 4.20 136"
W14X283 "83.3 16.7 16.1 1.29 2.07 3840 1440 542 459 6.79 274 179 4.17 104"
W14X257 "75.6 16.4 16.0 1.18 1.89 3400 1290 487 415 6.71 246 161 4.13 79.1"
W14X233 "68.5 16.0 15.9 1.07 1.72 3010 1150 436 375 6.63 221 145 4.10 59.5"
W14X211 "62.0 15.7 15.8 0.980 1.56 2660 1030 390 338 6.55 198 130 4.07 44.6"
W14X193 "56.8 15.5 15.7 0.890 1.44 2400 931 355 310 6.50 180 119 4.05 34.8"
W14X176 "51.8 15.2 15.7 0.830 1.31 2140 838 320 281 6.43 163 107 4.02 26.5"
W14X159 "46.7 15.0 15.6 0.745 1.19 1900 748 287 254 6.38 146 96.2 4.00 19.7"
W14X145 "42.7 14.8 15.5 0.680 1.09 1710 677 260 232 6.33 133 87.3 3.98 15.2"
W14X132 "38.8 14.7 14.7 0.645 1.03 1530 548 234 209 6.28 113 74.5 3.76 12.3"
W14X120 "35.3 14.5 14.7 0.590 0.940 1380 495 212 190 6.24 102 67.5 3.74 9.37"
W14X109 "32.0 14.3 14.6 0.525 0.860 1240 447 192 173 6.22 92.7 61.2 3.73 7.12"
W14X99 "29.1 14.2 14.6 0.485 0.780 1110 402 173 157 6.17 83.6 55.2 3.71 5.37"
W14X90 "26.5 14.0 14.5 0.440 0.710 999 362 157 143 6.14 75.6 49.9 3.70 4.06"
W14X82 "24.0 14.3 10.1 0.510 0.855 881 148 139 123 6.05 44.8 29.3 2.48 5.07"
W14X74 "21.8 14.2 10.1 0.450 0.785 795 134 126 112 6.04 40.5 26.6 2.48 3.87"
W14X68 "20.0 14.0 10.0 0.415 0.720 722 121 115 103 6.01 36.9 24.2 2.46 3.01"
W14X61 "17.9 13.9 10.0 0.375 0.645 640 107 102 92.1 5.98 32.8 21.5 2.45 2.19"
W14X53 "15.6 13.9 8.06 0.370 0.660 541 57.7 87.1 77.8 5.89 22.0 14.3 1.92 1.94"
W14X48 "14.1 13.8 8.03 0.340 0.595 484 51.4 78.4 70.2 5.85 19.6 12.8 1.91 1.45"
W14X43 "12.6 13.7 8.00 0.305 0.530 428 45.2 69.6 62.6 5.82 17.3 11.3 1.89 1.05"
W14X38 "11.2 14.1 6.77 0.310 0.515 385 26.7 61.5 54.6 5.87 12.1 7.88 1.55 0.798"
W14X34 "10.0 14.0 6.75 0.285 0.455 340 23.3 54.6 48.6 5.83 10.6 6.91 1.53 0.569"
W14X30 "8.85 13.8 6.73 0.270 0.385 291 19.6 47.3 42.0 5.73 8.99 5.82 1.49 0.380"
W14X26 "7.69 13.9 5.03 0.255 0.420 245 8.91 40.2 35.3 5.65 5.54 3.55 1.08 0.358"
W14X22 "6.49 13.7 5.00 0.230 0.335 199 7.00 33.2 29.0 5.54 4.39 2.80 1.04 0.208"
W12X336 "98.9 16.8 13.4 1.78 2.96 4060 1190 603 483 6.41 274 177 3.47 243"
W12X305 "89.5 16.3 13.2 1.63 2.71 3550 1050 537 435 6.29 244 159 3.42 185"
W12X279 "81.9 15.9 13.1 1.53 2.47 3110 937 481 393 6.16 220 143 3.38 143"
W12X252 "74.1 15.4 13.0 1.40 2.25 2720 828 428 353 6.06 196 127 3.34 108"
W12X230 "67.7 15.1 12.9 1.29 2.07 2420 742 386 321 5.97 177 115 3.31 83.8"
W12X210 "61.8 14.7 12.8 1.18 1.90 2140 664 348 292 5.89 159 104 3.28 64.7"
W12X190 "56.0 14.4 12.7 1.06 1.74 1890 589 311 263 5.82 143 93.0 3.25 48.8"
W12X170 "50.0 14.0 12.6 0.960 1.56 1650 517 275 235 5.74 126 82.3 3.22 35.6"
W12X152 "44.7 13.7 12.5 0.870 1.40 1430 454 243 209 5.66 111 72.8 3.19 25.8"
W12X136 "39.9 13.4 12.4 0.790 1.25 1240 398 214 186 5.58 98.0 64.2 3.16 18.5"
W12X120 "35.2 13.1 12.3 0.710 1.11 1070 345 186 163 5.51 85.4 56.0 3.13 12.9"
W12X106 "31.2 12.9 12.2 0.610 0.990 933 301 164 145 5.47 75.1 49.3 3.11 9.13"
W12X96 "28.2 12.7 12.2 0.550 0.900 833 270 147 131 5.44 67.5 44.4 3.09 6.85"
W12X87 "25.6 12.5 12.1 0.515 0.810 740 241 132 118 5.38 60.4 39.7 3.07 5.10"
W12X79 "23.2 12.4 12.1 0.470 0.735 662 216 119 107 5.34 54.3 35.8 3.05 3.84"
W12X72 "21.1 12.3 12.0 0.430 0.670 597 195 108 97.4 5.31 49.2 32.4 3.04 2.93"
W12X65 "19.1 12.1 12.0 0.390 0.605 533 174 96.8 87.9 5.28 44.1 29.1 3.02 2.18"
W12X58 "17.0 12.2 10.0 0.360 0.640 475 107 86.4 78.0 5.28 32.5 21.4 2.51 2.10"
W12X53 "15.6 12.1 10.0 0.345 0.575 425 95.8 77.9 70.6 5.23 29.1 19.2 2.48 1.58"
W12X50 "14.6 12.2 8.08 0.370 0.640 391 56.3 71.9 64.2 5.18 21.3 13.9 1.96 1.71"
W12X45 "13.1 12.1 8.05 0.335 0.575 348 50.0 64.2 57.7 5.15 19.0 12.4 1.95 1.26"
W12X40 "11.7 11.9 8.01 0.295 0.515 307 44.1 57.0 51.5 5.13 16.8 11.0 1.94 0.906"
W12X35 "10.3 12.5 6.56 0.300 0.520 285 24.5 51.2 45.6 5.25 11.5 7.47 1.54 0.741"
W12X30 "8.79 12.3 6.52 0.260 0.440 238 20.3 43.1 38.6 5.21 9.56 6.24 1.52 0.457"
W12X26 "7.65 12.2 6.49 0.230 0.380 204 17.3 37.2 33.4 5.17 8.17 5.34 1.51 0.300"
W12X22 "6.48 12.3 4.03 0.260 0.425 156 4.66 29.3 25.4 4.91 3.66 2.31 0.848 0.293"
W12X19 "5.57 12.2 4.01 0.235 0.350 130 3.76 24.7 21.3 4.82 2.98 1.88 0.822 0.180"
W12X16 "4.71 12.0 3.99 0.220 0.265 103 2.82 20.1 17.1 4.67 2.26 1.41 0.773 0.103"
W12X14 "4.16 11.9 3.97 0.200 0.225 88.6 2.36 17.4 14.9 4.62 1.90 1.19 0.753 0.0704"
W10X112 "32.9 11.4 10.4 0.755 1.25 716 236 147 126 4.66 69.2 45.3 2.68 15.1"
W10X100 "29.3 11.1 10.3 0.680 1.12 623 207 130 112 4.60 61.0 40.0 2.65 10.9"
W10X88 "26.0 10.8 10.3 0.605 0.990 534 179 113 98.5 4.54 53.1 34.8 2.63 7.53"
W10X77 "22.7 10.6 10.2 0.530 0.870 455 154 97.6 85.9 4.49 45.9 30.1 2.60 5.11"
W10X68 "19.9 10.4 10.1 0.470 0.770 394 134 85.3 75.7 4.44 40.1 26.4 2.59 3.56"
W10X60 "17.7 10.2 10.1 0.420 0.680 341 116 74.6 66.7 4.39 35.0 23.0 2.57 2.48"
W10X54 "15.8 10.1 10.0 0.370 0.615 303 103 66.6 60.0 4.37 31.3 20.6 2.56 1.82"
W10X49 "14.4 10.0 10.0 0.340 0.560 272 93.4 60.4 54.6 4.35 28.3 18.7 2.54 1.39"
W10X45 "13.3 10.1 8.02 0.350 0.620 248 53.4 54.9 49.1 4.32 20.3 13.3 2.01 1.51"
W10X39 "11.5 9.92 7.99 0.315 0.530 209 45.0 46.8 42.1 4.27 17.2 11.3 1.98 0.976"
W10X33 "9.71 9.73 7.96 0.290 0.435 171 36.6 38.8 35.0 4.19 14.0 9.20 1.94 0.583"
W10X30 "8.84 10.5 5.81 0.300 0.510 170 16.7 36.6 32.4 4.38 8.84 5.75 1.37 0.622"
W10X26 "7.61 10.3 5.77 0.260 0.440 144 14.1 31.3 27.9 4.35 7.50 4.89 1.36 0.402"
W10X22 "6.49 10.2 5.75 0.240 0.360 118 11.4 26.0 23.2 4.27 6.10 3.97 1.33 0.239"
W10X19 "5.62 10.2 4.02 0.250 0.395 96.3 4.29 21.6 18.8 4.14 3.35 2.14 0.874 0.233"
W10X17 "4.99 10.1 4.01 0.240 0.330 81.9 3.56 18.7 16.2 4.05 2.80 1.78 0.845 0.156"
W10X15 "4.41 9.99 4.00 0.230 0.270 68.9 2.89 16.0 13.8 3.95 2.30 1.45 0.810 0.104"
W10X12 "3.54 9.87 3.96 0.190 0.210 53.8 2.18 12.6 10.9 3.90 1.74 1.10 0.785 0.0547"
W8X67 "19.7 9.00 8.28 0.570 0.935 272 88.6 70.1 60.4 3.72 32.7 21.4 2.12 5.05"
W8X58 "17.1 8.75 8.22 0.510 0.810 228 75.1 59.8 52.0 3.65 27.9 18.3 2.10 3.33"
W8X48 "14.1 8.50 8.11 0.400 0.685 184 60.9 49.0 43.2 3.61 22.9 15.0 2.08 1.96"
W8X40 "11.7 8.25 8.07 0.360 0.560 146 49.1 39.8 35.5 3.53 18.5 12.2 2.04 1.12"
W8X35 "10.3 8.12 8.02 0.310 0.495 127 42.6 34.7 31.2 3.51 16.1 10.6 2.03 0.769"
W8X31 "9.13 8.00 8.00 0.285 0.435 110 37.1 30.4 27.5 3.47 14.1 9.27 2.02 0.536"
W8X28 "8.25 8.06 6.54 0.285 0.465 98.0 21.7 27.2 24.3 3.45 10.1 6.63 1.62 0.537"
W8X24 "7.08 7.93 6.50 0.245 0.400 82.7 18.3 23.1 20.9 3.42 8.57 5.63 1.61 0.346"
W8X21 "6.16 8.28 5.27 0.250 0.400 75.3 9.77 20.4 18.2 3.49 5.69 3.71 1.26 0.282"
W8X18 "5.26 8.14 5.25 0.230 0.330 61.9 7.97 17.0 15.2 3.43 4.66 3.04 1.23 0.172"
W8X15 "4.44 8.11 4.02 0.245 0.315 48.0 3.41 13.6 11.8 3.29 2.67 1.70 0.876 0.137"
W8X13 "3.84 7.99 4.00 0.230 0.255 39.6 2.73 11.4 9.91 3.21 2.15 1.37 0.843 0.0871"
W8X10 "2.96 7.89 3.94 0.170 0.205 30.8 2.09 8.87 7.81 3.22 1.66 1.06 0.841 0.0426"
W6X25 "7.34 6.38 6.08 0.320 0.455 53.4 17.1 18.9 16.7 2.70 8.56 5.61 1.52 0.461"
W6X20 "5.87 6.20 6.02 0.260 0.365 41.4 13.3 15.0 13.4 2.66 6.72 4.41 1.50 0.240"
W6X15 "4.43 5.99 5.99 0.230 0.260 29.1 9.32 10.8 9.72 2.56 4.75 3.11 1.45 0.101"
W6X16 "4.74 6.28 4.03 0.260 0.405 32.1 4.43 11.7 10.2 2.60 3.39 2.20 0.967 0.223"
W6X12 "3.55 6.03 4.00 0.230 0.280 22.1 2.99 8.30 7.31 2.49 2.32 1.50 0.918 0.0903"
W6X9  "2.68 5.90 3.94 0.170 0.215 16.4 2.20 6.23 5.56 2.47 1.72 1.11 0.905 0.0405"
W6X8.5 "2.52 5.83 3.94 0.170 0.195 14.9 1.99 5.73 5.10 2.43 1.56 1.01 0.890 0.0333"
W5X19 "5.56 5.15 5.03 0.270 0.430 26.3 9.13 11.6 10.2 2.17 5.53 3.63 1.28 0.316"
W5X16 "4.71 5.01 5.00 0.240 0.360 21.4 7.51 9.63 8.55 2.13 4.58 3.00 1.26 0.192"
W4X13 "3.83 4.16 4.06 0.280 0.345 11.3 3.86 6.28 5.46 1.72 2.92 1.90 1.00 0.151"
FMK1  "50. 20.23 24.43 0.1 0.1 562 562 562 55.6"
FMK2  "50. 144.51 83.56 0.002 0.002 2248 2248 31.1 31.1"
FMK3  "50. 254.09 96.93 0.0005 0.0005 2248 2248 17.7 17.7"
FMK4  "50. 299.73 2137.71 0.0002 0.0002 1686 1686 11.3 11.3"
FMK5  "50. 674.40 134.73 0.00003 0.00003 1686 1686 5.0 5.0"
}




# AISC_Manual_Label "W A h b tdes Ix Zx Sx rx Iy Zy Sy ry J C
array set HSSection {
HSS20X12X5/8	"127.37	35.0	18.3	10.3	0.581	1880	230	188	7.33	851	162	142	4.93	1890 257"
HSS20X12X1/2	"103.3	28.3	18.6	10.6	0.465	1550	188	155	7.39	705	132	117	4.99	1540 209"
HSS20X12X3/8	"78.52	21.5	19.0	11.0	0.349	1200	144	120	7.45	547	102	91.1	5.04	1180 160"
HSS20X12X5/16	"65.87	18.1	19.1	11.1	0.291	1010	122	101	7.48	464	85.8	77.3	5.07	997  134"
HSS20X8X5/8	"110.36	30.3	18.3	6.26	0.581	1440	185	144	6.89	338	96.4	84.6	3.34	916  167"
HSS20X8X1/2	"89.68	24.6	18.6	6.60	0.465	1190	152	119	6.96	283	79.5	70.8	3.39	757  137"
HSS20X8X3/8	"68.31	18.7	19.0	6.95	0.349	926	117	92.6	7.03	222	61.5	55.6	3.44	586  105"
HSS20X8X5/16	"57.36	15.7	19.1	7.13	0.291	786	98.6	78.6	7.07	189	52.0	47.4	3.47	496  88.3"
HSS20X4X1/2	"76.07	20.9	18.6	2.60	0.465	838	115	83.8	6.33	58.7	34.0	29.3	1.68	195  63.8"
HSS20X4X3/8	"58.1	16.0	19.0	2.95	0.349	657	89.3	65.7	6.42	47.6	26.8	23.8	1.73	156  49.9"
HSS20X4X5/16	"48.86	13.4	19.1	3.13	0.291	560	75.6	56.0	6.46	41.2	22.9	20.6	1.75	134  42.4"
HSS20X4X1/4	"39.43	10.8	19.3	3.30	0.233	458	61.5	45.8	6.50	34.3	18.7	17.1	1.78	111  34.7"
HSS18X6X5/8	"93.34	25.7	16.3	4.26	0.581	923	135	103	6.00	158	61.0	52.7	2.48	462  109"
HSS18X6X1/2	"76.07	20.9	16.6	4.61	0.465	770	112	85.6	6.07	134	50.7	44.6	2.53	387  89.9"
HSS18X6X3/8	"58.1	16.0	17.0	4.95	0.349	602	86.4	66.9	6.15	106	39.5	35.5	2.58	302  69.5"
HSS18X6X5/16	"48.86	13.4	17.1	5.13	0.291	513	73.1	57.0	6.18	91.3	33.5	30.4	2.61	257  58.7"
HSS18X6X1/4	"39.43	10.8	17.3	5.30	0.233	419	59.4	46.5	6.22	75.1	27.3	25.0	2.63	210  47.7"
HSS16X16X5/8	"127.37	35.0	14.3	14.3	0.581	1370	200	171	6.25	1370	200	171	6.25	2170 276"
HSS16X16X1/2	"103.3	28.3	14.6	14.6	0.465	1130	164	141	6.31	1130	164	141	6.31	1770 224"
HSS16X16X3/8	"78.52	21.5	15.0	15.0	0.349	873	126	109	6.37	873	126	109	6.37	1350 171"
HSS16X16X5/16	"65.87	18.1	15.1	15.1	0.291	739	106	92.3	6.39	739	106	92.3	6.39	1140 144"
HSS16X12X5/8	"110.36	30.3	14.3	10.3	0.581	1090	165	136	6.00	700	135	117	4.80	1370 204"
HSS16X12X1/2	"89.68	24.6	14.6	10.6	0.465	904	135	113	6.06	581	111	96.8	4.86	1120 166"
HSS16X12X3/8	"68.31	18.7	15.0	11.0	0.349	702	104	87.7	6.12	452	85.5	75.3	4.91	862  127"
HSS16X12X5/16	"57.36	15.7	15.1	11.1	0.291	595	87.7	74.4	6.15	384	72.2	64.0	4.94	727  107"
HSS16X8X5/8	"93.34	25.7	14.3	6.26	0.581	815	129	102	5.64	274	79.2	68.6	3.27	681  132"
HSS16X8X1/2	"76.07	20.9	14.6	6.60	0.465	679	106	84.9	5.70	230	65.5	57.6	3.32	563  108"
HSS16X8X3/8	"58.1	16.0	15.0	6.95	0.349	531	82.1	66.3	5.77	181	50.8	45.3	3.37	436  83.4"
HSS16X8X5/16	"48.86	13.4	15.1	7.13	0.291	451	69.4	56.4	5.80	155	43.0	38.7	3.40	369  70.4"
HSS16X8X1/4	"39.43	10.8	15.3	7.30	0.233	368	56.4	46.1	5.83	127	35.0	31.7	3.42	300  57.0"
HSS16X4X5/8	"76.33	21.0	14.3	2.26	0.581	539	92.9	67.3	5.06	54.1	32.5	27.0	1.60	174  60.5"
HSS16X4X1/2	"62.46	17.2	14.6	2.60	0.465	455	77.3	56.9	5.15	47.0	27.4	23.5	1.65	150  50.7"
HSS16X4X3/8	"47.9	13.2	15.0	2.95	0.349	360	60.2	45.0	5.23	38.3	21.7	19.1	1.71	120  39.7"
HSS16X4X5/16	"40.35	11.1	15.1	3.13	0.291	308	51.1	38.5	5.27	33.2	18.5	16.6	1.73	103  33.8"
HSS16X4X1/4	"32.63	8.96	15.3	3.30	0.233	253	41.7	31.6	5.31	27.7	15.2	13.8	1.76	85.2 27.6"
HSS16X4X3/16	"24.73	6.76	15.5	3.48	0.174	193	31.7	24.2	5.35	21.5	11.7	10.8	1.78	65.5 21.1"
HSS14X14X5/8	"110.36	30.3	12.3	12.3	0.581	897	151	128	5.44	897	151	128	5.44	1430 208"
HSS14X14X1/2	"89.68	24.6	12.6	12.6	0.465	743	124	106	5.49	743	124	106	5.49	1170 170"
HSS14X14X3/8	"68.31	18.7	13.0	13.0	0.349	577	95.4	82.5	5.55	577	95.4	82.5	5.55	900  130"
HSS14X14X5/16	"57.36	15.7	13.1	13.1	0.291	490	80.5	69.9	5.58	490	80.5	69.9	5.58	759  109"
HSS14X10X5/8	"93.34	25.7	12.3	8.26	0.581	687	120	98.2	5.17	407	95.1	81.5	3.98	832  146"
HSS14X10X1/2	"76.07	20.9	12.6	8.60	0.465	573	98.8	81.8	5.23	341	78.5	68.1	4.04	685  120"
HSS14X10X3/8	"58.1	16.0	13.0	8.95	0.349	447	76.3	63.9	5.29	267	60.7	53.4	4.09	528  91.8"
HSS14X10X5/16	"48.86	13.4	13.1	9.13	0.291	380	64.6	54.3	5.32	227	51.4	45.5	4.12	446  77.4"
HSS14X10X1/4	"39.43	10.8	13.3	9.30	0.233	310	52.4	44.3	5.35	186	41.8	37.2	4.14	362  62.6"
HSS14X6X5/8	"76.33	21.0	12.3	4.26	0.581	478	88.7	68.3	4.77	124	48.4	41.2	2.43	334  83.7"
HSS14X6X1/2	"62.46	17.2	12.6	4.61	0.465	402	73.6	57.4	4.84	105	40.4	35.1	2.48	279  69.3"
HSS14X6X3/8	"47.9	13.2	13.0	4.95	0.349	317	57.3	45.3	4.91	84.1	31.6	28.0	2.53	219  53.7"
HSS14X6X5/16	"40.35	11.1	13.1	5.13	0.291	271	48.6	38.7	4.94	72.3	26.9	24.1	2.55	186  45.5"
HSS14X6X1/4	"32.63	8.96	13.3	5.30	0.233	222	39.6	31.7	4.98	59.6	22.0	19.9	2.58	152  36.9"
HSS14X6X3/16	"24.73	6.76	13.5	5.48	0.174	170	30.1	24.3	5.01	45.9	16.7	15.3	2.61	116  28.0"
HSS14X4X5/8	"67.82	18.7	12.3	2.26	0.581	373	73.1	53.3	4.47	47.2	28.5	23.6	1.59	148  52.6"
HSS14X4X1/2	"55.66	15.3	12.6	2.60	0.465	317	61.0	45.3	4.55	41.2	24.1	20.6	1.64	127  44.1"
HSS14X4X3/8	"42.79	11.8	13.0	2.95	0.349	252	47.8	36.0	4.63	33.6	19.1	16.8	1.69	102  34.6"
HSS14X4X5/16	"36.1	9.92	13.1	3.13	0.291	216	40.6	30.9	4.67	29.2	16.4	14.6	1.72	87.7 29.5"
HSS14X4X1/4	"29.23	8.03	13.3	3.30	0.233	178	33.2	25.4	4.71	24.4	13.5	12.2	1.74	72.4 24.1"
HSS14X4X3/16	"22.18	6.06	13.5	3.48	0.174	137	25.3	19.5	4.74	19.0	10.3	9.48	1.77	55.8 18.4"
HSS12X12X5/8	"93.34	25.7	10.3	10.3	0.581	548	109	91.4	4.62	548	109	91.4	4.62	885  151"
HSS12X12X1/2	"76.07	20.9	10.6	10.6	0.465	457	89.6	76.2	4.68	457	89.6	76.2	4.68	728  123"
HSS12X12X3/8	"58.1	16.0	11.0	11.0	0.349	357	69.2	59.5	4.73	357	69.2	59.5	4.73	561  94.6"
HSS12X12X5/16	"48.86	13.4	11.1	11.1	0.291	304	58.6	50.7	4.76	304	58.6	50.7	4.76	474  79.7"
HSS12X12X1/4	"39.43	10.8	11.3	11.3	0.233	248	47.6	41.4	4.79	248	47.6	41.4	4.79	384  64.5"
HSS12X12X3/16	"29.84	8.15	11.5	11.5	0.174	189	36.0	31.5	4.82	189	36.0	31.5	4.82	290  48.6"
HSS12X10X1/2	"69.27	19.0	10.6	8.60	0.465	395	78.8	65.9	4.56	298	69.6	59.7	3.96	545  102"
HSS12X10X3/8	"53.0	14.6	11.0	8.95	0.349	310	61.1	51.6	4.61	234	54.0	46.9	4.01	421  78.3"
HSS12X10X5/16	"44.6	12.2	11.1	9.13	0.291	264	51.7	44.0	4.64	200	45.7	40.0	4.04	356  66.1"
HSS12X10X1/4	"36.03	9.90	11.3	9.30	0.233	216	42.1	36.0	4.67	164	37.2	32.7	4.07	289  53.5"
HSS12X8X5/8	"76.33	21.0	10.3	6.26	0.581	397	82.1	66.1	4.34	210	61.9	52.5	3.16	454  97.7"
HSS12X8X1/2	"62.46	17.2	10.6	6.60	0.465	333	68.1	55.6	4.41	178	51.5	44.4	3.21	377  80.4"
HSS12X8X3/8	"47.9	13.2	11.0	6.95	0.349	262	53.0	43.7	4.47	140	40.1	35.1	3.27	293  62.1"
HSS12X8X5/16	"40.35	11.1	11.1	7.13	0.291	224	44.9	37.4	4.50	120	34.1	30.1	3.29	248  52.4"
HSS12X8X1/4	"32.63	8.96	11.3	7.30	0.233	184	36.6	30.6	4.53	98.8	27.8	24.7	3.32	202  42.5"
HSS12X8X3/16	"24.73	6.76	11.5	7.48	0.174	140	27.8	23.4	4.56	75.7	21.1	18.9	3.35	153  32.2"
HSS12X6X5/8	"67.82	18.7	10.3	4.26	0.581	321	68.8	53.4	4.14	107	42.1	35.5	2.39	271  71.1"
HSS12X6X1/2	"55.66	15.3	10.6	4.61	0.465	271	57.4	45.2	4.21	91.1	35.2	30.4	2.44	227  59.0"
HSS12X6X3/8	"42.79	11.8	11.0	4.95	0.349	215	44.8	35.9	4.28	72.9	27.7	24.3	2.49	178  45.8"
HSS12X6X5/16	"36.1	9.92	11.1	5.13	0.291	184	38.1	30.7	4.31	62.8	23.6	20.9	2.52	152  38.8"
HSS12X6X1/4	"29.23	8.03	11.3	5.30	0.233	151	31.1	25.2	4.34	51.9	19.3	17.3	2.54	124  31.6"
HSS12X6X3/16	"22.18	6.06	11.5	5.48	0.174	116	23.7	19.4	4.38	40.0	14.7	13.3	2.57	94.6 24.0"
HSS12X4X5/8	"59.32	16.4	10.3	2.26	0.581	245	55.5	40.8	3.87	40.4	24.5	20.2	1.57	122  44.6"
HSS12X4X1/2	"48.85	13.5	10.6	2.60	0.465	210	46.7	34.9	3.95	35.3	20.9	17.7	1.62	105  37.5"
HSS12X4X3/8	"37.69	10.4	11.0	2.95	0.349	168	36.7	28.0	4.02	28.9	16.6	14.5	1.67	84.1 29.5"
HSS12X4X5/16	"31.84	8.76	11.1	3.13	0.291	144	31.3	24.1	4.06	25.2	14.2	12.6	1.70	72.4 25.2"
HSS12X4X1/4	"25.82	7.10	11.3	3.30	0.233	119	25.6	19.9	4.10	21.0	11.7	10.5	1.72	59.8 20.6"
HSS12X4X3/16	"19.63	5.37	11.5	3.48	0.174	91.8	19.6	15.3	4.13	16.4	9.00	8.20	1.75	46.1 15.7"
HSS12X3-1/2X3/8	"36.41	10.0	11.0	2.45	0.349	156	34.7	26.0	3.94	21.3	14.0	12.2	1.46	64.7 25.5"
HSS12X3-1/2X5/16	"30.78	8.46	11.1	2.63	0.291	134	29.6	22.4	3.98	18.6	12.1	10.6	1.48 56.0 21.8"
HSS12X3X5/16		"29.72	8.17	11.1	2.13	0.291	124	27.9	20.7	3.90	13.1	10.0	8.73	1.27 41.3 18.4"
HSS12X3X1/4		"24.12	6.63	11.3	2.30	0.233	103	22.9	17.2	3.94	11.1	8.28	7.38	1.29 34.5 15.1"
HSS12X3X3/16		"18.35	5.02	11.5	2.48	0.174	79.6	17.5	13.3	3.98	8.72	6.40	5.81	1.32 26.8 11.6"
HSS12X2X5/16		"27.59	7.59	11.1	1.13	0.291	104	24.5	17.4	3.71	5.10	6.05	5.10	0.820	  17.6 11.6"
HSS12X2X1/4		"22.42	6.17	11.3	1.30	0.233	86.9	20.1	14.5	3.75	4.41	5.08	4.41	0.845	  15.1 9.64"
HSS12X2X3/16		"17.08	4.67	11.5	1.48	0.174	67.4	15.5	11.2	3.80	3.55	3.97	3.55	0.872	  12.0 7.49"
HSS10X10X5/8		"76.33	21.0	8.26	8.26	0.581	304	73.2	60.8	3.80	304	73.2	60.8	3.80	  498  102"
HSS10X10X1/2		"62.46	17.2	8.60	8.60	0.465	256	60.7	51.2	3.86	256	60.7	51.2	3.86	  412  84.2"
HSS10X10X3/8		"47.9	13.2	8.95	8.95	0.349	202	47.2	40.4	3.92	202	47.2	40.4	3.92	  320  64.8"
HSS10X10X5/16		"40.35	11.1	9.13	9.13	0.291	172	40.1	34.5	3.94	172	40.1	34.5	3.94	  271  54.8"
HSS10X10X1/4		"32.63	8.96	9.30	9.30	0.233	141	32.7	28.3	3.97	141	32.7	28.3	3.97	  220  44.4"
HSS10X10X3/16		"24.73	6.76	9.48	9.48	0.174	108	24.8	21.6	4.00	108	24.8	21.6	4.00	  167  33.6"
HSS10X8X5/8		"67.82	18.7	8.26	6.26	0.581	253	62.2	50.5	3.68	178	53.3	44.5	3.09	  346  80.4"
HSS10X8X1/2		"55.66	15.3	8.60	6.60	0.465	214	51.9	42.7	3.73	151	44.5	37.8	3.14	  288  66.4"
HSS10X8X3/8		"42.79	11.8	8.95	6.95	0.349	169	40.5	33.9	3.79	120	34.8	30.0	3.19	  224  51.4"
HSS10X8X5/16		"36.1	9.92	9.13	7.13	0.291	145	34.4	29.0	3.82	103	29.6	25.7	3.22	  190  43.5"
HSS10X8X1/4		"29.23	8.03	9.30	7.30	0.233	119	28.1	23.8	3.85	84.7	24.2	21.2	3.25	  155  35.3"
HSS10X8X3/16		"22.18	6.06	9.48	7.48	0.174	91.4	21.4	18.3	3.88	65.1	18.4	16.3	3.28	  118  26.7"
HSS10X6X5/8		"59.32	16.4	8.26	4.26	0.581	201	51.3	40.2	3.50	89.4	35.8	29.8	2.34	  209  58.6"
HSS10X6X1/2		"48.85	13.5	8.60	4.61	0.465	171	43.0	34.3	3.57	76.8	30.1	25.6	2.39	  176  48.7"
HSS10X6X3/8		"37.69	10.4	8.95	4.95	0.349	137	33.8	27.4	3.63	61.8	23.7	20.6	2.44	  139  37.9"
HSS10X6X5/16		"31.84	8.76	9.13	5.13	0.291	118	28.8	23.5	3.66	53.3	20.2	17.8	2.47	  118  32.2"
HSS10X6X1/4		"25.82	7.10	9.30	5.30	0.233	96.9	23.6	19.4	3.69	44.1	16.6	14.7	2.49	  96.7 26.2"
HSS10X6X3/16		"19.63	5.37	9.48	5.48	0.174	74.6	18.0	14.9	3.73	34.1	12.7	11.4	2.52	  73.8 19.9"
HSS10X5X3/8		"35.13	9.67	8.95	3.95	0.349	120	30.4	24.1	3.53	40.6	18.7	16.2	2.05	  100  31.2"
HSS10X5X5/16		"29.72	8.17	9.13	4.13	0.291	104	26.0	20.8	3.56	35.2	16.0	14.1	2.07	  86.0 26.5"
HSS10X5X1/4		"24.12	6.63	9.30	4.30	0.233	85.8	21.3	17.2	3.60	29.3	13.2	11.7	2.10	  70.7 21.6"
HSS10X5X3/16		"18.35	5.02	9.48	4.48	0.174	66.2	16.3	13.2	3.63	22.7	10.1	9.09	2.13	  54.1 16.5"
HSS10X4X5/8		"50.81	14.0	8.26	2.26	0.581	149	40.3	29.9	3.26	33.5	20.6	16.8	1.54	  95.7 36.7"
HSS10X4X1/2		"42.05	11.6	8.60	2.60	0.465	129	34.1	25.8	3.34	29.5	17.6	14.7	1.59	  82.6 31.0"
HSS10X4X3/8		"32.58	8.97	8.95	2.95	0.349	104	27.0	20.8	3.41	24.3	14.0	12.1	1.64	  66.5 24.4"
HSS10X4X5/16		"27.59	7.59	9.13	3.13	0.291	90.1	23.1	18.0	3.44	21.2	12.1	10.6	1.67	  57.3 20.9"
HSS10X4X1/4		"22.42	6.17	9.30	3.30	0.233	74.7	19.0	14.9	3.48	17.7	10.0	8.87	1.70	  47.4 17.1"
HSS10X4X3/16		"17.08	4.67	9.48	3.48	0.174	57.8	14.6	11.6	3.52	13.9	7.66	6.93	1.72	  36.5 13.1"
HSS10X4X1/8		"11.56	3.16	9.65	3.65	0.116	39.8	10.0	7.97	3.55	9.65	5.26	4.83	1.75	  25.1 8.90"
HSS10X3-1/2X1/2		"40.34	11.1	8.60	2.10	0.465	118	31.9	23.7	3.26	21.4	14.7	12.2	1.39	  63.2 26.5"
HSS10X3-1/2X3/8		"31.31	8.62	8.95	2.45	0.349	96.1	25.3	19.2	3.34	17.8	11.8	10.2	1.44	  51.5 21.1"
HSS10X3-1/2X5/16	"26.53	7.30	9.13	2.63	0.291	83.2	21.7	16.6	3.38	15.6	10.2	8.92	1.46	  44.6 18.0"
HSS10X3-1/2X1/4		"21.57	5.93	9.30	2.80	0.233	69.1	17.9	13.8	3.41	13.1	8.45	7.51	1.49	  37.0 14.8"
HSS10X3-1/2X3/16	"16.44	4.50	9.48	2.98	0.174	53.6	13.7	10.7	3.45	10.3	6.52	5.89	1.51	  28.6 11.4"
HSS10X3-1/2X1/8		"11.13	3.04	9.65	3.15	0.116	37.0	9.37	7.40	3.49	7.22	4.48	4.12	1.54	  19.8 7.75"
HSS10X3X3/8		"30.03	8.27	8.95	1.95	0.349	88.0	23.7	17.6	3.26	12.4	9.73	8.28	1.22	  37.8 17.7"
HSS10X3X5/16		"25.46	7.01	9.13	2.13	0.291	76.3	20.3	15.3	3.30	11.0	8.42	7.30	1.25	  33.0 15.2"
HSS10X3X1/4		"20.72	5.70	9.30	2.30	0.233	63.6	16.7	12.7	3.34	9.28	6.99	6.19	1.28	  27.6 12.5"
HSS10X3X3/16		"15.8	4.32	9.48	2.48	0.174	49.4	12.8	9.87	3.38	7.33	5.41	4.89	1.30	  21.5 9.64"
HSS10X3X1/8		"10.71	2.93	9.65	2.65	0.116	34.2	8.80	6.83	3.42	5.16	3.74	3.44	1.33	  14.9 6.61"
HSS10X2X3/8		"27.48	7.58	8.95	0.953	0.349	71.7	20.3	14.3	3.08	4.70	5.76	4.70	0.787	  15.9 11.0"
HSS10X2X5/16		"23.34	6.43	9.13	1.13	0.291	62.6	17.5	12.5	3.12	4.24	5.06	4.24	0.812	  14.2 9.56"
HSS10X2X1/4		"19.02	5.24	9.30	1.30	0.233	52.5	14.4	10.5	3.17	3.67	4.26	3.67	0.838	  12.2 7.99"
HSS10X2X3/16		"14.53	3.98	9.48	1.48	0.174	41.0	11.1	8.19	3.21	2.97	3.34	2.97	0.864	  9.74 6.22"
HSS10X2X1/8		"9.86	2.70	9.65	1.65	0.116	28.5	7.65	5.70	3.25	2.14	2.33	2.14	0.890	  6.90 4.31"
HSS9X9X5/8		"67.82	18.7	7.26	7.26	0.581	216	58.1	47.9	3.40	216	58.1	47.9	3.40	  356  81.6"
HSS9X9X1/2		"55.66	15.3	7.60	7.60	0.465	183	48.4	40.6	3.45	183	48.4	40.6	3.45	  296  67.4"
HSS9X9X3/8		"42.79	11.8	7.95	7.95	0.349	145	37.8	32.2	3.51	145	37.8	32.2	3.51	  231  52.1"
HSS9X9X5/16		"36.1	9.92	8.13	8.13	0.291	124	32.1	27.6	3.54	124	32.1	27.6	3.54	  196  44.0"
HSS9X9X1/4		"29.23	8.03	8.30	8.30	0.233	102	26.2	22.7	3.56	102	26.2	22.7	3.56	  159  35.8"
HSS9X9X3/16		"22.18	6.06	8.48	8.48	0.174	78.2	20.0	17.4	3.59	78.2	20.0	17.4	3.59	  121  27.1"
HSS9X9X1/8		"14.96	4.09	8.65	8.65	0.116	53.5	13.6	11.9	3.62	53.5	13.6	11.9	3.62	  82.0 18.3"
HSS9X7X5/8		"59.32	16.4	7.26	5.26	0.581	174	48.3	38.7	3.26	117	40.5	33.5	2.68	  235  62.0"
HSS9X7X1/2		"48.85	13.5	7.60	5.60	0.465	149	40.5	33.0	3.32	100	34.0	28.7	2.73	  197  51.5"
HSS9X7X3/8		"37.69	10.4	7.95	5.95	0.349	119	31.8	26.4	3.38	80.4	26.7	23.0	2.78	  154  40.0"
HSS9X7X5/16		"31.84	8.76	8.13	6.13	0.291	102	27.1	22.6	3.41	69.2	22.8	19.8	2.81	  131  33.9"
HSS9X7X1/4		"25.82	7.10	8.30	6.30	0.233	84.1	22.2	18.7	3.44	57.2	18.7	16.3	2.84	  107  27.6"
HSS9X7X3/16		"19.63	5.37	8.48	6.48	0.174	64.7	16.9	14.4	3.47	44.1	14.3	12.6	2.87	  81.7 20.9"
HSS9X5X5/8		"50.81	14.0	7.26	3.26	0.581	133	38.5	29.6	3.08	52.0	25.3	20.8	1.92	  128  42.5"
HSS9X5X1/2		"42.05	11.6	7.60	3.60	0.465	115	32.5	25.5	3.14	45.2	21.5	18.1	1.97	  109  35.6"
HSS9X5X3/8		"32.58	8.97	7.95	3.95	0.349	92.5	25.7	20.5	3.21	36.8	17.1	14.7	2.03	  86.9 27.9"
HSS9X5X5/16		"27.59	7.59	8.13	4.13	0.291	79.8	22.0	17.7	3.24	32.0	14.6	12.8	2.05	  74.4 23.8"
HSS9X5X1/4		"22.42	6.17	8.30	4.30	0.233	66.1	18.1	14.7	3.27	26.6	12.0	10.6	2.08	  61.2 19.4"
HSS9X5X3/16		"17.08	4.67	8.48	4.48	0.174	51.1	13.8	11.4	3.31	20.7	9.25	8.28	2.10	  46.9 14.8"
HSS9X3X1/2		"35.24	9.74	7.60	1.60	0.465	80.8	24.6	18.0	2.88	13.2	10.8	8.81	1.17	  40.0 19.7"
HSS9X3X3/8		"27.48	7.58	7.95	1.95	0.349	66.3	19.7	14.7	2.96	11.2	8.80	7.45	1.21	  33.1 15.8"
HSS9X3X5/16		"23.34	6.43	8.13	2.13	0.291	57.7	16.9	12.8	3.00	9.88	7.63	6.59	1.24	  28.9 13.6"
HSS9X3X1/4		"19.02	5.24	8.30	2.30	0.233	48.2	14.0	10.7	3.04	8.38	6.35	5.59	1.27	  24.2 11.3"
HSS9X3X3/16		"14.53	3.98	8.48	2.48	0.174	37.6	10.8	8.35	3.07	6.64	4.92	4.42	1.29	  18.9 8.66"
HSS8X8X5/8		"59.32	16.4	6.26	6.26	0.581	146	44.7	36.5	2.99	146	44.7	36.5	2.99	  244  63.2"
HSS8X8X1/2		"48.85	13.5	6.60	6.60	0.465	125	37.5	31.2	3.04	125	37.5	31.2	3.04	  204  52.4"
HSS8X8X3/8		"37.69	10.4	6.95	6.95	0.349	100	29.4	24.9	3.10	100	29.4	24.9	3.10	  160  40.7"
HSS8X8X5/16		"31.84	8.76	7.13	7.13	0.291	85.6	25.1	21.4	3.13	85.6	25.1	21.4	3.13	  136  34.5"
HSS8X8X1/4		"25.82	7.10	7.30	7.30	0.233	70.7	20.5	17.7	3.15	70.7	20.5	17.7	3.15	  111  28.1"
HSS8X8X3/16		"19.63	5.37	7.48	7.48	0.174	54.4	15.7	13.6	3.18	54.4	15.7	13.6	3.18	  84.5 21.3"
HSS8X8X1/8		"13.26	3.62	7.65	7.65	0.116	37.4	10.7	9.34	3.21	37.4	10.7	9.34	3.21	  57.3 14.4"
HSS8X6X5/8		"50.81	14.0	6.26	4.26	0.581	114	36.1	28.5	2.85	72.3	29.5	24.1	2.27	  150  46.0"
HSS8X6X1/2		"42.05	11.6	6.60	4.61	0.465	98.2	30.5	24.6	2.91	62.5	24.9	20.8	2.32	  127  38.4"
HSS8X6X3/8		"32.58	8.97	6.95	4.95	0.349	79.1	24.1	19.8	2.97	50.6	19.8	16.9	2.38	  100  30.0"
HSS8X6X5/16		"27.59	7.59	7.13	5.13	0.291	68.3	20.6	17.1	3.00	43.8	16.9	14.6	2.40	  85.8 25.5"
HSS8X6X1/4		"22.42	6.17	7.30	5.30	0.233	56.6	16.9	14.2	3.03	36.4	13.9	12.1	2.43	  70.3 20.8"
HSS8X6X3/16		"17.08	4.67	7.48	5.48	0.174	43.7	13.0	10.9	3.06	28.2	10.7	9.39	2.46	  53.7 15.8"
HSS8X4X5/8		"42.3	11.7	6.26	2.26	0.581	82.0	27.4	20.5	2.64	26.6	16.6	13.3	1.51	  70.3 28.7"
HSS8X4X1/2		"35.24	9.74	6.60	2.60	0.465	71.8	23.5	17.9	2.71	23.6	14.3	11.8	1.56	  61.1 24.4"
HSS8X4X3/8		"27.48	7.58	6.95	2.95	0.349	58.7	18.8	14.7	2.78	19.6	11.5	9.80	1.61	  49.3 19.3"
HSS8X4X5/16		"23.34	6.43	7.13	3.13	0.291	51.0	16.1	12.8	2.82	17.2	9.91	8.58	1.63	  42.6 16.5"
HSS8X4X1/4		"19.02	5.24	7.30	3.30	0.233	42.5	13.3	10.6	2.85	14.4	8.20	7.21	1.66	  35.3 13.6"
HSS8X4X3/16		"14.53	3.98	7.48	3.48	0.174	33.1	10.2	8.27	2.88	11.3	6.33	5.65	1.69	  27.2 10.4"
HSS8X4X1/8		"9.86	2.70	7.65	3.65	0.116	22.9	7.02	5.73	2.92	7.90	4.36	3.95	1.71	  18.7 7.10"
HSS8X3X1/2		"31.84	8.81	6.60	1.60	0.465	58.6	20.0	14.6	2.58	11.7	9.64	7.81	1.15	  34.3 17.4"
HSS8X3X3/8		"24.93	6.88	6.95	1.95	0.349	48.5	16.1	12.1	2.65	10.0	7.88	6.63	1.20	  28.5 14.0"
HSS8X3X5/16		"21.21	5.85	7.13	2.13	0.291	42.4	13.9	10.6	2.69	8.81	6.84	5.87	1.23	  24.9 12.1"
HSS8X3X1/4		"17.32	4.77	7.30	2.30	0.233	35.5	11.5	8.88	2.73	7.49	5.70	4.99	1.25	  20.8 10.0"
HSS8X3X3/16		"13.25	3.63	7.48	2.48	0.174	27.8	8.87	6.94	2.77	5.94	4.43	3.96	1.28	  16.2 7.68"
HSS8X3X1/8		"9.01	2.46	7.65	2.65	0.116	19.3	6.11	4.83	2.80	4.20	3.07	2.80	1.31	  11.3 5.27"
HSS8X2X3/8		"22.37	6.18	6.95	0.953	0.349	38.2	13.4	9.56	2.49	3.73	4.61	3.73	0.777	  12.1 8.65"
HSS8X2X5/16		"19.08	5.26	7.13	1.13	0.291	33.7	11.6	8.43	2.53	3.38	4.06	3.38	0.802	  10.9 7.57"
HSS8X2X1/4		"15.62	4.30	7.30	1.30	0.233	28.5	9.68	7.12	2.57	2.94	3.43	2.94	0.827	  9.36 6.35"
HSS8X2X3/16		"11.97	3.28	7.48	1.48	0.174	22.4	7.51	5.61	2.61	2.39	2.70	2.39	0.853	  7.48 4.95"
HSS8X2X1/8		"8.16	2.23	7.65	1.65	0.116	15.7	5.19	3.93	2.65	1.72	1.90	1.72	0.879	  5.30 3.44"
HSS7X7X5/8		"50.81	14.0	5.26	5.26	0.581	93.4	33.1	26.7	2.58	93.4	33.1	26.7	2.58	  158  47.1"
HSS7X7X1/2		"42.05	11.6	5.60	5.60	0.465	80.5	27.9	23.0	2.63	80.5	27.9	23.0	2.63	  133  39.3"
HSS7X7X3/8		"32.58	8.97	5.95	5.95	0.349	65.0	22.1	18.6	2.69	65.0	22.1	18.6	2.69	  105  30.7"
HSS7X7X5/16		"27.59	7.59	6.13	6.13	0.291	56.1	18.9	16.0	2.72	56.1	18.9	16.0	2.72	  89.7 26.1"
HSS7X7X1/4		"22.42	6.17	6.30	6.30	0.233	46.5	15.5	13.3	2.75	46.5	15.5	13.3	2.75	  73.5 21.3"
HSS7X7X3/16		"17.08	4.67	6.48	6.48	0.174	36.0	11.9	10.3	2.77	36.0	11.9	10.3	2.77	  56.1 16.2"
HSS7X7X1/8		"11.56	3.16	6.65	6.65	0.116	24.8	8.13	7.09	2.80	24.8	8.13	7.09	2.80	  38.2 11.0"
HSS7X5X1/2		"35.24	9.74	5.60	3.60	0.465	60.6	21.9	17.3	2.50	35.6	17.3	14.2	1.91	  75.8 27.2"
HSS7X5X3/8		"27.48	7.58	5.95	3.95	0.349	49.5	17.5	14.1	2.56	29.3	13.8	11.7	1.97	  60.6 21.4"
HSS7X5X5/16		"23.34	6.43	6.13	4.13	0.291	43.0	15.0	12.3	2.59	25.5	11.9	10.2	1.99	  52.1 18.3"
HSS7X5X1/4		"19.02	5.24	6.30	4.30	0.233	35.9	12.4	10.2	2.62	21.3	9.83	8.53	2.02	  42.9 15.0"
HSS7X5X3/16		"14.53	3.98	6.48	4.48	0.174	27.9	9.52	7.96	2.65	16.6	7.57	6.65	2.05	  32.9 11.4"
HSS7X5X1/8		"9.86	2.70	6.65	4.65	0.116	19.3	6.53	5.52	2.68	11.6	5.20	4.63	2.07	  22.5 7.79"
HSS7X4X1/2		"31.84	8.81	5.60	2.60	0.465	50.7	18.8	14.5	2.40	20.7	12.6	10.4	1.53	  50.5 21.1"
HSS7X4X3/8		"24.93	6.88	5.95	2.95	0.349	41.8	15.1	11.9	2.46	17.3	10.2	8.63	1.58	  41.0 16.8"
HSS7X4X5/16		"21.21	5.85	6.13	3.13	0.291	36.5	13.1	10.4	2.50	15.2	8.83	7.58	1.61	  35.4 14.4"
HSS7X4X1/4		"17.32	4.77	6.30	3.30	0.233	30.5	10.8	8.72	2.53	12.8	7.33	6.38	1.64	  29.3 11.8"
HSS7X4X3/16		"13.25	3.63	6.48	3.48	0.174	23.8	8.33	6.81	2.56	10.0	5.67	5.02	1.66	  22.7 9.07"
HSS7X4X1/8		"9.01	2.46	6.65	3.65	0.116	16.6	5.73	4.73	2.59	7.03	3.91	3.51	1.69	  15.6 6.20"
HSS7X3X1/2		"28.43	7.88	5.60	1.60	0.465	40.7	15.8	11.6	2.27	10.2	8.46	6.80	1.14	  28.6 15.0"
HSS7X3X3/8		"22.37	6.18	5.95	1.95	0.349	34.1	12.8	9.73	2.35	8.71	6.95	5.81	1.19	  23.9 12.1"
HSS7X3X5/16		"19.08	5.26	6.13	2.13	0.291	29.9	11.1	8.54	2.38	7.74	6.05	5.16	1.21	  20.9 10.5"
HSS7X3X1/4		"15.62	4.30	6.30	2.30	0.233	25.2	9.22	7.19	2.42	6.60	5.06	4.40	1.24	  17.5 8.68"
HSS7X3X3/16		"11.97	3.28	6.48	2.48	0.174	19.8	7.14	5.65	2.45	5.24	3.94	3.50	1.26	  13.7 6.69"
HSS7X3X1/8		"8.16	2.23	6.65	2.65	0.116	13.8	4.93	3.95	2.49	3.71	2.73	2.48	1.29	  9.48 4.60"
HSS7X2X1/4		"13.91	3.84	6.30	1.30	0.233	19.8	7.64	5.67	2.27	2.58	3.02	2.58	0.819	  7.95 5.52"
HSS7X2X3/16		"10.7	2.93	6.48	1.48	0.174	15.7	5.95	4.49	2.31	2.10	2.39	2.10	0.845	  6.35 4.32"
HSS7X2X1/8		"7.31	2.00	6.65	1.65	0.116	11.1	4.13	3.16	2.35	1.52	1.68	1.52	0.871	  4.51 3.00"
HSS6X6X5/8		"42.3	11.7	4.26	4.26	0.581	55.2	23.2	18.4	2.17	55.2	23.2	18.4	2.17	  94.9 33.4"
HSS6X6X1/2		"35.24	9.74	4.61	4.61	0.465	48.3	19.8	16.1	2.23	48.3	19.8	16.1	2.23	  81.1 28.1"
HSS6X6X3/8		"27.48	7.58	4.95	4.95	0.349	39.5	15.8	13.2	2.28	39.5	15.8	13.2	2.28	  64.6 22.1"
HSS6X6X5/16		"23.34	6.43	5.13	5.13	0.291	34.3	13.6	11.4	2.31	34.3	13.6	11.4	2.31	  55.4 18.9"
HSS6X6X1/4		"19.02	5.24	5.30	5.30	0.233	28.6	11.2	9.54	2.34	28.6	11.2	9.54	2.34	  45.6 15.4"
HSS6X6X3/16		"14.53	3.98	5.48	5.48	0.174	22.3	8.63	7.42	2.37	22.3	8.63	7.42	2.37	  35.0 11.8"
HSS6X6X1/8		"9.86	2.70	5.65	5.65	0.116	15.5	5.92	5.15	2.39	15.5	5.92	5.15	2.39	  23.9 8.03"
HSS6X5X1/2		"31.84	8.81	4.61	3.60	0.465	41.1	17.2	13.7	2.16	30.8	15.2	12.3	1.87	  59.8 23.0"
HSS6X5X3/8		"24.93	6.88	4.95	3.95	0.349	33.9	13.8	11.3	2.22	25.5	12.2	10.2	1.92	  48.1 18.2"
HSS6X5X5/16		"21.21	5.85	5.13	4.13	0.291	29.6	11.9	9.85	2.25	22.3	10.5	8.91	1.95	  41.4 15.6"
HSS6X5X1/4		"17.32	4.77	5.30	4.30	0.233	24.7	9.87	8.25	2.28	18.7	8.72	7.47	1.98	  34.2 12.8"
HSS6X5X3/16		"13.25	3.63	5.48	4.48	0.174	19.3	7.62	6.44	2.31	14.6	6.73	5.84	2.01	  26.3 9.76"
HSS6X5X1/8		"9.01	2.46	5.65	4.65	0.116	13.4	5.24	4.48	2.34	10.2	4.63	4.07	2.03	  18.0 6.66"
HSS6X4X1/2		"28.43	7.88	4.61	2.60	0.465	34.0	14.6	11.3	2.08	17.8	11.0	8.89	1.50	  40.3 17.8"
HSS6X4X3/8		"22.37	6.18	4.95	2.95	0.349	28.3	11.9	9.43	2.14	14.9	8.94	7.47	1.55	  32.8 14.2"
HSS6X4X5/16		"19.08	5.26	5.13	3.13	0.291	24.8	10.3	8.27	2.17	13.2	7.75	6.58	1.58	  28.4 12.2"
HSS6X4X1/4		"15.62	4.30	5.30	3.30	0.233	20.9	8.53	6.96	2.20	11.1	6.45	5.56	1.61	  23.6 10.1"
HSS6X4X3/16		"11.97	3.28	5.48	3.48	0.174	16.4	6.60	5.46	2.23	8.76	5.00	4.38	1.63	  18.2 7.74"
HSS6X4X1/8		"8.16	2.23	5.65	3.65	0.116	11.4	4.56	3.81	2.26	6.15	3.46	3.08	1.66	  12.6 5.30"
HSS6X3X1/2		"25.03	6.95	4.61	1.60	0.465	26.8	12.1	8.95	1.97	8.69	7.28	5.79	1.12	  23.1 12.7"
HSS6X3X3/8		"19.82	5.48	4.95	1.95	0.349	22.7	9.90	7.57	2.04	7.48	6.03	4.99	1.17	  19.3 10.3"
HSS6X3X5/16		"16.96	4.68	5.13	2.13	0.291	20.1	8.61	6.69	2.07	6.67	5.27	4.45	1.19	  16.9 8.91"
HSS6X3X1/4		"13.91	3.84	5.30	2.30	0.233	17.0	7.19	5.66	2.10	5.70	4.41	3.80	1.22	  14.2 7.39"
HSS6X3X3/16		"10.7	2.93	5.48	2.48	0.174	13.4	5.59	4.47	2.14	4.55	3.45	3.03	1.25	  11.1 5.71"
HSS6X3X1/8		"7.31	2.00	5.65	2.65	0.116	9.43	3.87	3.14	2.17	3.23	2.40	2.15	1.27	  7.73 3.93"
HSS6X2X3/8		"17.27	4.78	4.95	0.953	0.349	17.1	7.93	5.71	1.89	2.77	3.46	2.77	0.760	  8.42 6.35"
HSS6X2X5/16		"14.83	4.10	5.13	1.13	0.291	15.3	6.95	5.11	1.93	2.52	3.07	2.52	0.785	  7.60 5.58"
HSS6X2X1/4		"12.21	3.37	5.30	1.30	0.233	13.1	5.84	4.37	1.97	2.21	2.61	2.21	0.810	  6.55 4.70"
HSS6X2X3/16		"9.42	2.58	5.48	1.48	0.174	10.5	4.58	3.49	2.01	1.80	2.07	1.80	0.836	  5.24 3.68"
HSS6X2X1/8		"6.46	1.77	5.65	1.65	0.116	7.42	3.19	2.47	2.05	1.31	1.46	1.31	0.861	  3.72 2.57"
HSS5-1/2X5-1/2X3/8	"24.93	6.88	4.45	4.45	0.349	29.7	13.1	10.8	2.08	29.7	13.1	10.8	2.08	  49.0 18.4"
HSS5-1/2X5-1/2X5/16	"21.21	5.85	4.63	4.63	0.291	25.9	11.3	9.43	2.11	25.9	11.3	9.43	2.11	  42.2 15.7"
HSS5-1/2X5-1/2X1/4	"17.32	4.77	4.80	4.80	0.233	21.7	9.32	7.90	2.13	21.7	9.32	7.90	2.13	  34.8 12.9"
HSS5-1/2X5-1/2X3/16	"13.25	3.63	4.98	4.98	0.174	17.0	7.19	6.17	2.16	17.0	7.19	6.17	2.16	  26.7 9.85"
HSS5-1/2X5-1/2X1/8	"9.01	2.46	5.15	5.15	0.116	11.8	4.95	4.30	2.19	11.8	4.95	4.30	2.19	  18.3 6.72"
HSS5X5X1/2		"28.43	7.88	3.60	3.60	0.465	26.0	13.1	10.4	1.82	26.0	13.1	10.4	1.82	  44.6 18.7"
HSS5X5X3/8		"22.37	6.18	3.95	3.95	0.349	21.7	10.6	8.68	1.87	21.7	10.6	8.68	1.87	  36.1 14.9"
HSS5X5X5/16		"19.08	5.26	4.13	4.13	0.291	19.0	9.16	7.62	1.90	19.0	9.16	7.62	1.90	  31.2 12.8"
HSS5X5X1/4		"15.62	4.30	4.30	4.30	0.233	16.0	7.61	6.41	1.93	16.0	7.61	6.41	1.93	  25.8 10.5"
HSS5X5X3/16		"11.97	3.28	4.48	4.48	0.174	12.6	5.89	5.03	1.96	12.6	5.89	5.03	1.96	  19.9 8.08"
HSS5X5X1/8		"8.16	2.23	4.65	4.65	0.116	8.80	4.07	3.52	1.99	8.80	4.07	3.52	1.99	  13.7 5.53"
HSS5X4X1/2		"25.03	6.95	3.60	2.60	0.465	21.2	10.9	8.49	1.75	14.9	9.35	7.43	1.46	  30.3 14.5"
HSS5X4X3/8		"19.82	5.48	3.95	2.95	0.349	17.9	8.96	7.17	1.81	12.6	7.67	6.30	1.52	  24.9 11.7"
HSS5X4X5/16		"16.96	4.68	4.13	3.13	0.291	15.8	7.79	6.32	1.84	11.1	6.67	5.57	1.54	  21.7 10.1"
HSS5X4X1/4		"13.91	3.84	4.30	3.30	0.233	13.4	6.49	5.35	1.87	9.46	5.57	4.73	1.57	  18.0 8.32"
HSS5X4X3/16		"10.7	2.93	4.48	3.48	0.174	10.6	5.05	4.22	1.90	7.48	4.34	3.74	1.60	  14.0 6.41"
HSS5X4X1/8		"7.31	2.00	4.65	3.65	0.116	7.42	3.50	2.97	1.93	5.27	3.01	2.64	1.62	  9.66 4.39"
HSS5X3X1/2		"21.63	6.02	3.60	1.60	0.465	16.4	8.83	6.57	1.65	7.18	6.10	4.78	1.09	  17.6 10.3"
HSS5X3X3/8		"17.27	4.78	3.95	1.95	0.349	14.1	7.34	5.65	1.72	6.25	5.10	4.16	1.14	  14.9 8.44"
HSS5X3X5/16		"14.83	4.10	4.13	2.13	0.291	12.6	6.42	5.03	1.75	5.60	4.48	3.73	1.17	  13.1 7.33"
HSS5X3X1/4		"12.21	3.37	4.30	2.30	0.233	10.7	5.38	4.29	1.78	4.81	3.77	3.21	1.19	  11.0 6.10"
HSS5X3X3/16		"9.42	2.58	4.48	2.48	0.174	8.53	4.21	3.41	1.82	3.85	2.96	2.57	1.22	  8.64 4.73"
HSS5X3X1/8		"6.46	1.77	4.65	2.65	0.116	6.03	2.93	2.41	1.85	2.75	2.07	1.83	1.25	  6.02 3.26"
HSS5X2-1/2X1/4		"11.36	3.14	4.30	1.80	0.233	9.40	4.83	3.76	1.73	3.13	2.95	2.50	0.999	  7.93 4.99"
HSS5X2-1/2X3/16		"8.78	2.41	4.48	1.98	0.174	7.51	3.79	3.01	1.77	2.53	2.33	2.03	1.02	  6.26 3.89"
HSS5X2-1/2X1/8		"6.03	1.65	4.65	2.15	0.116	5.34	2.65	2.14	1.80	1.82	1.64	1.46	1.05	  4.40 2.70"
HSS5X2X3/8		"14.72	4.09	3.95	0.953	0.349	10.4	5.71	4.14	1.59	2.28	2.88	2.28	0.748	  6.61 5.20"
HSS5X2X5/16		"12.7	3.52	4.13	1.13	0.291	9.35	5.05	3.74	1.63	2.10	2.57	2.10	0.772	  5.99 4.59"
HSS5X2X1/4		"10.51	2.91	4.30	1.30	0.233	8.08	4.27	3.23	1.67	1.84	2.20	1.84	0.797	  5.17 3.88"
HSS5X2X3/16		"8.15	2.24	4.48	1.48	0.174	6.50	3.37	2.60	1.70	1.51	1.75	1.51	0.823	  4.15 3.05"
HSS5X2X1/8		"5.61	1.54	4.65	1.65	0.116	4.65	2.37	1.86	1.74	1.10	1.24	1.10	0.848	  2.95 2.13"
HSS4-1/2X4-1/2X1/2	"25.03	6.95	3.10	3.10	0.465	18.1	10.2	8.03	1.61	18.1	10.2	8.03	1.61	  31.3 14.8"
HSS4-1/2X4-1/2X3/8	"19.82	5.48	3.45	3.45	0.349	15.3	8.36	6.79	1.67	15.3	8.36	6.79	1.67	  25.7 11.9"
HSS4-1/2X4-1/2X5/16	"16.96	4.68	3.63	3.63	0.291	13.5	7.27	6.00	1.70	13.5	7.27	6.00	1.70	  22.3 10.2"
HSS4-1/2X4-1/2X1/4	"13.91	3.84	3.80	3.80	0.233	11.4	6.06	5.08	1.73	11.4	6.06	5.08	1.73	  18.5 8.44"
HSS4-1/2X4-1/2X3/16	"10.7	2.93	3.98	3.98	0.174	9.02	4.71	4.01	1.75	9.02	4.71	4.01	1.75	  14.4 6.49"
HSS4-1/2X4-1/2X1/8	"7.31	2.00	4.15	4.15	0.116	6.35	3.27	2.82	1.78	6.35	3.27	2.82	1.78	  9.92 4.45"
HSS4X4X1/2		"21.63	6.02	2.60	2.60	0.465	11.9	7.70	5.97	1.41	11.9	7.70	5.97	1.41	  21.0 11.2"
HSS4X4X3/8		"17.27	4.78	2.95	2.95	0.349	10.3	6.39	5.13	1.47	10.3	6.39	5.13	1.47	  17.5 9.14"
HSS4X4X5/16		"14.83	4.10	3.13	3.13	0.291	9.14	5.59	4.57	1.49	9.14	5.59	4.57	1.49	  15.3 7.91"
HSS4X4X1/4		"12.21	3.37	3.30	3.30	0.233	7.80	4.69	3.90	1.52	7.80	4.69	3.90	1.52	  12.8 6.56"
HSS4X4X3/16		"9.42	2.58	3.48	3.48	0.174	6.21	3.67	3.10	1.55	6.21	3.67	3.10	1.55	  10.0 5.07"
HSS4X4X1/8		"6.46	1.77	3.65	3.65	0.116	4.40	2.56	2.20	1.58	4.40	2.56	2.20	1.58	  6.91 3.49"
HSS4X3X3/8		"14.72	4.09	2.95	1.95	0.349	7.93	5.12	3.97	1.39	5.01	4.18	3.34	1.11	  10.6 6.59"
HSS4X3X5/16		"12.7	3.52	3.13	2.13	0.291	7.14	4.51	3.57	1.42	4.52	3.69	3.02	1.13	  9.41 5.75"
HSS4X3X1/4		"10.51	2.91	3.30	2.30	0.233	6.15	3.81	3.07	1.45	3.91	3.12	2.61	1.16	  7.96 4.81"
HSS4X3X3/16		"8.15	2.24	3.48	2.48	0.174	4.93	3.00	2.47	1.49	3.16	2.46	2.10	1.19	  6.26 3.74"
HSS4X3X1/8		"5.61	1.54	3.65	2.65	0.116	3.52	2.11	1.76	1.52	2.27	1.73	1.51	1.21	  4.38 2.59"
HSS4X2-1/2X3/8		"13.44	3.74	2.95	1.45	0.349	6.77	4.48	3.38	1.35	3.17	3.20	2.54	0.922	  7.57 5.32"
HSS4X2-1/2X5/16		"11.64	3.23	3.13	1.63	0.291	6.13	3.97	3.07	1.38	2.89	2.85	2.32	0.947	  6.77 4.67"
HSS4X2-1/2X1/4		"9.66	2.67	3.30	1.80	0.233	5.32	3.38	2.66	1.41	2.53	2.43	2.02	0.973	  5.78 3.93"
HSS4X2-1/2X3/16		"7.51	2.06	3.48	1.98	0.174	4.30	2.67	2.15	1.44	2.06	1.93	1.65	0.999	  4.59 3.08"
HSS4X2-1/2X1/8		"5.18	1.42	3.65	2.15	0.116	3.09	1.88	1.54	1.47	1.49	1.36	1.19	1.03	  3.23 2.14"
HSS4X2X3/8		"12.17	3.39	2.95	0.953	0.349	5.60	3.84	2.80	1.29	1.80	2.31	1.80	0.729	  4.83 4.04"
HSS4X2X5/16		"10.58	2.94	3.13	1.13	0.291	5.13	3.43	2.56	1.32	1.67	2.08	1.67	0.754	  4.40 3.59"
HSS4X2X1/4		"8.81	2.44	3.30	1.30	0.233	4.49	2.94	2.25	1.36	1.48	1.79	1.48	0.779	  3.82 3.05"
HSS4X2X3/16		"6.87	1.89	3.48	1.48	0.174	3.66	2.34	1.83	1.39	1.22	1.43	1.22	0.804	  3.08 2.41"
HSS4X2X1/8		"4.75	1.30	3.65	1.65	0.116	2.65	1.66	1.32	1.43	0.898	1.02	0.898	0.830	  2.20 1.69"
HSS3-1/2X3-1/2X3/8	"14.72	4.09	2.45	2.45	0.349	6.49	4.69	3.71	1.26	6.49	4.69	3.71	1.26	  11.2 6.77"
HSS3-1/2X3-1/2X5/16	"12.7	3.52	2.63	2.63	0.291	5.84	4.14	3.34	1.29	5.84	4.14	3.34	1.29	  9.89 5.90"
HSS3-1/2X3-1/2X1/4	"10.51	2.91	2.80	2.80	0.233	5.04	3.50	2.88	1.32	5.04	3.50	2.88	1.32	  8.35 4.92"
HSS3-1/2X3-1/2X3/16	"8.15	2.24	2.98	2.98	0.174	4.05	2.76	2.31	1.35	4.05	2.76	2.31	1.35	  6.56 3.83"
HSS3-1/2X3-1/2X1/8	"5.61	1.54	3.15	3.15	0.116	2.90	1.93	1.66	1.37	2.90	1.93	1.66	1.37	  4.58 2.65"
HSS3-1/2X2-1/2X3/8	"12.17	3.39	2.45	1.45	0.349	4.75	3.59	2.72	1.18	2.77	2.82	2.21	0.904	  6.16 4.57"
HSS3-1/2X2-1/2X5/16	"10.58	2.94	2.63	1.63	0.291	4.34	3.20	2.48	1.22	2.54	2.52	2.03	0.930	  5.53 4.03"
HSS3-1/2X2-1/2X1/4	"8.81	2.44	2.80	1.80	0.233	3.79	2.74	2.17	1.25	2.23	2.16	1.78	0.956	  4.75 3.40"
HSS3-1/2X2-1/2X3/16	"6.87	1.89	2.98	1.98	0.174	3.09	2.18	1.76	1.28	1.82	1.72	1.46	0.983	  3.78 2.67"
HSS3-1/2X2-1/2X1/8	"4.75	1.30	3.15	2.15	0.116	2.23	1.54	1.28	1.31	1.33	1.22	1.06	1.01	  2.67 1.87"
HSS3-1/2X2X1/4		"7.96	2.21	2.80	1.30	0.233	3.17	2.36	1.81	1.20	1.30	1.58	1.30	0.766	  3.16 2.64"
HSS3-1/2X2X3/16		"6.23	1.71	2.98	1.48	0.174	2.61	1.89	1.49	1.23	1.08	1.27	1.08	0.792	  2.55 2.09"
HSS3-1/2X2X1/8		"4.33	1.19	3.15	1.65	0.116	1.90	1.34	1.09	1.27	0.795	0.912	0.795	0.818	  1.83 1.47"
HSS3-1/2X1-1/2X1/4	"7.11	1.97	2.80	0.801	0.233	2.55	1.98	1.46	1.14	0.638	1.06	0.851	0.569	  1.79 1.88"
HSS3-1/2X1-1/2X3/16	"5.59	1.54	2.98	0.978	0.174	2.12	1.60	1.21	1.17	0.544	0.867	0.725	0.594	  1.49 1.51"
HSS3-1/2X1-1/2X1/8	"3.90	1.07	3.15	1.15	0.116	1.57	1.15	0.896	1.21	0.411	0.630	0.548	0.619	  1.09 1.08"
HSS3X3X3/8		"12.17	3.39	1.95	1.95	0.349	3.78	3.25	2.52	1.06	3.78	3.25	2.52	1.06	  6.64 4.74"
HSS3X3X5/16		"10.58	2.94	2.13	2.13	0.291	3.45	2.90	2.30	1.08	3.45	2.90	2.30	1.08	  5.94 4.18"
HSS3X3X1/4		"8.81	2.44	2.30	2.30	0.233	3.02	2.48	2.01	1.11	3.02	2.48	2.01	1.11	  5.08 3.52"
HSS3X3X3/16		"6.87	1.89	2.48	2.48	0.174	2.46	1.97	1.64	1.14	2.46	1.97	1.64	1.14	  4.03 2.76"
HSS3X3X1/8		"4.75	1.30	2.65	2.65	0.116	1.78	1.40	1.19	1.17	1.78	1.40	1.19	1.17	  2.84 1.92"
HSS3X2-1/2X5/16		"9.51	2.64	2.13	1.63	0.291	2.92	2.51	1.94	1.05	2.18	2.20	1.74	0.908	  4.34 3.39"
HSS3X2-1/2X1/4		"7.96	2.21	2.30	1.80	0.233	2.57	2.16	1.72	1.08	1.93	1.90	1.54	0.935	  3.74 2.87"
HSS3X2-1/2X3/16		"6.23	1.71	2.48	1.98	0.174	2.11	1.73	1.41	1.11	1.59	1.52	1.27	0.963	  3.00 2.27"
HSS3X2-1/2X1/8		"4.33	1.19	2.65	2.15	0.116	1.54	1.23	1.03	1.14	1.16	1.09	0.931	0.990	  2.13 1.59"
HSS3X2X5/16		"8.45	2.35	2.13	1.13	0.291	2.38	2.11	1.59	1.01	1.24	1.58	1.24	0.725	  2.87 2.60"
HSS3X2X1/4		"7.11	1.97	2.30	1.30	0.233	2.13	1.83	1.42	1.04	1.11	1.38	1.11	0.751	  2.52 2.23"
HSS3X2X3/16		"5.59	1.54	2.48	1.48	0.174	1.77	1.48	1.18	1.07	0.932	1.12	0.932	0.778	  2.05 1.78"
HSS3X2X1/8		"3.90	1.07	2.65	1.65	0.116	1.30	1.06	0.867	1.10	0.692	0.803	0.692	0.804	  1.47 1.25"
HSS3X1-1/2X1/4		"6.26	1.74	2.30	0.801	0.233	1.68	1.51	1.12	0.982	0.543	0.911	0.725	0.559	  1.44 1.58"
HSS3X1-1/2X3/16		"4.96	1.37	2.48	0.978	0.174	1.42	1.24	0.945	1.02	0.467	0.752	0.622	0.584	  1.21 1.28"
HSS3X1-1/2X1/8		"3.48	0.956	2.65	1.15	0.116	1.06	0.895	0.706	1.05	0.355	0.550	0.474	0.610	  0.886 0.920"
HSS3X1X3/16		"4.32	1.19	2.48	0.478	0.174	1.07	0.989	0.713	0.947	0.173	0.432	0.345	0.380	  0.526	0.792"
HSS3X1X1/8		"3.05	0.840	2.65	0.652	0.116	0.817	0.728	0.545	0.987	0.138	0.325	0.276	0.405	  0.408	0.585"
HSS2-1/2X2-1/2X5/16	"8.45	2.35	1.63	1.63	0.291	1.82	1.88	1.46	0.880	1.82	1.88	1.46	0.880	  3.20	2.74"
HSS2-1/2X2-1/2X1/4	"7.11	1.97	1.80	1.80	0.233	1.63	1.63	1.30	0.908	1.63	1.63	1.30	0.908	  2.79	2.35"
HSS2-1/2X2-1/2X3/16	"5.59	1.54	1.98	1.98	0.174	1.35	1.32	1.08	0.937	1.35	1.32	1.08	0.937	  2.25	1.86"
HSS2-1/2X2-1/2X1/8	"3.90	1.07	2.15	2.15	0.116	0.998	0.947	0.799	0.965	0.998	0.947	0.799	0.965	  1.61	1.31"
HSS2-1/2X2X1/4		"6.26	1.74	1.80	1.30	0.233	1.33	1.37	1.06	0.874	0.930	1.17	0.930	0.731	  1.90	1.82"
HSS2-1/2X2X3/16		"4.96	1.37	1.98	1.48	0.174	1.12	1.12	0.894	0.904	0.786	0.956	0.786	0.758	  1.55	1.46"
HSS2-1/2X2X1/8		"3.48	0.956	2.15	1.65	0.116	0.833	0.809	0.667	0.934	0.589	0.694	0.589	0.785	  1.12	1.04"
HSS2-1/2X1-1/2X1/4	"5.41	1.51	1.80	0.801	0.233	1.03	1.11	0.822	0.826	0.449	0.764	0.599	0.546	  1.10	1.29"
HSS2-1/2X1-1/2X3/16	"4.32	1.19	1.98	0.978	0.174	0.882	0.915	0.705	0.860	0.390	0.636	0.520	0.572	  0.929	1.05"
HSS2-1/2X1-1/2X1/8	"3.05	0.840	2.15	1.15	0.116	0.668	0.671	0.535	0.892	0.300	0.469	0.399	0.597	  0.687	0.759"
HSS2-1/2X1X3/16		"3.68	1.02	1.98	0.478	0.174	0.646	0.713	0.517	0.796	0.143	0.360	0.285	0.374	  0.412	0.648"
HSS2-1/2X1X1/8		"2.63	0.724	2.15	0.652	0.116	0.503	0.532	0.403	0.834	0.115	0.274	0.230	0.399	  0.322	0.483"
HSS2-1/4X2-1/4X1/4	"6.26	1.74	1.55	1.55	0.233	1.13	1.28	1.01	0.806	1.13	1.28	1.01	0.806	  1.96	1.85"
HSS2-1/4X2-1/4X3/16	"4.96	1.37	1.73	1.73	0.174	0.953	1.04	0.847	0.835	0.953	1.04	0.847	0.835	  1.60	1.48"
HSS2-1/4X2-1/4X1/8	"3.48	0.956	1.90	1.90	0.116	0.712	0.755	0.633	0.863	0.712	0.755	0.633	0.863	  1.15	1.05"
HSS2-1/4X2X3/16		"4.64	1.28	1.73	1.48	0.174	0.859	0.952	0.764	0.819	0.713	0.877	0.713	0.747	  1.32	1.30"
HSS2-1/4X2X1/8		"3.27	0.898	1.90	1.65	0.116	0.646	0.693	0.574	0.848	0.538	0.639	0.538	0.774	  0.957	0.927"
HSS2X2X1/4		"5.41	1.51	1.30	1.30	0.233	0.747	0.964	0.747	0.704	0.747	0.964	0.747	0.704	  1.31	1.41"
HSS2X2X3/16		"4.32	1.19	1.48	1.48	0.174	0.641	0.797	0.641	0.733	0.641	0.797	0.641	0.733	  1.09	1.14"
HSS2X2X1/8		"3.05	0.840	1.65	1.65	0.116	0.486	0.584	0.486	0.761	0.486	0.584	0.486	0.761	  0.796	0.817"
HSS2X1-1/2X3/16		"3.68	1.02	1.48	0.978	0.174	0.495	0.639	0.495	0.697	0.313	0.521	0.417	0.554	  0.664	0.822"
HSS2X1-1/2X1/8		"2.63	0.724	1.65	1.15	0.116	0.383	0.475	0.383	0.728	0.244	0.389	0.325	0.581	  0.496	0.599"
HSS2X1X3/16		"3.04	0.845	1.48	0.478	0.174	0.350	0.480	0.350	0.643	0.112	0.288	0.225	0.365	  0.301	0.505"
HSS2X1X1/8		"2.20	0.608	1.65	0.652	0.116	0.280	0.366	0.280	0.679	0.0922	0.223	0.184	0.390	  0.238	0.380"
}


#AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
array set WSection {
W44X335 "98.5 44.0 15.9 1.03 1.77 31100 1200 1620 1410 17.8 236 150 3.49 74.7"
W44X290 "85.4 43.6 15.8 0.865 1.58 27000 1040 1410 1240 17.8 205 132 3.49 50.9"
W44X262 "77.2 43.3 15.8 0.785 1.42 24100 923 1270 1110 17.7 182 117 3.47 37.3"
W44X230 "67.8 42.9 15.8 0.710 1.22 20800 796 1100 971 17.5 157 101 3.43 24.9"
W40X593 "174 43.0 16.7 1.79 3.23 50400 2520 2760 2340 17.0 481 302 3.80 445"
W40X503 "148 42.1 16.4 1.54 2.76 41600 2040 2320 1980 16.8 394 249 3.72 277"
W40X431 "127 41.3 16.2 1.34 2.36 34800 1690 1960 1690 16.6 328 208 3.65 177"
W40X397 "117 41.0 16.1 1.22 2.20 32000 1540 1800 1560 16.6 300 191 3.64 142"
W40X372 "110 40.6 16.1 1.16 2.05 29600 1420 1680 1460 16.5 277 177 3.60 116"
W40X362 "106 40.6 16.0 1.12 2.01 28900 1380 1640 1420 16.5 270 173 3.60 109"
W40X324 "95.3 40.2 15.9 1.00 1.81 25600 1220 1460 1280 16.4 239 153 3.58 79.4"
W40X297 "87.3 39.8 15.8 0.930 1.65 23200 1090 1330 1170 16.3 215 138 3.54 61.2"
W40X277 "81.5 39.7 15.8 0.830 1.58 21900 1040 1250 1100 16.4 204 132 3.58 51.5"
W40X249 "73.5 39.4 15.8 0.750 1.42 19600 926 1120 993 16.3 182 118 3.55 38.1"
W40X215 "63.5 39.0 15.8 0.650 1.22 16700 803 964 859 16.2 156 101 3.54 24.8"
W40X199 "58.8 38.7 15.8 0.650 1.07 14900 695 869 770 16.0 137 88.2 3.45 18.3"
W40X392 "116 41.6 12.4 1.42 2.52 29900 803 1710 1440 16.1 212 130 2.64 172"
W40X331 "97.7 40.8 12.2 1.22 2.13 24700 644 1430 1210 15.9 172 106 2.57 105"
W40X327 "95.9 40.8 12.1 1.18 2.13 24500 640 1410 1200 16.0 170 105 2.58 103"
W40X294 "86.2 40.4 12.0 1.06 1.93 21900 562 1270 1080 15.9 150 93.5 2.55 76.6"
W40X278 "82.3 40.2 12.0 1.03 1.81 20500 521 1190 1020 15.8 140 87.1 2.52 65.0"
W40X264 "77.4 40.0 11.9 0.960 1.73 19400 493 1130 971 15.8 132 82.6 2.52 56.1"
W40X235 "69.1 39.7 11.9 0.830 1.58 17400 444 1010 875 15.9 118 74.6 2.54 41.3"
W40X211 "62.1 39.4 11.8 0.750 1.42 15500 390 906 786 15.8 105 66.1 2.51 30.4"
W40X183 "53.3 39.0 11.8 0.650 1.20 13200 331 774 675 15.7 88.3 56.0 2.49 19.3"
W40X167 "49.3 38.6 11.8 0.650 1.03 11600 283 693 600 15.3 76.0 47.9 2.40 14.0"
W40X149 "43.8 38.2 11.8 0.630 0.830 9800 229 598 513 15.0 62.2 38.8 2.29 9.36"
W36X652 "192 41.1 17.6 1.97 3.54 50600 3230 2910 2460 16.2 581 367 4.10 593"
W36X529 "156 39.8 17.2 1.61 2.91 39600 2490 2330 1990 16.0 454 289 4.00 327"
W36X487 "143 39.3 17.1 1.50 2.68 36000 2250 2130 1830 15.8 412 263 3.96 258"
W36X441 "130 38.9 17.0 1.36 2.44 32100 1990 1910 1650 15.7 368 235 3.92 194"
W36X395 "116 38.4 16.8 1.22 2.20 28500 1750 1710 1490 15.7 325 208 3.88 142"
W36X361 "106 38.0 16.7 1.12 2.01 25700 1570 1550 1350 15.6 293 188 3.85 109"
W36X330 "96.9 37.7 16.6 1.02 1.85 23300 1420 1410 1240 15.5 265 171 3.83 84.3"
W36X302 "89.0 37.3 16.7 0.945 1.68 21100 1300 1280 1130 15.4 241 156 3.82 64.3"
W36X282 "82.9 37.1 16.6 0.885 1.57 19600 1200 1190 1050 15.4 223 144 3.80 52.7"
W36X262 "77.2 36.9 16.6 0.840 1.44 17900 1090 1100 972 15.3 204 132 3.76 41.6"
W36X247 "72.5 36.7 16.5 0.800 1.35 16700 1010 1030 913 15.2 190 123 3.74 34.7"
W36X231 "68.2 36.5 16.5 0.760 1.26 15600 940 963 854 15.1 176 114 3.71 28.7"
W36X256 "75.3 37.4 12.2 0.960 1.73 16800 528 1040 895 14.9 137 86.5 2.65 52.9"
W36X232 "68.0 37.1 12.1 0.870 1.57 15000 468 936 809 14.8 122 77.2 2.62 39.6"
W36X210 "61.9 36.7 12.2 0.830 1.36 13200 411 833 719 14.6 107 67.5 2.58 28.0"
W36X194 "57.0 36.5 12.1 0.765 1.26 12100 375 767 664 14.6 97.7 61.9 2.56 22.2"
W36X182 "53.6 36.3 12.1 0.725 1.18 11300 347 718 623 14.5 90.7 57.6 2.55 18.5"
W36X170 "50.0 36.2 12.0 0.680 1.10 10500 320 668 581 14.5 83.8 53.2 2.53 15.1"
W36X160 "47.0 36.0 12.0 0.650 1.02 9760 295 624 542 14.4 77.3 49.1 2.50 12.4"
W36X150 "44.3 35.9 12.0 0.625 0.940 9040 270 581 504 14.3 70.9 45.1 2.47 10.1"
W36X135 "39.9 35.6 12.0 0.600 0.790 7800 225 509 439 14.0 59.7 37.7 2.38 7.00"
W33X387 "114 36.0 16.2 1.26 2.28 24300 1620 1560 1350 14.6 312 200 3.77 148"
W33X354 "104 35.6 16.1 1.16 2.09 22000 1460 1420 1240 14.5 282 181 3.74 115"
W33X318 "93.7 35.2 16.0 1.04 1.89 19500 1290 1270 1110 14.5 250 161 3.71 84.4"
W33X291 "85.6 34.8 15.9 0.960 1.73 17700 1160 1160 1020 14.4 226 146 3.68 65.1"
W33X263 "77.4 34.5 15.8 0.870 1.57 15900 1040 1040 919 14.3 202 131 3.66 48.7"
W33X241 "71.1 34.2 15.9 0.830 1.40 14200 933 940 831 14.1 182 118 3.62 36.2"
W33X221 "65.3 33.9 15.8 0.775 1.28 12900 840 857 759 14.1 164 106 3.59 27.8"
W33X201 "59.1 33.7 15.7 0.715 1.15 11600 749 773 686 14.0 147 95.2 3.56 20.8"
W33X169 "49.5 33.8 11.5 0.670 1.22 9290 310 629 549 13.7 84.4 53.9 2.50 17.7"
W33X152 "44.9 33.5 11.6 0.635 1.06 8160 273 559 487 13.5 73.9 47.2 2.47 12.4"
W33X141 "41.5 33.3 11.5 0.605 0.960 7450 246 514 448 13.4 66.9 42.7 2.43 9.70"
W33X130 "38.3 33.1 11.5 0.580 0.855 6710 218 467 406 13.2 59.5 37.9 2.39 7.37"
W33X118 "34.7 32.9 11.5 0.550 0.740 5900 187 415 359 13.0 51.3 32.6 2.32 5.30"
W30X391 "115 33.2 15.6 1.36 2.44 20700 1550 1450 1250 13.4 310 198 3.67 173"
W30X357 "105 32.8 15.5 1.24 2.24 18700 1390 1320 1140 13.3 279 179 3.64 134"
W30X326 "95.9 32.4 15.4 1.14 2.05 16800 1240 1190 1040 13.2 252 162 3.60 103"
W30X292 "86.0 32.0 15.3 1.02 1.85 14900 1100 1060 930 13.2 223 144 3.58 75.2"
W30X261 "77.0 31.6 15.2 0.930 1.65 13100 959 943 829 13.1 196 127 3.53 54.1"
W30X235 "69.3 31.3 15.1 0.830 1.50 11700 855 847 748 13.0 175 114 3.51 40.3"
W30X211 "62.3 30.9 15.1 0.775 1.32 10300 757 751 665 12.9 155 100 3.49 28.4"
W30X191 "56.1 30.7 15.0 0.710 1.19 9200 673 675 600 12.8 138 89.5 3.46 21.0"
W30X173 "50.9 30.4 15.0 0.655 1.07 8230 598 607 541 12.7 123 79.8 3.42 15.6"
W30X148 "43.6 30.7 10.5 0.650 1.18 6680 227 500 436 12.4 68.0 43.3 2.28 14.5"
W30X132 "38.8 30.3 10.5 0.615 1.00 5770 196 437 380 12.2 58.4 37.2 2.25 9.72"
W30X124 "36.5 30.2 10.5 0.585 0.930 5360 181 408 355 12.1 54.0 34.4 2.23 7.99"
W30X116 "34.2 30.0 10.5 0.565 0.850 4930 164 378 329 12.0 49.2 31.3 2.19 6.43"
W30X108 "31.7 29.8 10.5 0.545 0.760 4470 146 346 299 11.9 43.9 27.9 2.15 4.99"
W30X99  "29.0 29.7 10.5 0.520 0.670 3990 128 312 269 11.7 38.6 24.5 2.10 3.77"
W30X90  "26.3 29.5 10.4 0.470 0.610 3610 115 283 245 11.7 34.7 22.1 2.09 2.84"
W27X539 "159 32.5 15.3 1.97 3.54 25600 2110 1890 1570 12.7 437 277 3.65 496"
W27X368 "109 30.4 14.7 1.38 2.48 16200 1310 1240 1060 12.2 279 179 3.48 170"
W27X336 "99.2 30.0 14.6 1.26 2.28 14600 1180 1130 972 12.1 252 162 3.45 131"
W27X307 "90.2 29.6 14.4 1.16 2.09 13100 1050 1030 887 12.0 227 146 3.41 101"
W27X281 "83.1 29.3 14.4 1.06 1.93 11900 953 936 814 12.0 206 133 3.39 79.5"
W27X258 "76.1 29.0 14.3 0.980 1.77 10800 859 852 745 11.9 187 120 3.36 61.6"
W27X235 "69.4 28.7 14.2 0.910 1.61 9700 769 772 677 11.8 168 108 3.33 47.0"
W27X217 "63.9 28.4 14.1 0.830 1.50 8910 704 711 627 11.8 154 100 3.32 37.6"
W27X194 "57.1 28.1 14.0 0.750 1.34 7860 619 631 559 11.7 136 88.1 3.29 27.1"
W27X178 "52.5 27.8 14.1 0.725 1.19 7020 555 570 505 11.6 122 78.8 3.25 20.1"
W27X161 "47.6 27.6 14.0 0.660 1.08 6310 497 515 458 11.5 109 70.9 3.23 15.1"
W27X146 "43.2 27.4 14.0 0.605 0.975 5660 443 464 414 11.5 97.7 63.5 3.20 11.3"
W27X129 "37.8 27.6 10.0 0.610 1.10 4760 184 395 345 11.2 57.6 36.8 2.21 11.1"
W27X114 "33.6 27.3 10.1 0.570 0.930 4080 159 343 299 11.0 49.3 31.5 2.18 7.33"
W27X102 "30.0 27.1 10.0 0.515 0.830 3620 139 305 267 11.0 43.4 27.8 2.15 5.28"
W27X94  "27.6 26.9 10.0 0.490 0.745 3270 124 278 243 10.9 38.8 24.8 2.12 4.03"
W27X84  "24.7 26.7 10.0 0.460 0.640 2850 106 244 213 10.7 33.2 21.2 2.07 2.81"
W24X370 "109 28.0 13.7 1.52 2.72 13400 1160 1130 957 11.1 267 170 3.27 201"
W24X335 "98.3 27.5 13.5 1.38 2.48 11900 1030 1020 864 11.0 238 152 3.23 152"
W24X306 "89.7 27.1 13.4 1.26 2.28 10700 919 922 789 10.9 214 137 3.20 117"
W24X279 "81.9 26.7 13.3 1.16 2.09 9600 823 835 718 10.8 193 124 3.17 90.5"
W24X250 "73.5 26.3 13.2 1.04 1.89 8490 724 744 644 10.7 171 110 3.14 66.6"
W24X229 "67.2 26.0 13.1 0.960 1.73 7650 651 675 588 10.7 154 99.4 3.11 51.3"
W24X207 "60.7 25.7 13.0 0.870 1.57 6820 578 606 531 10.6 137 88.8 3.08 38.3"
W24X192 "56.5 25.5 13.0 0.810 1.46 6260 530 559 491 10.5 126 81.8 3.07 30.8"
W24X176 "51.7 25.2 12.9 0.750 1.34 5680 479 511 450 10.5 115 74.3 3.04 23.9"
W24X162 "47.8 25.0 13.0 0.705 1.22 5170 443 468 414 10.4 105 68.4 3.05 18.5"
W24X146 "43.0 24.7 12.9 0.650 1.09 4580 391 418 371 10.3 93.2 60.5 3.01 13.4"
W24X131 "38.6 24.5 12.9 0.605 0.960 4020 340 370 329 10.2 81.5 53.0 2.97 9.50"
W24X117 "34.4 24.3 12.8 0.550 0.850 3540 297 327 291 10.1 71.4 46.5 2.94 6.72"
W24X104 "30.7 24.1 12.8 0.500 0.750 3100 259 289 258 10.1 62.4 40.7 2.91 4.72"
W24X103 "30.3 24.5 9.00 0.550 0.980 3000 119 280 245 10.0 41.5 26.5 1.99 7.07"
W24X94 "27.7 24.3 9.07 0.515 0.875 2700 109 254 222 9.87 37.5 24.0 1.98 5.26"
W24X84 "24.7 24.1 9.02 0.470 0.770 2370 94.4 224 196 9.79 32.6 20.9 1.95 3.70"
W24X76 "22.4 23.9 8.99 0.440 0.680 2100 82.5 200 176 9.69 28.6 18.4 1.92 2.68"
W24X68 "20.1 23.7 8.97 0.415 0.585 1830 70.4 177 154 9.55 24.5 15.7 1.87 1.87"
W24X62 "18.2 23.7 7.04 0.430 0.590 1550 34.5 153 131 9.23 15.7 9.80 1.38 1.71"
W24X55 "16.2 23.6 7.01 0.395 0.505 1350 29.1 134 114 9.11 13.3 8.30 1.34 1.18"
W21X201 "59.3 23.0 12.6 0.910 1.63 5310 542 530 461 9.47 133 86.1 3.02 40.9"
W21X182 "53.6 22.7 12.5 0.830 1.48 4730 483 476 417 9.40 119 77.2 3.00 30.7"
W21X166 "48.8 22.5 12.4 0.750 1.36 4280 435 432 380 9.36 108 70.0 2.99 23.6"
W21X147 "43.2 22.1 12.5 0.720 1.15 3630 376 373 329 9.17 92.6 60.1 2.95 15.4"
W21X132 "38.8 21.8 12.4 0.650 1.04 3220 333 333 295 9.12 82.3 53.5 2.93 11.3"
W21X122 "35.9 21.7 12.4 0.600 0.960 2960 305 307 273 9.09 75.6 49.2 2.92 8.98"
W21X111 "32.6 21.5 12.3 0.550 0.875 2670 274 279 249 9.05 68.2 44.5 2.90 6.83"
W21X101 "29.8 21.4 12.3 0.500 0.800 2420 248 253 227 9.02 61.7 40.3 2.89 5.21"
W21X93 "27.3 21.6 8.42 0.580 0.930 2070 92.9 221 192 8.70 34.7 22.1 1.84 6.03"
W21X83 "24.4 21.4 8.36 0.515 0.835 1830 81.4 196 171 8.67 30.5 19.5 1.83 4.34"
W21X73 "21.5 21.2 8.30 0.455 0.740 1600 70.6 172 151 8.64 26.6 17.0 1.81 3.02"
W21X68 "20.0 21.1 8.27 0.430 0.685 1480 64.7 160 140 8.60 24.4 15.7 1.80 2.45"
W21X62 "18.3 21.0 8.24 0.400 0.615 1330 57.5 144 127 8.54 21.7 14.0 1.77 1.83"
W21X55 "16.2 20.8 8.22 0.375 0.522 1140 48.4 126 110 8.40 18.4 11.8 1.73 1.24"
W21X48 "14.1 20.6 8.14 0.350 0.430 959 38.7 107 93.0 8.24 14.9 9.52 1.66 0.803"
W21X57 "16.7 21.1 6.56 0.405 0.650 1170 30.6 129 111 8.36 14.8 9.35 1.35 1.77"
W21X50 "14.7 20.8 6.53 0.380 0.535 984 24.9 110 94.5 8.18 12.2 7.64 1.30 1.14"
W21X44 "13.0 20.7 6.50 0.350 0.450 843 20.7 95.4 81.6 8.06 10.2 6.37 1.26 0.770"
W18X311 "91.6 22.3 12.0 1.52 2.74 6970 795 754 624 8.72 207 132 2.95 176"
W18X283 "83.3 21.9 11.9 1.40 2.50 6170 704 676 565 8.61 185 118 2.91 134"
W18X258 "76.0 21.5 11.8 1.28 2.30 5510 628 611 514 8.53 166 107 2.88 103"
W18X234 "68.6 21.1 11.7 1.16 2.11 4900 558 549 466 8.44 149 95.8 2.85 78.7"
W18X211 "62.3 20.7 11.6 1.06 1.91 4330 493 490 419 8.35 132 85.3 2.82 58.6"
W18X192 "56.2 20.4 11.5 0.960 1.75 3870 440 442 380 8.28 119 76.8 2.79 44.7"
W18X175 "51.4 20.0 11.4 0.890 1.59 3450 391 398 344 8.20 106 68.8 2.76 33.8"
W18X158 "46.3 19.7 11.3 0.810 1.44 3060 347 356 310 8.12 94.8 61.4 2.74 25.2"
W18X143 "42.0 19.5 11.2 0.730 1.32 2750 311 322 282 8.09 85.4 55.5 2.72 19.2"
W18X130 "38.3 19.3 11.2 0.670 1.20 2460 278 290 256 8.03 76.7 49.9 2.70 14.5"
W18X119 "35.1 19.0 11.3 0.655 1.06 2190 253 262 231 7.90 69.1 44.9 2.69 10.6"
W18X106 "31.1 18.7 11.2 0.590 0.940 1910 220 230 204 7.84 60.5 39.4 2.66 7.48"
W18X97 "28.5 18.6 11.1 0.535 0.870 1750 201 211 188 7.82 55.3 36.1 2.65 5.86"
W18X86 "25.3 18.4 11.1 0.480 0.770 1530 175 186 166 7.77 48.4 31.6 2.63 4.10"
W18X76 "22.3 18.2 11.0 0.425 0.680 1330 152 163 146 7.73 42.2 27.6 2.61 2.83"
W18X71 "20.9 18.5 7.64 0.495 0.810 1170 60.3 146 127 7.50 24.7 15.8 1.70 3.49"
W18X65 "19.1 18.4 7.59 0.450 0.750 1070 54.8 133 117 7.49 22.5 14.4 1.69 2.73"
W18X60 "17.6 18.2 7.56 0.415 0.695 984 50.1 123 108 7.47 20.6 13.3 1.68 2.17"
W18X55 "16.2 18.1 7.53 0.390 0.630 890 44.9 112 98.3 7.41 18.5 11.9 1.67 1.66"
W18X50 "14.7 18.0 7.50 0.355 0.570 800 40.1 101 88.9 7.38 16.6 10.7 1.65 1.24"
W18X46 "13.5 18.1 6.06 0.360 0.605 712 22.5 90.7 78.8 7.25 11.7 7.43 1.29 1.22"
W18X40 "11.8 17.9 6.02 0.315 0.525 612 19.1 78.4 68.4 7.21 10.0 6.35 1.27 0.810"
W18X35 "10.3 17.7 6.00 0.300 0.425 510 15.3 66.5 57.6 7.04 8.06 5.12 1.22 0.506"
W16X100 "29.4 17.0 10.4 0.585 0.985 1490 186 198 175 7.10 54.9 35.7 2.51 7.73"
W16X89 "26.2 16.8 10.4 0.525 0.875 1300 163 175 155 7.05 48.1 31.4 2.49 5.45"
W16X77 "22.6 16.5 10.3 0.455 0.760 1110 138 150 134 7.00 41.1 26.9 2.47 3.57"
W16X67 "19.6 16.3 10.2 0.395 0.665 954 119 130 117 6.96 35.5 23.2 2.46 2.39"
W16X57 "16.8 16.4 7.12 0.430 0.715 758 43.1 105 92.2 6.72 18.9 12.1 1.60 2.22"
W16X50 "14.7 16.3 7.07 0.380 0.630 659 37.2 92.0 81.0 6.68 16.3 10.5 1.59 1.52"
W16X45 "13.3 16.1 7.04 0.345 0.565 586 32.8 82.3 72.7 6.65 14.5 9.34 1.57 1.11"
W16X40 "11.8 16.0 7.00 0.305 0.505 518 28.9 73.0 64.7 6.63 12.7 8.25 1.57 0.794"
W16X36 "10.6 15.9 6.99 0.295 0.430 448 24.5 64.0 56.5 6.51 10.8 7.00 1.52 0.545"
W16X31 "9.13 15.9 5.53 0.275 0.440 375 12.4 54.0 47.2 6.41 7.03 4.49 1.17 0.461"
W16X26 "7.68 15.7 5.50 0.250 0.345 301 9.59 44.2 38.4 6.26 5.48 3.49 1.12 0.262"
W14X730 "215 22.4 17.9 3.07 4.91 14300 4720 1660 1280 8.17 816 527 4.69 1450"
W14X665 "196 21.6 17.7 2.83 4.52 12400 4170 1480 1150 7.98 730 472 4.62 1120"
W14X605 "178 20.9 17.4 2.60 4.16 10800 3680 1320 1040 7.80 652 423 4.55 869"
W14X550 "162 20.2 17.2 2.38 3.82 9430 3250 1180 931 7.63 583 378 4.49 669"
W14X500 "147 19.6 17.0 2.19 3.50 8210 2880 1050 838 7.48 522 339 4.43 514"
W14X455 "134 19.0 16.8 2.02 3.21 7190 2560 936 756 7.33 468 304 4.38 395"
W14X426 "125 18.7 16.7 1.88 3.04 6600 2360 869 706 7.26 434 283 4.34 331"
W14X398 "117 18.3 16.6 1.77 2.85 6000 2170 801 656 7.16 402 262 4.31 273"
W14X370 "109 17.9 16.5 1.66 2.66 5440 1990 736 607 7.07 370 241 4.27 222"
W14X342 "101 17.5 16.4 1.54 2.47 4900 1810 672 558 6.98 338 221 4.24 178"
W14X311 "91.4 17.1 16.2 1.41 2.26 4330 1610 603 506 6.88 304 199 4.20 136"
W14X283 "83.3 16.7 16.1 1.29 2.07 3840 1440 542 459 6.79 274 179 4.17 104"
W14X257 "75.6 16.4 16.0 1.18 1.89 3400 1290 487 415 6.71 246 161 4.13 79.1"
W14X233 "68.5 16.0 15.9 1.07 1.72 3010 1150 436 375 6.63 221 145 4.10 59.5"
W14X211 "62.0 15.7 15.8 0.980 1.56 2660 1030 390 338 6.55 198 130 4.07 44.6"
W14X193 "56.8 15.5 15.7 0.890 1.44 2400 931 355 310 6.50 180 119 4.05 34.8"
W14X176 "51.8 15.2 15.7 0.830 1.31 2140 838 320 281 6.43 163 107 4.02 26.5"
W14X159 "46.7 15.0 15.6 0.745 1.19 1900 748 287 254 6.38 146 96.2 4.00 19.7"
W14X145 "42.7 14.8 15.5 0.680 1.09 1710 677 260 232 6.33 133 87.3 3.98 15.2"
W14X132 "38.8 14.7 14.7 0.645 1.03 1530 548 234 209 6.28 113 74.5 3.76 12.3"
W14X120 "35.3 14.5 14.7 0.590 0.940 1380 495 212 190 6.24 102 67.5 3.74 9.37"
W14X109 "32.0 14.3 14.6 0.525 0.860 1240 447 192 173 6.22 92.7 61.2 3.73 7.12"
W14X99 "29.1 14.2 14.6 0.485 0.780 1110 402 173 157 6.17 83.6 55.2 3.71 5.37"
W14X90 "26.5 14.0 14.5 0.440 0.710 999 362 157 143 6.14 75.6 49.9 3.70 4.06"
W14X82 "24.0 14.3 10.1 0.510 0.855 881 148 139 123 6.05 44.8 29.3 2.48 5.07"
W14X74 "21.8 14.2 10.1 0.450 0.785 795 134 126 112 6.04 40.5 26.6 2.48 3.87"
W14X68 "20.0 14.0 10.0 0.415 0.720 722 121 115 103 6.01 36.9 24.2 2.46 3.01"
W14X61 "17.9 13.9 10.0 0.375 0.645 640 107 102 92.1 5.98 32.8 21.5 2.45 2.19"
W14X53 "15.6 13.9 8.06 0.370 0.660 541 57.7 87.1 77.8 5.89 22.0 14.3 1.92 1.94"
W14X48 "14.1 13.8 8.03 0.340 0.595 484 51.4 78.4 70.2 5.85 19.6 12.8 1.91 1.45"
W14X43 "12.6 13.7 8.00 0.305 0.530 428 45.2 69.6 62.6 5.82 17.3 11.3 1.89 1.05"
W14X38 "11.2 14.1 6.77 0.310 0.515 385 26.7 61.5 54.6 5.87 12.1 7.88 1.55 0.798"
W14X34 "10.0 14.0 6.75 0.285 0.455 340 23.3 54.6 48.6 5.83 10.6 6.91 1.53 0.569"
W14X30 "8.85 13.8 6.73 0.270 0.385 291 19.6 47.3 42.0 5.73 8.99 5.82 1.49 0.380"
W14X26 "7.69 13.9 5.03 0.255 0.420 245 8.91 40.2 35.3 5.65 5.54 3.55 1.08 0.358"
W14X22 "6.49 13.7 5.00 0.230 0.335 199 7.00 33.2 29.0 5.54 4.39 2.80 1.04 0.208"
W12X336 "98.9 16.8 13.4 1.78 2.96 4060 1190 603 483 6.41 274 177 3.47 243"
W12X305 "89.5 16.3 13.2 1.63 2.71 3550 1050 537 435 6.29 244 159 3.42 185"
W12X279 "81.9 15.9 13.1 1.53 2.47 3110 937 481 393 6.16 220 143 3.38 143"
W12X252 "74.1 15.4 13.0 1.40 2.25 2720 828 428 353 6.06 196 127 3.34 108"
W12X230 "67.7 15.1 12.9 1.29 2.07 2420 742 386 321 5.97 177 115 3.31 83.8"
W12X210 "61.8 14.7 12.8 1.18 1.90 2140 664 348 292 5.89 159 104 3.28 64.7"
W12X190 "56.0 14.4 12.7 1.06 1.74 1890 589 311 263 5.82 143 93.0 3.25 48.8"
W12X170 "50.0 14.0 12.6 0.960 1.56 1650 517 275 235 5.74 126 82.3 3.22 35.6"
W12X152 "44.7 13.7 12.5 0.870 1.40 1430 454 243 209 5.66 111 72.8 3.19 25.8"
W12X136 "39.9 13.4 12.4 0.790 1.25 1240 398 214 186 5.58 98.0 64.2 3.16 18.5"
W12X120 "35.2 13.1 12.3 0.710 1.11 1070 345 186 163 5.51 85.4 56.0 3.13 12.9"
W12X106 "31.2 12.9 12.2 0.610 0.990 933 301 164 145 5.47 75.1 49.3 3.11 9.13"
W12X96 "28.2 12.7 12.2 0.550 0.900 833 270 147 131 5.44 67.5 44.4 3.09 6.85"
W12X87 "25.6 12.5 12.1 0.515 0.810 740 241 132 118 5.38 60.4 39.7 3.07 5.10"
W12X79 "23.2 12.4 12.1 0.470 0.735 662 216 119 107 5.34 54.3 35.8 3.05 3.84"
W12X72 "21.1 12.3 12.0 0.430 0.670 597 195 108 97.4 5.31 49.2 32.4 3.04 2.93"
W12X65 "19.1 12.1 12.0 0.390 0.605 533 174 96.8 87.9 5.28 44.1 29.1 3.02 2.18"
W12X58 "17.0 12.2 10.0 0.360 0.640 475 107 86.4 78.0 5.28 32.5 21.4 2.51 2.10"
W12X53 "15.6 12.1 10.0 0.345 0.575 425 95.8 77.9 70.6 5.23 29.1 19.2 2.48 1.58"
W12X50 "14.6 12.2 8.08 0.370 0.640 391 56.3 71.9 64.2 5.18 21.3 13.9 1.96 1.71"
W12X45 "13.1 12.1 8.05 0.335 0.575 348 50.0 64.2 57.7 5.15 19.0 12.4 1.95 1.26"
W12X40 "11.7 11.9 8.01 0.295 0.515 307 44.1 57.0 51.5 5.13 16.8 11.0 1.94 0.906"
W12X35 "10.3 12.5 6.56 0.300 0.520 285 24.5 51.2 45.6 5.25 11.5 7.47 1.54 0.741"
W12X30 "8.79 12.3 6.52 0.260 0.440 238 20.3 43.1 38.6 5.21 9.56 6.24 1.52 0.457"
W12X26 "7.65 12.2 6.49 0.230 0.380 204 17.3 37.2 33.4 5.17 8.17 5.34 1.51 0.300"
W12X22 "6.48 12.3 4.03 0.260 0.425 156 4.66 29.3 25.4 4.91 3.66 2.31 0.848 0.293"
W12X19 "5.57 12.2 4.01 0.235 0.350 130 3.76 24.7 21.3 4.82 2.98 1.88 0.822 0.180"
W12X16 "4.71 12.0 3.99 0.220 0.265 103 2.82 20.1 17.1 4.67 2.26 1.41 0.773 0.103"
W12X14 "4.16 11.9 3.97 0.200 0.225 88.6 2.36 17.4 14.9 4.62 1.90 1.19 0.753 0.0704"
W10X112 "32.9 11.4 10.4 0.755 1.25 716 236 147 126 4.66 69.2 45.3 2.68 15.1"
W10X100 "29.3 11.1 10.3 0.680 1.12 623 207 130 112 4.60 61.0 40.0 2.65 10.9"
W10X88 "26.0 10.8 10.3 0.605 0.990 534 179 113 98.5 4.54 53.1 34.8 2.63 7.53"
W10X77 "22.7 10.6 10.2 0.530 0.870 455 154 97.6 85.9 4.49 45.9 30.1 2.60 5.11"
W10X68 "19.9 10.4 10.1 0.470 0.770 394 134 85.3 75.7 4.44 40.1 26.4 2.59 3.56"
W10X60 "17.7 10.2 10.1 0.420 0.680 341 116 74.6 66.7 4.39 35.0 23.0 2.57 2.48"
W10X54 "15.8 10.1 10.0 0.370 0.615 303 103 66.6 60.0 4.37 31.3 20.6 2.56 1.82"
W10X49 "14.4 10.0 10.0 0.340 0.560 272 93.4 60.4 54.6 4.35 28.3 18.7 2.54 1.39"
W10X45 "13.3 10.1 8.02 0.350 0.620 248 53.4 54.9 49.1 4.32 20.3 13.3 2.01 1.51"
W10X39 "11.5 9.92 7.99 0.315 0.530 209 45.0 46.8 42.1 4.27 17.2 11.3 1.98 0.976"
W10X33 "9.71 9.73 7.96 0.290 0.435 171 36.6 38.8 35.0 4.19 14.0 9.20 1.94 0.583"
W10X30 "8.84 10.5 5.81 0.300 0.510 170 16.7 36.6 32.4 4.38 8.84 5.75 1.37 0.622"
W10X26 "7.61 10.3 5.77 0.260 0.440 144 14.1 31.3 27.9 4.35 7.50 4.89 1.36 0.402"
W10X22 "6.49 10.2 5.75 0.240 0.360 118 11.4 26.0 23.2 4.27 6.10 3.97 1.33 0.239"
W10X19 "5.62 10.2 4.02 0.250 0.395 96.3 4.29 21.6 18.8 4.14 3.35 2.14 0.874 0.233"
W10X17 "4.99 10.1 4.01 0.240 0.330 81.9 3.56 18.7 16.2 4.05 2.80 1.78 0.845 0.156"
W10X15 "4.41 9.99 4.00 0.230 0.270 68.9 2.89 16.0 13.8 3.95 2.30 1.45 0.810 0.104"
W10X12 "3.54 9.87 3.96 0.190 0.210 53.8 2.18 12.6 10.9 3.90 1.74 1.10 0.785 0.0547"
W8X67 "19.7 9.00 8.28 0.570 0.935 272 88.6 70.1 60.4 3.72 32.7 21.4 2.12 5.05"
W8X58 "17.1 8.75 8.22 0.510 0.810 228 75.1 59.8 52.0 3.65 27.9 18.3 2.10 3.33"
W8X48 "14.1 8.50 8.11 0.400 0.685 184 60.9 49.0 43.2 3.61 22.9 15.0 2.08 1.96"
W8X40 "11.7 8.25 8.07 0.360 0.560 146 49.1 39.8 35.5 3.53 18.5 12.2 2.04 1.12"
W8X35 "10.3 8.12 8.02 0.310 0.495 127 42.6 34.7 31.2 3.51 16.1 10.6 2.03 0.769"
W8X31 "9.13 8.00 8.00 0.285 0.435 110 37.1 30.4 27.5 3.47 14.1 9.27 2.02 0.536"
W8X28 "8.25 8.06 6.54 0.285 0.465 98.0 21.7 27.2 24.3 3.45 10.1 6.63 1.62 0.537"
W8X24 "7.08 7.93 6.50 0.245 0.400 82.7 18.3 23.1 20.9 3.42 8.57 5.63 1.61 0.346"
W8X21 "6.16 8.28 5.27 0.250 0.400 75.3 9.77 20.4 18.2 3.49 5.69 3.71 1.26 0.282"
W8X18 "5.26 8.14 5.25 0.230 0.330 61.9 7.97 17.0 15.2 3.43 4.66 3.04 1.23 0.172"
W8X15 "4.44 8.11 4.02 0.245 0.315 48.0 3.41 13.6 11.8 3.29 2.67 1.70 0.876 0.137"
W8X13 "3.84 7.99 4.00 0.230 0.255 39.6 2.73 11.4 9.91 3.21 2.15 1.37 0.843 0.0871"
W8X10 "2.96 7.89 3.94 0.170 0.205 30.8 2.09 8.87 7.81 3.22 1.66 1.06 0.841 0.0426"
W6X25 "7.34 6.38 6.08 0.320 0.455 53.4 17.1 18.9 16.7 2.70 8.56 5.61 1.52 0.461"
W6X20 "5.87 6.20 6.02 0.260 0.365 41.4 13.3 15.0 13.4 2.66 6.72 4.41 1.50 0.240"
W6X15 "4.43 5.99 5.99 0.230 0.260 29.1 9.32 10.8 9.72 2.56 4.75 3.11 1.45 0.101"
W6X16 "4.74 6.28 4.03 0.260 0.405 32.1 4.43 11.7 10.2 2.60 3.39 2.20 0.967 0.223"
W6X12 "3.55 6.03 4.00 0.230 0.280 22.1 2.99 8.30 7.31 2.49 2.32 1.50 0.918 0.0903"
W6X9  "2.68 5.90 3.94 0.170 0.215 16.4 2.20 6.23 5.56 2.47 1.72 1.11 0.905 0.0405"
W6X8.5 "2.52 5.83 3.94 0.170 0.195 14.9 1.99 5.73 5.10 2.43 1.56 1.01 0.890 0.0333"
W5X19 "5.56 5.15 5.03 0.270 0.430 26.3 9.13 11.6 10.2 2.17 5.53 3.63 1.28 0.316"
W5X16 "4.71 5.01 5.00 0.240 0.360 21.4 7.51 9.63 8.55 2.13 4.58 3.00 1.26 0.192"
W4X13 "3.83 4.16 4.06 0.280 0.345 11.3 3.86 6.28 5.46 1.72 2.92 1.90 1.00 0.151"
}
