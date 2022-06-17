foreach rigidConstraint {no yes} {

    puts "RIGID CONSTRAINT: $rigidConstraint"

    foreach matType {steel concrete} {
	foreach eleType {forceBeamColumn dispBeamColumn} {

	    wipe

	    # create model builder
	    model basic -ndm 2 -ndf 3
	    
	    set width    360
	    set height   144
	    
	    # create nodes
	    
	    node  1       0.0     0.0 
	    node  2    $width     0.0 
	    node  3       0.0 $height
	    node  4    $width $height
	    
	    # Fix supports at base of columns
	    #    tag   DX   DY   RZ
	    
	    fix   1     1    1    1
	    fix   2     1    1    1
	    
	    if {$rigidConstraint == "yes"} {
		equalDOF 3 4 1
	    }
	    
	    if {$matType == "concrete"} {
		# CONCRETE                  tag   f'c        ec0   f'cu        ecu
		# Core concrete (confined)
		uniaxialMaterial Concrete01  1  -6.0  -0.004   -5.0     -0.014
		
		# Cover concrete (unconfined)
		uniaxialMaterial Concrete01  2  -5.0   -0.002   0.0     -0.006
		
		# STEEL
		# Reinforcing steel 
		set fy 60.0;      # Yield stress
		set E 30000.0;    # Young's modulus
		#                        tag  fy E0    b
		uniaxialMaterial Steel01  3  $fy $E 0.01
		
		
	    } else {
		
		set fy 60.0;      # Yield stress
		set E 30000.0;    # Young's modulus
		#                        tag  fy E0    b
		uniaxialMaterial Steel01  1  $fy $E 0.01
		uniaxialMaterial Steel01  2  $fy $E 0.01
		uniaxialMaterial Steel01  3  $fy $E 0.01
	    }
	    
	    # Define cross-section for nonlinear columns
	    # ------------------------------------------
	    
	    # set some parameters
	    set colWidth 15
	    set colDepth 24 
	    
	    set cover  1.5
	    set As    0.60;     # area of no. 7 bars
	    
	    # some variables derived from the parameters
	    set y1 [expr $colDepth/2.0]
	    set z1 [expr $colWidth/2.0]
	    
	    section Fiber 1 {
		
		# Create the concrete core fibers
		patch rect 1 10 1 [expr $cover-$y1] [expr $cover-$z1] [expr $y1-$cover] [expr $z1-$cover]
		
		# Create the concrete cover fibers (top, bottom, left, right)
		patch rect 2 10 1  [expr -$y1] [expr $z1-$cover] $y1 $z1
		patch rect 2 10 1  [expr -$y1] [expr -$z1] $y1 [expr $cover-$z1]
		patch rect 2  2 1  [expr -$y1] [expr $cover-$z1] [expr $cover-$y1] [expr $z1-$cover]
		patch rect 2  2 1  [expr $y1-$cover] [expr $cover-$z1] $y1 [expr $z1-$cover]
		
		# Create the reinforcing fibers (left, middle, right)
		layer straight 3 3 $As [expr $y1-$cover] [expr $z1-$cover] [expr $y1-$cover] [expr $cover-$z1]
		layer straight 3 2 $As 0.0 [expr $z1-$cover] 0.0 [expr $cover-$z1]
		layer straight 3 3 $As [expr $cover-$y1] [expr $z1-$cover] [expr $cover-$y1] [expr $cover-$z1]	
	    }    
	    
	    
	    # Define elements
	    # ----------------------
	    
	    geomTransf Linear 1  
	    set np 5
	    

	    # Create the columns using Beam-column elements
	    #               e            tag ndI ndJ nsecs secID transfTag
	    element $eleType  1   1   3   $np    1       1 
	    element $eleType  2   2   4   $np    1       1 
	    element $eleType  3   3   4   $np    1       1
	    
	    # Define gravity loads
	    # --------------------
	    
	    # Set a parameter for the axial load
	    set P 180;                # 10% of axial capacity of columns
	    
	    # Create a Plain load pattern with a Linear TimeSeries
	    pattern Plain 1 "Linear" {
		# Create nodal loads at nodes 3 & 4
		#    nd    FX          FY  MZ 
		load  3   0.0  [expr -$P] 0.0
		load  4   0.0  [expr -$P] 0.0
	    }
	    
	    
	    system BandGeneral
	    constraints Transformation
	    numberer RCM

	    test NormDispIncr 1.0e-12  10 3
	    algorithm Newton
	    integrator LoadControl 0.1
	    analysis Static
	    
	    analyze 10
	    
	    loadConst -time 0.0
	    
	    pattern Plain 2 "Linear" {
		load  3   1.0  0.0 0.0
		load  4   1.0  0.0 0.0
	    }

	    integrator DisplacementControl 3 1 0.1
	    analyze 10
	    set strains [eleResponse 3 basicDeformation]
	    set strains1 [eleResponse 3 section 1 forces]
	    set forces [eleResponse 3 forces]
	    

	    puts "$strains1 $forces"
	    puts "eleType: $eleType matType $matType axialForce [lindex $forces 0] axialDeformation: [lindex $strains 0]"
	}
    }
}