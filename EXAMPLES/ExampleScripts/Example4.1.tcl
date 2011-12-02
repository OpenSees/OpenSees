# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# 2 Story Multi Bay Frame Example 1
# -----------------------------
#  Reinforced concrete multi-bay, two-story frame
#  Distributed vertical load on girder
# 
# Example Objectives
# -----------------
#  Nonlinear beam-column elements
#  Gravity load analysis followed by pushover analysis
#  Demonstrate scripting for the algorithmic level
#
# 
# Units: kips, in, sec
#
# Written: GLF/MHS/fmk
# Date: January 2001


# Parameter identifying the number of bays
set numBay          3

# ------------------------------
# Start of model generation
# ------------------------------

# Create ModelBuilder (with two-dimensions and 3 DOF/node)
model basic -ndm 2 -ndf 3

# Create nodes
# ------------

# Set parameters for overall model geometry
set bayWidth      288
set nodeID          1

# Define nodes
for {set i 0} {$i <= $numBay} {incr i 1} {
    
    set xDim [expr $i * $bayWidth]

    #             tag             X   Y
    node           $nodeID    $xDim  0
    node  [expr $nodeID+1]    $xDim 180
    node  [expr $nodeID+2]    $xDim 324

    incr nodeID 3
}

# Fix supports at base of columns
for {set i 0} {$i <= $numBay} {incr i 1} {
#                node  DX   DY   RZ
    fix [expr $i*3+1]   1    1    1
}

# Define materials for nonlinear columns
# ------------------------------------------

# CONCRETE
# Cover concrete
#                  tag -f'c  -epsco  -f'cu -epscu
# Core concrete (confined)
uniaxialMaterial Concrete01  1  -6.0  -0.004   -5.0     -0.014


# Cover concrete (unconfined)
uniaxialMaterial Concrete01  2  -5.0   -0.002   0.0     -0.006

# STEEL
# Reinforcing steel 
#                        tag fy   E0     b
uniaxialMaterial Steel01  3  60 30000 0.015


# Define cross-section for nonlinear columns
# ------------------------------------------

# Interior column section
section Fiber 1 {
   #           mat nfIJ nfJK   yI  zI    yK  zK  
   #           mat nfIJ nfJK   yI  zI    yJ  zJ    yK  zK    yL  zL
   patch quadr  2    1   12 -11.5  10 -11.5 -10  11.5 -10  11.5  10
   patch quadr  1    1   14 -13.5 -10 -13.5 -12  13.5 -12  13.5 -10
   patch quadr  1    1   14 -13.5  12 -13.5  10  13.5  10  13.5  12
   patch quadr  1    1    2 -13.5  10 -13.5 -10 -11.5 -10 -11.5  10
   patch quadr  1    1    2  11.5  10  11.5 -10  13.5 -10  13.5  10

   #              mat nBars area    yI zI    yF zF
   layer straight  3    6   1.56 -10.5  9 -10.5 -9
   layer straight  3    6   1.56  10.5  9  10.5 -9
}

# Exterior column section
section Fiber 2 {
   patch quadr 2 1 10 -10  10 -10 -10  10 -10  10  10
   patch quadr 1 1 12 -12 -10 -12 -12  12 -12  12 -10
   patch quadr 1 1 12 -12  12 -12  10  12  10  12  12
   patch quadr 1 1  2 -12  10 -12 -10 -10 -10 -10  10
   patch quadr 1 1  2  10  10  10 -10  12 -10  12  10
   layer straight 3 6 0.79 -9 9 -9 -9
   layer straight 3 6 0.79  9 9  9 -9
}

# Girder section
section Fiber 3 {
   patch quadr 1 1 12 -12 9 -12 -9 12 -9 12 9
   layer straight 3 4 1.00 -9 9 -9 -9
   layer straight 3 4 1.00  9 9  9 -9
}


# Define column elements
# ----------------------

# Number of integration points
set nP 4

# Geometric transformation
geomTransf Linear 1

set beamID          1

# Define elements
for {set i 0} {$i <= $numBay} {incr i 1} {

    # set some parameters
    set iNode [expr $i*3 + 1]
    set jNode [expr $i*3 + 2]

    for {set j 1} {$j < 3} {incr j 1} {    
	# add the column element (secId == 2 if external, 1 if internal column)
	if {$i == 0} {
	    element nonlinearBeamColumn  $beamID   $iNode $jNode    $nP   2      1
	} elseif {$i == $numBay} {
	    element nonlinearBeamColumn  $beamID   $iNode $jNode    $nP   2      1
	} else {
	    element nonlinearBeamColumn  $beamID   $iNode $jNode    $nP   1      1
	}

	# increment the parameters
	incr iNode  1
	incr jNode  1
	incr beamID 1
    }

}


# Define beam elements
# ----------------------

# Number of integration points
set nP 4

# Geometric transformation
geomTransf Linear 2

# Define elements
for {set j 1} {$j < 3} {incr j 1} {    
    # set some parameters
    set iNode [expr $j + 1]
    set jNode [expr $iNode + 3]


    for {set i 1} {$i <= $numBay} {incr i 1} {
	element nonlinearBeamColumn  $beamID   $iNode $jNode    $nP   3      2

	# increment the parameters
	incr iNode  3
	incr jNode  3
	incr beamID 1
    }
}

# Define gravity loads
# --------------------

# Constant gravity load
set P -192

# Create a Plain load pattern with a Linear TimeSeries
pattern Plain 1 Linear {

    # Create nodal loads at nodes 
    for {set i 0} {$i <= $numBay} {incr i 1} {

	# set some parameters
	set node1 [expr $i*3 + 2]
	set node2 [expr $node1 + 1]

	if {$i == 0} {
	    load $node1 0.0            $P  0.0 
	    load $node2 0.0 [expr $P/2.0]  0.0 
	} elseif {$i == $numBay} {
	    load $node1 0.0            $P  0.0 
	    load $node2 0.0 [expr $P/2.0]  0.0 
	} else {
	    load $node1 0.0 [expr 2.0*$P]  0.0 
	    load $node2 0.0            $P  0.0 
	}
    }
}

# ------------------------------
# End of model generation
# ------------------------------


# ------------------------------------------------
# Start of analysis generation for gravity analysis
# -------------------------------------------------

# Create the convergence test, the norm of the residual with a tolerance of 
# 1e-12 and a max number of iterations of 10
test NormDispIncr 1.0e-8  10 0

# Create the solution algorithm, a Newton-Raphson algorithm
algorithm Newton

# Create the integration scheme, the LoadControl scheme using steps of 0.1 
integrator LoadControl 0.1

# Create the system of equation, a SPD using a profile storage scheme
system BandGeneral

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer RCM

# Create the constraint handler, the transformation method
constraints Plain

# Create the analysis object
analysis Static


# ------------------------------------------------
# End of analysis generation for gravity analysis
# -------------------------------------------------


# ------------------------------
# Perform gravity load analysis
# ------------------------------

# initialize the model, done to set initial tangent
initialize

# perform the gravity load analysis, requires 10 steps to reach the load level
analyze 10


# set gravity loads to be const and set pseudo time to be 0.0
#  for start of lateral load analysis
loadConst -time 0.0


# ------------------------------
# Add lateral loads 
# ------------------------------

# Reference lateral load for pushover analysis
set H   10

# Reference lateral loads
# Create a Plain load pattern with a Linear TimeSeries
pattern Plain 2 Linear {

    load 2 [expr $H/2.0]  0.0  0.0
    load 3            $H  0.0  0.0
}

# ------------------------------
# Start of recorder generation
# ------------------------------

# Create a recorder which writes to Node.out and prints
# the current load factor (pseudo-time) and dof 1 displacements at node 2 & 3
recorder Node -file Node41.out -time -node 2 3 -dof 1 disp

# Source in some commands to display the model
# comment out one of lines
set displayMode "displayON"
#set displayMode "displayOFF"

if {$displayMode == "displayON"} {
    # a window to plot the nodal displacements versus load for node 3
    recorder plot Node41.out Node_3_Xdisp 10 340 300 300 -columns 3 1 -dT 0.1
}

# ------------------------------
# End of recorder generation
# ------------------------------


# ------------------------------
# Start of lateral load analysis
# ------------------------------

# Change the integrator to take a min and max load increment
integrator LoadControl 1.0 4 0.02 2.0

# Perform the analysis

# Perform the pushover analysis
# Set some parameters
set maxU 10.0;	        # Max displacement
set controlDisp 0.0;
set ok 0;

while {$controlDisp < $maxU && $ok == 0} {
    set ok [analyze 1]
    set controlDisp [nodeDisp 3 1]
    if {$ok != 0} {
	puts "... trying an initial tangent iteration"
	test NormDispIncr 1.0e-8  4000 0
 	algorithm ModifiedNewton -initial
	set ok [analyze 1]
	test NormDispIncr 1.0e-8  10 0
	algorithm Newton
    }
}

if {$ok != 0} {
    puts "Pushover analysis FAILED"
} else {
    puts "Pushover analysis completed SUCCESSFULLY"
}











