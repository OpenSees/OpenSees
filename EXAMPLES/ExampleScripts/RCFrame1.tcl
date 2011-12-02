# RCFrame1.tcl: R/C two story, two bay frame
# Units: kip, in
# MHS, Sept 1999
#   email: mhscott@ce.berkeley.edu
#
# Pushover analysis
# No TCL constructs
#
#    _________________________   _
#   |            |            |  
#   |            |            | 12'
#   |            |            |
#   |____________|____________|  _
#   |            |            |
#   |            |            |
#   |            |            | 15'
#   |            |            |
#   |            |            |  _
#  ===          ===          ===
#   |     24'    |     24'    |
#
#
# NOTE: to RUN this example, run the g3 interpreter and 
#       type the command: source RCFrame1.tcl
#
# $Revision: 1.1.1.1 $
# $Date: 2000-09-15 08:23:09 $
# $Source: /usr/local/cvs/OpenSees/EXAMPLES/ExampleScripts/RCFrame1.tcl,v $

# comment out one of lines
set displayMode "displayON"
#set displayMode "displayOFF"

# Define the model builder
model BasicBuilder -ndm 2 -ndf 3

# Define nodes
#    tag  X   Y
node  1   0   0
node  2   0 180
node  3   0 324
node  4 288   0
node  5 288 180
node  6 288 324
node  7 576   0
node  8 576 180
node  9 576 324

# Single point constraints
#   node DX DY RZ
fix   1   1  1  1
fix   4   1  1  1
fix   7   1  1  1

# Define materials
# Cover concrete
#                  tag -f'c  -epsco  -f'cu -epscu
uniaxialMaterial Concrete01 1 -4.00  -0.002    0.0 -0.006
# Core concrete
uniaxialMaterial Concrete01 2 -5.20  -0.005  -4.70  -0.02

# Steel model with softening
#                   tag   E      epsY
uniaxialMaterial ElasticPP   4  29500  0.002333
uniaxialMaterial ElasticPP   5   1100  0.03
uniaxialMaterial Elastic     6   -600
#                   tag mat1 mat2 mat3 ...
uniaxialMaterial Parallel    3    4    5    6

# Interior column section
section fiberSec 1 {
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
section fiberSec 2 {
   patch quadr 2 1 10 -10  10 -10 -10  10 -10  10  10
   patch quadr 1 1 12 -12 -10 -12 -12  12 -12  12 -10
   patch quadr 1 1 12 -12  12 -12  10  12  10  12  12
   patch quadr 1 1  2 -12  10 -12 -10 -10 -10 -10  10
   patch quadr 1 1  2  10  10  10 -10  12 -10  12  10
   layer straight 3 6 0.79 -9 9 -9 -9
   layer straight 3 6 0.79  9 9  9 -9
}

# Girder section
section fiberSec 3 {
   patch quadr 1 1 12 -12 9 -12 -9 12 -9 12 9
   layer straight 3 4 1.00 -9 9 -9 -9
   layer straight 3 4 1.00  9 9  9 -9
}

# Number of integration points
set nP 4

# Geometric transformation
geomTransf Linear 1

# Define elements
# Columns
#                           tag ndI ndJ  nPts secID transf
element nonlinearBeamColumn  1   1   2    $nP   2      1
element nonlinearBeamColumn  2   2   3    $nP   2      1
element nonlinearBeamColumn  3   4   5    $nP   1      1
element nonlinearBeamColumn  4   5   6    $nP   1      1
element nonlinearBeamColumn  5   7   8    $nP   2      1
element nonlinearBeamColumn  6   8   9    $nP   2      1

# Beams
element nonlinearBeamColumn  7   2   5    $nP   3      1
element nonlinearBeamColumn  8   5   8    $nP   3      1
element nonlinearBeamColumn  9   3   6    $nP   3      1
element nonlinearBeamColumn 10   6   9    $nP   3      1

# Constant gravity loads
pattern Plain 1 Linear {
   #    node  FX    FY   MZ
   load   2  0.0 -192.0 0.0 -const
   load   5  0.0 -384.0 0.0 -const
   load   8  0.0 -192.0 0.0 -const
   load   3  0.0  -96.0 0.0 -const
   load   6  0.0 -192.0 0.0 -const
   load   9  0.0  -96.0 0.0 -const
}

# Load control with variable load steps
#                     init Jd min max
integrator LoadControl  1   3  .2  1

# Convergence test
#                  tolerance maxIter displayCode
test NormDispIncr  1.0e-06     10         1

# Solution algorithm
algorithm Newton

# DOF numberer
numberer RCM

# Cosntraint handler
constraints Transformation

# System of equations solver
system SparseGeneral -piv

# Analysis for gravity load
analysis Static

# Perform the gravity load analysis
analyze 1


# Reference lateral loads
pattern Plain 2 Linear {
   #    node  FX   FY   MZ
   load   2 10.0  0.0  0.0
   load   3 10.0  0.0  0.0
}

# Create a recorder which writes to Node.out and prints
# the current load factor (pseudo-time) and dof 1 displacements at node 2 & 3
recorder Node Node.out disp -time -node 2 3 -dof 1

# Create another recorder which writes to Element.out and prints
# the section deformations at the base of the columns  (elements 1,3,5)
recorder Element 1 3 5 -file Element.out section 1 deformations

# Source in some g3 commands to display the model
if {$displayMode == "displayON"} {
    # a window to plot the nodal displacements versus load for node 3
    recorder plot Node.out Node3Xdisp 10 340 300 300 -columns 3 1
	
    # a window to show the displayed shape
    source RCFrameDisplay.tcl 
}

# Create an analysis object
analysis Static

# Perform 40 analysis steps
analyze 40

# Switch integrator to displacement control
integrator DisplacementControl 3 1 .001 3 .001 1.0

# perform 25 analysis steps
analyze 25

