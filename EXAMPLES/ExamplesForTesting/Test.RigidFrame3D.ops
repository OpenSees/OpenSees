# 3-D model of rigid floor diaphragm
#     Three stories
#     Elastic RC columns
#     Nonlinear RC beams
#     Non-linear dynamic analysis
#
# Written: MHS
#   email: mhscott@ce.berkeley.edu
# Date: 3 Dec 1999
#
# Units: kip, in
#
# NOTE: to RUN this example, run the g3 interpreter and 
#       type the command: source RigidFrame3D.ops


# Define the model builder
model BasicBuilder -ndm 3 -ndf 6

# Gravitational constant
set g 386.4

# Story height
set h 144

# Width in Y direction
set by 240

# Width in X direction
set bx 240

# Dead Load (ksi)
set dl [expr 0.3/144]

# Live Load (ksi)
set ll [expr 0.2/144]

# Define nodes
#                   X             Y Z 
node  1 [expr -$by/2] [expr  $bx/2] 0
node  2 [expr  $by/2] [expr  $bx/2] 0 
node  3 [expr  $by/2] [expr -$bx/2] 0 
node  4 [expr -$by/2] [expr -$bx/2] 0 

node  5 [expr -$by/2] [expr  $bx/2] $h 
node  6 [expr  $by/2] [expr  $bx/2] $h 
node  7 [expr  $by/2] [expr -$bx/2] $h 
node  8 [expr -$by/2] [expr -$bx/2] $h 

node 10 [expr -$by/2] [expr  $bx/2] [expr 2*$h]
node 11 [expr  $by/2] [expr  $bx/2] [expr 2*$h] 
node 12 [expr  $by/2] [expr -$bx/2] [expr 2*$h] 
node 13 [expr -$by/2] [expr -$bx/2] [expr 2*$h] 

node 15 [expr -$by/2] [expr  $bx/2] [expr 3*$h] 
node 16 [expr  $by/2] [expr  $bx/2] [expr 3*$h] 
node 17 [expr  $by/2] [expr -$bx/2] [expr 3*$h] 
node 18 [expr -$by/2] [expr -$bx/2] [expr 3*$h]

# Set base constraints
#     DX DY DZ RX RY RZ
fix 1  1  1  1  1  1  1
fix 2  1  1  1  1  1  1
fix 3  1  1  1  1  1  1
fix 4  1  1  1  1  1  1

# Mass lumped at retained node
set m [expr ($dl+$ll)*$bx*$by/$g]

# Rotary inertia of floor about retained node
set i [expr $m*($bx^2+$by^2)/12.0]

# Retained nodes for rigid diaphragm
node  9 0 0         $h  -mass $m $m 0 0 0 $i
node 14 0 0 [expr 2*$h] -mass $m $m 0 0 0 $i
node 19 0 0 [expr 3*$h] -mass $m $m 0 0 0 $i

# Define rigid diaphragm constraints
#               normalDir  retained     constrained
rigidDiaphragm     3         9      5  6  7  8
rigidDiaphragm     3        14     10 11 12 13
rigidDiaphragm     3        19     15 16 17 18

# Constraints for rigid diaphragm retained nodes
fix  9  0 0 1 1 1 0
fix 14  0 0 1 1 1 0
fix 19  0 0 1 1 1 0

# Define materials
# Core concrete
#                  tag -f'c -epsc0 -f'cu -epscu
uniaxialMaterial Concrete01 1  -5.0 -0.005  -3.5  -0.02

# Cover concrete
uniaxialMaterial Concrete01 2  -4.0 -0.002   0.0 -0.006

# Steel model
#                        tag fy   E     b
uniaxialMaterial Steel01  3  60 30000 0.02

# Torsional stiffness
set G 2000

# Column polar moment of inertia
set Jcol 5000

# Column cross-sectional area
set Acol 576

# Column moment of inertia
set Icol 28000

# Elastic modulus of concrete
set E 4000

# Geometric transformation
#                tag  vecxz
geomTransf Linear 1   1 0 0

# Define elements
# COLUMNS
#                         tag ndI ndJ  A    E  G   Jx    Iy    Iz  transf
element elasticBeamColumn  1   1   5 $Acol $E $G $Jcol $Icol $Icol    1
element elasticBeamColumn  2   2   6 $Acol $E $G $Jcol $Icol $Icol    1
element elasticBeamColumn  3   3   7 $Acol $E $G $Jcol $Icol $Icol    1
element elasticBeamColumn  4   4   8 $Acol $E $G $Jcol $Icol $Icol    1

element elasticBeamColumn  5   5  10 $Acol $E $G $Jcol $Icol $Icol    1
element elasticBeamColumn  6   6  11 $Acol $E $G $Jcol $Icol $Icol    1
element elasticBeamColumn  7   7  12 $Acol $E $G $Jcol $Icol $Icol    1
element elasticBeamColumn  8   8  13 $Acol $E $G $Jcol $Icol $Icol    1

element elasticBeamColumn  9  10  15 $Acol $E $G $Jcol $Icol $Icol    1
element elasticBeamColumn 10  11  16 $Acol $E $G $Jcol $Icol $Icol    1
element elasticBeamColumn 11  12  17 $Acol $E $G $Jcol $Icol $Icol    1
element elasticBeamColumn 12  13  18 $Acol $E $G $Jcol $Icol $Icol    1

# Number of Lobatto integration points
set np 4

# Beam polar moment of inertia
set Jbeam 4000

# Torsional stiffness
set GJ [expr $G*$Jbeam]

# Source in the procedure definition for generating an RC section
source Test.RCsection.ops

# Call the procedure to generate a column section
#          id  h  b cover core cover steel nBars barArea nfCoreY nfCoreZ nfCoverY nfCoverZ
RCsection   1 18 15   2.5    1     2     3     3    0.79       8       8       10       10

# Linear elastic torion
uniaxialMaterial Elastic 10 $GJ

# Attach torsion to the RC section
#                 tag uniTag uniCode       secTag
section Aggregator 2    10      T    -section 1

# Geometric transformation
#                tag  vecxz
geomTransf Linear 2   0 0 1

# BEAMS
#                           tag ndI ndJ nPts secID transf
element nonlinearBeamColumn  13  5   6   $np   2     2
element nonlinearBeamColumn  14  6   7   $np   2     2
element nonlinearBeamColumn  15  7   8   $np   2     2
element nonlinearBeamColumn  16  8   5   $np   2     2

element nonlinearBeamColumn  17 10  11   $np   2     2
element nonlinearBeamColumn  18 11  12   $np   2     2
element nonlinearBeamColumn  19 12  13   $np   2     2
element nonlinearBeamColumn  20 13  10   $np   2     2

element nonlinearBeamColumn  21 15  16   $np   2     2
element nonlinearBeamColumn  22 16  17   $np   2     2
element nonlinearBeamColumn  23 17  18   $np   2     2
element nonlinearBeamColumn  24 18  15   $np   2     2

# Gravity load applied at each corner node
set p [expr -$m*$g/4]

# Define gravity loads
pattern Plain 1 Linear {
   load  5  0.0 0.0 $p 0.0 0.0 0.0 -const
   load  6  0.0 0.0 $p 0.0 0.0 0.0 -const
   load  7  0.0 0.0 $p 0.0 0.0 0.0 -const
   load  8  0.0 0.0 $p 0.0 0.0 0.0 -const

   load 10  0.0 0.0 $p 0.0 0.0 0.0 -const
   load 11  0.0 0.0 $p 0.0 0.0 0.0 -const
   load 12  0.0 0.0 $p 0.0 0.0 0.0 -const
   load 13  0.0 0.0 $p 0.0 0.0 0.0 -const

   load 15  0.0 0.0 $p 0.0 0.0 0.0 -const
   load 16  0.0 0.0 $p 0.0 0.0 0.0 -const
   load 17  0.0 0.0 $p 0.0 0.0 0.0 -const
   load 18  0.0 0.0 $p 0.0 0.0 0.0 -const
}

set accelSeries "Path -filePath Test.tabasFN.txt -dt 0.02 -factor $g"

# Define the ground motion excitation using Tabas fault parallel and fault normal records
#                         tag dir         accel series args
pattern UniformExcitation  2   1  -accel    $accelSeries
pattern UniformExcitation  3   2  -accel    $accelSeries

# Record DOF 1 and 2 displacements at nodes 9, 14, and 19
recorder Node RigidFrame3D.out disp -time -node 9 14 19 -dof 1 2

# Source in commands to display the structure
#source RigidFrame3Ddisplay.ops

# Convergence test
#                tol   maxIter  printFlag
test EnergyIncr 1.0e-8   20         1

# Solution algorithm
algorithm Newton

# System of equations solver
system SparseGeneral -piv

# Transient integrator
#                   gamma beta
integrator Newmark   0.5  0.25

# DOF numberer
numberer RCM

# Constraint handler
#                      aSP    aMP
#constraints Lagrange   1.0e6    1.0e6
constraints Transformation

# Transient analysis
#                   dt    T
analysis Transient

# Perform the analysis
analyze 5 0.01

wipe





