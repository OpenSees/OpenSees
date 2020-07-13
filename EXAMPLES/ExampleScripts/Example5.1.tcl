# Units: kip, in, sec

# ----------------------------
# Start of model generation
# ----------------------------

# Create ModelBuilder with 3 dimensions and 6 DOF/node
model BasicBuilder -ndm 3 -ndf 6

# Define geometry
# ---------------

# Set parameters for model geometry
set h 144.0;       # Story height
set by 240.0;      # Bay width in Y-direction
set bx 240.0;      # Bay width in X-direction

# Create nodes
#    tag             X             Y           Z 
node  1  [expr -$bx/2] [expr  $by/2]           0
node  2  [expr  $bx/2] [expr  $by/2]           0
node  3  [expr  $bx/2] [expr -$by/2]           0 
node  4  [expr -$bx/2] [expr -$by/2]           0 

node  5  [expr -$bx/2] [expr  $by/2]          $h 
node  6  [expr  $bx/2] [expr  $by/2]          $h 
node  7  [expr  $bx/2] [expr -$by/2]          $h 
node  8  [expr -$bx/2] [expr -$by/2]          $h 

node 10  [expr -$bx/2] [expr  $by/2] [expr 2*$h]
node 11  [expr  $bx/2] [expr  $by/2] [expr 2*$h] 
node 12  [expr  $bx/2] [expr -$by/2] [expr 2*$h] 
node 13  [expr -$bx/2] [expr -$by/2] [expr 2*$h] 

node 15  [expr -$bx/2] [expr  $by/2] [expr 3*$h] 
node 16  [expr  $bx/2] [expr  $by/2] [expr 3*$h] 
node 17  [expr  $bx/2] [expr -$by/2] [expr 3*$h] 
node 18  [expr -$bx/2] [expr -$by/2] [expr 3*$h]

# Retained nodes for rigid diaphragm
#    tag X Y          Z 
node  9  0 0         $h 
node 14  0 0 [expr 2*$h]
node 19  0 0 [expr 3*$h]

# Set base constraints
#   tag DX DY DZ RX RY RZ
fix  1   1  1  1  1  1  1
fix  2   1  1  1  1  1  1
fix  3   1  1  1  1  1  1
fix  4   1  1  1  1  1  1

# Define rigid diaphragm multi-point constraints
#               normalDir  retained constrained
rigidDiaphragm     3          9     5  6  7  8
rigidDiaphragm     3         14    10 11 12 13
rigidDiaphragm     3         19    15 16 17 18

# Constraints for rigid diaphragm retained nodes
#   tag DX DY DZ RX RY RZ
fix  9   0  0  1  1  1  0
fix 14   0  0  1  1  1  0
fix 19   0  0  1  1  1  0

# Define materials for nonlinear columns
# --------------------------------------
# CONCRETE
set fc 4.0
set Ec [expr 57000.0*sqrt($fc*1000)/1000];
# Core concrete (confined)
#                           tag  f'c  epsc0  f'cu  epscu
uniaxialMaterial Concrete01  1  -5.0 -0.005  -3.5  -0.02

# Cover concrete (unconfined)
uniaxialMaterial Concrete01  2   -$fc -0.002   0.0 -0.006

# STEEL
set fy 60.0;       # Yield stress
set Es 30000.0;    # Young's modulus
# Reinforcing steel
#                        tag fy   E   b
uniaxialMaterial Steel01  3  $fy $Es 0.02

# Column width
set h 18.0
set GJ 1.0e10;
set colSec 1

# Source in a procedure for generating an RC fiber section
source RCsection.tcl

# Call the procedure to generate the column section
#               id  h  b cover core cover steel nBars barArea nfCoreY nfCoreZ nfCoverY nfCoverZ GJ
RCsection  $colSec $h $h   2.5    1     2     3     3    0.79       8       8       10       10 $GJ

# Concrete elastic stiffness

# Define column elements
# ----------------------
set PDelta "OFF"
#set PDelta "ON"

# Geometric transformation for columns
if {$PDelta == "OFF"} {
   #                           tag  vecxz
   geomTransf Linear  1   1 0 0
} else {
   geomTransf PDelta  1   1 0 0
}

# Number of column integration points (sections)
set np 4

# Create the nonlinear column elements
#                           tag ndI ndJ nPts   secID  transf
element nonlinearBeamColumn  1   1   5   $np  $colSec    1
element nonlinearBeamColumn  2   2   6   $np  $colSec    1
element nonlinearBeamColumn  3   3   7   $np  $colSec    1
element nonlinearBeamColumn  4   4   8   $np  $colSec    1

element nonlinearBeamColumn  5   5  10   $np  $colSec    1
element nonlinearBeamColumn  6   6  11   $np  $colSec    1
element nonlinearBeamColumn  7   7  12   $np  $colSec    1
element nonlinearBeamColumn  8   8  13   $np  $colSec    1

element nonlinearBeamColumn  9  10  15   $np  $colSec    1
element nonlinearBeamColumn 10  11  16   $np  $colSec    1
element nonlinearBeamColumn 11  12  17   $np  $colSec    1
element nonlinearBeamColumn 12  13  18   $np  $colSec    1

# Define beam elements
# --------------------

# Define material properties for elastic beams
# Using beam depth of 24 and width of 18
# --------------------------------------------
set Abeam [expr 18.0*24.0];
# "Cracked" second moments of area
set Ibeamzz [expr 0.5*1.0/12.0*18.0*pow(24.0,3)];
set Ibeamyy [expr 0.5*1.0/12.0*24.0*pow(18.0,3)];
set beamSec 2

# Define elastic section for beams
#                  tag    E    A      Iz       Iy     G   J
section Elastic $beamSec $Ec $Abeam $Ibeamzz $Ibeamyy $GJ 1.0

# Geometric transformation for beams
#                tag  vecxz
geomTransf Linear 2   1 1 0

# Number of beam integration points (sections)
set np 3

# Create the beam elements
#                           tag ndI ndJ nPts    secID   transf
element nonlinearBeamColumn  13  5   6   $np  $beamSec     2
element nonlinearBeamColumn  14  6   7   $np  $beamSec     2
element nonlinearBeamColumn  15  7   8   $np  $beamSec     2
element nonlinearBeamColumn  16  8   5   $np  $beamSec     2

element nonlinearBeamColumn  17 10  11   $np  $beamSec     2
element nonlinearBeamColumn  18 11  12   $np  $beamSec     2
element nonlinearBeamColumn  19 12  13   $np  $beamSec     2
element nonlinearBeamColumn  20 13  10   $np  $beamSec     2

element nonlinearBeamColumn  21 15  16   $np  $beamSec     2
element nonlinearBeamColumn  22 16  17   $np  $beamSec     2
element nonlinearBeamColumn  23 17  18   $np  $beamSec     2
element nonlinearBeamColumn  24 18  15   $np  $beamSec     2

# Define gravity loads
# --------------------
# Gravity load applied at each corner node
# 10% of column capacity
set p [expr 0.1*$fc*$h*$h]

# Mass lumped at retained nodes
set g 386.4;            # Gravitational constant
set m [expr (4.0*$p)/$g]

# Rotary inertia of floor about retained node
set i [expr $m*($bx*$bx+$by*$by)/12.0]

# Set mass at the retained nodes
#    tag MX MY MZ RX RY RZ
mass  9  $m $m  0  0  0 $i
mass 14  $m $m  0  0  0 $i
mass 19  $m $m  0  0  0 $i

# Define gravity loads
pattern Plain 1 Constant {
   foreach node {5 6 7 8  10 11 12 13  15 16 17 18} {
      load $node 0.0 0.0 -$p 0.0 0.0 0.0
   }
}

#set rayleigh damping factors
rayleigh 0.0 0.0 0.0 0.0018


# Define earthquake excitation
# ----------------------------
# Set up the acceleration records for Tabas fault normal and fault parallel
set tabasFN "Path -filePath tabasFN.txt -dt 0.02 -factor $g"
set tabasFP "Path -filePath tabasFP.txt -dt 0.02 -factor $g"

# Define the excitation using the Tabas ground motion records
#                         tag dir         accel series args
pattern UniformExcitation  2   1  -accel      $tabasFN
pattern UniformExcitation  3   2  -accel      $tabasFP


# ----------------------- 
# End of model generation
# -----------------------


# ----------------------------
# Start of analysis generation
# ----------------------------

# Create the convergence test
#                tol   maxIter  printFlag
test EnergyIncr 1.0e-8   20         3

# Create the solution algorithm
algorithm Newton

# Create the constraint handler
constraints Transformation

# Create the time integration scheme
#                   gamma beta
integrator Newmark   0.5  0.25 

# Create the system of equation storage and solver
system UmfPack

# Create the DOF numberer
numberer Plain

# Create the transient analysis
analysis Transient

# --------------------------
# End of analysis generation
# --------------------------
# ----------------------------

# Record DOF 1 and 2 displacements at nodes 9, 14, and 19
recorder Node -file Node51.out -time -node 9 14 19 -dof 1 2 disp
recorder plot Node51.out Node9_14_19_Xdisp 10 340 300 300 -columns 1 2 -columns 1 4 -columns 1 6  -dT 1.0

# --------------------------
# End of recorder generation
# --------------------------


# --------------------
# Perform the analysis
# --------------------

# record once at time 0
record

# Analysis duration of 20 seconds
#       numSteps  dt
set ok [analyze   2000   0.01]

if {$ok != 0} {
    puts "analysis FAILED"
} else {
    puts "analysis SUCCESSFUL"
}

wipe
