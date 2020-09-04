# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# 3 Story One-by-One Bay Frame Example 5.1
# ----------------------------------------
#  Reinforced concrete one-bay, three-story frame
#  Distributed vertical load on girder
# 
# Example Objectives
# ------------------
#  3D building with rigid diaphragms
#  Nonlinear beam-column elements
#  Gravity load analysis followed by transient analysis
#
# Units: kips, in, sec
#
# Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
# Date: September 2017

# import the OpenSees Python module
from opensees import *
import RCsection
import math

# ----------------------------
# Start of model generation
# ----------------------------

# remove existing model
wipe()

# create ModelBuilder (with three-dimensions and 6 DOF/node)
model("BasicBuilder", "-ndm",3, "-ndf",6)

# set default units
defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")

# Define geometry
# ---------------

# Set parameters for model geometry
h  = 144.0;      # Story height
by = 240.0;      # Bay width in Y-direction
bx = 240.0;      # Bay width in X-direction

# Create nodes
#       tag    X        Y        Z 
node( 1, -bx/2.0,  by/2.0, 0.0)
node( 2,  bx/2.0,  by/2.0, 0.0)
node( 3,  bx/2.0, -by/2.0, 0.0) 
node( 4, -bx/2.0, -by/2.0, 0.0) 

node( 5, -bx/2.0,  by/2.0, h) 
node( 6,  bx/2.0,  by/2.0, h) 
node( 7,  bx/2.0, -by/2.0, h) 
node( 8, -bx/2.0, -by/2.0, h) 

node(10, -bx/2.0,  by/2.0, 2.0*h)
node(11,  bx/2.0,  by/2.0, 2.0*h) 
node(12,  bx/2.0, -by/2.0, 2.0*h) 
node(13, -bx/2.0, -by/2.0, 2.0*h) 

node(15, -bx/2.0,  by/2.0, 3.0*h) 
node(16,  bx/2.0,  by/2.0, 3.0*h) 
node(17,  bx/2.0, -by/2.0, 3.0*h) 
node(18, -bx/2.0, -by/2.0, 3.0*h)

# Retained nodes for rigid diaphragm
#        tag   X    Y    Z 
node( 9,  0.0, 0.0,     h)
node(14,  0.0, 0.0, 2.0*h)
node(19,  0.0, 0.0, 3.0*h)

# Set base constraints
#      tag DX DY DZ RX RY RZ
fix(1, 1, 1, 1, 1, 1, 1)
fix(2, 1, 1, 1, 1, 1, 1)
fix(3, 1, 1, 1, 1, 1, 1)
fix(4, 1, 1, 1, 1, 1, 1)

# Define rigid diaphragm multi-point constraints
#              normalDir retained constrained
rigidDiaphragm(3,  9,  5,  6,  7,  8)
rigidDiaphragm(3, 14, 10, 11, 12, 13)
rigidDiaphragm(3, 19, 15, 16, 17, 18)

# Constraints for rigid diaphragm retained nodes
#      tag DX DY DZ RX RY RZ
fix( 9, 0, 0, 1, 1, 1, 0)
fix(14, 0, 0, 1, 1, 1, 0)
fix(19, 0, 0, 1, 1, 1, 0)

# Define materials for nonlinear columns
# --------------------------------------
# CONCRETE
fc = 4.0
Ec = 57000.0*math.sqrt(fc*1000.0)/1000.0;

# Core concrete (confined)
#                                 tag  f'c   epsc0  f'cu  epscu
uniaxialMaterial("Concrete01", 1, -5.0, -0.005, -3.5, -0.02)

# Cover concrete (unconfined)
#                                 tag  f'c   epsc0  f'cu  epscu
uniaxialMaterial("Concrete01", 2, -fc, -0.002, 0.0, -0.006)

# STEEL
fy = 60.0;       # Yield stress
Es = 30000.0;    # Young's modulus
# Reinforcing steel 
#                              tag fy  E0  b
uniaxialMaterial("Steel01", 3, fy, Es, 0.02)

# Column parameters
h = 18.0
GJ = 1.0E10
colSec = 1

# Call the RCsection procedure to generate the column section
#                        id  h  b cover core cover steel nBars barArea nfCoreY nfCoreZ nfCoverY nfCoverZ GJ
RCsection.create(colSec, h, h, 2.5, 1,    2,    3,    3,   0.79,     8,      8,      10,      10,   GJ)

# Define column elements
# ----------------------
PDelta = "OFF"
#PDelta = "ON"

# Geometric transformation for columns
if (PDelta == "OFF"):
   geomTransf("Linear", 1, 1.0, 0.0, 0.0)
else:
   geomTransf("PDelta", 1, 1.0, 0.0, 0.0)

# Number of column integration points (sections)
np = 4
beamIntegration("Lobatto", colSec, colSec, np)

# Create the nonlinear column elements
eleType = "forceBeamColumn"
#                   tag ndI ndJ transfTag integrationTag
element(eleType, 1, 1, 5, 1, colSec)
element(eleType, 2, 2, 6, 1, colSec)
element(eleType, 3, 3, 7, 1, colSec)
element(eleType, 4, 4, 8, 1, colSec)

element(eleType, 5, 5, 10, 1, colSec)
element(eleType, 6, 6, 11, 1, colSec)
element(eleType, 7, 7, 12, 1, colSec)
element(eleType, 8, 8, 13, 1, colSec)

element(eleType,  9, 10, 15, 1, colSec)
element(eleType, 10, 11, 16, 1, colSec)
element(eleType, 11, 12, 17, 1, colSec)
element(eleType, 12, 13, 18, 1, colSec)

# Define beam elements
# --------------------
# Define material properties for elastic beams
# Using beam depth of 24 and width of 18
Abeam = 18.0*24.0
# "Cracked" second moments of area
Ibeamzz = 0.5*1.0/12.0*18.0*pow(24.0,3)
Ibeamyy = 0.5*1.0/12.0*24.0*pow(18.0,3)
beamSec = 2

# Define elastic section for beams
#                       tag     E    A      Iz       Iy     G    J
section("Elastic", beamSec, Ec, Abeam, Ibeamzz, Ibeamyy, GJ, 1.0)

# Geometric transformation for beams
geomTransf("Linear", 2, 1.0, 1.0, 0.0)

# Number of beam integration points (sections)
np = 3
beamIntegration("Lobatto", beamSec, beamSec, np)

# Create the beam elements
eleType = "forceBeamColumn"
#                   tag ndI ndJ transfTag integrationTag
element(eleType, 13, 5, 6, 2, beamSec)
element(eleType, 14, 6, 7, 2, beamSec)
element(eleType, 15, 7, 8, 2, beamSec)
element(eleType, 16, 8, 5, 2, beamSec)

element(eleType, 17, 10, 11, 2, beamSec)
element(eleType, 18, 11, 12, 2, beamSec)
element(eleType, 19, 12, 13, 2, beamSec)
element(eleType, 20, 13, 10, 2, beamSec)

element(eleType, 21, 15, 16, 2, beamSec)
element(eleType, 22, 16, 17, 2, beamSec)
element(eleType, 23, 17, 18, 2, beamSec)
element(eleType, 24, 18, 15, 2, beamSec)

# Define gravity loads
# --------------------
# Gravity load applied at each corner node
# 10% of column capacity
p = 0.1*fc*h*h

# Mass lumped at retained nodes
m = (4.0*p)/g

# Rotary inertia of floor about retained node
i = m*(bx*bx + by*by)/12.0

# Set mass at the retained nodes
#        tag MX MY MZ   RX   RY   RZ
mass( 9, m, m, 0.0, 0.0, 0.0, i)
mass(14, m, m, 0.0, 0.0, 0.0, i)
mass(19, m, m, 0.0, 0.0, 0.0, i)

# Define gravity loads
# create a Constant TimeSeries
timeSeries("Constant", 1)
# create a Plain load pattern
pattern("Plain", 1, 1, "-fact", 1.0)

for i in [5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18]:
    load(i, 0.0, 0.0, -p, 0.0, 0.0, 0.0)

# set rayleigh damping factors
rayleigh(0.0, 0.0, 0.0, 0.0018)

# Define earthquake excitation
# ----------------------------
dt = 0.02
# Set up the acceleration records for Tabas fault normal and fault parallel
timeSeries("Path", 2, "-filePath", "tabasFN.txt", "-dt", dt, "-factor", g)
timeSeries("Path", 3, "-filePath", "tabasFP.txt", "-dt", dt, "-factor", g)

# Define the excitation using the Tabas ground motion records
#                         tag dir         accel series args
pattern("UniformExcitation", 2, 1, "-accel", 2)
pattern("UniformExcitation", 3, 2, "-accel", 3)

# print model
#printModel()
printModel("-JSON", "-file", "Example5.1.json")

# ----------------------- 
# End of model generation
# -----------------------


# ----------------------------
# Start of analysis generation
# ----------------------------

# create the system of equation
system("UmfPack")

# create the DOF numberer
numberer("Plain")

# create the constraint handler
constraints("Transformation")

# create the convergence test
test("EnergyIncr", 1.0E-8, 20)

# create the solution algorithm, a Newton-Raphson algorithm
algorithm("Newton")

# create the integration scheme, the Newmark with gamma=0.5 and beta=0.25
integrator("Newmark", 0.5, 0.25) 

# create the analysis object 
analysis("Transient")

# --------------------------
# End of analysis generation
# --------------------------


# ----------------------------
# Start of recorder generation
# ----------------------------

# Record DOF 1 and 2 displacements at nodes 9, 14, and 19
recorder("Node", "-file", "Node51.out", "-time", "-node", 9, 14, 19, "-dof", 1, 2, "disp")
#recorder("plot", "Node51.out", "Node9_14_19_Xdisp", 10, 340, 300, 300, "-columns", 1, 2, "-columns", 1, 4, "-columns", 1, 6, "-dT", 1.0)

# --------------------------
# End of recorder generation
# --------------------------


# --------------------
# Perform the analysis
# --------------------

# record once at time 0
record()

# Analysis duration of 20 seconds
#              numSteps dt
ok = analyze(2000, 0.01)

if (ok != 0):
    print("analysis FAILED")
else:
    print("analysis SUCCESSFUL")

wipe()
