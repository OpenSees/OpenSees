# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# 2 Story Multi Bay Frame Example 4.1
# -----------------------------------
#  Reinforced concrete multi-bay, two-story frame
#  Distributed vertical load on girder
# 
# Example Objectives
# ------------------
#  Nonlinear beam-column elements
#  Gravity load analysis followed by pushover analysis
#  Demonstrate scripting for the algorithmic level
#
# Units: kips, in, sec
#
# Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
# Date: August 2017

# import the OpenSees Python module
from opensees import *
import math

# ------------------------------
# Start of model generation
# ------------------------------

# Parameter identifying the number of bays
numBay = 3

# remove existing model
wipe()

# create ModelBuilder (with two-dimensions and 3 DOF/node)
model("BasicBuilder", "-ndm",2, "-ndf",3)

# set default units
defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")

# Create nodes
# ------------
# Set parameters for overall model geometry
bayWidth = 288.0
m = 0.1
nodeID = 1

# Define nodes
for i in range(numBay+1):
    xDim = i * bayWidth
    
    #         tag       X      Y
    node(nodeID,   xDim,   0.0)
    node(nodeID+1, xDim, 180.0, "-mass", m, m, 0.0)
    node(nodeID+2, xDim, 324.0, "-mass", m, m, 0.0)

    nodeID += 3

# Fix supports at base of columns
for i in range(numBay+1):
    #       node  DX DY RZ
    fix(i*3+1, 1, 1, 1)

# Define materials for nonlinear columns
# ------------------------------------------
# CONCRETE                        tag  f'c    ec0    f'cu   ecu
# Core concrete (confined)
uniaxialMaterial("Concrete01", 1, -6.0, -0.004, -5.0, -0.014)
# Cover concrete (unconfined)
uniaxialMaterial("Concrete01", 2, -5.0, -0.002, -0.0, -0.006)

# STEEL
# Reinforcing steel 
fy = 60.0;      # Yield stress
E = 30000.0;    # Young's modulus
#                              tag fy  E0  b
uniaxialMaterial("Steel01", 3, fy, E, 0.015)

# Define cross-section for nonlinear columns
# ------------------------------------------
# Interior column section
section("Fiber", 1)
#                mat nfIJ nfJK   yI     zI     yJ     zJ     yK     zK     yL     zL
patch("quad", 2,   1,  12, -11.5,  10.0, -11.5, -10.0,  11.5, -10.0,  11.5,  10.0)
patch("quad", 1,   1,  14, -13.5, -10.0, -13.5, -12.0,  13.5, -12.0,  13.5, -10.0)
patch("quad", 1,   1,  14, -13.5,  12.0, -13.5,  10.0,  13.5,  10.0,  13.5,  12.0)
patch("quad", 1,   1,   2, -13.5,  10.0, -13.5, -10.0, -11.5, -10.0, -11.5,  10.0)
patch("quad", 1,   1,   2,  11.5,  10.0,  11.5, -10.0,  13.5, -10.0,  13.5,  10.0)
#                    mat nBars area    yI    zI     yF    zF
layer("straight", 3,   6,  1.56, -10.5,  9.0, -10.5, -9.0)
layer("straight", 3,   6,  1.56,  10.5,  9.0,  10.5, -9.0)
# define beam integration
np = 4;  # number of integration points along length of element
beamIntegration("Lobatto", 1, 1, np)

# Exterior column section
section("Fiber", 2)
patch("quad", 2, 1, 10, -10.0,  10.0, -10.0, -10.0,  10.0, -10.0,  10.0,  10.0)
patch("quad", 1, 1, 12, -12.0, -10.0, -12.0, -12.0,  12.0, -12.0,  12.0, -10.0)
patch("quad", 1, 1, 12, -12.0,  12.0, -12.0,  10.0,  12.0,  10.0,  12.0,  12.0)
patch("quad", 1, 1,  2, -12.0,  10.0, -12.0, -10.0, -10.0, -10.0, -10.0,  10.0)
patch("quad", 1, 1,  2,  10.0,  10.0,  10.0, -10.0,  12.0, -10.0,  12.0,  10.0)
layer("straight", 3, 6, 0.79, -9.0, 9.0, -9.0, -9.0)
layer("straight", 3, 6, 0.79,  9.0, 9.0,  9.0, -9.0)
# define beam integration
beamIntegration("Lobatto", 2, 2, np)

# Girder section
section("Fiber", 3)
patch("quad", 1, 1, 12, -12.0, 9.0, -12.0, -9.0, 12.0, -9.0, 12.0, 9.0)
layer("straight", 3, 4, 1.0, -9.0, 9.0, -9.0, -9.0)
layer("straight", 3, 4, 1.0,  9.0, 9.0,  9.0, -9.0)
# define beam integration
beamIntegration("Lobatto", 3, 3, np)

# Define column elements
# ----------------------
# Geometric transformation
geomTransf("Linear", 1)

beamID = 1
eleType = "forceBeamColumn"

# Define elements
for i in range(numBay+1):
    # set some parameters
    iNode = i*3 + 1
    jNode = i*3 + 2

    for j in range(1, 3):
        # add the column element (secId == 2 if external, 1 if internal column)
        if (i == 0):
            element(eleType, beamID, iNode, jNode, 1, 2)
        elif (i == numBay):
            element(eleType, beamID, iNode, jNode, 1, 2)
        else:
            element(eleType, beamID, iNode, jNode, 1, 1)
        
        # increment the parameters
        iNode += 1
        jNode += 1
        beamID += 1

# Define beam elements
# ----------------------
# Geometric transformation
geomTransf("Linear", 2)

# Define elements
for j in range(1, 3):
    # set some parameters
    iNode = j + 1
    jNode = iNode + 3

    for i in range(1, numBay+1):
        element(eleType, beamID, iNode, jNode, 2, 3)
        
        # increment the parameters
        iNode += 3
        jNode += 3
        beamID += 1
 
# Define gravity loads
# --------------------
# Constant gravity load
P = -192.0

# create a Linear TimeSeries
timeSeries("Linear", 1)

# create a Plain load pattern
pattern("Plain", 1, 1, "-fact", 1.0)

# Create nodal loads at nodes 
for i in range(numBay+1):
    # set some parameters
    node1 = i*3 + 2
    node2 = node1 + 1
    
    if (i == 0):
        load(node1, 0.0, P,     0.0) 
        load(node2, 0.0, P/2.0, 0.0)
    elif (i == numBay):
        load(node1, 0.0, P,     0.0)
        load(node2, 0.0, P/2.0, 0.0)
    else:
        load(node1, 0.0, 2.0*P, 0.0)
        load(node2, 0.0, P,     0.0)

# print model
#printModel()
printModel("-JSON", "-file", "Example4.1.json")

# ------------------------------
# End of model generation
# ------------------------------


# --------------------------------------------------
# Start of analysis generation for gravity analysis
# --------------------------------------------------

# create the system of equation
system("BandGeneral")

# create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer("RCM")

# create the constraint handler, a Plain handler is used as homo constraints
constraints("Plain")

# Create the convergence test, the norm of the residual with a tolerance of 
# 1e-12 and a max number of iterations of 10
test("NormDispIncr", 1.0E-8, 10, 0)

# create the solution algorithm, a Newton-Raphson algorithm
algorithm("Newton")

# create the integration scheme, the LoadControl scheme using steps of 0.1
integrator("LoadControl", 0.1)

# create the analysis object 
analysis("Static")

# ------------------------------------------------
# End of analysis generation for gravity analysis
# ------------------------------------------------


# ------------------------------
# Perform gravity load analysis
# ------------------------------

# initialize the model, done to set initial tangent
initialize()

# perform the gravity load analysis, requires 10 steps to reach the load level
analyze(10)

print("Gravity load analysis completed\n")

# set gravity loads to be const and set pseudo time to be 0.0
# for start of lateral load analysis
loadConst("-time", 0.0)


# ------------------------------
# Add lateral loads 
# ------------------------------

# Reference lateral load for pushover analysis
H = 10.0

# Set lateral load pattern with a Linear TimeSeries
pattern("Plain", 2, 1, "-fact", 1.0)
load(2, H/2.0, 0.0, 0.0)
load(3, H,     0.0, 0.0)


# ------------------------------
# Start of recorder generation
# ------------------------------

# Create a recorder which writes to Node.out and prints
# the current load factor (pseudo-time) and dof 1 displacements at node 2 & 3
recorder("Node", "-file", "Node41.out", "-time", "-node", 2, 3, "-dof", 1, "disp")

# Source in some commands to display the model
# comment out one of lines
#displayMode = "displayON"
displayMode = "displayOFF"

if (displayMode == "displayON"):  
    # a window to plot the nodal displacements versus load for node 3  
    recorder("plot", "Node41.out", "Node_3_Xdisp", 10, 340, 300, 300, "-columns", 3, 1, "-dT", 0.1)

# ------------------------------
# End of recorder generation
# ------------------------------


# ------------------------------
# Start of lateral load analysis
# ------------------------------

# Perform an eigenvalue analysis
lam = eigen(2)
Tstart = 2.0*math.pi/math.sqrt(lam[0])
print("Fundamental period at start of pushover analysis: ", Tstart, "sec\n")

# Change the integrator to take a min and max load increment
integrator("LoadControl", 1.0, 4, 0.02, 2.0)

# record once at time 0
record()

# Perform the pushover analysis
# Set some parameters
maxU = 10.0;	        # Max displacement
controlDisp = 0.0
ok = analyze(1)
while ((ok == 0) and (controlDisp < maxU)):
    ok = analyze(1)
    controlDisp = nodeDisp(3, 1)
    if (ok != 0):
        print("... trying an initial tangent iteration")
        test("NormDispIncr", 1.0e-8, 4000, 0)
        algorithm("ModifiedNewton", "-initial")
        ok = analyze(1)
        test("NormDispIncr", 1.0e-8, 10, 0)
        algorithm("Newton")

# Print a message to indicate if analysis successful or not
if (ok == 0):
    print("\nPushover analysis completed SUCCESSFULLY\n")
else:
    print("\nPushover analysis FAILED\n")

# Print the state at node 3
printModel("node", 3)
wipe()
