# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# Basic Truss Example 1.1
# -----------------------
#  2d 3 Element Elastic Truss
#  Single Nodal Load, Static Analysis
# 
# Example Objectives
# ------------------
#  Simple Introduction to OpenSees
# 
# Units: kips, in, sec
#
# Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
# Date: June 2017

# import the OpenSees Python module
from opensees import *

# ------------------------------
# Start of model generation
# ------------------------------

# remove existing model
wipe()

# create ModelBuilder (with two-dimensions and 2 DOF/node)
model("BasicBuilder", "-ndm",2, "-ndf",2)

# set default units
defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")

# create nodes & add to Domain - command: node nodeId xCrd yCrd
#node(1, 0.0,    0.0, "-disp",0.0,0.0, "-vel", 0.0,0.0, "-mass", 0.0,0.0)
node(1, 0.0,    0.0)
node(2, 144.0,  0.0)
node(3, 168.0,  0.0)
node(4,  72.0, 96.0)

# set the boundary conditions - command: fix nodeID xRestrnt? yRestrnt?
fix(1, 1, 1)
fix(2, 1, 1)
fix(3, 1, 1)

# Define materials for truss elements
# -----------------------------------
# Create Elastic material prototype - command: uniaxialMaterial Elastic matID E
uniaxialMaterial("Elastic", 1, 3000.0)

# Define elements
# ---------------
# Create truss elements - command: element truss trussID node1 node2 A matID
element("truss", 1, 1, 4, 10.0, 1)
element("truss", 2, 2, 4,  5.0, 1)
element("truss", 3, 3, 4,  5.0, 1)

# Define loads
# ------------
# create a Linear TimeSeries (load factor varies linearly with time) - command: timeSeries Linear $tag
timeSeries("Linear", 1)

# create a Plain load pattern - command: pattern Plain $tag $timeSeriesTag { $loads }
pattern("Plain", 1, 1, "-fact", 1.0)
# create the nodal load - command: load nodeID xForce yForce
load(4, 100.0, -50.0)

# print model
#printModel()
printModel("-JSON", "-file", "Example1.1.json")

# ------------------------------
# End of model generation
# ------------------------------


# ------------------------------
# Start of analysis generation
# ------------------------------

# create the system of equation, a SPD using a band storage scheme
system("BandSPD")

# create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer("RCM")

# create the constraint handler, a Plain handler is used as homo constraints
constraints("Plain")

# create the solution algorithm, a Linear algorithm is created
algorithm("Linear")

# create the integration scheme, the LoadControl scheme using steps of 1.0
integrator("LoadControl", 1.0)

# create the analysis object 
analysis("Static")

# ------------------------------
# End of analysis generation
# ------------------------------


# ------------------------------
# Start of recorder generation
# ------------------------------

# create a Recorder object for the nodal displacements at node 4
recorder("Node", "-file", "example.out", "-time", "-node", 4, "-dof", 1, 2, "disp")

# create a recorder for element forces, one in global and the other local system
recorder("Element", "-file", "eleGlobal.out", "-time", "-ele", 1, 2, 3, "forces")
recorder("Element", "-file", "eleLocal.out", "-time", "-ele", 1, 2, 3, "basicForces")

# ------------------------------
# End of recorder generation
# ------------------------------


# ------------------------------
# Finally perform the analysis
# ------------------------------

# perform the analysis
analyze(1)


# ------------------------------
# Print Stuff to Screen
# ------------------------------

# print the current state at node 4 and at all elements
#print("node 4 displacement: ", nodeDisp(4))
printModel("node", 4)
printModel("ele")
wipe()
