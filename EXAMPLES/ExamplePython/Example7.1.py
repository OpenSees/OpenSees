# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# 3D Shell Structure Example 7.1
# ------------------------------
#  Shell roof modeled with three
#  dimensional linear shell elements
# 
# Example Objectives
# ------------------
#  test linear-elastic shell element
#  free vibration analysis starting from static deflection
#
# Units: kips, in, sec
#
# Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
# Date: September 2017

# import the OpenSees Python module
from openseespy.opensees import *

# ----------------------------
# Start of model generation
# ----------------------------
# remove existing model
wipe()

# create ModelBuilder (with three-dimensions and 6 DOF/node)
model("BasicBuilder", "-ndm",3, "-ndf",6)

# set default units
defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")

# Define the section
# ------------------
#                                       secTag  E     nu     h    rho
section("ElasticMembranePlateSection", 1, 3.0E3, 0.25, 1.175, 1.27)

# Define geometry
# ---------------
# these should both be even
nx = 8
ny = 2

# loaded nodes
mid   = int(((nx+1)*(ny+1) + 1)/2)
side1 = int((nx+2)/2) 
side2 = int((nx+1)*(ny+1) - side1 + 1)

# generate the nodes and elements
#          numX numY startNode startEle eleType eleArgs? coords?
block2D(nx, ny, 1, 1,
            "ShellMITC4", 1,
            1, -20.0,  0.0,  0.0,
            2, -20.0,  0.0, 40.0,
            3,  20.0,  0.0, 40.0,
            4,  20.0,  0.0,  0.0,
            5, -10.0, 10.0, 20.0, 
            7,  10.0, 10.0, 20.0,   
            9,   0.0, 10.0, 20.0)

# define the boundary conditions
# rotation free about x-axis (remember right-hand-rule)
fixZ( 0.0, 1, 1, 1, 0, 1, 1)
fixZ(40.0, 1, 1, 1, 0, 1, 1)  

# create a Linear time series
timeSeries("Linear", 1)
# add some loads
pattern("Plain", 1, 1, "-fact", 1.0)
load(mid  , 0.0, -0.50, 0.0, 0.0, 0.0, 0.0)
load(side1, 0.0, -0.25, 0.0, 0.0, 0.0, 0.0)
load(side2, 0.0, -0.25, 0.0, 0.0, 0.0, 0.0)

# print model
#printModel()
printModel("-JSON", "-file", "Example7.1.json")

# ----------------------- 
# End of model generation
# -----------------------


# ------------------------
# Start of static analysis
# ------------------------

# Load control with variable load steps
#                            init  Jd  min  max
integrator("LoadControl", 1.0, 1, 1.0, 10.0)

# Convergence test
#                  tolerance maxIter displayCode
test("EnergyIncr", 1.0E-10, 20, 0)

# Solution algorithm
algorithm("Newton")

# DOF numberer
numberer("RCM")

# Constraint handler
constraints("Plain") 

# System of equations solver
system("SparseGeneral", "-piv")
#system("ProfileSPD")

# Analysis for gravity load
analysis("Static") 

# Perform the gravity load analysis
analyze(5)

# --------------------------
# End of static analysis
# --------------------------


# ----------------------------
# Start of recorder generation
# ----------------------------

recorder("Node", "-file", "Node.out", "-time", "-node", mid, "-dof", 2, "disp")
#recorder("plot", "Node.out", "CenterNodeDisp", 625, 10, 625, 450, "-columns", 1, 2)

# create the display
#recorder("display", "shellDynamics", 10, 10, 600, 600, "-wipe")
#prp -0 0 1000
#vup 0 1 0 
#display 2 4 100

# --------------------------
# End of recorder generation
# --------------------------


# ---------------------------------------
# Create and Perform the dynamic analysis
# ---------------------------------------

# Remove the static analysis & reset the time to 0.0
wipeAnalysis()
setTime(0.0)

# Now remove the loads and let the beam vibrate
remove("loadPattern", 1)

# Create the transient analysis
test("EnergyIncr", 1.0E-10, 20, 0)
algorithm("Newton")
numberer("RCM")
constraints("Plain") 
system("SparseGeneral", "-piv")
integrator("Newmark", 0.50, 0.25)
analysis("Transient")

# record once at time 0
record()

# Perform the transient analysis (20 sec)
analyze(100, 0.2)

wipe()
