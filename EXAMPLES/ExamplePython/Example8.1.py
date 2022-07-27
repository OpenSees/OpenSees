# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# Cantilever Beam Example 8.1
# ---------------------------
#  Cantilever beam modeled with
#  three dimensional brick elements
# 
# Example Objectives
# ------------------
#  test different brick elements
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

# create ModelBuilder (with three-dimensions and 3 DOF/node)
model("BasicBuilder", "-ndm",3, "-ndf",3)

# set default units
defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")

# Define the material
# -------------------
#                               matTag  E     nu   rho
nDMaterial("ElasticIsotropic", 1, 100.0, 0.25, 1.27) 

# Define geometry
# ---------------
Brick = "stdBrick"
#Brick = "bbarBrick"
#Brick = "SSPbrick"

nz = 6
nx = 2 
ny = 2

nn = int((nz+1)*(nx+1)*(ny+1))

# mesh generation
#          numX numY numZ startNode startEle eleType eleArgs? coords?
block3D(nx, ny, nz, 1, 1,
            Brick, 1,
            1, -1.0, -1.0,  0.0,
            2,  1.0, -1.0,  0.0,
            3,  1.0,  1.0,  0.0,
            4, -1.0,  1.0,  0.0, 
            5, -1.0, -1.0, 10.0,
            6,  1.0, -1.0, 10.0,
            7,  1.0,  1.0, 10.0,
            8, -1.0,  1.0, 10.0)

# boundary conditions
fixZ(0.0, 1, 1, 1) 

# Define point load
# create a Linear time series
timeSeries("Linear", 1)
# create a Plain load pattern
p = 0.10
pattern("Plain", 1, 1, "-fact", 1.0)
load(nn, p, p, 0.0)

# print model
#printModel()
printModel("-JSON", "-file", "Example8.1.json")

# ----------------------- 
# End of model generation
# -----------------------


# ------------------------
# Start of static analysis
# ------------------------

# Load control with variable load steps
#                            init  Jd  min   max
integrator("LoadControl", 1.0, 1) 

# Convergence test
#                     tolerance maxIter displayCode
test("NormUnbalance", 1.0E-10, 20, 0)

# Solution algorithm
algorithm("Newton")

# DOF numberer
numberer("RCM")

# Constraint handler
constraints("Plain")

# System of equations solver
system("ProfileSPD")

# Analysis for gravity load
analysis("Static")

# Perform the analysis
analyze(5)

# --------------------------
# End of static analysis
# --------------------------


# ----------------------------
# Start of recorder generation
# ----------------------------

recorder("Node", "-file", "Node.out", "-time", "-node", nn, "-dof", 1, "disp")
recorder("Element", "-file", "Elem.out", "-time", "-eleRange", 1, 10, "material", "1", "strains")
#recorder("plot", "Node.out", "CenterNodeDisp", 625, 10, 625, 450, "-columns", 1, 2)

# create the display
#recorder("display", "VibratingBeam", 100, 40, 500, 500, "-wipe")
#prp -100 100 120.5
#vup 0 1 0 
#display 1 4 1 

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

# add some mass proportional damping
rayleigh(0.01, 0.0, 0.0, 0.0)

# Create the transient analysis
test("EnergyIncr", 1.0E-10, 20, 0)
algorithm("Newton")
numberer("RCM")
constraints("Plain")
system("ProfileSPD")
integrator("Newmark", 0.5, 0.25)
analysis("Transient")

# record once at time 0
record()

# Perform the transient analysis (20 sec)
#         numSteps dt
analyze(1000, 1.0)

wipe()
