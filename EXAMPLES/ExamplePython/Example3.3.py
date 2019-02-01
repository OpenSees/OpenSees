# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# Portal Frame Example 3.3
# ------------------------
#  Reinforced concrete one-bay, one-story frame
#  Distributed vertical load on girder
#  Uniform excitation acting at fixed nodes in horizontal direction
#  
# 
# Example Objectives
# -----------------
#  Nonlinear dynamic analysis using Portal Frame Example 1 as staring point
#  Using Tcl Procedures 
#
# 
# Units: kips, in, sec
#
# Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
# Date: June 2017

# import the OpenSees Python module
import opensees as ops
import math

# ----------------------------------------------------
# Start of Model Generation & Initial Gravity Analysis
# ----------------------------------------------------

# remove existing model
ops.wipe()

# create ModelBuilder (with two-dimensions and 3 DOF/node)
ops.model("BasicBuilder", "-ndm",2, "-ndf",3)

# set default units
ops.defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")

# Create nodes
# ------------
# Set parameters for overall model geometry
width = 360.0
height = 144.0

# create nodes & add to Domain - command: node nodeId xCrd yCrd
ops.node(1, 0.0,   0.0)
ops.node(2, width, 0.0)
ops.node(3, 0.0,   height)
ops.node(4, width, height)

# set the boundary conditions - command: fix nodeID uxRestrnt? uyRestrnt? rzRestrnt?
ops.fix(1, 1, 1, 1)
ops.fix(2, 1, 1, 1)

# Define materials for nonlinear columns
# ------------------------------------------
# CONCRETE                        tag  f'c    ec0    f'cu   ecu
# Core concrete (confined)
ops.uniaxialMaterial("Concrete01", 1, -6.0, -0.004, -5.0, -0.014)
# Cover concrete (unconfined)
ops.uniaxialMaterial("Concrete01", 2, -5.0, -0.002, -0.0, -0.006)

# STEEL
# Reinforcing steel 
fy = 60.0;      # Yield stress
E = 30000.0;    # Young's modulus
#                              tag fy  E0  b
ops.uniaxialMaterial("Steel01", 3, fy, E, 0.01)

# Define cross-section for nonlinear columns
# ------------------------------------------
# set some parameters
colWidth = 15.0
colDepth = 24.0 
cover = 1.5
As = 0.60;     # area of no. 7 bars

# some variables derived from the parameters
y1 = colDepth/2.0
z1 = colWidth/2.0

ops.section("Fiber", 1)
# Create the concrete core fibers
ops.patch("rect", 1, 10, 1, cover-y1, cover-z1, y1-cover, z1-cover)
# Create the concrete cover fibers (top, bottom, left, right)
ops.patch("rect", 2, 10, 1, -y1, z1-cover, y1, z1)
ops.patch("rect", 2, 10, 1, -y1, -z1, y1, cover-z1)
ops.patch("rect", 2,  2, 1, -y1, cover-z1, cover-y1, z1-cover)
ops.patch("rect", 2,  2, 1,  y1-cover, cover-z1, y1, z1-cover)
# Create the reinforcing fibers (left, middle, right)
ops.layer("straight", 3, 3, As, y1-cover, z1-cover, y1-cover, cover-z1)
ops.layer("straight", 3, 2, As, 0.0, z1-cover, 0.0, cover-z1)
ops.layer("straight", 3, 3, As, cover-y1, z1-cover, cover-y1, cover-z1)
# define beam integration
np = 5;  # number of integration points along length of element
ops.beamIntegration("Lobatto", 1, 1, np)

# Define column elements
# ----------------------
# Geometric transformation for columns
#                       tag 
ops.geomTransf("PDelta", 1)

# Create the column elements
eleType = "forceBeamColumn"
#                   tag ndI ndJ transfTag integrationTag
ops.element(eleType, 1, 1, 3, 1, 1)
ops.element(eleType, 2, 2, 4, 1, 1)

# Define beam element
# ------------------
# Geometric transformation for beams
#                tag 
ops.geomTransf("Linear", 2)

# Create the beam element
#                               tag ndI ndJ  A     E       Iz   transfTag
ops.element("elasticBeamColumn", 3, 3, 4, 360.0, 4030.0, 8640.0, 2)

# Define gravity loads
# --------------------
# Set a parameter for the axial load
P = 180.0;                # 10% of axial capacity of columns

# create a Linear TimeSeries (load factor varies linearly with time) - command: timeSeries Linear $tag
ops.timeSeries("Linear", 1)

# create a Plain load pattern - command: pattern Plain $tag $timeSeriesTag { $loads }
ops.pattern("Plain", 1, 1, "-fact", 1.0)

# create the nodal load - command: load nodeID xForce yForce zMoment
ops.load(3, 0.0, -P, 0.0)
ops.load(4, 0.0, -P, 0.0)

# ------------------------------
# End of model generation
# ------------------------------


# ------------------------------
# Start of analysis generation
# ------------------------------

# create the system of equation
ops.system("BandGeneral")

# create the DOF numberer, the reverse Cuthill-McKee algorithm
ops.numberer("RCM")

# create the constraint handler, a Plain handler is used as homo constraints
ops.constraints("Transformation")

# Create the convergence test, the norm of the residual with a tolerance of 
# 1e-12 and a max number of iterations of 10
ops.test("NormDispIncr", 1.0E-12, 10, 3)

# create the solution algorithm, a Newton-Raphson algorithm
ops.algorithm("Newton")

# create the integration scheme, the LoadControl scheme using steps of 0.1
ops.integrator("LoadControl", 0.1)

# create the analysis object 
ops.analysis("Static")

# ------------------------------
# End of analysis generation
# ------------------------------


# ------------------------------
# Finally perform the analysis
# ------------------------------

# perform the gravity load analysis, requires 10 steps to reach the load level
ops.analyze(10)

print("Gravity load analysis completed\n")

# Set the gravity loads to be constant & reset the time in the domain
ops.loadConst("-time", 0.0)

# ----------------------------------------------------
# End of Model Generation & Initial Gravity Analysis
# ----------------------------------------------------


# ----------------------------------------------------
# Start of additional modelling for dynamic loads
# ----------------------------------------------------

# Define nodal mass in terms of axial load on columns
m = P/g

#       tag MX MY RZ
ops.mass(3, m, m, 0.0)
ops.mass(4, m, m, 0.0)

# Define dynamic loads
# --------------------
# Set some parameters
gmFile = "ARL360.g3"
dt = 0.02

# Set time series to be passed to uniform excitation
ops.timeSeries("Path", 2, "-filePath", gmFile, "-dt", dt, "-factor", g)

# Create UniformExcitation load pattern
#                               tag dir        tsTag
ops.pattern("UniformExcitation", 2, 1, "-accel", 2)

# Set the rayleigh damping factors for nodes & elements
ops.rayleigh(0.0, 0.0, 0.0, 0.000625)

# ----------------------------------------------------
# End of additional modelling for dynamic loads
# ----------------------------------------------------


# ---------------------------------------------------------
# Start of modifications to analysis for transient analysis
# ---------------------------------------------------------

# delete the old analysis and all its component objects
ops.wipeAnalysis()

# create the system of equation, a banded general storage scheme
ops.system("BandGeneral")

# create the DOF numberer, the reverse Cuthill-McKee algorithm
ops.numberer("RCM")

# create the constraint handler, a plain handler as homogeneous boundary
ops.constraints("Plain")

# create the convergence test, the norm of the residual with a tolerance of 
# 1e-12 and a max number of iterations of 10
ops.test("NormDispIncr", 1.0E-12, 10)

# create the solution algorithm, a Newton-Raphson algorithm
ops.algorithm("Newton")

# create the integration scheme, the Newmark with gamma=0.5 and beta=0.25
ops.integrator("Newmark", 0.5, 0.25) 

# create the analysis object 
ops.analysis("Transient")

# ---------------------------------------------------------
# End of modifications to analysis for transient analysis
# ---------------------------------------------------------


# ------------------------------
# Start of recorder generation
# ------------------------------

# Create a recorder to monitor nodal displacements
ops.recorder("Node", "-time", "-file", "disp.out", "-node", 3, 4, "-dof", 1, 2, 3, "disp")
ops.recorder("Node", "-time", "-file", "accel.out", "-node", 3, 4, "-dof", 1, 2, 3, "accel")
ops.recorder("Node", "-time", "-file", "totAccel.out", "-timeSeries", 2, 0, 0, "-node", 3, 4, "-dof", 1, 2, 3, "accel")

# Create recorders to monitor section forces and deformations
# at the base of the left column
ops.recorder("Element", "-time", "-file", "ele1secForce.out", "-ele", 1, "section", 1, "force")
ops.recorder("Element", "-time", "-file", "ele1secDef.out", "-ele", 1, "section", 1, "deformation")

# --------------------------------
# End of recorder generation
# ---------------------------------


# ------------------------------
# Finally perform the analysis
# ------------------------------

# record once at time 0
ops.record()

# Perform an eigenvalue analysis
lam = ops.eigen(2)
Tstart = 2.0*math.pi/math.sqrt(lam[0])
print("Fundamental period at start of transient analysis: ", Tstart, "sec\n")

# set some variables
tFinal = 2000 * 0.01
tCurrent = ops.getTime()
ok = 0

# Perform the transient analysis
while ((ok == 0) and (tCurrent < tFinal)):
    
    ok = ops.analyze(1, 0.01)
    
    # if the analysis fails try initial tangent iteration
    if (ok != 0):
        print("regular newton failed .. lets try an initial stiffness for this step")
        ops.test("NormDispIncr", 1.0E-12, 100, 0)
        ops.algorithm("ModifiedNewton", "-initial")
        ok = ops.analyze(1, 0.01)
        if (ok == 0):
            print("that worked .. back to regular newton")
        ops.test("NormDispIncr", 1.0E-12, 10) 
        ops.algorithm("Newton")
    
    tCurrent = ops.getTime()

# Print a message to indicate if analysis successful or not
if (ok == 0):
    print("\nTransient analysis completed SUCCESSFULLY\n")
else:
    print("\nTransient analysis FAILED\n")

# Perform an eigenvalue analysis
lam = ops.eigen(2)
Tend = 2.0*math.pi/math.sqrt(lam[0])
print("Fundamental period at end of transient analysis: ", Tend, "sec")

# Print state of node 3
ops.printModel("node", 3)
ops.wipe()
