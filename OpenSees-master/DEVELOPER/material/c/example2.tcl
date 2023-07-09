# OpenSees Example: Simple supported beam modeled with 2D solid elements

# Units: kips, in, sec
# ----------------------------

wipe; # clear opensees model


#-----------------------------
# Define the model
# ----------------------------

# Create ModelBuilder with 2 dimensions and 2 DOF/node
model BasicBuilder -ndm 2 -ndf 2

# create the material
nDMaterial elasticPlaneStressC 1 30000 0.25

# set type of quadrilateral element (uncomment one of the three #options)
set Quad quad
#set Quad bbarQuad
#set Quad enhancedQuad

# set up the arguments for the three considered elements
if {$Quad == "enhancedQuad" } {
	set eleArgs "PlaneStrain2D 1"
}
if {$Quad == "quad" } {
	set eleArgs "1 PlaneStrain2D 1"
}
if {$Quad == "bbarQuad" } {
	set eleArgs "1"
}

#define nodes
node 1 	0.0 		0.0   	
node 2 	100.0 	0.0 		
node 3 	100.0 	100.0 	
node 4 	0.0 		100.0 	

#define element
element quad 1 1 2 3 4 1.0 "PlaneStress" 1 

# define boundary conditions
fix   1 1 1
fix   4 1 1

# define the recorder
#---------------------
recorder Node -file RCUniaxialCompression.out -time -node 2 -dof disp

# define load pattern
#---------------------
pattern Plain 1 Linear {
load 2 1.0 0.0 
load 3 1.0 0.0 
}
# --------------------------------------------------------------------
# Start of static analysis (creation of the analysis & analysis itself)
# --------------------------------------------------------------------

# Load control with variable load steps
#                      init Jd min max
integrator LoadControl 1.0   1 1.0 100.0

# Convergence test
#              tolerance maxIter displayCode
test EnergyIncr 1.0e-12    10         0

# Solution algorithm
algorithm Newton

# DOF numberer
numberer RCM

# Cosntraint handler
constraints Plain

# System of equations solver
system ProfileSPD

# Type of analysis analysis
analysis Static

# Perform the analysis
analyze 10

# --------------------------
# End of static analysis
# --------------------------

print ele