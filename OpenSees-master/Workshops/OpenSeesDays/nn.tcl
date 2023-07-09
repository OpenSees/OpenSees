# OpenSees Example: Simple supported beam modeled with 2D solid elements

# Units: kips, in, sec
# ----------------------------

wipe; # clear opensees model

# create data directory
file mkdir Data

#-----------------------------
# Define the model
# ----------------------------

# Create ModelBuilder with 2 dimensions and 2 DOF/node
model BasicBuilder -ndm 2 -ndf 2

# create the material
nDMaterial ElasticIsotropic 1 1000 0.25 3.0

# set type of quadrilateral element (uncomment one of the three options)
set Quad quad
#set Quad bbarQuad
#set Quad enhancedQuad

# set up the arguments for the three considered elements
if {$Quad == "enhancedQuad" } {
	set eleArgs "PlaneStress2D 1"
}
if {$Quad == "quad" } {
	set eleArgs "1 PlaneStress2D 1"
}
if {$Quad == "bbarQuad" } {
	set eleArgs "1"
}

# set up the number of elements in x (nx) and y (ny) direction
set nx 2; # NOTE: nx MUST BE EVEN FOR THIS EXAMPLE
set ny 1

# define numbering of node at the left support (bn), and the two nodes at load application (l1, l2)
set bn [expr $nx + 1]
set l1 [expr $nx/2 + 1]
set l2 [expr $l1 + $ny*($nx+1)]

# create the nodes and elements using the block2D command
block2D $nx $ny 1 1 $Quad $eleArgs {
	1   0   0
	2  40   0
	3  40  10
	4   0  10
}

# define boundary conditions
fix   1 1 1
fix $bn 0 1

# define the recorder
#---------------------
recorder Node -file Data/Node.out -time -node $l1 -dof 2 disp

# define load pattern
#---------------------
pattern Plain 1 Linear {
load $l1 0.0 -1.0
load $l2 0.0 -1.0
}
# --------------------------------------------------------------------
# Start of static analysis (creation of the analysis & analysis itself)
# --------------------------------------------------------------------

# Load control with variable load steps
#                      init Jd min max
integrator LoadControl 1.0   1 1.0 10.0

# Convergence test
#              tolerance maxIter displayCode
test EnergyIncr 1.0e-12    10         0

# Solution algorithm
algorithm Newton

# DOF numberer
numberer RCM

# Constraint handler
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

# -------------------------------------
# create display for transient analysis
#--------------------------------------
#                    $windowTitle       $xLoc $yLoc $xPixels $yPixels
recorder display "Simply Supported Beam" 10     10      800     200    -wipe  
prp 20 5.0 1.0;                                      # projection reference point (prp); defines the center of projection (viewer eye)
vup  0  1 0;                                         # view-up vector (vup) 
vpn  0  0 1;                                         # view-plane normal (vpn)     
viewWindow -30 30 -10 10;                            # coordinates of the window relative to prp  
display 10 0 5;                                      # the 1st arg. is the tag for display mode
                                                     # the 2nd arg. is magnification factor for nodes, the 3rd arg. is magnif. factor of deformed shape

# ---------------------------------------
# Create and Perform the dynamic analysis
# ---------------------------------------
#define damping
rayleigh 0. 0. 0. [expr 2*0.02/sqrt([eigen 1])];

# Remove the static analysis & reset the time to 0.0
wipeAnalysis
setTime 0.0

# Now remove the loads and let the beam vibrate
remove loadPattern 1

# Create the transient analysis
test EnergyIncr 1.0e-12 10 0
algorithm Newton
numberer RCM
constraints Plain
integrator Newmark 0.5 0.25
system BandGeneral
analysis Transient

# Perform the transient analysis (50 sec)
#analyze 1500 0.5

