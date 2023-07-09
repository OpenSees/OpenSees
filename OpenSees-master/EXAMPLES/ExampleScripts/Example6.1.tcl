# ----------------------------
# Start of model generation
# ----------------------------

# Create ModelBuilder with 3 dimensions and 6 DOF/node
model basic -ndm 2 -ndf 2

# create the material
nDMaterial ElasticIsotropic   1   1000   0.25  [expr 6.75/(12.0*32.174)]

# Define geometry
# ---------------

# define some  parameters
set Quad  quad
#set Quad SSPquad
#set Quad  bbarQuad
#set Quad  enhancedQuad

set thick 2.0
if {$Quad == "quad" } {
    set eleArgs "$thick PlaneStrain 1"
} 
if {$Quad == "SSPquad" } {
    set eleArgs "1 PlaneStrain $thick"
} 
if {$Quad == "bbarQuad" } {
    set eleArgs "$thick 1"
}
if {$Quad == "enhancedQuad" } {
    set eleArgs "$thick PlaneStrain 1"
} 

set nx 10; # NOTE: nx MUST BE EVEN FOR THIS EXAMPLE
set ny 4
set bn [expr $nx + 1 ] 
set l1 [expr $nx/2 + 1 ] 
set l2 [expr $l1 + $ny*($nx+1) ]

# now create the nodes and elements using the block2D command
block2D $nx $ny   1 1  $Quad  $eleArgs {
    1   0   0
    2  40   0
    3  40  10
    4   0  10
}

# Single point constraints
#   node   u1  u2    
fix    1    1    1   
fix  $bn    0    1   

# Gravity loads
pattern Plain 1 Linear {
    load $l1  0.0  -1.0
    load $l2  0.0  -1.0
}

# ----------------------------
# End of model generation
# ----------------------------


# --------------------------------------------------------------------
# Start of static analysis (creation of the analysis & analysis itself)
# --------------------------------------------------------------------

# Load control with variable load steps
#                       init   Jd  min   max
integrator LoadControl  1.0  1   1.0   10.0

# Convergence test
#                  tolerance maxIter displayCode
test EnergyIncr  1.0e-12    10         0

# Solution algorithm
algorithm Newton

# DOF numberer
numberer RCM

# Constraint handler
constraints Plain 

# System of equations solver
system ProfileSPD

# Analysis for gravity load
analysis Static

# Perform the analysis
analyze   10     

# --------------------------
# End of static analysis
# --------------------------


# ----------------------------
# Start of recorder generation
# ----------------------------

recorder Node -file Node.out -time -node $l1 -dof 2 disp
recorder plot Node.out CenterNodeDisp 10 210 800 400 -columns 1 2

# create the display
recorder display OpenSees 10 10 800 200 -wipe
prp 20 5.0 100.0
vup 0 1 0
viewWindow -30 30 -10 10
display 1 4 5

# --------------------------
# End of recorder generation
# --------------------------


# ---------------------------------------
# Create and Perform the dynamic analysis
# ---------------------------------------

# Remove the static analysis & reset the time to 0.0
wipeAnalysis
setTime 0.0

# Now remove the loads and let the beam vibrate
remove loadPattern 1

# Create the transient analysis
test EnergyIncr  1.0e-12    10         0
algorithm Newton
numberer RCM
constraints Plain 
integrator Newmark 0.5 0.25
system ProfileSPD
#integrator GeneralizedMidpoint 0.50
analysis Transient

record

# Perform the transient analysis (20 sec)
#       numSteps  dt
analyze  1000     0.05

print node $l1 $l2
