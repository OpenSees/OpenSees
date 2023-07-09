# ----------------------------
# Start of model generation
# ----------------------------

model basic -ndm 2 -ndf 2

nDMaterial ElasticIsotropic   1   1000   0.25  6.75 

# define some  parameters

set Quad  quad
set Quad  bbarQuad
set Quad  enhancedQuad

if {$Quad == "enhancedQuad" } {
    set eleArgs "PlaneStrain2D  1"
} 

if {$Quad == "quad" } {
    set eleArgs "1 PlaneStrain2D  1"
} 

if {$Quad == "bbarQuad" } {
    set eleArgs "1"
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

recorder Node -xml Node.out -time -node $l1 -dof 2 disp

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

# Perform the transient analysis (20 sec)
#       numSteps  dt
analyze  500     0.5

