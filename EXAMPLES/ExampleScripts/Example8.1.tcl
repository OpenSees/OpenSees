# ----------------------------
# Start of model generation
# ----------------------------

# Create ModelBuilder with 3 dimensions and 6 DOF/node
model basic -ndm 3 -ndf 3


# create the material
nDMaterial ElasticIsotropic   1   100   0.25  1.27

# Define geometry
# ---------------

# define some  parameters
set eleArgs "1" 

#set element stdBrick
set element bbarBrick

set nz 6
set nx 2 
set ny 2

set nn [expr ($nz+1)*($nx+1)*($ny+1) ]

# mesh generation
block3D $nx $ny $nz   1 1  $element  $eleArgs {
    1   -1     -1      0
    2    1     -1      0
    3    1      1      0
    4   -1      1      0 
    5   -1     -1     10
    6    1     -1     10
    7    1      1     10
    8   -1      1     10
}


set load 0.10

# Constant point load
pattern Plain 1 Linear {
   load $nn  $load  $load  0.0
}

# boundary conditions
fixZ 0.0   1 1 1 

# --------------------------------------------------------------------
# Start of static analysis (creation of the analysis & analysis itself)
# --------------------------------------------------------------------

# Load control with variable load steps
#                       init   Jd  min   max
integrator LoadControl  1.0  1 

# Convergence test
#                  tolerance maxIter displayCode
test NormUnbalance     1.0e-10    20     0

# Solution algorithm
algorithm Newton

# DOF numberer
numberer RCM

# Cosntraint handler
constraints Plain 

# System of equations solver
system ProfileSPD

# Analysis for gravity load
analysis Static 

# Perform the analysis
analyze 5


# --------------------------
# End of static analysis
# --------------------------

# ----------------------------
# Start of recorder generation
# ----------------------------

recorder Node -file Node.out -time -node $nn -dof 1 disp
#recorder plot Node.out CenterNodeDisp 625 10 625 450 -columns 1 2

recorder display ShakingBeam 100 40 500 500 -wipe
prp -100 100 120.5
vup 0 1 0 
display 1 4 1 

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

# add some mass proportional damping
rayleigh 0.01 0.0 0.0 0.0

# Create the transient analysis
test EnergyIncr     1.0e-10    20   0
algorithm Newton
numberer RCM
constraints Plain 
integrator Newmark 0.5 0.25
#integrator GeneralizedMidpoint 0.50
analysis Transient


# Perform the transient analysis (20 sec)
#       numSteps  dt
#analyze 1000 1.0
puts [eigen -standard 1]
puts [eigen -standard -fullGenLapack 1]




