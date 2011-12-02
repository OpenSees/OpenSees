# ----------------------------
# Start of model generation
# ----------------------------
model basic -ndm 3 -ndf 6

# create the material
section ElasticMembranePlateSection  1   3.0e3  0.25  1.175  1.27

# set some parameters for node and element generation
set Plate shell

set eleArgs "1"

#these should both be even
set nx 8
set ny 2

#loaded nodes
set mid [expr (  ($nx+1)*($ny+1)+1 ) / 2 ]
set side1 [expr ($nx + 2)/2 ] 
set side2 [expr ($nx+1)*($ny+1) - $side1 + 1 ]

# generate the nodes and elements
block2D $nx $ny 1 1 $Plate $eleArgs {
    1   -20    0     0
    2   -20    0    40
    3    20    0    40
    4    20    0     0
    5   -10   10    20 
    7    10   10    20   
    9    0    10    20 
} 

# add some loads
pattern Plain 1 Linear {
    load $mid    0.0  -0.5   0.0   0.0   0.0  0.0
    load $side1  0.0  -0.25  0.0   0.0   0.0  0.0
    load $side2  0.0  -0.25  0.0   0.0   0.0  0.0
}


# define the boundary conditions
# rotation free about x-axis (remember right-hand-rule)
fixZ 0.0   1 1 1  0 1 1
fixZ 40.0  1 1 1  0 1 1   

# Load control with variable load steps
#                       init   Jd  min   max
integrator LoadControl  1.0  1   1.0   10.0

# Convergence test
#                  tolerance maxIter displayCode
test EnergyIncr     1.0e-10    20       0

# Solution algorithm
algorithm Newton

# DOF numberer
numberer RCM

# Cosntraint handler
constraints Plain 

# System of equations solver
system SparseGeneral -piv
#system ProfileSPD

# Analysis for gravity load
#analysis Transient 
analysis Static 


# Perform the gravity load analysis
analyze 5

# --------------------------
# End of static analysis
# --------------------------

# ----------------------------
# Start of recorder generation
# ----------------------------

recorder Node -file Node.out -time -node $mid -dof 2 disp
recorder plot Node.out CenterNodeDisp 625 10 625 450 -columns 1 2

recorder display shellDynamics 10 10 600 600 -wipe
prp -100 20 30
vup 0 1 0 
display 2 4 100

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
test EnergyIncr     1.0e-10    20    0
algorithm Newton
numberer RCM
constraints Plain 
system SparseGeneral -piv
#integrator GeneralizedMidpoint 0.50
integrator Newmark 0.50 0.25
analysis Transient

# Perform the transient analysis (20 sec)
analyze 100 0.2
