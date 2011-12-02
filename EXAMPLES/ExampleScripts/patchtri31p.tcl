# ----------------------------
# Start of model generation
# ----------------------------
wipe 
file mkdir Output;

# Create ModelBuilder with 3 dimensions and 6 DOF/node
model basic -ndm 2 -ndf 2

# create the material
nDMaterial ElasticIsotropic   1   1e3   0.3  6.75 


# ################################ 
# build the model 
# ################################# 

node 1 0   0 
node 2 0.5 0 
node 3 1   0
node 4 0   0.5
node 5 0.4 0.6 
node 6 1   0.5
node 7 0   1
node 8 0.5 1
node 9 1   1

element tri31 1 1 2 4 1 PlaneStrain 1 
element tri31 2 2 5 4 1 PlaneStrain 1 
element tri31 3 2 3 5 1 PlaneStrain 1 
element tri31 4 3 6 5 1 PlaneStrain 1 
element tri31 5 4 5 7 1 PlaneStrain 1 
element tri31 6 5 8 7 1 PlaneStrain 1 
element tri31 7 5 6 8 1 PlaneStrain 1 
element tri31 8 6 9 8 1 PlaneStrain 1 

eval "recorder Node -file Output/disp.out -node 1 2 3 4 5 6 7 8 9 -dof 1 2 disp"

fix 1 1 1 
fix 3 0 1   

# point loads
pattern Plain 1 Constant {
  load 3  0.25 0.0 
  load 6  0.50 0.0
  load 9  0.25 0.0 
  load 4 -0.50 0.0
  load 7 -0.25 0.0  
}

# --------------------------------------------------------------------
# Start of static analysis (creation of the analysis & analysis itself)
# --------------------------------------------------------------------

# Load control with variable load steps
#                      init  Jd   min   max
#integrator LoadControl  0.01   1   1.0   10.0
integrator LoadControl  1   

# Convergence test
#                tolerance maxIter displayCode
test EnergyIncr  1.0e-6    100         5

# Solution algorithm
algorithm Newton

# DOF numberer
numberer RCM

# Cosntraint handler
constraints Plain 

# System of equations solver
system UmfPack

# Analysis for gravity load
analysis Static

# Perform the analysis
analyze   1 

# print ele
# print node

# --------------------------
# End of static analysis
# --------------------------

# ----------------------------
# Start of recorder generation
# ----------------------------

# recorder Node -file Node.out -time -node 4 -dof 2 disp
# recorder plot Node.out CenterNodeDisp 625 10 625 450 -columns 1 2

# create the display
recorder display g3 10 10 800 200 -wipe
prp 20 5.0 100.0
vup 0 1 0
viewWindow -2 2 -2 2
display 1 4 5


