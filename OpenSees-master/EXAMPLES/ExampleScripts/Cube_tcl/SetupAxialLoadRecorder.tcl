#================================================================
# SetupAxialLoadRecorder.tcl
# --  Set up recorders for axial loading analysis
# Zhaohui Yang UC Davis
# April 6, 2002
#================================================================

# Stop the old recorders by destroying them
remove recorders

#Set up recorder for axial loading stage

# Displacement recorder for a top node
recorder Node DisplL.out disp -time -node all -dof 1 2 3
#recorder plot Disp_L.out "disp load"  10 10 300 300 -columns 2 1

# Plastic info recorder
recorder Element all  -file PlastL.out plastic

# Stress recorder
recorder Element all -file  Stress_shear.out stress

