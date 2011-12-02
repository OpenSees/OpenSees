#================================================================
# SetupIsotropicCompressionRecorder.tcl
# --  Set up recorders for isotropic compression analysis
# Zhaohui Yang UC Davis
# April 6, 2002
#================================================================

#Set up recorder

# Disp recorder for all nodes
recorder Node DisplW.out disp -time -node all -dof 1 2 3
#recorder plot Disp_iso.out "disp load"  10 10 300 300 -columns 2 1

# Stress recorder for all element
recorder Element all -time -file SigmaW.out stress
#recorder plot stress_iso.out "time Stress"  10 400 300 300 -columns 3 1

# Plstic info recorder for all element
recorder Element all -file PlastW.out plastic

# Gauss point coordinate info recorder for all element at selfweight stage
recorder Element all -file GaussPoint.out gausspoint

## create the display
#recorder display g3 10 10 200 200 -wipe
#prp 20 5.0 100.0
#vup 0  1 0
#viewWindow -30 30 -10 10
#display 10 0 5

