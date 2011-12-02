#================================================================
# SetupPBowlAnalysisRecorder.tcl
# --Set up recorders for plastic bowl analysis
# Zhaohui Yang & Boris Jeremic UC Davis
# Nov. 1, 2002
#================================================================

remove recorders

# Displacement recorder for a top node
recorder Node $Dir/DisplL.out disp -time -node 109 -dof 1 2 3
recorder plot $Dir/DisplL.out "disp load"  10 10 300 300 -columns 2 1

# Plastic info recorder
#recorder Element all  -file PlastL.out plastic
