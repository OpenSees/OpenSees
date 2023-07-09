# [jeremic,zhyang]@ucdavis.edu
# static analysis -- Push-Over test on Single Pile in Sands
#
# Two load stages:
#   stage 1: Self Weight
#   stage 2: Push-Over at Pile Head
# Feb. 25, 2002

wipe

#================================================================
# Create the modelbuilder
#================================================================
model BasicBuilder -ndm 3 -ndf 3


#================================================================
# Define some parameters and Units
#================================================================
source Parameters.tcl


#================================================================
# build the model			      
#================================================================
source InputNodes.tcl


#================================================================
#Elastic-plastic model
#================================================================
source CreateMaterialModels.tcl
   

#================================================================
# Read element info
#================================================================
source InputElements.tcl


##================================================================
## Set up recorders for Self Weight analysis
##================================================================

# Disp recorder for a top node
recorder Node Disp_iso.out disp -time -node 807 810 -dof 1 2 3
#recorder plot Disp_iso.out "disp load"  10 10 300 300 -columns 2 1


# Stress recorder for all element close to bottom
recorder Element all -time -file stress_iso.out stress
#recorder plot stress_iso.out "time Stress"  10 400 300 300 -columns 3 1

# Plstic info recorder for all element
recorder Element all -file Plastic_iso.out plastic


# Gauss point coordinate info recorder for all element at selfweight stage
recorder Element all -file GaussPoint.out gausspoint

#================================================================
# stage 1 -- Self Weight analysis
#================================================================
source SelfWeight_Analysis.tcl



#================================================================
# stage 2 -- push-over analysis
#================================================================
#Clean selfweight analysis method
wipeAnalysis

#set previous load constant
# The following command is necessary in order to avoid seg fault
setTime 0.0
loadConst
#remove loadPattern 1

#================================================================
# Plug in pile elements
#================================================================
#source Plug_In_PileElements.tcl


source Apply_PushoverBC.tcl
source Apply_LateralLoad.tcl

#================================================================
# Set up recorders for Self Weight analysis
#================================================================
remove recorders

# Displacement recorder for a top node
recorder Node Disp_L.out disp -time -node 807 -dof 1 2 3
#recorder plot Disp_L.out "disp load"  10 10 300 300 -columns 2 1

# Plastic info recorder
recorder Element 528  -file Plastic_L.out plastic


#================================================================
# push-over analysis
#================================================================
source Pushover_Analysis.tcl

wipe
