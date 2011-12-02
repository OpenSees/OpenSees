#================================================================
#
# [jeremic,zhyang]@ucdavis.edu
# Dynamic analysis -- Plastic bowl loading
#
# Zhaohui Yang & Boris Jeremic UC Davis
# Nov. 1, 2002
#
#================================================================

wipe

#Directory to find all segments of the tcl files
set Dir DRMtcl;

#================================================================
# Create the modelbuilder
#================================================================
model BasicBuilder -ndm 3 -ndf 3


#================================================================
# Define some parameters and Units
#================================================================
source $Dir/Parameters.tcl


#================================================================
# build the model			      
#================================================================
source $Dir/InputNodes.tcl


#================================================================
#Elastic-plastic model
#================================================================
source $Dir/CreateMaterialModels.tcl
   

#================================================================
# Read element info
#================================================================
source $Dir/InputElements.tcl



#================================================================
# apply boundary conditions (fix bottom of the model)
#================================================================
source $Dir/Apply_BC.tcl


#================================================================
# transient analysis
#================================================================
source $Dir/PBowl_Analysis.tcl


#================================================================
# Setup Link File for Visualizer Joey3D
#================================================================
#source SetupLinkFileforJoey3D.tcl

wipe
