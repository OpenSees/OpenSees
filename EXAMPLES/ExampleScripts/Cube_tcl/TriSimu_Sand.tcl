#================================================================
# [jeremic,zhyang]@ucdavis.edu
# Tcl example file using Template3Dep material 
# 
# Triaxial Simulation with cubic block sample
#
# Two load stages:
#   stage 1: Applying Confining pressure
#   stage 2: Shearing
# Mar. 20, 2002 Zhaohui Yang UC Davis
#================================================================

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


#stage 1
#================================================================
# Isotropic compression analysis
#================================================================
source IsotropicCompression_Analysis.tcl


#stage 2

#===========================================================
# Apply friction(fixing the top and bottom plates with soil)
#================================================================
# You can test the differece between full frition or no friction 
# on the top and bottom platens by comment or de-comment the 
# following command
source Apply_FrictionBC.tcl


#================================================================
# Axial loading
#================================================================
source AxialLoad_Analysis.tcl

#================================================================
# Setup Link File for Visualizer Joey3D
#================================================================
source SetupLinkFileforJoey3D.tcl

wipe
