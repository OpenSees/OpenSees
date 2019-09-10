# ************************************************************************************************************
#              ***************************************************************************
#                                  *************************************
#				Mechanical Characterization of Integrally-Attached Timber Plate Structures
#                                   Aryan Rezaie Rad, EPFL, Suisse, 2018                  
#                                  *************************************
#              ***************************************************************************
# ************************************************************************************************************
#                           Section Name : In-Plane, Rectangular, Isotropic

# In-Plane Action of Forces: Timoshenko Formulation

wipe

model	basic	-ndm	3	-ndf	6;

source 01_Output_Folder.tcl
source 02_Parameters.tcl
source 021_DisplayModel3D.tcl
source 022_DisplayPlane.tcl
source 03_Materials.tcl
source 03_5_Section.tcl
source 04_Nodes.tcl
source 05_Elements.tcl
source 06_Recorders.tcl
# source 07_Load_Patterns_Displacement_Control.tcl
source 07_Load_Patterns_Load_Control.tcl
# source 08_Displacement_Control_Analysis.tcl
source 08_Load_Control_Analysis.tcl

puts "Analysis Done, Successfully"