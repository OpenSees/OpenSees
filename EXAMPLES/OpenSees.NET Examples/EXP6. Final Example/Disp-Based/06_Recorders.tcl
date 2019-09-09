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

# ------------------------------ Recorders -------------------------------------------------------------------

recorder Node -file	$folder/Displacements/RF1.out	-node 412 	-dof	1	disp
recorder Node -file	$folder/Displacements/RF2.out	-node 404 	-dof	1	disp
recorder Node -file	$folder/Displacements/RF3.out	-node 217 	-dof	1	disp
recorder Node -file	$folder/Displacements/RF4.out	-node 210 	-dof	1	disp
recorder Node -file	$folder/Displacements/RF5.out	-node 204 	-dof	1	disp


# Define DISPLAY -------------------------------------------------------------
set  xPixels 1900;	# height of graphical window in pixels
set  yPixels 1100;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set dAmp 20.;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D DeformedShape $dAmp $xLoc1 $yLoc1  $xPixels $yPixels


puts "The Recorders are defined"