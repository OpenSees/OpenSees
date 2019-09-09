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

recorder Node -file	$folder/Displacements/Disp_21.out	-node 21 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_23.out	-node 23 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_33.out	-node 33 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_13.out	-node 13 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_8.out	-node 8 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_7.out	-node 7 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_3.out	-node 3 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_10.out	-node 10 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_19.out	-node 19 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_18.out	-node 18 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_17.out	-node 17 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_20.out	-node 20 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_22.out	-node 22 -dof	1	2	3	4	5	6	disp
recorder Node -file	$folder/Displacements/Disp_28.out	-node 28 -dof	1	2	3	4	5	6	disp

recorder Node -file	$folder/Nodal_Reactions/Nodal_Reactions_27.out	-node 27 -dof	1	2	3	4	5	6	reaction
recorder Node -file	$folder/Nodal_Reactions/Nodal_Reactions_29.out	-node 29 -dof	1	2	3	4	5	6	reaction


recorder Element -file	$folder/Element_Forces/Element_Forces_1.out	-ele	1	force
recorder Element -file	$folder/Element_Forces/Element_Forces_2.out	-ele	2	force


# Define DISPLAY -------------------------------------------------------------
set  xPixels 1900;	# height of graphical window in pixels
set  yPixels 1100;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set dAmp 20.;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D DeformedShape $dAmp $xLoc1 $yLoc1  $xPixels $yPixels


puts "The Recorders are defined"