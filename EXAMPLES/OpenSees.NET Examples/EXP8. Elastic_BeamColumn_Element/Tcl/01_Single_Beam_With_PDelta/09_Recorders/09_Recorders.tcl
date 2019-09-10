# ************************************************************************************************************
#              ***************************************************************************
#                                  *************************************
#									Collaboration: Aryan, Petras, Iman                 
#                                  *************************************
#              ***************************************************************************
# ************************************************************************************************************





recorder Node -file	$folder/Displacements/Node_2_Disp.out	-node 2 	-dof	1	2	3	4	5	6	disp


# Define DISPLAY -------------------------------------------------------------
set  xPixels 1900;	# height of graphical window in pixels
set  yPixels 1100;	# height of graphical window in pixels
set  xLoc1 10;	# horizontal location of graphical window (0=upper left-most corner)
set  yLoc1 10;	# vertical location of graphical window (0=upper left-most corner)
set dAmp 50.;	# scaling factor for viewing deformed shape, it depends on the dimensions of the model
DisplayModel3D DeformedShape $dAmp $xLoc1 $yLoc1  $xPixels $yPixels


puts "The Recorders are defined"