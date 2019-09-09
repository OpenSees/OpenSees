# ************************************************************************************************************
#              ***************************************************************************
#                                  *************************************
#									Collaboration: Aryan, Petras, Iman                 
#                                  *************************************
#              ***************************************************************************
# ************************************************************************************************************





set 	Node_1_x			0.
set 	Node_1_y			0.
set 	Node_1_z			0.

set 	Node_2_x			100.0
set 	Node_2_y			-250.0
set 	Node_2_z			600.0


set		Vec_X_x				[expr	$Node_2_x-$Node_1_x]
set		Vec_X_y				[expr	$Node_2_y-$Node_1_y]	
set		Vec_X_z				[expr	$Node_2_z-$Node_1_z]

set		Vec_Y_x				[expr	-$Vec_X_y]
set		Vec_Y_y				[expr	$Vec_X_x]
set		Vec_Y_z				0.0

set		Vec_Z_x				[expr	$Vec_X_y*$Vec_Y_z-$Vec_X_z*$Vec_Y_y]
set		Vec_Z_y				[expr	$Vec_X_z*$Vec_Y_x-$Vec_X_x*$Vec_Y_z]
set		Vec_Z_z				[expr	$Vec_X_x*$Vec_Y_y-$Vec_X_y*$Vec_Y_x]


puts	$Vec_Z_x
puts	$Vec_Z_y
puts	$Vec_Z_z



set Trans_tag	2
geomTransf PDelta  $Trans_tag	$Vec_Z_x	$Vec_Z_y	$Vec_Z_z