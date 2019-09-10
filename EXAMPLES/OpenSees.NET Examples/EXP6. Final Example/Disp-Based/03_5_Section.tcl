# ************************************************************************************************************
#              ***************************************************************************
#                                  *************************************
#				Mechanical Characterization of Integrally-Attached Timber Plate Structures
#                                   Aryan Rezaie Rad, EPFL, Suisse, 2018                  
#                                  *************************************
#              ***************************************************************************
# ************************************************************************************************************
#                           Section Name : Out-of-Plane, Rectangular, Isotropic


# ------------------------------ Sections -------------------------------------------------------------------




set width_Perimeter_long_right_tag	1

section Fiber $width_Perimeter_long_right_tag {

patch quad $timbermattag_long	 20		5		[expr -$width_Perimeter_long_right/2.0]	 	[expr $tickness/2.0]		 [expr $width_Perimeter_long_right/2.0]	 [expr $tickness/2.0]	 [expr $width_Perimeter_long_right/2.0]	 [expr -$tickness/2.0] 	[expr -$width_Perimeter_long_right/2.0]		 [expr -$tickness/2.0]

}



set width_Perimeter_long_left_tag	2

section Fiber $width_Perimeter_long_left_tag {

patch quad $G_Timber_Tag_Long	 20		5		[expr -$width_Perimeter_long_left/2.0]	 	[expr $tickness/2.0]		 [expr $width_Perimeter_long_left/2.0]	 [expr $tickness/2.0]	 [expr $width_Perimeter_long_left/2.0]	 [expr -$tickness/2.0] 	[expr -$width_Perimeter_long_left/2.0]		 [expr -$tickness/2.0]

}


set width_Perimeter_trans_up_tag	3

section Fiber $width_Perimeter_trans_up_tag {

patch quad $timbermattag_trans	 5		20		[expr -$width_Perimeter_trans_up/2.0]	[expr -$tickness/2.0]		[expr -$width_Perimeter_trans_up/2.0]	[expr $tickness/2.0]		[expr $width_Perimeter_trans_up/2.0]		[expr $tickness/2.0]		[expr $width_Perimeter_trans_up/2.0]	[expr -$tickness/2.0]

}


set width_Perimeter_trans_down_tag	4

section Fiber $width_Perimeter_trans_down_tag {

patch quad $G_Timber_Tag_Trans	 5		20		[expr -$width_Perimeter_trans_down/2.0]	[expr -$tickness/2.0]		[expr -$width_Perimeter_trans_down/2.0]	[expr $tickness/2.0]		[expr $width_Perimeter_trans_down/2.0]		[expr $tickness/2.0]		[expr $width_Perimeter_trans_down/2.0]	[expr -$tickness/2.0]

}


set width_inside_long_tag	5

section Fiber $width_inside_long_tag {

patch quad $Long_Shear_Timoshenko_Parallel	 5		20		[expr -$width_inside_long/2.0]	 	[expr $tickness/2.0]		 [expr $width_inside_long/2.0]	 [expr $tickness/2.0]	 [expr $width_inside_long/2.0]	 [expr -$tickness/2.0] 	[expr -$width_inside_long/2.0]		 [expr -$tickness/2.0]

}



set width_inside_trans_tag	6

section Fiber $width_inside_trans_tag {

patch quad $Long_Shear_Timoshenko_Serial	 5		20		[expr -$width_inside_trans_Rec/2.0]	[expr -$tickness/2.0]		[expr -$width_inside_trans_Rec/2.0]	[expr $tickness/2.0]		[expr $width_inside_trans_Rec/2.0]		[expr $tickness/2.0]		[expr $width_inside_trans_Rec/2.0]	[expr -$tickness/2.0]

}






puts "The Sections are defined"