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

# ------------------------------ Parameters -------------------------------------------------------------------
# Lenght = mm
# Force = Newton
# Mass = Tonne
# Time = Second

set	length 1400.0;
set width_Timoshenko_long	1100.0;
set width_Timoshenko_trans	1400.0;
set width_inside_long	385.0
set width_inside_trans	465.0
set width_Perimeter_long_right	220.0;
set width_Perimeter_long_left	110.0;
set width_Perimeter_trans_up	250.0;
set width_Perimeter_trans_down	140.0;
set tickness	40.0;

set Tension	666.0
set IP	8590.0
set OOP	7177.0
set about_Tension	1000000000000000.0
set about_IP	1000000000000000.0
set about_OOP	93000.0

set	E_Timber_Long	13200.0;
set	E_Timber_Trans	2200.0;
set	poisson_ratio	0.4;
set	G_Timber	820.0

set Area_Timoshenko_long			[expr $width_Timoshenko_long*$tickness]
set Iz_Timoshenko_long				[expr ($tickness*$width_Timoshenko_long*$width_Timoshenko_long*$width_Timoshenko_long)/(12.0)];
set Iy_Timoshenko_long				[expr ($width_Timoshenko_long*$tickness*$tickness*$tickness)/(12.0)];
set Av_Timoshenko_long				[expr (5*$width_Timoshenko_long*$tickness)/(6)];
set Torsional_J_Timoshenko_long		[expr ($width_Timoshenko_long*$tickness*($width_Timoshenko_long*$width_Timoshenko_long+$tickness*$tickness))/12.0]

set Area_Timoshenko_trans			[expr $width_Timoshenko_trans*$tickness]
set Iz_Timoshenko_trans				[expr ($tickness*$width_Timoshenko_trans*$width_Timoshenko_trans*$width_Timoshenko_trans)/(12.0)];
set Iy_Timoshenko_trans				[expr ($width_Timoshenko_trans*$tickness*$tickness*$tickness)/(12.0)];
set Av_Timoshenko_trans				[expr (5*$width_Timoshenko_trans*$tickness)/(6)];
set Torsional_J_Timoshenko_trans	[expr ($width_Timoshenko_trans*$tickness*($width_Timoshenko_trans*$width_Timoshenko_trans+$tickness*$tickness))/12.0]

set Area_Perimeter_long_right				[expr $width_Perimeter_long_right*$tickness]
set Iz_Perimeter_long_right					[expr 1.*($tickness*$width_Perimeter_long_right*$width_Perimeter_long_right*$width_Perimeter_long_right)/(12.0)];
set Iy_Perimeter_long_right					[expr 1.*($width_Perimeter_long_right*$tickness*$tickness*$tickness)/(12.0)];
set Av_Perimeter_long_right					[expr (5*$width_Perimeter_long_right*$tickness)/(6)];
set Torsional_J_Perimeter_long_right		[expr ($width_Perimeter_long_right*$tickness*($width_Perimeter_long_right*$width_Perimeter_long_right+$tickness*$tickness))/12.0]

set Area_Perimeter_long_left				[expr $width_Perimeter_long_left*$tickness]
set Iz_Perimeter_long_left				[expr 1.*($tickness*$width_Perimeter_long_left*$width_Perimeter_long_left*$width_Perimeter_long_left)/(12.0)];
set Iy_Perimeter_long_left				[expr 1.*($width_Perimeter_long_left*$tickness*$tickness*$tickness)/(12.0)];
set Av_Perimeter_long_left				[expr (5*$width_Perimeter_long_left*$tickness)/(6)];
set Torsional_J_Perimeter_long_left		[expr ($width_Perimeter_long_left*$tickness*($width_Perimeter_long_left*$width_Perimeter_long_left+$tickness*$tickness))/12.0]

set Area_Perimeter_trans_up			[expr $width_Perimeter_trans_up*$tickness]
set Iz_Perimeter_trans_up			[expr 1.0*($tickness*$width_Perimeter_trans_up*$width_Perimeter_trans_up*$width_Perimeter_trans_up)/(12.0)];
set Iy_Perimeter_trans_up			[expr 1.0*($width_Perimeter_trans_up*$tickness*$tickness*$tickness)/(12.0)];
set Av_Perimeter_trans_up			[expr (5*$width_Perimeter_trans_up*$tickness)/(6)];
set Torsional_J_Perimeter_trans_up	[expr ($width_Perimeter_trans_up*$tickness*($width_Perimeter_trans_up*$width_Perimeter_trans_up+$tickness*$tickness))/12.0]

set Area_Perimeter_trans_down			[expr $width_Perimeter_trans_down*$tickness]
set Iz_Perimeter_trans_down				[expr 1.0*($tickness*$width_Perimeter_trans_down*$width_Perimeter_trans_down*$width_Perimeter_trans_down)/(12.0)];
set Iy_Perimeter_trans_down				[expr 1.0*($width_Perimeter_trans_down*$tickness*$tickness*$tickness)/(12.0)];
set Av_Perimeter_trans_down				[expr (5*$width_Perimeter_trans_down*$tickness)/(6)];
set Torsional_J_Perimeter_trans_down	[expr ($width_Perimeter_trans_down*$tickness*($width_Perimeter_trans_down*$width_Perimeter_trans_down+$tickness*$tickness))/12.0]

set Area_inside_long				[expr $width_inside_long*$tickness]
set Iz_inside_long					[expr ($tickness*$width_inside_long*$width_inside_long*$width_inside_long)/(12.0)];
set Iy_inside_long					[expr ($width_inside_long*$tickness*$tickness*$tickness)/(12.0)];
set Av_inside_long					[expr (5*$width_inside_long*$tickness)/(6)];
set Torsional_J_inside_long			[expr ($width_inside_long*$tickness*($width_inside_long*$width_inside_long+$tickness*$tickness))/12.0]

set Area_inside_trans				[expr $width_inside_trans*$tickness]
set Iz_inside_trans					[expr ($tickness*$width_inside_trans*$width_inside_trans*$width_inside_trans)/(12.0)];
set Iy_inside_trans					[expr ($width_inside_trans*$tickness*$tickness*$tickness)/(12.0)];
set Av_inside_trans					[expr (5*$width_inside_trans*$tickness)/(6)];
set Torsional_J_inside_trans		[expr ($width_inside_trans*$tickness*($width_inside_trans*$width_inside_trans+$tickness*$tickness))/12.0]

puts "The Parameters are defined"