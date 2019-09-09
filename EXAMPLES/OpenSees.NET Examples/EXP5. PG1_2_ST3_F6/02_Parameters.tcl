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

set	length 									40.
set width_inside_long						40.
set width_inside_trans_Rec					40.
set width_Perimeter_long_right				40.
set width_Perimeter_long_left				40.
set width_Perimeter_trans_up				40.
set width_Perimeter_trans_down				40.
set tickness								10.

set Tension									[expr 260000.]
set IP										[expr 260000.]
set OOP										[expr 260000.]
set about_Tension							[expr 260000.]
set about_IP								[expr 260000.]
set about_OOP								[expr 260000.]

set	E_Timber_Long							13200.0;
set	E_Timber_Long_mod						15200.0;
set	E_Timber_Trans							2200.0;
set	poisson_ratio							0.04;
set	G_Timber_Long							240.0;
set	G_Timber_Trans							220.0;




set Area_Perimeter_long_right				[expr $width_Perimeter_long_right*$tickness]
set Iz_Perimeter_long_right					[expr ($tickness*$width_Perimeter_long_right*$width_Perimeter_long_right*$width_Perimeter_long_right)/(12.0)];
set Iy_Perimeter_long_right					[expr 1.0*($width_Perimeter_long_right*$tickness*$tickness*$tickness)/(12.0)];
set Av_Perimeter_long_right					[expr (5*$width_Perimeter_long_right*$tickness)/(6)];
set Torsional_J_Perimeter_long_right		[expr ($width_Perimeter_long_right*$tickness*($width_Perimeter_long_right*$width_Perimeter_long_right+$tickness*$tickness))/12.0]
	
set Iy_Perimeter_long_right_mod_down		[expr 20.0*($width_Perimeter_long_right*$tickness*$tickness*$tickness)/(12.0)];
set Iy_Perimeter_long_right_mod_mid			[expr 3.0*($width_Perimeter_long_right*$tickness*$tickness*$tickness)/(12.0)];
set Iy_Perimeter_long_right_mod_up			[expr 0.1*($width_Perimeter_long_right*$tickness*$tickness*$tickness)/(12.0)];

set Iz_Perimeter_long_right_mod_down		[expr 48*100*($tickness*$width_Perimeter_long_right*$width_Perimeter_long_right*$width_Perimeter_long_right)/(12.0)];
set Iz_Perimeter_long_right_mod_mid			[expr 48*100*($tickness*$width_Perimeter_long_right*$width_Perimeter_long_right*$width_Perimeter_long_right)/(12.0)];
set Iz_Perimeter_long_right_mod_up			[expr 48*100*($tickness*$width_Perimeter_long_right*$width_Perimeter_long_right*$width_Perimeter_long_right)/(12.0)];


set Area_Perimeter_long_left				[expr $width_Perimeter_long_left*$tickness]
set Iz_Perimeter_long_left					[expr ($tickness*$width_Perimeter_long_left*$width_Perimeter_long_left*$width_Perimeter_long_left)/(12.0)];
set Iy_Perimeter_long_left					[expr ($width_Perimeter_long_left*$tickness*$tickness*$tickness)/(12.0)];
set Av_Perimeter_long_left					[expr (5*$width_Perimeter_long_left*$tickness)/(6)];
set Torsional_J_Perimeter_long_left			[expr ($width_Perimeter_long_left*$tickness*($width_Perimeter_long_left*$width_Perimeter_long_left+$tickness*$tickness))/12.0]

# set Iz_Perimeter_long_left_mod			[expr 0.1*($tickness*$width_Perimeter_long_left*$width_Perimeter_long_left*$width_Perimeter_long_left)/(12.0)];
set Iz_Perimeter_long_left_mod_down			[expr 0.25*0.01*($tickness*$width_Perimeter_long_left*$width_Perimeter_long_left*$width_Perimeter_long_left)/(12.0)];
set Iz_Perimeter_long_left_mod_mid			[expr 0.25*0.1*($tickness*$width_Perimeter_long_left*$width_Perimeter_long_left*$width_Perimeter_long_left)/(12.0)];
set Iz_Perimeter_long_left_mod_up			[expr 0.25*0.1*($tickness*$width_Perimeter_long_left*$width_Perimeter_long_left*$width_Perimeter_long_left)/(12.0)];

	
set Area_Perimeter_trans_up					[expr $width_Perimeter_trans_up*$tickness]
set Iz_Perimeter_trans_up					[expr ($tickness*$width_Perimeter_trans_up*$width_Perimeter_trans_up*$width_Perimeter_trans_up)/(12.0)];
set Iy_Perimeter_trans_up					[expr ($width_Perimeter_trans_up*$tickness*$tickness*$tickness)/(12.0)];
set Av_Perimeter_trans_up					[expr (5*$width_Perimeter_trans_up*$tickness)/(6)];
set Torsional_J_Perimeter_trans_up			[expr ($width_Perimeter_trans_up*$tickness*($width_Perimeter_trans_up*$width_Perimeter_trans_up+$tickness*$tickness))/12.0]
	
set Iy_Perimeter_trans_up_mod				[expr 10.0*($width_Perimeter_trans_up*$tickness*$tickness*$tickness)/(12.0)];
set Iy_Perimeter_trans_up_mod_right			[expr 0.5*($width_Perimeter_trans_up*$tickness*$tickness*$tickness)/(12.0)];
	
set Iz_Perimeter_trans_up_mod				[expr 2.*($tickness*$width_Perimeter_trans_up*$width_Perimeter_trans_up*$width_Perimeter_trans_up)/(12.0)];
	
	
set Area_Perimeter_trans_down				[expr $width_Perimeter_trans_down*$tickness]
set Iz_Perimeter_trans_down					[expr ($tickness*$width_Perimeter_trans_down*$width_Perimeter_trans_down*$width_Perimeter_trans_down)/(12.0)];
set Iy_Perimeter_trans_down					[expr ($width_Perimeter_trans_down*$tickness*$tickness*$tickness)/(12.0)];
set Av_Perimeter_trans_down					[expr (5*$width_Perimeter_trans_down*$tickness)/(6)];
set Torsional_J_Perimeter_trans_down		[expr ($width_Perimeter_trans_down*$tickness*($width_Perimeter_trans_down*$width_Perimeter_trans_down+$tickness*$tickness))/12.0]
	
set Iz_Perimeter_trans_down_mod				[expr 0.20*0.775*($tickness*$width_Perimeter_trans_down*$width_Perimeter_trans_down*$width_Perimeter_trans_down)/(12.0)];

	
set Area_inside_long						[expr $width_inside_long*$tickness]
set Iz_inside_long							[expr ($tickness*$width_inside_long*$width_inside_long*$width_inside_long)/(12.0)];
set Iy_inside_long							[expr 2.0*($width_inside_long*$tickness*$tickness*$tickness)/(12.0)];
set Av_inside_long							[expr (5*$width_inside_long*$tickness)/(6)];
set Torsional_J_inside_long					[expr ($width_inside_long*$tickness*($width_inside_long*$width_inside_long+$tickness*$tickness))/12.0]

set Iz_inside_long_mod						[expr 1.*($tickness*$width_inside_long*$width_inside_long*$width_inside_long)/(12.0)];



set Area_inside_trans_Rec					[expr $width_inside_trans_Rec*$tickness]
set Iz_inside_trans_Rec						[expr ($tickness*$width_inside_trans_Rec*$width_inside_trans_Rec*$width_inside_trans_Rec)/(12.0)];
set Iy_inside_trans_Rec						[expr 2.0*($width_inside_trans_Rec*$tickness*$tickness*$tickness)/(12.0)];
set Av_inside_trans_Rec						[expr (5*$width_inside_trans_Rec*$tickness)/(6)];
set Torsional_J_inside_trans_Rec			[expr ($width_inside_trans_Rec*$tickness*($width_inside_trans_Rec*$width_inside_trans_Rec+$tickness*$tickness))/12.0]

set Iz_inside_trans_Rec_mod					[expr 1.*($tickness*$width_inside_trans_Rec*$width_inside_trans_Rec*$width_inside_trans_Rec)/(12.0)];


puts "The Parameters are defined"