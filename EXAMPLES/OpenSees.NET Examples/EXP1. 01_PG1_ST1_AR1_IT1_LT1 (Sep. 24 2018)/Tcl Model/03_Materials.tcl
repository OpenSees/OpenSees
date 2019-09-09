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

# ------------------------------ Materials -------------------------------------------------------------------


set rigidmattag	1
set semi_rigidmattag	101

set Freemattag	2

set timbermattag_long	30
set timbermattag_trans	31

set Axial_perimeter_long_right	40
set Axial_perimeter_long_left	41
set Axial_perimeter_trans_up	42
set Axial_perimeter_trans_down	43
set Axial_inside_long	44
set Axial_inside_trans	45

set Shear_1_perimeter_long_right	50
set Shear_1_perimeter_long_left	51
set Shear_1_perimeter_trans_up	52
set Shear_1_perimeter_trans_down	53
set Shear_1_inside_long	54
set Shear_1_inside_trans	55

set Shear_2_perimeter_long	60
set Shear_2_perimeter_trans	61
set Shear_2_inside_long	62
set Shear_2_inside_trans	63

set Trans_Shear_Timoshenko_Eq	70
set Long_Shear_Timoshenko_Eq	71

set Stiffness_Tension	10
set Stiffness_IP	11
set Stiffness_OOP	12
set Stiffness_about_Tension	13
set Stiffness_about_IP	14
set Stiffness_about_OOP	15


uniaxialMaterial	Elastic	$rigidmattag  				[expr 1.0e+15]
uniaxialMaterial	Elastic	$semi_rigidmattag			[expr 1.0e+4]

uniaxialMaterial	Elastic	$Freemattag  				[expr 1.0e-1]

uniaxialMaterial	Elastic	$timbermattag_long			$E_Timber_Long
uniaxialMaterial	Elastic	$timbermattag_trans			$E_Timber_Trans

uniaxialMaterial	Elastic	$Axial_perimeter_long_right			[expr 100.0*$E_Timber_Long*$Area_Perimeter_long_right/$length]; 		# for LT2
uniaxialMaterial	Elastic	$Axial_perimeter_long_left			[expr 100.0*$E_Timber_Long*$Area_Perimeter_long_left/$length]; 		# for LT2
uniaxialMaterial	Elastic	$Axial_perimeter_trans_up			[expr 100.0*$E_Timber_Trans*$Area_Perimeter_trans_up/$length]; 		# for LT2
uniaxialMaterial	Elastic	$Axial_perimeter_trans_down			[expr 100.0*$E_Timber_Trans*$Area_Perimeter_trans_down/$length]; 		# for LT2
uniaxialMaterial	Elastic	$Axial_inside_long					[expr 1.*$E_Timber_Long*$Area_inside_long/$length]; 			# for LT2
uniaxialMaterial	Elastic	$Axial_inside_trans					[expr 1.*$E_Timber_Trans*$Area_inside_trans/$length]; 			# for LT2

uniaxialMaterial	Elastic	$Shear_1_perimeter_long_right		[expr 100.0*(5*$width_Perimeter_long_right*$tickness)/(6)];														
uniaxialMaterial	Elastic	$Shear_1_perimeter_long_left		[expr 100.0*(5*$width_Perimeter_long_left*$tickness)/(6)];														
uniaxialMaterial	Elastic	$Shear_1_perimeter_trans_up			[expr 100.0*(5*$width_Perimeter_trans_up*$tickness)/(6)];												 		
uniaxialMaterial	Elastic	$Shear_1_perimeter_trans_down		[expr 100.0*(5*$width_Perimeter_trans_down*$tickness)/(6)];												 		
uniaxialMaterial	Elastic	$Shear_1_inside_long				[expr 1.0*(5*$width_inside_long*$tickness)/(6)];			
uniaxialMaterial	Elastic	$Shear_1_inside_trans				[expr 1.0*(5*$width_inside_trans*$tickness)/(6)]; 			

uniaxialMaterial	Elastic	$Shear_2_perimeter_long		1.0e+15;	
uniaxialMaterial	Elastic	$Shear_2_perimeter_trans	1.0e+15;	
uniaxialMaterial	Elastic	$Shear_2_inside_long		1.0e+15;	
uniaxialMaterial	Elastic	$Shear_2_inside_trans		1.0e+15;	

uniaxialMaterial	Elastic	$Trans_Shear_Timoshenko_Eq	[expr 34787.878]
uniaxialMaterial	Elastic	$Long_Shear_Timoshenko_Eq	[expr 16081.971]

uniaxialMaterial	Elastic	$Stiffness_Tension			[expr	$Tension]
uniaxialMaterial	Elastic	$Stiffness_IP				[expr	$IP]
uniaxialMaterial	Elastic	$Stiffness_OOP				[expr	$OOP]
uniaxialMaterial	Elastic	$Stiffness_about_Tension	[expr	$about_Tension]
uniaxialMaterial	Elastic	$Stiffness_about_IP			[expr	$about_IP]
uniaxialMaterial	Elastic	$Stiffness_about_OOP		[expr	$about_OOP]


puts "The materials are defined"