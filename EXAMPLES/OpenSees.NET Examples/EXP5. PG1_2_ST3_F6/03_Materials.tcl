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


set rigidmattag						1
set semi_rigidmattag				101

set Freemattag						2

set timbermattag_long				30
set timbermattag_trans				31
set G_Timber_Tag_Long				32
set G_Timber_Tag_Trans				33

set timbermattag_Ortho_Long			21
set timbermattag_Ortho_Trans		22


set Stiffness_Tension				10
set Stiffness_IP					11
set Stiffness_OOP					12
set Stiffness_about_Tension			13
set Stiffness_about_IP				14
set Stiffness_about_OOP				15

set Trans_Shear_Timoshenko_Eq		70
set Long_Shear_Timoshenko_Eq		71


uniaxialMaterial	Elastic	$rigidmattag  						[expr 1.0e+15]
uniaxialMaterial	Elastic	$semi_rigidmattag					[expr 1.0e+4]

uniaxialMaterial	Elastic	$Freemattag  						[expr 1.0e-1]

uniaxialMaterial	Elastic	$timbermattag_long					$E_Timber_Long
uniaxialMaterial	Elastic	$timbermattag_trans					$E_Timber_Trans
uniaxialMaterial	Elastic	$G_Timber_Tag_Long					$G_Timber_Long
uniaxialMaterial	Elastic	$G_Timber_Tag_Trans					$G_Timber_Trans

uniaxialMaterial	Elastic	$Stiffness_Tension					[expr $Tension]
uniaxialMaterial	Elastic	$Stiffness_IP						[expr $IP]
uniaxialMaterial	Elastic	$Stiffness_OOP						[expr $OOP]
uniaxialMaterial	Elastic	$Stiffness_about_Tension			[expr $about_Tension]
uniaxialMaterial	Elastic	$Stiffness_about_IP					[expr $about_IP]
uniaxialMaterial	Elastic	$Stiffness_about_OOP				[expr $about_OOP]


uniaxialMaterial	Elastic	$Trans_Shear_Timoshenko_Eq	[expr 34787.878]
uniaxialMaterial	Elastic	$Long_Shear_Timoshenko_Eq	[expr 16081.971]


puts "The materials are defined"

puts "The materials are defined"