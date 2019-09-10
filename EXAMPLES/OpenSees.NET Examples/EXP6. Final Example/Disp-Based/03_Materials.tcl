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
set Long_Shear_Timoshenko_Parallel	72
set Long_Shear_Timoshenko_Serial	73


uniaxialMaterial	Elastic	$rigidmattag  						[expr 1.0e+15]
uniaxialMaterial	Elastic	$semi_rigidmattag					[expr 1.0e+4]

uniaxialMaterial	Elastic	$Freemattag  						[expr 1.0e-1]


uniaxialMaterial 	ElasticPPGap $timbermattag_long $E_Timber_Long 430. 20.
# uniaxialMaterial	Elastic	$timbermattag_long					$E_Timber_Long

uniaxialMaterial 	ElasticMultiLinear $timbermattag_trans -strain -100.	-50.	-20.	0.	5.	12.	47.	60.	120. -stress -100.	-50.	-20.	0.	50.	120.	4700.	600.	120.
# uniaxialMaterial	Elastic	$timbermattag_trans					$E_Timber_Trans

uniaxialMaterial 	ElasticBilin $G_Timber_Tag_Long 80. 10. 12.0 50. 20. 9.0
# uniaxialMaterial	Elastic	$G_Timber_Tag_Long					$G_Timber_Long

# uniaxialMaterial 	ENT $G_Timber_Tag_Trans $G_Timber_Trans
uniaxialMaterial	Elastic	$G_Timber_Tag_Trans					$G_Timber_Trans

uniaxialMaterial	Elastic	$Stiffness_Tension					[expr	$Tension]
uniaxialMaterial	Elastic	$Stiffness_IP						[expr	$IP]
uniaxialMaterial	Elastic	$Stiffness_OOP						[expr	$OOP]
uniaxialMaterial	Elastic	$Stiffness_about_Tension			[expr	$about_Tension]
uniaxialMaterial	Elastic	$Stiffness_about_IP					[expr	$about_IP]
uniaxialMaterial	Elastic	$Stiffness_about_OOP				[expr	$about_OOP]


uniaxialMaterial	Elastic	$Trans_Shear_Timoshenko_Eq	[expr 34.878]
uniaxialMaterial	Elastic	$Long_Shear_Timoshenko_Eq	[expr 11.971]


uniaxialMaterial Parallel $Long_Shear_Timoshenko_Parallel $Trans_Shear_Timoshenko_Eq $Long_Shear_Timoshenko_Eq
uniaxialMaterial Series $Long_Shear_Timoshenko_Serial $Trans_Shear_Timoshenko_Eq $Trans_Shear_Timoshenko_Eq

puts "The materials are defined"

puts "The materials are defined"