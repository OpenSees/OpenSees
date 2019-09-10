# ************************************************************************************************************
#              ***************************************************************************
#                                  *************************************
#				Mechanical Characterization of Integrally-Attached Timber Plate Structures
#                                   Aryan Rezaie Rad, EPFL, Suisse, 2018                  
#                                  *************************************
#              ***************************************************************************
# ************************************************************************************************************
#                           Section Name : Material Definition


set Stiffness_Tension				10
set Stiffness_IP					11
set Stiffness_OOP					12
set Stiffness_about_Tension			13
set Stiffness_about_IP				14
set Stiffness_about_OOP				15


uniaxialMaterial	Elastic	$Stiffness_Tension					[expr	$Tension]
uniaxialMaterial	Elastic	$Stiffness_IP						[expr	$IP]
uniaxialMaterial	Elastic	$Stiffness_OOP						[expr	$OOP]
uniaxialMaterial	Elastic	$Stiffness_about_Tension			[expr	$about_Tension]
uniaxialMaterial	Elastic	$Stiffness_about_IP					[expr	$about_IP]
uniaxialMaterial	Elastic	$Stiffness_about_OOP				[expr	$about_OOP]
