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

# ------------------------------ Elements -------------------------------------------------------------------

set lineartag_Timoshenko 1
geomTransf Linear $lineartag_Timoshenko	0	1	0

set lineartag_Vertical 2
geomTransf Linear $lineartag_Vertical	0	1	0

set lineartag_Horizontal 3
geomTransf Linear $lineartag_Horizontal	0	1	0


# ************************************************************************************************************
# Defining Vertical Perimeter Elements



element dispBeamColumn 		1		100		101			5 $width_Perimeter_long_left_tag $lineartag_Vertical

# element elasticBeamColumn 	1		100		101	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	2		101		102	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	3		102		103	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	4		103		104	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	5		104		105	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	6		105		106	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	7		106		107	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	8		107		108	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	9		108		109	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	129		109		110	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	11		110		111	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	12		111		112	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	13		112		113	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	14		113		114	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	15		114		115	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	16		115		116	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	17		116		117	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	18		117		118	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	19		118		119	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	130		119		120	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical


element forceBeamColumn 	20		200		201			5 $width_Perimeter_long_right_tag $lineartag_Vertical ;#******* iman modifications ***********
element dispBeamColumn 		21		201		202			5 $width_Perimeter_long_right_tag $lineartag_Vertical


#******* iman modifications *********** element elasticBeamColumn 	20		200		201	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
# element elasticBeamColumn 	21		201		202	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	22		202		203	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	23		203		204	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	24		204		205	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	25		205		206	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	26		206		207	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	27		207		208	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	28		208		209	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	29		209		210	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	30		210		211	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	31		211		212	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	32		212		213	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	33		213		214	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	34		214		215	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
element elasticBeamColumn 	35		215		216	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
element elasticBeamColumn 	36		216		217	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
element elasticBeamColumn 	37		217		218	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
element elasticBeamColumn 	38		218		219	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
element elasticBeamColumn 	131		219		220	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical







element forceBeamColumn 	200		10001		10101			5 $width_Perimeter_trans_up_tag $lineartag_Vertical ;#******* iman modifications ***********
element dispBeamColumn 		201		10101		10201			5 $width_Perimeter_trans_up_tag $lineartag_Vertical


#******* iman modifications *********** element elasticBeamColumn 	200		10001		10101	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
# element elasticBeamColumn 	201		10101		10201	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	202		10201		10301	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	203		10301		10401	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	204		10401		10501	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	205		10501		10601	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	206		10601		10701	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
element elasticBeamColumn 	207		10701		10801	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	208		10801		10901	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	209		10901		11001	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	210		11001		11101	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	211		11101		11201	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	212		11201		11301	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
element elasticBeamColumn 	213		11301		11401	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	214		11401		11501	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	215		11501		11601	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	216		11601		11701	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	217		11701		11801	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	218		11801		11901	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
element elasticBeamColumn 	219		11901		12001	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical


element elasticBeamColumn 	220		20001		20101	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	221		20101		20201	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	222		20201		20301	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	223		20301		20401	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	224		20401		20501	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	225		20501		20601	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	226		20601		20701	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
element elasticBeamColumn 	227		20701		20801	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	228		20801		20901	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	229		20901		21001	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	230		21001		21101	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	231		21101		21201	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	232		21201		21301	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	233		21301		21401	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
element elasticBeamColumn 	234		21401		21501	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
element elasticBeamColumn 	235		21501		21601	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
element elasticBeamColumn 	236		21601		21701	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
element elasticBeamColumn 	237		21701		21801	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
element elasticBeamColumn 	238		21801		21901	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
element elasticBeamColumn 	239		21901		22001	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical








puts "The vertical Perimeter elements are defined"
# ************************************************************************************************************

# Defining Horizontal Perimeter Elements

element forceBeamColumn 	39		200		301			5 $width_Perimeter_trans_down_tag $lineartag_Horizontal ;#******* iman modifications ***********
element dispBeamColumn 		40		301		302			5 $width_Perimeter_trans_down_tag $lineartag_Horizontal


#******* iman modifications *********** element elasticBeamColumn 	39		200		301			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
# element elasticBeamColumn 	40		301		302			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	41		302		303			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	42		303		304			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	43		304		305			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	44		305		306			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	45		306		307			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	46		307		308			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	47		308		309			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	48		309		310			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	49		310		311			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	50		311		312			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	51		312		313			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	52		313		314			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	53		314		315			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	54		315		100			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal

element elasticBeamColumn 	55		220		401			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	56		401		402			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	57		402		403			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	58		403		404			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	59		404		405			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	60		405		406			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	61		406		407			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	62		407		408			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	63		408		409			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	64		409		410			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	65		410		411			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	66		411		412			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	67		412		413			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	68		413		414			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	69		414		415			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	70		415		120			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal









element forceBeamColumn 	250		20001		30101			5 $width_inside_long_tag $lineartag_Horizontal ;#******* iman modifications ***********
element dispBeamColumn 		251		30101		30201			5 $width_inside_long_tag $lineartag_Horizontal



#******* iman modifications *********** element elasticBeamColumn 	250		20001		30101			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
# element elasticBeamColumn 	251		30101		30201			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	252		30201		30301			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	253		30301		30401			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	254		30401		30501			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	255		30501		30601			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	256		30601		30701			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	257		30701		30801			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	258		30801		30901			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	259		30901		31001			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	260		31001		31101			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	261		31101		31201			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	262		31201		31301			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	263		31301		31401			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	264		31401		31501			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
element elasticBeamColumn 	265		31501		10001			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal

element elasticBeamColumn 	266		22001		40101			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	267		40101		40201			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	268		40201		40301			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	269		40301		40401			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	270		40401		40501			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	271		40501		40601			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	272		40601		40701			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	273		40701		40801			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	274		40801		40901			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	275		40901		41001			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	276		41001		41101			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	277		41101		41201			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	278		41201		41301			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	279		41301		41401			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	280		41401		41501			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
element elasticBeamColumn 	281		41501		12001			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal


puts "The horizontal Perimeter elements are defined"
# ************************************************************************************************************

# Defining Vertical Inside Elements



element forceBeamColumn 	71	315		415			5 $width_inside_trans_tag $lineartag_Vertical ;#******* iman modifications ***********
element dispBeamColumn 		72	310		408			5 $width_inside_trans_tag $lineartag_Vertical



#******* iman modifications *********** element elasticBeamColumn 71	315		415		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
# element elasticBeamColumn 72	310		408		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
element elasticBeamColumn 73	308		405		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
element elasticBeamColumn 74	301		401		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical



element elasticBeamColumn 300	31501		41501		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
element elasticBeamColumn 301	31001		40801		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
element elasticBeamColumn 302	30801		40501		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
element elasticBeamColumn 303	30101		40101		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical


puts "The Vertical Inside elements are defined"
# ************************************************************************************************************

# Defining Horizontal Inside Elements

element elasticBeamColumn	75	201		106		$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
element elasticBeamColumn	76	207		107		$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
element elasticBeamColumn	77	213		113		$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
element elasticBeamColumn	78	214		119		$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal


element elasticBeamColumn	304	20101	10601	$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
element elasticBeamColumn	305	20701	10701	$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
element elasticBeamColumn	306	21301	11301	$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
element elasticBeamColumn	307	21401	11901	$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal


puts "The Horizontal Inside elements are defined"
# ************************************************************************************************************


# Defining Joint Elements

element twoNodeLink		79		2020		202		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		80		2030		203		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		81		2040		204		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		82		2050		205		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		83		2060		206		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		84		2080		208		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		85		2090		209		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		86		2100		210		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		87		2110		211		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		88		2120		212		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		89		2150		215		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		90		2160		216		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		91		2170		217		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		92		2180		218		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		93		2190		219		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1


element twoNodeLink		94		101			1010	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		95		102			1020	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		96		103			1030	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		97		104			1040	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		98		105			1050	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		99		108			1080	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		100		109			1090	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		101		110			1100	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		102		111			1110	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		103		112			1120	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		104		114			1140	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		105		115			1150	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		106		116			1160	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		107		117			1170	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink		108		118			1180	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1


element twoNodeLink		109		3020		302		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		110		3040		304		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		111		3050		305		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		112		3060		306		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		113		3070		307		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		114		3090		309		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		115		3110		311		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		116		3120		312		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		117		3130		313		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		118		3140		314		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
	
		
element twoNodeLink		119		402			4020	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		120		403			4030	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		121		404			4040	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		122		406			4060	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		123		407			4070	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		124		409			4090	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		125		411			4110	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		126		413			4130	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		127		414			4140	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink		128		412			4120	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0


puts "The Joint elements are defined"

