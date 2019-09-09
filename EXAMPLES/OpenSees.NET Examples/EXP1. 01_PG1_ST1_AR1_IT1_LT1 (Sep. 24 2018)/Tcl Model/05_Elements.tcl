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

# Defining Timoshenko Beam

# element ElasticTimoshenkoBeam	1	35	33	$E_Timber_Long		$G_Timber	$Area_Timoshenko_long	$Torsional_J_Timoshenko_long	$Iy_Timoshenko_long		$Iz_Timoshenko_long		$Av_Timoshenko_long		$Av_Timoshenko_long		$lineartag_Timoshenko
# element ElasticTimoshenkoBeam	58	1	4	$E_Timber_Trans		$G_Timber	$Area_Timoshenko_trans	$Torsional_J_Timoshenko_trans	$Iy_Timoshenko_trans	$Iz_Timoshenko_trans	$Av_Timoshenko_trans	$Av_Timoshenko_trans	$lineartag_Timoshenko

element twoNodeLink	1	35	33	-mat	$Long_Shear_Timoshenko_Eq	-dir	2	-orient	1	0	0
element twoNodeLink	58	1	4	-mat	$Trans_Shear_Timoshenko_Eq	-dir	2	-orient	0	0	1


puts "The Timoshenko beams are defined"
# ************************************************************************************************************

# Defining Link Elements Between Zero-Length Nodes

element twoNodeLink	2	16		1601	-mat	$Axial_perimeter_long_left		$Shear_1_perimeter_long_left		$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	0	0	1	1	0	0
element twoNodeLink	3	1602	16		-mat	$Axial_perimeter_trans_down		$Shear_1_perimeter_trans_down		$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	-1	0	0	0	0	1

element twoNodeLink	4	1301	13		-mat	$Axial_perimeter_long_left		$Shear_1_perimeter_long_left		$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	0	0	1	1	0	0
element twoNodeLink	5	1302	13		-mat	$Axial_perimeter_trans_up		$Shear_1_perimeter_trans_up			$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	-1	0	0	0	0	1

element twoNodeLink	6	11		1101	-mat	$Axial_perimeter_long_right		$Shear_1_perimeter_long_right		$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	0	0	1	1	0	0
element twoNodeLink	7	11		1102	-mat	$Axial_perimeter_trans_down		$Shear_1_perimeter_trans_down		$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	-1	0	0	0	0	1

element twoNodeLink	8	801		8		-mat	$Axial_perimeter_long_right		$Shear_1_perimeter_long_right		$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	0	0	1	1	0	0
element twoNodeLink	9	8		802		-mat	$Axial_perimeter_trans_up		$Shear_1_perimeter_trans_up			$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	-1	0	0	0	0	1

element twoNodeLink	10	2		202		-mat	$Axial_inside_trans		$Shear_1_inside_trans	$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	-1	0	0	0	0	1
element twoNodeLink	11	9		902		-mat	$Axial_inside_trans		$Shear_1_inside_trans	$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	-1	0	0	0	0	1
element twoNodeLink	12	502		5		-mat	$Axial_inside_trans		$Shear_1_inside_trans	$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	-1	0	0	0	0	1
element twoNodeLink	13	1402	14		-mat	$Axial_inside_trans		$Shear_1_inside_trans	$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	-1	0	0	0	0	1

element twoNodeLink	14	36		3601	-mat	$Axial_inside_long		$Shear_1_inside_long	$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	0	0	1	1	0	0
element twoNodeLink	15	34		3401	-mat	$Axial_inside_long		$Shear_1_inside_long	$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	0	0	1	1	0	0
element twoNodeLink	16	3201	32		-mat	$Axial_inside_long		$Shear_1_inside_long	$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	0	0	1	1	0	0
element twoNodeLink	17	3101	31		-mat	$Axial_inside_long		$Shear_1_inside_long	$rigidmattag	$Freemattag	$Freemattag	$Freemattag	-dir	1	2	3	4	5	6 -orient	0	0	1	1	0	0

puts "The Link elements for zero-length nodes are defined"
# ************************************************************************************************************

# Defining Vertical Perimeter Elements

element ElasticTimoshenkoBeam	18	1601	15		$E_Timber_Long $G_Timber $Area_Perimeter_long_left $Torsional_J_Perimeter_long_left $Iy_Perimeter_long_left $Iz_Perimeter_long_left $Av_Perimeter_long_left $Av_Perimeter_long_left $lineartag_Vertical
element ElasticTimoshenkoBeam	19	15		14 		$E_Timber_Long $G_Timber $Area_Perimeter_long_left $Torsional_J_Perimeter_long_left $Iy_Perimeter_long_left $Iz_Perimeter_long_left $Av_Perimeter_long_left $Av_Perimeter_long_left $lineartag_Vertical
element ElasticTimoshenkoBeam	20	14		6 		$E_Timber_Long $G_Timber $Area_Perimeter_long_left $Torsional_J_Perimeter_long_left $Iy_Perimeter_long_left $Iz_Perimeter_long_left $Av_Perimeter_long_left $Av_Perimeter_long_left $lineartag_Vertical
element ElasticTimoshenkoBeam	21	6		4 		$E_Timber_Long $G_Timber $Area_Perimeter_long_left $Torsional_J_Perimeter_long_left $Iy_Perimeter_long_left $Iz_Perimeter_long_left $Av_Perimeter_long_left $Av_Perimeter_long_left $lineartag_Vertical
element ElasticTimoshenkoBeam	22	4		5 		$E_Timber_Long $G_Timber $Area_Perimeter_long_left $Torsional_J_Perimeter_long_left $Iy_Perimeter_long_left $Iz_Perimeter_long_left $Av_Perimeter_long_left $Av_Perimeter_long_left $lineartag_Vertical
element ElasticTimoshenkoBeam	23	5		12 		$E_Timber_Long $G_Timber $Area_Perimeter_long_left $Torsional_J_Perimeter_long_left $Iy_Perimeter_long_left $Iz_Perimeter_long_left $Av_Perimeter_long_left $Av_Perimeter_long_left $lineartag_Vertical
element ElasticTimoshenkoBeam	24	12		1301 	$E_Timber_Long $G_Timber $Area_Perimeter_long_left $Torsional_J_Perimeter_long_left $Iy_Perimeter_long_left $Iz_Perimeter_long_left $Av_Perimeter_long_left $Av_Perimeter_long_left $lineartag_Vertical

element ElasticTimoshenkoBeam	25	1101	10		$E_Timber_Long $G_Timber $Area_Perimeter_long_right $Torsional_J_Perimeter_long_right $Iy_Perimeter_long_right $Iz_Perimeter_long_right $Av_Perimeter_long_right $Av_Perimeter_long_right $lineartag_Vertical
element ElasticTimoshenkoBeam	26	10		2 		$E_Timber_Long $G_Timber $Area_Perimeter_long_right $Torsional_J_Perimeter_long_right $Iy_Perimeter_long_right $Iz_Perimeter_long_right $Av_Perimeter_long_right $Av_Perimeter_long_right $lineartag_Vertical
element ElasticTimoshenkoBeam	27	2		1 		$E_Timber_Long $G_Timber $Area_Perimeter_long_right $Torsional_J_Perimeter_long_right $Iy_Perimeter_long_right $Iz_Perimeter_long_right $Av_Perimeter_long_right $Av_Perimeter_long_right $lineartag_Vertical
element ElasticTimoshenkoBeam	28	1		3 		$E_Timber_Long $G_Timber $Area_Perimeter_long_right $Torsional_J_Perimeter_long_right $Iy_Perimeter_long_right $Iz_Perimeter_long_right $Av_Perimeter_long_right $Av_Perimeter_long_right $lineartag_Vertical
element ElasticTimoshenkoBeam	29	3		9 		$E_Timber_Long $G_Timber $Area_Perimeter_long_right $Torsional_J_Perimeter_long_right $Iy_Perimeter_long_right $Iz_Perimeter_long_right $Av_Perimeter_long_right $Av_Perimeter_long_right $lineartag_Vertical
element ElasticTimoshenkoBeam	30	9		7 		$E_Timber_Long $G_Timber $Area_Perimeter_long_right $Torsional_J_Perimeter_long_right $Iy_Perimeter_long_right $Iz_Perimeter_long_right $Av_Perimeter_long_right $Av_Perimeter_long_right $lineartag_Vertical
element ElasticTimoshenkoBeam	31	7		801 	$E_Timber_Long $G_Timber $Area_Perimeter_long_right $Torsional_J_Perimeter_long_right $Iy_Perimeter_long_right $Iz_Perimeter_long_right $Av_Perimeter_long_right $Av_Perimeter_long_right $lineartag_Vertical


puts "The vertical Perimeter elements are defined"
# ************************************************************************************************************

# Defining Horizontal Perimeter Elements

element ElasticTimoshenkoBeam	32	1102	30		$E_Timber_Trans $G_Timber $Area_Perimeter_trans_down $Torsional_J_Perimeter_trans_down $Iy_Perimeter_trans_down $Iz_Perimeter_trans_down $Av_Perimeter_trans_down $Av_Perimeter_trans_down $lineartag_Horizontal
element ElasticTimoshenkoBeam	33	30		34		$E_Timber_Trans $G_Timber $Area_Perimeter_trans_down $Torsional_J_Perimeter_trans_down $Iy_Perimeter_trans_down $Iz_Perimeter_trans_down $Av_Perimeter_trans_down $Av_Perimeter_trans_down $lineartag_Horizontal
element ElasticTimoshenkoBeam	34	34		35		$E_Timber_Trans $G_Timber $Area_Perimeter_trans_down $Torsional_J_Perimeter_trans_down $Iy_Perimeter_trans_down $Iz_Perimeter_trans_down $Av_Perimeter_trans_down $Av_Perimeter_trans_down $lineartag_Horizontal
element ElasticTimoshenkoBeam	35	35		28		$E_Timber_Trans $G_Timber $Area_Perimeter_trans_down $Torsional_J_Perimeter_trans_down $Iy_Perimeter_trans_down $Iz_Perimeter_trans_down $Av_Perimeter_trans_down $Av_Perimeter_trans_down $lineartag_Horizontal
element ElasticTimoshenkoBeam	36	28		36		$E_Timber_Trans $G_Timber $Area_Perimeter_trans_down $Torsional_J_Perimeter_trans_down $Iy_Perimeter_trans_down $Iz_Perimeter_trans_down $Av_Perimeter_trans_down $Av_Perimeter_trans_down $lineartag_Horizontal
element ElasticTimoshenkoBeam	37	36		1602	$E_Timber_Trans $G_Timber $Area_Perimeter_trans_down $Torsional_J_Perimeter_trans_down $Iy_Perimeter_trans_down $Iz_Perimeter_trans_down $Av_Perimeter_trans_down $Av_Perimeter_trans_down $lineartag_Horizontal


element ElasticTimoshenkoBeam	38	802		31		$E_Timber_Trans $G_Timber $Area_Perimeter_trans_up $Torsional_J_Perimeter_trans_up $Iy_Perimeter_trans_up $Iz_Perimeter_trans_up $Av_Perimeter_trans_up $Av_Perimeter_trans_up $lineartag_Horizontal
element ElasticTimoshenkoBeam	39	31		23		$E_Timber_Trans $G_Timber $Area_Perimeter_trans_up $Torsional_J_Perimeter_trans_up $Iy_Perimeter_trans_up $Iz_Perimeter_trans_up $Av_Perimeter_trans_up $Av_Perimeter_trans_up $lineartag_Horizontal
element ElasticTimoshenkoBeam	40	23		33		$E_Timber_Trans $G_Timber $Area_Perimeter_trans_up $Torsional_J_Perimeter_trans_up $Iy_Perimeter_trans_up $Iz_Perimeter_trans_up $Av_Perimeter_trans_up $Av_Perimeter_trans_up $lineartag_Horizontal
element ElasticTimoshenkoBeam	41	33		32		$E_Timber_Trans $G_Timber $Area_Perimeter_trans_up $Torsional_J_Perimeter_trans_up $Iy_Perimeter_trans_up $Iz_Perimeter_trans_up $Av_Perimeter_trans_up $Av_Perimeter_trans_up $lineartag_Horizontal
element ElasticTimoshenkoBeam	42	32		21		$E_Timber_Trans $G_Timber $Area_Perimeter_trans_up $Torsional_J_Perimeter_trans_up $Iy_Perimeter_trans_up $Iz_Perimeter_trans_up $Av_Perimeter_trans_up $Av_Perimeter_trans_up $lineartag_Horizontal
element ElasticTimoshenkoBeam	43	21		1302	$E_Timber_Trans $G_Timber $Area_Perimeter_trans_up $Torsional_J_Perimeter_trans_up $Iy_Perimeter_trans_up $Iz_Perimeter_trans_up $Av_Perimeter_trans_up $Av_Perimeter_trans_up $lineartag_Horizontal


puts "The horizontal Perimeter elements are defined"
# ************************************************************************************************************

# Defining Vertical Inside Elements

element ElasticTimoshenkoBeam	44	3601	3201	$E_Timber_Long $G_Timber $Area_inside_long $Torsional_J_inside_long $Iy_inside_long $Iz_inside_long $Av_inside_long $Av_inside_long $lineartag_Vertical
element ElasticTimoshenkoBeam	45	3401	3101	$E_Timber_Long $G_Timber $Area_inside_long $Torsional_J_inside_long $Iy_inside_long $Iz_inside_long $Av_inside_long $Av_inside_long $lineartag_Vertical


puts "The Vertical Inside elements are defined"
# ************************************************************************************************************

# Defining Horizontal Inside Elements

element ElasticTimoshenkoBeam	46	202	1402	$E_Timber_Trans $G_Timber $Area_inside_trans $Torsional_J_inside_trans $Iy_inside_trans $Iz_inside_trans $Av_inside_trans $Av_inside_trans $lineartag_Horizontal
element ElasticTimoshenkoBeam	47	902	502		$E_Timber_Trans $G_Timber $Area_inside_trans $Torsional_J_inside_trans $Iy_inside_trans $Iz_inside_trans $Av_inside_trans $Av_inside_trans $lineartag_Horizontal

puts "The Horizontal Inside elements are defined"
# ************************************************************************************************************


# Defining Joint Elements

element twoNodeLink	48	27	28	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink	49	29	30	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0

element twoNodeLink	50	23	22	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
element twoNodeLink	51	21	20	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0


element twoNodeLink	52	19	7	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink	53	18	3	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink	54	17	10	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1

element twoNodeLink	55	12	24	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink	56	6	25	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
element twoNodeLink	57	15	26	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1


puts "The Joint elements are defined"

