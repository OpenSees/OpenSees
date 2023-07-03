/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Written by: Rasool Ghorbani, University of Texas at San Antonio
// Created: Jun 2023
// Revision: A

// Description: This file contains the class implementation for GMG_CMAC2D Model. 
// GMG_CMAC2D model is intended for use within a two-dimensional lumped-plasticity 
// framework and aimed at capturing the lateral cyclic behavior of non-retrofitted 
// and retrofitted reinforced concrete columns subjected to seismic motions, up to 
// complete loss of lateral strength.

// The material is based on newly developed computational tools for decision-oriented reinforced concrete column simulations 
// [Ghorbani2022] Ghorbani, R. (2022). "Computational Framework for Decision-Oriented Reinforced Concrete Column Simulation Capabilities". PhD Dissertation, The University of Texas at San Antonio. <https://www.proquest.com/docview/2702490424?pq-origsite=gscholar&fromopenview=true> 

#include <iostream>
#include <math.h>
#include <string.h>
#include <array>
#include <algorithm>

#include <elementAPI.h>
#include <GMG_CMAC2D.h>
#include <Vector.h>
#include <Channel.h>
#include <MaterialResponse.h>
#include <vector>
#include <numeric>


#include <OPS_Globals.h>

Matrix GMG_CMAC2D::Damage_Data(1000000, 16);

// Static Variables
int GMG_CMAC2D::Counter_MACE = 0;

// Static Variables
static int numGMG_CMAC2D = 0;

#ifndef fmax
#define max(a,b) (((a)>(b)) ? (a) : (b))
#endif

#ifndef fmin
#define min(a,b) (((a)<(b)) ? (a) : (b))
#endif

void *
OPS_GMG_CMAC2D()
{
	if (numGMG_CMAC2D == 0) {
		numGMG_CMAC2D++;
		opserr << "GMG_CMAC2D Material Model - Simulating cyclic behavior of non-retrofitted and retrofitted reinforced concrete columns subjected to seismic motions, up to complete loss of lateral strength\n";
		opserr << "Written by: Rasool Ghorbani, University of Texas at San Antonio - 2023\n";
	}

	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial = 0;
	int argc = OPS_GetNumRemainingInputArgs();

	if (!(argc == 24)) {
		opserr << "WARNING GMG_CMAC2D -- insufficient arguments\n";
		opserr << "For direct input, the material needs:\n";
		opserr << "UniaxialMaterial GMG_CMAC2D\n";
		opserr << "(1) Material Tag?\n";
		opserr << "Section Geometry; Section height, section width, and clear cover:\n";
		opserr << "(3) H B C_C\n";
		opserr << "Column Geometry; Column height and shear span\n";
		opserr << "(2) L La\n";
		opserr << "Concrete Property; Compressive strength of concrete:\n";
		opserr << "(1) f'c\n";
		opserr << "Longitudinal Bars; Yielding strength of the longitudinal bars, diameter of the longitudinal bar, number of longitudinal bars along the width of the section (only count one side), and number of longitudinal bars along the height of the section (only count one side of the section)\n";
		opserr << "(4) fyL dbL nL_B nL_H\n";
		opserr << "Transverse Bars; Yield strength of the transverse reinforcement, \n";
		opserr << "(4) fyT dbT n_leg S\n";
		opserr << "Splice\n";
		opserr << "(1) lb\n";
		opserr << "Axial\n";
		opserr << "(1) P_Axial\n";
		opserr << "Interaction Flag for moment and shear strength\n";
		opserr << "(1) PVM_Flag\n";
		opserr << "Number of Averaged Steps\n";
		opserr << "(2) number_of_averaged_steps_Elastic, number_of_averaged_steps_Hardening\n";
		opserr << "Retrofit Flag\n";
		opserr << "(1) Retrofit_flag\n";
		opserr << "Retrofit Properties\n";
		opserr << "(3) Strength thickness bj\n";
		return 0;
	}

	int iTagData[1];

	double Sec_prop[3];
	double Col_geom[2];
	double Conc_prop[1];
	double L_bar[4];
	double T_bar[4];
	double Splice[1];
	double dMonoData_Axial[1];
	int PVM_Flag[1];
	int number_of_averaged_elements[2];
	int Retrofit_flag[1];
	double Retrofit[3];

	int numData;

	numData = 1;
	if (OPS_GetIntInput(&numData, iTagData) != 0) {
		opserr << "WARNING GMG_CMAC2D -- invalid uniaxialMaterial matTag\n";
		return 0;
	}

	numData = 3;
	if (OPS_GetDoubleInput(&numData, Sec_prop)) {
		opserr << "WARNING GMG_CMAC2D -- invalid section properties; section height, section width, clear cover\n";
		opserr << "H? B? C_C?\n";
		return 0;
	}

	numData = 2;
	if (OPS_GetDoubleInput(&numData, Col_geom)) {
		opserr << "WARNING GMG_CMAC2D -- invalid column geometry; total height, shear span\n";
		opserr << "L? La?\n";
		return 0;
	}

	numData = 1;
	if (OPS_GetDoubleInput(&numData, Conc_prop)) {
		opserr << "WARNING GMG_CMAC2D -- invalid compressive strength\n";
		opserr << "f'c?\n";
		return 0;
	}

	numData = 4;
	if (OPS_GetDoubleInput(&numData, L_bar)) {
		opserr << "WARNING GMG_CMAC2D -- invalid section properties - longitudinal bar\n";
		opserr << "fyL? dbL? nL_B? nL_H?\n";
		return 0;
	}

	numData = 4;
	if (OPS_GetDoubleInput(&numData, T_bar)) {
		opserr << "WARNING GMG_CMAC2D -- invalid section properties - transverse reinforcement\n";
		opserr << "fyT? dbT? n_leg? S?\n";
		return 0;
	}

	numData = 1;
	if (OPS_GetDoubleInput(&numData, Splice)) {
		opserr << "WARNING GMG_CMAC2D -- invalid section properties - Provided lap splice\n";
		opserr << "l_spl?\n";
		return 0;
	}
	numData = 1;
	if (OPS_GetDoubleInput(&numData, dMonoData_Axial)) {
		opserr << "WARNING GMG_CMAC2D -- invalid initial axial force\n";
		opserr << "P_Axial?\n";
		return 0;
	}

	numData = 1;
	if (OPS_GetIntInput(&numData, PVM_Flag) != 0) {
		opserr << "WARNING GMG_CMACA2D -- invalid axial/shear/Moment flag\n";
		opserr << "PVM_Flag?\n";
		return 0;
	}

	numData = 2;
	if (OPS_GetIntInput(&numData, number_of_averaged_elements) != 0) {
		opserr << "WARNING GMG_CMAC2D -- invalid number of prior converged steps considered to obtain the average axial demand that is used to calculate yield and capping transition surfaces\n";
		opserr << "Number_of_Step_ElasticBranch? Number_of_Step_HardeningBranch?\n";
		return 0;
	}

	numData = 1;
	if (OPS_GetIntInput(&numData, Retrofit_flag) != 0) {
		opserr << "WARNING GMG_CMAC2D -- invalid rerofit flag; 0 = Non-rerofitted; 1 = FRP jacketed; 2 = Steel jacketed\n";
		opserr << "Retrofit_flag?\n";
		return 0;
	}

	numData = 3;
	if (OPS_GetDoubleInput(&numData, Retrofit)) {
		opserr << "WARNING GMG_CMAC2D -- invalid retrofit properties\n";
		opserr << "f_je? t_j? bj?\n";
		return 0;
	}

	// Parsing was successful, allocate the material
	theMaterial = new GMG_CMAC2D(iTagData[0],

		Sec_prop[0], Sec_prop[1], Sec_prop[2],
		Col_geom[0], Col_geom[1],
		Conc_prop[0],
		L_bar[0], L_bar[1], L_bar[2], L_bar[3],
		T_bar[0], T_bar[1], T_bar[2], T_bar[3],
		Splice[0],
		dMonoData_Axial[0],
		PVM_Flag[0],
		number_of_averaged_elements[0], number_of_averaged_elements[1],
		Retrofit_flag[0],
		Retrofit[0], Retrofit[1], Retrofit[2]);
	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type EnhancedSpring Material\n";
		return 0;
	}
	return theMaterial;
}

// Constructors
GMG_CMAC2D::GMG_CMAC2D(int tag,

	double H, double B, double C_C,
	double L, double La,
	double fc,
	double fyL, double dbL, double nL_EW, double nL_NS,
	double fyT, double dbT, double n_leg, double S,
	double lb,
	double P_Axial,
	int PVM_Flag,
	int number_of_averaged_elements_Elastic,
	int number_of_averaged_elements_Hardening,
	int Retrofit_flag,
	double fje,
	double tj,
	double bj)

	:UniaxialMaterial(tag, MAT_TAG_GMG_CMAC2D),
	Strain_vector(3),
	Tang_vector(3),
	Stress_vector(3),
	Damp_vector(3),

	H(H), B(B), C_C(C_C),
	L(L), La(La),
	fc(fc),
	fyL(fyL), dbL(dbL), nL_EW(nL_EW), nL_NS(nL_NS),
	fyT(fyT), dbT(dbT), n_leg(n_leg), S(S),
	lb(lb),
	P_Axial(P_Axial),
	PVM_Flag(PVM_Flag),
	number_of_averaged_elements_Elastic(number_of_averaged_elements_Elastic),
	number_of_averaged_elements_Hardening(number_of_averaged_elements_Hardening),
	Retrofit_flag(Retrofit_flag),
	fje(fje),
	tj(tj),
	bj(bj)
{
	// Initialize Variables by Calling revertToStart function
	this->DamageOutput();
	this->P_M_yield(-1.0*P_Axial);
	this->P_M_capping(-1.0*P_Axial);
	this->P_M_splice(-1.0*P_Axial);
	this->Shear_strength(-1.0*P_Axial);
	this->revertToStart();
	//revertToStart();
}

GMG_CMAC2D::GMG_CMAC2D()
	:UniaxialMaterial(0, MAT_TAG_GMG_CMAC2D),
	Strain_vector(3),
	Tang_vector(3),
	Stress_vector(3),
	Damp_vector(3),

	H(0.0), B(0.0), C_C(0.0),
	L(0.0), La(0.0),
	fc(0.0),
	fyL(0.0), dbL(0.0), nL_EW(0.0), nL_NS(0.0),
	fyT(0.0), dbT(0.0), n_leg(0.0), S(0.0),
	lb(0.0),
	P_Axial(0.0),
	Kepos_Rot(0.0), Keneg_Rot(0.0), fypos_Rot(0.0), fyneg_Rot(0.0), fcappos_Rot(0.0),
	fcapneg_Rot(0.0), dcappos_Rot(0.0), dcapneg_Rot(0.0),
	Kdegpos_Rot(0.0), Kdegneg_rot(0.0),
	frespos_Rot(0.0), fresneg_Rot(0.0),
	delUpos_Rot(0.0), delUneg_Rot(0.0),
	alpha_Er_Rot_Post_Yielding(0.0), beta_Er_Rot_Post_Yielding(0.0),
	alpha_Er_Rot_Post_Capping(0.0), beta_Er_Rot_Post_Capping(0.0),
	Er_Rot_Post_Yielding(0.0), Er_Rot_Post_Capping(0.0),
	Kun_Rot_Post_Yielding(0.0), Kun_Rot_Post_Capping(0.0),
	Kr_Rot_Post_Yielding(0.0), Kr_Rot_Post_Capping(0.0),
	delta_ratio_max_hard_Rot(0.0), Ref_Energy_Coe_Rot(0.0),
	C1_Rot(0.0), C2_Rot(0.0), C3_Rot(0.0), solpe_post_yielding_Rot(0.0), solpe_post_capping_Rot(0.0),

	fypos_splice_Rot(0.0), fyneg_splice_Rot(0.0), fypos_splice_deg_Rot(0.0), fyneg_splice_deg_Rot(0.0),
	fcappos_splice_Rot(0.0), fcapneg_splice_Rot(0.0), dcappos_splice_Rot(0.0), dcapneg_splice_Rot(0.0),
	Kdegpos_splice_Rot(0.0), Kdegneg_splice_rot(0.0),
	frespos_splice_Rot(0.0), fresneg_splice_Rot(0.0),
	delUpos_splice_Rot(0.0), delUneg_splice_Rot(0.0),
	alpha_Er_Rot_Post_Capping_splice(0.0), beta_Er_Rot_Post_Capping_splice(0.0),
	Er_Rot_Post_Capping_splice(0.0),
	Kun_Rot_Post_Capping_splice(0.0),
	Kr_Rot_Post_Capping_splice(0.0),
	C1_splice_Rot(0.0), C2_splice_Rot(0.0), C3_splice_Rot(0.0), solpe_post_capping_splice_Rot(0.0),

	dcappos_FlexShr_Rot(0.0), dcapneg_FlexShr_Rot(0.0),
	Kdegpos_FlexShr_Rot(0.0), Kdegneg_FlexShr_rot(0.0),
	frespos_FlexShr_Rot(0.0), fresneg_FlexShr_Rot(0.0),
	delUpos_FlexShr_Rot(0.0), delUneg_FlexShr_Rot(0.0),
	alpha_Er_Rot_Post_Capping_FlexShr(0.0), beta_Er_Rot_Post_Capping_FlexShr(0.0),
	Er_Rot_Post_Capping_FlexShr(0.0),
	Kun_Rot_Post_Capping_FlexShr(0.0),
	Kr_Rot_Post_Capping_FlexShr(0.0),
	C1_FlexShr_Rot(0.0), C2_FlexShr_Rot(0.0), C3_FlexShr_Rot(0.0), solpe_post_capping_FlexShr_Rot(0.0),

	Kepos_Shr(0.0), Keneg_Shr(0.0), fypos_Shr(0.0), fyneg_Shr(0.0), fcappos_Shr(0.0),
	fcapneg_Shr(0.0), dcappos_Shr(0.0), dcapneg_Shr(0.0), Kdegpos_Shr(0.0), Kdegneg_Shr(0.0), frespos_Shr(0.0), fresneg_Shr(0.0),
	delUpos_Shr(0.0), delUneg_Shr(0.0),
	Er_Shr_Post_Yielding(0.0), Er_Shr_Post_Capping(0.0),
	Kun_Shr_Post_Yielding(0.0), Kun_Shr_Post_Capping(0.0),
	Kr_Shr_Post_Yielding(0.0), Kr_Shr_Post_Capping(0.0),
	delta_ratio_max_hard_Shr(0.0), Ref_Energy_Coe_Shr(0.0),
	C1_Shr(0.0), C2_Shr(0.0), C3_Shr(0.0), solpe_post_yielding_Shr(0.0), solpe_post_capping_Shr(0.0),

	Kepos_Axil(0.0), Keneg_Axil(0.0),
	PVM_Flag(0),
	number_of_averaged_elements_Elastic(0),
	number_of_averaged_elements_Hardening(0),
	Retrofit_flag(0),
	fje(0.0),
	tj(0.0),
	bj(0.0)


{
	// Initialize Variables by Calling revertToStart function
	this->DamageOutput();
	this->P_M_yield(-1.0*P_Axial);
	this->P_M_capping(-1.0*P_Axial);
	this->P_M_splice(-1.0*P_Axial);
	this->Shear_strength(-1.0*P_Axial);
	this->revertToStart();
	//revertToStart();
}

GMG_CMAC2D::~GMG_CMAC2D()
{
	// does nothing
	ofstream MyFile("Damage Outputs.txt");

	// Write to the file
	MyFile << "Row MaterialTag Flex_Yield Shr_Yield %My %Vy %Ms %Flx %Shr %Flx-Shr %Flx-Splice %Splice FailureType Def_Failure %Strength_Loss Def_Residual" << endln;
	int counter = 1;
	for (int count = 0; count < Counter_MACE; count++) {
		//MyFile << " "<<Damage_Data(count) << "           "<< Damage_Data2(count) << endln;
		if (count % 2 == 0) {
			MyFile << counter << " " << Damage_Data(count, 1) << " " << Damage_Data(count, 2) << " " << Damage_Data(count, 3) << " " << Damage_Data(count, 4) << " " << Damage_Data(count, 5) << " " << Damage_Data(count, 6) << " " << Damage_Data(count, 7) << " " << Damage_Data(count, 8) << " " << Damage_Data(count, 9) << " " << Damage_Data(count, 10) << " " << Damage_Data(count, 11) << " " << Damage_Data(count, 12) << " " << Damage_Data(count, 13) << " " << Damage_Data(count, 14) << " " << Damage_Data(count, 15) << endln;
			counter++;
		}
		//opserr << "Deconstructor " << count << "  " << Damage_Data(count) << endln;
	}

	// Close the file
	MyFile.close();
}

int GMG_CMAC2D::setTrialStrain(double strain, double strainRate)
{
	//does nothing
	//this->revertToStart();
	return 0;
}

int
GMG_CMAC2D::setTrialStrain(const Vector strain_from_element)
{
	// all variables to the last commit state
	//******************************************Adding My Code****************************************
	this->revertToLastCommit();

	/* ********************************************************************************************** **
	************************************************************************************************* **
	**                                                                                                **
	**                                   Dealing with strain vector                                   **
	**                                                                                                **
	************************************************************************************************* **
	** ********************************************************************************************** */

	//double d_Axial, d_Shr, d_Rot, d;
	Strain_vector = strain_from_element;
	d_Shr = strain_from_element(0);
	Strain_vector(0) = d_Shr;

	d_Axial = strain_from_element(1);
	Strain_vector(1) = d_Axial;

	d_Rot = strain_from_element(2);
	//d = d_Rot;
	Strain_vector(2) = d_Rot;

	/* ********************************************************************************************** **
	************************************************************************************************* **
	**                                                                                                **
	**                                    Which Material Goes First                                   **
	**                                                                                                **
	************************************************************************************************* **
	** ********************************************************************************************** */
	int MatCount_Rot = 2; int MatCount_Shr = 3; int MatCount_Axl = 1;

	if (flag_entering_hardening_Flex_Rot != 1 && flag_entering_hardening_Shr != 1 && flag_entering_hardening_Splice_Rot != 1) {
		if (d_Shr > 0 && d_Rot > 0) {
			if (fmax((d_Rot * Kepos_Rot / fypos_Rot), (d_Rot * Kepos_Rot / fypos_splice_Rot)) > (d_Shr * Kepos_Shr / fypos_Shr)) {
				MatCount_Rot = 2;
				MatCount_Shr = 3;
			}
			else {
				MatCount_Rot = 3;
				MatCount_Shr = 2;
			}
		}

		if (d_Shr > 0 && d_Rot < 0) {
			if (fmax((d_Rot * Keneg_Rot / fyneg_Rot), (d_Rot * Keneg_Rot / fyneg_splice_Rot)) > (d_Shr * Kepos_Shr / fypos_Shr)) {
				MatCount_Rot = 2;
				MatCount_Shr = 3;
			}
			else {
				MatCount_Rot = 3;
				MatCount_Shr = 2;
			}
		}

		if (d_Shr < 0 && d_Rot > 0) {
			if (fmax((d_Rot * Kepos_Rot / fypos_Rot), (d_Rot * Kepos_Rot / fypos_splice_Rot)) > (d_Shr * Keneg_Shr / fyneg_Shr)) {
				MatCount_Rot = 2;
				MatCount_Shr = 3;
			}
			else {
				MatCount_Rot = 3;
				MatCount_Shr = 2;
			}
		}

		if (d_Shr < 0 && d_Rot < 0) {
			if (fmax((d_Rot * Keneg_Rot / fyneg_Rot), (d_Rot * Keneg_Rot / fyneg_splice_Rot)) > (d_Shr * Keneg_Shr / fyneg_Shr)) {
				MatCount_Rot = 2;
				MatCount_Shr = 3;
			}
			else {
				MatCount_Rot = 3;
				MatCount_Shr = 2;
			}
		}
	}

	if ((flag_entering_hardening_Flex_Rot == 1 || flag_entering_hardening_Splice_Rot == 1) && flag_entering_hardening_Shr != 1) {
		MatCount_Rot = 2;
		MatCount_Shr = 3;
	}

	if ((flag_entering_hardening_Flex_Rot != 1 || flag_entering_hardening_Splice_Rot != 1) && flag_entering_hardening_Shr == 1) {
		MatCount_Rot = 3;
		MatCount_Shr = 2;
	}
	/* ********************************************************************************************** **
	************************************************************************************************* **
	**                                                                                                **
	**                            End of Which Material Goes First Section                            **
	**                                                                                                **
	************************************************************************************************* **
	** ********************************************************************************************** */



	/* ********************************************************************************************** **
	************************************************************************************************* **
	************************************************************************************************* **
	************************************************************************************************* **
	************************************************************************************************* **
	************************************************************************************************* **
	**                                                                                                **
	**                                                                                                **
	**                                                                                                **
	**                                                                                                **
	**                            Calculating the Stiffness and stress                                **
	**                                                                                                **
	**                                                                                                **
	**                                                                                                **
	**                                                                                                **
	************************************************************************************************* **
	************************************************************************************************* **
	************************************************************************************************* **
	************************************************************************************************* **
	************************************************************************************************* **
	** ********************************************************************************************** */


	for (int MatCount = 1; MatCount <= 3; MatCount++) {
		//this->revertToLastCommit();
		/* ********************************************************************************************** **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		**                                                                                                **
		**                                       Rotational Material                                      **
		**                                                                                                **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		** ********************************************************************************************** */
		if (MatCount == MatCount_Rot) {
			Tdu_Rot = d_Rot - Cstrain_Rot;


			/* ********************************************************************************************** **
			************************************************************************************************* **
			**                                                                                                **
			**                              Check if Splice failure has happend                               **
			**                                                                                                **
			************************************************************************************************* **
			** ********************************************************************************************** */

			if (flag_entering_hardening_Splice_Rot == 1 && flag_Splice_Filure_Rot != 1) {

				if ((d_Rot >= dcappos_splice_Rot && Tdu_Rot >= 0) || (d_Rot <= dcapneg_splice_Rot && Tdu_Rot <= 0)) {
					flag_Splice_Filure_Rot = 1;
					Calibration();
					defineBackbone();
				}
			}

			/* ********************************************************************************************** **
			************************************************************************************************* **
			**                                                                                                **
			**                        Check for Flexure-Splice or Flexure-Shear failure                       **
			**                                     deformation criterion                                      **
			**                                                                                                **
			************************************************************************************************* **
			** ********************************************************************************************** */
			if (flag_entering_hardening_Shr != 1 && flag_entering_hardening_Flex_Rot == 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1 && flag_entering_hardening_Splice_Rot != 1) {

				if (d_Rot >= 0) {

					if (dcappos_Rot <= dcappos_FlexShr_Rot && dcappos_Rot <= dcappos_splice_Rot) {
						if (d_Rot >= dcappos_Rot && Tdu_Rot >= 0) {
							flag_Flx_Filure_Rot = 1;
							Calibration();
							defineBackbone();
						}
					}

					if (dcappos_FlexShr_Rot < dcappos_Rot && dcappos_FlexShr_Rot <= dcappos_splice_Rot) {
						if (d_Rot >= dcappos_FlexShr_Rot && Tdu_Rot >= 0) {
							flag_FlexShr_Filure_Rot = 1;
							Calibration();
							defineBackbone();
						}
					}

					if (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL && dcappos_splice_Rot < dcappos_Rot && dcappos_splice_Rot < dcappos_FlexShr_Rot) {

						if (d_Rot >= dcappos_splice_Rot && Tdu_Rot >= 0) {
							flag_FlexSplice_Filure_Rot = 1;
							Calibration();
							defineBackbone();
						}
					}
				}

				if (d_Rot < 0) {

					if (dcapneg_Rot >= dcapneg_FlexShr_Rot && dcapneg_Rot >= dcapneg_splice_Rot) {
						if (d_Rot <= dcapneg_Rot && Tdu_Rot <= 0) {
							flag_Flx_Filure_Rot = 1;
							Calibration();
							defineBackbone();
						}
					}

					if (dcapneg_FlexShr_Rot > dcapneg_Rot && dcapneg_FlexShr_Rot >= dcapneg_splice_Rot) {
						if (d_Rot <= dcapneg_FlexShr_Rot && Tdu_Rot <= 0) {
							flag_FlexShr_Filure_Rot = 1;
							Calibration();
							defineBackbone();
						}
					}

					if (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL && dcapneg_splice_Rot > dcapneg_Rot && dcapneg_splice_Rot > dcapneg_FlexShr_Rot) {
						if (d_Rot <= dcapneg_splice_Rot && Tdu_Rot <= 0) {
							flag_FlexSplice_Filure_Rot = 1;
							Calibration();
							defineBackbone();
						}
					}
				}
			}

			/* ********************************************************************************************** **
			************************************************************************************************* **
			**                                                                                                **
			**                                     End of the failure Checks                                  **
			**                                                                                                **
			************************************************************************************************* **
			** ********************************************************************************************** */
			if (TstateFlag == 0 || flag_entering_hardening_Shr == 1)
			{
				flagdmg = 0;
				if (d_Rot >= 0.0) {
					F_Rot = Kepos_Rot * d_Rot;
					ek_Rot = Kepos_Rot;
					if (d_Rot >= dpeakmax && flag_entering_hardening_Shr != 1) {
						TstateFlag = 12;
						if (dpeakmax == dypos) {
							flag_entering_hardening_Flex_Rot = 1;
							flag_entering_hardening_Flex_pos_Rot = 1;
							Calibration();
						}
						if (dpeakmax == dypos_splice) {
							flag_entering_hardening_Splice_Rot = 1;
							flag_entering_hardening_Splice_pos_Rot = 1;
							Calibration();
						}
						defineBackbone();
						F_Rot = slope_pos * d_Rot + Intcpt_slope_pos;
						ek_Rot = slope_pos;

					}

				}
				else if (d_Rot < 0.0) {
					F_Rot = Keneg_Rot * d_Rot;
					ek_Rot = Keneg_Rot;
					if (d_Rot <= dpeakmin && flag_entering_hardening_Shr != 1) {  // if flag_entering_hardening_Shr = 1 means Shear spring is in the hardening
						TstateFlag = -12;
						if (dpeakmin == dyneg) {
							flag_entering_hardening_Flex_Rot = 1;
							flag_entering_hardening_Flex_neg_Rot = 1;
							Calibration();
						}
						if (dpeakmin == dyneg_splice) {
							flag_entering_hardening_Splice_Rot = 1;
							flag_entering_hardening_Splice_neg_Rot = 1;
							Calibration();
						}
						defineBackbone();
						F_Rot = -(slope_neg * fabs(d_Rot) - Intcpt_slope_neg);
						ek_Rot = slope_neg;

					}
				}

			}

			else
			{
				if (flagFlur_Shr == 1) {
					if ((flagFlur_Shr == 1 && CflagFlur_Shr != 1) && ((CstateFlag_Shr != 2 && TstateFlag_Shr == 2) || (CstateFlag_Shr != -2 && TstateFlag_Shr == -2))) {
						d_Bench_Rot = Cstrain;
						f_Bench_Rot = Cstress_Rot;
					}
					F_Rot = Kepos_Rot * (d_Rot - d_Bench_Rot) + f_Bench_Rot;
					ek_Rot = Kepos_Rot;
					if (abs(F_Rot) > abs(f_Bench_Rot)) {
						F_Rot = f_Bench_Rot;
					}
				}

				else {
					TstateFlag = getStateFlag();
					switch (TstateFlag)
					{
					case 1:	F_Rot = slope_pos * d_Rot + Intcpt_slope_pos;
						ek_Rot = slope_pos;
						break;
					case -1: F_Rot = -(slope_neg * fabs(d_Rot) - Intcpt_slope_neg);
						ek_Rot = slope_neg;
						break;
					case 12:
						flag_entering_hardening_Flex_pos_Rot = 1;
						flag_entering_hardening_Splice_pos_Rot = 1;
						define_peak();
						flagdmg = 0;
						flagdmg_Hardening_strength = 0;
						flagdmg_Hardening = 1;
						if (BenMark > 0) {
							ek_Rot = slope_pos;
							F_Rot = slope_pos * d_Rot + Intcpt_slope_pos;
							d12 = d_Rot;
							f12 = F_Rot;
						}
						else if (BenMark < 0) {
							ek_Rot = slope_pos;
							F_Rot = slope_pos * fabs(d_Rot) + Intcpt_slope_pos;
							d_12 = d_Rot;
							f_12 = F_Rot;
						}
						dpeakmax = d_Rot;
						if (flagdmg_Hardening_strength == 1)
							update_damage();
						if (flagdmg_Hardening == 1)
							update_damage_hardeingin();
						break;
					case -12:
						flag_entering_hardening_Flex_neg_Rot = 1;
						flag_entering_hardening_Splice_neg_Rot = 1;
						define_peak();
						flagdmg = 0;
						flagdmg_Hardening_strength = 0;
						flagdmg_Hardening = 1;
						if (BenMark > 0) {
							ek_Rot = slope_neg;
							F_Rot = -(slope_neg * fabs(d_Rot) - Intcpt_slope_neg);
							d12neg = d_Rot;
							f12neg = F_Rot;
						}
						else if (BenMark < 0) {
							ek_Rot = slope_neg;
							F_Rot = -(slope_neg * fabs(d_Rot) - Intcpt_slope_neg);
							d_12neg = d_Rot;
							f_12neg = F_Rot;
						}
						dpeakmin = d_Rot;
						if (flagdmg_Hardening_strength == 1)
							update_damage();
						if (flagdmg_Hardening == 1)
							update_damage_hardeingin();
						break;

					case 2:
						//---------Damage will start after entering case 2 or case -2 for the first time---------
						flag_entering_hardening_Flex_pos_Rot = 0;
						flag_entering_hardening_Splice_pos_Rot = 0;
						flag_entering_hardening_Flex_neg_Rot = 0;
						flag_entering_hardening_Splice_neg_Rot = 0;
						flagdmg = 1;
						flagdmg_Hardening_strength = 0;
						flagdmg_Hardening = 0;
						//---------------------------------------------------------------------------------------
						if (BenMark > 0) {
							ek_Rot = R_Kdegpos;
							F_Rot = R_Kdegpos * fabs(d_Rot) + Intcpt_deg_pos;
							d2 = d_Rot;
							f2 = F_Rot;
						}
						else if (BenMark < 0) {
							ek_Rot = R_Kdegpos;
							F_Rot = R_Kdegpos * fabs(d_Rot) + Intcpt_deg_pos;
							d_2 = d_Rot;
							f_2 = F_Rot;
						}
						if (flagdmg == 1)
							update_damage();
						break;
					case -2:
						//---------Damage will start after entering case 2 or case -2 for the first time---------
						flag_entering_hardening_Flex_pos_Rot = 0;
						flag_entering_hardening_Splice_neg_Rot = 0;
						flag_entering_hardening_Splice_pos_Rot = 0;
						flag_entering_hardening_Flex_neg_Rot = 0;
						flagdmg = 1;
						flagdmg_Hardening_strength = 0;
						flagdmg_Hardening = 0;
						//---------------------------------------------------------------------------------------
						if (BenMark > 0) {
							ek_Rot = R_Kdegneg;
							F_Rot = -(R_Kdegneg * fabs(d_Rot) + Intcpt_deg_neg);
							d2neg = d_Rot;
							f2neg = F_Rot;
						}
						else if (BenMark < 0) {
							ek_Rot = R_Kdegneg;
							F_Rot = -(R_Kdegneg * fabs(d_Rot) + Intcpt_deg_neg);
							d_2neg = d_Rot;
							f_2neg = F_Rot;

						}
						if (flagdmg == 1)
							update_damage();
						break;
					case 3:
						flag_entering_hardening_Flex_pos_Rot = 0;
						flag_entering_hardening_Splice_pos_Rot = 0;
						flag_entering_hardening_Flex_neg_Rot = 0;
						flag_entering_hardening_Splice_neg_Rot = 0;
						flag_entering_residual_Flex_Rot = 1;
						flagdmg = 0;
						flagdmg_Hardening_strength = 0;
						flagdmg_Hardening = 0;

						if (BenMark > 0) {
							ek_Rot = 0.0001;
							F_Rot = R_frespos;
							d3 = d_Rot;
							f3 = F_Rot;
						}
						else if (BenMark < 0) {
							ek_Rot = 0.0001;
							F_Rot = R_frespos;
							d_3 = d_Rot;
							f_3 = F_Rot;
						}
						break;
					case -3:
						flag_entering_hardening_Flex_pos_Rot = 0;
						flag_entering_hardening_Splice_neg_Rot = 0;
						flag_entering_hardening_Flex_neg_Rot = 0;
						flag_entering_hardening_Splice_pos_Rot = 0;
						flag_entering_residual_Flex_Rot = 1;
						flagdmg = 0;
						flagdmg_Hardening_strength = 0;
						flagdmg_Hardening = 0;

						if (BenMark > 0) {
							ek_Rot = 0.0001;
							F_Rot = R_fresneg;
							d3neg = d_Rot;
							f3neg = F_Rot;
						}
						else if (BenMark < 0) {
							ek_Rot = 0.0001;
							F_Rot = R_fresneg;
							d_3neg = d_Rot;
							f_3neg = F_Rot;
						}
						break;
					case 30:
						flagdmg = 0;
						flagdmg_Hardening_strength = 0;
						flagdmg_Hardening = 0;
						if (BenMark > 0) {
							ek_Rot = R_Kdegpos;
							F_Rot = R_Kdegpos * fabs(d_Rot) + Intcpt_res_pos;
						}
						else if (BenMark < 0) {
							ek_Rot = R_Kdegpos;
							F_Rot = R_Kdegpos * fabs(d_Rot) + Intcpt_res_pos;
						}
						break;

					case -30:
						flagdmg = 0;
						flagdmg_Hardening_strength = 0;
						flagdmg_Hardening = 0;

						if (BenMark > 0) {
							ek_Rot = R_Kdegneg;
							F_Rot = -(R_Kdegneg * fabs(d_Rot) + Intcpt_res_neg);
						}
						else if (BenMark < 0) {
							ek_Rot = R_Kdegneg;
							F_Rot = -(R_Kdegneg * fabs(d_Rot) + Intcpt_res_neg);
						}
						break;

					case 40:
						flagdmg = 0;
						flagdmg_Hardening_strength = 0;
						flagdmg_Hardening = 0;

						if (BenMark > 0) {
							ek_Rot = -0.0001;
							F_Rot = 0.0;
						}
						else if (BenMark < 0) {
							ek_Rot = -0.0001;
							F_Rot = 0.0;
						}
						break;

					case -40:

						flagdmg = 0;
						flagdmg_Hardening_strength = 0;
						flagdmg_Hardening = 0;
						if (BenMark > 0) {
							ek_Rot = -0.0001;
							F_Rot = 0.0;
						}
						else if (BenMark < 0) {
							ek_Rot = -0.0001;
							F_Rot = 0.0;
						}
						break;

					case 31:
						if (CstateFlag != 31)
							define_peak();
						ek_Rot = (ffmax - ffmin) / (dpeakmax - dpeakmin);
						F_Rot = ek_Rot * (d_Rot - dpeakmax) + ffmax;

						checkEnvelope();
						break;

					case 41:

						ek_Rot = 0.0001;
						F_Rot = 0.0;
						break;

					case 4:

						define_peak();
						dirtag = 1; // Positive to Negative Direction

						MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
						splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);

						Krel = fmin(pow(MtoRref / ((Kepos_Rot + Keneg_Rot) / 2), KrR) * (Kepos_Rot + Keneg_Rot) / 2, (Kepos_Rot + Keneg_Rot) / 2);
						Kun = fmin(KuR * MtoRref, (Kepos_Rot + Keneg_Rot) / 2);
						E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Kepos_Rot);  // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G

						spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot);
						F_Rot = fspl;
						ek_Rot = ek;

						if (flagdmg == 1)
							update_damage();
						if (flagdmg_Hardening == 1)
							update_damage_hardeingin();
						if (flagdmg_Hardening_strength == 1)
							update_damage();
						dpeakmin_inner = d_Rot;
						ffmin_inner = F_Rot;
						break;

					case 5:
						define_peak();
						spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot);      //Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						F_Rot = fspl;
						ek_Rot = ek;


						if (flagdmg == 1)
							update_damage();
						if (flagdmg_Hardening_strength == 1)
							update_damage();
						if (flagdmg_Hardening == 1)
							update_damage_hardeingin();
						checkEnvelope();
						break;

					case 6:
						define_peak();

						// Because of the damage model the targeted peak points may be updated so we need to re-define spline parametere again 
						if (CstateFlag != 6) {
							MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
							splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);
							Krel = fmin(pow(MtoRref / ((Kepos_Rot + Keneg_Rot) / 2), KrR) * (Kepos_Rot + Keneg_Rot) / 2, (Kepos_Rot + Keneg_Rot) / 2);
							Kun = fmin(KuR * MtoRref, (Kepos_Rot + Keneg_Rot) / 2);
							E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Kepos_Rot);  // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						}
						//------------------------------------------------------------------------------------------------------

						// -----Check if the current point is above or under the new lower limit------
						if (CstateFlag != 6) {

							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain_Rot); // NP
							Benchmark56_up = fspl;

							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain_Rot); // PN
							Benchmark56_down = fspl;
						}
						//---------------------------------------------------------------------------------------

						if ((ffmin_inner <= Benchmark56_up && ffmin_inner >= Benchmark56_down)) {

							if (CstateFlag != 6) {

								x = Cstrain_Rot;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
								fouterPN_min = fspl;
								ekouterPN = ek;

								x = Cstrain_Rot;
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
								fouterNP = fspl;
								ekouterNP = ek;

								InCycFac = (fouterNP - ffmin_inner) / (fouterNP - fouterPN_min);

							}

							x = d_Rot;
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
							fouterNP = fspl;
							ekouterNP = ek;
							// update inner point
							x = dpeakmin + (d_Rot - dpeakmin_inner);
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
							finner = fouterNP - InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)));
							ekinner = (ek - ekouterNP) / ((fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) * (InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) + ekouterNP;
							F_Rot = fmin(finner, fouterNP);
							ek_Rot = ekinner;
							if (F_Rot == fouterNP) {
								ek_Rot = ekouterNP;
								TstateFlag = -5;
							}
						}

						if (ffmin_inner < Benchmark56_down) {

							// Before reaching the lower limit it should be a line
							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot);
							fouterPN = fspl;
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot);
							fouterNP = fspl;
							ekouterNP = ek;
							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN
							F_Rot = ffmin_inner + ek * (d_Rot - dpeakmin_inner);
							ek_Rot = ek;

							F_Rot = fmin(F_Rot, fouterNP);
							if (F_Rot == fouterNP) {
								ek_Rot = ekouterNP;
								TstateFlag = -5;
							}

							// Check if it has reached the lower limit or not, because after it passes the lower limit the coppy-pasting process of the spline curve should get initiated
							if (F_Rot >= fouterPN && Cstress_Rot < fouterPN &&  F_Rot < fouterNP) {
								dpeakmin_inner = d_Rot;
								ffmin_inner = F_Rot;
								x = d_Rot;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
								fouterPN_min = fspl;
								ekouterPN = ek;

								x = d_Rot;
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
								fouterNP = fspl;
								ekouterNP = ek;

								InCycFac = (fouterNP - ffmin_inner) / (fouterNP - fouterPN_min);
							}

							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmin); // NP
							ek_Rot = ek;

							if (F_Rot >= fouterPN && Cstress_Rot > fouterPN && F_Rot < fouterNP) {
								x = d_Rot;
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
								fouterNP = fspl;
								ekouterNP = ek;
								x = dpeakmin + (d_Rot - dpeakmin_inner);
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
								finner = fouterNP - InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)));
								ekinner = (ek - ekouterNP) / ((fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) * (InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) + ekouterNP;
								ek_Rot = ekinner;
								F_Rot = fmin(finner, fouterNP);
								if (F_Rot == fouterNP) {
									ek_Rot = ekouterNP;
									TstateFlag = -5;
								}
							}
						}

						if (ffmin_inner > Benchmark56_up) {

							if (CstateFlag != 6) {
								InCycFac = (ffmin_inner - Benchmark56_up);
							}

							x = d_Rot;
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP 
							F_Rot = fspl + (InCycFac + (-InCycFac / (dpeakmax - dpeakmin_inner))*(d_Rot - dpeakmin_inner));
							ek_Rot = ek - (InCycFac / (dpeakmax - dpeakmin_inner));
						}

						if (flagdmg == 1)
							update_damage();
						if (flagdmg_Hardening_strength == 1)
							update_damage();
						if (flagdmg_Hardening == 1)
							update_damage_hardeingin();
						checkEnvelope();
						break;

					case 7:
						define_peak();

						// Because of the damage model the targeted peak points may be updated so spline parametere should be redefined again 
						if (CstateFlag != 7) {
							MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
							splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);
							Krel = fmin(pow(MtoRref / ((Kepos_Rot + Keneg_Rot) / 2), KrR) * (Kepos_Rot + Keneg_Rot) / 2, (Kepos_Rot + Keneg_Rot) / 2);
							Kun = fmin(KuR * MtoRref, (Kepos_Rot + Keneg_Rot) / 2);
							E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Kepos_Rot);  // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						}
						//------------------------------------------------------------------------------------------------------

						// -----Check if the current point is above or under the new lower limit------
						if (CstateFlag != 7) {
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain_Rot); // NP
							Benchmark67_up = fspl;

							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain_Rot); // PN
							Benchmark67_down = fspl;
						}
						//---------------------------------------------------------------------------------------

						if ((ffmax_inner_inner <= Benchmark67_up && ffmax_inner_inner >= Benchmark67_down)) {

							if (CstateFlag != 7) {

								x = Cstrain_Rot;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
								fouterPN = fspl;
								ekouterPN = ek;

								x = Cstrain_Rot;
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
								fouterNP_max = fspl;
								ekouterNP = ek;

								InCycFac = (ffmax_inner_inner - fouterPN) / (fouterNP_max - fouterPN);

							}

							x = d_Rot;
							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
							fouterPN = fspl;
							ekouterPN = ek;
							// update inner point
							x = -fabs(dpeakmax_inner_inner - d_Rot) + dpeakmax;
							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN

							finner = fouterPN + InCycFac * ((fouterNP_max - fabs(ffmax - fspl)) - fouterPN);
							ekinner = (ek - ekouterPN) / (((fouterNP_max - fabs(ffmax - fspl)) - fouterPN))*(InCycFac*((fouterNP_max - fabs(ffmax - fspl)) - fouterPN)) + ekouterPN;

							ek_Rot = ekinner;
							F_Rot = fmax(finner, fouterPN);
							if (F_Rot == fouterPN) {
								ek_Rot = ekouterPN;
								TstateFlag = 5;
							}
						}

						if (ffmax_inner_inner > Benchmark67_up) {

							// Before reaching the upper limit it should be a line
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot); //NP 
							fouterNP = fspl;
							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot); //PN 
							fouterPN = fspl;
							ekouterPN = ek;

							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN
							F_Rot = ffmax_inner_inner - ek * (dpeakmax_inner_inner - d_Rot);
							ek_Rot = ek;

							F_Rot = fmax(F_Rot, fouterPN);
							if (F_Rot == fouterPN) {
								ek_Rot = ekouterPN;
								TstateFlag = 5;
							}
							// Check if it has reached the upper limit or not, because after it passes the upper limit the coppy-pasting process of the spline curve should get started
							if (F_Rot <= fouterNP && Cstress_Rot > fouterNP && F_Rot > fouterPN) {

								dpeakmax_inner_inner = d_Rot;
								ffmax_inner_inner = F_Rot;

								x = d_Rot;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
								fouterPN = fspl;
								ekouterPN = ek;

								x = d_Rot;
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
								fouterNP_max = fspl;
								ekouterNP = ek;

								InCycFac = (ffmax_inner_inner - fouterPN) / (fouterNP_max - fouterPN);
							}

							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN
							ek_Rot = ek;

							if (F_Rot <= fouterNP && Cstress_Rot < fouterNP && F_Rot > fouterPN) {
								x = d_Rot;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
								fouterPN = fspl;
								ekouterPN = ek;
								// update inner point
								x = -fabs(dpeakmax_inner_inner - d_Rot) + dpeakmax;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
								finner = fouterPN + InCycFac * ((fouterNP_max - fabs(ffmax - fspl)) - fouterPN);
								ekinner = (ek - ekouterPN) / (((fouterNP_max - fabs(ffmax - fspl)) - fouterPN))*(InCycFac*((fouterNP_max - fabs(ffmax - fspl)) - fouterPN)) + ekouterPN;
								ek_Rot = ekinner;
								F_Rot = fmax(finner, fouterPN);
								if (F_Rot == fouterPN) {
									ek_Rot = ekouterPN;
									TstateFlag = 5;
								}
							}
						}

						if (ffmax_inner_inner < Benchmark67_down) {

							if (CstateFlag != 7) {
								InCycFac = (Benchmark67_down - ffmax_inner_inner);
							}

							x = d_Rot;
							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN 
							F_Rot = fspl - (InCycFac + (-InCycFac / (dpeakmax_inner_inner - dpeakmin))*(dpeakmax_inner_inner - d_Rot));

							ek_Rot = ek - (InCycFac / (dpeakmax_inner_inner - dpeakmin));
						}

						if (flagdmg == 1)
							update_damage();
						if (flagdmg_Hardening_strength == 1)
							update_damage();
						if (flagdmg_Hardening == 1)
							update_damage_hardeingin();
						checkEnvelope();
						break;

					case -4:

						define_peak();
						dirtag = 0; // Positive to Negative Direction

						MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
						splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);

						Krel = fmin(pow(MtoRref / ((Kepos_Rot + Keneg_Rot) / 2), KrR) * (Kepos_Rot + Keneg_Rot) / 2, (Kepos_Rot + Keneg_Rot) / 2);
						Kun = fmin(KuR * MtoRref, (Kepos_Rot + Keneg_Rot) / 2);
						E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Kepos_Rot);  // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G

						spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot);
						F_Rot = fspl;
						ek_Rot = ek;

						if (flagdmg == 1)
							update_damage();
						if (flagdmg_Hardening_strength == 1)
							update_damage();
						if (flagdmg_Hardening == 1)
							update_damage_hardeingin();
						dpeakmax_inner = d_Rot;
						ffmax_inner = F_Rot;
						break;

					case -5:
						define_peak();
						spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot);      //Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						F_Rot = fspl;
						ek_Rot = ek;

						if (flagdmg == 1)
							update_damage();
						if (flagdmg_Hardening_strength == 1)
							update_damage();
						if (flagdmg_Hardening == 1)
							update_damage_hardeingin();
						checkEnvelope();
						break;

					case -6:

						define_peak();
						// Because of the damage model the targeted peak points may be updated so the spline parametere should get redefined again 
						if (CstateFlag != -6) {
							MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
							splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);

							Krel = fmin(pow(MtoRref / ((Kepos_Rot + Keneg_Rot) / 2), KrR) * (Kepos_Rot + Keneg_Rot) / 2, (Kepos_Rot + Keneg_Rot) / 2);
							Kun = fmin(KuR * MtoRref, (Kepos_Rot + Keneg_Rot) / 2);
							E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Kepos_Rot);  // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						}
						//------------------------------------------------------------------------------------------------------

						// -----Check if my current point is above or under the new lower limit------
						if (CstateFlag != -6) {

							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain_Rot); // NP
							Benchmark_neg_56_up = fspl;

							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain_Rot); // PN
							Benchmark_neg_56_down = fspl;
						}
						//---------------------------------------------------------------------------------------

						if ((ffmax_inner <= Benchmark_neg_56_up && ffmax_inner >= Benchmark_neg_56_down)) {

							if (CstateFlag != -6) {

								x = Cstrain_Rot;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
								fouterPN = fspl;
								ekouterPN = ek;

								x = Cstrain_Rot;
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
								fouterNP_max = fspl;
								ekouterNP = ek;

								InCycFac = (ffmax_inner - fouterPN) / (fouterNP_max - fouterPN);
							}

							x = d_Rot;
							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
							fouterPN = fspl;
							ekouterPN = ek;
							// update inner point
							x = -fabs(dpeakmax_inner - d_Rot) + dpeakmax;
							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
							finner = fouterPN + InCycFac * ((fouterNP_max - fabs(ffmax - fspl)) - fouterPN);
							ekinner = (ek - ekouterPN) / (((fouterNP_max - fabs(ffmax - fspl)) - fouterPN))*(InCycFac*((fouterNP_max - fabs(ffmax - fspl)) - fouterPN)) + ekouterPN;
							F_Rot = max(finner, fouterPN);
							ek_Rot = ekinner;
							if (F_Rot == fouterPN) {
								ek_Rot = ekouterPN;
								TstateFlag = 5;
							}
						}

						if (ffmax_inner > Benchmark_neg_56_up) {

							// Before reaching the upper limit it should be a line
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot); //NP 
							fouterNP = fspl;
							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot); //PN 
							fouterPN = fspl;
							ekouterPN = ek;

							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN
							F_Rot = ffmax_inner - ek * (dpeakmax_inner - d_Rot);
							ek_Rot = ek;

							F_Rot = fmax(F_Rot, fouterPN);
							if (F_Rot == fouterPN) {
								ek_Rot = ekouterPN;
								TstateFlag = 5;
							}

							// Check if it has reached the upper limit or not, because after it passes the upper limit the coppy-pasting process of the spline curve should get started

							if (F_Rot <= fouterNP && Cstress_Rot > fouterNP && F_Rot > fouterPN) {
								dpeakmax_inner = d_Rot;
								ffmax_inner = F_Rot;

								x = d_Rot;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
								fouterPN = fspl;
								ekouterPN = ek;

								x = d_Rot;
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
								fouterNP_max = fspl;
								ekouterNP = ek;

								InCycFac = (ffmax_inner - fouterPN) / (fouterNP_max - fouterPN);
							}

							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN
							ek_Rot = ek;

							if (F_Rot <= fouterNP && Cstress_Rot < fouterNP && F_Rot > fouterPN) {
								x = d_Rot;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
								fouterPN = fspl;
								ekouterPN = ek;
								// update inner point
								x = -fabs(dpeakmax_inner - d_Rot) + dpeakmax;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN

								finner = fouterPN + InCycFac * ((fouterNP_max - fabs(ffmax - fspl)) - fouterPN);
								ekinner = (ek - ekouterPN) / (((fouterNP_max - fabs(ffmax - fspl)) - fouterPN))*(InCycFac*((fouterNP_max - fabs(ffmax - fspl)) - fouterPN)) + ekouterPN;
								ek_Rot = ekinner;
								F_Rot = fmax(finner, fouterPN);

								if (F_Rot == fouterPN) {
									ek_Rot = ekouterPN;
									TstateFlag = 5;
								}
							}
						}

						if (ffmax_inner < Benchmark_neg_56_down) {

							if (CstateFlag != -6) {
								InCycFac = (Benchmark_neg_56_down - ffmax_inner);
							}

							x = d_Rot;
							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN 
							F_Rot = fspl - (InCycFac + (-InCycFac / (dpeakmax_inner - dpeakmin))*(dpeakmax_inner - d_Rot));
							ek_Rot = ek - (InCycFac / (dpeakmax_inner - dpeakmin));
						}

						if (flagdmg == 1)
							update_damage();
						if (flagdmg_Hardening_strength == 1)
							update_damage();
						if (flagdmg_Hardening == 1)
							update_damage_hardeingin();

						checkEnvelope();
						break;

					case -7:
						define_peak();

						// Because of the damage model the targeted peak points may be updated so the spline parametere should be redefined again 
						if (CstateFlag != -7) {
							MtoRref = (ffmax - ffmin) / (dpeakmax - dpeakmin);
							splineparam(MtoRref, dpeakmax, R_dcappos, dpeakmin, R_dcapneg);

							Krel = fmin(pow(MtoRref / ((Kepos_Rot + Keneg_Rot) / 2), KrR) * (Kepos_Rot + Keneg_Rot) / 2, (Kepos_Rot + Keneg_Rot) / 2);
							Kun = fmin(KuR * MtoRref, (Kepos_Rot + Keneg_Rot) / 2);
							E = ER * (ffmax - ffmin) * ((dpeakmax - dpeakmin) - (ffmax - ffmin) / Kepos_Rot);  // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						}
						//------------------------------------------------------------------------------------------------------

						// -----Check if the current point is above or under the new lower limit------
						if (CstateFlag != -7) {
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain_Rot); // NP
							Benchmark_neg_67_up = fspl;

							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, Cstrain_Rot); // PN
							Benchmark_neg_67_down = fspl;
						}
						//---------------------------------------------------------------------------------------
						if ((ffmin_inner_inner <= Benchmark_neg_67_up && ffmin_inner_inner >= Benchmark_neg_67_down)) {

							if (CstateFlag != -7) {

								x = Cstrain_Rot;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
								fouterPN_min = fspl;
								ekouterPN = ek;

								x = Cstrain_Rot;
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
								fouterNP = fspl;
								ekouterNP = ek;

								InCycFac = (fouterNP - ffmin_inner_inner) / (fouterNP - fouterPN_min);
							}


							x = d_Rot;
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
							fouterNP = fspl;
							ekouterNP = ek;

							x = dpeakmin + (d_Rot - dpeakmin_inner_inner);
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP

							finner = fouterNP - InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)));
							ekinner = (ek - ekouterNP) / ((fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) * (InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) + ekouterNP;

							ek_Rot = ekinner;
							F_Rot = fmin(finner, fouterNP);
							if (F_Rot == fouterNP) {
								ek_Rot = ekouterNP;
								TstateFlag = -5;
							}
						}


						if (ffmin_inner_inner < Benchmark_neg_67_down) {

							// Before reaching the lower limit it should be a line
							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot);
							fouterPN = fspl;
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, d_Rot);
							fouterNP = fspl;
							ek = ekouterNP;

							spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmax); // PN
							F_Rot = ffmin_inner_inner + ek * (d_Rot - dpeakmin_inner_inner);
							ek_Rot = ek;

							// Check if it has reached the lower limit or not, because after it passes the lower limit the coppy-pasting process of spline curve should be started

							if (F_Rot >= fouterPN && Cstress_Rot < fouterPN && F_Rot < fouterNP) {

								dpeakmin_inner_inner = d_Rot;
								ffmin_inner_inner = F_Rot;
								x = d_Rot;
								spline_curve(1, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // PN
								fouterPN_min = fspl;
								ekouterPN = ek;

								x = d_Rot;
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
								fouterNP = fspl;
								ekouterNP = ek;

								InCycFac = (fouterNP - ffmin_inner_inner) / (fouterNP - fouterPN_min);
							}

							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, dpeakmin); // NP
							ek_Rot = ek;

							if (F_Rot >= fouterPN && Cstress_Rot > fouterPN && F_Rot < fouterNP) {
								x = d_Rot;
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP
								fouterNP = fspl;
								ekouterNP = ek;

								x = dpeakmin + (d_Rot - dpeakmin_inner_inner);
								spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP

								finner = fouterNP - InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)));
								ekinner = (ek - ekouterNP) / ((fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) * (InCycFac * (fouterNP - (fouterPN_min + fabs(ffmin - fspl)))) + ekouterNP;

								ek_Rot = ekinner;
								F_Rot = fmin(finner, fouterNP);
								if (F_Rot == fouterNP) {
									ek_Rot = ekouterNP;
									TstateFlag = -5;
								}
							}
						}

						if (ffmin_inner_inner > Benchmark_neg_67_up) {

							if (CstateFlag != -7) {
								InCycFac = (ffmin_inner_inner - Benchmark_neg_67_up);
							}

							x = d_Rot;
							spline_curve(0, dpeakmin, ffmin, dpeakmax, ffmax, Krel, Kun, E, x); // NP 
							F_Rot = fspl + (InCycFac + (-InCycFac / (dpeakmax - dpeakmin_inner_inner))*(d_Rot - dpeakmin_inner_inner));
							ek_Rot = ek - (InCycFac / (dpeakmax - dpeakmin_inner_inner));
						}

						if (flagdmg == 1)
							update_damage();
						if (flagdmg_Hardening_strength == 1)
							update_damage();
						if (flagdmg_Hardening == 1)
							update_damage_hardeingin();

						checkEnvelope();
						break;
					}

				}

			}

		} //First if

		  /* ********************************************************************************************** **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  **                                                                                                **
		  **                                         Shear Material                                         **
		  **                                                                                                **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ************************************************************************************************* **
		  ** ********************************************************************************************** */
		if (MatCount == MatCount_Shr) {

			Tdu_Shr = d_Shr - Cstrain_Shr;

			if (TstateFlag_Shr == 0 || flag_entering_hardening_Flex_Rot == 1 || flag_entering_hardening_Splice_Rot == 1)
			{

				flagdmg_Shr = 0;

				if (d_Shr >= 0.0) {

					F_Shr = Kepos_Shr * d_Shr;
					ek_Shr = Kepos_Shr;
					if (d_Shr >= dpeakmax_Shr && flag_entering_hardening_Flex_Rot != 1 && flag_entering_hardening_Splice_Rot != 1) {
						Vcap_for_ER = fabs(F_Shr);

						TstateFlag_Shr = 12;
						flag_entering_hardening_Shr = 1;
						flag_entering_hardening_Shr_pos = 1;
						Calibration_Shr();
						defineBackbone_Shr();

						F_Shr = slope_pos_Shr * d_Shr + Intcpt_slope_pos_Shr;
						ek_Shr = slope_pos_Shr;
					}
				}
				else if (d_Shr < 0.0) {

					F_Shr = Keneg_Shr * d_Shr;
					ek_Shr = Keneg_Shr;
					if (d_Shr <= dpeakmin_Shr && flag_entering_hardening_Flex_Rot != 1 && flag_entering_hardening_Splice_Rot != 1) {
						Vcap_for_ER = fabs(F_Shr);

						TstateFlag_Shr = -12;
						flag_entering_hardening_Shr = 1;
						flag_entering_hardening_Shr_neg = 1;
						Calibration_Shr();
						defineBackbone_Shr();
						F_Shr = -(slope_neg_Shr * fabs(d_Shr) - Intcpt_slope_neg_Shr);
						ek_Shr = slope_neg_Shr;
					}
				}
			}

			else
			{
				if (flagFlur_Rot == 1) {
					if ((flagFlur_Rot == 1 && CflagFlur_Rot != 1) && ((CstateFlag != 2 && TstateFlag == 2) || (CstateFlag != -2 && TstateFlag == -2))) {
						d_Bench_Shr = Cstrain_Shr;
						f_Bench_Shr = Cstress_Shr;
					}
					F_Shr = Kepos_Shr * (d_Shr - d_Bench_Shr) + f_Bench_Shr;
					ek_Shr = Kepos_Shr;
					if (abs(F_Shr) > abs(f_Bench_Shr)) {
						F_Shr = f_Bench_Shr;
					}
					break;
				}

				else {
					TstateFlag_Shr = getStateFlag_Shr();
					switch (TstateFlag_Shr)
					{
					case 1:	F_Shr = slope_pos_Shr * d_Shr + Intcpt_slope_pos_Shr;
						ek_Shr = slope_pos_Shr;

						break;
					case -1: F_Shr = -(slope_neg_Shr * fabs(d_Shr) - Intcpt_slope_neg_Shr);
						ek_Shr = slope_neg_Shr;

						break;

					case 12:
						flag_entering_hardening_Shr = 1;
						flag_entering_hardening_Shr_pos = 1;
						define_peak_Shr();
						flagdmg_Shr = 0;
						flagdmg_Hardening_Shr = 1;
						flagdmg_Hardening_strength_Shr = 0;

						if (BenMark_Shr > 0) {

							ek_Shr = slope_pos_Shr;
							F_Shr = slope_pos_Shr * d_Shr + Intcpt_slope_pos_Shr;
							d12 = d_Shr;
							f12 = F_Shr;
						}
						else if (BenMark_Shr < 0) {
							ek_Shr = slope_pos_Shr;
							F_Shr = slope_pos_Shr * fabs(d_Shr) + Intcpt_slope_pos_Shr;
							d_12 = d_Shr;
							f_12 = F_Shr;
						}
						dpeakmax_Shr = d_Shr;
						if (flagdmg_Hardening_Shr == 1)
							update_damage_hardeingin_Shr();

						if (flagdmg_Hardening_strength_Shr == 1)
							update_damage_Shr();
						break;

					case -12:
						flag_entering_hardening_Shr = 1;
						flag_entering_hardening_Shr_neg = 1;
						define_peak_Shr();
						flagdmg_Shr = 0;
						flagdmg_Hardening_Shr = 1;
						flagdmg_Hardening_strength_Shr = 0;

						if (BenMark_Shr > 0) {
							ek_Shr = slope_neg_Shr;
							F_Shr = -(slope_neg_Shr * fabs(d_Shr) - Intcpt_slope_neg_Shr);
							d12neg = d_Shr;
							f12neg = F_Shr;
						}
						else if (BenMark_Shr < 0) {
							ek_Shr = slope_neg_Shr;
							F_Shr = -(slope_neg_Shr * fabs(d_Shr) - Intcpt_slope_neg_Shr);
							d_12neg = d_Shr;
							f_12neg = F_Shr;
						}

						dpeakmin_Shr = d_Shr;

						if (flagdmg_Hardening_strength_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_Shr == 1)
							update_damage_hardeingin_Shr();

						break;

					case 2:
						//---------Damage will start after entering case 2 or case -2 for the first time---------
						flagdmg_Shr = 1;
						flag_entering_hardening_Shr_pos = 0;
						flag_entering_hardening_Shr_neg = 0;
						flagdmg_Hardening_Shr = 0;
						flagdmg_Hardening_strength_Shr = 0;
						flag_Filure_Shr = 1;
						//---------------------------------------------------------------------------------------
						if (BenMark_Shr > 0) {
							ek_Shr = Kdegpos_Shr;
							F_Shr = Kdegpos_Shr * fabs(d_Shr) + Intcpt_deg_pos_Shr;
							d2 = d_Shr;
							f2 = F_Shr;
						}
						else if (BenMark_Shr < 0) {
							ek_Shr = Kdegpos_Shr;
							F_Shr = Kdegpos_Shr * fabs(d_Shr) + Intcpt_deg_pos_Shr;
							d_2 = d_Shr;
							f_2 = F_Shr;
						}
						if (flagdmg_Shr == 1)
							update_damage_Shr();
						break;

					case -2:
						//---------Damage will start after entering case 2 or case -2 for the first time---------
						flagdmg_Shr = 1;
						flag_entering_hardening_Shr_pos = 0;
						flag_entering_hardening_Shr_neg = 0;
						flagdmg_Hardening_Shr = 0;
						flagdmg_Hardening_strength_Shr = 0;
						flag_Filure_Shr = 1;
						//---------------------------------------------------------------------------------------
						if (BenMark_Shr > 0) {
							ek_Shr = Kdegneg_Shr;
							F_Shr = -(Kdegneg_Shr * fabs(d_Shr) + Intcpt_deg_neg_Shr);
							d2neg = d_Shr;
							f2neg = F_Shr;
						}
						else if (BenMark_Shr < 0) {
							ek_Shr = Kdegneg_Shr;
							F_Shr = -(Kdegneg_Shr * fabs(d_Shr) + Intcpt_deg_neg_Shr);
							d_2neg = d_Shr;
							f_2neg = F_Shr;

						}
						if (flagdmg_Shr == 1)
							update_damage_Shr();
						break;

					case 3:
						//---------Damage will start after entering case 2 or case -2 for the first time---------
						flagdmg_Shr = 0;
						flag_entering_hardening_Shr_pos = 0;
						flag_entering_hardening_Shr_neg = 0;
						flagdmg_Hardening_strength_Shr = 0;
						flagdmg_Hardening_Shr = 0;
						flag_entering_residual_Shr = 1;
						//---------------------------------------------------------------------------------------
						if (BenMark_Shr > 0) {
							ek_Shr = -0.0001;
							F_Shr = R_frespos_Shr;
							d3 = d_Shr;
							f3 = F_Shr;
						}
						else if (BenMark_Shr < 0) {
							ek_Shr = -0.0001;
							F_Shr = R_frespos_Shr;
							d_3 = d_Shr;
							f_3 = F_Shr;
						}
						break;

					case -3:
						//---------Damage will start after entering case 2 or case -2 for the first time---------
						flagdmg_Shr = 0;
						flag_entering_hardening_Shr_pos = 0;
						flag_entering_hardening_Shr_neg = 0;
						flagdmg_Hardening_strength_Shr = 0;
						flagdmg_Hardening_Shr = 0;
						flag_entering_residual_Shr = 1;
						//---------------------------------------------------------------------------------------
						if (BenMark_Shr > 0) {
							ek_Shr = -0.0001;
							F_Shr = R_fresneg_Shr;
							d3neg = d_Shr;
							f3neg = F_Shr;
						}
						else if (BenMark_Shr < 0) {
							ek_Shr = 0.0001;
							F_Shr = R_fresneg_Shr;
							d_3neg = d_Shr;
							f_3neg = F_Shr;
						}
						break;

					case 30:
						//---------Damage will start after entering case 2 or case -2 for the first time---------
						flagdmg_Shr = 0;
						flag_entering_hardening_Shr_pos = 0;
						flag_entering_hardening_Shr_neg = 0;
						flagdmg_Hardening_strength_Shr = 0;
						flagdmg_Hardening_Shr = 0;
						//---------------------------------------------------------------------------------------
						if (BenMark_Shr > 0) {
							ek_Shr = Kdegpos_Shr;
							F_Shr = Kdegpos_Shr * fabs(d_Shr) + Intcpt_res_pos_Shr;
						}
						else if (BenMark_Shr < 0) {
							ek_Shr = Kdegpos_Shr;
							F_Shr = Kdegpos_Shr * fabs(d_Shr) + Intcpt_res_pos_Shr;
						}
						break;


					case -30:
						//---------Damage will start after entering case 2 or case -2 for the first time---------
						flagdmg_Shr = 0;
						flag_entering_hardening_Shr_pos = 0;
						flag_entering_hardening_Shr_neg = 0;
						flagdmg_Hardening_strength_Shr = 0;
						flagdmg_Hardening_Shr = 0;
						//---------------------------------------------------------------------------------------
						if (BenMark_Shr > 0) {
							ek_Shr = Kdegneg_Shr;
							F_Shr = -(Kdegneg_Shr * fabs(d_Shr) + Intcpt_res_neg_Shr);
						}
						else if (BenMark_Shr < 0) {
							ek_Shr = Kdegneg_Shr;
							F_Shr = -(Kdegneg_Shr * fabs(d_Shr) + Intcpt_res_neg_Shr);
						}
						break;

					case 40:

						flagdmg_Shr = 0;
						flagdmg_Hardening_Shr = 0;
						flagdmg_Hardening_strength_Shr = 0;
						//---------------------------------------------------------------------------------------
						if (BenMark_Shr > 0) {
							ek_Shr = -0.0001;
							F_Shr = 0.0;
						}
						else if (BenMark_Shr < 0) {
							ek_Shr = -0.0001;
							F_Shr = 0.0;
						}
						break;

					case -40:

						flagdmg_Shr = 0;
						flagdmg_Hardening_strength_Shr = 0;
						flagdmg_Hardening_Shr = 0;
						//---------------------------------------------------------------------------------------
						if (BenMark_Shr > 0) {
							ek_Shr = -0.0001;
							F_Shr = 0.0;
						}
						else if (BenMark_Shr < 0) {
							ek_Shr = -0.0001;
							F_Shr = 0.0;
						}
						break;

					case 31:
						if (CstateFlag_Shr != 31)
							define_peak_Shr();
						ek_Shr = (ffmax_Shr - ffmin_Shr) / (dpeakmax_Shr - dpeakmin_Shr);
						F_Shr = ek_Shr * (d_Shr - dpeakmax_Shr) + ffmax_Shr;

						checkEnvelope_Shr();
						break;

					case 41:

						ek_Shr = 0.0001;
						F_Shr = 0.0;
						break;

					case 4:

						define_peak_Shr();
						dirtag = 1; // Positive to Negative Direction
						MtoRref_Shr = (ffmax_Shr - ffmin_Shr) / (dpeakmax_Shr - dpeakmin_Shr);   // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						splineparam_Shr(MtoRref_Shr, dpeakmax_Shr, R_dcappos_Shr, dpeakmin_Shr, R_dcapneg_Shr);

						Krel_Shr = fmin(pow(MtoRref_Shr / ((Kepos_Shr + Keneg_Shr) / 2), KrR_Shr) * (Kepos_Shr + Keneg_Shr) / 2, (Kepos_Shr + Keneg_Shr) / 2);
						Kun_Shr = fmin(KuR_Shr * MtoRref_Shr, (Kepos_Shr + Keneg_Shr) / 2);
						E_Shr = ER_Shr * (ffmax_Shr - ffmin_Shr) * ((dpeakmax_Shr - dpeakmin_Shr) - (ffmax_Shr - ffmin_Shr) / Kepos_Shr);

						spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr);      //Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						F_Shr = fspl;
						ek_Shr = ek;

						if (flagdmg_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_strength_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_Shr == 1)
							update_damage_hardeingin_Shr();
						dpeakmin_inner_Shr = d_Shr;
						ffmin_inner_Shr = F_Shr;
						break;

					case 5:

						define_peak_Shr();
						spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr);      //Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						F_Shr = fspl;
						ek_Shr = ek;
						if (flagdmg_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_strength_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_Shr == 1)
							update_damage_hardeingin_Shr();
						checkEnvelope_Shr();

						break;

					case 6:

						define_peak_Shr();
						// Because of the damage model the targeted peak points may be updated so the spline parametere should be redefined again 
						if (CstateFlag_Shr != 6) {
							MtoRref_Shr = (ffmax_Shr - ffmin_Shr) / (dpeakmax_Shr - dpeakmin_Shr);   // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
							splineparam_Shr(MtoRref_Shr, dpeakmax_Shr, R_dcappos_Shr, dpeakmin_Shr, R_dcapneg_Shr);

							Krel_Shr = fmin(pow(MtoRref_Shr / ((Kepos_Shr + Keneg_Shr) / 2), KrR_Shr) * (Kepos_Shr + Keneg_Shr) / 2, (Kepos_Shr + Keneg_Shr) / 2);
							Kun_Shr = fmin(KuR_Shr * MtoRref_Shr, (Kepos_Shr + Keneg_Shr) / 2);
							E_Shr = ER_Shr * (ffmax_Shr - ffmin_Shr) * ((dpeakmax_Shr - dpeakmin_Shr) - (ffmax_Shr - ffmin_Shr) / Kepos_Shr);
						}
						//------------------------------------------------------------------------------------------------------

						// -----I needed to check if the current point is above or under the new lower limit------
						if (CstateFlag_Shr != 6) {

							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, Cstrain_Shr); // NP
							Benchmark56_up_Shr = fspl;

							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, Cstrain_Shr); // PN
							Benchmark56_down_Shr = fspl;
						}
						//---------------------------------------------------------------------------------------

						if ((ffmin_inner_Shr <= Benchmark56_up_Shr && ffmin_inner_Shr >= Benchmark56_down_Shr)) {

							if (CstateFlag_Shr != 6) {

								x = Cstrain_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
								fouterPN_min_Shr = fspl;
								ekouterPN_Shr = ek;

								x = Cstrain_Shr;
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
								fouterNP_Shr = fspl;
								ekouterNP_Shr = ek;

								InCycFac_Shr = (fouterNP_Shr - ffmin_inner_Shr) / (fouterNP_Shr - fouterPN_min_Shr);
							}


							x = d_Shr;
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
							fouterNP_Shr = fspl;
							ekouterNP_Shr = ek;
							x = dpeakmin_Shr + (d_Shr - dpeakmin_inner_Shr);
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP

							finner_Shr = fouterNP_Shr - InCycFac_Shr * (fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)));
							ekinner_Shr = (ek - ekouterNP_Shr) / ((fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)))) * (InCycFac_Shr * (fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)))) + ekouterNP_Shr;

							F_Shr = fmin(finner_Shr, fouterNP_Shr);
							ek_Shr = ekinner_Shr;
							if (F_Shr == fouterNP_Shr) {
								ek_Shr = ekouterNP_Shr;
								TstateFlag_Shr = -5;
							}
						}

						if (ffmin_inner_Shr < Benchmark56_down_Shr) {

							// Before reaching the lower limit it should be a line
							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr);
							fouterPN_Shr = fspl;
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr);
							fouterNP_Shr = fspl;
							ekouterNP_Shr = ek;

							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, dpeakmax_Shr); // PN
							F_Shr = ffmin_inner_Shr + ek * (d_Shr - dpeakmin_inner_Shr);
							ek_Shr = ek;

							F_Shr = fmin(F_Shr, fouterNP_Shr);
							if (F_Shr == fouterNP_Shr) {
								ek = ekouterNP_Shr;
								TstateFlag_Shr = -5;
							}

							if (F_Shr >= fouterPN_Shr && Cstress_Shr < fouterPN_Shr && F_Shr < fouterNP_Shr) {
								dpeakmin_inner_Shr = d_Shr;
								ffmin_inner_Shr = F_Shr;
								x = d_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
								fouterPN_min_Shr = fspl;
								ekouterPN_Shr = ek;

								x = d_Shr;
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
								fouterNP_Shr = fspl;
								ekouterNP_Shr = ek;

								InCycFac_Shr = (fouterNP_Shr - ffmin_inner_Shr) / (fouterNP_Shr - fouterPN_min_Shr);
							}

							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, dpeakmin_Shr); // NP
							ek_Shr = ek;

							if (F_Shr >= fouterPN_Shr && Cstress_Shr > fouterPN_Shr && F_Shr < fouterNP_Shr) {
								x = d_Shr;
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
								fouterNP_Shr = fspl;
								ekouterNP_Shr = ek;
								x = dpeakmin_Shr + (d_Shr - dpeakmin_inner_Shr);
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
								finner_Shr = fouterNP_Shr - InCycFac_Shr * (fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)));
								ekinner_Shr = (ek - ekouterNP_Shr) / ((fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)))) * (InCycFac_Shr * (fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)))) + ekouterNP_Shr;
								ek_Shr = ekinner_Shr;
								F_Shr = fmin(finner_Shr, fouterNP_Shr);
								if (F_Shr == fouterNP_Shr) {
									ek_Shr = ekouterNP_Shr;
									TstateFlag_Shr = -5;
								}
							}
						}

						if (ffmin_inner_Shr > Benchmark56_up_Shr) {

							if (CstateFlag_Shr != 6) {
								InCycFac_Shr = (ffmin_inner_Shr - Benchmark56_up_Shr);
							}

							x = d_Shr;
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP 
							F_Shr = fspl + (InCycFac_Shr + (-InCycFac_Shr / (dpeakmax_Shr - dpeakmin_inner_Shr))*(d_Shr - dpeakmin_inner_Shr));
							ek_Shr = ek - (InCycFac_Shr / (dpeakmax_Shr - dpeakmin_inner_Shr));
						}

						if (flagdmg_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_strength_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_Shr == 1)
							update_damage_hardeingin_Shr();
						checkEnvelope_Shr();
						break;

					case 7:
						define_peak_Shr();

						// Because of the damage model the targeted peak points may be updated so the spline parametere should be redefined again 
						if (CstateFlag_Shr != 7) {
							MtoRref_Shr = (ffmax_Shr - ffmin_Shr) / (dpeakmax_Shr - dpeakmin_Shr);   // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
							splineparam_Shr(MtoRref_Shr, dpeakmax_Shr, R_dcappos_Shr, dpeakmin_Shr, R_dcapneg_Shr);

							Krel_Shr = fmin(pow(MtoRref_Shr / ((Kepos_Shr + Keneg_Shr) / 2), KrR_Shr) * (Kepos_Shr + Keneg_Shr) / 2, (Kepos_Shr + Keneg_Shr) / 2);
							Kun_Shr = fmin(KuR_Shr * MtoRref_Shr, (Kepos_Shr + Keneg_Shr) / 2);
							E_Shr = ER_Shr * (ffmax_Shr - ffmin_Shr) * ((dpeakmax_Shr - dpeakmin_Shr) - (ffmax_Shr - ffmin_Shr) / Kepos_Shr);
						}
						//------------------------------------------------------------------------------------------------------

						// -----Check if the current point is above or under the new lower limit------
						if (CstateFlag_Shr != 7) {
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, Cstrain_Shr); // NP
							Benchmark67_up_Shr = fspl;

							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, Cstrain_Shr); // PN
							Benchmark67_down_Shr = fspl;
						}
						//---------------------------------------------------------------------------------------

						if ((ffmax_inner_inner_Shr <= Benchmark67_up_Shr && ffmax_inner_inner_Shr >= Benchmark67_down_Shr)) {

							if (CstateFlag_Shr != 7) {

								x = Cstrain_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
								fouterPN_Shr = fspl;
								ekouterPN_Shr = ek;

								x = Cstrain_Shr;
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
								fouterNP_max_Shr = fspl;
								ekouterNP_Shr = ek;

								InCycFac_Shr = (ffmax_inner_inner_Shr - fouterPN_Shr) / (fouterNP_max_Shr - fouterPN_Shr);
							}
							x = d_Shr;
							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
							fouterPN_Shr = fspl;
							ekouterPN_Shr = ek;
							// update inner point
							x = -fabs(dpeakmax_inner_inner_Shr - d_Shr) + dpeakmax_Shr;
							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN

							finner_Shr = fouterPN_Shr + InCycFac_Shr * ((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr);
							ekinner_Shr = (ek - ekouterPN_Shr) / (((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr))*(InCycFac_Shr*((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr)) + ekouterPN_Shr;
							ek_Shr = ekinner_Shr;
							F_Shr = fmax(finner_Shr, fouterPN_Shr);
							if (F_Shr == fouterPN_Shr) {
								ek_Shr = ekouterPN_Shr;
								TstateFlag_Shr = 5;
							}
						}

						if (ffmax_inner_inner_Shr > Benchmark67_up_Shr) {

							// Before reaching the upper limit it should be a line
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr); //NP 
							fouterNP_Shr = fspl;
							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr); //PN 
							fouterPN_Shr = fspl;
							ekouterPN_Shr = ek;

							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, dpeakmax_Shr); // PN
							F_Shr = ffmax_inner_inner_Shr - ek * (dpeakmax_inner_inner_Shr - d_Shr);
							ek_Shr = ek;

							F_Shr = fmax(F_Shr, fouterPN_Shr);
							if (F_Shr == fouterPN_Shr) {
								ek_Shr = ekouterPN_Shr;
								TstateFlag_Shr = 5;
							}

							// Check if it has reached the upper limit or not, because after it passes the upper limit the coppy-pasting process of the sline curve should be started
							if (F_Shr <= fouterNP_Shr && Cstress_Shr > fouterNP_Shr && F_Shr > fouterPN_Shr) {

								dpeakmax_inner_inner_Shr = d_Shr;
								ffmax_inner_inner_Shr = F_Shr;

								x = d_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
								fouterPN_Shr = fspl;
								ekouterPN_Shr = ek;

								x = d_Shr;
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
								fouterNP_max_Shr = fspl;
								ekouterNP_Shr = ek;

								InCycFac_Shr = (ffmax_inner_inner_Shr - fouterPN_Shr) / (fouterNP_max_Shr - fouterPN_Shr);
							}

							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, dpeakmax_Shr); // PN
							ek_Shr = ek;

							if (F_Shr <= fouterNP_Shr && Cstress_Shr < fouterNP_Shr && F_Shr > fouterPN_Shr) {
								x = d_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
								fouterPN_Shr = fspl;
								ekouterPN_Shr = ek;
								// update inner point
								x = -fabs(dpeakmax_inner_inner_Shr - d_Shr) + dpeakmax_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN

								finner_Shr = fouterPN_Shr + InCycFac_Shr * ((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr);
								ekinner_Shr = (ek - ekouterPN_Shr) / (((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr))*(InCycFac_Shr * ((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr)) + ekouterPN_Shr;
								ek_Shr = ekinner_Shr;
								F_Shr = fmax(finner_Shr, fouterPN_Shr);
								if (F_Shr == fouterPN_Shr) {
									ek_Shr = ekouterPN_Shr;
									TstateFlag_Shr = 5;
								}
							}
						}

						if (ffmax_inner_inner_Shr < Benchmark67_down_Shr) {

							if (CstateFlag_Shr != 7) {
								InCycFac_Shr = (Benchmark67_down_Shr - ffmax_inner_inner_Shr);
							}

							x = d_Shr;
							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN 
							F_Shr = fspl - (InCycFac_Shr + (-InCycFac_Shr / (dpeakmax_inner_inner_Shr - dpeakmin_Shr))*(dpeakmax_inner_inner_Shr - d_Shr));

							ek_Shr = ek - (InCycFac_Shr / (dpeakmax_inner_inner_Shr - dpeakmin_Shr));
						}

						if (flagdmg_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_strength_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_Shr == 1)
							update_damage_hardeingin_Shr();

						checkEnvelope_Shr();
						break;

					case -4:

						define_peak_Shr();
						dirtag = 0; // Positive to Negative Direction
						// calculate interpolation factor
						MtoRref_Shr = (ffmax_Shr - ffmin_Shr) / (dpeakmax_Shr - dpeakmin_Shr);   // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						splineparam_Shr(MtoRref_Shr, dpeakmax_Shr, R_dcappos_Shr, dpeakmin_Shr, R_dcapneg_Shr);

						Krel_Shr = fmin(pow(MtoRref_Shr / ((Kepos_Shr + Keneg_Shr) / 2), KrR_Shr) * (Kepos_Shr + Keneg_Shr) / 2, (Kepos_Shr + Keneg_Shr) / 2);
						Kun_Shr = fmin(KuR_Shr * MtoRref_Shr, (Kepos_Shr + Keneg_Shr) / 2);
						E_Shr = ER_Shr * (ffmax_Shr - ffmin_Shr) * ((dpeakmax_Shr - dpeakmin_Shr) - (ffmax_Shr - ffmin_Shr) / Kepos_Shr);

						spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr);      //Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						F_Shr = fspl;
						ek_Shr = ek;

						if (flagdmg_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_strength_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_Shr == 1)
							update_damage_hardeingin_Shr();

						dpeakmax_inner_Shr = d_Shr;
						ffmax_inner_Shr = F_Shr;
						break;

					case -5:
						define_peak_Shr();
						spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr);      //Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
						F_Shr = fspl;
						ek_Shr = ek;
						if (flagdmg_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_strength_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_Shr == 1)
							update_damage_hardeingin_Shr();

						checkEnvelope_Shr();
						break;

					case -6:

						define_peak_Shr();

						// Because of the damage model the targeted peak points may be updated so the spline parametere should be redefined again 
						if (CstateFlag_Shr != -6) {
							MtoRref_Shr = (ffmax_Shr - ffmin_Shr) / (dpeakmax_Shr - dpeakmin_Shr);   // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
							splineparam_Shr(MtoRref_Shr, dpeakmax_Shr, R_dcappos_Shr, dpeakmin_Shr, R_dcapneg_Shr);

							Krel_Shr = fmin(pow(MtoRref_Shr / ((Kepos_Shr + Keneg_Shr) / 2), KrR_Shr) * (Kepos_Shr + Keneg_Shr) / 2, (Kepos_Shr + Keneg_Shr) / 2);
							Kun_Shr = fmin(KuR_Shr * MtoRref_Shr, (Kepos_Shr + Keneg_Shr) / 2);
							E_Shr = ER_Shr * (ffmax_Shr - ffmin_Shr) * ((dpeakmax_Shr - dpeakmin_Shr) - (ffmax_Shr - ffmin_Shr) / Kepos_Shr);
						}
						//------------------------------------------------------------------------------------------------------

						// -----Check if my current point is above or under the new lower limit------
						if (CstateFlag_Shr != -6) {

							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, Cstrain_Shr); // NP
							Benchmark_neg_56_up_Shr = fspl;

							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, Cstrain_Shr); // PN
							Benchmark_neg_56_down_Shr = fspl;
						}
						//---------------------------------------------------------------------------------------

						if ((ffmax_inner_Shr <= Benchmark_neg_56_up_Shr && ffmax_inner_Shr >= Benchmark_neg_56_down_Shr)) {
							if (CstateFlag_Shr != -6) {

								x = Cstrain_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
								fouterPN_Shr = fspl;
								ekouterPN_Shr = ek;

								x = Cstrain_Shr;
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
								fouterNP_max_Shr = fspl;
								ekouterNP_Shr = ek;

								InCycFac_Shr = (ffmax_inner_Shr - fouterPN_Shr) / (fouterNP_max_Shr - fouterPN_Shr);
							}

							x = d_Shr;
							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
							fouterPN_Shr = fspl;
							ekouterPN_Shr = ek;

							x = -fabs(dpeakmax_inner_Shr - d_Shr) + dpeakmax_Shr;
							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN

							finner_Shr = fouterPN_Shr + InCycFac_Shr * ((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr);
							ekinner_Shr = (ek - ekouterPN_Shr) / (((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr))*(InCycFac_Shr*((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr)) + ekouterPN_Shr;
							ek_Shr = ekinner_Shr;
							F_Shr = fmax(finner_Shr, fouterPN_Shr);
							if (F_Shr == fouterPN_Shr) {
								ek_Shr = ekouterPN_Shr;
								TstateFlag_Shr = 5;
							}
						}

						if (ffmax_inner_Shr > Benchmark_neg_56_up_Shr) {

							// Before reaching the upper limit it should be a line
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr); //NP 
							fouterNP_Shr = fspl;
							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr); //PN 
							fouterPN_Shr = fspl;
							ekouterPN_Shr = ek;

							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, dpeakmax_Shr); // PN
							F_Shr = ffmax_inner_Shr - ek * (dpeakmax_inner_Shr - d_Shr);
							ek_Shr = ek;

							F_Shr = fmax(F_Shr, fouterPN_Shr);
							if (F_Shr == fouterPN_Shr) {
								ek_Shr = ekouterPN_Shr;
								TstateFlag_Shr = 5;
							}
							// Check if it has reached the upper limit or not, because after it passes the upper limit the coppy-pasting process of the spline curve should be started
							if (F_Shr <= fouterNP_Shr && Cstress_Shr > fouterNP_Shr && F_Shr > fouterPN_Shr) {

								dpeakmax_inner_Shr = d_Shr;
								ffmax_inner_Shr = F_Shr;

								x = d_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
								fouterPN_Shr = fspl;
								ekouterPN_Shr = ek;

								x = d_Shr;
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
								fouterNP_max_Shr = fspl;
								ekouterNP_Shr = ek;

								InCycFac_Shr = (ffmax_inner_Shr - fouterPN_Shr) / (fouterNP_max_Shr - fouterPN_Shr);
							}

							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, dpeakmax_Shr); // PN
							ek_Shr = ek;

							if (F_Shr <= fouterNP_Shr && Cstress_Shr < fouterNP_Shr && F_Shr > fouterPN_Shr) {
								x = d_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
								fouterPN_Shr = fspl;
								ekouterPN_Shr = ek;

								x = -fabs(dpeakmax_inner_Shr - d_Shr) + dpeakmax_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN

								finner_Shr = fouterPN_Shr + InCycFac_Shr * ((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr);
								ekinner_Shr = (ek - ekouterPN_Shr) / (((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr))*(InCycFac_Shr*((fouterNP_max_Shr - fabs(ffmax_Shr - fspl)) - fouterPN_Shr)) + ekouterPN_Shr;
								ek_Shr = ekinner_Shr;
								F_Shr = fmax(finner_Shr, fouterPN_Shr);
								if (F_Shr == fouterPN_Shr) {
									ek_Shr = ekouterPN_Shr;
									TstateFlag_Shr = 5;
								}

							}
						}

						if (ffmax_inner_Shr < Benchmark_neg_56_down_Shr) {

							if (CstateFlag_Shr != -6) {
								InCycFac_Shr = (Benchmark_neg_56_down_Shr - ffmax_inner_Shr);
							}

							x = d_Shr;
							spline_curve(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN 
							F_Shr = fspl - (InCycFac_Shr + (-InCycFac_Shr / (dpeakmax_inner_Shr - dpeakmin_Shr))*(dpeakmax_inner_Shr - d_Shr));
							ek_Shr = ek - (InCycFac_Shr / (dpeakmax_inner_Shr - dpeakmin_Shr));
						}

						if (flagdmg_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_strength_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_Shr == 1)
							update_damage_hardeingin_Shr();

						checkEnvelope_Shr();
						break;

					case -7:

						define_peak_Shr();

						// Because of the damage model the targeted peak points may be updated so the spline parametere should be redefined again 
						if (CstateFlag_Shr != -7) {
							MtoRref_Shr = (ffmax_Shr - ffmin_Shr) / (dpeakmax_Shr - dpeakmin_Shr);   // Changes:dpeakmax->Tstrain12G, ffmax->Tstress12G
							splineparam_Shr(MtoRref_Shr, dpeakmax_Shr, R_dcappos_Shr, dpeakmin_Shr, R_dcapneg_Shr);

							Krel_Shr = fmin(pow(MtoRref_Shr / ((Kepos_Shr + Keneg_Shr) / 2), KrR_Shr) * (Kepos_Shr + Keneg_Shr) / 2, (Kepos_Shr + Keneg_Shr) / 2);
							Kun_Shr = fmin(KuR_Shr * MtoRref_Shr, (Kepos_Shr + Keneg_Shr) / 2);
							E_Shr = ER_Shr * (ffmax_Shr - ffmin_Shr) * ((dpeakmax_Shr - dpeakmin_Shr) - (ffmax_Shr - ffmin_Shr) / Kepos_Shr);
						}
						//------------------------------------------------------------------------------------------------------

						// -----Check if my current point is above or under the new lower limit------
						if (CstateFlag_Shr != -7) {
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, Cstrain_Shr); // NP
							Benchmark_neg_67_up_Shr = fspl;

							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, Cstrain_Shr); // PN
							Benchmark_neg_67_down_Shr = fspl;
						}
						//---------------------------------------------------------------------------------------
						if ((ffmin_inner_inner_Shr <= Benchmark_neg_67_up_Shr && ffmin_inner_inner_Shr >= Benchmark_neg_67_down_Shr)) {

							if (CstateFlag_Shr != -7) {

								x = Cstrain_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
								fouterPN_min_Shr = fspl;
								ekouterPN_Shr = ek;

								x = Cstrain_Shr;
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
								fouterNP_Shr = fspl;
								ekouterNP_Shr = ek;

								InCycFac_Shr = (fouterNP_Shr - ffmin_inner_inner_Shr) / (fouterNP_Shr - fouterPN_min_Shr);
							}

							x = d_Shr;
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
							fouterNP_Shr = fspl;
							ekouterNP_Shr = ek;

							x = dpeakmin_Shr + (d_Shr - dpeakmin_inner_inner_Shr);
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP

							finner_Shr = fouterNP_Shr - InCycFac_Shr * (fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)));
							ekinner_Shr = (ek - ekouterNP_Shr) / ((fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)))) * (InCycFac_Shr * (fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)))) + ekouterNP_Shr;

							ek_Shr = ekinner_Shr;
							F_Shr = fmin(finner_Shr, fouterNP_Shr);
							if (F_Shr == fouterNP_Shr) {
								ek_Shr = ekouterNP_Shr;
								TstateFlag_Shr = -5;
							}
						}

						if (ffmin_inner_inner_Shr < Benchmark_neg_67_down_Shr) {

							// Before reaching the lower limit it should be a line
							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr);
							fouterPN_Shr = fspl;
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, d_Shr);
							fouterNP_Shr = fspl;
							ekouterNP = ek;

							spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, dpeakmax_Shr); // PN
							F_Shr = ffmin_inner_inner_Shr + ek * (d_Shr - dpeakmin_inner_inner_Shr);
							ek_Shr = ek;

							// Check if it has reached the lower limit or not, because after it passes the lower the coppy-pasting process of the spline curve should be started

							if (F_Shr >= fouterPN_Shr && Cstress_Shr < fouterPN_Shr && F_Shr < fouterNP_Shr) {

								dpeakmin_inner_inner_Shr = d_Shr;
								ffmin_inner_inner_Shr = F_Shr;
								x = d_Shr;
								spline_curve_Shr(1, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // PN
								fouterPN_min_Shr = fspl;
								ekouterPN_Shr = ek;

								x = d_Shr;
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
								fouterNP_Shr = fspl;
								ekouterNP_Shr = ek;

								InCycFac_Shr = (fouterNP_Shr - ffmin_inner_inner_Shr) / (fouterNP_Shr - fouterPN_min_Shr);

							}

							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, dpeakmin_Shr); // NP
							ek_Shr = ek;

							if (F_Shr >= fouterPN_Shr && Cstress_Shr > fouterPN_Shr && F_Shr < fouterNP_Shr) {
								x = d_Shr;
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP
								fouterNP_Shr = fspl;
								ekouterNP_Shr = ek;

								x = dpeakmin_Shr + (d_Shr - dpeakmin_inner_inner_Shr);
								spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP

								finner_Shr = fouterNP_Shr - InCycFac_Shr * (fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)));
								ekinner_Shr = (ek - ekouterNP_Shr) / ((fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)))) * (InCycFac_Shr * (fouterNP_Shr - (fouterPN_min_Shr + fabs(ffmin_Shr - fspl)))) + ekouterNP_Shr;

								ek_Shr = ekinner_Shr;
								F_Shr = fmin(finner_Shr, fouterNP_Shr);
								if (F_Shr == fouterNP_Shr) {
									ek_Shr = ekouterNP_Shr;
									TstateFlag_Shr = -5;
								}

							}
						}

						if (ffmin_inner_inner_Shr > Benchmark_neg_67_up_Shr) {

							if (CstateFlag_Shr != -7) {
								InCycFac_Shr = (ffmin_inner_inner_Shr - Benchmark_neg_67_up_Shr);
							}

							x = d_Shr;
							spline_curve_Shr(0, dpeakmin_Shr, ffmin_Shr, dpeakmax_Shr, ffmax_Shr, Krel_Shr, Kun_Shr, E_Shr, x); // NP 

							F_Shr = fspl + (InCycFac_Shr + (-InCycFac_Shr / (dpeakmax_Shr - dpeakmin_inner_inner_Shr))*(d_Shr - dpeakmin_inner_inner_Shr));
							ek_Shr = ek - (InCycFac_Shr / (dpeakmax_Shr - dpeakmin_inner_inner_Shr));
						}

						if (flagdmg_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_strength_Shr == 1)
							update_damage_Shr();
						if (flagdmg_Hardening_Shr == 1)
							update_damage_hardeingin_Shr();

						checkEnvelope_Shr();
						break;

					}
				}

			}

		}

		/* ********************************************************************************************** **
		************************************************************************************************* **
		**                                                                                                **
		**                                         Axial Material                                         **
		**                                                                                                **
		************************************************************************************************* **
		** ********************************************************************************************** */
		if (MatCount == 1) {
			if (d_Axial > 0.0) {
				F_Axial = d_Axial * Kepos_Axil;
				Stress_vector(1) = Kepos_Axil * d_Axial;
				Tang_vector(1) = Kepos_Axil;
			}
			else {
				F_Axial = d_Axial * Keneg_Axil;
				Stress_vector(1) = Keneg_Axil * d_Axial;
				Tang_vector(1) = Keneg_Axil;
			}

			Damp_vector(1) = Kepos_Axil;
			Needed_Inpt_to_get_failure_mode();
		}

		if (MatCount == MatCount_Rot) {
			Stress_vector(2) = F_Rot;
			Tang_vector(2) = ek_Rot;
			Damp_vector(2) = Cek_Rot;
		}
		if (MatCount == MatCount_Shr) {
			Stress_vector(0) = F_Shr;
			Tang_vector(0) = ek_Shr;
			Damp_vector(0) = ek_Shr;
		}

	} // End of the main loop

	return 0;
}

double GMG_CMAC2D::getStrain(void)
{
	d = getStrain_Multi();
	return d;
}

double GMG_CMAC2D::getStress(void)
{
	f = getStress_Multi();
	return f;
}

double GMG_CMAC2D::getTangent(void)
{
	ek = getTangent_Multi();
	return ek;
}

double GMG_CMAC2D::getInitialTangent(void)
{
	double InitialTangent = getInitialTangent_Multi();
	return (InitialTangent);
}

double GMG_CMAC2D::getDampTangent(void)
{
	double DampTangent = getDampTangent_Multi();
	return 0;
}

double GMG_CMAC2D::getStrainRate(void)
{
	return 0;
}

int GMG_CMAC2D::commitState(void)
{

	if (TstateFlag == 12) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmax = d_Rot;
		ffmax = F_Rot;
	}

	if (TstateFlag == -12) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmin = d_Rot;
		ffmin = F_Rot;
	}

	if (TstateFlag == 2) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmax = d_Rot;
		ffmax = F_Rot;
	}

	if (TstateFlag == -2) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmin = d_Rot;
		ffmin = F_Rot;
	}

	if (TstateFlag == 5) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmin_inner = d_Rot;
		ffmin_inner = F_Rot;
	}

	if (TstateFlag == -5) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmax_inner = d_Rot;
		ffmax_inner = F_Rot;
	}

	if (TstateFlag == 6) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmax_inner_inner = d_Rot;
		ffmax_inner_inner = F_Rot;
	}

	if (TstateFlag == -6) {

		//-------------------------Case -7 Calculation---------------------
		dpeakmin_inner_inner = d_Rot;
		ffmin_inner_inner = F_Rot;
	}

	if (TstateFlag == 7) {

		//-------------------------Case -7 Calculation---------------------
		dpeakmin_inner_inner = d_Rot;
		ffmin_inner_inner = F_Rot;
	}

	if (TstateFlag == -7) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmax_inner_inner = d_Rot;
		ffmax_inner_inner = F_Rot;
	}
	//-----------------------------------------End of part 1-------------------------------------------

	   /* ********************************************************************************************** **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		**                                                                                                **
		**                                         Shear Material                                         **
		**                                                                                                **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		************************************************************************************************* **
		** ********************************************************************************************** */
	if (TstateFlag_Shr == 12) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmax_Shr = d_Shr;
		ffmax_Shr = F_Shr;
	}

	if (TstateFlag_Shr == -12) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmin_Shr = d_Shr;
		ffmin_Shr = F_Shr;
	}

	if (TstateFlag_Shr == 2) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmax_Shr = d_Shr;
		ffmax_Shr = F_Shr;
	}

	if (TstateFlag_Shr == -2) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmin_Shr = d_Shr;
		ffmin_Shr = F_Shr;
	}


	if (TstateFlag_Shr == 5) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmin_inner_Shr = d_Shr;
		ffmin_inner_Shr = F_Shr;
	}

	if (TstateFlag_Shr == -5) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmax_inner_Shr = d_Shr;
		ffmax_inner_Shr = F_Shr;
	}


	if (TstateFlag_Shr == 6) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmax_inner_inner_Shr = d_Shr;
		ffmax_inner_inner_Shr = F_Shr;
	}

	if (TstateFlag_Shr == -6) {

		//-------------------------Case -7 Calculation---------------------
		dpeakmin_inner_inner_Shr = d_Shr;
		ffmin_inner_inner_Shr = F_Shr;
	}

	if (TstateFlag_Shr == 7) {

		//-------------------------Case -7 Calculation---------------------
		dpeakmin_inner_inner_Shr = d_Shr;
		ffmin_inner_inner_Shr = F_Shr;
	}

	if (TstateFlag_Shr == -7) {

		//-------------------------Case 7 Calculation---------------------
		dpeakmax_inner_inner_Shr = d_Shr;
		ffmax_inner_inner_Shr = F_Shr;
	}


	//*********************************************************************************************

	Cstrain = d;
	Cstrain_Rot = d_Rot;
	Cstrain_Shr = d_Shr;
	Cstrain_Axial = d_Axial;

	Cstress = f;
	Cstress_Rot = F_Rot;
	Cstress_Shr = F_Shr;
	Cstress_Axial = F_Axial;

	if (d_Axial > 0) {
		if (flag_entering_hardening_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {

			if ((TstateFlag == 12 && CstateFlag == 0) || (TstateFlag == -12 && CstateFlag == 0))
			{
				for (int i = 0; i <= number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_average;
				}

				for (int i = 0; i < number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_Vector(i + 1);
				}
				if (Kepos_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 2) < 0) {
					if (TstateFlag != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = max(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) - 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Kepos_Axil * d_Axial;
				}

				else {
					if (TstateFlag != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = min(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) + 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Kepos_Axil * d_Axial;
				}
				P_Axial_sum = 0.0;
				for (int i = 0; i <= number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_sum += P_Axial_Vector(i);
				}

			}
			else
			{
				for (int i = 0; i < number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_Vector(i + 1);
				}
				if (Kepos_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 2) < 0) {
					if (TstateFlag != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = max(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) - 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Kepos_Axil * d_Axial;
				}

				else {
					if (TstateFlag != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = min(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) + 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Kepos_Axil * d_Axial;
				}

				P_Axial_sum = 0.0;
				if (TstateFlag != 0) {
					for (int i = 0; i < number_of_averaged_elements_Hardening; i++) {
						P_Axial_sum += P_Axial_Vector(i);
					}
				}
				else {
					for (int i = number_of_averaged_elements_Hardening - number_of_averaged_elements_Elastic; i < number_of_averaged_elements_Hardening; i++) {
						P_Axial_sum += P_Axial_Vector(i);
					}
				}
			}
		}



		else if (flag_Filure_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {

			if ((TstateFlag_Shr == 12 && CstateFlag == 0) || (TstateFlag_Shr == -12 && CstateFlag == 0))
			{
				for (int i = 0; i <= number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_average;
				}

				for (int i = 0; i < number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_Vector(i + 1);
				}
				if (Kepos_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 2) < 0) {
					if (TstateFlag_Shr != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = max(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) - 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Kepos_Axil * d_Axial;
				}

				else {
					if (TstateFlag_Shr != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = min(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) + 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Kepos_Axil * d_Axial;
				}
				P_Axial_sum = 0.0;
				for (int i = 0; i <= number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_sum += P_Axial_Vector(i);
				}

			}
			else
			{
				for (int i = 0; i < number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_Vector(i + 1);
				}
				if (Kepos_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 2) < 0) {
					if (TstateFlag_Shr != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = max(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) - 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Kepos_Axil * d_Axial;
				}

				else {
					if (TstateFlag_Shr != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = min(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) + 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Kepos_Axil * d_Axial;
				}

				P_Axial_sum = 0.0;
				if (TstateFlag_Shr != 0) {
					for (int i = 0; i < number_of_averaged_elements_Hardening; i++) {
						P_Axial_sum += P_Axial_Vector(i);
					}
				}
				else {
					for (int i = number_of_averaged_elements_Hardening - number_of_averaged_elements_Elastic; i < number_of_averaged_elements_Hardening; i++) {
						P_Axial_sum += P_Axial_Vector(i);
					}
				}
			}
		}
	}


	else if (d_Axial <= 0) {
		if (flag_entering_hardening_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {

			if ((TstateFlag == 12 && CstateFlag == 0) || (TstateFlag == -12 && CstateFlag == 0))
			{
				for (int i = 0; i <= number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_average;
				}

				for (int i = 0; i < number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_Vector(i + 1);
				}
				if (Keneg_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 2) < 0) {
					if (TstateFlag != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = max(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) - 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Keneg_Axil * d_Axial;
				}

				else {
					if (TstateFlag != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = min(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) + 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Keneg_Axil * d_Axial;
				}
				P_Axial_sum = 0.0;
				for (int i = 0; i <= number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_sum += P_Axial_Vector(i);
				}

			}
			else
			{
				for (int i = 0; i < number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_Vector(i + 1);
				}
				if (Keneg_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 2) < 0) {
					if (TstateFlag != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = max(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) - 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Keneg_Axil * d_Axial;
				}

				else {
					if (TstateFlag != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = min(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) + 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Keneg_Axil * d_Axial;
				}

				P_Axial_sum = 0.0;
				if (TstateFlag != 0) {
					for (int i = 0; i < number_of_averaged_elements_Hardening; i++) {
						P_Axial_sum += P_Axial_Vector(i);
					}
				}
				else {
					for (int i = number_of_averaged_elements_Hardening - number_of_averaged_elements_Elastic; i < number_of_averaged_elements_Hardening; i++) {
						P_Axial_sum += P_Axial_Vector(i);
					}
				}
			}
		}



		else if (flag_Filure_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {

			if ((TstateFlag_Shr == 12 && CstateFlag == 0) || (TstateFlag_Shr == -12 && CstateFlag == 0))
			{
				for (int i = 0; i <= number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_average;
				}

				for (int i = 0; i < number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_Vector(i + 1);
				}
				if (Keneg_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 2) < 0) {
					if (TstateFlag_Shr != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = max(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) - 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Keneg_Axil * d_Axial;
				}

				else {
					if (TstateFlag_Shr != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = min(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) + 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Keneg_Axil * d_Axial;
				}
				P_Axial_sum = 0.0;
				for (int i = 0; i <= number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_sum += P_Axial_Vector(i);
				}

			}
			else
			{
				for (int i = 0; i < number_of_averaged_elements_Hardening - 1; i++) {
					P_Axial_Vector(i) = P_Axial_Vector(i + 1);
				}
				if (Keneg_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 2) < 0) {
					if (TstateFlag_Shr != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = max(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) - 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Keneg_Axil * d_Axial;
				}

				else {
					if (TstateFlag_Shr != 0)
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = min(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 2) + 0.01*fc*B*H);
					else
						P_Axial_Vector(number_of_averaged_elements_Hardening - 1) = Keneg_Axil * d_Axial;
				}

				P_Axial_sum = 0.0;
				if (TstateFlag_Shr != 0) {
					for (int i = 0; i < number_of_averaged_elements_Hardening; i++) {
						P_Axial_sum += P_Axial_Vector(i);
					}
				}
				else {
					for (int i = number_of_averaged_elements_Hardening - number_of_averaged_elements_Elastic; i < number_of_averaged_elements_Hardening; i++) {
						P_Axial_sum += P_Axial_Vector(i);
					}
				}
			}
		}
	}

	// Ratio of VyE/Vcol at flexural yielding for defferentiating between flexure and flexure-shear
	if (flag_entering_hardening_Flex_Rot == 1 && Cflag_entering_hardening_Flex_Rot != 1) {
		VyE_to_Vcol_Ratio = fypos_Rot / La / Vcol;
	}

	if (flag_entering_hardening_Shr == 1 && Cflag_entering_hardening_Shr != 1) {
		Vcap_for_ER = fabs(F_Shr);
	}

	// Damage outputs
	if ((flag_entering_hardening_Flex_Rot == 1 || flag_entering_hardening_Shr == 1 || flag_entering_hardening_Splice_Rot == 1) && (flag_Flx_Filure_Rot != 1 && flag_Filure_Shr != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1)) {
		DamageOutput_Hardening();
	}

	if ((flag_entering_hardening_Flex_Rot == 1 && Cflag_entering_hardening_Flex_Rot != 1) || (flag_entering_hardening_Splice_Rot == 1 && Cflag_entering_hardening_Splice_Rot != 1) || (flag_entering_hardening_Shr == 1 && Cflag_entering_hardening_Shr != 1) || (flag_Flx_Filure_Rot == 1 && Cflag_Flx_Filure_Rot != 1) || (flag_Filure_Shr == 1 && Cflag_Filure_Shr != 1) || (flag_FlexShr_Filure_Rot == 1 && Cflag_FlexShr_Filure_Rot != 1) || (flag_Splice_Filure_Rot == 1 && Cflag_Splice_Filure_Rot != 1) || (flag_FlexSplice_Filure_Rot == 1 && Cflag_FlexSplice_Filure_Rot != 1)) {
		DamageOutput_FailureType();
	}

	if ((flag_Flx_Filure_Rot == 1 || flag_Filure_Shr == 1 || flag_FlexShr_Filure_Rot == 1 || flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) && flag_entering_residual_Flex_Rot != 1 && flag_entering_residual_Shr != 1) {
		DamageOutput_PostCapping();
	}

	if ((flag_entering_residual_Flex_Rot == 1 && Cflag_entering_residual_Flex_Rot != 1) || (flag_entering_residual_Shr == 1 && Cflag_entering_residual_Shr != 1)) {
		DamageOutput_PostResidual();
	}

	Cek = ek;
	Cek_Rot = ek_Rot;
	Cek_Shr = ek_Shr;
	Cek_Axi = ek_Axi;

	CInCycFac = InCycFac;
	CInCycFac_Shr = InCycFac_Shr;

	CVyE_to_Vcol_Ratio = VyE_to_Vcol_Ratio;

	Cflagdmg = flagdmg;
	Cflagdmg_Shr = flagdmg_Shr;

	Cflagdmg_Hardening_strength = flagdmg_Hardening_strength;
	Cflagdmg_Hardening_strength_Shr = flagdmg_Hardening_strength_Shr;

	Cflagdmg_Hardening = flagdmg_Hardening;
	Cflagdmg_Hardening_Shr = flagdmg_Hardening_Shr;

	Cflag_entering_hardening_Flex_Rot = flag_entering_hardening_Flex_Rot;
	Cflag_entering_residual_Flex_Rot = flag_entering_residual_Flex_Rot;
	Cflag_entering_residual_Shr = flag_entering_residual_Shr;
	Cflag_entering_hardening_Flex_pos_Rot = flag_entering_hardening_Flex_pos_Rot;
	Cflag_entering_hardening_Flex_neg_Rot = flag_entering_hardening_Flex_neg_Rot;

	Cflag_entering_hardening_Splice_pos_Rot = flag_entering_hardening_Splice_pos_Rot;
	Cflag_entering_hardening_Splice_neg_Rot = flag_entering_hardening_Splice_neg_Rot;

	Cflag_entering_hardening_Shr_pos = flag_entering_hardening_Shr_pos;
	Cflag_entering_hardening_Shr_neg = flag_entering_hardening_Shr_neg;

	Cflag_entering_hardening_Splice_Rot = flag_entering_hardening_Splice_Rot;
	Cflag_entering_hardening_Shr = flag_entering_hardening_Shr;

	Cflag_FlexShr_Filure_Rot = flag_FlexShr_Filure_Rot;
	Cflag_Splice_Filure_Rot = flag_Splice_Filure_Rot;
	Cflag_FlexSplice_Filure_Rot = flag_FlexSplice_Filure_Rot;
	Cflag_Flx_Filure_Rot = flag_Flx_Filure_Rot;
	Cflag_Filure_Shr = flag_Filure_Shr;

	CflagFlur_Rot = flagFlur_Rot;
	CflagFlur_Shr = flagFlur_Shr;

	Cf_Bench_Rot = f_Bench_Rot;
	Cf_Bench_Shr = f_Bench_Shr;

	Cd_Bench_Rot = d_Bench_Rot;
	Cd_Bench_Shr = d_Bench_Shr;

	Cdpeakmax = dpeakmax;
	Cdpeakmax_Shr = dpeakmax_Shr;

	Cdpeakmax_bench = dpeakmax_bench;
	Cdpeakmax_bench_Shr = dpeakmax_bench_Shr;

	Cdpeakmax_inner = dpeakmax_inner;
	Cdpeakmax_inner_Shr = dpeakmax_inner_Shr;

	Cdpeakmax_inner_inner = dpeakmax_inner_inner;
	Cdpeakmax_inner_inner_Shr = dpeakmax_inner_inner_Shr;

	Cdpeakmin = dpeakmin;
	Cdpeakmin_Shr = dpeakmin_Shr;

	Cdpeakmin_bench = dpeakmin_bench;
	Cdpeakmin_bench_Shr = dpeakmin_bench_Shr;

	Cdpeakmin_inner = dpeakmin_inner;
	Cdpeakmin_inner_Shr = dpeakmin_inner_Shr;

	Cdpeakmin_inner_inner = dpeakmin_inner_inner;
	Cdpeakmin_inner_inner_Shr = dpeakmin_inner_inner_Shr;

	Cffmax = ffmax;
	Cffmax_Shr = ffmax_Shr;

	Cffmax_inner = ffmax_inner;
	Cffmax_inner_Shr = ffmax_inner_Shr;

	Cffmax_inner_inner = ffmax_inner_inner;
	Cffmax_inner_inner_Shr = ffmax_inner_inner_Shr;

	CfouterNP_max = fouterNP_max;
	CfouterNP_max_Shr = fouterNP_max_Shr;

	Cffmin = ffmin;
	Cffmin_Shr = ffmin_Shr;

	Cffmin_inner = ffmin_inner;
	Cffmin_inner_Shr = ffmin_inner_Shr;

	Cffmin_inner_inner = ffmin_inner_inner;
	Cffmin_inner_inner_Shr = ffmin_inner_inner_Shr;

	CfouterPN_min = fouterPN_min;
	CfouterPN_min_Shr = fouterPN_min_Shr;

	Cspos = spos;
	Cspos_Shr = spos_Shr;

	Csneg = sneg;
	Csneg_Shr = sneg_Shr;

	CdmgSpos = dmgSpos;
	CdmgSpos_Shr = dmgSpos_Shr;

	CdmgSneg = dmgSneg;
	CdmgSneg_Shr = dmgSneg_Shr;

	CEt = Et;
	CEt_Shr = Et_Shr;

	CEc = Ec;
	CEc_Shr = Ec_Shr;

	CEposnorm = Eposnorm;
	CEposnorm_Shr = Eposnorm_Shr;

	CEnegnorm = Enegnorm;
	CEnegnorm_Shr = Enegnorm_Shr;

	CMtoRref = MtoRref;
	CMtoRref_Shr = MtoRref_Shr;

	CKrR = KrR;
	CKrR_Shr = KrR_Shr;

	CKrel = Krel;
	CKrel_Shr = Krel_Shr;

	CKuR = KuR;
	CKuR_Shr = KuR_Shr;

	CKun = Kun;
	CKun_Shr = Kun_Shr;

	CER = ER;
	CER_Shr = ER_Shr;

	CE = E;
	CE_Shr = E_Shr;

	CTdu = Tdu;
	CTdu_Rot = Tdu_Rot;
	CTdu_Shr = Tdu_Shr;

	CstrainRate = TstrainRate;
	CstrainMax = TstrainMax;
	CstrainMin = TstrainMin;

	CstateFlag = TstateFlag;
	CstateFlag_Shr = TstateFlag_Shr;

	CVcol = Vcol;

	//************************************************
	Cd12 = d12;
	Cf12 = f12;
	Cd_12 = d_12;
	Cd12neg = d12neg;
	Cf12neg = f12neg;
	Cd_12neg = d_12neg;
	Cf_12neg = f_12neg;
	Cd2 = d2;
	Cf2 = f2;
	Cf_12 = f_12;
	Cd_2 = d_2;
	Cf_2 = f_2;
	Cd2neg = d2neg;
	Cf2neg = f2neg;
	Cd_2neg = d_2neg;
	Cf_2neg = f_2neg;
	Cd3 = d3;
	Cf3 = f3;
	Cd_3 = d_3;
	Cf_3 = f_3;
	Cd3neg = d3neg;
	Cf3neg = f3neg;
	Cd_3neg = d_3neg;
	Cf_3neg = f_3neg;
	//************************************************

	Cratio = ratio;
	Cratio_Shr = ratio_Shr;

	CdmgCounter = dmgCounter;
	CdmgCounter_Shr = dmgCounter_Shr;

	CR_dcapneg = R_dcapneg;
	CR_dcapneg_Shr = R_dcapneg_Shr;

	CR_dcappos = R_dcappos;
	CR_dcappos_Shr = R_dcappos_Shr;

	CIntcpt_slope_neg = Intcpt_slope_neg;
	CIntcpt_slope_neg_Shr = Intcpt_slope_neg_Shr;

	CR_delUneg = R_delUneg;

	CR_fresneg = R_fresneg;
	CR_fresneg_Shr = R_fresneg_Shr;

	Cslope_neg = slope_neg;
	Cslope_neg_Shr = slope_neg_Shr;

	Cdpneg = dpneg;
	Cdpneg_Shr = dpneg_Shr;

	CIntcpt_slope_pos = Intcpt_slope_pos;
	CIntcpt_slope_pos_Shr = Intcpt_slope_pos_Shr;

	CR_delUpos = R_delUpos;
	CR_delUpos_Shr = R_delUpos_Shr;
	CR_delUneg_Shr = R_delUneg_Shr;

	CR_frespos = R_frespos;
	CR_frespos_Shr = R_frespos_Shr;

	Cslope_pos = slope_pos;
	Cslope_pos_Shr = slope_pos_Shr;

	Cdppos = dppos;
	Cdppos_Shr = dppos_Shr;

	CR_Kdegneg = R_Kdegneg;
	CR_Kdegneg_Shr = R_Kdegneg_Shr;

	CR_fcapneg = R_fcapneg;
	CR_fcapneg_Shr = R_fcapneg_Shr;

	CR_dyneg = R_dyneg;
	CR_dyneg_Shr = R_dyneg_Shr;

	CR_fyneg = R_fyneg;
	CR_fyneg_Shr = R_fyneg_Shr;

	CR_Kdegpos = R_Kdegpos;
	CR_Kdegpos_Shr = R_Kdegpos_Shr;

	CR_fcappos = R_fcappos;
	CR_fcappos_Shr = R_fcappos_Shr;

	CR_fypos = R_fypos;
	CR_fypos_Shr = R_fypos_Shr;

	CR_dypos = R_dypos;
	CR_dypos_Shr = R_dypos_Shr;

	CR_dresneg = R_dresneg;
	CR_dresneg_Shr = R_dresneg_Shr;

	CR_drespos = R_drespos;
	CR_drespos_Shr = R_drespos_Shr;

	CIntcpt_deg_pos = Intcpt_deg_pos;
	CIntcpt_deg_pos_Shr = Intcpt_deg_pos_Shr;

	CIntcpt_Xaxis_pos_Shr = Intcpt_Xaxis_pos_Shr;
	CIntcpt_Xaxis_neg_Shr = Intcpt_Xaxis_neg_Shr;

	CIntcpt_Xaxis_pos = Intcpt_Xaxis_pos;
	CIntcpt_Xaxis_neg = Intcpt_Xaxis_neg;

	CIntcpt_deg_neg = Intcpt_deg_neg;
	CIntcpt_deg_neg_Shr = Intcpt_deg_neg_Shr;

	CIntcpt_res_pos = Intcpt_res_pos;
	CIntcpt_res_pos_Shr = Intcpt_res_pos_Shr;

	CIntcpt_res_neg = Intcpt_res_neg;
	CIntcpt_res_neg_Shr = Intcpt_res_neg_Shr;


	CInCycFac_6 = InCycFac_6;
	CInCycFac_6_Shr = InCycFac_6_Shr;

	CInCycFac_7 = InCycFac_7;
	CInCycFac_7_Shr = InCycFac_7_Shr;

	CInCycFac_neg_6 = InCycFac_neg_6;
	CInCycFac_neg_6_Shr = InCycFac_neg_6_Shr;

	CInCycFac_neg_7 = InCycFac_neg_7;
	CInCycFac_neg_7_Shr = InCycFac_neg_7_Shr;

	CBenchmark56_up = Benchmark56_up;
	CBenchmark56_up_Shr = Benchmark56_up_Shr;

	CBenchmark56_down = Benchmark56_down;
	CBenchmark56_down_Shr = Benchmark56_down_Shr;


	CBenchmark_neg_56_up = Benchmark_neg_56_up;
	CBenchmark_neg_56_up_Shr = Benchmark_neg_56_up_Shr;

	CBenchmark_neg_56_down = Benchmark_neg_56_down;
	CBenchmark_neg_56_down_Shr = Benchmark_neg_56_down_Shr;

	CBenchmark67_up = Benchmark67_up;
	CBenchmark67_up_Shr = Benchmark67_up_Shr;

	CBenchmark67_down = Benchmark67_down;
	CBenchmark67_down_Shr = Benchmark67_down_Shr;

	CBenchmark_neg_67_up = Benchmark_neg_67_up;
	CBenchmark_neg_67_up_Shr = Benchmark_neg_67_up_Shr;

	CBenchmark_neg_67_down = Benchmark_neg_67_down;
	CBenchmark_neg_67_down_Shr = Benchmark_neg_67_down_Shr;

	Cspos_pos = spos_pos;
	Cspos_pos_Shr = spos_pos_Shr;

	Csneg_pos = sneg_pos;
	Csneg_pos_Shr = sneg_pos_Shr;

	CEt_pos = Et_pos;
	CEt_pos_Shr = Et_pos_Shr;

	CEc_pos = Ec_pos;
	CEc_pos_Shr = Ec_pos_Shr;

	Cspos_neg = spos_neg;
	Cspos_neg_Shr = spos_neg_Shr;

	Csneg_neg = sneg_neg;
	Csneg_neg_Shr = sneg_neg_Shr;

	CEt_neg = Et_neg;
	CEt_neg_Shr = Et_neg_Shr;

	CEc_neg = Ec_neg;
	CEc_neg_Shr = Ec_neg_Shr;

	CEpos_pos = Epos_pos;
	CEpos_pos_Shr = Epos_pos_Shr;

	CEneg_neg = Eneg_neg;
	CEneg_neg_Shr = Eneg_neg_Shr;

	CT_area = T_area;

	CEt_hard = Et_hard;
	CEc_hard = Ec_hard;
	Cspos_hard = spos_hard;
	Csneg_hard = sneg_hard;
	Cspos_pos_hard = spos_pos_hard;
	Csneg_pos_hard = sneg_pos_hard;
	CEt_pos_hard = Et_pos_hard;
	CEc_pos_hard = Ec_pos_hard;
	Cspos_neg_hard = spos_neg_hard;
	Csneg_neg_hard = sneg_neg_hard;
	CEt_neg_hard = Et_neg_hard;
	CEc_neg_hard = Ec_neg_hard;
	CEpos_pos_hard = Epos_pos_hard;
	CEneg_neg_hard = Eneg_neg_hard;
	CEpos_hard = Epos_hard;
	CEneg_hard = Eneg_hard;
	CEposnorm_hard = Eposnorm_hard;
	CEnegnorm_hard = Enegnorm_hard;
	Cdelta_pos_hard = delta_pos_hard;
	Cdelta_pos_max_hard = delta_pos_max_hard;
	Cdelta_neg_hard = delta_neg_hard;
	Cdelta_neg_max_hard = delta_neg_max_hard;
	Calpha_pos = alpha_pos;
	Calpha_neg = alpha_neg;

	//Hardening damage model Shear
	CEt_hard_Shr = Et_hard_Shr;
	CEc_hard_Shr = Ec_hard_Shr;
	Cspos_hard_Shr = spos_hard_Shr;
	Csneg_hard_Shr = sneg_hard_Shr;
	Cspos_pos_hard_Shr = spos_pos_hard_Shr;
	Csneg_pos_hard_Shr = sneg_pos_hard_Shr;
	CEt_pos_hard_Shr = Et_pos_hard_Shr;
	CEc_pos_hard_Shr = Ec_pos_hard_Shr;
	Cspos_neg_hard_Shr = spos_neg_hard_Shr;
	Csneg_neg_hard_Shr = sneg_neg_hard_Shr;
	CEt_neg_hard_Shr = Et_neg_hard_Shr;
	CEc_neg_hard_Shr = Ec_neg_hard_Shr;
	CEpos_pos_hard_Shr = Epos_pos_hard_Shr;
	CEneg_neg_hard_Shr = Eneg_neg_hard_Shr;
	CEpos_hard_Shr = Epos_hard_Shr;
	CEneg_hard_Shr = Eneg_hard_Shr;
	CEposnorm_hard_Shr = Eposnorm_hard_Shr;
	CEnegnorm_hard_Shr = Enegnorm_hard_Shr;
	Cdelta_pos_hard_Shr = delta_pos_hard_Shr;
	Cdelta_pos_max_hard_Shr = delta_pos_max_hard_Shr;
	Cdelta_neg_hard_Shr = delta_neg_hard_Shr;
	Cdelta_neg_max_hard_Shr = delta_neg_max_hard_Shr;
	Calpha_pos_Shr = alpha_pos_Shr;
	Calpha_neg_Shr = alpha_neg_Shr;

	return 0;
}

int GMG_CMAC2D::revertToLastCommit(void)
{
	d = Cstrain;
	d_Rot = Cstrain_Rot;
	d_Shr = Cstrain_Shr;
	d_Axial = Cstrain_Axial;

	f = Cstress;
	F_Rot = Cstress_Rot;
	F_Shr = Cstress_Shr;
	F_Axial = Cstress_Axial;

	ek = Cek;
	ek_Rot = Cek_Rot;
	ek_Shr = Cek_Shr;
	ek_Axi = Cek_Axi;

	InCycFac = CInCycFac;
	InCycFac_Shr = CInCycFac_Shr;

	VyE_to_Vcol_Ratio = CVyE_to_Vcol_Ratio;
	Vcap_for_ER = CVcap_for_ER;

	flagdmg = Cflagdmg;
	flagdmg_Shr = Cflagdmg_Shr;

	flagdmg_Hardening_strength = Cflagdmg_Hardening_strength;
	flagdmg_Hardening_strength_Shr = Cflagdmg_Hardening_strength_Shr;

	flagdmg_Hardening = Cflagdmg_Hardening;
	flagdmg_Hardening_Shr = Cflagdmg_Hardening_Shr;

	flag_entering_hardening_Flex_Rot = Cflag_entering_hardening_Flex_Rot;
	flag_entering_residual_Flex_Rot = Cflag_entering_residual_Flex_Rot;
	flag_entering_residual_Shr = Cflag_entering_residual_Shr;
	flag_entering_hardening_Flex_pos_Rot = Cflag_entering_hardening_Flex_pos_Rot;
	flag_entering_hardening_Flex_neg_Rot = Cflag_entering_hardening_Flex_neg_Rot;

	flag_entering_hardening_Splice_pos_Rot = Cflag_entering_hardening_Splice_pos_Rot;
	flag_entering_hardening_Splice_neg_Rot = Cflag_entering_hardening_Splice_neg_Rot;

	flag_entering_hardening_Shr_pos = Cflag_entering_hardening_Shr_pos;
	flag_entering_hardening_Shr_neg = Cflag_entering_hardening_Shr_neg;

	flag_entering_hardening_Splice_Rot = Cflag_entering_hardening_Splice_Rot;
	flag_entering_hardening_Shr = Cflag_entering_hardening_Shr;

	flag_FlexShr_Filure_Rot = Cflag_FlexShr_Filure_Rot;
	flag_Splice_Filure_Rot = Cflag_Splice_Filure_Rot;
	flag_FlexSplice_Filure_Rot = Cflag_FlexSplice_Filure_Rot;
	flag_Flx_Filure_Rot = Cflag_Flx_Filure_Rot;
	flag_Filure_Shr = Cflag_Filure_Shr;

	flagFlur_Rot = CflagFlur_Rot;
	flagFlur_Shr = CflagFlur_Shr;

	f_Bench_Rot = Cf_Bench_Rot;
	f_Bench_Shr = Cf_Bench_Shr;

	d_Bench_Rot = Cd_Bench_Rot;
	d_Bench_Shr = Cd_Bench_Shr;

	dpeakmax = Cdpeakmax;
	dpeakmax_Shr = Cdpeakmax_Shr;

	dpeakmax_inner = Cdpeakmax_inner;
	dpeakmax_inner_Shr = Cdpeakmax_inner_Shr;

	dpeakmax_bench = Cdpeakmax_bench;
	dpeakmax_bench_Shr = Cdpeakmax_bench_Shr;

	dpeakmax_inner_inner = Cdpeakmax_inner_inner;
	dpeakmax_inner_inner_Shr = Cdpeakmax_inner_inner_Shr;

	dpeakmin = Cdpeakmin;
	dpeakmin_Shr = Cdpeakmin_Shr;

	dpeakmin_bench = Cdpeakmin_bench;
	dpeakmin_bench_Shr = Cdpeakmin_bench_Shr;

	dpeakmin_inner = Cdpeakmin_inner;
	dpeakmin_inner_Shr = Cdpeakmin_inner_Shr;

	dpeakmin_inner_inner = Cdpeakmin_inner_inner;
	dpeakmin_inner_inner_Shr = Cdpeakmin_inner_inner_Shr;

	ffmax = Cffmax;
	ffmax_Shr = Cffmax_Shr;

	ffmax_inner = Cffmax_inner;
	ffmax_inner_Shr = Cffmax_inner_Shr;

	ffmax_inner_inner = Cffmax_inner_inner;
	ffmax_inner_inner_Shr = Cffmax_inner_inner_Shr;

	fouterNP_max = CfouterNP_max;
	fouterNP_max_Shr = CfouterNP_max_Shr;

	ffmin = Cffmin;
	ffmin_Shr = Cffmin_Shr;

	ffmin_inner = Cffmin_inner;
	ffmin_inner_Shr = Cffmin_inner_Shr;

	ffmin_inner_inner = Cffmin_inner_inner;
	ffmin_inner_inner_Shr = Cffmin_inner_inner_Shr;

	fouterPN_min = CfouterPN_min;
	fouterPN_min_Shr = CfouterPN_min_Shr;

	spos = Cspos;
	spos_Shr = Cspos_Shr;

	sneg = Csneg;
	sneg_Shr = Csneg_Shr;

	dmgSpos = CdmgSpos;
	dmgSpos_Shr = CdmgSpos_Shr;

	dmgSneg = CdmgSneg;
	dmgSneg_Shr = CdmgSneg_Shr;

	Et = CEt;
	Et_Shr = CEt_Shr;

	Ec = CEc;
	Ec_Shr = CEc_Shr;

	Eposnorm = CEposnorm;
	Eposnorm_Shr = CEposnorm_Shr;

	Enegnorm = CEnegnorm;
	Enegnorm_Shr = CEnegnorm_Shr;

	MtoRref = CMtoRref;
	MtoRref_Shr = CMtoRref_Shr;

	KrR = CKrR;
	KrR_Shr = CKrR_Shr;

	Krel = CKrel;
	Krel_Shr = CKrel_Shr;

	KuR = CKuR;
	KuR_Shr = CKuR_Shr;

	Kun = CKun;
	Kun_Shr = CKun_Shr;

	ER = CER;
	ER_Shr = CER_Shr;

	E = CE;
	E_Shr = CE_Shr;

	//----------------------My Parameters-------------------------

	Tdu = CTdu;
	Tdu_Rot = CTdu_Rot;
	Tdu_Shr = CTdu_Shr;

	TstrainRate = CstrainRate;
	TstrainMax = CstrainMax;
	TstrainMin = CstrainMin;

	TstateFlag = CstateFlag;
	TstateFlag_Shr = CstateFlag_Shr;

	Vcol = CVcol;

	d12 = Cd12;
	f12 = Cf12;
	d_12 = Cd_12;
	d12neg = Cd12neg;
	f12neg = Cf12neg;
	d_12neg = Cd_12neg;
	f_12neg = Cf_12neg;
	d2 = Cd2;
	f2 = Cf2;
	f_12 = Cf_12;
	d_2 = Cd_2;
	f_2 = Cf_2;
	d2neg = Cd2neg;
	f2neg = Cf2neg;
	d_2neg = Cd_2neg;
	f_2neg = Cf_2neg;
	d3 = Cd3;
	f3 = Cf3;
	d_3 = Cd_3;
	f_3 = Cf_3;
	d3neg = Cd3neg;
	f3neg = Cf3neg;
	d_3neg = Cd_3neg;
	f_3neg = Cf_3neg;
	dmgCounter = CdmgCounter;
	R_dcapneg = CR_dcapneg;
	R_dcapneg_Shr = CR_dcapneg_Shr;

	R_dcappos = CR_dcappos;
	R_dcappos_Shr = CR_dcappos_Shr;

	Intcpt_slope_neg = CIntcpt_slope_neg;
	Intcpt_slope_neg_Shr = CIntcpt_slope_neg_Shr;

	R_delUneg = CR_delUneg;

	R_fresneg = CR_fresneg;
	R_fresneg_Shr = CR_fresneg_Shr;

	slope_neg = Cslope_neg;
	slope_neg_Shr = Cslope_neg_Shr;

	dpneg = Cdpneg;
	dpneg_Shr = Cdpneg_Shr;

	Intcpt_slope_pos = CIntcpt_slope_pos;
	Intcpt_slope_pos_Shr = CIntcpt_slope_pos_Shr;

	R_delUpos = CR_delUpos;
	R_delUpos_Shr = CR_delUpos_Shr;
	R_delUneg_Shr = CR_delUneg_Shr;

	R_frespos = CR_frespos;
	R_frespos_Shr = CR_frespos_Shr;

	slope_pos = Cslope_pos;
	slope_pos_Shr = Cslope_pos_Shr;

	dppos = Cdppos;
	dppos_Shr = Cdppos_Shr;

	R_Kdegneg = CR_Kdegneg;
	R_Kdegneg_Shr = CR_Kdegneg_Shr;

	R_fcapneg = CR_fcapneg;
	R_fcapneg_Shr = CR_fcapneg_Shr;

	R_dyneg = CR_dyneg;
	R_dyneg_Shr = CR_dyneg_Shr;

	R_fyneg = CR_fyneg;
	R_fyneg_Shr = CR_fyneg_Shr;

	R_Kdegpos = CR_Kdegpos;
	R_Kdegpos_Shr = CR_Kdegpos_Shr;

	R_fcappos = CR_fcappos;
	R_fcappos_Shr = CR_fcappos_Shr;

	R_fypos = CR_fypos;
	R_fypos_Shr = CR_fypos_Shr;

	R_dypos = CR_dypos;
	R_dypos_Shr = CR_dypos_Shr;

	R_dresneg = CR_dresneg;
	R_dresneg_Shr = CR_dresneg_Shr;

	R_drespos = CR_drespos;
	R_drespos_Shr = CR_drespos_Shr;

	Intcpt_deg_pos = CIntcpt_deg_pos;
	Intcpt_deg_pos_Shr = CIntcpt_deg_pos_Shr;

	Intcpt_Xaxis_pos_Shr = CIntcpt_Xaxis_pos_Shr;
	Intcpt_Xaxis_neg_Shr = CIntcpt_Xaxis_neg_Shr;

	Intcpt_Xaxis_pos = CIntcpt_Xaxis_pos;
	Intcpt_Xaxis_neg = CIntcpt_Xaxis_neg;

	Intcpt_deg_neg = CIntcpt_deg_neg;
	Intcpt_deg_neg_Shr = CIntcpt_deg_neg_Shr;

	Intcpt_res_pos = CIntcpt_res_pos;
	Intcpt_res_pos_Shr = CIntcpt_res_pos_Shr;

	Intcpt_res_neg = CIntcpt_res_neg;
	Intcpt_res_neg_Shr = CIntcpt_res_neg_Shr;

	InCycFac_6 = CInCycFac_6;
	InCycFac_6_Shr = CInCycFac_6_Shr;

	InCycFac_7 = CInCycFac_7;
	InCycFac_7_Shr = CInCycFac_7_Shr;

	InCycFac_neg_6 = CInCycFac_neg_6;
	InCycFac_neg_6_Shr = CInCycFac_neg_6_Shr;

	InCycFac_neg_7 = CInCycFac_neg_7;
	InCycFac_neg_7_Shr = CInCycFac_neg_7_Shr;

	Benchmark56_up = CBenchmark56_up;
	Benchmark56_up_Shr = CBenchmark56_up_Shr;

	Benchmark56_down = CBenchmark56_down;
	Benchmark56_down_Shr = CBenchmark56_down_Shr;

	Benchmark_neg_56_up = CBenchmark_neg_56_up;
	Benchmark_neg_56_down_Shr = CBenchmark_neg_56_down_Shr;

	Benchmark67_up = CBenchmark67_up;
	Benchmark67_up_Shr = CBenchmark67_up_Shr;

	Benchmark67_down = CBenchmark67_down;
	Benchmark67_down_Shr = CBenchmark67_down_Shr;

	Benchmark_neg_67_up = CBenchmark_neg_67_up;
	Benchmark_neg_67_up_Shr = CBenchmark_neg_67_up_Shr;

	Benchmark_neg_67_down = CBenchmark_neg_67_down;
	Benchmark_neg_67_down_Shr = CBenchmark_neg_67_down_Shr;

	spos_pos = Cspos_pos;
	spos_pos_Shr = Cspos_pos_Shr;

	sneg_pos = Csneg_pos;
	sneg_pos_Shr = Csneg_pos_Shr;

	Et_pos = CEt_pos;
	Et_pos_Shr = CEt_pos_Shr;

	Ec_pos = CEc_pos;
	Ec_pos_Shr = CEc_pos_Shr;

	spos_neg = Cspos_neg;
	spos_neg_Shr = Cspos_neg_Shr;

	sneg_neg = Csneg_neg;
	sneg_neg_Shr = Csneg_neg_Shr;

	Et_neg = CEt_neg;
	Et_neg_Shr = CEt_neg_Shr;

	Ec_neg = CEc_neg;
	Ec_neg_Shr = CEc_neg_Shr;

	Epos_pos = CEpos_pos;
	Epos_pos_Shr = CEpos_pos_Shr;

	Eneg_neg = CEneg_neg;
	Eneg_neg_Shr = CEneg_neg_Shr;

	T_area = CT_area;


	Et_hard = CEt_hard;
	Ec_hard = CEc_hard;
	spos_hard = Cspos_hard;
	sneg_hard = Csneg_hard;
	spos_pos_hard = Cspos_pos_hard;
	sneg_pos_hard = Csneg_pos_hard;
	Et_pos_hard = CEt_pos_hard;
	Ec_pos_hard = CEc_pos_hard;
	spos_neg_hard = Cspos_neg_hard;
	sneg_neg_hard = Csneg_neg_hard;
	Et_neg_hard = CEt_neg_hard;
	Ec_neg_hard = CEc_neg_hard;
	Epos_pos_hard = CEpos_pos_hard;
	Eneg_neg_hard = CEneg_neg_hard;
	Epos_hard = CEpos_hard;
	Eneg_hard = CEneg_hard;
	Eposnorm_hard = CEposnorm_hard;
	Enegnorm_hard = CEnegnorm_hard;
	delta_pos_hard = Cdelta_pos_hard;
	delta_pos_max_hard = Cdelta_pos_max_hard;
	delta_neg_hard = Cdelta_neg_hard;
	delta_neg_max_hard = Cdelta_neg_max_hard;
	alpha_pos = Calpha_pos;
	alpha_neg = Calpha_neg;

	//Hardening damage model Shear
	Et_hard_Shr = CEt_hard_Shr;
	Ec_hard_Shr = CEc_hard_Shr;
	spos_hard_Shr = Cspos_hard_Shr;
	sneg_hard_Shr = Csneg_hard_Shr;
	spos_pos_hard_Shr = Cspos_pos_hard_Shr;
	sneg_pos_hard_Shr = Csneg_pos_hard_Shr;
	Et_pos_hard_Shr = CEt_pos_hard_Shr;
	Ec_pos_hard_Shr = CEc_pos_hard_Shr;
	spos_neg_hard_Shr = Cspos_neg_hard_Shr;
	sneg_neg_hard_Shr = Csneg_neg_hard_Shr;
	Et_neg_hard_Shr = CEt_neg_hard_Shr;
	Ec_neg_hard_Shr = CEc_neg_hard_Shr;
	Epos_pos_hard_Shr = CEpos_pos_hard_Shr;
	Eneg_neg_hard_Shr = CEneg_neg_hard_Shr;
	Epos_hard_Shr = CEpos_hard_Shr;
	Eneg_hard_Shr = CEneg_hard_Shr;
	Eposnorm_hard_Shr = CEposnorm_hard_Shr;
	Enegnorm_hard_Shr = CEnegnorm_hard_Shr;
	delta_pos_hard_Shr = Cdelta_pos_hard_Shr;
	delta_pos_max_hard_Shr = Cdelta_pos_max_hard_Shr;
	delta_neg_hard_Shr = Cdelta_neg_hard_Shr;
	delta_neg_max_hard_Shr = Cdelta_neg_max_hard_Shr;
	alpha_pos_Shr = Calpha_pos_Shr;
	alpha_neg_Shr = Calpha_neg_Shr;

	return 0;
}

int GMG_CMAC2D::revertToStart(void)
{
	// Initialize state variables
	// Initial value:
	const double PI = 3.141592653589793238463;

	ek = Cek = ek_Rot = Cek_Rot = Kepos_Rot = Keneg_Rot;
	ek_Shr = Cek_Shr = Kepos_Shr = Keneg_Shr;
	ek_Axi = Cek_Axi = Keneg_Axil;

	P_Axial_average = 0.0;
	P_Axial_sum = 0.0;
	for (int i = 0; i < number_of_averaged_elements_Hardening; i++) {
		P_Axial_Vector[i] = 0.0;
	}


	d = Cstrain = 0.0;
	d_Rot = Cstrain_Rot = 0.0;
	d_Shr = Cstrain_Shr = 0.0;
	d_Axial = Cstrain_Axial = 0.0;

	f = Cstress = 0.0;
	F_Rot = Cstress_Rot = 0.0;
	F_Shr = Cstress_Shr = 0.0;
	F_Axial = Cstress_Axial = 0.0;

	n = 10.0;
	Lamda_KsecS = 0.0;
	Lamda_KsecG = 0.0;
	Lamda_KrelG = 0.0;
	Lamda_KunG = 0.0;
	Lamda_KrelS = 0.0;
	Lamda_KunS = 0.0;

	dypos = fypos_Rot / Kepos_Rot;

	dypos_Shr = fypos_Shr / Kepos_Shr;

	dyneg = fyneg_Rot / Keneg_Rot;

	dyneg_Shr = fyneg_Shr / Keneg_Shr;

	dpeakmax = Cdpeakmax = dpeakmax_bench = Cdpeakmax_bench = fmin(dypos, dypos_splice);
	dpeakmax_Shr = Cdpeakmax_Shr = dpeakmax_bench_Shr = Cdpeakmax_bench_Shr = dypos_Shr;

	dpeakmin = Cdpeakmin = dpeakmin_bench = Cdpeakmin_bench = fmax(dyneg, dyneg_splice);
	dpeakmin_Shr = Cdpeakmin_Shr = dpeakmin_bench_Shr = Cdpeakmin_bench_Shr = dyneg_Shr;

	ffmax = Cffmax = fypos_Rot;
	ffmax_Shr = Cffmax_Shr = fypos_Shr;

	ffmin = Cffmin = fyneg_Rot;
	ffmin_Shr = Cffmin_Shr = fyneg_Shr;


	CMtoRref = MtoRref = 0.0;
	CMtoRref_Shr = MtoRref_Shr = 0.0;

	CKrR = KrR = 0.0;
	CKrR_Shr = KrR_Shr = 0.0;

	CKrel = Krel = 0.0;
	CKrel_Shr = Krel_Shr = 0.0;

	CKuR = KuR = 0.0;
	CKuR_Shr = KuR_Shr = 0.0;

	CKun = Kun = 0.0;
	CKun_Shr = Kun_Shr = 0.0;

	CER = ER = 0.0;
	CER_Shr = ER_Shr = 0.0;

	CE = E = 0.0;
	CE_Shr = E_Shr = 0.0;

	Ed0pos = 0.0;
	Ed0pos_Shr = 0.0;

	Ed0pos_hard = 0.0;
	Ed0pos_hard_Shr = 0.0;

	Ed0neg_hard = 0.0;
	Ed0neg_hard_Shr = 0.0;

	dmgSpos = CdmgSpos = 0.0;
	dmgSpos_Shr = CdmgSpos_Shr = 0.0;

	spos = Cspos = 0.0;
	spos_Shr = Cspos_Shr = 0.0;

	Et = CEt = 0.0;
	Et_Shr = CEt_Shr = 0.0;

	Eposnorm = CEposnorm = 0.0;
	Eposnorm_Shr = CEposnorm_Shr = 0.0;

	Ed0neg = 0.0;
	Ed0neg_Shr = 0.0;

	dmgSneg = CdmgSneg = 0.0;
	dmgSneg_Shr = CdmgSneg_Shr = 0.0;

	sneg = Csneg = 0.0;
	sneg_Shr = Csneg_Shr = 0.0;

	Ec = CEc = 0.0;
	Ec_Shr = CEc_Shr = 0.0;

	Enegnorm = CEnegnorm = 0.0;

	CstateFlag = 0;
	CstateFlag_Shr = 0;

	VyE_to_Vcol_Ratio = CVyE_to_Vcol_Ratio = 0.0;
	Vcap_for_ER = CVcap_for_ER = 0.0;

	flagdmg = 0;
	flagdmg_Shr = 0;

	flagdmg_Hardening_strength = 0;
	flagdmg_Hardening_strength_Shr = 0;

	flagdmg_Hardening = 0;
	flagdmg_Hardening_Shr = 0;

	flag_entering_hardening_Flex_Rot = 0;
	flag_entering_residual_Flex_Rot = 0;
	flag_entering_residual_Shr = 0;
	flag_entering_hardening_Flex_pos_Rot = 0;
	flag_entering_hardening_Flex_neg_Rot = 0;
	flag_entering_hardening_Shr_pos = 0;
	flag_entering_hardening_Shr_neg = 0;
	flag_entering_hardening_Splice_Rot = 0;
	flag_entering_hardening_Shr = 0;

	flag_entering_hardening_Splice_pos_Rot = 0;
	flag_entering_hardening_Splice_neg_Rot = 0;

	flag_FlexShr_Filure_Rot = 0;
	flag_Splice_Filure_Rot = 0;
	flag_FlexSplice_Filure_Rot = 0;
	flag_Flx_Filure_Rot = 0;
	flag_Filure_Shr = 0;

	flagFlur_Rot = 0;
	flagFlur_Shr = 0;


	Et = 0.0;
	Et_Shr = 0.0;

	Ec = 0.0;
	Ec_Shr = 0.0;

	dmgCounter = CdmgCounter = 0;
	dmgCounter_Shr = CdmgCounter_Shr = 0;

	spos_pos = Cspos_pos = 0.0;
	spos_pos_Shr = Cspos_pos_Shr = 0.0;

	sneg_pos = Csneg_pos = 0.0;
	sneg_pos_Shr = Csneg_pos_Shr = 0.0;

	Et_pos = CEt_pos = 0.0;
	Et_pos_Shr = CEt_pos_Shr = 0.0;

	Ec_pos = CEc_pos = 0.0;
	Ec_pos_Shr = CEc_pos_Shr = 0.0;

	spos_neg = Cspos_neg = 0.0;
	spos_neg_Shr = Cspos_neg_Shr = 0.0;

	sneg_neg = Csneg_neg = 0.0;
	sneg_neg_Shr = Csneg_neg_Shr = 0.0;

	Et_neg = CEt_neg = 0.0;
	Et_neg_Shr = CEt_neg_Shr = 0.0;

	Ec_neg = CEc_neg = 0.0;
	Ec_neg_Shr = CEc_neg_Shr = 0.0;

	Epos_pos = CEpos_pos = 0.0;
	Epos_pos_Shr = CEpos_pos_Shr = 0.0;

	Eneg_neg = CEneg_neg = 0.0;
	Eneg_neg_Shr = CEneg_neg_Shr = 0.0;

	T_area = CT_area = 0.0;
	T_area_Shr = CT_area_Shr = 0.0;


	Et_hard = CEt_hard = 0;
	Ec_hard = CEc_hard = 0;
	spos_hard = Cspos_hard = 0;
	sneg_hard = Csneg_hard = 0;
	spos_pos_hard = Cspos_pos_hard = 0;
	sneg_pos_hard = Csneg_pos_hard = 0;
	Et_pos_hard = CEt_pos_hard = 0;
	Ec_pos_hard = CEc_pos_hard = 0;
	spos_neg_hard = Cspos_neg_hard = 0;
	sneg_neg_hard = Csneg_neg_hard = 0;
	Et_neg_hard = CEt_neg_hard = 0;
	Ec_neg_hard = CEc_neg_hard = 0;
	Epos_pos_hard = CEpos_pos_hard = 0;
	Eneg_neg_hard = CEneg_neg_hard = 0;
	Epos_hard = CEpos_hard = 0;
	Eneg_hard = CEneg_hard = 0;
	Eposnorm_hard = CEposnorm_hard = 0;
	Enegnorm_hard = CEnegnorm_hard = 0;
	delta_pos_hard = Cdelta_pos_hard = 0;
	delta_pos_max_hard = Cdelta_pos_max_hard = 0;
	delta_neg_hard = Cdelta_neg_hard = 0;
	delta_neg_max_hard = Cdelta_neg_max_hard = 0;
	alpha_pos = Calpha_pos = 0;
	alpha_neg = Calpha_neg = 0;

	//Hardening damage model Shear
	Et_hard_Shr = CEt_hard_Shr = 0;
	Ec_hard_Shr = CEc_hard_Shr = 0;
	spos_hard_Shr = Cspos_hard_Shr = 0;
	sneg_hard_Shr = Csneg_hard_Shr = 0;
	spos_pos_hard_Shr = Cspos_pos_hard_Shr = 0;
	sneg_pos_hard_Shr = Csneg_pos_hard_Shr = 0;
	Et_pos_hard_Shr = CEt_pos_hard_Shr = 0;
	Ec_pos_hard_Shr = CEc_pos_hard_Shr = 0;
	spos_neg_hard_Shr = Cspos_neg_hard_Shr = 0;
	sneg_neg_hard_Shr = Csneg_neg_hard_Shr = 0;
	Et_neg_hard_Shr = CEt_neg_hard_Shr = 0;
	Ec_neg_hard_Shr = CEc_neg_hard_Shr = 0;
	Epos_pos_hard_Shr = CEpos_pos_hard_Shr = 0;
	Eneg_neg_hard_Shr = CEneg_neg_hard_Shr = 0;
	Epos_hard_Shr = CEpos_hard_Shr = 0;
	Eneg_hard_Shr = CEneg_hard_Shr = 0;
	Eposnorm_hard_Shr = CEposnorm_hard_Shr = 0;
	Enegnorm_hard_Shr = CEnegnorm_hard_Shr = 0;
	delta_pos_hard_Shr = Cdelta_pos_hard_Shr = 0;
	delta_pos_max_hard_Shr = Cdelta_pos_max_hard_Shr = 0;
	delta_neg_hard_Shr = Cdelta_neg_hard_Shr = 0;
	delta_neg_max_hard_Shr = Cdelta_neg_max_hard_Shr = 0;
	alpha_pos_Shr = Calpha_pos_Shr = 0;
	alpha_neg_Shr = Calpha_neg_Shr = 0;

	return 0;
}

UniaxialMaterial *
GMG_CMAC2D::getCopy(void)
{
	GMG_CMAC2D *theCopy = new GMG_CMAC2D(this->getTag(),

		H, B, C_C,
		L, La,
		fc,
		fyL, dbL, nL_EW, nL_NS,
		fyT, dbT, n_leg, S,
		lb,
		P_Axial,
		PVM_Flag,
		number_of_averaged_elements_Elastic,
		number_of_averaged_elements_Hardening,
		Retrofit_flag,
		fje,
		tj,
		bj);

	// Converged state variables
	theCopy->CstateFlag = CstateFlag;
	theCopy->CstateFlag_Shr = CstateFlag_Shr;

	theCopy->Cstrain = Cstrain;
	theCopy->Cstrain_Rot = Cstrain_Rot;
	theCopy->Cstrain_Shr = Cstrain_Shr;
	theCopy->Cstrain_Axial = Cstrain_Axial;

	theCopy->Cstress = Cstress;
	theCopy->Cstress_Rot = Cstress_Rot;
	theCopy->Cstress_Shr = Cstress_Shr;
	theCopy->Cstress_Axial = Cstress_Axial;

	theCopy->Cek = Cek;
	theCopy->Cek_Rot = Cek_Rot;
	theCopy->Cek_Shr = Cek_Shr;
	theCopy->Cek_Axi = Cek_Axi;


	theCopy->CInCycFac = CInCycFac;
	theCopy->CInCycFac_Shr = CInCycFac_Shr;

	theCopy->CInCycFac_6 = CInCycFac_6;
	theCopy->CInCycFac_6_Shr = CInCycFac_6_Shr;

	theCopy->CInCycFac_7 = CInCycFac_7;
	theCopy->CInCycFac_7_Shr = CInCycFac_7_Shr;

	theCopy->CInCycFac_neg_6 = CInCycFac_neg_6;
	theCopy->CInCycFac_neg_6_Shr = CInCycFac_neg_6_Shr;

	theCopy->CInCycFac_neg_7 = CInCycFac_neg_7;
	theCopy->CInCycFac_neg_7_Shr = CInCycFac_neg_7_Shr;

	theCopy->Cflagstop = Cflagstop;

	theCopy->Cflagdmg = Cflagdmg;
	theCopy->Cflagdmg_Shr = Cflagdmg_Shr;

	theCopy->Cflagdmg_Hardening_strength = Cflagdmg_Hardening_strength;
	theCopy->Cflagdmg_Hardening_strength_Shr = Cflagdmg_Hardening_strength_Shr;

	theCopy->Cflagdmg_Hardening = Cflagdmg_Hardening;
	theCopy->Cflagdmg_Hardening_Shr = Cflagdmg_Hardening_Shr;

	theCopy->CVyE_to_Vcol_Ratio = CVyE_to_Vcol_Ratio;
	theCopy->CVcap_for_ER = CVcap_for_ER;

	theCopy->Cflag_entering_hardening_Flex_Rot = Cflag_entering_hardening_Flex_Rot;
	theCopy->Cflag_entering_residual_Flex_Rot = Cflag_entering_residual_Flex_Rot;
	theCopy->Cflag_entering_residual_Shr = Cflag_entering_residual_Shr;

	theCopy->flag_entering_hardening_Flex_pos_Rot = flag_entering_hardening_Flex_pos_Rot;
	theCopy->Cflag_entering_hardening_Flex_pos_Rot = Cflag_entering_hardening_Flex_pos_Rot;
	theCopy->flag_entering_hardening_Flex_neg_Rot = flag_entering_hardening_Flex_neg_Rot;
	theCopy->Cflag_entering_hardening_Flex_neg_Rot = Cflag_entering_hardening_Flex_neg_Rot;

	theCopy->flag_entering_hardening_Splice_pos_Rot = flag_entering_hardening_Splice_pos_Rot;
	theCopy->Cflag_entering_hardening_Splice_pos_Rot = Cflag_entering_hardening_Splice_pos_Rot;
	theCopy->flag_entering_hardening_Splice_neg_Rot = flag_entering_hardening_Splice_neg_Rot;
	theCopy->Cflag_entering_hardening_Splice_neg_Rot = Cflag_entering_hardening_Splice_neg_Rot;

	theCopy->Cflag_entering_hardening_Shr_pos = Cflag_entering_hardening_Shr_pos;
	theCopy->flag_entering_hardening_Shr_pos = flag_entering_hardening_Shr_pos;
	theCopy->Cflag_entering_hardening_Shr_neg = Cflag_entering_hardening_Shr_neg;
	theCopy->flag_entering_hardening_Shr_neg = flag_entering_hardening_Shr_neg;

	theCopy->Cflag_entering_hardening_Splice_Rot = Cflag_entering_hardening_Splice_Rot;
	theCopy->Cflag_entering_hardening_Shr = Cflag_entering_hardening_Shr;

	theCopy->Cflag_FlexShr_Filure_Rot = Cflag_FlexShr_Filure_Rot;
	theCopy->Cflag_Splice_Filure_Rot = Cflag_Splice_Filure_Rot;
	theCopy->Cflag_FlexSplice_Filure_Rot = Cflag_FlexSplice_Filure_Rot;
	theCopy->Cflag_Flx_Filure_Rot = Cflag_Flx_Filure_Rot;
	theCopy->Cflag_Filure_Shr = Cflag_Filure_Shr;

	theCopy->CflagFlur_Rot = CflagFlur_Rot;
	theCopy->CflagFlur_Shr = CflagFlur_Shr;

	theCopy->Cf_Bench_Rot = Cf_Bench_Rot;
	theCopy->Cf_Bench_Shr = Cf_Bench_Shr;

	theCopy->Cd_Bench_Rot = Cd_Bench_Rot;
	theCopy->Cd_Bench_Shr = Cd_Bench_Shr;

	theCopy->Cdpeakmax = Cdpeakmax;
	theCopy->Cdpeakmax_Shr = Cdpeakmax_Shr;

	theCopy->Cdpeakmax_bench = Cdpeakmax_bench;
	theCopy->Cdpeakmax_bench_Shr = Cdpeakmax_bench_Shr;

	theCopy->Cdpeakmax_inner = Cdpeakmax_inner;
	theCopy->Cdpeakmax_inner_Shr = Cdpeakmax_inner_Shr;

	theCopy->Cdpeakmax_inner_inner = Cdpeakmax_inner_inner;
	theCopy->Cdpeakmax_inner_inner_Shr = Cdpeakmax_inner_inner_Shr;

	theCopy->Cdpeakmin = Cdpeakmin;
	theCopy->Cdpeakmin_Shr = Cdpeakmin_Shr;

	theCopy->Cdpeakmin_bench = Cdpeakmin_bench;
	theCopy->Cdpeakmin_bench_Shr = Cdpeakmin_bench_Shr;

	theCopy->Cdpeakmin_inner = Cdpeakmin_inner;
	theCopy->Cdpeakmin_inner_Shr = Cdpeakmin_inner_Shr;

	theCopy->Cdpeakmin_inner_inner = Cdpeakmin_inner_inner;
	theCopy->Cdpeakmin_inner_inner_Shr = Cdpeakmin_inner_inner_Shr;

	theCopy->Cffmax = Cffmax;
	theCopy->Cffmax_Shr = Cffmax_Shr;

	theCopy->Cffmax_inner = Cffmax_inner;
	theCopy->Cffmax_inner_Shr = Cffmax_inner_Shr;

	theCopy->Cffmax_inner_inner = Cffmax_inner_inner;
	theCopy->Cffmax_inner_inner_Shr = Cffmax_inner_inner_Shr;

	theCopy->CfouterNP_max = CfouterNP_max;
	theCopy->CfouterNP_max_Shr = CfouterNP_max_Shr;

	theCopy->Cffmin = Cffmin;
	theCopy->Cffmin_Shr = Cffmin_Shr;

	theCopy->Cffmin_inner = Cffmin_inner;
	theCopy->Cffmin_inner_Shr = Cffmin_inner_Shr;

	theCopy->Cffmin_inner_inner = Cffmin_inner_inner;
	theCopy->Cffmin_inner_inner_Shr = Cffmin_inner_inner_Shr;

	theCopy->CfouterPN_min = CfouterPN_min;
	theCopy->CfouterPN_min_Shr = CfouterPN_min_Shr;

	theCopy->Cspos = Cspos;
	theCopy->Cspos_Shr = Cspos_Shr;

	theCopy->Csneg = Csneg;
	theCopy->Csneg_Shr = Csneg_Shr;

	theCopy->CdmgSpos = CdmgSpos;
	theCopy->CdmgSpos_Shr = CdmgSpos_Shr;

	theCopy->CdmgSneg = CdmgSneg;
	theCopy->CdmgSneg_Shr = CdmgSneg_Shr;

	theCopy->CEt = CEt;
	theCopy->CEt_Shr = CEt_Shr;

	theCopy->CEc = CEc;
	theCopy->CEc_Shr = CEc_Shr;

	theCopy->CEposnorm = CEposnorm;
	theCopy->CEposnorm_Shr = CEposnorm_Shr;

	theCopy->CEnegnorm = CEnegnorm;
	theCopy->CEnegnorm_Shr = CEnegnorm_Shr;

	theCopy->CBenchmark56_up = CBenchmark56_up;
	theCopy->CBenchmark56_up_Shr = CBenchmark56_up_Shr;

	theCopy->CBenchmark56_down = CBenchmark56_down;
	theCopy->CBenchmark56_down_Shr = CBenchmark56_down_Shr;

	theCopy->CBenchmark_neg_56_up = CBenchmark_neg_56_up;
	theCopy->CBenchmark_neg_56_up_Shr = CBenchmark_neg_56_up_Shr;

	theCopy->CBenchmark_neg_56_down = CBenchmark_neg_56_down;
	theCopy->CBenchmark_neg_56_down_Shr = CBenchmark_neg_56_down_Shr;

	theCopy->CBenchmark67_up = CBenchmark67_up;
	theCopy->CBenchmark67_up_Shr = CBenchmark67_up_Shr;

	theCopy->CBenchmark67_down = CBenchmark67_down;
	theCopy->CBenchmark67_down_Shr = CBenchmark67_down_Shr;

	theCopy->CBenchmark_neg_67_up = CBenchmark_neg_67_up;
	theCopy->CBenchmark_neg_67_up_Shr = CBenchmark_neg_67_up_Shr;

	theCopy->CBenchmark_neg_67_down = CBenchmark_neg_67_down;
	theCopy->CBenchmark_neg_67_down_Shr = CBenchmark_neg_67_down_Shr;

	theCopy->Cspos_pos = Cspos_pos;
	theCopy->Cspos_pos_Shr = Cspos_pos_Shr;

	theCopy->Csneg_pos = Csneg_pos;
	theCopy->Csneg_pos_Shr = Csneg_pos_Shr;

	theCopy->CEt_pos = CEt_pos;
	theCopy->CEt_pos_Shr = CEt_pos_Shr;

	theCopy->CEc_pos = CEc_pos;
	theCopy->CEc_pos_Shr = CEc_pos_Shr;

	theCopy->Cspos_neg = Cspos_neg;
	theCopy->Cspos_neg_Shr = Cspos_neg_Shr;

	theCopy->Csneg_neg = Csneg_neg;
	theCopy->Csneg_neg_Shr = Csneg_neg_Shr;

	theCopy->CEt_neg = CEt_neg;
	theCopy->CEt_neg_Shr = CEt_neg_Shr;

	theCopy->CEt_hard = CEt_hard;
	theCopy->CEc_hard = CEc_hard;
	theCopy->Cspos_hard = Cspos_hard;
	theCopy->Csneg_hard = Csneg_hard;
	theCopy->Cspos_pos_hard = Cspos_pos_hard;
	theCopy->Csneg_pos_hard = Csneg_pos_hard;
	theCopy->CEt_pos_hard = CEt_pos_hard;
	theCopy->CEc_pos_hard = CEc_pos_hard;
	theCopy->Cspos_neg_hard = Cspos_neg_hard;
	theCopy->Csneg_neg_hard = Csneg_neg_hard;
	theCopy->CEt_neg_hard = CEt_neg_hard;
	theCopy->CEc_neg_hard = CEc_neg_hard;
	theCopy->CEpos_pos_hard = CEpos_pos_hard;
	theCopy->CEneg_neg_hard = CEneg_neg_hard;
	theCopy->CEpos_hard = CEpos_hard;
	theCopy->CEneg_hard = CEneg_hard;
	theCopy->CEposnorm_hard = CEposnorm_hard;
	theCopy->CEnegnorm_hard = CEnegnorm_hard;
	theCopy->Cdelta_pos_hard = Cdelta_pos_hard;
	theCopy->Cdelta_pos_max_hard = Cdelta_pos_max_hard;
	theCopy->Cdelta_neg_hard = Cdelta_neg_hard;
	theCopy->Cdelta_neg_max_hard = Cdelta_neg_max_hard;
	theCopy->Calpha_pos = Calpha_pos;
	theCopy->Calpha_neg = Calpha_neg;

	//Hardening damage model Shear
	theCopy->CEt_hard_Shr = CEt_hard_Shr;
	theCopy->CEc_hard_Shr = CEc_hard_Shr;
	theCopy->Cspos_hard_Shr = Cspos_hard_Shr;
	theCopy->Csneg_hard_Shr = Csneg_hard_Shr;
	theCopy->Cspos_pos_hard_Shr = Cspos_pos_hard_Shr;
	theCopy->Csneg_pos_hard_Shr = Csneg_pos_hard_Shr;
	theCopy->CEt_pos_hard_Shr = CEt_pos_hard_Shr;
	theCopy->CEc_pos_hard_Shr = CEc_pos_hard_Shr;
	theCopy->Cspos_neg_hard_Shr = Cspos_neg_hard_Shr;
	theCopy->Csneg_neg_hard_Shr = Csneg_neg_hard_Shr;
	theCopy->CEt_neg_hard_Shr = CEt_neg_hard_Shr;
	theCopy->CEc_neg_hard_Shr = CEc_neg_hard_Shr;
	theCopy->CEpos_pos_hard_Shr = CEpos_pos_hard_Shr;
	theCopy->CEneg_neg_hard_Shr = CEneg_neg_hard_Shr;
	theCopy->CEpos_hard_Shr = CEpos_hard_Shr;
	theCopy->CEneg_hard_Shr = CEneg_hard_Shr;
	theCopy->CEposnorm_hard_Shr = CEposnorm_hard_Shr;
	theCopy->CEnegnorm_hard_Shr = CEnegnorm_hard_Shr;
	theCopy->Cdelta_pos_hard_Shr = Cdelta_pos_hard_Shr;
	theCopy->Cdelta_pos_max_hard_Shr = Cdelta_pos_max_hard_Shr;
	theCopy->Cdelta_neg_hard_Shr = Cdelta_neg_hard_Shr;
	theCopy->Cdelta_neg_max_hard_Shr = Cdelta_neg_max_hard_Shr;
	theCopy->Calpha_pos_Shr = Calpha_pos_Shr;
	theCopy->Calpha_neg_Shr = Calpha_neg_Shr;

	theCopy->CEc_neg = CEc_neg;
	theCopy->CEc_neg_Shr = CEc_neg_Shr;

	theCopy->CEpos_pos = CEpos_pos;
	theCopy->CEpos_pos_Shr = CEpos_pos_Shr;

	theCopy->CEneg_neg = CEneg_neg;
	theCopy->CEneg_neg_Shr = CEneg_neg_Shr;

	theCopy->CT_area = CT_area;

	theCopy->CMtoRref = CMtoRref;
	theCopy->CMtoRref_Shr = CMtoRref_Shr;

	theCopy->Ckd = Ckd;

	theCopy->CKrR = CKrR;
	theCopy->CKrR_Shr = CKrR_Shr;

	theCopy->CKrel = CKrel;
	theCopy->CKrel_Shr = CKrel_Shr;

	theCopy->CKuR = CKuR;
	theCopy->CKuR_Shr = CKuR_Shr;

	theCopy->CKun = CKun;
	theCopy->CKun_Shr = CKun_Shr;

	theCopy->CER = CER;
	theCopy->CER_Shr = CER_Shr;

	theCopy->CE = CE;
	theCopy->CE_Shr = CE_Shr;

	theCopy->CIntcpt_res_pos = CIntcpt_res_pos;
	theCopy->CIntcpt_res_pos_Shr = CIntcpt_res_pos_Shr;

	theCopy->CIntcpt_res_neg = CIntcpt_res_neg;
	theCopy->CIntcpt_res_neg_Shr = CIntcpt_res_neg_Shr;

	theCopy->CIntcpt_deg_neg = CIntcpt_deg_neg;
	theCopy->CIntcpt_deg_neg_Shr = CIntcpt_deg_neg_Shr;

	theCopy->CIntcpt_deg_pos = CIntcpt_deg_pos;
	theCopy->CIntcpt_deg_pos_Shr = CIntcpt_deg_pos_Shr;

	theCopy->CIntcpt_Xaxis_pos = CIntcpt_Xaxis_pos;
	theCopy->CIntcpt_Xaxis_neg = CIntcpt_Xaxis_neg;

	theCopy->CIntcpt_Xaxis_pos_Shr = CIntcpt_Xaxis_pos_Shr;
	theCopy->CIntcpt_Xaxis_neg_Shr = CIntcpt_Xaxis_neg_Shr;


	// Trial state variables
	theCopy->TstateFlag = TstateFlag;
	theCopy->TstateFlag_Shr = TstateFlag_Shr;

	theCopy->d = d;
	theCopy->d_Rot = d_Rot;
	theCopy->d_Shr = d_Shr;

	theCopy->f = f;
	theCopy->F_Rot = F_Rot;
	theCopy->F_Shr = F_Shr;

	theCopy->ek = ek;
	theCopy->ek_Rot = ek_Rot;
	theCopy->ek_Shr = ek_Shr;
	theCopy->ek_Axi = ek_Axi;

	theCopy->P_Axial_average = P_Axial_average;

	theCopy->InCycFac = InCycFac;
	theCopy->InCycFac_Shr = InCycFac_Shr;

	theCopy->InCycFac_6 = InCycFac_6;
	theCopy->InCycFac_6_Shr = InCycFac_6_Shr;

	theCopy->InCycFac_7 = InCycFac_7;
	theCopy->InCycFac_7_Shr = InCycFac_7_Shr;

	theCopy->InCycFac_neg_6 = InCycFac_neg_6;
	theCopy->InCycFac_neg_6_Shr = InCycFac_neg_6_Shr;

	theCopy->InCycFac_neg_7 = InCycFac_neg_7;
	theCopy->InCycFac_neg_7_Shr = InCycFac_neg_7_Shr;

	theCopy->flagstop = flagstop;

	theCopy->flagdmg = flagdmg;
	theCopy->flagdmg_Shr = flagdmg_Shr;

	theCopy->flagdmg_Hardening_strength = flagdmg_Hardening_strength;
	theCopy->flagdmg_Hardening_strength_Shr = flagdmg_Hardening_strength_Shr;

	theCopy->flagdmg_Hardening = flagdmg_Hardening;
	theCopy->flagdmg_Hardening_Shr = flagdmg_Hardening_Shr;

	theCopy->VyE_to_Vcol_Ratio = VyE_to_Vcol_Ratio;
	theCopy->Vcap_for_ER = Vcap_for_ER;

	theCopy->flag_entering_hardening_Flex_Rot = flag_entering_hardening_Flex_Rot;
	theCopy->flag_entering_residual_Flex_Rot = flag_entering_residual_Flex_Rot;
	theCopy->flag_entering_residual_Shr = flag_entering_residual_Shr;

	theCopy->flag_entering_hardening_Splice_Rot = flag_entering_hardening_Splice_Rot;
	theCopy->flag_entering_hardening_Shr = flag_entering_hardening_Shr;

	theCopy->flag_FlexShr_Filure_Rot = flag_FlexShr_Filure_Rot;
	theCopy->flag_Splice_Filure_Rot = flag_Splice_Filure_Rot;
	theCopy->flag_FlexSplice_Filure_Rot = flag_FlexSplice_Filure_Rot;
	theCopy->flag_Flx_Filure_Rot = flag_Flx_Filure_Rot;
	theCopy->flag_Filure_Shr = flag_Filure_Shr;

	theCopy->flagFlur_Rot = flagFlur_Rot;
	theCopy->flagFlur_Shr = flagFlur_Shr;

	theCopy->f_Bench_Rot = f_Bench_Rot;
	theCopy->f_Bench_Shr = f_Bench_Shr;

	theCopy->d_Bench_Rot = d_Bench_Rot;
	theCopy->d_Bench_Shr = d_Bench_Shr;

	theCopy->dpeakmax = dpeakmax;
	theCopy->dpeakmax_Shr = dpeakmax_Shr;

	theCopy->dpeakmax_bench = dpeakmax_bench;
	theCopy->dpeakmax_bench_Shr = dpeakmax_bench_Shr;

	theCopy->dpeakmax_inner = dpeakmax_inner;
	theCopy->dpeakmax_inner_Shr = dpeakmax_inner_Shr;

	theCopy->dpeakmin = dpeakmin;
	theCopy->dpeakmin_Shr = dpeakmin_Shr;

	theCopy->dpeakmin_bench = dpeakmin_bench;
	theCopy->dpeakmin_bench_Shr = dpeakmin_bench_Shr;

	theCopy->dpeakmin_inner = dpeakmin_inner;
	theCopy->dpeakmin_inner_Shr = dpeakmin_inner_Shr;

	theCopy->dpeakmax_inner_inner = dpeakmax_inner_inner;
	theCopy->dpeakmax_inner_inner_Shr = dpeakmax_inner_inner_Shr;

	theCopy->dpeakmin_inner_inner = dpeakmin_inner_inner;
	theCopy->dpeakmin_inner_inner_Shr = dpeakmin_inner_inner_Shr;

	theCopy->ffmax = ffmax;
	theCopy->ffmax_Shr = ffmax_Shr;

	theCopy->ffmin = ffmin;
	theCopy->ffmin_Shr = ffmin_Shr;

	theCopy->ffmin_inner = ffmin_inner;
	theCopy->ffmin_inner_Shr = ffmin_inner_Shr;

	theCopy->ffmax_inner = ffmax_inner;
	theCopy->ffmax_inner_Shr = ffmax_inner_Shr;

	theCopy->ffmax_inner_inner = ffmax_inner_inner;
	theCopy->ffmax_inner_inner_Shr = ffmax_inner_inner_Shr;

	theCopy->fouterNP_max = fouterNP_max;
	theCopy->fouterNP_max_Shr = fouterNP_max_Shr;

	theCopy->ffmin_inner_inner = ffmin_inner_inner;
	theCopy->ffmin_inner_inner_Shr = ffmin_inner_inner_Shr;

	theCopy->fouterPN_min = fouterPN_min;
	theCopy->fouterPN_min_Shr = fouterPN_min_Shr;

	theCopy->spos = spos;
	theCopy->spos_Shr = spos_Shr;

	theCopy->sneg = sneg;
	theCopy->sneg_Shr = sneg_Shr;

	theCopy->dmgSpos = dmgSpos;
	theCopy->dmgSpos_Shr = dmgSpos_Shr;

	theCopy->dmgSneg = dmgSneg;
	theCopy->dmgSneg_Shr = dmgSneg_Shr;

	theCopy->Et = Et;
	theCopy->Et_Shr = Et_Shr;

	theCopy->Ec = Ec;
	theCopy->Ec_Shr = Ec_Shr;

	theCopy->Eposnorm = Eposnorm;
	theCopy->Eposnorm_Shr = Eposnorm_Shr;

	theCopy->Enegnorm = Enegnorm;
	theCopy->Enegnorm_Shr = Enegnorm_Shr;

	theCopy->dmgCounter = dmgCounter;

	theCopy->R_dcapneg = R_dcapneg;
	theCopy->R_dcapneg_Shr = R_dcapneg_Shr;

	theCopy->R_dcappos = R_dcappos;
	theCopy->R_dcappos_Shr = R_dcappos_Shr;

	theCopy->Intcpt_slope_pos = Intcpt_slope_pos;
	theCopy->Intcpt_slope_pos_Shr = Intcpt_slope_pos_Shr;

	theCopy->CIntcpt_slope_pos = CIntcpt_slope_pos;
	theCopy->CIntcpt_slope_pos_Shr = CIntcpt_slope_pos_Shr;

	theCopy->R_fcappos = R_fcappos;
	theCopy->R_fcappos_Shr = R_fcappos_Shr;

	theCopy->CR_fcappos = CR_fcappos;
	theCopy->CR_fcappos_Shr = CR_fcappos_Shr;

	theCopy->R_fcapneg = R_fcapneg;
	theCopy->R_fcapneg_Shr = R_fcapneg_Shr;

	theCopy->CR_fcapneg = CR_fcapneg;
	theCopy->CR_fcapneg_Shr = CR_fcapneg_Shr;

	theCopy->R_dypos = R_dypos;
	theCopy->R_dypos_Shr = R_dypos_Shr;

	theCopy->fs = fs;
	theCopy->lb_deg = lb_deg;
	theCopy->fs_deg = fs_deg;

	theCopy->CR_dypos = CR_dypos;
	theCopy->CR_dypos_Shr = CR_dypos_Shr;

	theCopy->Intcpt_slope_neg = Intcpt_slope_neg;
	theCopy->Intcpt_slope_neg_Shr = Intcpt_slope_neg_Shr;

	theCopy->CIntcpt_slope_neg = CIntcpt_slope_neg;
	theCopy->CIntcpt_slope_neg_Shr = CIntcpt_slope_neg_Shr;

	theCopy->R_delUneg = R_delUneg;
	theCopy->CR_delUneg = CR_delUneg;

	theCopy->R_delUpos = R_delUpos;

	theCopy->CR_delUpos = CR_delUpos;

	theCopy->R_delUpos_Shr = R_delUpos_Shr;
	theCopy->R_delUneg_Shr = R_delUneg_Shr;

	theCopy->CR_delUpos_Shr = CR_delUpos_Shr;
	theCopy->CR_delUneg_Shr = CR_delUneg_Shr;

	theCopy->R_frespos = R_frespos;
	theCopy->R_frespos_Shr = R_frespos_Shr;

	theCopy->CR_frespos = CR_frespos;
	theCopy->CR_frespos_Shr = CR_frespos_Shr;

	theCopy->R_fresneg = R_fresneg;
	theCopy->R_fresneg_Shr = R_fresneg_Shr;

	theCopy->CR_fresneg = CR_fresneg;
	theCopy->CR_fresneg_Shr = CR_fresneg_Shr;

	theCopy->slope_pos = slope_pos;
	theCopy->slope_pos_Shr = slope_pos_Shr;

	theCopy->Cslope_pos = Cslope_pos;
	theCopy->Cslope_pos_Shr = Cslope_pos_Shr;

	theCopy->slope_neg = slope_neg;
	theCopy->slope_neg_Shr = slope_neg_Shr;

	theCopy->Cslope_neg = Cslope_neg;
	theCopy->Cslope_neg_Shr = Cslope_neg_Shr;

	theCopy->dppos = dppos;
	theCopy->dppos_Shr = dppos_Shr;

	theCopy->Cdppos = Cdppos;
	theCopy->Cdppos_Shr = Cdppos_Shr;

	theCopy->dpneg = dpneg;
	theCopy->dpneg_Shr = dpneg_Shr;

	theCopy->Cdpneg = Cdpneg;
	theCopy->Cdpneg_Shr = Cdpneg_Shr;

	theCopy->R_Kdegneg = R_Kdegneg;
	theCopy->R_Kdegneg_Shr = R_Kdegneg_Shr;

	theCopy->CR_Kdegneg = CR_Kdegneg;
	theCopy->CR_Kdegneg_Shr = CR_Kdegneg_Shr;

	theCopy->R_dyneg = R_dyneg;
	theCopy->R_dyneg_Shr = R_dyneg_Shr;

	theCopy->CR_dyneg = CR_dyneg;
	theCopy->CR_dyneg_Shr = CR_dyneg_Shr;

	theCopy->R_Kdegpos = R_Kdegpos;
	theCopy->R_Kdegpos_Shr = R_Kdegpos_Shr;

	theCopy->CR_Kdegpos = CR_Kdegpos;
	theCopy->CR_Kdegpos_Shr = CR_Kdegpos_Shr;

	theCopy->R_fypos = R_fypos;
	theCopy->R_fypos_Shr = R_fypos_Shr;

	theCopy->CR_fypos = CR_fypos;
	theCopy->CR_fypos_Shr = CR_fypos_Shr;

	theCopy->R_fyneg = R_fyneg;
	theCopy->R_fyneg_Shr = R_fyneg_Shr;

	theCopy->CR_fyneg = CR_fyneg;
	theCopy->CR_fyneg_Shr = CR_fyneg_Shr;

	theCopy->CR_dcapneg = CR_dcapneg;
	theCopy->CR_dcapneg_Shr = CR_dcapneg_Shr;

	theCopy->CR_dcappos = CR_dcappos;
	theCopy->CR_dcappos_Shr = CR_dcappos_Shr;

	theCopy->R_dresneg = R_dresneg;
	theCopy->R_dresneg_Shr = R_dresneg_Shr;

	theCopy->CR_dresneg = CR_dresneg;
	theCopy->CR_dresneg_Shr = CR_dresneg_Shr;

	theCopy->R_drespos = R_drespos;
	theCopy->R_drespos_Shr = R_drespos_Shr;

	theCopy->CR_drespos = CR_drespos;
	theCopy->CR_drespos_Shr = CR_drespos_Shr;

	theCopy->Intcpt_deg_pos = Intcpt_deg_pos;
	theCopy->Intcpt_deg_pos_Shr = Intcpt_deg_pos_Shr;

	theCopy->Intcpt_Xaxis_pos_Shr = Intcpt_Xaxis_pos_Shr;
	theCopy->Intcpt_Xaxis_neg_Shr = Intcpt_Xaxis_neg_Shr;

	theCopy->Intcpt_Xaxis_pos = Intcpt_Xaxis_pos;
	theCopy->Intcpt_Xaxis_neg = Intcpt_Xaxis_neg;

	theCopy->Intcpt_deg_neg = Intcpt_deg_neg;
	theCopy->Intcpt_deg_neg_Shr = Intcpt_deg_neg_Shr;

	theCopy->Intcpt_res_pos = Intcpt_res_pos;
	theCopy->Intcpt_res_pos_Shr = Intcpt_res_pos_Shr;

	theCopy->Intcpt_res_neg = Intcpt_res_neg;
	theCopy->Intcpt_res_neg_Shr = Intcpt_res_neg_Shr;

	theCopy->Benchmark56_up = Benchmark56_up;
	theCopy->Benchmark56_up_Shr = Benchmark56_up_Shr;

	theCopy->Benchmark56_down = Benchmark56_down;
	theCopy->Benchmark56_down_Shr = Benchmark56_down_Shr;

	theCopy->Benchmark_neg_56_up = Benchmark_neg_56_up;
	theCopy->Benchmark_neg_56_up_Shr = Benchmark_neg_56_up_Shr;

	theCopy->Benchmark_neg_56_down = Benchmark_neg_56_down;
	theCopy->Benchmark_neg_56_down_Shr = Benchmark_neg_56_down_Shr;

	theCopy->Benchmark67_up = Benchmark67_up;
	theCopy->Benchmark67_up_Shr = Benchmark67_up_Shr;

	theCopy->Benchmark67_down = Benchmark67_down;
	theCopy->Benchmark67_down_Shr = Benchmark67_down_Shr;

	theCopy->Benchmark_neg_67_up = Benchmark_neg_67_up;
	theCopy->Benchmark_neg_67_up_Shr = Benchmark_neg_67_up_Shr;

	theCopy->Benchmark_neg_67_down = Benchmark_neg_67_down;
	theCopy->Benchmark_neg_67_down_Shr = Benchmark_neg_67_down_Shr;

	theCopy->spos_pos = spos_pos;
	theCopy->spos_pos_Shr = spos_pos_Shr;

	theCopy->sneg_pos = sneg_pos;
	theCopy->sneg_pos_Shr = sneg_pos_Shr;

	theCopy->Et_pos = Et_pos;
	theCopy->Et_pos_Shr = Et_pos_Shr;

	theCopy->Ec_pos = Ec_pos;
	theCopy->Ec_pos_Shr = Ec_pos_Shr;

	theCopy->spos_neg = spos_neg;
	theCopy->spos_neg_Shr = spos_neg_Shr;

	theCopy->sneg_neg = sneg_neg;
	theCopy->sneg_neg_Shr = sneg_neg_Shr;

	theCopy->Et_neg = Et_neg;
	theCopy->Et_neg_Shr = Et_neg_Shr;

	theCopy->CEt_hard = CEt_hard;
	theCopy->CEc_hard = CEc_hard;
	theCopy->Cspos_hard = Cspos_hard;
	theCopy->Csneg_hard = Csneg_hard;
	theCopy->Cspos_pos_hard = Cspos_pos_hard;
	theCopy->Csneg_pos_hard = Csneg_pos_hard;
	theCopy->CEt_pos_hard = CEt_pos_hard;
	theCopy->CEc_pos_hard = CEc_pos_hard;
	theCopy->Cspos_neg_hard = Cspos_neg_hard;
	theCopy->Csneg_neg_hard = Csneg_neg_hard;
	theCopy->CEt_neg_hard = CEt_neg_hard;
	theCopy->CEc_neg_hard = CEc_neg_hard;
	theCopy->CEpos_pos_hard = CEpos_pos_hard;
	theCopy->CEneg_neg_hard = CEneg_neg_hard;
	theCopy->CEpos_hard = CEpos_hard;
	theCopy->CEneg_hard = CEneg_hard;
	theCopy->CEposnorm_hard = CEposnorm_hard;
	theCopy->CEnegnorm_hard = CEnegnorm_hard;
	theCopy->Cdelta_pos_hard = Cdelta_pos_hard;
	theCopy->Cdelta_pos_max_hard = Cdelta_pos_max_hard;
	theCopy->Cdelta_neg_hard = Cdelta_neg_hard;
	theCopy->Cdelta_neg_max_hard = Cdelta_neg_max_hard;
	theCopy->Calpha_pos = Calpha_pos;
	theCopy->Calpha_neg = Calpha_neg;

	//Hardening damage model Shear
	theCopy->Et_hard_Shr = Et_hard_Shr;
	theCopy->Ec_hard_Shr = Ec_hard_Shr;
	theCopy->spos_hard_Shr = spos_hard_Shr;
	theCopy->sneg_hard_Shr = sneg_hard_Shr;
	theCopy->spos_pos_hard_Shr = spos_pos_hard_Shr;
	theCopy->sneg_pos_hard_Shr = sneg_pos_hard_Shr;
	theCopy->Et_pos_hard_Shr = Et_pos_hard_Shr;
	theCopy->Ec_pos_hard_Shr = Ec_pos_hard_Shr;
	theCopy->spos_neg_hard_Shr = spos_neg_hard_Shr;
	theCopy->sneg_neg_hard_Shr = sneg_neg_hard_Shr;
	theCopy->Et_neg_hard_Shr = Et_neg_hard_Shr;
	theCopy->Ec_neg_hard_Shr = Ec_neg_hard_Shr;
	theCopy->Epos_pos_hard_Shr = Epos_pos_hard_Shr;
	theCopy->Eneg_neg_hard_Shr = Eneg_neg_hard_Shr;
	theCopy->Epos_hard_Shr = Epos_hard_Shr;
	theCopy->Eneg_hard_Shr = Eneg_hard_Shr;
	theCopy->Eposnorm_hard_Shr = Eposnorm_hard_Shr;
	theCopy->Enegnorm_hard_Shr = Enegnorm_hard_Shr;
	theCopy->delta_pos_hard_Shr = delta_pos_hard_Shr;
	theCopy->delta_pos_max_hard_Shr = delta_pos_max_hard_Shr;
	theCopy->delta_neg_hard_Shr = delta_neg_hard_Shr;
	theCopy->delta_neg_max_hard_Shr = delta_neg_max_hard_Shr;
	theCopy->alpha_pos_Shr = alpha_pos_Shr;
	theCopy->alpha_neg_Shr = alpha_neg_Shr;

	theCopy->Ec_neg = Ec_neg;
	theCopy->Ec_neg_Shr = Ec_neg_Shr;

	theCopy->Epos_pos = Epos_pos;
	theCopy->Epos_pos_Shr = Epos_pos_Shr;

	theCopy->Eneg_neg = Eneg_neg;
	theCopy->Eneg_neg_Shr = Eneg_neg_Shr;

	theCopy->T_area = T_area;

	theCopy->MtoRref = MtoRref;
	theCopy->MtoRref_Shr = MtoRref_Shr;

	theCopy->KrR = KrR;
	theCopy->KrR_Shr = KrR_Shr;

	theCopy->Krel = Krel;
	theCopy->Krel_Shr = Krel_Shr;

	theCopy->KuR = KuR;
	theCopy->KuR_Shr = KuR_Shr;

	theCopy->Kun = Kun;
	theCopy->Kun_Shr = Kun_Shr;

	theCopy->ER = ER;
	theCopy->ER_Shr = ER_Shr;

	theCopy->E = E;
	theCopy->E_Shr = E_Shr;

	theCopy->ld = ld;

	theCopy->Vcol = Vcol;
	theCopy->CVcol = CVcol;

	return theCopy;
}

int
GMG_CMAC2D::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	static Vector data(541);
	data(0) = this->getTag();

	// Monotonic Backbone Flexure Critical
	data(1) = Kepos_Rot;
	data(2) = Keneg_Rot;
	data(3) = fypos_Rot;
	data(4) = fyneg_Rot;
	data(5) = fcappos_Rot;
	data(6) = fcapneg_Rot;
	data(7) = dcappos_Rot;
	data(8) = dcapneg_Rot;
	data(9) = Kdegpos_Rot;
	data(10) = Kdegneg_rot;
	data(11) = frespos_Rot;
	data(12) = fresneg_Rot;
	data(13) = delUpos_Rot;
	data(14) = delUneg_Rot;
	data(15) = alpha_Er_Rot_Post_Yielding;
	data(16) = beta_Er_Rot_Post_Yielding;
	data(17) = alpha_Er_Rot_Post_Capping;
	data(18) = beta_Er_Rot_Post_Capping;
	data(19) = Er_Rot_Post_Yielding;
	data(20) = Er_Rot_Post_Capping;
	data(21) = Kun_Rot_Post_Yielding;
	data(22) = Kun_Rot_Post_Capping;
	data(23) = Kr_Rot_Post_Yielding;
	data(24) = Kr_Rot_Post_Capping;

	// Monotonic Backbone Splice Critical
	data(25) = fypos_splice_Rot;
	data(26) = fyneg_splice_Rot;
	data(27) = fypos_splice_deg_Rot;
	data(28) = fyneg_splice_deg_Rot;
	data(29) = fcappos_splice_Rot;
	data(30) = fcapneg_splice_Rot;
	data(31) = dcappos_splice_Rot;
	data(32) = dcapneg_splice_Rot;
	data(33) = Kdegpos_splice_Rot;
	data(34) = Kdegneg_splice_rot;
	data(35) = frespos_splice_Rot;
	data(36) = fresneg_splice_Rot;
	data(37) = delUpos_splice_Rot;
	data(38) = delUneg_splice_Rot;
	data(39) = alpha_Er_Rot_Post_Capping_splice;
	data(40) = beta_Er_Rot_Post_Capping_splice;
	data(41) = Er_Rot_Post_Capping_splice;
	data(42) = Kun_Rot_Post_Capping_splice;
	data(43) = Kr_Rot_Post_Capping_splice;
	data(44) = fs;
	data(45) = lb_deg;
	data(46) = fs_deg;

	// Monotonic Backbone Flexure-Shear Critical
	data(47) = dcappos_FlexShr_Rot;
	data(48) = dcapneg_FlexShr_Rot;
	data(49) = Kdegpos_FlexShr_Rot;
	data(50) = Kdegneg_FlexShr_rot;
	data(51) = frespos_FlexShr_Rot;
	data(52) = fresneg_FlexShr_Rot;
	data(53) = delUpos_FlexShr_Rot;
	data(54) = delUneg_FlexShr_Rot;
	data(55) = alpha_Er_Rot_Post_Capping_FlexShr;
	data(56) = beta_Er_Rot_Post_Capping_FlexShr;
	data(57) = Er_Rot_Post_Capping_FlexShr;
	data(58) = Kun_Rot_Post_Capping_FlexShr;
	data(59) = Kr_Rot_Post_Capping_FlexShr;

	// Damage Properties Rotation Flexure Critical
	data(60) = C1_Rot;
	data(61) = C2_Rot;
	data(62) = C3_Rot;
	data(63) = solpe_post_yielding_Rot;
	data(64) = solpe_post_capping_Rot;
	data(65) = delta_ratio_max_hard_Rot;
	data(66) = Ref_Energy_Coe_Rot;

	// Damage Properties Rotation Splice Critical
	data(67) = C1_splice_Rot;
	data(68) = C2_splice_Rot;
	data(69) = C3_splice_Rot;
	data(70) = solpe_post_capping_splice_Rot;

	// Damage Properties Rotation Flexure-Shear Critical
	data(71) = C1_FlexShr_Rot;
	data(72) = C2_FlexShr_Rot;
	data(73) = C3_FlexShr_Rot;
	data(74) = solpe_post_capping_FlexShr_Rot;

	// Converged Variables
	data(75) = Cstrain;
	data(76) = Cstress;
	data(77) = Cek;
	data(78) = CInCycFac;
	data(79) = Cflagdmg;
	data(80) = flagdmg_Hardening_strength;
	data(81) = Cflagdmg_Hardening_strength;
	data(82) = flagdmg_Hardening_strength_Shr;
	data(83) = Cflagdmg_Hardening_strength_Shr;
	data(84) = flagdmg_Hardening;
	data(85) = Cflagdmg_Hardening;
	data(86) = flagdmg_Hardening_Shr;
	data(87) = Cflagdmg_Hardening_Shr;
	data(88) = flag_entering_hardening_Flex_Rot;
	data(89) = Cflag_entering_hardening_Flex_Rot;
	data(90) = flag_entering_residual_Flex_Rot;
	data(91) = Cflag_entering_residual_Flex_Rot;
	data(92) = flag_entering_residual_Shr;
	data(93) = Cflag_entering_residual_Shr;

	data(94) = flag_entering_hardening_Flex_pos_Rot;
	data(95) = Cflag_entering_hardening_Flex_pos_Rot;
	data(96) = flag_entering_hardening_Flex_neg_Rot;
	data(97) = Cflag_entering_hardening_Flex_neg_Rot;

	data(98) = flag_entering_hardening_Splice_pos_Rot;
	data(99) = Cflag_entering_hardening_Splice_pos_Rot;
	data(100) = flag_entering_hardening_Splice_neg_Rot;
	data(101) = Cflag_entering_hardening_Splice_neg_Rot;

	data(102) = flag_entering_hardening_Shr_pos;
	data(103) = Cflag_entering_hardening_Shr_pos;
	data(104) = flag_entering_hardening_Shr_neg;
	data(105) = Cflag_entering_hardening_Shr_neg;

	data(106) = flag_entering_hardening_Splice_Rot;
	data(107) = Cflag_entering_hardening_Splice_Rot;
	data(108) = flag_entering_hardening_Shr;
	data(109) = Cflag_entering_hardening_Shr;

	data(110) = flag_FlexShr_Filure_Rot;
	data(111) = Cflag_FlexShr_Filure_Rot;
	data(112) = flag_Splice_Filure_Rot;
	data(113) = Cflag_Splice_Filure_Rot;
	data(114) = flag_FlexSplice_Filure_Rot;
	data(115) = Cflag_FlexSplice_Filure_Rot;
	data(116) = flag_Flx_Filure_Rot;
	data(117) = Cflag_Flx_Filure_Rot;
	data(118) = flag_Filure_Shr;
	data(119) = Cflag_Filure_Shr;

	data(120) = Cdpeakmax;
	data(121) = Cdpeakmin;
	data(122) = Cffmax;
	data(123) = Cffmin;

	data(124) = Cspos;
	data(125) = Csneg;

	data(126) = CdmgSpos;
	data(127) = CdmgSneg;
	data(128) = CEt;
	data(129) = CEc;
	data(130) = CEposnorm;
	data(131) = CEnegnorm;
	data(132) = Cdpeakmin_inner;
	data(133) = Cffmin_inner;
	data(134) = Cffmax_inner;
	data(135) = Cdpeakmax_inner;
	data(136) = Cdpeakmax_inner_inner;
	data(137) = Cffmax_inner_inner;
	data(138) = Cdpeakmin_inner_inner;
	data(139) = Cffmin_inner_inner;
	data(140) = CfouterNP_max;
	data(141) = CfouterPN_min;
	data(142) = CInCycFac_6;
	data(143) = CInCycFac_7;
	data(144) = CInCycFac_neg_6;
	data(145) = CInCycFac_neg_7;

	data(146) = CBenchmark56_up;
	data(147) = Benchmark56_up;
	data(148) = CBenchmark56_down;
	data(149) = Benchmark56_down;

	data(150) = CBenchmark_neg_56_up;
	data(151) = Benchmark_neg_56_up;
	data(152) = CBenchmark_neg_56_down;
	data(153) = Benchmark_neg_56_down;

	data(154) = CBenchmark67_up;
	data(155) = CBenchmark67_down;
	data(156) = CBenchmark_neg_67_up;
	data(157) = CBenchmark_neg_67_down;
	data(158) = Cspos_pos;
	data(159) = Csneg_pos;
	data(160) = CEt_pos;
	data(161) = CEc_pos;
	data(162) = Cspos_neg;
	data(163) = Csneg_neg;
	data(164) = CEt_neg;
	data(165) = CEc_neg;
	data(166) = CEpos_pos;
	data(167) = CEneg_neg;
	data(168) = CT_area;
	data(169) = Cek_Rot;
	data(170) = Cstress_Rot;

	data(171) = Cflagdmg_Shr;
	data(172) = flagdmg;

	data(173) = Kepos_Shr;
	data(174) = Keneg_Shr;
	data(175) = fypos_Shr;
	data(176) = fyneg_Shr;
	data(177) = fcappos_Shr;
	data(178) = fcapneg_Shr;
	data(179) = dcappos_Shr;
	data(180) = dcapneg_Shr;
	data(181) = Kdegpos_Shr;
	data(182) = Kdegneg_Shr;
	data(183) = frespos_Shr;
	data(184) = fresneg_Shr;
	data(184) = delUpos_Shr;
	data(185) = delUneg_Shr;

	// Damage Properties Shear
	data(187) = C1_Shr;
	data(188) = C2_Shr;
	data(189) = C3_Shr;
	data(190) = solpe_post_yielding_Shr;
	data(191) = solpe_post_capping_Shr;
	data(192) = delta_ratio_max_hard_Shr;
	data(193) = Ref_Energy_Coe_Shr;

	// Cyclic propertie Shear
	data(194) = Kun_Shr_Post_Yielding;
	data(195) = Kun_Shr_Post_Capping;
	data(196) = Kr_Shr_Post_Yielding;
	data(197) = Kr_Shr_Post_Capping;
	data(198) = Er_Shr_Post_Yielding;
	data(199) = Er_Shr_Post_Capping;
	data(200) = alpha_Er_Shr_Post_Yielding;
	data(201) = beta_Er_Shr_Post_Yielding;
	data(202) = alpha_Er_Shr_Post_Capping;
	data(203) = beta_Er_Shr_Post_Capping;

	// Axial parameter
	data(204) = Cstress_Axial;
	data(205) = Cstrain_Axial;

	// Converged Variables
	data(206) = Cstrain_Shr;

	data(207) = Cstress_Shr;
	data(208) = Cek_Shr;
	data(209) = CInCycFac_Shr;

	data(210) = Cflagdmg_Shr;
	data(211) = flagdmg_Shr;

	data(212) = Cdpeakmax_Shr;
	data(212) = Cdpeakmin_Shr;
	data(214) = Cffmax_Shr;
	data(215) = Cffmin_Shr;

	data(216) = Cspos_Shr;
	data(217) = Csneg_Shr;
	data(218) = CdmgSpos_Shr;
	data(219) = CdmgSneg_Shr;
	data(220) = CEt_Shr;
	data(221) = CEc_Shr;
	data(222) = CEposnorm_Shr;
	data(223) = CEnegnorm_Shr;
	data(224) = Cdpeakmin_inner_Shr;
	data(225) = Cffmin_inner_Shr;
	data(226) = Cffmax_inner_Shr;
	data(227) = Cdpeakmax_inner_Shr;
	data(228) = Cdpeakmax_inner_inner_Shr;
	data(229) = Cffmax_inner_inner_Shr;
	data(230) = Cdpeakmin_inner_inner_Shr;
	data(231) = Cffmin_inner_inner_Shr;
	data(232) = CfouterNP_max_Shr;
	data(233) = CfouterPN_min_Shr;
	data(234) = CInCycFac_6_Shr;
	data(235) = CInCycFac_7_Shr;
	data(236) = CInCycFac_neg_6_Shr;
	data(237) = CInCycFac_neg_7_Shr;

	data(238) = CBenchmark56_up_Shr;
	data(239) = Benchmark56_up_Shr;
	data(240) = CBenchmark56_down_Shr;
	data(241) = Benchmark56_down_Shr;

	data(242) = CBenchmark_neg_56_up_Shr;
	data(243) = Benchmark_neg_56_up_Shr;
	data(244) = CBenchmark_neg_56_down_Shr;
	data(245) = Benchmark_neg_56_down_Shr;

	data(246) = CBenchmark67_up_Shr;
	data(247) = CBenchmark67_down_Shr;
	data(248) = CBenchmark_neg_67_up_Shr;
	data(249) = CBenchmark_neg_67_down_Shr;
	data(250) = Cspos_pos_Shr;
	data(251) = Csneg_pos_Shr;
	data(252) = CEt_pos_Shr;
	data(253) = CEc_pos_Shr;
	data(254) = Cspos_neg_Shr;
	data(255) = Csneg_neg_Shr;
	data(256) = CEt_neg_Shr;
	data(257) = CEc_neg_Shr;
	data(258) = CEpos_pos_Shr;
	data(259) = CEneg_neg_Shr;
	data(260) = CT_area_Shr;

	data(261) = MtoRref;
	data(262) = MtoRref_Shr;
	data(263) = KrR;
	data(264) = Krel_Shr;
	data(265) = KuR;
	data(266) = KuR_Shr;
	data(267) = Kun;
	data(268) = KuR_Shr;
	data(269) = ER;
	data(270) = ER_Shr;
	data(271) = E;
	data(272) = E_Shr;
	data(273) = CMtoRref;
	data(274) = CMtoRref_Shr;
	data(275) = CKrR;
	data(276) = CKrel_Shr;
	data(277) = CKuR;
	data(278) = CKuR_Shr;
	data(279) = CKun;
	data(280) = CKuR_Shr;
	data(281) = CER;
	data(282) = CER_Shr;
	data(283) = CE;
	data(284) = CE_Shr;
	data(285) = KrR_Shr;
	data(286) = Krel;
	data(287) = CKrR_Shr;
	data(288) = CKrel;
	data(289) = CflagFlur_Rot;
	data(290) = CflagFlur_Shr;
	data(291) = flagFlur_Rot;
	data(292) = flagFlur_Shr;

	data(293) = Cf_Bench_Rot;
	data(294) = Cf_Bench_Shr;
	data(295) = f_Bench_Rot;
	data(296) = f_Bench_Shr;
	data(297) = Intcpt_res_pos;
	data(298) = CIntcpt_res_pos;
	data(299) = Intcpt_res_pos_Shr;
	data(300) = CIntcpt_res_pos_Shr;
	data(301) = Intcpt_deg_pos;
	data(302) = CIntcpt_deg_pos;
	data(303) = Intcpt_deg_pos_Shr;
	data(304) = CIntcpt_deg_pos_Shr;
	data(305) = Intcpt_res_neg;
	data(306) = CIntcpt_res_neg;
	data(307) = Intcpt_res_neg_Shr;
	data(308) = CIntcpt_res_neg_Shr;
	data(309) = Intcpt_Xaxis_pos;
	data(310) = CIntcpt_Xaxis_pos;
	data(311) = Intcpt_Xaxis_neg;
	data(312) = CIntcpt_Xaxis_neg;
	data(313) = Intcpt_Xaxis_pos_Shr;
	data(314) = CIntcpt_Xaxis_pos_Shr;
	data(315) = Intcpt_Xaxis_neg_Shr;
	data(316) = CIntcpt_Xaxis_neg_Shr;

	data(317) = Et_hard;
	data(318) = CEt_hard;
	data(319) = Ec_hard;
	data(320) = CEc_hard;
	data(321) = spos_hard;
	data(322) = Cspos_hard;
	data(323) = sneg_hard;
	data(324) = Csneg_hard;
	data(325) = spos_pos_hard;
	data(326) = Cspos_pos_hard;
	data(327) = sneg_pos_hard;
	data(328) = Csneg_pos_hard;
	data(329) = Et_pos_hard;
	data(330) = CEt_pos_hard;
	data(331) = Ec_pos_hard;
	data(332) = CEc_pos_hard;
	data(333) = spos_neg_hard;
	data(334) = Cspos_neg_hard;
	data(335) = sneg_neg_hard;
	data(336) = Csneg_neg_hard;
	data(337) = Et_neg_hard;
	data(338) = CEt_neg_hard;
	data(339) = Ec_neg_hard;
	data(340) = CEc_neg_hard;
	data(341) = Epos_pos_hard;
	data(342) = CEpos_pos_hard;
	data(343) = Eneg_neg_hard;
	data(344) = CEneg_neg_hard;
	data(345) = Epos_hard;
	data(346) = CEpos_hard;
	data(347) = Eneg_hard;
	data(348) = CEneg_hard;
	data(349) = Eposnorm_hard;
	data(350) = CEposnorm_hard;
	data(351) = Enegnorm_hard;
	data(352) = CEnegnorm_hard;
	data(353) = delta_pos_hard;
	data(354) = Cdelta_pos_hard;
	data(355) = delta_pos_max_hard;
	data(356) = Cdelta_pos_max_hard;
	data(357) = delta_neg_hard;
	data(358) = Cdelta_neg_hard;
	data(359) = delta_neg_max_hard;
	data(360) = Cdelta_neg_max_hard;
	data(361) = alpha_pos;
	data(362) = Calpha_pos;
	data(363) = alpha_neg;
	data(364) = Calpha_neg;

	//Hardening damage model Shear
	data(365) = Et_hard_Shr;
	data(366) = CEt_hard_Shr;
	data(367) = Ec_hard_Shr;
	data(368) = CEc_hard_Shr;
	data(369) = spos_hard_Shr;
	data(370) = Cspos_hard_Shr;
	data(371) = sneg_hard_Shr;
	data(372) = Csneg_hard_Shr;
	data(373) = spos_pos_hard_Shr;
	data(374) = Cspos_pos_hard_Shr;
	data(375) = sneg_pos_hard_Shr;
	data(376) = Csneg_pos_hard_Shr;
	data(377) = Et_pos_hard_Shr;
	data(378) = CEt_pos_hard_Shr;
	data(379) = Ec_pos_hard_Shr;
	data(380) = CEc_pos_hard_Shr;
	data(381) = spos_neg_hard_Shr;
	data(382) = Cspos_neg_hard_Shr;
	data(383) = sneg_neg_hard_Shr;
	data(384) = Csneg_neg_hard_Shr;
	data(385) = Et_neg_hard_Shr;
	data(386) = CEt_neg_hard_Shr;
	data(387) = Ec_neg_hard_Shr;
	data(388) = CEc_neg_hard_Shr;
	data(389) = Epos_pos_hard_Shr;
	data(390) = CEpos_pos_hard_Shr;
	data(391) = Eneg_neg_hard_Shr;
	data(392) = CEneg_neg_hard_Shr;
	data(393) = Epos_hard_Shr;
	data(394) = CEpos_hard_Shr;
	data(395) = Eneg_hard_Shr;
	data(396) = CEneg_hard_Shr;
	data(397) = Eposnorm_hard_Shr;
	data(398) = CEposnorm_hard_Shr;
	data(399) = Enegnorm_hard_Shr;
	data(400) = CEnegnorm_hard_Shr;
	data(401) = delta_pos_hard_Shr;
	data(402) = Cdelta_pos_hard_Shr;
	data(403) = delta_pos_max_hard_Shr;
	data(404) = Cdelta_pos_max_hard_Shr;
	data(405) = delta_neg_hard_Shr;
	data(406) = Cdelta_neg_hard_Shr;
	data(407) = delta_neg_max_hard_Shr;
	data(408) = Cdelta_neg_max_hard_Shr;
	data(409) = alpha_pos_Shr;
	data(410) = Calpha_pos_Shr;
	data(411) = alpha_neg_Shr;
	data(412) = Calpha_neg_Shr;
	data(413) = dpeakmax_bench;
	data(414) = Cdpeakmax_bench;
	data(415) = dpeakmax_bench_Shr;
	data(416) = Cdpeakmax_bench_Shr;
	data(417) = dpeakmin_bench;
	data(418) = Cdpeakmin_bench;
	data(419) = dpeakmin_bench_Shr;
	data(420) = Cdpeakmin_bench_Shr;

	data(421) = R_dcapneg;
	data(422) = CR_dcapneg;
	data(423) = R_dcapneg_Shr;
	data(424) = CR_dcapneg_Shr;
	data(425) = R_dcappos;
	data(426) = CR_dcappos;
	data(427) = R_dcappos_Shr;
	data(428) = CR_dcappos_Shr;
	data(429) = R_dypos;
	data(430) = R_dypos_Shr;
	data(431) = R_fypos;
	data(432) = R_fypos_Shr;
	data(433) = R_Kdegpos;
	data(434) = R_Kdegpos_Shr;
	data(435) = R_fyneg;
	data(436) = R_fyneg_Shr;
	data(437) = R_dyneg;
	data(438) = R_dyneg_Shr;
	data(439) = R_fcappos;
	data(440) = R_fcappos_Shr;
	data(441) = R_fcapneg;
	data(442) = R_fcapneg_Shr;
	data(443) = R_Kdegneg;
	data(444) = R_Kdegneg_Shr;
	data(445) = dppos;
	data(446) = dppos_Shr;
	data(447) = slope_pos;
	data(448) = slope_pos_Shr;
	data(449) = R_frespos;
	data(450) = R_frespos_Shr;
	data(451) = R_drespos;
	data(452) = R_drespos_Shr;
	data(453) = R_delUpos;
	data(454) = dpneg;
	data(455) = dpneg_Shr;
	data(456) = slope_neg;
	data(457) = slope_neg_Shr;
	data(458) = R_fresneg;
	data(459) = R_fresneg_Shr;
	data(460) = CR_drespos;
	data(461) = CR_drespos_Shr;
	data(462) = R_dresneg;
	data(463) = R_dresneg_Shr;
	data(464) = CR_dresneg;
	data(465) = CR_dresneg_Shr;
	data(466) = R_delUneg;
	data(467) = Intcpt_slope_neg;
	data(468) = Intcpt_slope_neg_Shr;
	data(469) = Intcpt_deg_neg;
	data(470) = CIntcpt_deg_neg;
	data(471) = Intcpt_deg_neg_Shr;
	data(472) = CIntcpt_deg_neg_Shr;
	data(473) = CR_dypos;
	data(474) = CR_dypos_Shr;
	data(475) = CR_fypos;
	data(476) = CR_fypos_Shr;
	data(477) = CR_fcappos;
	data(478) = CR_fcappos_Shr;
	data(479) = CR_Kdegpos;
	data(480) = CR_Kdegpos_Shr;
	data(481) = CR_fyneg;
	data(482) = CR_fyneg_Shr;
	data(483) = R_dyneg;
	data(484) = R_dyneg_Shr;
	data(485) = CR_fcapneg;
	data(486) = CR_fcapneg_Shr;
	data(487) = CR_Kdegneg;
	data(488) = CR_Kdegneg_Shr;
	data(489) = Cdppos;
	data(490) = Cdppos_Shr;
	data(491) = Cslope_pos;
	data(492) = Cslope_pos_Shr;
	data(493) = CR_frespos;
	data(494) = CR_frespos_Shr;
	data(495) = CR_delUpos;
	data(496) = Intcpt_slope_pos;
	data(497) = Intcpt_slope_pos_Shr;
	data(498) = CIntcpt_slope_pos;
	data(499) = CIntcpt_slope_pos_Shr;
	data(500) = Cdpneg;
	data(501) = Cdpneg_Shr;
	data(502) = Cslope_neg;
	data(503) = Cslope_neg_Shr;
	data(504) = CR_fresneg;
	data(505) = CR_fresneg_Shr;
	data(506) = CR_delUneg;
	data(507) = CIntcpt_slope_neg;
	data(508) = CIntcpt_slope_neg_Shr;
	data(509) = ek_Axi;
	data(510) = Cek_Axi;
	data(511) = H;
	data(512) = B;
	data(513) = C_C;
	data(514) = L;
	data(515) = La;
	data(516) = fc;
	data(517) = fyL;
	data(518) = dbL;
	data(519) = nL_EW;
	data(520) = nL_NS;
	data(521) = fyT;
	data(522) = dbT;
	data(523) = n_leg;
	data(524) = S;
	data(525) = lb;
	data(526) = ld;
	data(527) = R_delUpos_Shr;
	data(528) = CR_delUpos_Shr;
	data(529) = R_delUneg_Shr;
	data(530) = CR_delUneg_Shr;
	data(531) = Ckd;
	data(532) = Vcol;
	data(533) = CVcol;
	data(534) = P_Axial_average;
	data(535) = VyE_to_Vcol_Ratio;
	data(536) = CVyE_to_Vcol_Ratio;
	data(537) = Vcap_for_ER;
	data(538) = CVcap_for_ER;

	data(539) = PVM_Flag;
	data(540) = Retrofit_flag;

	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "GMG_CMAC2D::sendSelf() - failed to send data\n";
	return res;
}

int
GMG_CMAC2D::recvSelf(int cTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector data(541);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);

	if (res < 0) {
		opserr << "GMG_CMAC2D::recvSelf() - failed to receive data\n";
		this->setTag(0);
	}
	else {
		this->setTag((int)data(0));
		Kepos_Rot = data(1);
		Keneg_Rot = data(2);
		fypos_Rot = data(3);
		fyneg_Rot = data(4);
		fcappos_Rot = data(5);
		fcapneg_Rot = data(6);
		dcappos_Rot = data(7);
		dcapneg_Rot = data(8);
		Kdegpos_Rot = data(9);
		Kdegneg_rot = data(10);
		frespos_Rot = data(11);
		fresneg_Rot = data(12);
		delUpos_Rot = data(13);
		delUneg_Rot = data(14);
		alpha_Er_Rot_Post_Yielding = data(15);
		beta_Er_Rot_Post_Yielding = data(16);
		alpha_Er_Rot_Post_Capping = data(17);
		beta_Er_Rot_Post_Capping = data(18);
		Er_Rot_Post_Yielding = data(19);
		Er_Rot_Post_Capping = data(20);
		Kun_Rot_Post_Yielding = data(21);
		Kun_Rot_Post_Capping = data(22);
		Kr_Rot_Post_Yielding = data(23);
		Kr_Rot_Post_Capping = data(24);

		//	Monotonic	Backbone
		fypos_splice_Rot = data(25);
		fyneg_splice_Rot = data(26);
		fypos_splice_deg_Rot = data(27);
		fyneg_splice_deg_Rot = data(28);
		fcappos_splice_Rot = data(29);
		fcapneg_splice_Rot = data(30);
		dcappos_splice_Rot = data(31);
		dcapneg_splice_Rot = data(32);
		Kdegpos_splice_Rot = data(33);
		Kdegneg_splice_rot = data(34);
		frespos_splice_Rot = data(35);
		fresneg_splice_Rot = data(36);
		delUpos_splice_Rot = data(37);
		delUneg_splice_Rot = data(38);
		alpha_Er_Rot_Post_Capping_splice = data(39);
		beta_Er_Rot_Post_Capping_splice = data(40);
		Er_Rot_Post_Capping_splice = data(41);
		Kun_Rot_Post_Capping_splice = data(42);
		Kr_Rot_Post_Capping_splice = data(43);
		fs = data(44);
		lb_deg = data(45);
		fs_deg = data(46);

		//	Monotonic	Backbone
		dcappos_FlexShr_Rot = data(47);
		dcapneg_FlexShr_Rot = data(48);
		Kdegpos_FlexShr_Rot = data(49);
		Kdegneg_FlexShr_rot = data(50);
		frespos_FlexShr_Rot = data(51);
		fresneg_FlexShr_Rot = data(52);
		delUpos_FlexShr_Rot = data(53);
		delUneg_FlexShr_Rot = data(54);
		alpha_Er_Rot_Post_Capping_FlexShr = data(55);
		beta_Er_Rot_Post_Capping_FlexShr = data(56);
		Er_Rot_Post_Capping_FlexShr = data(57);
		Kun_Rot_Post_Capping_FlexShr = data(58);
		Kr_Rot_Post_Capping_FlexShr = data(59);

		//	Damage	Properties
		C1_Rot = data(60);
		C2_Rot = data(61);
		C3_Rot = data(62);
		solpe_post_yielding_Rot = data(63);
		solpe_post_capping_Rot = data(64);
		delta_ratio_max_hard_Rot = data(65);
		Ref_Energy_Coe_Rot = data(66);

		//	Damage	Properties
		C1_splice_Rot = data(67);
		C2_splice_Rot = data(68);
		C3_splice_Rot = data(69);
		solpe_post_capping_splice_Rot = data(70);

		//	Damage	Properties
		C1_FlexShr_Rot = data(71);
		C2_FlexShr_Rot = data(72);
		C3_FlexShr_Rot = data(73);
		solpe_post_capping_FlexShr_Rot = data(74);

		//	Converged	Variables
		Cstrain = data(75);
		Cstress = data(76);
		Cek = data(77);
		CInCycFac = data(78);
		Cflagdmg = data(79);
		flagdmg_Hardening_strength = data(80);
		Cflagdmg_Hardening_strength = data(81);
		flagdmg_Hardening_strength_Shr = data(82);
		Cflagdmg_Hardening_strength_Shr = data(83);
		flagdmg_Hardening = data(84);
		Cflagdmg_Hardening = data(85);
		flagdmg_Hardening_Shr = data(86);
		Cflagdmg_Hardening_Shr = data(87);
		flag_entering_hardening_Flex_Rot = data(88);
		Cflag_entering_hardening_Flex_Rot = data(89);
		flag_entering_residual_Flex_Rot = data(90);
		Cflag_entering_residual_Flex_Rot = data(91);
		flag_entering_residual_Shr = data(92);
		Cflag_entering_residual_Shr = data(93);

		flag_entering_hardening_Flex_pos_Rot = data(94);
		Cflag_entering_hardening_Flex_pos_Rot = data(95);
		flag_entering_hardening_Flex_neg_Rot = data(96);
		Cflag_entering_hardening_Flex_neg_Rot = data(97);

		flag_entering_hardening_Splice_pos_Rot = data(98);
		Cflag_entering_hardening_Splice_pos_Rot = data(99);
		flag_entering_hardening_Splice_neg_Rot = data(100);
		Cflag_entering_hardening_Splice_neg_Rot = data(101);

		flag_entering_hardening_Shr_pos = data(102);
		Cflag_entering_hardening_Shr_pos = data(103);
		flag_entering_hardening_Shr_neg = data(104);
		Cflag_entering_hardening_Shr_neg = data(105);

		flag_entering_hardening_Splice_Rot = data(106);
		Cflag_entering_hardening_Splice_Rot = data(107);
		flag_entering_hardening_Shr = data(108);
		Cflag_entering_hardening_Shr = data(109);

		flag_FlexShr_Filure_Rot = data(110);
		Cflag_FlexShr_Filure_Rot = data(111);
		flag_Splice_Filure_Rot = data(112);
		Cflag_Splice_Filure_Rot = data(113);
		flag_FlexSplice_Filure_Rot = data(114);
		Cflag_FlexSplice_Filure_Rot = data(115);
		flag_Flx_Filure_Rot = data(116);
		Cflag_Flx_Filure_Rot = data(117);
		flag_Filure_Shr = data(118);
		Cflag_Filure_Shr = data(119);

		Cdpeakmax = data(120);
		Cdpeakmin = data(121);
		Cffmax = data(122);
		Cffmin = data(123);

		Cspos = data(124);
		Csneg = data(125);

		CdmgSpos = data(126);
		CdmgSneg = data(127);
		CEt = data(128);
		CEc = data(129);
		CEposnorm = data(130);
		CEnegnorm = data(131);
		Cdpeakmin_inner = data(132);
		Cffmin_inner = data(133);
		Cffmax_inner = data(134);
		Cdpeakmax_inner = data(135);
		Cdpeakmax_inner_inner = data(136);
		Cffmax_inner_inner = data(137);
		Cdpeakmin_inner_inner = data(138);
		Cffmin_inner_inner = data(139);
		CfouterNP_max = data(140);
		CfouterPN_min = data(141);
		CInCycFac_6 = data(142);
		CInCycFac_7 = data(143);
		CInCycFac_neg_6 = data(144);
		CInCycFac_neg_7 = data(145);

		CBenchmark56_up = data(146);
		Benchmark56_up = data(147);
		CBenchmark56_down = data(148);
		Benchmark56_down = data(149);

		CBenchmark_neg_56_up = data(150);
		Benchmark_neg_56_up = data(151);
		CBenchmark_neg_56_down = data(152);
		Benchmark_neg_56_down = data(153);

		CBenchmark67_up = data(154);
		CBenchmark67_down = data(155);
		CBenchmark_neg_67_up = data(156);
		CBenchmark_neg_67_down = data(157);
		Cspos_pos = data(158);
		Csneg_pos = data(159);
		CEt_pos = data(160);
		CEc_pos = data(161);
		Cspos_neg = data(162);
		Csneg_neg = data(163);
		CEt_neg = data(164);
		CEc_neg = data(165);
		CEpos_pos = data(166);
		CEneg_neg = data(167);
		CT_area = data(168);
		Cek_Rot = data(169);
		Cstress_Rot = data(170);

		Cflagdmg_Shr = data(171);
		flagdmg = data(172);

		Kepos_Shr = data(173);
		Keneg_Shr = data(174);
		fypos_Shr = data(175);
		fyneg_Shr = data(176);
		fcappos_Shr = data(177);
		fcapneg_Shr = data(178);
		dcappos_Shr = data(179);
		dcapneg_Shr = data(180);
		Kdegpos_Shr = data(181);
		Kdegneg_Shr = data(182);
		frespos_Shr = data(183);
		fresneg_Shr = data(184);
		delUpos_Shr = data(184);
		delUneg_Shr = data(185);

		//	Damage	Properties
		C1_Shr = data(187);
		C2_Shr = data(188);
		C3_Shr = data(189);
		solpe_post_yielding_Shr = data(190);
		solpe_post_capping_Shr = data(191);
		delta_ratio_max_hard_Shr = data(192);
		Ref_Energy_Coe_Shr = data(193);

		//	Cyclic	propertie
		Kun_Shr_Post_Yielding = data(194);
		Kun_Shr_Post_Capping = data(195);
		Kr_Shr_Post_Yielding = data(196);
		Kr_Shr_Post_Capping = data(197);
		Er_Shr_Post_Yielding = data(198);
		Er_Shr_Post_Capping = data(199);
		alpha_Er_Shr_Post_Yielding = data(200);
		beta_Er_Shr_Post_Yielding = data(201);
		alpha_Er_Shr_Post_Capping = data(202);
		beta_Er_Shr_Post_Capping = data(203);

		//	Axial parameter
		Cstress_Axial = data(204);
		Cstrain_Axial = data(205);

		// Converged Variables
		Cstrain_Shr = data(206);

		Cstress_Shr = data(207);
		Cek_Shr = data(208);
		CInCycFac_Shr = data(209);

		Cflagdmg_Shr = data(210);
		flagdmg_Shr = data(211);

		Cdpeakmax_Shr = data(212);
		Cdpeakmin_Shr = data(212);
		Cffmax_Shr = data(214);
		Cffmin_Shr = data(215);

		Cspos_Shr = data(216);
		Csneg_Shr = data(217);
		CdmgSpos_Shr = data(218);
		CdmgSneg_Shr = data(219);
		CEt_Shr = data(220);
		CEc_Shr = data(221);
		CEposnorm_Shr = data(222);
		CEnegnorm_Shr = data(223);
		Cdpeakmin_inner_Shr = data(224);
		Cffmin_inner_Shr = data(225);
		Cffmax_inner_Shr = data(226);
		Cdpeakmax_inner_Shr = data(227);
		Cdpeakmax_inner_inner_Shr = data(228);
		Cffmax_inner_inner_Shr = data(229);
		Cdpeakmin_inner_inner_Shr = data(230);
		Cffmin_inner_inner_Shr = data(231);
		CfouterNP_max_Shr = data(232);
		CfouterPN_min_Shr = data(233);
		CInCycFac_6_Shr = data(234);
		CInCycFac_7_Shr = data(235);
		CInCycFac_neg_6_Shr = data(236);
		CInCycFac_neg_7_Shr = data(237);

		CBenchmark56_up_Shr = data(238);
		Benchmark56_up_Shr = data(239);
		CBenchmark56_down_Shr = data(240);
		Benchmark56_down_Shr = data(241);

		CBenchmark_neg_56_up_Shr = data(242);
		Benchmark_neg_56_up_Shr = data(243);
		CBenchmark_neg_56_down_Shr = data(244);
		Benchmark_neg_56_down_Shr = data(245);

		CBenchmark67_up_Shr = data(246);
		CBenchmark67_down_Shr = data(247);
		CBenchmark_neg_67_up_Shr = data(248);
		CBenchmark_neg_67_down_Shr = data(249);
		Cspos_pos_Shr = data(250);
		Csneg_pos_Shr = data(251);
		CEt_pos_Shr = data(252);
		CEc_pos_Shr = data(253);
		Cspos_neg_Shr = data(254);
		Csneg_neg_Shr = data(255);
		CEt_neg_Shr = data(256);
		CEc_neg_Shr = data(257);
		CEpos_pos_Shr = data(258);
		CEneg_neg_Shr = data(259);
		CT_area_Shr = data(260);

		MtoRref = data(261);
		MtoRref_Shr = data(262);
		KrR = data(263);
		Krel_Shr = data(264);
		KuR = data(265);
		KuR_Shr = data(266);
		Kun = data(267);
		KuR_Shr = data(268);
		ER = data(269);
		ER_Shr = data(270);
		E = data(271);
		E_Shr = data(272);
		CMtoRref = data(273);
		CMtoRref_Shr = data(274);
		CKrR = data(275);
		CKrel_Shr = data(276);
		CKuR = data(277);
		CKuR_Shr = data(278);
		CKun = data(279);
		CKuR_Shr = data(280);
		CER = data(281);
		CER_Shr = data(282);
		CE = data(283);
		CE_Shr = data(284);
		KrR_Shr = data(285);
		Krel = data(286);
		CKrR_Shr = data(287);
		CKrel = data(288);
		CflagFlur_Rot = data(289);
		CflagFlur_Shr = data(290);
		flagFlur_Rot = data(291);
		flagFlur_Shr = data(292);

		Cf_Bench_Rot = data(293);
		Cf_Bench_Shr = data(294);
		f_Bench_Rot = data(295);
		f_Bench_Shr = data(296);
		Intcpt_res_pos = data(297);
		CIntcpt_res_pos = data(298);
		Intcpt_res_pos_Shr = data(299);
		CIntcpt_res_pos_Shr = data(300);
		Intcpt_deg_pos = data(301);
		CIntcpt_deg_pos = data(302);
		Intcpt_deg_pos_Shr = data(303);
		CIntcpt_deg_pos_Shr = data(304);
		Intcpt_res_neg = data(305);
		CIntcpt_res_neg = data(306);
		Intcpt_res_neg_Shr = data(307);
		CIntcpt_res_neg_Shr = data(308);
		Intcpt_Xaxis_pos = data(309);
		CIntcpt_Xaxis_pos = data(310);
		Intcpt_Xaxis_neg = data(311);
		CIntcpt_Xaxis_neg = data(312);
		Intcpt_Xaxis_pos_Shr = data(313);
		CIntcpt_Xaxis_pos_Shr = data(314);
		Intcpt_Xaxis_neg_Shr = data(315);
		CIntcpt_Xaxis_neg_Shr = data(316);

		Et_hard = data(317);
		CEt_hard = data(318);
		Ec_hard = data(319);
		CEc_hard = data(320);
		spos_hard = data(321);
		Cspos_hard = data(322);
		sneg_hard = data(323);
		Csneg_hard = data(324);
		spos_pos_hard = data(325);
		Cspos_pos_hard = data(326);
		sneg_pos_hard = data(327);
		Csneg_pos_hard = data(328);
		Et_pos_hard = data(329);
		CEt_pos_hard = data(330);
		Ec_pos_hard = data(331);
		CEc_pos_hard = data(332);
		spos_neg_hard = data(333);
		Cspos_neg_hard = data(334);
		sneg_neg_hard = data(335);
		Csneg_neg_hard = data(336);
		Et_neg_hard = data(337);
		CEt_neg_hard = data(338);
		Ec_neg_hard = data(339);
		CEc_neg_hard = data(340);
		Epos_pos_hard = data(341);
		CEpos_pos_hard = data(342);
		Eneg_neg_hard = data(343);
		CEneg_neg_hard = data(344);
		Epos_hard = data(345);
		CEpos_hard = data(346);
		Eneg_hard = data(347);
		CEneg_hard = data(348);
		Eposnorm_hard = data(349);
		CEposnorm_hard = data(350);
		Enegnorm_hard = data(351);
		CEnegnorm_hard = data(352);
		delta_pos_hard = data(353);
		Cdelta_pos_hard = data(354);
		delta_pos_max_hard = data(355);
		Cdelta_pos_max_hard = data(356);
		delta_neg_hard = data(357);
		Cdelta_neg_hard = data(358);
		delta_neg_max_hard = data(359);
		Cdelta_neg_max_hard = data(360);
		alpha_pos = data(361);
		Calpha_pos = data(362);
		alpha_neg = data(363);
		Calpha_neg = data(364);

		//Hardening	damage model
		Et_hard_Shr = data(365);
		CEt_hard_Shr = data(366);
		Ec_hard_Shr = data(367);
		CEc_hard_Shr = data(368);
		spos_hard_Shr = data(369);
		Cspos_hard_Shr = data(370);
		sneg_hard_Shr = data(371);
		Csneg_hard_Shr = data(372);
		spos_pos_hard_Shr = data(373);
		Cspos_pos_hard_Shr = data(374);
		sneg_pos_hard_Shr = data(375);
		Csneg_pos_hard_Shr = data(376);
		Et_pos_hard_Shr = data(377);
		CEt_pos_hard_Shr = data(378);
		Ec_pos_hard_Shr = data(379);
		CEc_pos_hard_Shr = data(380);
		spos_neg_hard_Shr = data(381);
		Cspos_neg_hard_Shr = data(382);
		sneg_neg_hard_Shr = data(383);
		Csneg_neg_hard_Shr = data(384);
		Et_neg_hard_Shr = data(385);
		CEt_neg_hard_Shr = data(386);
		Ec_neg_hard_Shr = data(387);
		CEc_neg_hard_Shr = data(388);
		Epos_pos_hard_Shr = data(389);
		CEpos_pos_hard_Shr = data(390);
		Eneg_neg_hard_Shr = data(391);
		CEneg_neg_hard_Shr = data(392);
		Epos_hard_Shr = data(393);
		CEpos_hard_Shr = data(394);
		Eneg_hard_Shr = data(395);
		CEneg_hard_Shr = data(396);
		Eposnorm_hard_Shr = data(397);
		CEposnorm_hard_Shr = data(398);
		Enegnorm_hard_Shr = data(399);
		CEnegnorm_hard_Shr = data(400);
		delta_pos_hard_Shr = data(401);
		Cdelta_pos_hard_Shr = data(402);
		delta_pos_max_hard_Shr = data(403);
		Cdelta_pos_max_hard_Shr = data(404);
		delta_neg_hard_Shr = data(405);
		Cdelta_neg_hard_Shr = data(406);
		delta_neg_max_hard_Shr = data(407);
		Cdelta_neg_max_hard_Shr = data(408);
		alpha_pos_Shr = data(409);
		Calpha_pos_Shr = data(410);
		alpha_neg_Shr = data(411);
		Calpha_neg_Shr = data(412);
		dpeakmax_bench = data(413);
		Cdpeakmax_bench = data(414);
		dpeakmax_bench_Shr = data(415);
		Cdpeakmax_bench_Shr = data(416);
		dpeakmin_bench = data(417);
		Cdpeakmin_bench = data(418);
		dpeakmin_bench_Shr = data(419);
		Cdpeakmin_bench_Shr = data(420);

		R_dcapneg = data(421);
		CR_dcapneg = data(422);
		R_dcapneg_Shr = data(423);
		CR_dcapneg_Shr = data(424);
		R_dcappos = data(425);
		CR_dcappos = data(426);
		R_dcappos_Shr = data(427);
		CR_dcappos_Shr = data(428);
		R_dypos = data(429);
		R_dypos_Shr = data(430);
		R_fypos = data(431);
		R_fypos_Shr = data(432);
		R_Kdegpos = data(433);
		R_Kdegpos_Shr = data(434);
		R_fyneg = data(435);
		R_fyneg_Shr = data(436);
		R_dyneg = data(437);
		R_dyneg_Shr = data(438);
		R_fcappos = data(439);
		R_fcappos_Shr = data(440);
		R_fcapneg = data(441);
		R_fcapneg_Shr = data(442);
		R_Kdegneg = data(443);
		R_Kdegneg_Shr = data(444);
		dppos = data(445);
		dppos_Shr = data(446);
		slope_pos = data(447);
		slope_pos_Shr = data(448);
		R_frespos = data(449);
		R_frespos_Shr = data(450);
		R_drespos = data(451);
		R_drespos_Shr = data(452);
		R_delUpos = data(453);
		dpneg = data(454);
		dpneg_Shr = data(455);
		slope_neg = data(456);
		slope_neg_Shr = data(457);
		R_fresneg = data(458);
		R_fresneg_Shr = data(459);
		CR_drespos = data(460);
		CR_drespos_Shr = data(461);
		R_dresneg = data(462);
		R_dresneg_Shr = data(463);
		CR_dresneg = data(464);
		CR_dresneg_Shr = data(465);
		R_delUneg = data(466);
		Intcpt_slope_neg = data(467);
		Intcpt_slope_neg_Shr = data(468);
		Intcpt_deg_neg = data(469);
		CIntcpt_deg_neg = data(470);
		Intcpt_deg_neg_Shr = data(471);
		CIntcpt_deg_neg_Shr = data(472);
		CR_dypos = data(473);
		CR_dypos_Shr = data(474);
		CR_fypos = data(475);
		CR_fypos_Shr = data(476);
		CR_fcappos = data(477);
		CR_fcappos_Shr = data(478);
		CR_Kdegpos = data(479);
		CR_Kdegpos_Shr = data(480);
		CR_fyneg = data(481);
		CR_fyneg_Shr = data(482);
		R_dyneg = data(483);
		R_dyneg_Shr = data(484);
		CR_fcapneg = data(485);
		CR_fcapneg_Shr = data(486);
		CR_Kdegneg = data(487);
		CR_Kdegneg_Shr = data(488);
		Cdppos = data(489);
		Cdppos_Shr = data(490);
		Cslope_pos = data(491);
		Cslope_pos_Shr = data(492);
		CR_frespos = data(493);
		CR_frespos_Shr = data(494);
		CR_delUpos = data(495);
		Intcpt_slope_pos = data(496);
		Intcpt_slope_pos_Shr = data(497);
		CIntcpt_slope_pos = data(498);
		CIntcpt_slope_pos_Shr = data(499);
		Cdpneg = data(500);
		Cdpneg_Shr = data(501);
		Cslope_neg = data(502);
		Cslope_neg_Shr = data(503);
		CR_fresneg = data(504);
		CR_fresneg_Shr = data(505);
		CR_delUneg = data(506);
		CIntcpt_slope_neg = data(507);
		CIntcpt_slope_neg_Shr = data(508);
		ek_Axi = data(509);
		Cek_Axi = data(510);
		H = data(511);
		B = data(512);
		C_C = data(513);
		L = data(514);
		La = data(515);
		fc = data(516);
		fyL = data(517);
		dbL = data(518);
		nL_EW = data(519);
		nL_NS = data(520);
		fyT = data(521);
		dbT = data(522);
		n_leg = data(523);
		S = data(524);
		lb = data(525);
		ld = data(526);
		R_delUpos_Shr = data(527);
		CR_delUpos_Shr = data(528);
		R_delUneg_Shr = data(529);
		CR_delUneg_Shr = data(530);
		Ckd = data(531);
		Vcol = data(532);
		CVcol = data(533);
		P_Axial_average = data(534);
		VyE_to_Vcol_Ratio = data(535);
		CVyE_to_Vcol_Ratio = data(536);
		Vcap_for_ER = data(537);
		CVcap_for_ER = data(538);

		PVM_Flag = data(539);
		Retrofit_flag = data(540);
	}

	return res;
}

void
GMG_CMAC2D::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "GMG_CMAC2D tag: " << this->getTag() << endln;
		s << "  Kepos_Rot: " << Kepos_Rot << endln;
		s << "  Keneg_Rot: " << Keneg_Rot << endln;
		s << "  fypos_Rot: " << fypos_Rot << endln;
		s << "  fyneg_Rot: " << fyneg_Rot << endln;
		s << "  fcappos_Rot: " << fcappos_Rot << endln;
		s << "  fcapneg_Rot: " << fcapneg_Rot << endln;
		s << "  dcappos_Rot: " << dcappos_Rot << endln;
		s << "  dcapneg_Rot: " << dcapneg_Rot << endln;
		s << "  Kdegpos_Rot: " << Kdegpos_Rot << endln;
		s << "  Kdegneg_rot: " << Kdegneg_rot << endln;
		s << "  frespos_Rot: " << frespos_Rot << endln;
		s << "  fresneg_Rot: " << fresneg_Rot << endln;
		s << "  delUpos_Rot: " << delUpos_Rot << endln;
		s << "  delUneg_Rot: " << delUneg_Rot << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"GMG_CMAC2D\", ";
		s << "\"Kepos_Rot\": " << Kepos_Rot << ", ";
		s << "\"Keneg_Rot\": " << Keneg_Rot << ", ";
		s << "\"fypos_Rot\": " << fypos_Rot << ", ";
		s << "\"fyneg_Rot\": " << fyneg_Rot << ", ";
		s << "\"fcappos_Rot\": " << fcappos_Rot << ", ";
		s << "\"fcapneg_Rot\": " << fcapneg_Rot << ", ";
		s << "\"dcappos_Rot\": " << dcappos_Rot << ", ";
		s << "\"dcapneg_Rot\": " << dcapneg_Rot << ", ";
		s << "\"Kdegpos_Rot\": " << Kdegpos_Rot << ", ";
		s << "\"Kdegneg_rot\": " << Kdegneg_rot << ", ";
		s << "\"frespos_Rot\": " << frespos_Rot << ", ";
		s << "\"fresneg_Rot\": " << fresneg_Rot << ", ";
		s << "\"delUpos_Rot\": " << delUpos_Rot << ", ";
		s << "\"delUneg_Rot\": " << delUneg_Rot << ", ";
	}
}


//-----------------------------My subrotin----------------------------------
/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Rotational Material                                    **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */

void
GMG_CMAC2D::Needed_Inpt_to_get_failure_mode(void) { //Capping, 


	const double PI = 3.141592653589793238463;

	//Transverse reinf. calculation
	double AsV = PI * pow(dbT, 2.0) / 4.0;
	double roT = (n_leg*AsV) / (B*S);

	//Longitudinal reinf. calculation
	double AsL = PI * pow(dbL, 2.0) / 4.0;
	double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);

	// Geometry
	double deff_EW = 0.8*H;

	if (d_Axial > 0) {
		if (flag_entering_hardening_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {
			double F_Axial_mod;
			if (Kepos_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 1) < 0) {
				if (TstateFlag != 0)
					F_Axial_mod = max(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 1) - 0.01*fc*B*H);
				else
					F_Axial_mod = Kepos_Axil * d_Axial;
			}

			else {
				if (TstateFlag != 0)
					F_Axial_mod = min(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 1) + 0.01*fc*B*H);
				else
					F_Axial_mod = Kepos_Axil * d_Axial;
			}

			if (TstateFlag != 0) {
				P_Axial_average = ((P_Axial_sum / (double)number_of_averaged_elements_Hardening)*((double)number_of_averaged_elements_Hardening - 1.0) + F_Axial_mod) / (double)number_of_averaged_elements_Hardening;
			}
			else {
				P_Axial_average = ((P_Axial_sum / (double)number_of_averaged_elements_Elastic)*((double)number_of_averaged_elements_Elastic - 1.0) + F_Axial_mod) / (double)number_of_averaged_elements_Elastic;
			}

		}

		else if (flag_Filure_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {
			double F_Axial_mod;
			if (Kepos_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 1) < 0) {
				if (TstateFlag_Shr != 0)
					F_Axial_mod = max(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 1) - 0.01*fc*B*H);
				else
					F_Axial_mod = Kepos_Axil * d_Axial;
			}

			else {
				if (TstateFlag_Shr != 0)
					F_Axial_mod = min(Kepos_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 1) + 0.01*fc*B*H);
				else
					F_Axial_mod = Kepos_Axil * d_Axial;
			}

			if (TstateFlag_Shr != 0) {
				P_Axial_average = ((P_Axial_sum / (double)number_of_averaged_elements_Hardening)*((double)number_of_averaged_elements_Hardening - 1.0) + F_Axial_mod) / (double)number_of_averaged_elements_Hardening;
			}
			else {
				P_Axial_average = ((P_Axial_sum / (double)number_of_averaged_elements_Elastic)*((double)number_of_averaged_elements_Elastic - 1.0) + F_Axial_mod) / (double)number_of_averaged_elements_Elastic;
			}

		}
	}


	else if (d_Axial <= 0) {
		if (flag_entering_hardening_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {
			double F_Axial_mod;
			if (Keneg_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 1) < 0) {
				if (TstateFlag != 0)
					F_Axial_mod = max(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 1) - 0.01*fc*B*H);
				else
					F_Axial_mod = Keneg_Axil * d_Axial;
			}

			else {
				if (TstateFlag != 0)
					F_Axial_mod = min(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 1) + 0.01*fc*B*H);
				else
					F_Axial_mod = Keneg_Axil * d_Axial;
			}

			if (TstateFlag != 0) {
				P_Axial_average = ((P_Axial_sum / (double)number_of_averaged_elements_Hardening)*((double)number_of_averaged_elements_Hardening - 1.0) + F_Axial_mod) / (double)number_of_averaged_elements_Hardening;
			}
			else {
				P_Axial_average = ((P_Axial_sum / (double)number_of_averaged_elements_Elastic)*((double)number_of_averaged_elements_Elastic - 1.0) + F_Axial_mod) / (double)number_of_averaged_elements_Elastic;
			}

		}

		else if (flag_Filure_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {
			double F_Axial_mod;
			if (Keneg_Axil * d_Axial - P_Axial_Vector(number_of_averaged_elements_Hardening - 1) < 0) {
				if (TstateFlag_Shr != 0)
					F_Axial_mod = max(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 1) - 0.01*fc*B*H);
				else
					F_Axial_mod = Keneg_Axil * d_Axial;
			}

			else {
				if (TstateFlag_Shr != 0)
					F_Axial_mod = min(Keneg_Axil * d_Axial, P_Axial_Vector(number_of_averaged_elements_Hardening - 1) + 0.01*fc*B*H);
				else
					F_Axial_mod = Keneg_Axil * d_Axial;
			}

			if (TstateFlag_Shr != 0) {
				P_Axial_average = ((P_Axial_sum / (double)number_of_averaged_elements_Hardening)*((double)number_of_averaged_elements_Hardening - 1.0) + F_Axial_mod) / (double)number_of_averaged_elements_Hardening;
			}
			else {
				P_Axial_average = ((P_Axial_sum / (double)number_of_averaged_elements_Elastic)*((double)number_of_averaged_elements_Elastic - 1.0) + F_Axial_mod) / (double)number_of_averaged_elements_Elastic;
			}

		}
	}

	if (Retrofit_flag == 0) {
		if (flag_Filure_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {

			if (PVM_Flag == 1) {
				P_M_yield_Opt(P_Axial_average);
				P_M_capping_Opt(P_Axial_average);
				P_M_splice_Opt(P_Axial_average);
				Shear_strength_Opt(P_Axial_average);
			}

			if (PVM_Flag != 1) {
				P_M_yield_Opt(-1.0*P_Axial);
				P_M_capping_Opt(-1.0*P_Axial);
				P_M_splice_Opt(-1.0*P_Axial);
				Shear_strength_Opt(-1.0*P_Axial);
			}

			if (VyE_to_Vcol_Ratio == 0.0) {
				//Flexure-Critical Check============================================================================================ 
				dypos = fypos_Rot / Kepos_Rot;
				dyneg = fyneg_Rot / Keneg_Rot;
				dcappos_Rot = fypos_Rot / Kepos_Rot + fmax(-0.02 + 0.054*((roT*fyT / fc) / (0.2 - (8.0 / 7.0)*(P_Axial_average / (B*H*fc)))) + 0.053 / pow(S / fmin(0.125*B, 0.125*H), 0.5), 0.1*dypos); // I had to use different sign for the Axial load here
				dcapneg_Rot = -dcappos_Rot;

				//Flexure-Shear Critical Check======================================================================================
				dcappos_FlexShr_Rot = fypos_Rot / Kepos_Rot + fmax(0.008 + 0.036*P_Axial_average / (B*H*fc) + 2.05*roT + 0.0055*pow(fypos_Rot / La / Vcol, -0.9), 0.1*dypos);
				dcapneg_FlexShr_Rot = -dcappos_FlexShr_Rot;
			}
			else {
				if (VyE_to_Vcol_Ratio <= 0.6) {
					//Flexure-Critical Check============================================================================================ 
					dypos = fypos_Rot / Kepos_Rot;
					dyneg = fyneg_Rot / Keneg_Rot;
					//dcappos_Rot = fypos_Rot / Kepos_Rot + fmax((0.025 / pow(S / fmin(B, H), 0.5) + 0.0035*(P_Axial_average) / (B*H*fc) - 0.02), 0.1*dypos); // I had to use different sign for the Axial load here
					dcappos_Rot = fypos_Rot / Kepos_Rot + fmax(-0.02 + 0.054*((roT*fyT / fc) / (0.2 - (8.0 / 7.0)*(P_Axial_average / (B*H*fc)))) + 0.053 / pow(S / fmin(0.125*B, 0.125*H), 0.5), 0.1*dypos); // I had to use different sign for the Axial load here
					dcapneg_Rot = -dcappos_Rot;

					//Removing Flexure-Shear Critical Check======================================================================================
					dcappos_FlexShr_Rot = 1000.0*(fypos_Rot / Kepos_Rot + fmax(0.008 + 0.036*P_Axial_average / (B*H*fc) + 2.05*roT + 0.0055*pow(fypos_Rot / La / Vcol, -0.9), 0.1*dypos));
					dcapneg_FlexShr_Rot = -dcappos_FlexShr_Rot;
				}
				else {
					//Removing Flexure-Critical Check============================================================================================ 
					dypos = fypos_Rot / Kepos_Rot;
					dyneg = fyneg_Rot / Keneg_Rot;
					dcappos_Rot = fypos_Rot / Kepos_Rot + (fmax(-0.02 + 0.054*((roT*fyT / fc) / (0.2 - (8.0 / 7.0)*(P_Axial_average / (B*H*fc)))) + 0.053 / pow(S / fmin(0.125*B, 0.125*H), 0.5), 0.1*dypos)); // I had to use different sign for the Axial load here

					dcapneg_Rot = -dcappos_Rot;
					fcappos_Rot = fcappos_Rot;
					fcapneg_Rot = fcapneg_Rot;

					//Flexure-Shear Critical Check======================================================================================
					dcappos_FlexShr_Rot = fypos_Rot / Kepos_Rot + fmax(0.008 + 0.036*P_Axial_average / (B*H*fc) + 2.05*roT + 0.0055*pow(fypos_Rot / La / Vcol, -0.9), 0.1*dypos);
					dcapneg_FlexShr_Rot = -dcappos_FlexShr_Rot;
				}
			}
			// Flexure-splice Check=============================================================================================
			if (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL) {
				dypos_splice = 1.000001*fypos_Rot / Kepos_Rot;
				dyneg_splice = 1.000001*fyneg_Rot / Keneg_Rot;
				dcappos_splice_Rot = fypos_Rot / Kepos_Rot + fmax(fmin(0.016 + 0.03*roT*fyT / (roL*fyL), 0.03), 0.1*dypos);
				dcapneg_splice_Rot = -dcappos_splice_Rot;
			}
			else {
				dypos_splice = 1000.0 * fypos_Rot / Kepos_Rot;
				dyneg_splice = 1000.0 * fyneg_Rot / Keneg_Rot;
				dcappos_splice_Rot = 1000.0 * fmax(dcappos_Rot, dcappos_FlexShr_Rot);
				dcapneg_splice_Rot = -dcappos_splice_Rot;
			}

			// Shear-Critical Check============================================================================================
			fypos_Shr = Vcol;
			dypos_Shr = Vcol / Kepos_Shr;
			dcappos_Shr = dypos_Shr + max((0.0025 + 1.4*roT)*La, /*0.001*/ 0.1*dypos_Shr); // I had to use different sign for the Axial load here
			fyneg_Shr = -Vcol;
			dyneg_Shr = -Vcol / Kepos_Shr;
			dcapneg_Shr = -dcappos_Shr;

			//Splice check=====================================================================================================
			if (lb != 0.0 && lb < ld && fs < fyL) {
				dypos_splice = fypos_splice_Rot / Kepos_Rot;
				dyneg_splice = fyneg_splice_Rot / Keneg_Rot;
				dcappos_splice_Rot = fypos_splice_Rot / Kepos_Rot + fmax(fmin(0.016 + 0.03*roT*fyT / (roL*fyL), 0.03), 0.1*dypos_splice);
				dcapneg_splice_Rot = -dcappos_splice_Rot;
			}

			defineBackbone();
			defineBackbone_Shr();
		}
	}


	if (Retrofit_flag == 1) { // Fiber jacket
		if (flag_Filure_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {

			double roj = (2.0*(B + H)*tj) / (B*H);

			double psi_t = 1.0;
			double psi_e = 1.0;
			double psi_s;
			if (dbL <= 0.8)
				psi_s = 0.8;
			else
				psi_s = 1.0;

			double psi_g;
			if (fyL <= 80.0)                       // For grade 60
				psi_g = 1.0;
			else if (fyL <= 100.0 && fyL > 80.0)   // For grade 80
				psi_g = 1.15;
			else                                   // For grade 100
				psi_g = 1.3;

			double cb = fmin(C_C + dbT + dbL / 2.0, (B - 2.0 * (C_C + dbT + dbL / 2.0) / (nL_EW - 1.0)) / 2.0);

			double AsV = PI * pow(dbT, 2.0) / 4.0;
			double Atr = n_leg * AsV;
			double Ktr = 40.0*(Atr) / (S*nL_EW);

			ld = (3.0 / 40.0) * (fyL * 1000.0 * psi_t * psi_e * psi_s * psi_g / (pow(fc*1000.0, 0.5) * fmin((cb + Ktr) / dbL, 2.5)))*dbL;

			if (PVM_Flag == 1) {
				P_M_yield_Opt(P_Axial_average);
				P_M_capping_Opt(P_Axial_average);
				P_M_splice_Opt(P_Axial_average);
				Shear_strength_Opt(P_Axial_average);
			}

			if (PVM_Flag != 1) {
				P_M_yield_Opt(-1.0*P_Axial);
				P_M_capping_Opt(-1.0*P_Axial);
				P_M_splice_Opt(-1.0*P_Axial);
				Shear_strength_Opt(-1.0*P_Axial);
			}

			//Flexure-Critical (Retrofitted)============================================================================================ 
			dypos = fypos_Rot / Kepos_Rot;
			dyneg = fyneg_Rot / Keneg_Rot;
			if (lb == 0.0)
				dcappos_Rot = fypos_Rot / Kepos_Rot + 0.03 + 0.014*((roT*fyT + roj * fje*(fmin(16 * H / B, 16)) / bj) / (fc*(0.2 - (P_Axial_average / (B*H*fc)))));
			else
				dcappos_Rot = fypos_Rot / Kepos_Rot + 0.03 + 0.014*((roT*fyT + roj * fje*(fmin(16 * H / B, 16)) / bj) / (fc*(0.2 - (P_Axial_average / (B*H*fc)))))*(fmin(lb / ld, 1.0));
			dcapneg_Rot = -dcappos_Rot;

			// Removing Flexure-Shear Critical Check======================================================================================
			dcappos_FlexShr_Rot = 1000.0*dcappos_Rot;
			dcapneg_FlexShr_Rot = -dcappos_FlexShr_Rot;

			// Removing Flexure-splice Check=============================================================================================
			dypos_splice = 1000.0 * fypos_Rot / Kepos_Rot;
			dyneg_splice = 1000.0 * fyneg_Rot / Keneg_Rot;
			dcappos_splice_Rot = 1000.0 * dcappos_Rot;
			dcapneg_splice_Rot = -dcappos_splice_Rot;

			// Removing Shear-Critical Check============================================================================================
			fypos_Shr = 1000000.0 * Vcol;
			dypos_Shr = fypos_Shr / Kepos_Shr;
			dcappos_Shr = dypos_Shr + 100000.0 * max((0.0025 + 1.4*roT)*La, /*0.001*/ 0.1*dypos_Shr); // I had to use different sign for the Axial load here
			fyneg_Shr = -1000000.0*Vcol;
			dyneg_Shr = -fypos_Shr / Kepos_Shr;
			dcapneg_Shr = -dcappos_Shr;

			defineBackbone();
			defineBackbone_Shr();
		}
	}

	if (Retrofit_flag == 2) {  // Steel jacket
		if (flag_Filure_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {

			double roj = (2.0*(B + H)*tj) / (B*H);

			double psi_t = 1.0;
			double psi_e = 1.0;
			double psi_s;
			if (dbL <= 0.8)
				psi_s = 0.8;
			else
				psi_s = 1.0;

			double psi_g;
			if (fyL <= 80.0)                       // For grade 60
				psi_g = 1.0;
			else if (fyL <= 100.0 && fyL > 80.0)   // For grade 80
				psi_g = 1.15;
			else                                   // For grade 100
				psi_g = 1.3;

			double cb = fmin(C_C + dbT + dbL / 2.0, (B - 2.0 * (C_C + dbT + dbL / 2.0) / (nL_EW - 1.0)) / 2.0);

			double AsV = PI * pow(dbT, 2.0) / 4.0;
			double Atr = n_leg * AsV;
			double Ktr = 40.0*(Atr) / (S*nL_EW);

			ld = (3.0 / 40.0) * (fyL * 1000.0 * psi_t * psi_e * psi_s * psi_g / (pow(fc*1000.0, 0.5) * fmin((cb + Ktr) / dbL, 2.5)))*dbL;

			if (PVM_Flag == 1) {
				P_M_yield_Opt(P_Axial_average);
				P_M_capping_Opt(P_Axial_average);
				P_M_splice_Opt(P_Axial_average);
				Shear_strength_Opt(P_Axial_average);
			}

			if (PVM_Flag != 1) {
				P_M_yield_Opt(-1.0*P_Axial);
				P_M_capping_Opt(-1.0*P_Axial);
				P_M_splice_Opt(-1.0*P_Axial);
				Shear_strength_Opt(-1.0*P_Axial);
			}

			//Flexure-Critical (Retrofitted)============================================================================================ 
			dypos = fypos_Rot / Kepos_Rot;
			dyneg = fyneg_Rot / Keneg_Rot;
			if (lb == 0.0)
				dcappos_Rot = fypos_Rot / Kepos_Rot + 0.03 + 0.014*((roT*fyT + roj * fje*(fmin(16 * H / B, 16)) / bj) / (fc*(0.2 - (P_Axial_average / (B*H*fc)))));
			else
				dcappos_Rot = fypos_Rot / Kepos_Rot + 0.03 + 0.014*((roT*fyT + roj * fje*(fmin(16 * H / B, 16)) / bj) / (fc*(0.2 - (P_Axial_average / (B*H*fc)))))*(fmin(lb / ld, 1.0));
			dcapneg_Rot = -dcappos_Rot;

			// Removing Flexure-Shear Critical Check======================================================================================
			dcappos_FlexShr_Rot = 1000.0*dcappos_Rot;
			dcapneg_FlexShr_Rot = -dcappos_FlexShr_Rot;

			// Removing Flexure-splice Check=============================================================================================
			dypos_splice = 1000.0 * fypos_Rot / Kepos_Rot;
			dyneg_splice = 1000.0 * fyneg_Rot / Keneg_Rot;
			dcappos_splice_Rot = 1000.0 * dcappos_Rot;
			dcapneg_splice_Rot = -dcappos_splice_Rot;

			// Removing Shear-Critical Check============================================================================================
			fypos_Shr = 1000000.0 * Vcol;
			dypos_Shr = fypos_Shr / Kepos_Shr;
			dcappos_Shr = dypos_Shr + 100000.0 * max((0.0025 + 1.4*roT)*La, /*0.001*/ 0.1*dypos_Shr); // I had to use different sign for the Axial load here
			fyneg_Shr = -1000000.0*Vcol;
			dyneg_Shr = -fypos_Shr / Kepos_Shr;
			dcapneg_Shr = -dcappos_Shr;

			defineBackbone();
			defineBackbone_Shr();
		}
	}

	if (flag_entering_hardening_Flex_pos_Rot != 1 && flag_entering_hardening_Splice_pos_Rot != 1 && flag_entering_hardening_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {
		dpeakmax = dpeakmax_bench = fmin(dypos, dypos_splice);
		ffmax = slope_pos * dpeakmax + Intcpt_slope_pos;
	}
	if (flag_entering_hardening_Flex_neg_Rot != 1 && flag_entering_hardening_Splice_neg_Rot != 1 && flag_entering_hardening_Shr != 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {
		dpeakmin = dpeakmin_bench = fmax(dyneg, dyneg_splice);
		ffmin = -(slope_neg * fabs(dpeakmin) - Intcpt_slope_neg);
	}

	if (flag_entering_hardening_Shr_pos != 1 && flag_Filure_Shr != 1 && flag_entering_hardening_Flex_Rot != 1 && flag_entering_hardening_Splice_Rot != 1) {
		dpeakmax_Shr = dpeakmax_bench_Shr = dypos_Shr;
		ffmax_Shr = slope_pos_Shr * dpeakmax_Shr + Intcpt_slope_pos_Shr;
	}

	if (flag_entering_hardening_Shr_neg != 1 && flag_Filure_Shr != 1 && flag_entering_hardening_Flex_Rot != 1 && flag_entering_hardening_Splice_Rot != 1) {
		dpeakmin_Shr = dpeakmin_bench_Shr = dyneg_Shr;
		ffmin_Shr = -(slope_neg_Shr * fabs(dpeakmin_Shr) - Intcpt_slope_neg_Shr);
	}
}

void
GMG_CMAC2D::Calibration(void)
{
	if (Retrofit_flag == 0) {

		if (flag_entering_hardening_Flex_Rot == 1 && flag_FlexShr_Filure_Rot != 1) {

			const double PI = 3.141592653589793238463;

			//Transverse reinf. calculation
			double AsV = PI * pow(dbT, 2.0) / 4.0;
			double roT = (n_leg*AsV) / (B*S);
			//Longitudinal reinf. calculation
			double AsL = PI * pow(dbL, 2.0) / 4.0;
			double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);
			// Geometry
			double deff_EW = 0.8*H;

			// Cyclic Properties=============================================
			// (i) ER Calculations

			// Before Capping
			Er_Rot_Post_Yielding = 15.44*roL + 0.23;

			if (Er_Rot_Post_Yielding < 0.2)
				Er_Rot_Post_Yielding = 0.2;
			else if (Er_Rot_Post_Yielding > 0.6)
				Er_Rot_Post_Yielding = 0.6;

			alpha_Er_Rot_Post_Yielding = 0.2;
			beta_Er_Rot_Post_Yielding = -0.2 - 3.0*roL;

			// After Capping
			Er_Rot_Post_Capping = 15.44*roL + 0.23;
			if (Er_Rot_Post_Capping < 0.2)
				Er_Rot_Post_Capping = 0.2;
			else if (Er_Rot_Post_Capping > 0.6)
				Er_Rot_Post_Capping = 0.6;

			alpha_Er_Rot_Post_Capping = 0.2;
			beta_Er_Rot_Post_Capping = -0.2 - 3.0*roL;

			// (ii) Kun :
			// Before Capping
			if (Er_Rot_Post_Yielding <= 0.25)
				Kun_Rot_Post_Yielding = 4.0;
			else
				Kun_Rot_Post_Yielding = 65.7*Er_Rot_Post_Yielding - 12.4;

			// After Capping
			if (Er_Rot_Post_Capping <= 0.25)
				Kun_Rot_Post_Capping = 4.0;
			else
				Kun_Rot_Post_Capping = 65.7*Er_Rot_Post_Capping - 12.4;

			// Before Capping
			Kr_Rot_Post_Yielding = 1.25;
			// After Capping
			Kr_Rot_Post_Capping = 1.2;

			// Damage Properties=============================================
			delta_ratio_max_hard_Rot = 0.0;
			Ref_Energy_Coe_Rot = 20.0;
			solpe_post_yielding_Rot = 0.0;
			solpe_post_capping_Rot = 0.13;
			C1_Rot = 0.2;
			C2_Rot = 0.8;
			C3_Rot = 1.0;
		}

		if (flag_entering_hardening_Flex_Rot == 1 && flag_FlexShr_Filure_Rot == 1) {

			const double PI = 3.141592653589793238463;
			//Transverse reinf. calculation
			double AsV = PI * pow(dbT, 2.0) / 4.0;
			double roT = (n_leg*AsV) / (B*S);

			//Longitudinal reinf. calculation
			double AsL = PI * pow(dbL, 2.0) / 4.0;
			double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);

			// Geometry
			double deff_EW = 0.8*H;


			// Cyclic Properties=============================================
			// (i) ER:

			//Er_Rot_Post_Capping_FlexShr = 0.6;
			Er_Rot_Post_Capping_FlexShr = 13.2*roL;
			if (Er_Rot_Post_Capping_FlexShr < 0.2)
				Er_Rot_Post_Capping_FlexShr = 0.2;
			else if (Er_Rot_Post_Capping_FlexShr > 0.6)
				Er_Rot_Post_Capping_FlexShr = 0.6;

			alpha_Er_Rot_Post_Capping_FlexShr = 0.2;
			beta_Er_Rot_Post_Capping_FlexShr = -0.3 - 3.0*roL;

			// (ii) Kun:

			if (Er_Rot_Post_Capping_FlexShr <= 0.25)
				Kun_Rot_Post_Capping_FlexShr = 4.0;
			else
				Kun_Rot_Post_Capping_FlexShr = 65.7*Er_Rot_Post_Capping - 12.4;

			// (iii) Krl
			Kr_Rot_Post_Capping_FlexShr = 1.0;

			// Damage Properties=============================================
			solpe_post_capping_FlexShr_Rot = 0.51;
			C1_FlexShr_Rot = 0.2;
			C2_FlexShr_Rot = 0.8;
			C3_FlexShr_Rot = 1.0;

		}


		if (flag_entering_hardening_Splice_Rot == 1) {

			const double PI = 3.141592653589793238463;

			//Transverse reinf. calculation
			double AsV = PI * pow(dbT, 2.0) / 4.0;
			double roT = (n_leg*AsV) / (B*S);

			//Longitudinal reinf. calculation
			double AsL = PI * pow(dbL, 2.0) / 4.0;
			double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);

			// Geometry
			double deff_EW = 0.8*H;

			// Cyclic Properties=============================================
			// (i) ER:
			Er_Rot_Post_Yielding = -0.025 + 0.16 * lb / H + 0.026*La / deff_EW;

			if (Er_Rot_Post_Yielding < 0.2)
				Er_Rot_Post_Yielding = 0.2;
			else if (Er_Rot_Post_Yielding > 0.6)
				Er_Rot_Post_Yielding = 0.6;

			alpha_Er_Rot_Post_Yielding = 0.2;
			beta_Er_Rot_Post_Yielding = -0.3 - 3.0*roL;

			Er_Rot_Post_Capping_splice = -0.025 + 0.16 * lb / H + 0.026*La / deff_EW;

			if (Er_Rot_Post_Capping_splice < 0.2)
				Er_Rot_Post_Capping_splice = 0.2;
			else if (Er_Rot_Post_Capping_splice > 0.6)
				Er_Rot_Post_Capping_splice = 0.6;

			alpha_Er_Rot_Post_Capping_splice = 0.2;
			beta_Er_Rot_Post_Capping_splice = -0.3 - 3.0*roL;

			// (ii) Kun and Krel:
			if (Er_Rot_Post_Yielding <= 0.25)
				Kun_Rot_Post_Yielding = 4.0;
			else
				Kun_Rot_Post_Yielding = 65.7*Er_Rot_Post_Yielding - 12.4;

			if (Er_Rot_Post_Capping_splice <= 0.25)
				Kun_Rot_Post_Capping_splice = 4.0;
			else
				Kun_Rot_Post_Capping_splice = 65.7*Er_Rot_Post_Capping_splice - 12.4;

			// (iii) Krel:
			Kr_Rot_Post_Yielding = 1.25;
			Kr_Rot_Post_Capping_splice = 1.3;

			// Damage Properties=============================================
			delta_ratio_max_hard_Rot = 0.0;
			Ref_Energy_Coe_Rot = 20.0;
			solpe_post_yielding_Rot = 0.0;
			solpe_post_capping_splice_Rot = 0.19;
			C1_splice_Rot = 0.2;
			C2_splice_Rot = 0.8;
			C3_splice_Rot = 1.0;
			C1_Rot = 0.2;
			C2_Rot = 0.8;
			C3_Rot = 1.0;

		}


		if (flag_FlexSplice_Filure_Rot == 1) {

			const double PI = 3.141592653589793238463;

			//Transverse reinf. calculation
			double AsV = PI * pow(dbT, 2.0) / 4.0;
			double roT = (n_leg*AsV) / (B*S);

			//Longitudinal reinf. calculation
			double AsL = PI * pow(dbL, 2.0) / 4.0;
			double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);

			// Geometry
			double deff_EW = 0.8*H;

			// Cyclic Properties=============================================
			// (i) ER:
			Er_Rot_Post_Yielding = 15.44 * roL + 0.23;
			if (Er_Rot_Post_Yielding < 0.2)
				Er_Rot_Post_Yielding = 0.2;
			else if (Er_Rot_Post_Yielding > 0.6)
				Er_Rot_Post_Yielding = 0.6;

			alpha_Er_Rot_Post_Yielding = 0.2;
			beta_Er_Rot_Post_Yielding = -0.3 - 3.0*roL;


			Er_Rot_Post_Capping_splice = -0.025 + 0.16 * lb / H + 0.026*La / deff_EW;
			if (Er_Rot_Post_Capping_splice < 0.2)
				Er_Rot_Post_Capping_splice = 0.2;
			else if (Er_Rot_Post_Capping_splice > 0.6)
				Er_Rot_Post_Capping_splice = 0.6;

			alpha_Er_Rot_Post_Capping_splice = 0.2;
			beta_Er_Rot_Post_Capping_splice = -0.3 - 3.0*roL;

			// Kun and Krel :
			if (Er_Rot_Post_Yielding <= 0.25)
				Kun_Rot_Post_Yielding = 4.0;
			else
				Kun_Rot_Post_Yielding = 65.7*Er_Rot_Post_Yielding - 12.4;

			if (Er_Rot_Post_Capping_splice <= 0.25)
				Kun_Rot_Post_Capping_splice = 4.0;
			else
				Kun_Rot_Post_Capping_splice = 65.7*Er_Rot_Post_Capping_splice - 12.4;

			Kr_Rot_Post_Yielding = 1.25;
			Kr_Rot_Post_Capping_splice = 1.3;
			// Damage Properties=============================================
			delta_ratio_max_hard_Rot = 0.0;
			Ref_Energy_Coe_Rot = 20.0;
			solpe_post_yielding_Rot = 0.0;
			solpe_post_capping_splice_Rot = 0.19;
			C1_splice_Rot = 0.2;
			C2_splice_Rot = 0.8;
			C3_splice_Rot = 1.0;
		}
	}

	else if (Retrofit_flag == 1) { // FRP jacket

		const double PI = 3.141592653589793238463;

		//Transverse reinf. calculation
		double AsV = PI * pow(dbT, 2.0) / 4.0;
		double roT = (n_leg*AsV) / (B*S);
		//Longitudinal reinf. calculation
		double AsL = PI * pow(dbL, 2.0) / 4.0;
		double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);
		// Geometry
		double deff_EW = 0.8*H;

		// Cyclic Properties=============================================
		// (i) ER Calculations

		// Before Capping
		Er_Rot_Post_Yielding = 0.4;

		if (Er_Rot_Post_Yielding < 0.2)
			Er_Rot_Post_Yielding = 0.2;
		else if (Er_Rot_Post_Yielding > 0.6)
			Er_Rot_Post_Yielding = 0.6;

		alpha_Er_Rot_Post_Yielding = 0.2;
		beta_Er_Rot_Post_Yielding = -0.3;

		// After Capping
		Er_Rot_Post_Capping = 0.4;
		if (Er_Rot_Post_Capping < 0.2)
			Er_Rot_Post_Capping = 0.2;
		else if (Er_Rot_Post_Capping > 0.6)
			Er_Rot_Post_Capping = 0.6;

		alpha_Er_Rot_Post_Capping = 0.2;
		beta_Er_Rot_Post_Capping = -0.3;

		// (ii) Kun :
		// Before Capping
		if (Er_Rot_Post_Yielding <= 0.25)
			Kun_Rot_Post_Yielding = 4.0;
		else
			Kun_Rot_Post_Yielding = 65.7*Er_Rot_Post_Yielding - 12.4;

		// After Capping
		if (Er_Rot_Post_Capping <= 0.25)
			Kun_Rot_Post_Capping = 4.0;
		else
			Kun_Rot_Post_Capping = 65.7*Er_Rot_Post_Capping - 12.4;

		// (iii) Krel:
		// Before Capping
		Kr_Rot_Post_Yielding = 1.2;
		// After Capping
		Kr_Rot_Post_Capping = 1.1;


		// Damage Properties=============================================
		delta_ratio_max_hard_Rot = 0.0;
		Ref_Energy_Coe_Rot = 20.0;
		solpe_post_yielding_Rot = 0.0;
		solpe_post_capping_Rot = 0.13;
		C1_Rot = 0.2;
		C2_Rot = 0.8;
		C3_Rot = 1.0;

	}

	else { // Steel jacket

		const double PI = 3.141592653589793238463;

		//Transverse reinf. calculation
		double AsV = PI * pow(dbT, 2.0) / 4.0;
		double roT = (n_leg*AsV) / (B*S);
		//Longitudinal reinf. calculation
		double AsL = PI * pow(dbL, 2.0) / 4.0;
		double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);
		// Geometry
		double deff_EW = 0.8*H;

		// Cyclic Properties=============================================
		// (i) ER Calculations

		// Before Capping
		Er_Rot_Post_Yielding = 0.4;

		if (Er_Rot_Post_Yielding < 0.2)
			Er_Rot_Post_Yielding = 0.2;
		else if (Er_Rot_Post_Yielding > 0.6)
			Er_Rot_Post_Yielding = 0.6;

		alpha_Er_Rot_Post_Yielding = 0.2;
		beta_Er_Rot_Post_Yielding = -0.3;

		// After Capping
		Er_Rot_Post_Capping = 15.44*roL + 0.23;
		if (Er_Rot_Post_Capping < 0.2)
			Er_Rot_Post_Capping = 0.2;
		else if (Er_Rot_Post_Capping > 0.6)
			Er_Rot_Post_Capping = 0.6;

		alpha_Er_Rot_Post_Capping = 0.2;
		beta_Er_Rot_Post_Capping = -0.3;

		// (ii) Kun :
		// Before Capping
		if (Er_Rot_Post_Yielding <= 0.25)
			Kun_Rot_Post_Yielding = 4.0;
		else
			Kun_Rot_Post_Yielding = 65.7*Er_Rot_Post_Yielding - 12.4;

		// After Capping
		if (Er_Rot_Post_Capping <= 0.25)
			Kun_Rot_Post_Capping = 4.0;
		else
			Kun_Rot_Post_Capping = 65.7*Er_Rot_Post_Capping - 12.4;

		// (iii) Krel:
		// Before Capping
		Kr_Rot_Post_Yielding = 1.2;
		// After Capping
		Kr_Rot_Post_Capping = 1.1;


		// Damage Properties=============================================
		delta_ratio_max_hard_Rot = 0.0;
		Ref_Energy_Coe_Rot = 20.0;
		solpe_post_yielding_Rot = 0.0;
		solpe_post_capping_Rot = 0.13;
		C1_Rot = 0.2;
		C2_Rot = 0.8;
		C3_Rot = 1.0;

	}
}



void
GMG_CMAC2D::Calibration_Shr(void)
{

	const double PI = 3.141592653589793238463;

	//Transverse reinf. calculation
	double AsV = PI * pow(dbT, 2.0) / 4.0;
	double roT = (n_leg*AsV) / (B*S);

	//Longitudinal reinf. calculation
	double AsL = PI * pow(dbL, 2.0) / 4.0;
	double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);

	// Geometry
	double deff_EW = 0.8*H;

	// Cyclic Properties===========================================
	// (i) ER:
	// Before Capping
	Er_Shr_Post_Yielding = fmax(fmin(0.25 + 0.014*Vcap_for_ER*1000.0 / (B*deff_EW*pow(fc*1000.0, 0.5)), 0.6), 0.2);
	alpha_Er_Shr_Post_Yielding = 0.1;
	beta_Er_Shr_Post_Yielding = -0.3;
	// After Capping
	Er_Shr_Post_Capping = fmax(fmin(0.25 + 0.014*Vcap_for_ER*1000.0 / (B*deff_EW*pow(fc*1000.0, 0.5)), 0.6), 0.2);
	alpha_Er_Shr_Post_Capping = 0.1;
	beta_Er_Shr_Post_Capping = -0.3;

	// (ii) Kun:
	// Before Capping
	if (Er_Shr_Post_Yielding <= 0.25)
		Kun_Shr_Post_Yielding = 4.0;
	else
		Kun_Shr_Post_Yielding = 65.7*Er_Shr_Post_Yielding - 12.4;
	// After Capping
	if (Er_Shr_Post_Capping <= 0.25)
		Kun_Shr_Post_Capping = 4.0;
	else
		Kun_Shr_Post_Capping = 65.7*Er_Shr_Post_Capping - 12.4;

	// (iii) Krel:
	// Before Capping
	Kr_Shr_Post_Yielding = 1.25;
	// After Capping
	Kr_Shr_Post_Capping = 1.0;


	// Damage Properties=============================================
	delta_ratio_max_hard_Shr = 0.0;
	Ref_Energy_Coe_Shr = 20.0;
	solpe_post_yielding_Shr = 0.0;
	solpe_post_capping_Shr = 1.01;
	C1_Shr = 0.2;
	C2_Shr = 0.8;
	C3_Shr = 1.0;

}




void
GMG_CMAC2D::defineBackbone(void)
{
	const double PI = 3.141592653589793238463;

	//Transverse reinf. calculation
	double n = 10.0;
	double AsV = PI * pow(dbT, 2.0) / 4.0;
	double roT = (n_leg*AsV) / (B*S);

	//Longitudinal reinf. calculation
	double AsL = PI * pow(dbL, 2.0) / 4.0;
	double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);

	// Geometry
	double deff_EW = 0.8*H;

	BenMark = d_Rot;
	if (flag_entering_hardening_Flex_Rot == 1 && flag_entering_hardening_Splice_Rot != 1/* && flag_FlexShr_Filure_Rot != 1 */) {
		if (BenMark >= 0) {
			R_dypos = fypos_Rot / Kepos_Rot;//ok ok
			R_fypos = fypos_Rot;//ok ok
			R_dcappos = dcappos_Rot;//ok ok
			R_fcappos = fcappos_Rot;//ok ok

			// Converting global equation to loacal for k_deg
			double Lamda_deg_Spr_Flx = (fmin(-0.105 + 0.15 * (roT*fyT / fc) / (0.2 - (8.0 / 7.0)*P_Axial_average / (B*H*fc)) + 0.06 / pow(S / fmin(0.125*B, 0.125*H), 0.5), -0.005)) / (1.0 + n * (1.0 - (fmin(-0.105 + 0.15 * (roT*fyT / fc) / (0.2 - (8.0 / 7.0)*P_Axial_average / (B*H*fc)) + 0.06 / pow(S / fmin(0.125*B, 0.125*H), 0.5), -0.005))));
			R_Kdegpos = fmax(Lamda_deg_Spr_Flx, -0.015) * Kepos_Rot;


			R_fyneg = fyneg_Rot; //ok ok
			R_dyneg = fyneg_Rot / Keneg_Rot; //ok ok
			R_fcapneg = fcapneg_Rot;//ok ok
			R_dcapneg = dcapneg_Rot;
			R_Kdegneg = R_Kdegpos;
			dppos = R_dcappos - R_dypos; //ok ok
			slope_pos = (R_fcappos - fypos_Rot) / dppos; //ok ok

			if (d_Axial > 0)
				R_frespos = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos;//ok ok
			else
				R_frespos = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos;//ok ok

			R_drespos = R_dcappos + (R_frespos - R_fcappos) / R_Kdegpos; //ok ok
			R_delUpos = 3.0 * R_drespos; //ok ok
			Intcpt_slope_pos = R_fcappos - (slope_pos * R_dcappos); //ok ok
			Intcpt_deg_pos = fabs(R_fcappos - R_Kdegpos * R_dcappos); //ok ok
			Intcpt_res_pos = fabs(R_frespos - R_Kdegpos * R_delUpos); //ok ok
			Intcpt_Xaxis_pos = (R_Kdegpos * R_delUpos - R_frespos) / R_Kdegpos; //ok ok

			dpneg = R_dcapneg - R_dyneg; //ok ok
			slope_neg = (R_fcapneg - R_fyneg) / dpneg; //ok ok

			if (d_Axial > 0)
				R_fresneg = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg; //ok ok
			else
				R_fresneg = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg; //ok ok
			R_dresneg = R_dcapneg + (R_fresneg - R_fcapneg) / R_Kdegneg; //ok ok
			R_delUneg = 3.0 * R_dresneg; //ok ok
			Intcpt_slope_neg = R_fcapneg - (slope_neg * R_dcapneg); //ok ok
			Intcpt_deg_neg = fabs(R_fcapneg - R_Kdegneg * R_dcapneg); //Ok ok
			Intcpt_res_neg = fabs(R_fresneg - R_Kdegneg * R_delUneg); //ok ok
			Intcpt_Xaxis_neg = (R_Kdegneg * R_delUneg - R_fresneg) / R_Kdegneg; //ok ok

		}

		else if (BenMark < 0) {
			R_dyneg = fyneg_Rot / Keneg_Rot;
			R_fyneg = fyneg_Rot;
			R_dcapneg = dcapneg_Rot;
			R_fcapneg = fcapneg_Rot;
			// Converting global equation to loacal for k_deg
			double Lamda_deg_Spr_Flx = (fmin(-0.105 + 0.15 * (roT*fyT / fc) / (0.2 - (8.0 / 7.0)*P_Axial_average / (B*H*fc)) + 0.06 / pow(S / fmin(0.125*B, 0.125*H), 0.5), -0.005)) / (1.0 + n * (1.0 - (fmin(-0.105 + 0.15 * (roT*fyT / fc) / (0.2 - (8.0 / 7.0)*P_Axial_average / (B*H*fc)) + 0.06 / pow(S / fmin(0.125*B, 0.125*H), 0.5), -0.005))));
			R_Kdegneg = fmax(Lamda_deg_Spr_Flx, -0.015) * Keneg_Rot;

			R_fypos = fypos_Rot;
			R_dypos = fypos_Rot / Kepos_Rot;
			R_dcappos = dcappos_Rot;
			R_fcappos = fcappos_Rot;
			R_Kdegpos = R_Kdegneg;
			dpneg = R_dcapneg - R_dyneg;
			slope_neg = (R_fcapneg - R_fyneg) / dpneg;

			if (d_Axial > 0)
				R_fresneg = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg; //ok ok
			else
				R_fresneg = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg; //ok ok

			R_dresneg = R_dcapneg + (R_fresneg - R_fcapneg) / R_Kdegneg;
			R_delUneg = 3.0 * R_dresneg;
			Intcpt_slope_neg = R_fcapneg - (slope_neg * R_dcapneg);
			Intcpt_deg_neg = fabs(R_fcapneg - R_Kdegneg * R_dcapneg);
			Intcpt_res_neg = fabs(R_fresneg - R_Kdegneg * R_delUneg);
			Intcpt_Xaxis_neg = (R_Kdegneg * R_delUneg - R_fresneg) / R_Kdegneg;

			dppos = R_dcappos - R_dypos;
			slope_pos = (R_fcappos - R_fypos) / dppos;

			if (d_Axial > 0)
				R_frespos = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos;//ok ok
			else
				R_frespos = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos;//ok ok

			R_drespos = R_dcappos + (R_frespos - R_fcappos) / R_Kdegpos;
			R_delUpos = 3.0 * R_drespos;
			Intcpt_slope_pos = R_fcappos - (slope_pos * R_dcappos);
			Intcpt_deg_pos = fabs(R_fcappos - R_Kdegpos * R_dcappos);
			Intcpt_res_pos = fabs(R_frespos - R_Kdegpos * R_delUpos);
			Intcpt_Xaxis_pos = (R_Kdegpos * R_delUpos - R_frespos) / R_Kdegpos;
		}

	}


	if (flag_entering_hardening_Flex_Rot == 1 && flag_FlexShr_Filure_Rot == 1) {
		if (BenMark >= 0) {

			R_dcappos = fmin(dcappos_FlexShr_Rot, d_Rot);
			R_fcappos = slope_pos * R_dcappos + Intcpt_slope_pos;
			// Converting global equation to loacal for k_deg
			double Lamda_deg_Spr_Flx_Shr = fmin((-0.003 * pow(roT*fyT / (roL * fyL), -1.5)), -0.005) / (1.0 + n * (1.0 - fmin((-0.003 * pow(roT*fyT / (roL * fyL), -1.5)), -0.005)));
			R_Kdegpos = fmax(Lamda_deg_Spr_Flx_Shr, -0.015) * Kepos_Rot;

			ratio = (R_dcappos - R_dypos) / (dcappos_Rot - R_dypos);
			R_dcapneg = fmax(((dcapneg_Rot - R_dyneg) * ratio + R_dyneg), dcapneg_FlexShr_Rot);
			R_fcapneg = -(slope_neg * fabs(R_dcapneg) - Intcpt_slope_neg);

			//R_Kdegneg = -0.01*Keneg_Rot;
			R_Kdegneg = R_Kdegpos;

			if (d_Axial > 0)
				R_frespos = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos;//ok ok
			else
				R_frespos = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos;//ok ok

			R_drespos = R_dcappos + (R_frespos - R_fcappos) / R_Kdegpos;
			R_delUpos = 3.0 * R_drespos;

			Intcpt_deg_pos = fabs(R_fcappos - R_Kdegpos * R_dcappos);
			Intcpt_res_pos = fabs(R_frespos - R_Kdegpos * R_delUpos);
			Intcpt_Xaxis_pos = (R_Kdegpos * R_delUpos - R_frespos) / R_Kdegpos;

			if (d_Axial > 0)
				R_fresneg = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg; //ok ok
			else
				R_fresneg = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg; //ok ok

			R_dresneg = R_dcapneg + (R_fresneg - R_fcapneg) / R_Kdegneg;
			R_delUneg = 3.0 * R_dresneg;

			Intcpt_deg_neg = fabs(R_fcapneg - R_Kdegneg * R_dcapneg);
			Intcpt_res_neg = fabs(R_fresneg - R_Kdegneg * R_delUneg);
			Intcpt_Xaxis_neg = (R_Kdegneg * R_delUneg - R_fresneg) / R_Kdegneg;
		}

		else if (BenMark < 0) {

			R_dcapneg = fmax(dcapneg_FlexShr_Rot, d_Rot);
			R_fcapneg = -(slope_neg * fabs(R_dcapneg) - Intcpt_slope_neg);

			// Converting global equation to loacal for k_deg

			double Lamda_deg_Spr_Flx_Shr = fmin((-0.003 * pow(roT*fyT / (roL * fyL), -1.5)), -0.005) / (1.0 + n * (1.0 - fmin((-0.003 * pow(roT*fyT / (roL * fyL), -1.5)), -0.005)));
			R_Kdegneg = fmax(Lamda_deg_Spr_Flx_Shr, -0.015) * Keneg_Rot;

			ratio = (R_dcapneg - R_dyneg) / (dcapneg_Rot - R_dyneg);
			R_dcappos = fmin(((dcappos_Rot - R_dypos) * ratio + R_dypos), dcappos_FlexShr_Rot);
			R_fcappos = slope_pos * R_dcappos + Intcpt_slope_pos;
			//R_Kdegpos = -0.01*Kepos_Rot;
			R_Kdegpos = R_Kdegneg;

			if (d_Axial > 0)
				R_fresneg = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg; //ok ok
			else
				R_fresneg = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg; //ok ok

			R_dresneg = R_dcapneg + (R_fresneg - R_fcapneg) / R_Kdegneg;
			R_delUneg = 3.0 * R_dresneg;

			Intcpt_deg_neg = fabs(R_fcapneg - R_Kdegneg * R_dcapneg);
			Intcpt_res_neg = fabs(R_fresneg - R_Kdegneg * R_delUneg);
			Intcpt_Xaxis_neg = (R_Kdegneg * R_delUneg - R_fresneg) / R_Kdegneg;

			if (d_Axial > 0)
				R_frespos = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos;//ok ok
			else
				R_frespos = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos;//ok ok

			R_drespos = R_dcappos + (R_frespos - R_fcappos) / R_Kdegpos;
			R_delUpos = 3.0 * R_drespos;

			Intcpt_deg_pos = fabs(R_fcappos - R_Kdegpos * R_dcappos);
			Intcpt_res_pos = fabs(R_frespos - R_Kdegpos * R_delUpos);
			Intcpt_Xaxis_pos = (R_Kdegpos * R_delUpos - R_frespos) / R_Kdegpos;
		}
	}


	if (flag_entering_hardening_Splice_Rot == 1) {
		if (BenMark >= 0) {
			R_fypos = fypos_splice_Rot;
			R_dypos = R_fypos / Kepos_Rot;
			R_dcappos = dcappos_splice_Rot;
			R_fcappos = fcappos_splice_Rot;

			double Lamda_deg_Spr_Spl = (fmin(-0.13 + 0.5 * roT*fyT / (roL * fyL), -0.005)) / (1.0 + n * (1.0 - (fmin(-0.13 + 0.5 * roT*fyT / (roL * fyL), -0.005))));
			R_Kdegpos = fmax(Lamda_deg_Spr_Spl, -0.015) * Kepos_Rot;

			R_fyneg = fyneg_splice_Rot;
			R_dyneg = R_fyneg / Keneg_Rot;
			R_fcapneg = fcapneg_splice_Rot;
			R_dcapneg = dcapneg_splice_Rot;
			R_Kdegneg = R_Kdegpos;

			dppos = R_dcappos - R_dypos;
			slope_pos = (R_fcappos - R_fypos) / dppos;

			R_frespos = fmax(fmin(0.15 + 36.0 *fmax(roT, 0.05), 0.4), 0.15) * R_fypos;

			R_drespos = R_dcappos + (R_frespos - R_fcappos) / R_Kdegpos;
			R_delUpos = 3.0 * R_drespos;
			Intcpt_slope_pos = R_fcappos - (slope_pos * R_dcappos);
			Intcpt_deg_pos = fabs(R_fcappos - R_Kdegpos * R_dcappos);
			Intcpt_res_pos = fabs(R_frespos - R_Kdegpos * R_delUpos);
			Intcpt_Xaxis_pos = (R_Kdegpos * R_delUpos - R_frespos) / R_Kdegpos;

			dpneg = R_dcapneg - R_dyneg;
			slope_neg = (R_fcapneg - R_fyneg) / dpneg;

			R_fresneg = fmax(fmin(0.15 + 36.0 *fmax(roT, 0.05), 0.4), 0.15) * R_fyneg;

			R_dresneg = R_dcapneg + (R_fresneg - R_fcapneg) / R_Kdegneg;
			R_delUneg = 3.0 * R_dresneg;
			Intcpt_slope_neg = R_fcapneg - (slope_neg * R_dcapneg);
			Intcpt_deg_neg = fabs(R_fcapneg - R_Kdegneg * R_dcapneg);
			Intcpt_res_neg = fabs(R_fresneg - R_Kdegneg * R_delUneg);
			Intcpt_Xaxis_neg = (R_Kdegneg * R_delUneg - R_fresneg) / R_Kdegneg;
		}

		else if (BenMark < 0) {
			R_fyneg = fyneg_splice_Rot;
			R_dyneg = R_fyneg / Keneg_Rot;
			R_fcapneg = fcapneg_splice_Rot;
			R_dcapneg = dcapneg_splice_Rot;

			double Lamda_deg_Spr_Spl = (fmin(-0.13 + 0.5 * roT*fyT / (roL * fyL), -0.005)) / (1.0 + n * (1.0 - (fmin(-0.13 + 0.5 * roT*fyT / (roL * fyL), -0.005))));
			R_Kdegneg = fmax(Lamda_deg_Spr_Spl, -0.015) * Keneg_Rot;

			R_fypos = fypos_splice_Rot;
			R_dypos = R_fypos / Kepos_Rot;
			R_dcappos = dcappos_splice_Rot;
			R_fcappos = fcappos_splice_Rot;
			R_Kdegpos = R_Kdegneg;

			dpneg = R_dcapneg - R_dyneg;
			slope_neg = (R_fcapneg - R_fyneg) / dpneg;

			R_fresneg = fmax(fmin(0.15 + 36.0 *fmax(roT, 0.05), 0.4), 0.15) * R_fyneg;

			R_dresneg = R_dcapneg + (R_fresneg - R_fcapneg) / R_Kdegneg;
			R_delUneg = 3.0 * R_dresneg;
			Intcpt_slope_neg = R_fcapneg - (slope_neg * R_dcapneg);
			Intcpt_deg_neg = fabs(R_fcapneg - R_Kdegneg * R_dcapneg);
			Intcpt_res_neg = fabs(R_fresneg - R_Kdegneg * R_delUneg);
			Intcpt_Xaxis_neg = (R_Kdegneg * R_delUneg - R_fresneg) / R_Kdegneg;

			dppos = R_dcappos - R_dypos;
			slope_pos = (R_fcappos - R_fypos) / dppos;

			R_frespos = fmax(fmin(0.15 + 36.0 *fmax(roT, 0.05), 0.4), 0.15) * R_fypos;

			R_drespos = R_dcappos + (R_frespos - R_fcappos) / R_Kdegpos;
			R_delUpos = 3.0 * R_drespos;
			Intcpt_slope_pos = R_fcappos - (slope_pos * R_dcappos);
			Intcpt_deg_pos = fabs(R_fcappos - R_Kdegpos * R_dcappos);
			Intcpt_res_pos = fabs(R_frespos - R_Kdegpos * R_delUpos);
			Intcpt_Xaxis_pos = (R_Kdegpos * R_delUpos - R_frespos) / R_Kdegpos;

		}
	}



	if (flag_FlexSplice_Filure_Rot == 1) {
		if (BenMark >= 0) {

			R_dcappos = dcappos_splice_Rot;
			R_fcappos = slope_pos * R_dcappos + Intcpt_slope_pos;

			double Lamda_deg_Spr_Spl = (fmin(-0.13 + 0.5 * roT*fyT / (roT * fyL), -0.005)) / (1.0 + n * (1.0 - (fmin(-0.13 + 0.5 * roT*fyT / (roT * fyL), -0.005))));
			R_Kdegpos = fmax(Lamda_deg_Spr_Spl, -0.015) * Kepos_Rot;

			ratio = (R_dcappos - R_dypos) / (dcappos_Rot - R_dypos);
			R_dcapneg = fmax(((dcapneg_Rot - R_dyneg) * ratio + R_dyneg), dcapneg_splice_Rot);
			R_fcapneg = -(slope_neg * fabs(R_dcapneg) - Intcpt_slope_neg);
			R_Kdegneg = R_Kdegpos;

			R_frespos = fmax(fmin(0.15 + 36.0 *fmax(roT, 0.05), 0.4), 0.15) * R_fypos;

			R_drespos = R_dcappos + (R_frespos - R_fcappos) / R_Kdegpos;
			R_delUpos = 3.0 * R_drespos;

			Intcpt_deg_pos = fabs(R_fcappos - R_Kdegpos * R_dcappos);
			Intcpt_res_pos = fabs(R_frespos - R_Kdegpos * R_delUpos);
			Intcpt_Xaxis_pos = (R_Kdegpos * R_delUpos - R_frespos) / R_Kdegpos;

			R_fresneg = fmax(fmin(0.15 + 36.0 *fmax(roT, 0.05), 0.4), 0.15) * R_fyneg;

			R_dresneg = R_dcapneg + (R_fresneg - R_fcapneg) / R_Kdegneg;
			R_delUneg = 3.0 * R_dresneg;

			Intcpt_deg_neg = fabs(R_fcapneg - R_Kdegneg * R_dcapneg);
			Intcpt_res_neg = fabs(R_fresneg - R_Kdegneg * R_delUneg);
			Intcpt_Xaxis_neg = (R_Kdegneg * R_delUneg - R_fresneg) / R_Kdegneg;
		}

		else if (BenMark < 0) {

			R_dcapneg = dcapneg_splice_Rot;
			R_fcapneg = -(slope_neg * fabs(R_dcapneg) - Intcpt_slope_neg);

			// Converting global equation to loacal for k_deg

			double Lamda_deg_Spr_Spl = (fmin(-0.13 + 0.5 * roT*fyT / (roT * fyL), -0.005)) / (1.0 + n * (1.0 - (fmin(-0.13 + 0.5 * roT*fyT / (roT * fyL), -0.005))));
			R_Kdegneg = fmax(Lamda_deg_Spr_Spl, -0.015) * Keneg_Rot;

			ratio = (R_dcapneg - R_dyneg) / (dcapneg_Rot - R_dyneg);
			R_dcappos = fmin(((dcappos_Rot - R_dypos) * ratio + R_dypos), dcappos_splice_Rot);
			R_fcappos = slope_pos * R_dcappos + Intcpt_slope_pos;
			R_Kdegpos = R_Kdegneg;

			R_fresneg = fmax(fmin(0.15 + 36.0 *fmax(roT, 0.05), 0.4), 0.15) * R_fyneg;

			R_dresneg = R_dcapneg + (R_fresneg - R_fcapneg) / R_Kdegneg;
			R_delUneg = 3.0 * R_dresneg;

			Intcpt_deg_neg = fabs(R_fcapneg - R_Kdegneg * R_dcapneg);
			Intcpt_res_neg = fabs(R_fresneg - R_Kdegneg * R_delUneg);
			Intcpt_Xaxis_neg = (R_Kdegneg * R_delUneg - R_fresneg) / R_Kdegneg;

			R_frespos = fmax(fmin(0.15 + 36.0 *fmax(roT, 0.05), 0.4), 0.15) * R_fypos;

			R_drespos = R_dcappos + (R_frespos - R_fcappos) / R_Kdegpos;
			R_delUpos = 3.0 * R_drespos;

			Intcpt_deg_pos = fabs(R_fcappos - R_Kdegpos * R_dcappos);
			Intcpt_res_pos = fabs(R_frespos - R_Kdegpos * R_delUpos);
			Intcpt_Xaxis_pos = (R_Kdegpos * R_delUpos - R_frespos) / R_Kdegpos;

		}
	}



	/* ********************************************************************************************** **
	************************************************************************************************* **
	************************************************************************************************* **
	**                                                                                                **
	**                               Dealing with the referece erergy                                 **
	**                                   Update Damage Function                                       **
	**                                                                                                **
	************************************************************************************************* **
	************************************************************************************************* **
	** ********************************************************************************************** */
	Ed1pos = Ed0pos + fabs(0.5 * R_fypos * R_dypos) + (0.5 * fabs(R_fypos + R_fcappos) * fabs(R_dcappos - R_dypos)) + (0.5 * fabs(R_fcappos + R_frespos) * fabs(R_drespos - R_dcappos)) + 0.5 * fabs(R_frespos) * fabs(R_frespos / R_Kdegpos);
	Ed1neg = Ed0neg + fabs(0.5 * R_fyneg * R_dyneg) + (0.5 * fabs(R_fyneg + R_fcapneg) * fabs(R_dcapneg - R_dyneg)) + (0.5 * fabs(R_fcapneg + R_fresneg) * fabs(R_dresneg - R_dcapneg)) + 0.5 * fabs(R_fresneg) * fabs(R_fresneg / R_Kdegneg);
}


/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Shear Material                                         **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */

void
GMG_CMAC2D::defineBackbone_Shr(void)
{

	const double PI = 3.141592653589793238463;
	double n = 10.0;
	//Transverse reinf. calculation
	double AsV = PI * pow(dbT, 2.0) / 4.0;
	double roT = (n_leg*AsV) / (B*S);

	//Longitudinal reinf. calculation
	double AsL = PI * pow(dbL, 2.0) / 4.0;
	double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);

	// Geometry
	double deff_EW = 0.8*H;

	BenMark_Shr = d_Shr;
	if (BenMark_Shr >= 0) {
		R_dypos_Shr = fypos_Shr / Kepos_Shr;
		R_fypos_Shr = fypos_Shr;
		R_dcappos_Shr = dcappos_Shr;
		fcappos_Shr = 1.05 * R_fypos_Shr;
		R_fcappos_Shr = fcappos_Shr;

		// Converting global equation to loacal for k_deg

		double Lamda_deg_Spr_Shr = (fmin(0.024 - 0.033 * (fypos_Rot / La / Vcol), -0.005)) / (1.0 + n * (1.0 - (fmin(0.024 - 0.033 * (fypos_Rot / La / Vcol), -0.005))));
		Kdegpos_Shr = fmax(Lamda_deg_Spr_Shr, -0.015) * Kepos_Shr;
		R_Kdegpos_Shr = Kdegpos_Shr;

		R_fyneg_Shr = fyneg_Shr;
		R_dyneg_Shr = fyneg_Shr / Keneg_Shr;
		fcapneg_Shr = 1.05 * R_fyneg_Shr;
		R_fcapneg_Shr = fcapneg_Shr;

		R_dcapneg_Shr = dcapneg_Shr;

		Kdegneg_Shr = Kdegpos_Shr;
		R_Kdegneg_Shr = R_Kdegpos_Shr;

		dppos_Shr = R_dcappos_Shr - R_dypos_Shr;
		slope_pos_Shr = (R_fcappos_Shr - R_fypos_Shr) / dppos_Shr;

		if (d_Axial > 0)
			frespos_Shr = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos_Shr;
		else
			frespos_Shr = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos_Shr;
		R_frespos_Shr = frespos_Shr;

		drespos_Shr = R_dcappos_Shr + (R_frespos_Shr - R_fcappos_Shr) / R_Kdegpos_Shr;
		R_drespos_Shr = drespos_Shr;

		delUpos_Shr = 3.0*R_drespos_Shr;
		R_delUpos_Shr = delUpos_Shr;

		Intcpt_slope_pos_Shr = R_fcappos_Shr - (slope_pos_Shr * R_dcappos_Shr);
		Intcpt_deg_pos_Shr = fabs(R_fcappos_Shr - R_Kdegpos_Shr * R_dcappos_Shr);
		Intcpt_res_pos_Shr = fabs(R_frespos_Shr - R_Kdegpos_Shr * R_delUpos_Shr);
		Intcpt_Xaxis_pos_Shr = (R_Kdegpos_Shr * R_delUpos_Shr - R_frespos_Shr) / R_Kdegpos_Shr;

		dpneg_Shr = R_dcapneg_Shr - R_dyneg_Shr;
		slope_neg_Shr = (R_fcapneg_Shr - R_fyneg_Shr) / dpneg_Shr;

		if (d_Axial > 0)
			fresneg_Shr = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg_Shr;
		else
			fresneg_Shr = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg_Shr;

		R_fresneg_Shr = fresneg_Shr;

		dresneg_Shr = R_dcapneg_Shr + (R_fresneg_Shr - R_fcapneg_Shr) / R_Kdegneg_Shr;
		R_dresneg_Shr = dresneg_Shr;

		delUneg_Shr = 3.0*R_dresneg_Shr;
		R_delUneg_Shr = delUneg_Shr;

		Intcpt_slope_neg_Shr = R_fcapneg_Shr - (slope_neg_Shr * R_dcapneg_Shr);
		Intcpt_deg_neg_Shr = fabs(R_fcapneg_Shr - R_Kdegneg_Shr * R_dcapneg_Shr);
		Intcpt_res_neg_Shr = fabs(R_fresneg_Shr - R_Kdegneg_Shr * R_delUneg_Shr);
		Intcpt_Xaxis_neg_Shr = (R_Kdegneg_Shr * R_delUneg_Shr - R_fresneg_Shr) / R_Kdegneg_Shr;

	}

	else if (BenMark_Shr < 0) {
		R_dyneg_Shr = fyneg_Shr / Keneg_Shr;
		R_fyneg_Shr = fyneg_Shr;
		R_dcapneg_Shr = dcapneg_Shr;
		fcapneg_Shr = 1.05 * R_fyneg_Shr;
		R_fcapneg_Shr = fcapneg_Shr;

		// Converting global equation to loacal for k_deg

		double Lamda_deg_Spr_Shr = (fmin(0.024 - 0.033 * (fypos_Rot / La / Vcol), -0.005)) / (1.0 + n * (1.0 - (fmin(0.024 - 0.033 * (fypos_Rot / La / Vcol), -0.005))));
		Kdegneg_Shr = fmax(Lamda_deg_Spr_Shr, -0.015) * Keneg_Shr;
		R_Kdegneg_Shr = Kdegneg_Shr;

		R_fypos_Shr = fypos_Shr;
		R_dypos_Shr = R_fypos_Shr / Kepos_Shr;
		R_dcappos_Shr = dcappos_Shr;
		fcappos_Shr = 1.05 * R_fypos_Shr;
		R_fcappos_Shr = fcappos_Shr;

		Kdegpos_Shr = Kdegneg_Shr;
		R_Kdegpos_Shr = R_Kdegneg_Shr;

		dpneg_Shr = R_dcapneg_Shr - R_dyneg_Shr;
		slope_neg_Shr = (R_fcapneg_Shr - R_fyneg_Shr) / dpneg_Shr;

		if (d_Axial > 0)
			fresneg_Shr = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg_Shr;
		else
			fresneg_Shr = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fyneg_Shr;

		R_fresneg_Shr = fresneg_Shr;

		dresneg_Shr = R_dcapneg_Shr + (R_fresneg_Shr - R_fcapneg_Shr) / R_Kdegneg_Shr;
		R_dresneg_Shr = dresneg_Shr;

		delUneg_Shr = 3.0*R_dresneg_Shr;
		R_delUneg_Shr = delUneg_Shr;

		Intcpt_slope_neg_Shr = R_fcapneg_Shr - (slope_neg_Shr * R_dcapneg_Shr);
		Intcpt_deg_neg_Shr = fabs(R_fcapneg_Shr - R_Kdegneg_Shr * R_dcapneg_Shr);
		Intcpt_res_neg_Shr = fabs(R_fresneg_Shr - R_Kdegneg_Shr * R_delUneg_Shr);
		Intcpt_Xaxis_neg_Shr = (R_Kdegneg_Shr * R_delUneg_Shr - R_fresneg_Shr) / R_Kdegneg_Shr;

		dppos_Shr = R_dcappos_Shr - R_dypos_Shr;
		slope_pos_Shr = (R_fcappos_Shr - R_fypos_Shr) / dppos_Shr;

		if (d_Axial > 0)
			frespos_Shr = fmax(0.24 + 0.4*fmin((Kepos_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos_Shr;
		else
			frespos_Shr = fmax(0.24 + 0.4*fmin((Keneg_Axil * d_Axial) / (B*H*fc), -0.1), 0.15) * R_fypos_Shr;

		R_frespos_Shr = frespos_Shr;

		drespos_Shr = R_dcappos_Shr + (frespos_Shr - fcappos_Shr) / R_Kdegpos_Shr;
		R_drespos_Shr = drespos_Shr;

		delUpos_Shr = 3.0*R_drespos_Shr;
		R_delUpos_Shr = delUpos_Shr;

		Intcpt_slope_pos_Shr = R_fcappos_Shr - (slope_pos_Shr * R_dcappos_Shr);
		Intcpt_deg_pos_Shr = fabs(R_fcappos_Shr - R_Kdegpos_Shr * R_dcappos_Shr);
		Intcpt_res_pos_Shr = fabs(R_frespos_Shr - R_Kdegpos_Shr * R_delUpos_Shr);
		Intcpt_Xaxis_pos_Shr = (R_Kdegpos_Shr * R_delUpos_Shr - R_frespos_Shr) / R_Kdegpos_Shr;

	}

	Ed1pos_Shr = Ed0pos_Shr + fabs(0.5 * R_fypos_Shr * R_dypos_Shr) + (0.5 * fabs(R_fypos_Shr + R_fcappos_Shr) * fabs(R_dcappos_Shr - R_dypos_Shr)) + (0.5 * fabs(R_fcappos_Shr + R_frespos_Shr) * fabs(R_drespos_Shr - R_dcappos_Shr)) + 0.5 * fabs(R_frespos_Shr) * fabs(R_frespos_Shr / R_Kdegpos_Shr);
	Ed1neg_Shr = Ed0neg_Shr + fabs(0.5 * R_fyneg_Shr * R_dyneg_Shr) + (0.5 * fabs(R_fyneg_Shr + R_fcapneg_Shr) * fabs(dpneg_Shr)) + (0.5 * fabs(R_fcapneg_Shr + R_fresneg_Shr) * fabs(R_dresneg_Shr - R_dcapneg_Shr)) + 0.5 * fabs(R_fresneg_Shr) * fabs(R_fresneg_Shr / R_Kdegneg_Shr);
}


/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Rotational Material                                    **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */
int
GMG_CMAC2D::getStateFlag(void)
{

	if ((CstateFlag == 0 || CstateFlag == 12 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu_Rot > 0.0 && d_Rot >= dpeakmax && d_Rot < R_dcappos && BenMark > 0)
		return 12;
	else if ((CstateFlag == 0 || CstateFlag == 12 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu_Rot > 0.0 && d_Rot >= dpeakmax && d_Rot < R_dcappos && BenMark < 0)
		return 12;
	else if ((CstateFlag == 12 || CstateFlag == 2 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu_Rot > 0.0 && d_Rot >= R_dcappos && d_Rot >= dpeakmax && d_Rot < R_drespos && BenMark > 0)
		return 2;
	else if ((CstateFlag == 12 || CstateFlag == 2 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu_Rot > 0.0 && d_Rot >= R_dcappos && d_Rot >= dpeakmax && d_Rot < R_drespos && BenMark < 0)
		return 2;
	else if ((CstateFlag == 12 || CstateFlag == 2 || CstateFlag == 3 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu_Rot > 0.0 && d_Rot >= R_drespos && d_Rot < R_delUpos && d_Rot >= dpeakmax && BenMark > 0)
		return 3;
	else if ((CstateFlag == 12 || CstateFlag == 2 || CstateFlag == 3 || CstateFlag == -5 || CstateFlag == 6 || CstateFlag == -7 || CstateFlag == 5 || CstateFlag == 4) && Tdu_Rot > 0.0 && d_Rot >= R_drespos && d_Rot < R_delUpos && d_Rot >= dpeakmax && BenMark < 0)
		return 3;
	else if ((CstateFlag == 3 || CstateFlag == 30 || CstateFlag == 31) && Tdu_Rot > 0.0 && d_Rot >= R_delUpos && d_Rot <= Intcpt_Xaxis_pos && d_Rot >= dpeakmax && BenMark > 0)
		return 30;
	else if ((CstateFlag == 3 || CstateFlag == 30 || CstateFlag == 31) && Tdu_Rot > 0.0 && d_Rot >= R_delUpos && d_Rot <= Intcpt_Xaxis_pos && d_Rot >= dpeakmax && BenMark < 0)
		return 30;
	else if ((CstateFlag == 30 || CstateFlag == 40) && Tdu_Rot > 0.0 && d_Rot >= Intcpt_Xaxis_pos && d_Rot >= dpeakmax && BenMark > 0)
		return 40;
	else if ((CstateFlag == 30 || CstateFlag == 40) && Tdu_Rot > 0.0 && d_Rot >= Intcpt_Xaxis_pos && d_Rot >= dpeakmax && BenMark < 0)
		return 40;
	else if ((CstateFlag == 30 || CstateFlag == 31) && (Tdu_Rot > 0.0 || Tdu_Rot <= 0.0) && (d_Rot > dpeakmin || d_Rot < dpeakmax))
		return 31;
	else if ((CstateFlag == 40 || CstateFlag == 41) && (Tdu_Rot > 0.0 || Tdu_Rot <= 0.0))
		return 41;

	else if ((CstateFlag == 12 || CstateFlag == 2 || CstateFlag == 3) && Tdu_Rot < 0.0)
		return 4;
	else if ((CstateFlag == 4 || CstateFlag == 5) && Tdu_Rot < 0.0 && d_Rot > dpeakmin)
		return 5;

	else if ((CstateFlag == 4 || CstateFlag == 5 || CstateFlag == 6) && Tdu_Rot > 0.0 && d_Rot < dpeakmax)
		return 6;
	else if ((CstateFlag == 6 || CstateFlag == 7 || CstateFlag == -7) && Tdu_Rot < 0.0 && d_Rot > dpeakmin)
		return 7;

	else if ((CstateFlag == 0 || CstateFlag == -12 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu_Rot < 0.0 && d_Rot <= dpeakmin && d_Rot > R_dcapneg && BenMark > 0)
		return -12;
	else if ((CstateFlag == 0 || CstateFlag == -12 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu_Rot < 0.0 && d_Rot <= dpeakmin && d_Rot > R_dcapneg && BenMark < 0)
		return -12;
	else if ((CstateFlag == -12 || CstateFlag == -2 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu_Rot < 0.0 && d_Rot <= R_dcapneg && d_Rot <= dpeakmin && d_Rot > R_dresneg && BenMark > 0)
		return -2;
	else if ((CstateFlag == -12 || CstateFlag == -2 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu_Rot < 0.0 && d_Rot <= R_dcapneg && d_Rot <= dpeakmin && d_Rot > R_dresneg && BenMark < 0)
		return -2;
	else if ((CstateFlag == -12 || CstateFlag == -2 || CstateFlag == -3 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu_Rot < 0.0 && d_Rot <= R_dresneg && d_Rot > R_delUneg && d_Rot <= dpeakmin && BenMark > 0)
		return -3;
	else if ((CstateFlag == -12 || CstateFlag == -2 || CstateFlag == -3 || CstateFlag == 5 || CstateFlag == -6 || CstateFlag == 7 || CstateFlag == -5 || CstateFlag == -4) && Tdu_Rot < 0.0 && d_Rot <= R_dresneg && d_Rot > R_delUneg && d_Rot <= dpeakmin && BenMark < 0)
		return -3;
	else if ((CstateFlag == -3 || CstateFlag == -30 || CstateFlag == 31) && Tdu_Rot < 0.0 && d_Rot <= R_delUneg && d_Rot >= Intcpt_Xaxis_neg && d_Rot <= dpeakmin && BenMark > 0)
		return -30;
	else if ((CstateFlag == -3 || CstateFlag == -30 || CstateFlag == 31) && Tdu_Rot < 0.0 && d_Rot <= R_delUneg && d_Rot >= Intcpt_Xaxis_neg && d_Rot <= dpeakmin && BenMark < 0)
		return -30;
	else if ((CstateFlag == -30 || CstateFlag == -40) && Tdu_Rot < 0.0 && d_Rot < Intcpt_Xaxis_neg && d_Rot <= dpeakmin && BenMark > 0)
		return -40;
	else if ((CstateFlag == -30 || CstateFlag == -40) && Tdu_Rot < 0.0 && d_Rot < Intcpt_Xaxis_neg && d_Rot <= dpeakmin && BenMark < 0)
		return -40;
	else if ((CstateFlag == -30 || CstateFlag == 31) && (Tdu_Rot > 0.0 || Tdu_Rot < 0.0) && (d_Rot > dpeakmin || d_Rot < dpeakmax))
		return 31;
	else if ((CstateFlag == -40 || CstateFlag == 41) && (Tdu_Rot > 0.0 || Tdu_Rot < 0.0))
		return 41;

	else if ((CstateFlag == -12 || CstateFlag == -2 || CstateFlag == -3) && Tdu_Rot > 0.0)
		return -4;
	else if ((CstateFlag == -4 || CstateFlag == -5) && Tdu_Rot > 0.0 && d_Rot < dpeakmax)
		return -5;

	else if ((CstateFlag == -4 || CstateFlag == -5 || CstateFlag == -6) && Tdu_Rot < 0.0 && d_Rot > dpeakmin)
		return -6;
	else if ((CstateFlag == -6 || CstateFlag == -7 || CstateFlag == 7) && Tdu_Rot > 0.0 && d_Rot < dpeakmax)
		return -7;
	else
		return 999;

}

/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Shear Material                                         **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */

int
GMG_CMAC2D::getStateFlag_Shr(void)
{

	if ((CstateFlag_Shr == 0 || CstateFlag_Shr == 12 || CstateFlag_Shr == -5 || CstateFlag_Shr == 6 || CstateFlag_Shr == -7 || CstateFlag_Shr == 5 || CstateFlag_Shr == 4) && Tdu_Shr > 0.0 && d_Shr >= dypos_Shr && d_Shr >= dpeakmax_Shr && d_Shr < R_dcappos_Shr && BenMark_Shr > 0)
		return 12;
	else if ((CstateFlag_Shr == 0 || CstateFlag_Shr == 12 || CstateFlag_Shr == -5 || CstateFlag_Shr == 6 || CstateFlag_Shr == -7 || CstateFlag_Shr == 5 || CstateFlag_Shr == 4) && Tdu_Shr > 0.0 && d_Shr >= dypos_Shr && d_Shr >= dpeakmax_Shr && d_Shr < R_dcappos_Shr && BenMark_Shr < 0)
		return 12;
	else if ((CstateFlag_Shr == 12 || CstateFlag_Shr == 2 || CstateFlag_Shr == -5 || CstateFlag_Shr == 6 || CstateFlag_Shr == -7 || CstateFlag_Shr == 5 || CstateFlag_Shr == 4) && Tdu_Shr > 0.0 && d_Shr >= R_dcappos_Shr && d_Shr >= dpeakmax_Shr && d_Shr < R_drespos_Shr && BenMark_Shr > 0)
		return 2;
	else if ((CstateFlag_Shr == 12 || CstateFlag_Shr == 2 || CstateFlag_Shr == -5 || CstateFlag_Shr == 6 || CstateFlag_Shr == -7 || CstateFlag_Shr == 5 || CstateFlag_Shr == 4) && Tdu_Shr > 0.0 && d_Shr >= R_dcappos_Shr && d_Shr >= dpeakmax_Shr && d_Shr < R_drespos_Shr && BenMark_Shr < 0)
		return 2;
	else if ((CstateFlag_Shr == 12 || CstateFlag_Shr == 2 || CstateFlag_Shr == 3 || CstateFlag_Shr == -5 || CstateFlag_Shr == 6 || CstateFlag_Shr == -7 || CstateFlag_Shr == 5 || CstateFlag_Shr == 4) && Tdu_Shr > 0.0 && d_Shr >= R_drespos_Shr && d_Shr < delUpos_Shr && d_Shr >= dpeakmax_Shr && BenMark_Shr > 0)
		return 3;
	else if ((CstateFlag_Shr == 12 || CstateFlag_Shr == 2 || CstateFlag_Shr == 3 || CstateFlag_Shr == -5 || CstateFlag_Shr == 6 || CstateFlag_Shr == -7 || CstateFlag_Shr == 5 || CstateFlag_Shr == 4) && Tdu_Shr > 0.0 && d_Shr >= R_drespos_Shr && d_Shr < delUpos_Shr && d_Shr >= dpeakmax_Shr && BenMark_Shr < 0)
		return 3;

	else if ((CstateFlag_Shr == 3 || CstateFlag_Shr == 30 || CstateFlag_Shr == 31) && Tdu_Shr > 0.0 && d_Shr >= delUpos_Shr && d_Shr <= Intcpt_Xaxis_pos_Shr && d_Shr >= dpeakmax_Shr && BenMark_Shr > 0)
		return 30;
	else if ((CstateFlag_Shr == 3 || CstateFlag_Shr == 30 || CstateFlag_Shr == 31) && Tdu_Shr > 0.0 && d_Shr >= delUpos_Shr && d_Shr <= Intcpt_Xaxis_pos_Shr && d_Shr >= dpeakmax_Shr && BenMark_Shr < 0)
		return 30;
	else if ((CstateFlag_Shr == 30 || CstateFlag_Shr == 40) && Tdu_Shr > 0.0 && d_Shr >= Intcpt_Xaxis_pos_Shr && d_Shr >= dpeakmax_Shr && BenMark_Shr > 0)
		return 40;
	else if ((CstateFlag_Shr == 30 || CstateFlag_Shr == 40) && Tdu_Shr > 0.0 && d_Shr >= Intcpt_Xaxis_pos_Shr && d_Shr >= dpeakmax_Shr && BenMark_Shr < 0)
		return 40;
	else if ((CstateFlag_Shr == 30 || CstateFlag_Shr == 31) && (Tdu_Shr > 0.0 || Tdu_Shr < 0.0) && (d_Shr > dpeakmin_Shr || d_Shr < dpeakmax_Shr))
		return 31;
	else if ((CstateFlag_Shr == 40 || CstateFlag_Shr == 41) && (Tdu_Shr > 0.0 || Tdu_Shr < 0.0))
		return 41;


	else if ((CstateFlag_Shr == 12 || CstateFlag_Shr == 2 || CstateFlag_Shr == 3) && Tdu_Shr < 0.0)
		return 4;
	else if ((CstateFlag_Shr == 4 || CstateFlag_Shr == 5) && Tdu_Shr < 0.0 && d_Shr > dpeakmin_Shr)
		return 5;
	else if ((CstateFlag_Shr == 4 || CstateFlag_Shr == 5 || CstateFlag_Shr == 6) && Tdu_Shr > 0.0 && d_Shr < dpeakmax_Shr)
		return 6;
	else if ((CstateFlag_Shr == 6 || CstateFlag_Shr == 7 || CstateFlag_Shr == -7) && Tdu_Shr < 0.0 && d_Shr > dpeakmin_Shr)
		return 7;

	else if ((CstateFlag_Shr == 0 || CstateFlag_Shr == -12 || CstateFlag_Shr == 5 || CstateFlag_Shr == -6 || CstateFlag_Shr == 7 || CstateFlag_Shr == -5 || CstateFlag_Shr == -4) && Tdu_Shr < 0.0 && d_Shr <= dpeakmin_Shr && d_Shr > R_dcapneg_Shr && BenMark_Shr > 0)
		return -12;
	else if ((CstateFlag_Shr == 0 || CstateFlag_Shr == -12 || CstateFlag_Shr == 5 || CstateFlag_Shr == -6 || CstateFlag_Shr == 7 || CstateFlag_Shr == -5 || CstateFlag_Shr == -4) && Tdu_Shr < 0.0 && d_Shr <= dpeakmin_Shr && d_Shr > R_dcapneg_Shr && BenMark_Shr < 0)
		return -12;
	else if ((CstateFlag_Shr == -12 || CstateFlag_Shr == -2 || CstateFlag_Shr == 5 || CstateFlag_Shr == -6 || CstateFlag_Shr == 7 || CstateFlag_Shr == -5 || CstateFlag_Shr == -4) && Tdu_Shr < 0.0 && d_Shr <= R_dcapneg_Shr && d_Shr <= dpeakmin_Shr && d_Shr > R_dresneg_Shr && BenMark_Shr > 0)
		return -2;
	else if ((CstateFlag_Shr == -12 || CstateFlag_Shr == -2 || CstateFlag_Shr == 5 || CstateFlag_Shr == -6 || CstateFlag_Shr == 7 || CstateFlag_Shr == -5 || CstateFlag_Shr == -4) && Tdu_Shr < 0.0 && d_Shr <= R_dcapneg_Shr && d_Shr <= dpeakmin_Shr && d_Shr > R_dresneg_Shr && BenMark_Shr < 0)
		return -2;
	else if ((CstateFlag_Shr == -12 || CstateFlag_Shr == -2 || CstateFlag_Shr == -3 || CstateFlag_Shr == 5 || CstateFlag_Shr == -6 || CstateFlag_Shr == 7 || CstateFlag_Shr == -5 || CstateFlag_Shr == -4) && Tdu_Shr < 0.0 && d_Shr <= R_dresneg_Shr && d_Shr > delUneg_Shr && d_Shr <= dpeakmin_Shr && BenMark_Shr > 0)
		return -3;
	else if ((CstateFlag_Shr == -12 || CstateFlag_Shr == -2 || CstateFlag_Shr == -3 || CstateFlag_Shr == 5 || CstateFlag_Shr == -6 || CstateFlag_Shr == 7 || CstateFlag_Shr == -5 || CstateFlag_Shr == -4) && Tdu_Shr < 0.0 && d_Shr <= R_dresneg_Shr && d_Shr > delUneg_Shr && d_Shr <= dpeakmin_Shr && BenMark_Shr < 0)
		return -3;
	else if ((CstateFlag_Shr == -3 || CstateFlag_Shr == -30 || CstateFlag_Shr == 31) && Tdu_Shr < 0.0 && d_Shr <= delUneg_Shr && d_Shr >= Intcpt_Xaxis_neg_Shr && d_Shr <= dpeakmin_Shr && BenMark_Shr > 0)
		return -30;
	else if ((CstateFlag_Shr == -3 || CstateFlag_Shr == -30 || CstateFlag_Shr == 31) && Tdu_Shr < 0.0 && d_Shr <= delUneg_Shr && d_Shr >= Intcpt_Xaxis_neg_Shr && d_Shr <= dpeakmin_Shr && BenMark_Shr < 0)
		return -30;
	else if ((CstateFlag_Shr == -30 || CstateFlag_Shr == -40) && Tdu_Shr < 0.0 && d_Shr < Intcpt_Xaxis_neg_Shr && d_Shr <= dpeakmin_Shr && BenMark_Shr > 0)
		return -40;
	else if ((CstateFlag_Shr == -30 || CstateFlag_Shr == -40) && Tdu_Shr < 0.0 && d_Shr < Intcpt_Xaxis_neg_Shr && d_Shr <= dpeakmin_Shr && BenMark_Shr < 0)
		return -40;
	else if ((CstateFlag_Shr == -30 || CstateFlag_Shr == 31) && (Tdu_Shr > 0.0 || Tdu_Shr < 0.0) && (d_Shr > dpeakmin_Shr || d_Shr < dpeakmax_Shr))
		return 31;
	else if ((CstateFlag_Shr == -40 || CstateFlag_Shr == 41) && (Tdu_Shr > 0.0 || Tdu_Shr < 0.0))
		return 41;

	else if ((CstateFlag_Shr == -12 || CstateFlag_Shr == -2 || CstateFlag_Shr == -3) && Tdu_Shr > 0.0)
		return -4;
	else if ((CstateFlag_Shr == -4 || CstateFlag_Shr == -5) && Tdu_Shr > 0.0 && d_Shr < dpeakmax_Shr)
		return -5;

	else if ((CstateFlag_Shr == -4 || CstateFlag_Shr == -5 || CstateFlag_Shr == -6) && Tdu_Shr < 0.0 && d_Shr > dpeakmin_Shr)
		return -6;
	else if ((CstateFlag_Shr == -6 || CstateFlag_Shr == -7 || CstateFlag_Shr == 7) && Tdu_Shr > 0.0 && d_Shr < dpeakmax_Shr)
		return -7;
	else
		return 999;

}


/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Rotational Material                                    **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */

//Defining the peaks of the splin curve
void
GMG_CMAC2D::define_peak(void)
{
	double ffmax_effective, ffmin_effective;
	double delta_R_fcappos, delta_R_fcapneg;

	/* ************************************************************ **
	* ************************************************************ **

			  Post Yielding (before capping) Damage Mode
				  First Part, Backbone to Backbone

	* ************************************************************ **
	* ************************************************************ **/

	if (CstateFlag == 12 && TstateFlag == 4) {
		dpeakmax = Cstrain_Rot;
		ffmax = Cstress_Rot;

		if (flagdmg_Hardening_strength == 1) {

			delta_R_fcapneg = dmgSneg * fabs(R_fcapneg);
			R_dcapneg = CR_dcapneg + (fabs(R_Kdegneg) / (fabs(slope_neg) + fabs(R_Kdegneg))) * (fabs(delta_R_fcapneg) / fabs(R_Kdegneg));  //Updating dcapneg
			R_fcapneg = -(slope_neg * fabs(R_dcapneg) - Intcpt_slope_neg);
			R_dresneg = CR_dresneg + fabs(delta_R_fcapneg) / fabs(R_Kdegneg); //Updating dresneg
			Intcpt_deg_neg = CIntcpt_deg_neg - fabs(delta_R_fcapneg);
		}

		if (flagdmg_Hardening == 1) {

			dpeakmax_bench = dpeakmax * (1 + delta_ratio_max_hard_Rot);

			if (dpeakmin > dpeakmin_bench) {
				dpeakmin = fmax((1 + alpha_neg) * dpeakmin, dpeakmin_bench);
				ffmin = -(slope_neg * fabs(dpeakmin) - Intcpt_slope_neg);
				if (dpeakmin < R_dcapneg) {
					ffmin = -(R_Kdegneg * fabs(dpeakmin) + Intcpt_deg_neg);
				}
			}
			else if (dpeakmin <= dpeakmin_bench) {
				dpeakmin = dpeakmin_bench;
				ffmin = -(slope_neg * fabs(dpeakmin) - Intcpt_slope_neg);
				if (dpeakmin < R_dcapneg) {
					ffmin = -(R_Kdegneg * fabs(dpeakmin) + Intcpt_deg_neg);
				}
			}
		}

	}
	if (CstateFlag == -12 && TstateFlag == -4) {
		dpeakmin = Cstrain_Rot;
		ffmin = Cstress_Rot;

		if (flagdmg_Hardening_strength == 1) {

			delta_R_fcappos = dmgSpos * R_fcappos;
			R_dcappos = CR_dcappos - fabs(R_Kdegpos) / (fabs(slope_pos) + fabs(R_Kdegpos))*(fabs(delta_R_fcappos) / fabs(R_Kdegpos));
			R_fcappos = slope_pos * fabs(R_dcappos) + Intcpt_slope_pos;
			R_drespos = CR_drespos - fabs(delta_R_fcappos) / fabs(R_Kdegpos);
			Intcpt_deg_pos = CIntcpt_deg_pos - fabs(delta_R_fcappos);
		}

		if (flagdmg_Hardening == 1) {
			dpeakmin_bench = dpeakmin * (1 + delta_ratio_max_hard_Rot);

			if ((dpeakmax) < dpeakmax_bench) {
				dpeakmax = fmin((1 + alpha_pos) * dpeakmax, dpeakmax_bench);
				ffmax = slope_pos * dpeakmax + Intcpt_slope_pos;
				if (dpeakmax > R_dcappos) {
					ffmax = R_Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
				}
			}
			else if ((dpeakmax) >= dpeakmax_bench) {
				dpeakmax = dpeakmax_bench;
				ffmax = slope_pos * dpeakmax + Intcpt_slope_pos;
				if (dpeakmax > R_dcappos) {
					ffmax = R_Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
				}
			}
		}

	}

	/* ************************************************************ **
	 * ************************************************************ **

				  Post Yielding before capping Damage Mode
						Second Part, inner cyclic

	 * ************************************************************ **
	 * ************************************************************ **/
	if (flagdmg_Hardening_strength == 1) {
		if ((CstateFlag == -4 && TstateFlag == -6) || (CstateFlag == -5 && TstateFlag == -6) || (CstateFlag == 6 && TstateFlag == 7) || (CstateFlag == -7 && TstateFlag == 7)) {

			delta_R_fcapneg = dmgSneg * fabs(R_fcapneg);
			R_dcapneg = CR_dcapneg + (fabs(R_Kdegneg) / (fabs(slope_neg) + fabs(R_Kdegneg))) * (fabs(delta_R_fcapneg) / fabs(R_Kdegneg));  //Updating dcapneg
			R_fcapneg = -(slope_neg * fabs(R_dcapneg) - Intcpt_slope_neg);
			R_dresneg = CR_dresneg + fabs(delta_R_fcapneg) / fabs(R_Kdegneg); //Updating dresneg
			Intcpt_deg_neg = CIntcpt_deg_neg - fabs(delta_R_fcapneg);
			if (dpeakmin < R_dcapneg) {
				ffmin = -(R_Kdegneg * fabs(dpeakmin) + Intcpt_deg_neg);
			}
		}
	}


	if (flagdmg_Hardening_strength == 1) {

		if ((CstateFlag == 4 && TstateFlag == 6) || (CstateFlag == 5 && TstateFlag == 6) || (CstateFlag == -6 && TstateFlag == -7) || (CstateFlag == 7 && TstateFlag == -7)) {

			delta_R_fcappos = dmgSpos * R_fcappos;
			R_dcappos = CR_dcappos - fabs(R_Kdegpos) / (fabs(slope_pos) + fabs(R_Kdegpos))*(fabs(delta_R_fcappos) / fabs(R_Kdegpos));
			R_fcappos = slope_pos * fabs(R_dcappos) + Intcpt_slope_pos;
			R_drespos = CR_drespos - fabs(delta_R_fcappos) / fabs(R_Kdegpos);
			Intcpt_deg_pos = CIntcpt_deg_pos - fabs(delta_R_fcappos);
			if (dpeakmax > R_dcappos) {
				ffmax = R_Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
			}

		}
	}

	if (flagdmg_Hardening == 1 && (dpeakmax < dpeakmax_bench)) {

		if ((CstateFlag == 4 && TstateFlag == 6) || (CstateFlag == 5 && TstateFlag == 6) || (CstateFlag == -6 && TstateFlag == -7) || (CstateFlag == 7 && TstateFlag == -7)) {
			dpeakmax = fmin((1 + alpha_pos) * dpeakmax, dpeakmax_bench);
			ffmax = slope_pos * dpeakmax + Intcpt_slope_pos;
			if (dpeakmax > R_dcappos) {
				ffmax = R_Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
			}
		}
	}

	if (flagdmg_Hardening == 1 && (dpeakmin > dpeakmin_bench)) {
		if ((CstateFlag == -4 && TstateFlag == -6) || (CstateFlag == -5 && TstateFlag == -6) || (CstateFlag == 6 && TstateFlag == 7) || (CstateFlag == -7 && TstateFlag == 7)) {
			dpeakmin = fmax((1 + alpha_neg) * dpeakmin, dpeakmin_bench);
			ffmin = -(slope_neg * fabs(dpeakmin) - Intcpt_slope_neg);
			if (dpeakmin < R_dcapneg) {
				ffmin = -(R_Kdegneg * fabs(dpeakmin) + Intcpt_deg_neg);
			}
		}
	}

	/* ************************************************************ **
			 Dealing with PM rule in the post yielding region,
					   updating ffmax and ffmin
	* ************************************************************ **/
	if (flag_entering_hardening_Flex_pos_Rot == 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {
		ffmax = slope_pos * dpeakmax + Intcpt_slope_pos;
	}
	if (flag_entering_hardening_Flex_neg_Rot == 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {
		ffmin = -(slope_neg * fabs(dpeakmin) - Intcpt_slope_neg);
	}

	if (flag_entering_hardening_Splice_pos_Rot == 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {
		ffmax = slope_pos * dpeakmax + Intcpt_slope_pos;
	}
	if (flag_entering_hardening_Splice_neg_Rot == 1 && flag_Flx_Filure_Rot != 1 && flag_FlexShr_Filure_Rot != 1 && flag_Splice_Filure_Rot != 1 && flag_FlexSplice_Filure_Rot != 1) {
		ffmin = -(slope_neg * fabs(dpeakmin) - Intcpt_slope_neg);
	}

	/* ************************************************************ **
	 * ************************************************************ **

					 End of Post Yielding Damage Mode

	 * ************************************************************ **
	 * ************************************************************ **/

	 /* ************************************************************ **
	 * ************************************************************ **

				  Post Capping (before residual) Damage Mode
					  First Part, Backbone to Backbone

	 * ************************************************************ **
	 * ************************************************************ **/

	if (CstateFlag == 2 && TstateFlag == 4) {
		dpeakmax = Cstrain_Rot;
		ffmax = Cstress_Rot;

		ratio = (dpeakmax - R_dcappos) / (R_drespos - R_dcappos);
		dpeakmin = (R_dresneg - R_dcapneg) * ratio + R_dcapneg;

		ffmin_effective = -(R_Kdegneg * fabs(dpeakmin) + Intcpt_deg_neg);
		ffmin = ffmin_effective;

		if (dmgSneg > 0.0) {

			delta_R_fcapneg = dmgSneg * fabs(R_fcapneg);
			ffmin = fmin((ffmin_effective + delta_R_fcapneg), R_fresneg);
			R_dcapneg = CR_dcapneg + (fabs(R_Kdegneg) / (fabs(slope_neg) + fabs(R_Kdegneg))) * (fabs(delta_R_fcapneg) / fabs(R_Kdegneg));  //Updating dcapneg
			R_dresneg = CR_dresneg + fabs(delta_R_fcapneg) / fabs(R_Kdegneg); //Updating dresneg
			Intcpt_deg_neg = CIntcpt_deg_neg - fabs(delta_R_fcapneg);
		}

	}
	if (CstateFlag == -2 && TstateFlag == -4) {
		dpeakmin = Cstrain_Rot;
		ffmin = Cstress_Rot;

		ratio = (dpeakmin - R_dcapneg) / (R_dresneg - R_dcapneg);
		dpeakmax = (R_drespos - R_dcappos) * ratio + R_dcappos;

		ffmax_effective = R_Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
		ffmax = ffmax_effective;

		if (dmgSpos > 0.0) {

			delta_R_fcappos = dmgSpos * R_fcappos;
			ffmax = fmax(ffmax_effective - delta_R_fcappos, R_frespos);
			R_drespos = CR_drespos - fabs(delta_R_fcappos) / fabs(R_Kdegpos);
			R_dcappos = CR_dcappos - fabs(R_Kdegpos) / (fabs(slope_pos) + fabs(R_Kdegpos))*(fabs(delta_R_fcappos) / fabs(R_Kdegpos));
			Intcpt_deg_pos = CIntcpt_deg_pos - fabs(delta_R_fcappos);

		}
	}

	/* ************************************************************ **
	 * ************************************************************ **

			 Post Capping (before residual) Damage Mode
					Second Part, Inner cyclic

	 * ************************************************************ **
	 * ************************************************************ **/

	if (flagdmg == 1) {
		if ((CstateFlag == 4 && TstateFlag == 6) || (CstateFlag == 5 && TstateFlag == 6) || (CstateFlag == -6 && TstateFlag == -7) || (CstateFlag == 7 && TstateFlag == -7)) {

			ffmax_effective = R_Kdegpos * fabs(dpeakmax) + Intcpt_deg_pos;
			ffmax = ffmax_effective;
			if (dmgSpos > 0.0) {

				delta_R_fcappos = dmgSpos * R_fcappos;
				ffmax = fmax(ffmax_effective - delta_R_fcappos, R_frespos);
				R_drespos = CR_drespos - fabs(delta_R_fcappos) / fabs(R_Kdegpos);
				R_dcappos = CR_dcappos - fabs(R_Kdegpos) / (fabs(slope_pos) + fabs(R_Kdegpos))*(fabs(delta_R_fcappos) / fabs(R_Kdegpos));
				Intcpt_deg_pos = CIntcpt_deg_pos - fabs(delta_R_fcappos);
			}
		}
	}


	if (flagdmg == 1) {
		if ((CstateFlag == -4 && TstateFlag == -6) || (CstateFlag == -5 && TstateFlag == -6) || (CstateFlag == 6 && TstateFlag == 7) || (CstateFlag == -7 && TstateFlag == 7)) {


			ffmin_effective = -(R_Kdegneg * fabs(dpeakmin) + Intcpt_deg_neg);
			ffmin = ffmin_effective;
			if (dmgSneg > 0.0) {

				delta_R_fcapneg = dmgSneg * fabs(R_fcapneg);
				ffmin = fmin((ffmin_effective + delta_R_fcapneg), R_fresneg);
				R_dcapneg = CR_dcapneg + (fabs(R_Kdegneg) / (fabs(slope_neg) + fabs(R_Kdegneg))) * (fabs(delta_R_fcapneg) / fabs(R_Kdegneg));  //Updating dcapneg
				R_dresneg = CR_dresneg + fabs(delta_R_fcapneg) / fabs(R_Kdegneg); //Updating dresneg
				Intcpt_deg_neg = CIntcpt_deg_neg - fabs(delta_R_fcapneg);

			}

		}
	}

	/* ************************************************************ **
	 * ************************************************************ **

					 End of Post Capping Damage Mode

	 * ************************************************************ **
	 * ************************************************************ **/

	 /* ************************************************************ **
	  * ************************************************************ **

					Post Capping (After residual)

	 * ************************************************************ **
	  * ************************************************************ **/

	if (CstateFlag == 3 && TstateFlag == 4) {
		dpeakmax = Cstrain_Rot;
		ffmax = Cstress_Rot;
		ratio = (dpeakmax - R_drespos) / (R_delUpos - R_drespos);
		dpeakmin = (R_delUneg - R_dresneg) * ratio + R_dresneg;
		ffmin = R_fresneg;
	}

	if (CstateFlag == -3 && TstateFlag == -4) {
		dpeakmin = Cstrain_Rot;
		ffmin = Cstress_Rot;

		ratio = (dpeakmin - R_dresneg) / (R_delUneg - R_dresneg);
		dpeakmax = (R_delUpos - R_drespos) * ratio + R_drespos;
		ffmax = R_frespos;
	}

	if (CstateFlag == 30 && TstateFlag == 31) {
		dpeakmax = Cstrain_Rot;
		ffmax = Cstress_Rot;
		ratio = (dpeakmax - R_delUpos) / (Intcpt_Xaxis_pos - R_delUpos);
		dpeakmin = (Intcpt_Xaxis_neg - R_delUneg) * ratio + R_delUneg;
		ffmin = -(R_Kdegneg * fabs(dpeakmin) + Intcpt_res_neg);
	}

	if (CstateFlag == -30 && TstateFlag == 31) {
		dpeakmin = Cstrain_Rot;
		ffmin = Cstress_Rot;
		ratio = (dpeakmin - R_delUneg) / (Intcpt_Xaxis_neg - R_delUneg);
		dpeakmax = (Intcpt_Xaxis_pos - R_delUpos) * ratio + R_delUpos;
		ffmax = R_Kdegpos * fabs(dpeakmax) + Intcpt_res_pos;
	}
	//--------------------------------------End of Part 2-----------------------------------------------
}
//-----------------------------My subrotin----------------------------------

/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Shear Material                                         **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */


void
GMG_CMAC2D::define_peak_Shr(void)
{
	double ffmax_effective, ffmin_effective;
	double delta_R_fcappos_Shr, delta_R_fcapneg_Shr;

	/* ************************************************************ **
	* ************************************************************ **

			 Post Yielding (before capping) Damage Mode
				  First Part, Backbone to Backbone

	* ************************************************************ **
	* ************************************************************ **/

	if (CstateFlag_Shr == 12 && TstateFlag_Shr == 4) {
		dpeakmax_Shr = Cstrain_Shr;
		ffmax_Shr = Cstress_Shr;

		if (flagdmg_Hardening_strength_Shr == 1) {

			delta_R_fcapneg_Shr = dmgSneg_Shr * fabs(fcapneg_Shr);
			R_dcapneg_Shr = CR_dcapneg_Shr + fabs(delta_R_fcapneg_Shr) / (fabs(slope_neg_Shr) + fabs(R_Kdegneg_Shr));
			R_fcapneg_Shr = -(slope_neg_Shr * fabs(R_dcapneg_Shr) - Intcpt_slope_neg_Shr);
			R_dresneg_Shr = CR_dresneg_Shr + fabs(delta_R_fcapneg_Shr) / fabs(R_Kdegneg_Shr); //Updating dresneg
			Intcpt_deg_neg_Shr = CIntcpt_deg_neg_Shr - fabs(delta_R_fcapneg_Shr);

		}

		if (flagdmg_Hardening_Shr == 1) {
			dpeakmax_bench_Shr = dpeakmax_Shr * (1 + delta_ratio_max_hard_Shr);

			if (dpeakmin_Shr > dpeakmin_bench_Shr) {
				dpeakmin_Shr = fmax((1 + alpha_neg_Shr) * dpeakmin_Shr, dpeakmin_bench_Shr);
				ffmin_Shr = -(slope_neg_Shr * fabs(dpeakmin_Shr) - Intcpt_slope_neg_Shr);
				if (dpeakmin_Shr < R_dcapneg_Shr) {
					ffmin_Shr = -(R_Kdegneg_Shr * fabs(dpeakmin_Shr) + Intcpt_deg_neg_Shr);
				}
			}
			else if (dpeakmin_Shr <= dpeakmin_bench_Shr) {
				dpeakmin_Shr = dpeakmin_bench_Shr;
				ffmin_Shr = -(slope_neg_Shr * fabs(dpeakmin_Shr) - Intcpt_slope_neg_Shr);
				if (dpeakmin_Shr < R_dcapneg_Shr) {
					ffmin_Shr = -(R_Kdegneg_Shr * fabs(dpeakmin_Shr) + Intcpt_deg_neg_Shr);
				}
			}
		}
	}
	if (CstateFlag_Shr == -12 && TstateFlag_Shr == -4) {
		dpeakmin_Shr = Cstrain_Shr;
		ffmin_Shr = Cstress_Shr;


		if (flagdmg_Hardening_strength_Shr == 1) {

			delta_R_fcappos_Shr = dmgSpos_Shr * fcappos_Shr;
			R_dcappos_Shr = CR_dcappos_Shr - fabs(delta_R_fcappos_Shr) / (fabs(slope_pos_Shr) + fabs(Kdegpos_Shr));
			R_fcappos_Shr = slope_pos_Shr * fabs(R_dcappos_Shr) + Intcpt_slope_pos_Shr;
			R_drespos_Shr = CR_drespos_Shr - fabs(delta_R_fcappos_Shr) / fabs(Kdegpos_Shr);
			Intcpt_deg_pos_Shr = CIntcpt_deg_pos_Shr - fabs(delta_R_fcappos_Shr);

		}

		if (flagdmg_Hardening_Shr == 1) {
			dpeakmin_bench_Shr = dpeakmin_Shr * (1 + delta_ratio_max_hard_Shr);

			if ((dpeakmax_Shr) < dpeakmax_bench_Shr) {
				dpeakmax_Shr = fmin((1 + alpha_pos_Shr) * dpeakmax_Shr, dpeakmax_bench_Shr);
				ffmax_Shr = slope_pos_Shr * dpeakmax_Shr + Intcpt_slope_pos_Shr;
				if (dpeakmax_Shr > R_dcappos_Shr) {
					ffmax_Shr = Kdegpos_Shr * fabs(dpeakmax_Shr) + Intcpt_deg_pos_Shr;
				}
			}
			else if ((dpeakmax_Shr) >= dpeakmax_bench_Shr) {
				dpeakmax_Shr = dpeakmax_bench_Shr;
				ffmax_Shr = slope_pos_Shr * dpeakmax_Shr + Intcpt_slope_pos_Shr;
				if (dpeakmax_Shr > R_dcappos_Shr) {
					ffmax_Shr = Kdegpos_Shr * fabs(dpeakmax_Shr) + Intcpt_deg_pos_Shr;
				}
			}
		}

	}


	/* ************************************************************ **
	 * ************************************************************ **

				  Post Yielding before capping Damage Mode
						Second Part, inner cyclic

	 * ************************************************************ **
	 * ************************************************************ **/
	if (flagdmg_Hardening_strength_Shr == 1) {
		if ((CstateFlag_Shr == -4 && TstateFlag_Shr == -6) || (CstateFlag_Shr == -5 && TstateFlag_Shr == -6) || (CstateFlag_Shr == 6 && TstateFlag_Shr == 7) || (CstateFlag_Shr == -7 && TstateFlag_Shr == 7)) {

			delta_R_fcapneg_Shr = dmgSneg_Shr * fabs(fcapneg_Shr);
			R_dcapneg_Shr = CR_dcapneg_Shr + fabs(delta_R_fcapneg_Shr) / (fabs(slope_neg_Shr) + fabs(R_Kdegneg_Shr));
			R_fcapneg_Shr = -(slope_neg_Shr * fabs(R_dcapneg_Shr) - Intcpt_slope_neg_Shr);
			R_dresneg_Shr = CR_dresneg_Shr + fabs(delta_R_fcapneg_Shr) / fabs(R_Kdegneg_Shr); //Updating dresneg
			Intcpt_deg_neg_Shr = CIntcpt_deg_neg_Shr - fabs(delta_R_fcapneg_Shr);
			if (dpeakmin_Shr < R_dcapneg_Shr) {
				ffmin_Shr = -(R_Kdegneg_Shr * fabs(dpeakmin_Shr) + Intcpt_deg_neg_Shr);
			}
		}
	}

	if (flagdmg_Hardening_strength_Shr == 1) {
		if ((CstateFlag_Shr == 4 && TstateFlag_Shr == 6) || (CstateFlag_Shr == 5 && TstateFlag_Shr == 6) || (CstateFlag_Shr == -6 && TstateFlag_Shr == -7) || (CstateFlag_Shr == 7 && TstateFlag_Shr == -7)) {

			delta_R_fcappos_Shr = dmgSpos_Shr * fcappos_Shr;
			R_dcappos_Shr = CR_dcappos_Shr - fabs(delta_R_fcappos_Shr) / (fabs(slope_pos_Shr) + fabs(Kdegpos_Shr));
			R_fcappos_Shr = slope_pos_Shr * fabs(R_dcappos_Shr) + Intcpt_slope_pos_Shr;
			R_drespos_Shr = CR_drespos_Shr - fabs(delta_R_fcappos_Shr) / fabs(Kdegpos_Shr);
			Intcpt_deg_pos_Shr = CIntcpt_deg_pos_Shr - fabs(delta_R_fcappos_Shr);
			if (dpeakmax_Shr > R_dcappos_Shr) {
				ffmax_Shr = Kdegpos_Shr * fabs(dpeakmax_Shr) + Intcpt_deg_pos_Shr;
			}
		}
	}

	if (flagdmg_Hardening_Shr == 1 && (dpeakmin_Shr > dpeakmin_bench_Shr)) {
		if ((CstateFlag_Shr == -4 && TstateFlag_Shr == -6) || (CstateFlag_Shr == -5 && TstateFlag_Shr == -6) || (CstateFlag_Shr == 6 && TstateFlag_Shr == 7) || (CstateFlag_Shr == -7 && TstateFlag_Shr == 7)) {
			dpeakmin_Shr = fmax((1 + alpha_neg_Shr) * dpeakmin_Shr, dpeakmin_bench_Shr);
			ffmin_Shr = -(slope_neg_Shr * fabs(dpeakmin_Shr) - Intcpt_slope_neg_Shr);
			if (dpeakmin_Shr < R_dcapneg_Shr) {
				ffmin_Shr = -(R_Kdegneg_Shr * fabs(dpeakmin_Shr) + Intcpt_deg_neg_Shr);
			}
		}
	}

	if (flagdmg_Hardening_Shr == 1 && (dpeakmax_Shr < dpeakmax_bench_Shr)) {
		if ((CstateFlag_Shr == 4 && TstateFlag_Shr == 6) || (CstateFlag_Shr == 5 && TstateFlag_Shr == 6) || (CstateFlag_Shr == -6 && TstateFlag_Shr == -7) || (CstateFlag_Shr == 7 && TstateFlag_Shr == -7)) {
			dpeakmax_Shr = fmin((1 + alpha_pos_Shr) * dpeakmax_Shr, dpeakmax_bench_Shr);
			ffmax_Shr = slope_pos_Shr * dpeakmax_Shr + Intcpt_slope_pos_Shr;
			if (dpeakmax_Shr > R_dcappos_Shr) {
				ffmax_Shr = Kdegpos_Shr * fabs(dpeakmax_Shr) + Intcpt_deg_pos_Shr;
			}
		}
	}

	/* ************************************************************ **
		 Dealing with PM rule in the post yielding region,
				   updating ffmax and ffmin
* ************************************************************ **/
	if (flag_entering_hardening_Shr_pos == 1 && flag_Filure_Shr != 1) {
		ffmax_Shr = slope_pos_Shr * dpeakmax_Shr + Intcpt_slope_pos_Shr;
	}
	if (flag_entering_hardening_Shr_neg == 1 && flag_Filure_Shr != 1) {
		ffmin_Shr = -(slope_neg_Shr * fabs(dpeakmin_Shr) - Intcpt_slope_neg_Shr);
	}
	/* ************************************************************ **
	 * ************************************************************ **

					 End of Post Yielding Damage Mode

	 * ************************************************************ **
	 * ************************************************************ **/

	if (CstateFlag_Shr == 2 && TstateFlag_Shr == 4) {
		dpeakmax_Shr = Cstrain_Shr;
		ffmax_Shr = Cstress_Shr;
		//---------------------Changing this part to see if it is gonna correct the difference between the dpeakmin and depeakmax----------------------
		ratio_Shr = (dpeakmax_Shr - R_dcappos_Shr) / (R_drespos_Shr - R_dcappos_Shr);
		dpeakmin_Shr = (R_dresneg_Shr - R_dcapneg_Shr) * ratio_Shr + R_dcapneg_Shr;
		//---------------------------------------------------------------------------------------------------------------------------------------------
		ffmin_effective = -(R_Kdegneg_Shr*fabs(dpeakmin_Shr) + Intcpt_deg_neg_Shr);
		ffmin_Shr = ffmin_effective;

		if (dmgSneg_Shr > 0.0) {

			delta_R_fcapneg_Shr = dmgSneg_Shr * fabs(fcapneg_Shr);
			ffmin_Shr = fmin((ffmin_effective + delta_R_fcapneg_Shr), R_fresneg_Shr);
			R_dcapneg_Shr = CR_dcapneg_Shr + fabs(delta_R_fcapneg_Shr) / (fabs(slope_neg_Shr) + fabs(R_Kdegneg_Shr));
			R_dresneg_Shr = CR_dresneg_Shr + fabs(delta_R_fcapneg_Shr) / fabs(R_Kdegneg_Shr); //Updating dresneg
			Intcpt_deg_neg_Shr = CIntcpt_deg_neg_Shr - fabs(delta_R_fcapneg_Shr);

		}
	}

	if (CstateFlag_Shr == -2 && TstateFlag_Shr == -4) {
		dpeakmin_Shr = Cstrain_Shr;
		ffmin_Shr = Cstress_Shr;
		//---------------------Changing this part to see if it is gonna correct the difference between the dpeakmin and depeakmax----------------------

		ratio_Shr = (dpeakmin_Shr - R_dcapneg_Shr) / (R_dresneg_Shr - R_dcapneg_Shr);
		dpeakmax_Shr = (R_drespos_Shr - R_dcappos_Shr) * ratio_Shr + R_dcappos_Shr;
		//---------------------------------------------------------------------------------------------------------------------------------------------
		ffmax_effective = Kdegpos_Shr * fabs(dpeakmax_Shr) + Intcpt_deg_pos_Shr;
		ffmax_Shr = ffmax_effective;

		if (dmgSpos_Shr > 0.0) {

			delta_R_fcappos_Shr = dmgSpos_Shr * fcappos_Shr;
			ffmax_Shr = fmax(ffmax_effective - delta_R_fcappos_Shr, frespos_Shr);
			R_drespos_Shr = CR_drespos_Shr - fabs(delta_R_fcappos_Shr) / fabs(Kdegpos_Shr);
			R_dcappos_Shr = CR_dcappos_Shr - fabs(delta_R_fcappos_Shr) / (fabs(slope_pos_Shr) + fabs(Kdegpos_Shr));
			Intcpt_deg_pos_Shr = CIntcpt_deg_pos_Shr - fabs(delta_R_fcappos_Shr);

		}
	}

	/* ************************************************************ **
	 * ************************************************************ **

					  Post Capping Damage Mode

	 * ************************************************************ **
	 * ************************************************************ **/

	if (flagdmg_Shr == 1) {

		if ((CstateFlag_Shr == 4 && TstateFlag_Shr == 6) || (CstateFlag_Shr == 5 && TstateFlag_Shr == 6) || (CstateFlag_Shr == -6 && TstateFlag_Shr == -7) || (CstateFlag_Shr == 7 && TstateFlag_Shr == -7)) {

			ffmax_effective = Kdegpos_Shr * fabs(dpeakmax_Shr) + Intcpt_deg_pos_Shr;
			ffmax_Shr = ffmax_effective;

			if (dmgSpos_Shr > 0.0) {

				delta_R_fcappos_Shr = dmgSpos_Shr * fcappos_Shr;
				ffmax_Shr = fmax(ffmax_effective - delta_R_fcappos_Shr, frespos_Shr);
				R_drespos_Shr = CR_drespos_Shr - fabs(delta_R_fcappos_Shr) / fabs(Kdegpos_Shr);
				R_dcappos_Shr = CR_dcappos_Shr - fabs(delta_R_fcappos_Shr) / (fabs(slope_pos_Shr) + fabs(Kdegpos_Shr));
				Intcpt_deg_pos_Shr = CIntcpt_deg_pos_Shr - fabs(delta_R_fcappos_Shr);

			}
		}
	}
	if (flagdmg_Shr == 1) {
		if ((CstateFlag_Shr == -4 && TstateFlag_Shr == -6) || (CstateFlag_Shr == -5 && TstateFlag_Shr == -6) || (CstateFlag_Shr == 6 && TstateFlag_Shr == 7) || (CstateFlag_Shr == -7 && TstateFlag_Shr == 7)) {

			ffmin_effective = -(R_Kdegneg_Shr * fabs(dpeakmin_Shr) + Intcpt_deg_neg_Shr);
			ffmin_Shr = ffmin_effective;

			if (dmgSneg_Shr > 0.0) {

				delta_R_fcapneg_Shr = dmgSneg_Shr * fabs(fcapneg_Shr);
				ffmin_Shr = fmin((ffmin_effective + delta_R_fcapneg_Shr), R_fresneg_Shr);
				R_dcapneg_Shr = CR_dcapneg_Shr + fabs(delta_R_fcapneg_Shr) / (fabs(slope_neg_Shr) + fabs(R_Kdegneg_Shr));
				R_dresneg_Shr = CR_dresneg_Shr + fabs(delta_R_fcapneg_Shr) / fabs(R_Kdegneg_Shr); //Updating dresneg
				Intcpt_deg_neg_Shr = CIntcpt_deg_neg_Shr - fabs(delta_R_fcapneg_Shr);
			}
		}
	}

	/* ************************************************************ **
	 * ************************************************************ **

					End of Post Capping Damage Mode

	 * ************************************************************ **
	 * ************************************************************ **/

	if (CstateFlag_Shr == 3 && TstateFlag_Shr == 4) {
		dpeakmax_Shr = Cstrain_Shr;
		ffmax_Shr = Cstress_Shr;
		ratio_Shr = (dpeakmax_Shr - R_drespos_Shr) / (delUpos_Shr - R_drespos_Shr);
		dpeakmin_Shr = (delUneg_Shr - R_dresneg_Shr) * ratio_Shr + R_dresneg_Shr;
		ffmin_Shr = fresneg_Shr;
	}

	if (CstateFlag_Shr == -3 && TstateFlag_Shr == -4) {
		dpeakmin_Shr = Cstrain_Shr;
		ffmin_Shr = Cstress_Shr;
		ratio_Shr = (dpeakmin_Shr - R_dresneg_Shr) / (delUneg_Shr - R_dresneg_Shr);
		dpeakmax_Shr = (delUpos_Shr - R_drespos_Shr) * ratio_Shr + R_drespos_Shr;
		ffmax_Shr = frespos_Shr;
	}

	if (CstateFlag_Shr == 30 && TstateFlag_Shr == 31) {
		dpeakmax_Shr = Cstrain_Shr;
		ffmax_Shr = Cstress_Shr;
		ratio_Shr = (dpeakmax_Shr - delUpos_Shr) / (Intcpt_Xaxis_pos_Shr - delUpos_Shr);
		dpeakmin_Shr = (Intcpt_Xaxis_neg_Shr - delUneg_Shr) * ratio_Shr + delUneg_Shr;
		ffmin_Shr = -(R_Kdegneg_Shr * fabs(dpeakmin_Shr) + Intcpt_res_neg_Shr);
	}

	if (CstateFlag_Shr == -30 && TstateFlag_Shr == 31) {
		dpeakmin_Shr = Cstrain_Shr;
		ffmin_Shr = Cstress_Shr;
		ratio_Shr = (dpeakmin_Shr - delUneg_Shr) / (Intcpt_Xaxis_neg_Shr - delUneg_Shr);
		dpeakmax_Shr = (Intcpt_Xaxis_pos_Shr - delUpos_Shr) * ratio_Shr + delUpos_Shr;
		ffmax_Shr = Kdegpos_Shr * fabs(dpeakmax_Shr) + Intcpt_res_pos_Shr;
	}
	//--------------------------------------End of Part 2-----------------------------------------------
}

// Spline parameter subroutine

void
GMG_CMAC2D::splineparam(double MtoRref, double dpeakpos, double dcappose, double dpeakneg, double dcapneg)
{

	n = 10.0;
	Lamda_KsecS = (ffmax - ffmin) / (dpeakmax - dpeakmin) / ((Kepos_Rot + Keneg_Rot) / 2.0);
	Lamda_KsecG = Lamda_KsecS * (1.0 + n) / (1.0 + n * Lamda_KsecS);
	if ((fabs(dpeakpos) < fabs(dcappose)) && (fabs(dpeakneg) < fabs(dcapneg))) {

		ER = fmin(alpha_Er_Rot_Post_Yielding * pow(Lamda_KsecG, beta_Er_Rot_Post_Yielding), Er_Rot_Post_Yielding);
		KrR = Kr_Rot_Post_Yielding;
		KuR = Kun_Rot_Post_Yielding;
	}
	else {

		ER = fmin(alpha_Er_Rot_Post_Capping * pow(Lamda_KsecG, beta_Er_Rot_Post_Capping), Er_Rot_Post_Capping);
		KrR = Kr_Rot_Post_Capping;
		KuR = Kun_Rot_Post_Capping;
		if (flag_FlexShr_Filure_Rot == 1) {

			ER = fmin(alpha_Er_Rot_Post_Capping_FlexShr * pow(Lamda_KsecG, beta_Er_Rot_Post_Capping_FlexShr), Er_Rot_Post_Capping_FlexShr);
			KrR = Kr_Rot_Post_Capping_FlexShr;
			KuR = Kun_Rot_Post_Capping_FlexShr;
		}
		if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {

			ER = fmin(alpha_Er_Rot_Post_Capping_splice * pow(Lamda_KsecG, beta_Er_Rot_Post_Capping_splice), Er_Rot_Post_Capping_splice);
			KrR = Kr_Rot_Post_Capping_splice;
			KuR = Kun_Rot_Post_Capping_splice;
		}

	}

}


void
GMG_CMAC2D::splineparam_Shr(double MtoRref, double dpeakpos, double dcappose, double dpeakneg, double dcapneg)
{

	double n = 10.0;
	double Lamda_KsecS_Shr = (ffmax_Shr - ffmin_Shr) / (dpeakmax_Shr - dpeakmin_Shr) / ((Kepos_Shr + Keneg_Shr) / 2.0);
	double Lamda_KsecG_Shr = Lamda_KsecS_Shr * (1.0 + n) / (1.0 + n * Lamda_KsecS_Shr);

	if ((fabs(dpeakpos) < fabs(dcappose)) && (fabs(dpeakneg) < fabs(dcapneg))) {

		ER_Shr = fmin(alpha_Er_Shr_Post_Yielding * pow(Lamda_KsecG_Shr, beta_Er_Shr_Post_Yielding), Er_Shr_Post_Yielding);
		KrR_Shr = Kr_Shr_Post_Yielding;
		KuR_Shr = Kun_Shr_Post_Yielding;
	}
	else {

		ER_Shr = fmin(alpha_Er_Shr_Post_Capping * pow(Lamda_KsecG_Shr, beta_Er_Shr_Post_Capping), Er_Shr_Post_Capping);
		KrR_Shr = Kr_Shr_Post_Capping;
		KuR_Shr = Kun_Shr_Post_Capping;
	}
}

/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Rotational Material                                    **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */
void
GMG_CMAC2D::spline_curve(double dirtag, double P0x, double P0y, double P4x, double P4y, double Krel, double Kun, double E, double x)
{
	//opserr << "Rasoolllll " << " dirtag " << dirtag << endln;
	double delx, C, P1x, P2x, P3x, P1y, P2y, P3y;
	double a, b, c, Q, R, theta, u;
	double B11, B21, B31, B41, B51;
	double B12, B22, B32, B42, B52;
	const double PI = 3.141592653589793238463;

	if (dirtag == 0) { // NP dir
		x = P4x - (x - P0x);

	}

	delx = (P4x - P0x) / 4;
	if (P4x * 100 > 4) {
		C = 0.1;
	}
	else {
		C = 0.1;
	}
	P1x = P0x + C / 4 * (P4x - P0x);
	P2x = P0x + 2 * delx;
	P3x = P4x - C / 4 * (P4x - P0x);

	P3y = P4y - Kun * (C*delx);
	P1y = P0y + Krel * (C*delx);
	P2y = 0.5 * C*delx*(2.0 * (Kun - Krel)) - 0.5 * E / delx + .5 * (P0y + P4y);

	// mapping x to u
	a = (P0x - P4x)*(9.0 / 2.0 * C - 3.0) / ((P0x - P4x)*(2.0 - 3.0 * C));
	b = 3.0 / 2.0 * C*(P4x - P0x) / ((P0x - P4x)*(2.0 - 3.0 * C));
	c = (P0x - x) / ((P0x - P4x)*(2.0 - 3.0 * C));
	Q = (pow(a, 2.0) - 3.0 * b) / 9.0;
	R = (2.0 * pow(a, 3.0) - 9.0 * a*b + 27.0 * c) / 54.0;
	theta = acos(R / sqrt(pow(Q, 3.0)));
	u = -2.0 * sqrt(Q)*cos((theta - 2.0 * PI) / 3.0) - a / 3.0;
	// solve basis function and its derivative
	if (u < 0.5) {
		B11 = -8.0 * pow((u - 0.5), 3.0);
		B21 = 14.0 * u*(0.428571 - 1.28571*u + pow(u, 2.0));
		B31 = -8.0 * pow(u, 2.0) * (u - 0.75);
		B41 = 2.0 * pow(u, 3.0);
		B51 = 0.0;
		B12 = -24.0 * pow((u - 0.5), 2.0);
		B22 = 6.0 * (1.0 - 6.0 * u + 7.0 * pow(u, 2.0));
		B32 = 12.0 * u*(1.0 - 2.0 * u);
		B42 = 6.0 * pow(u, 2.0);
		B52 = 0.0;
	}
	else {
		B11 = 0.0;
		B21 = -2.0 * pow((u - 1.0), 3.0);
		B31 = 8.0 * pow((1.0 - u), 2.0) * (u - 0.25);
		B41 = 2.0 - 12.0 * u + 24.0 * pow(u, 2.0) - 14.0 * pow(u, 3.0);
		B51 = 8.0 * pow((u - 0.5), 3.0);
		B12 = 0.0;
		B22 = -6.0 * pow((1.0 - u), 2.0);
		B32 = 12.0 * (1.0 - 3.0 * u + 2.0 * pow(u, 2.0));
		B42 = -12.0 + 48.0 * u - 42.0 * pow(u, 2.0);
		B52 = 24.0 * pow((u - 0.5), 2.0);
	}

	fspl = P0y * B11 + P1y * B21 + P2y * B31 + P3y * B41 + P4y * B51;
	fderspl = P0y * B12 + P1y * B22 + P2y * B32 + P3y * B42 + P4y * B52;
	dderspl = P0x * B12 + P1x * B22 + P2x * B32 + P3x * B42 + P4x * B52;
	ek = fderspl / dderspl;

	if (dirtag == 0) { // NP dir

		fspl = -fspl + P0y + P4y;
	}

}

/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Shear Material                                         **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */

void
GMG_CMAC2D::spline_curve_Shr(double dirtag, double P0x, double P0y, double P4x, double P4y, double Krel, double Kun, double E, double x)
{
	//opserr << "Rasoolllll " << " dirtag " << dirtag << endln;
	double delx, C, P1x, P2x, P3x, P1y, P2y, P3y;
	double a, b, c, Q, R, theta, u;
	double B11, B21, B31, B41, B51;
	double B12, B22, B32, B42, B52;
	const double PI = 3.141592653589793238463;

	if (dirtag == 0) { // NP dir
		x = P4x - (x - P0x);
	}

	delx = (P4x - P0x) / 4.0;
	if (P4x * 100 > 4) {
		C = 0.1;
	}
	else {
		C = 0.1;
	}
	P1x = P0x + C / 4.0 * (P4x - P0x);
	P2x = P0x + 2.0 * delx;
	P3x = P4x - C / 4.0 * (P4x - P0x);

	P3y = P4y - Kun * (C*delx);
	P1y = P0y + Krel * (C*delx);
	P2y = 0.5 * C*delx*(2.0 * (Kun - Krel)) - 0.5 * E / delx + 0.5 * (P0y + P4y);

	// mapping x to u
	a = (P0x - P4x)*(9.0 / 2.0 * C - 3.0) / ((P0x - P4x)*(2.0 - 3.0 * C));
	b = 3.0 / 2.0 * C*(P4x - P0x) / ((P0x - P4x)*(2.0 - 3.0 * C));
	c = (P0x - x) / ((P0x - P4x)*(2.0 - 3.0 * C));
	Q = (pow(a, 2.0) - 3.0 * b) / 9.0;
	R = (2.0 * pow(a, 3.0) - 9.0 * a*b + 27.0 * c) / 54.0;
	theta = acos(R / sqrt(pow(Q, 3.0)));
	u = -2.0 * sqrt(Q)*cos((theta - 2.0 * PI) / 3.0) - a / 3.0;
	// solve basis function and its derivative
	if (u < 0.5) {
		B11 = -8.0 * pow((u - 0.5), 3.0);
		B21 = 14.0 * u*(0.428571 - 1.28571*u + pow(u, 2.0));
		B31 = -8.0 * pow(u, 2.0) * (u - 0.75);
		B41 = 2.0 * pow(u, 3.0);
		B51 = 0.0;
		B12 = -24.0 * pow((u - 0.5), 2.0);
		B22 = 6.0 * (1.0 - 6.0 * u + 7.0 * pow(u, 2.0));
		B32 = 12.0 * u*(1.0 - 2.0 * u);
		B42 = 6.0 * pow(u, 2.0);
		B52 = 0.0;
	}
	else {
		B11 = 0.0;
		B21 = -2.0 * pow((u - 1.0), 3.0);
		B31 = 8.0 * pow((1.0 - u), 2.0) * (u - 0.25);
		B41 = 2.0 - 12.0 * u + 24.0 * pow(u, 2.0) - 14.0 * pow(u, 3.0);
		B51 = 8.0 * pow((u - 0.5), 3.0);
		B12 = 0.0;
		B22 = -6.0 * pow((1.0 - u), 2.0);
		B32 = 12.0 * (1.0 - 3.0 * u + 2.0 * pow(u, 2.0));
		B42 = -12.0 + 48.0 * u - 42.0 * pow(u, 2.0);
		B52 = 24.0 * pow((u - 0.5), 2.0);
	}

	fspl = P0y * B11 + P1y * B21 + P2y * B31 + P3y * B41 + P4y * B51;
	fderspl = P0y * B12 + P1y * B22 + P2y * B32 + P3y * B42 + P4y * B52;
	dderspl = P0x * B12 + P1x * B22 + P2x * B32 + P3x * B42 + P4x * B52;
	ek = fderspl / dderspl;

	if (dirtag == 0) { // NP dir

		fspl = -fspl + P0y + P4y;
	}

}


/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Rotational Material                                    **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */

//This section has been added to not go beyond the backbone during unloading or reloading
void
GMG_CMAC2D::checkEnvelope(void)
{

	if (F_Rot >= 0.0 && d_Rot >= 0.0) {

		if (d_Rot >= dpeakmax && d_Rot <= R_dcappos)
		{
			TstateFlag = 12;
			ek_Rot = slope_pos;

			F_Rot = ffmax;

		}

		if (d_Rot >= dpeakmax && d_Rot >= R_dcappos && d_Rot <= R_drespos)
		{
			TstateFlag = 2;
			ek_Rot = R_Kdegpos;

			F_Rot = ffmax;

		}

		if (d_Rot >= R_dcappos && d_Rot <= R_drespos && d_Rot < dpeakmax && F_Rot >= R_Kdegpos * fabs(d_Rot) + Intcpt_deg_pos) {
			TstateFlag = 2;
			ek_Rot = R_Kdegpos;
			F_Rot = R_Kdegpos * fabs(d_Rot) + Intcpt_deg_pos;
			dpeakmax = d_Rot;
		}

		if (d_Rot >= dpeakmax && d_Rot >= R_drespos && d_Rot < R_delUpos)
		{
			TstateFlag = 3;
			ek_Rot = 0.0001;

			F_Rot = R_frespos;
		}

		if (d_Rot >= dpeakmax && d_Rot >= R_delUpos && d_Rot <= Intcpt_Xaxis_pos)
		{
			TstateFlag = 30;
			ek_Rot = R_Kdegpos;
			F_Rot = ffmax;
		}
	}

	else if (F_Rot < 0.0 && d_Rot < 0.0)
	{
		if (d_Rot <= dpeakmin && d_Rot >= R_dcapneg)
		{
			TstateFlag = -12;
			ek_Rot = slope_neg;
			F_Rot = ffmin;
		}
		if (d_Rot <= dpeakmin && d_Rot <= R_dcapneg && d_Rot >= R_dresneg)
		{
			TstateFlag = -2;
			ek_Rot = R_Kdegneg;

			F_Rot = -(R_Kdegneg * fabs(d_Rot) + Intcpt_deg_neg);

		}

		// For Hardening
		if (d_Rot > dpeakmin && d_Rot <= R_dcapneg && d_Rot >= R_dresneg && F_Rot <= -(R_Kdegneg * fabs(d_Rot) + Intcpt_deg_neg))
		{
			TstateFlag = -2;
			ek_Rot = R_Kdegneg;
			F_Rot = -(R_Kdegneg * fabs(d_Rot) + Intcpt_deg_neg);
			dpeakmin = d_Rot;
		}

		if (d_Rot <= dpeakmin && d_Rot <= R_dresneg && d_Rot > R_delUneg)
		{
			TstateFlag = -3;
			ek_Rot = 0.0001;
			F_Rot = R_fresneg;
		}

		if (d_Rot <= dpeakmin && d_Rot <= R_delUneg && d_Rot >= Intcpt_Xaxis_neg)
		{
			TstateFlag = -30;
			ek_Rot = R_Kdegneg;
			F_Rot = ffmin;

		}
	}
}
//************************************************************************************************************


/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Shear Material                                         **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */


void
GMG_CMAC2D::checkEnvelope_Shr(void)
{

	if (F_Shr >= 0.0 && d_Shr >= 0.0) {

		if (d_Shr >= dpeakmax_Shr && d_Shr <= R_dcappos_Shr)
		{
			TstateFlag_Shr = 12;
			ek_Shr = slope_pos_Shr;

			F_Shr = ffmax_Shr;
		}

		if (d_Shr >= dpeakmax_Shr && d_Shr >= R_dcappos_Shr && d_Shr <= R_drespos_Shr)
		{
			TstateFlag_Shr = 2;
			ek_Shr = Kdegpos_Shr;

			F_Shr = ffmax_Shr;
		}

		//For hardening cases
		if (d_Shr >= R_dcappos_Shr && d_Shr <= R_drespos_Shr && d_Shr < dpeakmax_Shr && F_Shr >= Kdegpos_Shr * fabs(d_Shr) + Intcpt_deg_pos_Shr) {
			TstateFlag_Shr = 2;
			ek_Shr = Kdegpos_Shr;
			F_Shr = Kdegpos_Shr * fabs(d_Shr) + Intcpt_deg_pos_Shr;
			dpeakmax_Shr = d_Shr;
		}

		if (d_Shr >= dpeakmax_Shr && d_Shr >= R_drespos_Shr && d_Shr < delUpos_Shr)
		{
			TstateFlag_Shr = 3;
			ek_Shr = 0.0001;
			F_Shr = frespos_Shr;
		}
		if (d_Shr >= dpeakmax_Shr && d_Shr >= delUpos_Shr && d_Shr <= Intcpt_Xaxis_pos_Shr)
		{
			TstateFlag_Shr = 30;
			ek_Shr = Kdegpos_Shr;

			F_Shr = Kdegpos_Shr * fabs(d_Shr) + Intcpt_res_pos_Shr;
		}
	}

	else if (F_Shr < 0.0 && d_Shr < 0.0)
	{

		if (d_Shr <= dpeakmin_Shr && d_Shr >= R_dcapneg_Shr)
		{
			TstateFlag_Shr = -12;
			ek_Shr = slope_neg_Shr;
			F_Shr = ffmin_Shr;
		}
		if (d_Shr < dpeakmin_Shr && d_Shr <= R_dcapneg_Shr && d_Shr >= R_dresneg_Shr)
		{
			TstateFlag_Shr = -2;
			ek_Shr = R_Kdegneg_Shr;
			F_Shr = ffmin_Shr;
		}

		// For Hardening
		if (d_Shr > dpeakmin_Shr && d_Shr <= R_dcapneg_Shr && d_Shr >= R_dresneg_Shr && F_Shr <= -(R_Kdegneg_Shr * fabs(d_Shr) + Intcpt_deg_neg_Shr))
		{
			TstateFlag_Shr = -2;
			ek_Shr = R_Kdegneg_Shr;
			F_Shr = -(R_Kdegneg_Shr * fabs(d_Shr) + Intcpt_deg_neg_Shr);
			dpeakmin_Shr = d_Shr;
		}

		if (d_Shr <= dpeakmin_Shr && d_Shr <= R_dresneg_Shr && d_Shr > delUneg_Shr)
		{
			TstateFlag_Shr = -3;
			ek_Shr = 0.0001;
			F_Shr = fresneg_Shr;
		}

		if (d_Shr <= dpeakmin_Shr && d_Shr <= delUneg_Shr && d_Shr >= Intcpt_Xaxis_neg_Shr)
		{
			TstateFlag_Shr = -30;
			ek_Shr = R_Kdegneg_Shr;
			F_Shr = ffmin_Shr;

		}
	}
}

/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Rotational Material                                    **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */

//This section has been added to update the damage after entering case 2 or -2 for the first time
void
GMG_CMAC2D::update_damage(void)
{
	//opserr << " update_damage" << "  " << " Hi " << endln;

	double Energy_Coe_Rot = 0.0;

	spos = (F_Rot + fabs(F_Rot)) / 2;
	sneg = (F_Rot - fabs(F_Rot)) / 2;

	/* ********************************************************************************************** **
	************************************************************************************************* **
	**                                                                                                **
	**                             Check which mode of failure is active                              **
	**                                                                                                **
	************************************************************************************************* **
	** ********************************************************************************************** */

	// Dealing with coefficient
	if (d_Rot >= dpeakmax || d_Rot <= dpeakmin) {
		Energy_Coe_Rot = C3_Rot;
		if (flag_FlexShr_Filure_Rot == 1) {
			Energy_Coe_Rot = C3_FlexShr_Rot;
		}
		if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
			Energy_Coe_Rot = C3_splice_Rot;
		}
	}
	else if (Tdu_Rot <= 0 && F_Rot >= 0) {
		Energy_Coe_Rot = C1_Rot;
		if (flag_FlexShr_Filure_Rot == 1) {
			Energy_Coe_Rot = C1_FlexShr_Rot;
		}
		if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
			Energy_Coe_Rot = C1_splice_Rot;
		}
	}
	else if (Tdu_Rot <= 0 && F_Rot < 0) {
		Energy_Coe_Rot = C2_Rot;
		if (flag_FlexShr_Filure_Rot == 1) {
			Energy_Coe_Rot = C2_FlexShr_Rot;
		}
		if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
			Energy_Coe_Rot = C2_splice_Rot;
		}
	}
	else if (Tdu_Rot >= 0 && F_Rot >= 0) {
		Energy_Coe_Rot = C2_Rot;
		if (flag_FlexShr_Filure_Rot == 1) {
			Energy_Coe_Rot = C2_FlexShr_Rot;
		}
		if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
			Energy_Coe_Rot = C2_splice_Rot;
		}
	}
	else if (Tdu_Rot >= 0 && F_Rot < 0) {
		Energy_Coe_Rot = C1_Rot;
		if (flag_FlexShr_Filure_Rot == 1) {
			Energy_Coe_Rot = C1_FlexShr_Rot;
		}
		if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
			Energy_Coe_Rot = C1_splice_Rot;
		}
	}
	Et = CEt + Energy_Coe_Rot * fabs(Cspos + spos) * fabs(Tdu_Rot) / 2;
	Ec = CEc + Energy_Coe_Rot * fabs(Csneg + sneg) * fabs(Tdu_Rot) / 2;
	Epos = Et + Ec;
	Eneg = Et + Ec;

	if (Tdu_Rot < 0) {

		if ((CstateFlag == -5 && TstateFlag == -6) || (CstateFlag == 6 && TstateFlag == 7) || (CstateFlag == -7 && TstateFlag == 7) || (CstateFlag == 2 && TstateFlag == 4) || (CstateFlag == 12 && TstateFlag == 4)) {

			spos_neg = Cspos_neg = sneg_neg = Csneg_neg = Et_neg = CEt_neg = Ec_neg = CEc_neg = Enegnorm = dmgSneg = 0.0;

		}

		if (Epos >= Ed0pos) {

			spos_pos = (F_Rot + fabs(F_Rot)) / 2;
			sneg_pos = (F_Rot - fabs(F_Rot)) / 2;
			if (d_Rot >= dpeakmax || d_Rot <= dpeakmin) {
				Energy_Coe_Rot = C3_Rot;
				if (flag_FlexShr_Filure_Rot == 1) {
					Energy_Coe_Rot = C3_FlexShr_Rot;
				}
				if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
					Energy_Coe_Rot = C3_splice_Rot;
				}
			}
			else if (Tdu_Rot <= 0 && F_Rot >= 0) {
				Energy_Coe_Rot = C1_Rot;
				if (flag_FlexShr_Filure_Rot == 1) {
					Energy_Coe_Rot = C1_FlexShr_Rot;
				}
				if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
					Energy_Coe_Rot = C1_splice_Rot;
				}
			}
			else if (Tdu_Rot <= 0 && F_Rot < 0) {
				Energy_Coe_Rot = C2_Rot;
				if (flag_FlexShr_Filure_Rot == 1) {
					Energy_Coe_Rot = C2_FlexShr_Rot;
				}
				if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
					Energy_Coe_Rot = C2_splice_Rot;
				}
			}

			Et_pos = CEt_pos + Energy_Coe_Rot * fabs(Cspos_pos + spos_pos) * fabs(Tdu_Rot) / 2;
			Ec_pos = CEc_pos + Energy_Coe_Rot * fabs(Csneg_pos + sneg_pos) * fabs(Tdu_Rot) / 2;
			Epos_pos = Et_pos + Ec_pos;
			Eposnorm = fmax(0.0, (Epos_pos) / (Ed1pos - Ed0pos));

			if (flagdmg == 1) {
				dmgSpos = fmin(solpe_post_capping_Rot * Eposnorm, 1); //cdf(d, Eposnorm);
				if (flag_FlexShr_Filure_Rot == 1) {
					dmgSpos = fmin(solpe_post_capping_FlexShr_Rot * Eposnorm, 1); //cdf(d, Eposnorm);
				}
				if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
					dmgSpos = fmin(solpe_post_capping_splice_Rot * Eposnorm, 1); //cdf(d, Eposnorm);
				}
			}

			if (flagdmg_Hardening_strength == 1) {
				dmgSpos = fmin(solpe_post_yielding_Rot * Eposnorm, 1); //cdf(d, Eposnorm);
			}

		}
	}

	if (Tdu_Rot > 0) {

		if ((CstateFlag == 5 && TstateFlag == 6) || (CstateFlag == -6 && TstateFlag == -7) || (CstateFlag == 7 && TstateFlag == -7) || (CstateFlag == -2 && TstateFlag == -4) || (CstateFlag == -12 && TstateFlag == -4)) {

			spos_pos = Cspos_pos = sneg_pos = Csneg_pos = Et_pos = CEt_pos = Ec_pos = CEc_pos = Eposnorm = dmgSpos = 0.0;

		}

		if (Epos >= Ed0pos) {

			spos_neg = (F_Rot + fabs(F_Rot)) / 2;
			sneg_neg = (F_Rot - fabs(F_Rot)) / 2;
			if (d_Rot >= dpeakmax || d_Rot <= dpeakmin) {
				Energy_Coe_Rot = C3_Rot;
				if (flag_FlexShr_Filure_Rot == 1) {
					Energy_Coe_Rot = C3_FlexShr_Rot;
				}
				if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
					Energy_Coe_Rot = C3_splice_Rot;
				}
			}
			else if (Tdu_Rot >= 0 && F_Rot < 0) {
				Energy_Coe_Rot = C1_Rot;
				if (flag_FlexShr_Filure_Rot == 1) {
					Energy_Coe_Rot = C1_FlexShr_Rot;
				}
				if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
					Energy_Coe_Rot = C1_splice_Rot;
				}
			}
			else if (Tdu_Rot >= 0 && F_Rot >= 0) {
				Energy_Coe_Rot = C2_Rot;
				if (flag_FlexShr_Filure_Rot == 1) {
					Energy_Coe_Rot = C2_FlexShr_Rot;
				}
				if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
					Energy_Coe_Rot = C2_splice_Rot;
				}
			}

			Et_neg = CEt_neg + Energy_Coe_Rot * fabs(Cspos_neg + spos_neg) * fabs(Tdu_Rot) / 2;
			Ec_neg = CEc_neg + Energy_Coe_Rot * fabs(Csneg_neg + sneg_neg) * fabs(Tdu_Rot) / 2;
			Eneg_neg = Et_neg + Ec_neg;
			Enegnorm = fmax(0.0, (Eneg_neg) / (Ed1neg - Ed0neg));

			if (flagdmg == 1) {
				dmgSneg = fmin(solpe_post_capping_Rot * Enegnorm, 1); //cdf(d, Enegnorm);
				if (flag_FlexShr_Filure_Rot == 1) {
					dmgSneg = fmin(solpe_post_capping_FlexShr_Rot * Enegnorm, 1); //cdf(d, Enegnorm);
				}
				if (flag_Splice_Filure_Rot == 1 || flag_FlexSplice_Filure_Rot == 1) {
					dmgSneg = fmin(solpe_post_capping_splice_Rot * Enegnorm, 1); //cdf(d, Enegnorm);
				}
			}

			if (flagdmg_Hardening_strength == 1) {
				dmgSneg = fmin(solpe_post_yielding_Rot * Enegnorm, 1); //cdf(d, Enegnorm);
			}

		}
	}
}
//************************************************************************************************************


/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Shear Material                                         **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */
void
GMG_CMAC2D::update_damage_Shr(void)
{
	double Energy_Coe_Shr;
	spos_Shr = (F_Shr + fabs(F_Shr)) / 2;
	sneg_Shr = (F_Shr - fabs(F_Shr)) / 2;
	if (d_Shr >= dpeakmax_Shr || d_Shr <= dpeakmin_Shr) {
		Energy_Coe_Shr = C3_Shr;
	}
	else if (Tdu_Shr <= 0 && F_Shr >= 0) {
		Energy_Coe_Shr = C1_Shr;
	}
	else if (Tdu_Shr <= 0 && F_Shr < 0) {
		Energy_Coe_Shr = C2_Shr;
	}
	else if (Tdu_Shr >= 0 && F_Shr >= 0) {
		Energy_Coe_Shr = C2_Shr;
	}
	else if (Tdu_Shr >= 0 && F_Shr < 0) {
		Energy_Coe_Shr = C1_Shr;
	}
	Et_Shr = CEt_Shr + Energy_Coe_Shr * fabs(Cspos_Shr + spos_Shr) * fabs(Tdu_Shr) / 2;
	Ec_Shr = CEt_Shr + Energy_Coe_Shr * fabs(Csneg_Shr + sneg_Shr) * fabs(Tdu_Shr) / 2;
	Epos_Shr = Et_Shr + Ec_Shr;
	Eneg_Shr = Et_Shr + Ec_Shr;

	if (Tdu_Shr < 0) {

		if ((CstateFlag_Shr == -5 && TstateFlag_Shr == -6) || (CstateFlag_Shr == 6 && TstateFlag_Shr == 7) || (CstateFlag_Shr == -7 && TstateFlag_Shr == 7) || (CstateFlag_Shr == 2 && TstateFlag_Shr == 4) || (CstateFlag_Shr == 12 && TstateFlag_Shr == 4)) {

			spos_neg_Shr = Cspos_neg_Shr = sneg_neg_Shr = Csneg_neg_Shr = Et_neg_Shr = CEt_neg_Shr = Ec_neg_Shr = CEc_neg_Shr = Eneg_neg_Shr = Enegnorm_Shr = dmgSneg_Shr = 0.0;

		}

		if (Epos_Shr >= Ed0pos_Shr) {

			spos_pos_Shr = (F_Shr + fabs(F_Shr)) / 2;
			sneg_pos_Shr = (F_Shr - fabs(F_Shr)) / 2;
			if (d_Shr >= dpeakmax_Shr || d_Shr <= dpeakmin_Shr) {
				Energy_Coe_Shr = C3_Shr;
			}
			else if (Tdu_Shr <= 0 && F_Shr >= 0) {
				Energy_Coe_Shr = C1_Shr;
			}
			else if (Tdu_Shr <= 0 && F_Shr < 0) {
				Energy_Coe_Shr = C2_Shr;
			}
			Et_pos_Shr = CEt_pos_Shr + Energy_Coe_Shr * fabs(Cspos_pos_Shr + spos_pos_Shr) * fabs(Tdu_Shr) / 2;
			Ec_pos_Shr = CEc_pos_Shr + Energy_Coe_Shr * fabs(Csneg_pos_Shr + sneg_pos_Shr) * fabs(Tdu_Shr) / 2;
			Epos_pos_Shr = Et_pos_Shr + Ec_pos_Shr;
			Eposnorm_Shr = fmax(0.0, (Epos_pos_Shr) / (Ed1pos_Shr - Ed0pos_Shr));

			if (flagdmg_Shr == 1) {
				dmgSpos_Shr = fmin(solpe_post_capping_Shr * Eposnorm_Shr, 1); //cdf(d, Enegnorm);
			}

			if (flagdmg_Hardening_strength_Shr == 1) {
				dmgSpos_Shr = fmin(solpe_post_yielding_Shr * Eposnorm_Shr, 1); //cdf(d, Enegnorm);
			}
		}
	}

	if (Tdu_Shr > 0) {

		if ((CstateFlag_Shr == 5 && TstateFlag_Shr == 6) || (CstateFlag_Shr == -6 && TstateFlag_Shr == -7) || (CstateFlag_Shr == 7 && TstateFlag_Shr == -7) || (CstateFlag_Shr == -2 && TstateFlag_Shr == -4) || (CstateFlag_Shr == -12 && TstateFlag_Shr == -4)) {

			spos_pos_Shr = Cspos_pos_Shr = sneg_pos_Shr = Csneg_pos_Shr = Et_pos_Shr = CEt_pos_Shr = Ec_pos_Shr = CEc_pos_Shr = Epos_pos_Shr = Eposnorm_Shr = dmgSpos_Shr = 0.0;

		}

		if (Eneg_Shr >= Ed0neg_Shr) {

			spos_neg_Shr = (F_Shr + fabs(F_Shr)) / 2;
			sneg_neg_Shr = (F_Shr - fabs(F_Shr)) / 2;
			if (d_Shr >= dpeakmax_Shr || d_Shr <= dpeakmin_Shr) {
				Energy_Coe_Shr = C3_Shr;
			}
			else if (Tdu_Shr >= 0 && F_Shr >= 0) {
				Energy_Coe_Shr = C2_Shr;
			}
			else if (Tdu_Shr >= 0 && F_Shr < 0) {
				Energy_Coe_Shr = C1_Shr;
			}
			Et_neg_Shr = CEt_neg_Shr + Energy_Coe_Shr * fabs(Cspos_neg_Shr + spos_neg_Shr) * fabs(Tdu_Shr) / 2;
			Ec_neg_Shr = CEc_neg_Shr + Energy_Coe_Shr * fabs(Csneg_neg_Shr + sneg_neg_Shr) * fabs(Tdu_Shr) / 2;
			Eneg_neg_Shr = Et_neg_Shr + Ec_neg_Shr;
			Enegnorm_Shr = fmax(0.0, (Eneg_neg_Shr) / (Ed1neg_Shr - Ed0neg_Shr));

			if (flagdmg_Shr == 1) {
				dmgSneg_Shr = fmin(solpe_post_capping_Shr * Enegnorm_Shr, 1); //cdf(d, Enegnorm);
			}

			if (flagdmg_Hardening_strength_Shr == 1) {
				dmgSneg_Shr = fmin(solpe_post_yielding_Shr * Enegnorm_Shr, 1); //cdf(d, Enegnorm);
			}
		}
	}
}

/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Rotational Material                                    **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */


void
GMG_CMAC2D::update_damage_hardeingin(void)
{
	double Energy_Coe_Stiff_Rot;
	spos_hard = (F_Rot + fabs(F_Rot)) / 2;
	sneg_hard = (F_Rot - fabs(F_Rot)) / 2;
	Energy_Coe_Stiff_Rot = 1.0;

	Et_hard = CEt_hard + Energy_Coe_Stiff_Rot * fabs(Cspos_hard + spos_hard) * fabs(Tdu_Rot) / 2;
	Ec_hard = CEc_hard + Energy_Coe_Stiff_Rot * fabs(Csneg_hard + sneg_hard) * fabs(Tdu_Rot) / 2;

	Epos_hard = Et_hard + Ec_hard;
	Eneg_hard = Et_hard + Ec_hard;

	if (Tdu_Rot < 0 && dpeakmax > R_dypos) {

		if ((CstateFlag == -5 && TstateFlag == -6) || (CstateFlag == 6 && TstateFlag == 7) || (CstateFlag == -7 && TstateFlag == 7) || (CstateFlag == 12 && TstateFlag == 4)) {

			spos_neg_hard = Cspos_neg_hard = sneg_neg_hard = Csneg_neg_hard = Et_neg_hard = CEt_neg_hard = Ec_neg_hard = CEc_neg_hard = 0.0;
			Et_neg_hard = CEt_neg_hard = Ec_neg_hard = CEc_neg_hard = Eneg_neg_hard = Enegnorm_hard = 0.0;
		}

		if (Epos_hard >= Ed0pos_hard) {

			spos_pos_hard = (F_Rot + fabs(F_Rot)) / 2;
			sneg_pos_hard = (F_Rot - fabs(F_Rot)) / 2;
			Energy_Coe_Stiff_Rot = 1.0;

			Et_pos_hard = CEt_pos_hard + Energy_Coe_Stiff_Rot * fabs(Cspos_pos_hard + spos_pos_hard) * fabs(Tdu_Rot) / 2;
			Ec_pos_hard = CEc_pos_hard + Energy_Coe_Stiff_Rot * fabs(Csneg_pos_hard + sneg_pos_hard) * fabs(Tdu_Rot) / 2;
			Epos_pos_hard = Et_pos_hard + Ec_pos_hard;
			Eposnorm_hard = fmax(0.0, (Epos_pos_hard) / (Ed1pos_hard - Ed0pos_hard));
			delta_pos_hard = solpe_post_capping_Rot * Eposnorm_hard /* + Cdelta_pos_hard*/;

			alpha_pos = delta_ratio_max_hard_Rot * Eposnorm_hard;

		}
	}

	if (Tdu_Rot > 0 && dpeakmin < R_dyneg) {

		if ((CstateFlag == 5 && TstateFlag == 6) || (CstateFlag == -6 && TstateFlag == -7) || (CstateFlag == 7 && TstateFlag == -7) || (CstateFlag == -12 && TstateFlag == -4)) {

			spos_pos_hard = Cspos_pos_hard = sneg_pos_hard = Csneg_pos_hard = Et_pos_hard = CEt_pos_hard = Ec_pos_hard = CEc_pos_hard = 0.0;
			Et_pos_hard = CEt_pos_hard = Ec_pos_hard = CEc_pos_hard = Epos_pos_hard = Eposnorm_hard = 0.0;
		}


		if (Epos_hard >= Ed0pos_hard) {

			spos_neg_hard = (F_Rot + fabs(F_Rot)) / 2;
			sneg_neg_hard = (F_Rot - fabs(F_Rot)) / 2;
			Energy_Coe_Stiff_Rot = 1.0;

			Et_neg_hard = CEt_neg_hard + Energy_Coe_Stiff_Rot * fabs(Cspos_neg_hard + spos_neg_hard) * fabs(Tdu_Rot) / 2;
			Ec_neg_hard = CEc_neg_hard + Energy_Coe_Stiff_Rot * fabs(Csneg_neg_hard + sneg_neg_hard) * fabs(Tdu_Rot) / 2;
			Eneg_neg_hard = Et_neg_hard + Ec_neg_hard;
			Enegnorm_hard = fmax(0.0, (Eneg_neg_hard) / (Ed1neg_hard - Ed0neg_hard));
			delta_neg_hard = solpe_post_capping_Rot * Enegnorm_hard /* + Cdelta_neg_hard */;

			alpha_neg = delta_ratio_max_hard_Rot * Enegnorm_hard;

		}
	}
}

/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                         Shear Material                                         **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */

void
GMG_CMAC2D::update_damage_hardeingin_Shr(void)
{
	double Energy_Coe_Stiff_Shr;
	spos_hard_Shr = (F_Shr + fabs(F_Shr)) / 2;
	sneg_hard_Shr = (F_Shr - fabs(F_Shr)) / 2;
	Energy_Coe_Stiff_Shr = 1.0;

	Et_hard_Shr = CEt_hard_Shr + Energy_Coe_Stiff_Shr * fabs(Cspos_hard_Shr + spos_hard_Shr) * fabs(Tdu_Shr) / 2;
	Ec_hard_Shr = CEc_hard_Shr + Energy_Coe_Stiff_Shr * fabs(Csneg_hard_Shr + sneg_hard_Shr) * fabs(Tdu_Shr) / 2;

	Epos_hard_Shr = Et_hard_Shr + Ec_hard_Shr;
	Eneg_hard_Shr = Et_hard_Shr + Ec_hard_Shr;

	if (Tdu_Shr < 0 && dpeakmax_Shr > dypos_Shr) {

		if ((CstateFlag_Shr == -5 && TstateFlag_Shr == -6) || (CstateFlag_Shr == 6 && TstateFlag_Shr == 7) || (CstateFlag_Shr == -7 && TstateFlag_Shr == 7) || (CstateFlag_Shr == 12 && TstateFlag_Shr == 4)) {

			spos_neg_hard_Shr = Cspos_neg_hard_Shr = sneg_neg_hard_Shr = Csneg_neg_hard_Shr = Et_neg_hard_Shr = CEt_neg_hard_Shr = Ec_neg_hard_Shr = CEc_neg_hard_Shr = 0.0;
			Et_neg_hard_Shr = CEt_neg_hard_Shr = Ec_neg_hard_Shr = CEc_neg_hard_Shr = Eneg_neg_hard_Shr = Enegnorm_hard_Shr = 0.0;
		}

		if (Epos_hard_Shr >= Ed0pos_hard_Shr) {

			spos_pos_hard_Shr = (F_Shr + fabs(F_Shr)) / 2;
			sneg_pos_hard_Shr = (F_Shr - fabs(F_Shr)) / 2;
			Energy_Coe_Stiff_Shr = 1.0;

			Et_pos_hard_Shr = CEt_pos_hard_Shr + Energy_Coe_Stiff_Shr * fabs(Cspos_pos_hard_Shr + spos_pos_hard_Shr) * fabs(Tdu_Shr) / 2;
			Ec_pos_hard_Shr = CEc_pos_hard_Shr + Energy_Coe_Stiff_Shr * fabs(Csneg_pos_hard_Shr + sneg_pos_hard_Shr) * fabs(Tdu_Shr) / 2;
			Epos_pos_hard_Shr = Et_pos_hard_Shr + Ec_pos_hard_Shr;
			Eposnorm_hard_Shr = fmax(0.0, (Epos_pos_hard_Shr) / (Ed1pos_hard_Shr - Ed0pos_hard_Shr));
			delta_pos_hard_Shr = solpe_post_capping_Shr * Eposnorm_hard_Shr /* + Cdelta_pos_hard*/;

			alpha_pos_Shr = delta_ratio_max_hard_Shr * Eposnorm_hard_Shr;

		}
	}

	if (Tdu_Shr > 0 && dpeakmin_Shr < dyneg_Shr) {

		if ((CstateFlag_Shr == 5 && TstateFlag_Shr == 6) || (CstateFlag_Shr == -6 && TstateFlag_Shr == -7) || (CstateFlag_Shr == 7 && TstateFlag_Shr == -7) || (CstateFlag_Shr == -12 && TstateFlag_Shr == -4)) {

			spos_pos_hard_Shr = Cspos_pos_hard_Shr = sneg_pos_hard_Shr = Csneg_pos_hard_Shr = Et_pos_hard_Shr = CEt_pos_hard_Shr = Ec_pos_hard_Shr = CEc_pos_hard_Shr = 0.0;
			Et_pos_hard_Shr = CEt_pos_hard_Shr = Ec_pos_hard_Shr = CEc_pos_hard_Shr = Epos_pos_hard_Shr = Eposnorm_hard_Shr = 0.0;
		}

		if (Epos_hard_Shr >= Ed0pos_hard_Shr) {

			spos_neg_hard_Shr = (F_Shr + fabs(F_Shr)) / 2;
			sneg_neg_hard_Shr = (F_Shr - fabs(F_Shr)) / 2;
			Energy_Coe_Stiff_Shr = 1.0;

			Et_neg_hard_Shr = CEt_neg_hard_Shr + Energy_Coe_Stiff_Shr * fabs(Cspos_neg_hard_Shr + spos_neg_hard_Shr) * fabs(Tdu_Shr) / 2;
			Ec_neg_hard_Shr = CEc_neg_hard_Shr + Energy_Coe_Stiff_Shr * fabs(Csneg_neg_hard_Shr + sneg_neg_hard_Shr) * fabs(Tdu_Shr) / 2;
			Eneg_neg_hard_Shr = Et_neg_hard_Shr + Ec_neg_hard_Shr;
			Enegnorm_hard_Shr = fmax(0.0, (Eneg_neg_hard_Shr) / (Ed1neg_hard_Shr - Ed0neg_hard_Shr));
			delta_neg_hard_Shr = solpe_post_capping_Shr * Enegnorm_hard_Shr /* + Cdelta_neg_hard */;

			alpha_neg_Shr = delta_ratio_max_hard_Shr * Enegnorm_hard_Shr;

		}
	}
}

void
GMG_CMAC2D::P_M_yield(double Axial_Force) {

	//opserr << "P_M_initials " << "  " << " H " << "  " << H << "  " << " B " << "  " << B << "  " << " C_C " << "  " << C_C << "  " << " L " << "  " << L << "  " << " La " << "  " << La << "  " << " fc " << "  " << fc << "  " << " fyL " << "  " << fyL << "  " << " dbL " << "  " << dbL << "  " << " nL_EW " << "  " << nL_EW << "  " << " nL_NS " << "  " << nL_NS << "  " << " fyT " << "  " << fyT << "  " << " dbT " << "  " << dbT << "  " << " S " << "  " << S << "  " << " n_leg " << "  " << n_leg << "  " << " lb " << "  " << lb << "  " << " P_Axial " << "  " << P_Axial << endln;
	// void P_M(double H, double B, double C_C, double fc, double fyL, double dbL, double dbT, double nL_EW, double nL_NS);
	/* ********************************************************************************************** **
	************************************************************************************************* **
	************************************************************************************************* **
	**                                                                                                **
	**                                    Moment Interaction Diagram                                  **
	**                                                                                                **
	************************************************************************************************* **
	************************************************************************************************* **
	** ********************************************************************************************** */

	// Using ACI 318 Assumptions, Rasool Ghorbani
		//           B
		//    < ----------->
		//  ^ | ---------- |
		//  | |  @   @   @ | --->First layer, d1 = deff, n_barT bars
		//  | |            | .
		//  | |            | .
		// H| |  @       @ | . --->n_layM layers in the middle, n_barM bars per each
		//  | |            | .
		//  | |            | .
		//  | |  @   @   @ | --->last layer, d(end), n_barT bars
		//  ^ | ---------- |

		//function[TrueMoment, Force_section, Moment_Section] = PM_INTERACTION_DIAGRAM(H, B, fc, fyL, deff, n_layM, n_barT, n_barM, dbL, P)
		//------------------------ - Concrete properties----------------------------

	const double PI = 3.141592653589793238463;
	double AsL = PI * pow(dbL, 2.0) / 4.0;
	double AsL_Total = (2.0 * nL_EW*AsL + 2.0 * (nL_NS - 2.0)*AsL);

	double deff = 0.8*H;
	Matrix d(1, nL_NS);
	for (int i = 0; i <= nL_NS - 1; i++) {
		if (i == 0)
			d(0, i) = H - (H - C_C - dbT - dbL / 2);
		else if (i > 0) {
			d(0, i) = (i)*(H - 2.0 * d(0, 0)) / (nL_NS - 1) + d(0, 0);
		}
	}

	double f_prime_c = fc * 1000.0; //psi
	double b = B;
	double h = H;
	double Ag = b * h;

	// ------------------ - Steel Properties(elasto plastic)------------
	double Es = 29000.0;
	double fy = fyL;
	double e_y = fy / Es;

	Matrix As(1, nL_NS);
	for (int i = 0; i <= nL_NS - 1; i++) {
		if (i == 0)
			As(0, i) = nL_EW * PI *pow(dbL, 2) / 4;
		else if (i >= 1 && i < nL_NS - 1)
			As(0, i) = 2.0 * PI *pow(dbL, 2) / 4;
		else if (i == nL_NS - 1)
			As(0, i) = nL_EW * PI *pow(dbL, 2) / 4;
	}
	double As_total = 0.0;
	for (int i = 0; i <= nL_NS - 1; i++) {
		As_total = As_total + As(0, i);
	}
	double Plastic_slope_Coeficient = 0.0;
	double E1 = Plastic_slope_Coeficient * Es; // Hardening slope

		// ------------------------Interaction Diagram------------------------ -
	double counter = 1;
	for (double z = 0.5; z >= -5.0; z = z - 0.1) {
		counter += 1;
	}

	Matrix Force_section(1, counter + 2);
	Matrix Moment_Section(1, counter + 2);
	Matrix Moment_concrete_total_Comp(1, counter + 1);
	Matrix Moment_Steel_total_layer(nL_NS, counter + 1);
	Matrix Fc_total_concrete_comp(1, counter + 1);
	Matrix Fs_total_steel_layer(nL_NS, counter + 1);
	Matrix kd(1, counter + 1);

	//---------------------- ACI 318-19 --------------------------
	double alphaa = 0.85;
	double betaa = 0;
	if (f_prime_c <= 4000.0)
		betaa = 0.85;
	else if (f_prime_c > 4000.0 && f_prime_c <= 8000.0)
		betaa = 0.85 - 0.05*(f_prime_c - 4000) / 1000;
	else
		betaa = 0.65;

	double e_u = 0.003;
	//------------------------Initializing of parameters for their first value------------------------ -
	int i = 0;
	Force_section(0, i) = (alphaa*f_prime_c*(Ag - As_total) / 1000 + As_total * fy); // Total force in section when bending is zero
	Moment_Section(0, i) = 0.0;                                       // Total moment of the section
	for (double z = 0.5; z >= -5.0; z = z - 0.1) {
		i = i + 1;
		//------------------------------------ - Finding Neutral Axis------------------------------------------
		kd(0, i) = (e_u / (e_u - z * e_y))*d(0, (nL_NS - 1)); //%++++++++++++++++++++++++++++++++++++++++++++++
	//------------------------Initializing of parameters for their ssecond value------------------------ -
		Moment_concrete_total_Comp(0, i) = 0;
		for (int j = 0; j <= nL_NS - 1; j++)
			Moment_Steel_total_layer(j, i) = 0;
		Fc_total_concrete_comp(1, i) = 0;
		for (int j = 0; j <= nL_NS - 1; j++)
			Fs_total_steel_layer(j, i) = 0;
		//Finding Force and moment that come from concrete
		double a = betaa * kd(0, i);
		if (a < h) {
			Fc_total_concrete_comp(0, i) = alphaa * f_prime_c*a*b / 1000.0;
			Moment_concrete_total_Comp(0, i) = Fc_total_concrete_comp(0, i)*(h / 2 - a / 2);
		}
		else {
			Fc_total_concrete_comp(0, i) = alphaa * f_prime_c*h*b / 1000.0;
			Moment_concrete_total_Comp(0, i) = 0;
		}

		Matrix e_s_layer(1, nL_NS);
		for (int j = 0; j <= nL_NS - 1; j++) {
			// Finding Force and moment that come from steel Layers
			if (kd(0, i) >= d(0, j)) {
				e_s_layer(0, j) = abs((kd(0, i) - d(0, j)) / kd(0, i) * e_u);
				if (abs(e_s_layer(0, j)) <= e_y) {
					if (d(0, j) < a)
						Fs_total_steel_layer(j, i) = e_s_layer(0, j) * Es*As(0, j) - 0.85*f_prime_c*As(0, j) / 1000;
					else
						Fs_total_steel_layer(j, i) = e_s_layer(0, j) * Es*As(0, j);
				}
				else {//entering plastic region
					if (d(0, j) < a)
						Fs_total_steel_layer(j, i) = (E1*e_s_layer(0, j) + (fy - E1 * e_y))*As(0, j) - 0.85*f_prime_c*As(0, j) / 1000.0;
					else
						Fs_total_steel_layer(j, i) = (E1*e_s_layer(0, j) + (fy - E1 * e_y))*As(0, j);
				}
				Moment_Steel_total_layer(j, i) = Fs_total_steel_layer(j, i)*(h / 2 - d(0, j));
			}
			else {
				e_s_layer(0, j) = (kd(0, i) - d(0, j)) / kd(0, i)*e_u;
				if (abs(e_s_layer(0, j)) <= e_y)
					Fs_total_steel_layer(j, i) = e_s_layer(0, j)*Es*As(0, j);
				else //entering plastic region
					Fs_total_steel_layer(j, i) = (E1*e_s_layer(0, j) - (fy - E1 * e_y))*As(0, j);

				Moment_Steel_total_layer(j, i) = Fs_total_steel_layer(j, i)*(h / 2 - d(0, j));
			}
		}
		double SUM_Moment_Steel = 0.0;
		for (int jj = 0; jj <= nL_NS - 1; jj++)
			SUM_Moment_Steel = SUM_Moment_Steel + Moment_Steel_total_layer(jj, i);
		Moment_Section(0, i) = Moment_concrete_total_Comp(0, i) + SUM_Moment_Steel;
		double SUM_FS_Steel = 0.0;
		for (int jj = 0; jj <= nL_NS - 1; jj++)
			SUM_FS_Steel = SUM_FS_Steel + Fs_total_steel_layer(jj, i);
		Force_section(0, i) = Fc_total_concrete_comp(0, i) + SUM_FS_Steel;
	}
	Force_section(0, i + 1) = (-As_total * fy); // Total force in section when bending is zero
	Moment_Section(0, i + 1) = 0;               // Total moment of the section
/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                    Needed Input Parameters                                     **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */
	double P = -1.0*Axial_Force;
	//Stiffness====================================================================================================
	double n = 10.0;
	double Iz_mod;
	double E_M = 57.0 * pow(fc * 1000.0, 0.5);

	Iz_mod = (fmax(fmin(0.15 + 0.2*(P_Axial / (B*H*fc)*(La / deff)), 1.0), 0.2)) * (1.0 / 12.0) * B * pow(H, 3.0);

	Kepos_Rot = Keneg_Rot = 6.0 * n * E_M * Iz_mod / La;
	Kepos_Shr = Keneg_Shr = 0.4 * E_M*B*H / (La);
	Kepos_Axil = 29000.0 * AsL_Total / La;
	Keneg_Axil = E_M * (B*H + (Es / E_M - 1.0)*AsL_Total) / La;

	// (Yield Strength)===============================================================================================
	int zz = i + 1;

	if (P >= 0.80*Force_section(0, 0))
		P = 0.80*Force_section(0, 0);
	else if (P <= 0.20*Force_section(0, i + 1))
		P = 0.20*Force_section(0, i + 1);

	for (int k = i + 1; k >= 0; k--) {
		if (P > Force_section(0, k)) {
			zz--;
			continue;
		}
		else
			break;

	}
	fypos_Rot = (P - Force_section(0, zz + 1)) * (Moment_Section(0, zz) - Moment_Section(0, zz + 1)) / (Force_section(0, zz) - Force_section(0, zz + 1)) + Moment_Section(0, zz + 1);
	fyneg_Rot = -1.0 * fypos_Rot;

	Force_section_yield = Force_section;
	Moment_Section_yield = Moment_Section;

}




void
GMG_CMAC2D::P_M_yield_Opt(double Axial_Force) {

	int i = 0;
	for (double z = 0.5; z >= -5.0; z = z - 0.1) {
		i = i + 1;
	}
	double P = -1.0*Axial_Force;

	// (Yield Strength)===============================================================================================
	int zz = i + 1;

	if (P >= 0.80*Force_section_yield(0, 0))
		P = 0.80*Force_section_yield(0, 0);
	else if (P <= 0.20*Force_section_yield(0, i + 1))
		P = 0.20*Force_section_yield(0, i + 1);

	for (int k = i + 1; k >= 0; k--) {
		if (P > Force_section_yield(0, k)) {
			zz--;
			continue;
		}
		else
			break;

	}

	fypos_Rot = (P - Force_section_yield(0, zz + 1)) * (Moment_Section_yield(0, zz) - Moment_Section_yield(0, zz + 1)) / (Force_section_yield(0, zz) - Force_section_yield(0, zz + 1)) + Moment_Section_yield(0, zz + 1);
	fyneg_Rot = -1.0 * fypos_Rot;

}


void
GMG_CMAC2D::P_M_capping(double Axial_Force) {
	// void P_M(double H, double B, double C_C, double fc, double fyL, double dbL, double dbT, double nL_EW, double nL_NS);
	/* ********************************************************************************************** **
	************************************************************************************************* **
	************************************************************************************************* **
	**                                                                                                **
	**                                    Moment Interaction Diagram                                  **
	**                                                                                                **
	************************************************************************************************* **
	************************************************************************************************* **
	** ********************************************************************************************** */

	// Using ACI 318 Assumptions, Rasool Ghorbani
		//           B
		//    < ----------->
		//  ^ | ---------- |
		//  | |  @   @   @ | --->First layer, d1 = deff, n_barT bars
		//  | |            | .
		//  | |            | .
		// H| |  @       @ | . --->n_layM layers in the middle, n_barM bars per each
		//  | |            | .
		//  | |            | .
		//  | |  @   @   @ | --->last layer, d(end), n_barT bars
		//  ^ | ---------- |

		//function[TrueMoment, Force_section, Moment_Section] = PM_INTERACTION_DIAGRAM(H, B, fc, fyL, deff, n_layM, n_barT, n_barM, dbL, P)
		//------------------------ - Concrete properties----------------------------

	const double PI = 3.141592653589793238463;
	double deff = 0.8*H;
	Matrix d(1, nL_NS);
	for (int i = 0; i <= nL_NS - 1; i++) {
		if (i == 0)
			d(0, i) = H - (H - C_C - dbT - dbL / 2);
		else if (i > 0) {
			d(0, i) = (i)*(H - 2.0 * d(0, 0)) / (nL_NS - 1) + d(0, 0);
		}
	}

	double f_prime_c = fc * 1000.0; //psi
	double b = B;
	double h = H;
	double Ag = b * h;

	// ------------------ - Steel Properties(elasto plastic)------------
	double Es = 29000;
	double fy = 1.25*fyL;
	if (Retrofit_flag == 0)                           // Non retrofitted
		fy = 1.25*fyL;
	else if (Retrofit_flag == 1) {                    // FRP
		if (H / B >= 1 && H / B <= 1.5)
			fy = 1.15*fyL + 0.2*(1.5 - H / B)*fyL;
		else if (H / B > 1.5)
			fy = 1.15*fyL;
	}
	else
		fy = 1.3*fyL;                                 // Steel

	double e_y = fy / Es;

	Matrix As(1, nL_NS);
	for (int i = 0; i <= nL_NS - 1; i++) {
		if (i == 0)
			As(0, i) = nL_EW * PI *pow(dbL, 2) / 4;
		else if (i >= 1 && i < nL_NS - 1)
			As(0, i) = 2.0 * PI *pow(dbL, 2) / 4;
		else if (i == nL_NS - 1)
			As(0, i) = nL_EW * PI *pow(dbL, 2) / 4;
	}

	double As_total = 0.0;
	for (int i = 0; i <= nL_NS - 1; i++) {
		As_total = As_total + As(0, i);
	}
	double Plastic_slope_Coeficient = 0.0;
	double E1 = Plastic_slope_Coeficient * Es; // Hardening slope

		// ------------------------Interaction Diagram------------------------ -

	double counter = 1;
	for (double z = 0.5; z >= -5.0; z = z - 0.1) {
		counter += 1;
	}

	Matrix Force_section(1, counter + 2);
	Matrix Moment_Section(1, counter + 2);
	Matrix Moment_concrete_total_Comp(1, counter + 1);
	Matrix Moment_Steel_total_layer(nL_NS, counter + 1);
	Matrix Fc_total_concrete_comp(1, counter + 1);
	Matrix Fs_total_steel_layer(nL_NS, counter + 1);
	Matrix kd(1, counter + 1);

	//---------------------- ACI 318-19 --------------------------
	double alphaa = 0.85;
	double betaa = 0;
	if (f_prime_c <= 4000.0)
		betaa = 0.85;
	else if (f_prime_c > 4000.0 && f_prime_c <= 8000.0)
		betaa = 0.85 - 0.05*(f_prime_c - 4000) / 1000;
	else
		betaa = 0.65;

	double e_u = 0.004;
	if (Retrofit_flag == 0)
		e_u = 0.004;
	else if (Retrofit_flag == 1 || Retrofit_flag == 2)
		e_u = 0.0047;

	//------------------------Initializing of parameters for their first value------------------------ -
	int i = 0;
	Force_section(0, i) = (alphaa*f_prime_c*(Ag - As_total) / 1000 + As_total * fy); // Total force in section when bending is zero
	Moment_Section(0, i) = 0.0;                                       // Total moment of the section
	for (double z = 0.5; z >= -5.0; z = z - 0.1) {
		i = i + 1;
		//------------------------------------ - Finding Neutral Axis------------------------------------------
		kd(0, i) = (e_u / (e_u - z * e_y))*d(0, (nL_NS - 1)); //%++++++++++++++++++++++++++++++++++++++++++++++
	//------------------------Initializing of parameters for their ssecond value------------------------ -
		Moment_concrete_total_Comp(0, i) = 0;
		for (int j = 0; j <= nL_NS - 1; j++)
			Moment_Steel_total_layer(j, i) = 0;
		Fc_total_concrete_comp(1, i) = 0;
		for (int j = 0; j <= nL_NS - 1; j++)
			Fs_total_steel_layer(j, i) = 0;
		//Finding Force and moment that come from concrete
		double a = betaa * kd(0, i);
		if (a < h) {
			Fc_total_concrete_comp(0, i) = alphaa * f_prime_c*a*b / 1000.0;
			Moment_concrete_total_Comp(0, i) = Fc_total_concrete_comp(0, i)*(h / 2 - a / 2);
		}
		else {
			Fc_total_concrete_comp(0, i) = alphaa * f_prime_c*h*b / 1000.0;
			Moment_concrete_total_Comp(0, i) = 0;
		}

		Matrix e_s_layer(1, nL_NS);
		for (int j = 0; j <= nL_NS - 1; j++) {
			// Finding Force and moment that come from steel Layers
			if (kd(0, i) >= d(0, j)) {
				e_s_layer(0, j) = abs((kd(0, i) - d(0, j)) / kd(0, i) * e_u);
				if (abs(e_s_layer(0, j)) <= e_y) {
					if (d(0, j) < a)
						Fs_total_steel_layer(j, i) = e_s_layer(0, j) * Es*As(0, j) - 0.85*f_prime_c*As(0, j) / 1000;
					else
						Fs_total_steel_layer(j, i) = e_s_layer(0, j) * Es*As(0, j);
				}
				else {//entering plastic region
					if (d(0, j) < a)
						Fs_total_steel_layer(j, i) = (E1*e_s_layer(0, j) + (fy - E1 * e_y))*As(0, j) - 0.85*f_prime_c*As(0, j) / 1000.0;
					else
						Fs_total_steel_layer(j, i) = (E1*e_s_layer(0, j) + (fy - E1 * e_y))*As(0, j);
				}
				Moment_Steel_total_layer(j, i) = Fs_total_steel_layer(j, i)*(h / 2 - d(0, j));
			}
			else {
				e_s_layer(0, j) = (kd(0, i) - d(0, j)) / kd(0, i)*e_u;
				if (abs(e_s_layer(0, j)) <= e_y)
					Fs_total_steel_layer(j, i) = e_s_layer(0, j)*Es*As(0, j);
				else //entering plastic region
					Fs_total_steel_layer(j, i) = (E1*e_s_layer(0, j) - (fy - E1 * e_y))*As(0, j);

				Moment_Steel_total_layer(j, i) = Fs_total_steel_layer(j, i)*(h / 2 - d(0, j));
			}
		}
		double SUM_Moment_Steel = 0.0;
		for (int jj = 0; jj <= nL_NS - 1; jj++)
			SUM_Moment_Steel = SUM_Moment_Steel + Moment_Steel_total_layer(jj, i);
		Moment_Section(0, i) = Moment_concrete_total_Comp(0, i) + SUM_Moment_Steel;
		double SUM_FS_Steel = 0.0;
		for (int jj = 0; jj <= nL_NS - 1; jj++)
			SUM_FS_Steel = SUM_FS_Steel + Fs_total_steel_layer(jj, i);
		Force_section(0, i) = Fc_total_concrete_comp(0, i) + SUM_FS_Steel;
		//opserr << "Rasoool " << "  " << " Force_section " << Force_section(0, i) << "  " << " Moment_Section " << Moment_Section(0, i) << "  " << " kd " << kd(0, i) << endln;
	}
	Force_section(0, i + 1) = (-As_total * fy); // Total force in section when bending is zero
	Moment_Section(0, i + 1) = 0;               // Total moment of the section
/* ********************************************************************************************** **
************************************************************************************************* **
************************************************************************************************* **
**                                                                                                **
**                                    Needed Input Parameters                                     **
**                                                                                                **
************************************************************************************************* **
************************************************************************************************* **
** ********************************************************************************************** */
	double P = -1.0*Axial_Force;
	// (Yield Strength)===============================================================================================
	int zz = i + 1;

	if (P >= 0.80*Force_section(0, 0))
		P = 0.80*Force_section(0, 0);
	else if (P <= 0.20*Force_section(0, i + 1))
		P = 0.20*Force_section(0, i + 1);

	for (int k = i + 1; k >= 0; k--) {
		if (P > Force_section(0, k)) {
			zz--;
			continue;
		}
		else
			break;

	}
	fcappos_Rot = (P - Force_section(0, zz + 1)) * (Moment_Section(0, zz) - Moment_Section(0, zz + 1)) / (Force_section(0, zz) - Force_section(0, zz + 1)) + Moment_Section(0, zz + 1);
	fcapneg_Rot = -1.0 * fcappos_Rot;

	Force_section_capping = Force_section;
	Moment_Section_capping = Moment_Section;

}

void
GMG_CMAC2D::P_M_capping_Opt(double Axial_Force) {

	int i = 0;
	for (double z = 0.5; z >= -5.0; z = z - 0.1) {
		i = i + 1;
	}
	double P = -1.0*Axial_Force;
	// (Yield Strength)===============================================================================================
	int zz = i + 1;

	if (P >= 0.80*Force_section_capping(0, 0))
		P = 0.80*Force_section_capping(0, 0);
	else if (P <= 0.20*Force_section_capping(0, i + 1))
		P = 0.20*Force_section_capping(0, i + 1);

	for (int k = i + 1; k >= 0; k--) {
		if (P > Force_section_capping(0, k)) {
			zz--;
			continue;
		}
		else
			break;

	}
	fcappos_Rot = (P - Force_section_capping(0, zz + 1)) * (Moment_Section_capping(0, zz) - Moment_Section_capping(0, zz + 1)) / (Force_section_capping(0, zz) - Force_section_capping(0, zz + 1)) + Moment_Section_capping(0, zz + 1);
	fcapneg_Rot = -1.0 * fcappos_Rot;

}

void
GMG_CMAC2D::P_M_splice(double Axial_Force) {

	// Final needed value (ld)==========================================================================
	const double PI = 3.141592653589793238463;
	double psi_t = 1.0;
	double psi_e = 1.0;
	double psi_s;
	if (dbL <= 0.8)
		psi_s = 0.8;
	else
		psi_s = 1.0;

	double psi_g;
	if (fyL <= 80.0)                       // For grade 60
		psi_g = 1.0;
	else if (fyL <= 100.0 && fyL > 80.0)   // For grade 80
		psi_g = 1.15;
	else                                   // For grade 100
		psi_g = 1.3;

	double cb = fmin(C_C + dbT + dbL / 2.0, (B - 2.0 * (C_C + dbT + dbL / 2.0) / (nL_EW - 1.0)) / 2.0);

	double AsV = PI * pow(dbT, 2.0) / 4.0;
	double Atr = n_leg * AsV;
	double Ktr = 40.0*(Atr) / (S*nL_EW);

	ld = (3.0 / 40.0) * (fyL * 1000.0 * psi_t * psi_e * psi_s * psi_g / (pow(fc*1000.0, 0.5) * fmin((cb + Ktr) / dbL, 2.5)))*dbL;

	// Final needed value (fs & fs_deg)==========================================================================
	fs = fmin(1.25*pow((lb / ld), (2.0 / 3.0))*fyL, fyL); // ksi *********************** ask about fy and fyl / E
	double deff = 0.8*H;
	lb_deg = lb - 2.0 / 3.0 * (0.8*H)/*deff*/;
	fs_deg = 1.25 * pow((lb_deg / ld), (2.0 / 3.0)) * fyL; // ksi *********************** ask about fy and fyl / E


	if (lb != 0.0 && lb < ld && fs < fyL) {
		/* ********************************************************************************************** **
		************************************************************************************************* **
		************************************************************************************************* **
		**                                                                                                **
		**                                    Moment Interaction Diagram                                  **
		**                                                                                                **
		************************************************************************************************* **
		************************************************************************************************* **
		** ********************************************************************************************** */

		// Using ACI 318 Assumptions, Rasool Ghorbani
			//           B
			//    < ----------->
			//  ^ | ---------- |
			//  | |  @   @   @ | --->First layer, d1 = deff, n_barT bars
			//  | |            | .
			//  | |            | .
			// H| |  @       @ | . --->n_layM layers in the middle, n_barM bars per each
			//  | |            | .
			//  | |            | .
			//  | |  @   @   @ | --->last layer, d(end), n_barT bars
			//  ^ | ---------- |

			//function[TrueMoment, Force_section, Moment_Section] = PM_INTERACTION_DIAGRAM(H, B, fc, fyL, deff, n_layM, n_barT, n_barM, dbL, P)
			//------------------------ - Concrete properties----------------------------

		const double PI = 3.141592653589793238463;
		double deff = 0.8*H;
		Matrix d(1, nL_NS);
		for (int i = 0; i <= nL_NS - 1; i++) {
			if (i == 0)
				d(0, i) = H - (H - C_C - dbT - dbL / 2);
			else if (i > 0) {
				d(0, i) = (i)*(H - 2.0 * d(0, 0)) / (nL_NS - 1) + d(0, 0);
			}
		}

		double f_prime_c = fc * 1000.0; //psi
		double b = B;
		double h = H;
		double Ag = b * h;

		// ------------------ - Steel Properties(elasto plastic)------------
		double Es = 29000;
		double fy = fs;
		double e_y = fy / Es;

		Matrix As(1, nL_NS);
		for (int i = 0; i <= nL_NS - 1; i++) {
			if (i == 0)
				As(0, i) = nL_EW * PI *pow(dbL, 2) / 4;
			else if (i >= 1 && i < nL_NS - 1)
				As(0, i) = 2.0 * PI *pow(dbL, 2) / 4;
			else if (i == nL_NS - 1)
				As(0, i) = nL_EW * PI *pow(dbL, 2) / 4;
		}
		double As_total = 0.0;
		for (int i = 0; i <= nL_NS - 1; i++) {
			As_total = As_total + As(0, i);
		}

		double Plastic_slope_Coeficient = 0.0;
		double E1 = Plastic_slope_Coeficient * Es; // Hardening slope

			// ------------------------Interaction Diagram------------------------ -

		double counter = 1;
		for (double z = 0.5; z >= -5.0; z = z - 0.1) {
			counter += 1;
		}

		Matrix Force_section(1, counter + 2);
		Matrix Moment_Section(1, counter + 2);
		Matrix Moment_concrete_total_Comp(1, counter + 1);
		Matrix Moment_Steel_total_layer(nL_NS, counter + 1);
		Matrix Fc_total_concrete_comp(1, counter + 1);
		Matrix Fs_total_steel_layer(nL_NS, counter + 1);
		Matrix kd(1, counter + 1);

		//---------------------- ACI 318-19 --------------------------
		double alphaa = 0.85;
		double betaa = 0;
		if (f_prime_c <= 4000.0)
			betaa = 0.85;
		else if (f_prime_c > 4000.0 && f_prime_c <= 8000.0)
			betaa = 0.85 - 0.05*(f_prime_c - 4000) / 1000;
		else
			betaa = 0.65;

		double e_u = 0.003;
		//------------------------Initializing of parameters for their first value------------------------ -
		int i = 0;
		Force_section(0, i) = (alphaa*f_prime_c*(Ag - As_total) / 1000 + As_total * fy); // Total force in section when bending is zero
		Moment_Section(0, i) = 0.0;                                       // Total moment of the section
		for (double z = 0.5; z >= -5.0; z = z - 0.1) {
			i = i + 1;
			//------------------------------------ - Finding Neutral Axis------------------------------------------
			kd(0, i) = (e_u / (e_u - z * e_y))*d(0, (nL_NS - 1)); //%++++++++++++++++++++++++++++++++++++++++++++++

		//------------------------Initializing of parameters for their ssecond value------------------------ -
			Moment_concrete_total_Comp(0, i) = 0;
			for (int j = 0; j <= nL_NS - 1; j++)
				Moment_Steel_total_layer(j, i) = 0;
			Fc_total_concrete_comp(1, i) = 0;
			for (int j = 0; j <= nL_NS - 1; j++)
				Fs_total_steel_layer(j, i) = 0;
			//Finding Force and moment that come from concrete
			double a = betaa * kd(0, i);
			if (a < h) {
				Fc_total_concrete_comp(0, i) = alphaa * f_prime_c*a*b / 1000.0;
				Moment_concrete_total_Comp(0, i) = Fc_total_concrete_comp(0, i)*(h / 2 - a / 2);
			}
			else {
				Fc_total_concrete_comp(0, i) = alphaa * f_prime_c*h*b / 1000.0;
				Moment_concrete_total_Comp(0, i) = 0;
			}

			Matrix e_s_layer(1, nL_NS);
			for (int j = 0; j <= nL_NS - 1; j++) {
				// Finding Force and moment that come from steel Layers
				if (kd(0, i) >= d(0, j)) {
					e_s_layer(0, j) = abs((kd(0, i) - d(0, j)) / kd(0, i) * e_u);
					if (abs(e_s_layer(0, j)) <= e_y) {
						if (d(0, j) < a)
							Fs_total_steel_layer(j, i) = e_s_layer(0, j) * Es*As(0, j) - 0.85*f_prime_c*As(0, j) / 1000;
						else
							Fs_total_steel_layer(j, i) = e_s_layer(0, j) * Es*As(0, j);
					}
					else {//entering plastic region
						if (d(0, j) < a)
							Fs_total_steel_layer(j, i) = (E1*e_s_layer(0, j) + (fy - E1 * e_y))*As(0, j) - 0.85*f_prime_c*As(0, j) / 1000.0;
						else
							Fs_total_steel_layer(j, i) = (E1*e_s_layer(0, j) + (fy - E1 * e_y))*As(0, j);
					}
					Moment_Steel_total_layer(j, i) = Fs_total_steel_layer(j, i)*(h / 2 - d(0, j));
				}
				else {
					e_s_layer(0, j) = (kd(0, i) - d(0, j)) / kd(0, i)*e_u;
					if (abs(e_s_layer(0, j)) <= e_y)
						Fs_total_steel_layer(j, i) = e_s_layer(0, j)*Es*As(0, j);
					else //entering plastic region
						Fs_total_steel_layer(j, i) = (E1*e_s_layer(0, j) - (fy - E1 * e_y))*As(0, j);

					Moment_Steel_total_layer(j, i) = Fs_total_steel_layer(j, i)*(h / 2 - d(0, j));
				}
			}
			double SUM_Moment_Steel = 0.0;
			for (int jj = 0; jj <= nL_NS - 1; jj++)
				SUM_Moment_Steel = SUM_Moment_Steel + Moment_Steel_total_layer(jj, i);
			Moment_Section(0, i) = Moment_concrete_total_Comp(0, i) + SUM_Moment_Steel;
			double SUM_FS_Steel = 0.0;
			for (int jj = 0; jj <= nL_NS - 1; jj++)
				SUM_FS_Steel = SUM_FS_Steel + Fs_total_steel_layer(jj, i);
			Force_section(0, i) = Fc_total_concrete_comp(0, i) + SUM_FS_Steel;
		}
		Force_section(0, i + 1) = (-As_total * fy); // Total force in section when bending is zero
		Moment_Section(0, i + 1) = 0;               // Total moment of the section
	/* ********************************************************************************************** **
	************************************************************************************************* **
	************************************************************************************************* **
	**                                                                                                **
	**                                    Needed Input Parameters                                     **
	**                                                                                                **
	************************************************************************************************* **
	************************************************************************************************* **
	** ********************************************************************************************** */
		double P = -1.0*Axial_Force;
		// (Yield Strength)===============================================================================================
		int zz = i + 1;
		//opserr << "Rasoool2 " << "  " << " zz" << "  " << zz << endln;

		if (P >= 0.80*Force_section(0, 0))
			P = 0.80*Force_section(0, 0);
		else if (P <= 0.20*Force_section(0, i + 1))
			P = 0.20*Force_section(0, i + 1);

		for (int k = i + 1; k >= 0; k--) {
			if (P > Force_section(0, k)) {
				zz--;
				continue;
			}
			else
				break;

		}
		fypos_splice_Rot = (P - Force_section(0, zz + 1)) * (Moment_Section(0, zz) - Moment_Section(0, zz + 1)) / (Force_section(0, zz) - Force_section(0, zz + 1)) + Moment_Section(0, zz + 1);
		fyneg_splice_Rot = -1.0 * fypos_splice_Rot;
		dypos_splice = fypos_splice_Rot / Kepos_Rot;
		dyneg_splice = -1.0*dypos_splice;
		fcappos_splice_Rot = 1.05 * fypos_splice_Rot;
		fcapneg_splice_Rot = 1.05 * fyneg_splice_Rot;

		Force_section_splice = Force_section;
		Moment_Section_splice = Moment_Section;
	}


	else if (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL) {
		dypos_splice = fypos_Rot / Kepos_Rot;
		dyneg_splice = fyneg_Rot / Keneg_Rot;
	}


	else { //(lb == 0.0 || lb >= ld) {
		fypos_splice_Rot = 2.0*fypos_Rot;
		dypos_splice = fypos_splice_Rot / Kepos_Rot;
		fyneg_splice_Rot = -1.0*fypos_splice_Rot;
		dyneg_splice = -1.0*dypos_splice;
		fypos_splice_deg_Rot = 2.0*fypos_splice_Rot;
		fyneg_splice_deg_Rot = -1.0*fypos_splice_deg_Rot;
		fcappos_splice_Rot = 1.2*fypos_splice_Rot;
		fcapneg_splice_Rot = -1.0*fcappos_splice_Rot;
	}
}

void
GMG_CMAC2D::P_M_splice_Opt(double Axial_Force) {
	if (lb != 0.0 && lb < ld && fs < fyL) {
		int i = 0;
		for (double z = 0.5; z >= -5.0; z = z - 0.1) {
			i = i + 1;
		}
		double P = -1.0*Axial_Force;
		// (Yield Strength)===============================================================================================
		int zz = i + 1;

		if (P >= 80 * Force_section_splice(0, 0))
			P = 0.80*Force_section_splice(0, 0);
		else if (P <= 0.20*Force_section_splice(0, i + 1))
			P = 0.20*Force_section_splice(0, i + 1);

		for (int k = i + 1; k >= 0; k--) {
			if (P > Force_section_splice(0, k)) {
				zz--;
				continue;
			}
			else
				break;

		}
		fypos_splice_Rot = (P - Force_section_splice(0, zz + 1)) * (Moment_Section_splice(0, zz) - Moment_Section_splice(0, zz + 1)) / (Force_section_splice(0, zz) - Force_section_splice(0, zz + 1)) + Moment_Section_splice(0, zz + 1);
		fyneg_splice_Rot = -1.0 * fypos_splice_Rot;
		dypos_splice = fypos_splice_Rot / Kepos_Rot;
		dyneg_splice = -1.0*dypos_splice;
		fcappos_splice_Rot = 1.05 * fypos_splice_Rot;
		fcapneg_splice_Rot = 1.05 * fyneg_splice_Rot;

	}


	else if (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL) {
		dypos_splice = fypos_Rot / Kepos_Rot;
		dyneg_splice = fyneg_Rot / Keneg_Rot;
	}


	else { //(lb == 0.0 || lb >= ld) {
		fypos_splice_Rot = 2.0*fypos_Rot;
		dypos_splice = fypos_splice_Rot / Kepos_Rot;
		fyneg_splice_Rot = -1.0*fypos_splice_Rot;
		dyneg_splice = -1.0*dypos_splice;
		fypos_splice_deg_Rot = 2.0*fypos_splice_Rot;
		fyneg_splice_deg_Rot = -1.0*fypos_splice_deg_Rot;
		fcappos_splice_Rot = 1.2*fypos_splice_Rot;
		fcapneg_splice_Rot = -1.0*fcappos_splice_Rot;
	}
}

void
GMG_CMAC2D::Shear_strength(double Axial_Force) {
	//Flexure-Shear================================================================================================
	//Calculating Vcol============================================================================
	double P = -1.0*Axial_Force;
	const double PI = 3.141592653589793238463;
	//Transverse reinf. calculation
	double AsV = PI * pow(dbT, 2.0) / 4.0;
	double roT = (n_leg*AsV) / (B*S);

	//Longitudinal reinf. calculation
	double AsL = PI * pow(dbL, 2.0) / 4.0;
	double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);

	// Geometry
	double deff_EW = 0.8*H;

	Vcol = 0.0;
	double alpha;
	double zeta;
	double kd;
	double A_str;
	double teta_str;

	double teta_cr = 40.0;

	// alpha calculatoin :
	if ((S / (deff_EW*(1.0 / tan(PI*teta_cr / 180.0)))) <= 1.0)
		alpha = 1.0;
	else if ((S / (deff_EW*(1.0 / tan(PI*teta_cr / 180.0)))) >= 1.25) {
		alpha = 0;
	}
	else/* if ((S / (deff_EW*(1.0 / tan(PI*teta_cr[i] / 180.0)))) > 0.75 && (S / (deff_EW*(1.0 / tan(PI*teta_cr[i] / 180.0)))) < 1.0)*/ {
		alpha = 1.0 + (-1.0 / 0.25)* (S / (deff_EW*(1.0 / tan(PI*teta_cr / 180.0))) - 1.0);
	}

	if (La / deff_EW < 2.0) {
		zeta = fmin(40.3 / pow(fc*1000.0, 0.5), 0.52);
		kd = (0.25 + fmax(0.85 * P / (B * H * fc), 0.0)) * H;                          // opposit sign for axial force
		A_str = kd * B;
		teta_str = atan((H - 2.0 * kd / 3.0) / L);
		Vcol = fmax((1.85 * zeta * fc * 1000.0 * A_str * sin(teta_str)) / 1000.0, alpha * (n_leg*AsV*fyT * deff_EW*(1.0 / tan(PI* teta_cr / 180.0)) / S));
	}

	else {

		Vcol = fmax((fmin(((7.5 * pow(fc*1000.0, 0.5) / (La / deff_EW)) * pow(fmax((1.0 + P * 1000.0 / (7.5 * sqrt(fc * 1000.0)*B*H)), 0.0), 0.5) * 0.8*B*H), 4.0*pow(fc*1000.0, 0.5)*B*H) + alpha * (n_leg*AsV*fyT * 1000.0 * deff_EW*(1.0 / tan(PI* teta_cr / 180.0)) / S)), 0.0) / 1000.0;

	}

	fypos_Shr = Vcol;
	dypos_Shr = Vcol / Kepos_Shr;
	dcappos_Shr = dypos_Shr + max((0.0025 + 1.4*roT)*La, 0.1*dypos_Shr); // I had to use different sign for the Axial load here
	fyneg_Shr = -Vcol;
	dyneg_Shr = -Vcol / Kepos_Shr;
	dcapneg_Shr = -dcappos_Shr;

}


void
GMG_CMAC2D::Shear_strength_Opt(double Axial_Force) {
	//Flexure-Shear================================================================================================
	//Calculating Vcol============================================================================
	double P = -1.0*Axial_Force;
	const double PI = 3.141592653589793238463;
	//Transverse reinf. calculation
	double AsV = PI * pow(dbT, 2.0) / 4.0;
	double roT = (n_leg*AsV) / (B*S);

	//Longitudinal reinf. calculation
	double AsL = PI * pow(dbL, 2.0) / 4.0;
	double roL = (2.0 * nL_EW*AsL + 2 * (nL_NS - 2.0)*AsL) / (B * H);

	// Geometry
	double deff_EW = 0.8*H;

	Vcol = 0.0;
	double alpha;
	double zeta;
	double kd;
	double A_str;
	double teta_str;
	double teta_cr = 40.0;

	// alpha calculatoin :
	if ((S / (deff_EW*(1.0 / tan(PI*teta_cr / 180.0)))) <= 1.0)
		alpha = 1.0;
	else if ((S / (deff_EW*(1.0 / tan(PI*teta_cr / 180.0)))) >= 1.25) {
		alpha = 0;
	}
	else/* if ((S / (deff_EW*(1.0 / tan(PI*teta_cr[i] / 180.0)))) > 0.75 && (S / (deff_EW*(1.0 / tan(PI*teta_cr[i] / 180.0)))) < 1.0)*/ {
		alpha = 1.0 + (-1.0 / 0.25)* (S / (deff_EW*(1.0 / tan(PI*teta_cr / 180.0))) - 1.0);
	}

	if (La / deff_EW < 2.0) {
		zeta = fmin(40.3 / pow(fc*1000.0, 0.5), 0.52);
		kd = (0.25 + fmax(0.85 * P / (B * H * fc), 0.0)) * H;                          // opposit sign for axial force
		A_str = kd * B;
		teta_str = atan((H - 2.0 * kd / 3.0) / L);
		Vcol = fmax((1.85 * zeta * fc * 1000.0 * A_str * sin(teta_str)) / 1000.0, alpha * (n_leg*AsV*fyT * deff_EW*(1.0 / tan(PI* teta_cr / 180.0)) / S));
	}

	else {
		Vcol = fmax((fmin(((7.5 * pow(fc*1000.0, 0.5) / (La / deff_EW)) * pow(fmax((1.0 + P * 1000.0 / (7.5 * sqrt(fc * 1000.0)*B*H)), 0.0), 0.5) * 0.8*B*H), 4.0*pow(fc*1000.0, 0.5)*B*H) + alpha * (n_leg*AsV*fyT * 1000.0 * deff_EW*(1.0 / tan(PI* teta_cr / 180.0)) / S)), 0.0) / 1000.0;
	}

	fypos_Shr = Vcol;
	dypos_Shr = Vcol / Kepos_Shr;
	dcappos_Shr = dypos_Shr + max((0.0025 + 1.4*roT)*La, 0.1*dypos_Shr); // I had to use different sign for the Axial load here
	fyneg_Shr = -Vcol;
	dyneg_Shr = -Vcol / Kepos_Shr;
	dcapneg_Shr = -dcappos_Shr;
}


double GMG_CMAC2D::getStress_Multi(void)
{

	static int Counter_multi_Strs = -1;
	Counter_multi_Strs++;
	if (Counter_multi_Strs == 3)
		Counter_multi_Strs = 0;
	return Stress_vector(Counter_multi_Strs);
}

double GMG_CMAC2D::getTangent_Multi(void)
{
	static int Counter_multi_Tan = -1;
	Counter_multi_Tan++;
	if (Counter_multi_Tan == 3)
		Counter_multi_Tan = 0;
	return Tang_vector(Counter_multi_Tan);
}

double GMG_CMAC2D::getInitialTangent_Multi(void)
{
	static int Counter_multi_Tan = -1;
	Counter_multi_Tan++;
	if (Counter_multi_Tan == 3)
		Counter_multi_Tan = 0;
	return Tang_vector(Counter_multi_Tan);
}

double GMG_CMAC2D::getDampTangent_Multi(void)
{
	static int Counter_multi_Damp = -1;
	Counter_multi_Damp++;
	if (Counter_multi_Damp == 3)
		Counter_multi_Damp = 0;
	return Damp_vector(Counter_multi_Damp);
}

double GMG_CMAC2D::getStrain_Multi(void)
{
	static int Counter_multi_Strain = -1;
	Counter_multi_Strain++;
	if (Counter_multi_Strain == 3)
		Counter_multi_Strain = 0;

	return Strain_vector(Counter_multi_Strain);
}

void GMG_CMAC2D::DamageOutput(void)
{
	Counter_MACE++;

	Damage_Data(Counter_MACE - 1, 0) = Counter_MACE;
	Damage_Data(Counter_MACE - 1, 1) = getTag();
}

void GMG_CMAC2D::DamageOutput_Hardening(void)
{
	for (int count = 0; count < Counter_MACE; count++) {
		if (count % 2 == 0 && Damage_Data(count, 1) == getTag()) {
			if (flag_entering_hardening_Flex_Rot == 1) {
				if (d_Rot >= 0) {
					if (VyE_to_Vcol_Ratio <= 0.6) {
						Damage_Data(count, 7) = fmin(dpeakmax / dcappos_Rot * 100.0, 100.0); // Flx
						Damage_Data(count, 9) = 0; // Flx-Shr
					}
					else {
						//Damage_Data(count, 7) = 0; // Flx
						Damage_Data(count, 7) = fmin(dpeakmax / dcappos_Rot * 100.0, 100.0); // Flx
						Damage_Data(count, 9) = fmin(dpeakmax / dcappos_FlexShr_Rot * 100.0, 100.0); // Flx-Shr
					}
					if (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL)
						Damage_Data(count, 10) = fmin(dpeakmax / dcappos_splice_Rot * 100.0, 100.0); // Flx-splice
					else
						Damage_Data(count, 10) = 0;
				}
				if (d_Rot < 0) {
					if (VyE_to_Vcol_Ratio <= 0.6) {
						Damage_Data(count, 7) = fmin(dpeakmin / dcapneg_Rot * 100.0, 100.0); // Flx
						Damage_Data(count, 9) = 0; // Flx-Shr
					}
					else {
						//Damage_Data(count, 7) = 0; // Flx
						Damage_Data(count, 7) = fmin(dpeakmin / dcapneg_Rot * 100.0, 100.0); // Flx
						Damage_Data(count, 9) = fmin(dpeakmin / dcapneg_FlexShr_Rot * 100.0, 100.0); // Flx-Shr
					}
					if (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL)
						Damage_Data(count, 10) = fmin(dpeakmin / dcapneg_splice_Rot * 100.0, 100.0); // Flx-splice
					else
						Damage_Data(count, 10) = 0;
				}
			}

			if (flag_entering_hardening_Shr == 1) {
				if (d_Shr >= 0)
					Damage_Data(count, 8) = fmin(dpeakmax_Shr / dcappos_Shr * 100.0, 100.0);
				else
					Damage_Data(count, 8) = fmin(dpeakmin_Shr / dcapneg_Shr * 100.0, 100.0);
			}

			if (flag_entering_hardening_Splice_Rot == 1) {
				if (d_Rot >= 0)
					Damage_Data(count, 11) = fmin(dpeakmax / dcappos_splice_Rot * 100.0, 100.0);
				else
					Damage_Data(count, 11) = fmin(dpeakmin / dcapneg_splice_Rot * 100.0, 100.0);
			}
		}
	}
}

void GMG_CMAC2D::DamageOutput_FailureType(void)
{
	for (int count = 0; count < Counter_MACE; count++) {
		if (count % 2 == 0 && Damage_Data(count, 1) == getTag()) {


			if (flag_entering_hardening_Flex_Rot == 1 && Cflag_entering_hardening_Flex_Rot != 1) {
				Damage_Data(count, 2) = d_Rot;
				Damage_Data(count, 4) = 100.0;
				Damage_Data(count, 5) = fmin(fabs(F_Shr / fypos_Shr)*100.0, 100.0);
				if (lb != 0.0 && lb < ld && fs < fyL) {
					Damage_Data(count, 6) = fmin(fabs(F_Rot / fypos_splice_Rot)*100.0, 100.0);
				}
			}

			if (flag_entering_hardening_Splice_Rot == 1 && Cflag_entering_hardening_Splice_Rot != 1) {
				Damage_Data(count, 2) = d_Rot;
				Damage_Data(count, 4) = fmin(fabs(F_Rot / fypos_Rot)*100.0, 100.0);
				Damage_Data(count, 5) = fmin(fabs(F_Shr / fypos_Shr)*100.0, 100.0);
				Damage_Data(count, 6) = 100.0;
			}

			if (flag_entering_hardening_Shr == 1 && Cflag_entering_hardening_Shr != 1) {
				Damage_Data(count, 3) = d_Shr;
				Damage_Data(count, 4) = fmin(fabs(F_Rot / fypos_Rot)*100.0, 100.0);
				Damage_Data(count, 5) = 100.0;
				if ((lb != 0.0 && lb < ld && fs < fyL) /*|| (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL)*/) {
					Damage_Data(count, 6) = fmin(fabs(F_Rot / fypos_splice_Rot)*100.0, 100.0);
				}
			}

			if (flag_Flx_Filure_Rot == 1 && Cflag_Flx_Filure_Rot != 1) {

				Damage_Data(count, 12) = 1;

				Damage_Data(count, 13) = d_Rot;
				if (d_Rot >= 0) {

					Damage_Data(count, 7) = 100;  // For flexure

					if (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL) {

						Damage_Data(count, 10) = fmin(d_Rot / dcappos_splice_Rot * 100.0, 100.0);  // For Flx-Splice
					}
					else {

						Damage_Data(count, 10) = 0;
					}
				}
				if (d_Rot < 0) {

					Damage_Data(count, 7) = 100;  // For flexure

					if (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL) {

						Damage_Data(count, 10) = fmin(d_Rot / dcapneg_splice_Rot * 100.0, 100.0);  // For Flx-Splice
					}
					else {

						Damage_Data(count, 10) = 0;
					}
				}
			}


			if (flag_Filure_Shr == 1 && Cflag_Filure_Shr != 1) {

				Damage_Data(count, 8) = 100;  // For Shr
				Damage_Data(count, 12) = 2;
				Damage_Data(count, 13) = d_Shr;

			}


			if (flag_FlexShr_Filure_Rot == 1 && Cflag_FlexShr_Filure_Rot != 1) {

				Damage_Data(count, 12) = 3;

				Damage_Data(count, 13) = d_Rot;
				if (d_Rot >= 0) {

					Damage_Data(count, 9) = 100;  // For Flx-Shr
					if (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL) {

						Damage_Data(count, 10) = fmin(d_Rot / dcappos_splice_Rot * 100.0, 100.0);  // For Flx-Splice
					}
					else {

						Damage_Data(count, 10) = 0;
					}
				}
				if (d_Rot < 0) {

					Damage_Data(count, 9) = 100;  // For Flx-Shr
					if (lb != 0.0 && lb < ld && fs == fyL && fs_deg < fyL) {

						Damage_Data(count, 10) = fmin(d_Rot / dcapneg_splice_Rot * 100.0, 100.0);  // For Flx-Splice
					}
					else {

						Damage_Data(count, 10) = 0;
					}
				}
			}


			if (flag_Splice_Filure_Rot == 1 && Cflag_Splice_Filure_Rot != 1) {

				Damage_Data(count, 12) = 4;

				Damage_Data(count, 13) = d_Rot;
				if (d_Rot >= 0) {

					Damage_Data(count, 11) = 100;  // For Splice
				}
				else {

					Damage_Data(count, 11) = 100;  // For Splice
				}
			}


			if (flag_FlexSplice_Filure_Rot == 1 && Cflag_FlexSplice_Filure_Rot != 1) {

				Damage_Data(count, 12) = 5;

				Damage_Data(count, 13) = d_Rot;
				if (d_Rot >= 0) {

					if (VyE_to_Vcol_Ratio <= 0.6) {
						Damage_Data(count, 7) = fmin(d_Rot / dcappos_Rot * 100.0, 100.0);  // For flexure

						Damage_Data(count, 9) = 0;  // For Flx-Shr
					}
					else {

						Damage_Data(count, 7) = fmin(d_Rot / dcappos_Rot * 100.0, 100.0);  // For flexure

						Damage_Data(count, 9) = fmin(d_Rot / dcappos_FlexShr_Rot * 100.0, 100.0);  // For Flx-Shr
					}

					Damage_Data(count, 10) = 100.0;  // For Flx-Splice
				}
				if (d_Rot < 0) {

					if (VyE_to_Vcol_Ratio <= 0.6) {
						Damage_Data(count, 7) = fmin(d_Rot / dcapneg_Rot * 100.0, 100.0);  // For flexure

						Damage_Data(count, 9) = 0;  // For Flx-Shr
					}
					else {

						Damage_Data(count, 7) = fmin(d_Rot / dcapneg_Rot * 100.0, 100.0);  // For flexure

						Damage_Data(count, 9) = fmin(d_Rot / dcapneg_FlexShr_Rot * 100.0, 100.0);  // For Flx-Shr
					}

					Damage_Data(count, 10) = 100;  // For Flx-Splice
				}
			}

		}
	}
}


void GMG_CMAC2D::DamageOutput_PostCapping(void)
{
	for (int count = 0; count < Counter_MACE; count++) {
		if (count % 2 == 0 && Damage_Data(count, 1) == getTag()) {
			if (flag_entering_hardening_Flex_Rot == 1) {
				Damage_Data(count, 14) = fmax(fabs(ffmax - R_fcappos) / fabs(R_frespos - R_fcappos) * 100.0, fabs(ffmin - R_fcapneg) / fabs(R_fresneg - R_fcapneg) * 100.0);
			}

			if (flag_entering_hardening_Shr == 1) {
				Damage_Data(count, 14) = fmax(fabs(ffmax_Shr - R_fcappos_Shr) / fabs(R_frespos_Shr - R_fcappos_Shr) * 100.0, fabs(ffmin_Shr - R_fcapneg_Shr) / fabs(R_fresneg_Shr - R_fcapneg_Shr) * 100.0);
			}

			if (flag_entering_hardening_Splice_Rot == 1) {
				Damage_Data(count, 14) = fmax(fabs(ffmax - R_fcappos) / fabs(R_frespos - R_fcappos) * 100.0, fabs(ffmin - R_fcapneg) / fabs(R_fresneg - R_fcapneg) * 100.0);
			}
		}
	}
}


void GMG_CMAC2D::DamageOutput_PostResidual(void)
{
	for (int count = 0; count < Counter_MACE; count++) {
		if (count % 2 == 0 && Damage_Data(count, 1) == getTag()) {

			Damage_Data(count, 14) = 100;

			if (flag_entering_residual_Flex_Rot == 1 && Cflag_entering_residual_Flex_Rot != 1) {

				Damage_Data(count, 15) = d_Rot;
			}

			if (flag_entering_residual_Shr == 1 && Cflag_entering_residual_Shr != 1) {

				Damage_Data(count, 15) = d_Shr;
			}
		}
	}
}