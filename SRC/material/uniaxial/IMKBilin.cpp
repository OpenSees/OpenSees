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

#include <math.h>
#include <IMKBilin.h>
#include <elementAPI.h>
#include <Vector.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <algorithm>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using namespace std;

static int numIMKBilinMaterials = 0;

void *
OPS_IMKBilin(void)
{
	if (numIMKBilinMaterials == 0) {
		numIMKBilinMaterials++;
		OPS_Error("Mod. IMK Model with Bilinear Hysteretic Response - Code by A. ELKADY\n", 1);
	}

	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial = 0;

	int    iData[1];
	double dData[21];
	int numData = 1;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial IMKBilin tag" << endln;
		return 0;
	}

	numData = 21;

	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid Args want: uniaxialMaterial IMKBilin tag? Ke? ";
		opserr << "Theta_p_pos? Theta_pc_pos? Theta_u_pos? Mpe_pos? MmaxMpe_pos? ResM_pos? ";
		opserr << "Theta_p_neg? Theta_pc_neg? Theta_u_neg? Mpe_neg? MmaxMpe_neg? ResM_neg? ";
		opserr << "LamdaS?  LamdaC? LamdaK? Cs? Cc? Ck? D_pos? D_neg? ";
		return 0;
	}


	// Parsing was successful, allocate the material
	theMaterial = new IMKBilin(iData[0],
		dData[0],
		dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
		dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
		dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20]);

	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type IMKBilin Material\n";
		return 0;
	}

	return theMaterial;
}

IMKBilin::IMKBilin(int tag, double p_Ke,
	double p_Theta_p_pos0, double p_Theta_pc_pos0, double p_Theta_u_pos0, double p_Mpe_pos0, double p_MmaxMpe_pos0, double p_ResM_pos0,
	double p_Theta_p_neg0, double p_Theta_pc_neg0, double p_Theta_u_neg0, double p_Mpe_neg0, double p_MmaxMpe_neg0, double p_ResM_neg0,
	double p_LAMBDA_S, double p_LAMBDA_C, double p_LAMBDA_K, double p_c_S, double p_c_C, double p_c_K, double p_D_pos, double p_D_neg)
	:UniaxialMaterial(tag, 0), Ke(p_Ke),
	Theta_p_pos0(p_Theta_p_pos0), Theta_pc_pos0(p_Theta_pc_pos0), Theta_u_pos0(p_Theta_u_pos0), Mpe_pos0(p_Mpe_pos0), MmaxMpe_pos0(p_MmaxMpe_pos0), ResM_pos0(p_ResM_pos0),
	Theta_p_neg0(p_Theta_p_neg0), Theta_pc_neg0(p_Theta_pc_neg0), Theta_u_neg0(p_Theta_u_neg0), Mpe_neg0(p_Mpe_neg0), MmaxMpe_neg0(p_MmaxMpe_neg0), ResM_neg0(p_ResM_neg0),
	LAMBDA_S(p_LAMBDA_S), LAMBDA_C(p_LAMBDA_C), LAMBDA_K(p_LAMBDA_K), c_S(p_c_S), c_C(p_c_C), c_K(p_c_K), D_pos(p_D_pos), D_neg(p_D_neg)
{
	this->revertToStart();
}

IMKBilin::IMKBilin()
	:UniaxialMaterial(0, 0), Ke(0),
	Theta_p_pos0(0), Theta_pc_pos0(0), Theta_u_pos0(0), Mpe_pos0(0), MmaxMpe_pos0(0), ResM_pos0(0),
	Theta_p_neg0(0), Theta_pc_neg0(0), Theta_u_neg0(0), Mpe_neg0(0), MmaxMpe_neg0(0), ResM_neg0(0),
	LAMBDA_S(0), LAMBDA_C(0), LAMBDA_K(0), c_S(0), c_C(0), c_K(0), D_pos(0), D_neg(0)
{
	this->revertToStart();
}

IMKBilin::~IMKBilin()
{
	// does nothing
}

int IMKBilin::setTrialStrain(double strain, double strainRate)
{
	//all variables to the last commit
	this->revertToLastCommit();

	//state determination algorithm: defines the current force and tangent stiffness
	U = strain; //set trial displacement
	Ri_1 = Ri;
	Mi_1 = Mi;
	Di_1 = Di;
	Ri = U;

	//cout << " Ri_1=" << Ri_1 << " Ri=" << Ri << " Ex_Flag=" << Excursion_Flag << " Rev_Flag=" << Reversal_Flag << " Mr+_Flag=" << Mrpos_Flag << " Mr-_Flag=" << Mrneg_Flag << " En_Flag=" << Energy_Flag << endln;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Theta_y_pos0 = Mpe_pos0 / Ke;
	Theta_max_pos0 = Theta_y_pos0 + Theta_p_pos0;
	slope_p_pos0 = (Mmax_pos0 - Mpe_pos0) / (Theta_p_pos0);
	slope_pc_pos0 = Mmax_pos0 / (Theta_pc_pos0);
	Mpe_pos0 = Mpe_pos0;
	Mmax_pos0 = MmaxMpe_pos0 * Mpe_pos0;
	MpeProject_pos0 = Mmax_pos0 - slope_p_pos0  * Theta_max_pos0;
	MmaxProject_pos0 = Mmax_pos0 + slope_pc_pos0 * Theta_max_pos0;

	Theta_y_neg0 = Mpe_neg0 / Ke;
	Theta_max_neg0 = Theta_y_neg0 + Theta_p_neg0;
	slope_p_neg0 = (Mmax_neg0 - Mpe_neg0) / (Theta_p_neg0);
	slope_pc_neg0 = Mmax_neg0 / (Theta_pc_neg0);
	Mpe_neg0 = Mpe_neg0;
	Mmax_neg0 = MmaxMpe_neg0 * Mpe_neg0;
	MpeProject_neg0 = Mmax_neg0 - slope_p_neg0  * Theta_max_neg0;
	MmaxProject_neg0 = Mmax_neg0 + slope_pc_neg0 * Theta_max_neg0;

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%% INITIALIZE CURRENT BACKBONE VALUES AS PREVIOUS %%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	double 		  Mmax_pos_j = Mmax_pos_j_1;
	double 		   Mpe_pos_j = Mpe_pos_j_1;
	double  MpeProject_pos_j = MpeProject_pos_j_1;
	double MmaxProject_pos_j = MmaxProject_pos_j_1;
	double 	   slope_p_pos_j = slope_p_pos_j_1;
	double 	  slope_pc_pos_j = slope_pc_pos_j_1;
	double 	   Theta_y_pos_j = Theta_y_pos_j_1;
	double 	 Theta_max_pos_j = Theta_max_pos_j_1;

	double 		  Mmax_neg_j = Mmax_neg_j_1;
	double 		   Mpe_neg_j = Mpe_neg_j_1;
	double  MpeProject_neg_j = MpeProject_neg_j_1;
	double MmaxProject_neg_j = MmaxProject_neg_j_1;
	double 	   slope_p_neg_j = slope_p_neg_j_1;
	double 	  slope_pc_neg_j = slope_pc_neg_j_1;
	double	   Theta_y_neg_j = Theta_y_neg_j_1;
	double 	 Theta_max_neg_j = Theta_max_neg_j_1;

	double  Mi_boundary_Pos = 0.0;
	double  Mi_boundary_Neg = 0.0;

	double QuarterFlag, Rintrsct_K, DISP_Rev;
	double beta_S_j, beta_C_j, beta_K_j, beta_F_j;
	double K_j;
	double Ki, slope_pi, slope_pci, Theta_yi, Theta_maxi, Mmaxi, MpeProjecti, MmaxProjecti, Theta_pci;
	double BackboneRangeFlag;

	QuarterFlag = 0;
	Reversal_Flag = 0;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%  CORE CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// Find the direction of the current increment "Di": 1:if Ri is moving right    and -1: if Ri is moving left
	if (Ri >= Ri_1) {
		Di = 1;
	}
	else {
		Di = -1;
	}

	//  Simple Notation for current step parameters
	Mi = Mi_1 + K_j_1*(Ri - Ri_1);

	//  Check for Fail Flag
	if ((Ri >= Theta_u_pos0)) {
		Fail_FlagPos = 1.0;
	}
	if ((Ri <= -Theta_u_neg0)) {
		Fail_FlagNeg = 1.0;
	}

	//  Get Information before first Yield
	if ((Mi > Mpe_pos0) && (Yield_Flag == 0.0)) {
		Yield_Flag = 1.0;
	}

	//  Check if previous point was a reversal point
	if (Di_1 / Di < 0) {
		Reversal_Flag = 1;
		Rreversal = Ri_1;
		Mreversal = Mi_1;
	}

	if (Reversal_Flag == 1) {
		Rintrsct_K = Rreversal - Mreversal / K_j_1;
		DISP_Rev = Energy_total - Energy_Excrsni_1 + 0.5*Mreversal *(Rintrsct_K - Rreversal);
		// Update Loading/Unloading Stiffness
		beta_K_j = pow((DISP_Rev / (2 * Ref_Energy_K - Energy_total + 0.5*Mreversal * (Rintrsct_K - Rreversal))), c_K);
		K_j = K_j_1*(1 - beta_K_j);
		if ((Mrpos_Flag == 1) || (Mrneg_Flag == 1)) {
			K_j = 0.5*Ke;
		}
/* 		// Update Paramters for Smooth Transition Option
		beta_F_j = pow((DISP_Rev / (2 * Ref_Energy_F - Energy_total + 0.5 * Mreversal * (Rintrsct_K - Rreversal))), c_F);
		if (Di > 0) {
			beta_Sx = pow(((DISP_Rev + 0.5*Mreversal* (Rintrsct_K - Rreversal)) / (Ref_Energy_S - Energy_total)), c_S);
			slope_p_x = slope_p_pos_j_1 * (1.0 - beta_Sx * D_pos);
			Mpe_x = Mpe_pos_j_1 	* (1.0 - beta_Sx * D_pos);
			Theta_y_x = Mpe_x / K_j;
			MpeProject_x = Mpe_x - slope_p_x * Theta_y_x;
		}
		else {
			beta_Sx = pow(((DISP_Rev + 0.5*Mreversal * (Rintrsct_K - Rreversal)) / (Ref_Energy_S - Energy_total)), c_S);
			slope_p_x = slope_p_neg_j_1 * (1.0 - beta_Sx * D_neg);
			Mpe_x = Mpe_neg_j_1 	* (1.0 - beta_Sx * D_neg);
			Theta_y_x = Mpe_x / K_j;
			MpeProject_x = Mpe_x - slope_p_x * Theta_y_x;
		} */
	}
	else {
		// beta_F_j = beta_F_j_1; 
		beta_K_j = beta_K_j_1;
		K_j = K_j_1;
	}

	//  Calculate Backbone parameters at current excursion based on Energy Dissipated in the previous Excursion
	if (Excursion_Flag == 1.0) {
		beta_S_j = pow((Energy_Excrsn / (Ref_Energy_S - Energy_total)), c_S);
		beta_C_j = pow((Energy_Excrsn / (Ref_Energy_C - Energy_total)), c_C);

		if (Ri - Ri_1 >= 0.0) {
			//  Update Mpe, Mmax Projection, slope_p, and slope_pc for Current Step
			Mpe_pos_j = Mpe_pos_j_1 		  * (1.0 - beta_S_j * D_pos);
			MmaxProject_pos_j = MmaxProject_pos_j_1 * (1.0 - beta_C_j * D_pos);
			slope_p_pos_j = slope_p_pos_j_1 	  * (1.0 - beta_S_j * D_pos);
			if (Mr_pos0 == 0.0) {
				slope_pc_pos_j = slope_pc_pos0 * (MmaxProject_pos_j - Mr_pos0) / MmaxProject_pos_j;
			}
			else {
				slope_pc_pos_j = slope_pc_pos0 * (Mpe_pos_j - Mr_pos0) / (Mpe_pos0 - Mr_pos0);
			}

			//  Calculate Rotation at Capping Point (Intersection of the two Slopes)
			Theta_y_pos_j = Mpe_pos_j / K_j;
			MpeProject_pos_j = Mpe_pos_j - slope_p_pos_j * Theta_y_pos_j;
			Theta_max_pos_j = fabs((MmaxProject_pos_j - MpeProject_pos_j) / (slope_pc_pos_j + slope_p_pos_j));
			Mmax_pos_j = MpeProject_pos_j + Theta_max_pos_j * slope_p_pos_j;
		}
		else {
			//  Update Mpe, Mmax Projection, slope_p, and slope_pc for Current Step
			Mpe_neg_j = Mpe_neg_j_1 		  * (1.0 - beta_S_j * D_neg);
			MmaxProject_neg_j = MmaxProject_neg_j_1 * (1.0 - beta_C_j * D_neg);
			slope_p_neg_j = slope_p_neg_j_1 	  * (1.0 - beta_S_j * D_neg);
			if (Mr_neg0 == 0.0) {
				slope_pc_neg_j = slope_pc_neg0 * (MmaxProject_neg_j - Mr_neg0) / MmaxProject_neg_j;
			}
			else {
				slope_pc_neg_j = slope_pc_neg0 * (Mpe_neg_j - Mr_neg0) / (Mpe_neg0 - Mr_neg0);
			}

			//  Calculate Rotation at Capping Point (Intersection of the two Slopes)
			Theta_y_neg_j = Mpe_neg_j / K_j;
			MpeProject_neg_j = Mpe_neg_j - slope_p_neg_j * Theta_y_neg_j;
			Theta_max_neg_j = fabs((MmaxProject_neg_j - MpeProject_neg_j) / (slope_pc_neg_j + slope_p_neg_j));
			Mmax_neg_j = MpeProject_neg_j + Theta_max_neg_j * slope_p_neg_j;
		}

	}
	else {

		beta_S_j = beta_S_j_1;
		beta_C_j = beta_C_j_1;

		if (Di >= 0.0) {
			Mpe_pos_j = Mpe_pos_j_1;
			MpeProject_pos_j = MpeProject_pos_j_1;
			MmaxProject_pos_j = MmaxProject_pos_j_1;
			Mmax_pos_j = Mmax_pos_j_1;
			slope_p_pos_j = slope_p_pos_j_1;
			slope_pc_pos_j = slope_pc_pos_j_1;
			Theta_y_pos_j = Theta_y_pos_j_1;
			Theta_max_pos_j = Theta_max_pos_j_1;
		}
		else {
			Mpe_neg_j = Mpe_neg_j_1;
			MpeProject_neg_j = MpeProject_neg_j_1;
			MmaxProject_neg_j = MmaxProject_neg_j_1;
			Mmax_neg_j = Mmax_neg_j_1;
			slope_p_neg_j = slope_p_neg_j_1;
			slope_pc_neg_j = slope_pc_neg_j_1;
			Theta_y_neg_j = Theta_y_neg_j_1;
			Theta_max_neg_j = Theta_max_neg_j_1;
		}
	}

	// If the residual moment is reached in a given direction, Ovverride the values of Mmax, Theta_max, slope_p and slope_pc
	if (Di >= 0.0) {
		if (Mmax_pos_j < Mr_pos0) {
			Mmax_pos_j = Mr_pos0;
			slope_pc_pos_j = pow(10., -6);
			slope_p_pos_j = pow(10., -6);
			Theta_max_pos_j = pow(10., -6);
		}
	}
	else {
		if (Mmax_neg_j < Mr_neg0) {
			Mmax_neg_j = Mr_neg0;
			slope_pc_neg_j = pow(10., -6);
			slope_p_neg_j = pow(10., -6);
			Theta_max_neg_j = pow(10., -6);
		}
	}

	//  Simple and unified notation for current bacbone parameters
	if (Di >= 0.0) {
		Ki = K_j;
		slope_pi = slope_p_pos_j;
		slope_pci = slope_pc_pos_j;
		Theta_yi = Theta_y_pos_j;
		Theta_maxi = Theta_max_pos_j;
		Mmaxi = Mmax_pos_j;
		MpeProjecti = MpeProject_pos_j;
		MmaxProjecti = MmaxProject_pos_j;
		Theta_pci = Mmaxi / slope_pci;
	}
	else {
		Ki = K_j;
		slope_pi = slope_p_neg_j;
		slope_pci = slope_pc_neg_j;
		Theta_yi = Theta_y_neg_j;
		Theta_maxi = Theta_max_neg_j;
		Mmaxi = Mmax_neg_j;
		MpeProjecti = MpeProject_neg_j;
		MmaxProjecti = MmaxProject_neg_j;
		Theta_pci = Mmaxi / slope_pci;
	}


	//  Moment Calculation Based on unloading/reloading stiffeness
	Mi = Mi_1 + Ki * (Ri - Ri_1);

	//  Location Flags
	if ((Ri >= 0.0) && (Mi > 0.0)) {
		QuarterFlag = 1;
	}
	else if ((Ri >= 0.0) && (Mi < 0.0)) {
		QuarterFlag = 2;
	}
	else if ((Ri <= 0.0) && (Mi < 0.0)) {
		QuarterFlag = 3;
	}
	else if ((Ri <= 0.0) && (Mi > 0.0)) {
		QuarterFlag = 4;
	}

	//  Get Boundary Moment at Current Step Based on Current BackBone Curve
	if ((Ri >= 0.0) && (fabs(Ri) <= Theta_maxi)) {
		Mi_boundary_Pos = MpeProjecti + slope_pi * Ri;
		Mi_boundary_Neg = -MpeProjecti + slope_pi * Ri;
	}
	else if ((Ri >= 0.0) && (fabs(Ri) > Theta_maxi)) {
		Mi_boundary_Pos = MmaxProjecti - slope_pci * Ri;
		Mi_boundary_Neg = -MpeProjecti + slope_pi  * Ri;
	}
	else if ((Ri <= 0.0) && (fabs(Ri) <= Theta_maxi)) {
		Mi_boundary_Pos = MpeProjecti + slope_pi * Ri;
		Mi_boundary_Neg = -MpeProjecti + slope_pi * Ri;
	}
	else if ((Ri <= 0.0) && (fabs(Ri) > Theta_maxi)) {
		Mi_boundary_Pos = MpeProjecti + slope_pi  * Ri;
		Mi_boundary_Neg = -MmaxProjecti - slope_pci * Ri;
	}

	//  Check and Modify Boundary Moment if it Exceeds Mr
	if ((QuarterFlag == 3) && (Mi_boundary_Neg >= 0.0)) {
		Mi_boundary_Neg = 0.0;
	}

	if ((QuarterFlag == 3) && (fabs(Mi_boundary_Neg) <= Mr_neg0))
	{
		Mi_boundary_Neg = -Mr_neg0;
		Mrneg_Flag = 1;
	}

	if ((QuarterFlag == 1) && (Mi_boundary_Pos <= 0.0)) {
		Mi_boundary_Pos = 0.0;
	}

	if ((QuarterFlag == 1) && (fabs(Mi_boundary_Pos) <= Mr_pos0))
	{
		Mi_boundary_Pos = Mr_pos0;
		Mrpos_Flag = 1;
	}

//	///////////////////////////////////////////////////////////////////////
//	///////////////////// CODE FOR SMOOTH SPLINE //////////////////////////
//	///////////////////////////////////////////////////////////////////////
//	//  Get Boundary Moment in the Transition Region at Current Step Based on Spline Curve
//
//	double Dx[4];
//	double Rxx[1];
//	double Mxx[1];
//	double Rintrsct, RintrsctL, MintrsctL, RintrsctR, MintrsctR;
//	double b1, b2;
//	double RintrsctXaxis;
//
//	if ((Di < 0)) {
//		b1 = Mreversal - Ki*Rreversal;
//		b2 = -MpeProject_x;
//		RintrsctXaxis = -b1 / Ki;
//		Rintrsct = -(b2 - b1) / (slope_p_x - Ki);
//
//		RintrsctL = Rintrsct - Roffset / (1 - beta_F_j);
//		if ((RintrsctL > Theta_maxi)) {
//			RintrsctL = Theta_maxi;
//		}
//		MintrsctL = b2 + slope_pi*RintrsctL;
//
//		RintrsctR = Rintrsct + n*Roffset / (1 - beta_F_j);
//		if (RintrsctR > Rreversal) {
//			RintrsctR = Rreversal;
//		}
//		MintrsctR = b1 + Ki*RintrsctR;
//
//		if ((Ri > RintrsctL) && (Ri < RintrsctR)) {
//			// R and M vectors for PCHIP interpolation
//			double Rx[4] = { RintrsctL - 0.01,   RintrsctL,   RintrsctR,          RintrsctR + 0.01 };
//			double Mx[4] = { b2 + slope_pi*(RintrsctL - 0.01),  MintrsctL,   MintrsctR,   b1 + Ki*(RintrsctR + 0.01) };
//			spline_pchip_set(4, Rx, Mx, Dx);
//			Rxx[0] = Ri;
//			spline_pchip_set(4, Rx, Mx, Dx);
//			spline_pchip_val(4, Rx, Mx, Dx, 1, Rxx, Mxx);
//			Mi_boundary_Pos = Mxx[0];
//			Mi_boundary_Neg = Mxx[0];
//			//cout <<  " R_rev=" << Rreversal  <<" R_L=" << RintrsctL  <<" R_R=" << RintrsctR << " Mbound-=" << Mi_boundary_Neg << " Mbound+=" << Mi_boundary_Neg<< " Mi="<< Mi << " beta="<<beta_F_j <<endl;
//		}
//	}
//
//	if ((Di >= 0)) {
//		b1 = Mreversal - Ki*Rreversal;
//		b2 = MpeProject_x;
//		RintrsctXaxis = -b1 / Ki;
//		Rintrsct = -(b2 - b1) / (slope_p_x - Ki);
//
//		RintrsctL = Rintrsct - n*Roffset / (1 - beta_F_j);
//		if (RintrsctL < Rreversal) {
//			RintrsctL = Rreversal;
//		}
//		MintrsctL = b1 + Ki*RintrsctL;
//
//		RintrsctR = Rintrsct + Roffset / (1 - beta_F_j);
//		if ((RintrsctR > Theta_maxi)) {
//			RintrsctR = Theta_maxi;
//		}
//		MintrsctR = b2 + slope_pi*RintrsctR;
//
//		if ((Ri > RintrsctL) && (Ri < RintrsctR)) {
//			// R and M vectors for PCHIP interpolation
//			double Rx[4] = { RintrsctL - 0.01,   RintrsctL,   RintrsctR,                RintrsctR + 0.01 };
//			double Mx[4] = { b1 + Ki*(RintrsctL - 0.01),  MintrsctL,   MintrsctR,   b2 + slope_pi*(RintrsctR + 0.01) };
//			spline_pchip_set(4, Rx, Mx, Dx);
//			Rxx[0] = Ri;
//			spline_pchip_set(4, Rx, Mx, Dx);
//			spline_pchip_val(4, Rx, Mx, Dx, 1, Rxx, Mxx);
//			Mi_boundary_Pos = Mxx[0];
//			Mi_boundary_Neg = Mxx[0];
//		}
//	}
//
//	///////////////////////////////////////////////////////////////////////
//	///////////////////////////////////////////////////////////////////////


	// If Failure took place in a given direction (Fail_Flag_dir=1), Set the Boundary Moment in the opposite direction to Mr
	if ((Ri < 0) && (Di > 0) && (Fail_FlagNeg == 1)) {
		Mi_boundary_Pos = Mr_pos0;
	}
	if ((Ri > 0) && (Di < 0) && (Fail_FlagPos == 1)) {
		Mi_boundary_Neg = -Mr_neg0;
	}

	// If the residual moment is reached in a given direction (Mrdir_Flag=1), Set the Boundary Moment in both directions to Mr
	if ((Mrpos_Flag == 1) || (Mrneg_Flag == 1)) {
		Mi_boundary_Pos = Mr_pos0;
		Mi_boundary_Neg = -Mr_neg0;
	}

	// %%%%%%% Current Step Moment Calculation %%%%%%%
	// If current moment based on unloading/reloading Ki is larger than the boundary moment, set it equal to the boundary moment
	if (Mi > Mi_boundary_Pos) {
		Mi = Mi_boundary_Pos;
	}
	if (Mi < Mi_boundary_Neg) {
		Mi = Mi_boundary_Neg;
	}

	// if fail flag is reached in a given loading direction, set current moment equal to zero
	if ((Ri > 0.0) && (Fail_FlagPos == 1.0)) {
		Mi = 0.0;
	}
	if ((Ri < 0.0) && (Fail_FlagNeg == 1.0)) {
		Mi = 0.0;
	}

	if (Energy_Flag == 1) {
		Mi = 0.0;
	}

	// %%%%%%%%%%%%% Energy Calculation %%%%%%%%%%%%%
	Energy_total = Energy_total + (Mi + Mi_1) * 0.5 * (Ri - Ri_1); //  total energy dissipated till current incremental step

	// Energy calculation at each new excursion
	if (Mi / Mi_1 < 0.0) {
		Energy_Excrsn = Energy_total - Energy_Excrsni_1;	// total energy dissipated in current excursion
		Energy_Excrsni_1 = Energy_total;					// total energy dissipated in previous excursion
		Excursion_Flag = 1;
	}
	else {
		Excursion_Flag = 0.0;
	}

	// Check if the Component inheret Reference Energy is Consumed
	if (Excursion_Flag == 1) {
		if ((Energy_total >= Ref_Energy_S) || (Energy_total >= Ref_Energy_C)) {
			Energy_Flag = 1;
		}
		if ((beta_S_j > 1) || (beta_C_j > 1)) {
			Energy_Flag = 1;
		}
	}
	else if (Reversal_Flag == 1) {
		if (Energy_total >= Ref_Energy_K) {
			Energy_Flag = 1;
		}
		if (beta_K_j > 1) {
			Energy_Flag = 1;
		}
	}

	// %%%%%%%%%% PREPARE RETURN VALUES %%%%%%%%%%%%%

	if (Ri >= Ri_1) {
		K_j_1 = K_j;
		Theta_y_pos_j_1 = Theta_y_pos_j;
		Theta_max_pos_j_1 = Theta_max_pos_j;
		slope_p_pos_j_1 = slope_p_pos_j;
		slope_pc_pos_j_1 = slope_pc_pos_j;
		Mpe_pos_j_1 = Mpe_pos_j;
		MpeProject_pos_j_1 = MpeProject_pos_j;
		Mmax_pos_j_1 = Mmax_pos_j;
		MmaxProject_pos_j_1 = MmaxProject_pos_j;
		Theta_u_pos0 = Theta_u_pos0;
	}
	else {
		K_j_1 = K_j;
		Theta_y_neg_j_1 = Theta_y_neg_j;
		Theta_max_neg_j_1 = Theta_max_neg_j;
		slope_p_neg_j_1 = slope_p_neg_j;
		slope_pc_neg_j_1 = slope_pc_neg_j;
		Mpe_neg_j_1 = Mpe_neg_j;
		MpeProject_neg_j_1 = MpeProject_neg_j;
		Mmax_neg_j_1 = Mmax_neg_j;
		MmaxProject_neg_j_1 = MmaxProject_neg_j;
		Theta_u_neg0 = Theta_u_neg0;
	}

	beta_S_j_1 = beta_S_j;
	beta_C_j_1 = beta_C_j;
	beta_K_j_1 = beta_K_j;
	//beta_F_j_1 = beta_F_j;

	// Tangent Stiffeness Calculation
	if (Mi == Mr_pos0) {
		TangentK = pow(10., -6);
	}

	if (Ri == Ri_1) {
		TangentK = Ke;
		Mi = Mi_1;
	}
	else {
		TangentK = (Mi - Mi_1) / (Ri - Ri_1);
	}


	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%% END OF MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	return 0;
}

double IMKBilin::getStress(void)
{
	//cout << " getStress" << endln;
	return (Mi);
}

double IMKBilin::getTangent(void)
{
	//cout << " getTangent" << endln;
	return (TangentK);
}

double IMKBilin::getInitialTangent(void)
{
	//cout << " getInitialTangent" << endln;
	return (Ke);
}

double IMKBilin::getStrain(void)
{
	//cout << " getStrain" << endln;
	return (U);
}

int IMKBilin::commitState(void)
{
	//cout << " commitState" << endln;

	//commit trial  variables

	cU = U;

	cRi = Ri;
	cMi = Mi;
	cDi = Di;
	cRi_1 = Ri_1;
	cMi_1 = Mi_1;
	cDi_1 = Di_1;

	cRreversal = Rreversal;
	cMreversal = Mreversal;
	cTangentK = TangentK;

	cK_j_1 = K_j_1;
	cTheta_y_pos_j_1 = Theta_y_pos_j_1;
	cTheta_max_pos_j_1 = Theta_max_pos_j_1;
	cslope_p_pos_j_1 = slope_p_pos_j_1;
	cslope_pc_pos_j_1 = slope_pc_pos_j_1;
	cMpe_pos_j_1 = Mpe_pos_j_1;
	cMpeProject_pos_j_1 = MpeProject_pos_j_1;
	cMmax_pos_j_1 = Mmax_pos_j_1;
	cMmaxProject_pos_j_1 = MmaxProject_pos_j_1;

	cTheta_y_neg_j_1 = Theta_y_neg_j_1;
	cTheta_max_neg_j_1 = Theta_max_neg_j_1;
	cslope_p_neg_j_1 = slope_p_neg_j_1;
	cslope_pc_neg_j_1 = slope_pc_neg_j_1;
	cMpe_neg_j_1 = Mpe_neg_j_1;
	cMpeProject_neg_j_1 = MpeProject_neg_j_1;
	cMmax_neg_j_1 = Mmax_neg_j_1;
	cMmaxProject_neg_j_1 = MmaxProject_neg_j_1;

	cbeta_S_j_1 = beta_S_j_1;
	cbeta_C_j_1 = beta_C_j_1;
	cbeta_K_j_1 = beta_K_j_1;
	cbeta_F_j_1 = beta_F_j_1;

	cExcursion_Flag = Excursion_Flag;
	cReversal_Flag = Reversal_Flag;
	cYield_Flag = Yield_Flag;
	cFail_FlagPos = Fail_FlagPos;
	cFail_FlagNeg = Fail_FlagNeg;
	cMrpos_Flag = Mrpos_Flag;
	cMrneg_Flag = Mrneg_Flag;
	cEnergy_Flag = Energy_Flag;

	cEnergy_Excrsni_1 = Energy_Excrsni_1;
	cEnergy_Excrsn = Energy_Excrsn;
	cEnergy_Rev = Energy_Rev;
	cEnergy_total = Energy_total;

	return 0;
}

int IMKBilin::revertToLastCommit(void)
{
	//cout << " revertToLastCommit" << endln;

	//the opposite of commit trial history variables
	U = cU;

	Ri = cRi;
	Mi = cMi;
	Di = cDi;
	Ri_1 = cRi_1;
	Mi_1 = cMi_1;
	Di_1 = cDi_1;

	Rreversal = cRreversal;
	Mreversal = cMreversal;
	TangentK = cTangentK;

	K_j_1 = cK_j_1;
	Theta_y_pos_j_1 = cTheta_y_pos_j_1;
	Theta_max_pos_j_1 = cTheta_max_pos_j_1;
	slope_p_pos_j_1 = cslope_p_pos_j_1;
	slope_pc_pos_j_1 = cslope_pc_pos_j_1;
	Mpe_pos_j_1 = cMpe_pos_j_1;
	MpeProject_pos_j_1 = cMpeProject_pos_j_1;
	Mmax_pos_j_1 = cMmax_pos_j_1;
	MmaxProject_pos_j_1 = cMmaxProject_pos_j_1;

	Theta_y_neg_j_1 = cTheta_y_neg_j_1;
	Theta_max_neg_j_1 = cTheta_max_neg_j_1;
	slope_p_neg_j_1 = cslope_p_neg_j_1;
	slope_pc_neg_j_1 = cslope_pc_neg_j_1;
	Mpe_neg_j_1 = cMpe_neg_j_1;
	MpeProject_neg_j_1 = cMpeProject_neg_j_1;
	Mmax_neg_j_1 = cMmax_neg_j_1;
	MmaxProject_neg_j_1 = cMmaxProject_neg_j_1;

	beta_S_j_1 = cbeta_S_j_1;
	beta_C_j_1 = cbeta_C_j_1;
	beta_K_j_1 = cbeta_K_j_1;
	beta_F_j_1 = cbeta_F_j_1;

	Excursion_Flag = cExcursion_Flag;
	Reversal_Flag = cReversal_Flag;
	Yield_Flag = cYield_Flag;
	Fail_FlagPos = cFail_FlagPos;
	Fail_FlagNeg = cFail_FlagNeg;
	Mrpos_Flag = cMrpos_Flag;
	Mrneg_Flag = cMrneg_Flag;
	Energy_Flag = cEnergy_Flag;

	Energy_Excrsni_1 = cEnergy_Excrsni_1;
	Energy_Excrsn = cEnergy_Excrsn;
	Energy_Rev = cEnergy_Rev;
	Energy_total = cEnergy_total;

	return 0;
}

int IMKBilin::revertToStart(void)
{
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ONE TIME CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

	Theta_y_pos0 = Mpe_pos0 / Ke;
	Theta_max_pos0 = Theta_y_pos0 + Theta_p_pos0;
	slope_p_pos0 = (Mmax_pos0 - Mpe_pos0) / (Theta_p_pos0);
	slope_pc_pos0 = Mmax_pos0 / (Theta_pc_pos0);
	Mpe_pos0 = Mpe_pos0;
	Mmax_pos0 = MmaxMpe_pos0 * Mpe_pos0;
	MpeProject_pos0 = Mmax_pos0 - slope_p_pos0  * Theta_max_pos0;
	MmaxProject_pos0 = Mmax_pos0 + slope_pc_pos0 * Theta_max_pos0;

	Theta_y_neg0 = Mpe_neg0 / Ke;
	Theta_max_neg0 = Theta_y_neg0 + Theta_p_neg0;
	slope_p_neg0 = (Mmax_neg0 - Mpe_neg0) / (Theta_p_neg0);
	slope_pc_neg0 = Mmax_neg0 / (Theta_pc_neg0);
	Mpe_neg0 = Mpe_neg0;
	Mmax_neg0 = MmaxMpe_neg0 * Mpe_neg0;
	MpeProject_neg0 = Mmax_neg0 - slope_p_neg0  * Theta_max_neg0;
	MmaxProject_neg0 = Mmax_neg0 + slope_pc_neg0 * Theta_max_neg0;

	Mr_pos0 = ResM_pos0*Mpe_pos0;
	Mr_neg0 = ResM_neg0*Mpe_neg0;

	Ref_Energy_S = LAMBDA_S*Mpe_pos0;
	Ref_Energy_C = LAMBDA_C*Mpe_pos0;
	Ref_Energy_K = LAMBDA_K*Mpe_pos0;
	Ref_Energy_F = LAMBDA_F*Mpe_pos0;

	K_j_1 = Ke;
	Theta_y_pos_j_1 = Mpe_pos0 / Ke;
	Theta_max_pos_j_1 = Theta_y_pos0 + Theta_p_pos0;
	slope_p_pos_j_1 = (Mmax_pos0 - Mpe_pos0) / (Theta_p_pos0);
	slope_pc_pos_j_1 = Mmax_pos0 / (Theta_pc_pos0);
	Mpe_pos_j_1 = Mpe_pos0;
	Mmax_pos_j_1 = MmaxMpe_pos0*Mpe_pos0;
	MpeProject_pos_j_1 = Mmax_pos_j_1 - slope_p_pos_j_1 *Theta_max_pos_j_1;
	MmaxProject_pos_j_1 = Mmax_pos_j_1 + slope_pc_pos_j_1*Theta_max_pos_j_1;

	Theta_y_neg_j_1 = Mpe_neg0 / Ke;
	Theta_max_neg_j_1 = Theta_y_neg0 + Theta_p_neg0;
	slope_p_neg_j_1 = (Mmax_neg0 - Mpe_neg0) / (Theta_p_neg0);
	slope_pc_neg_j_1 = Mmax_neg0 / (Theta_pc_neg0);
	Mpe_neg_j_1 = Mpe_neg0;
	Mmax_neg_j_1 = MmaxMpe_neg0*Mpe_neg0;
	MpeProject_neg_j_1 = Mmax_neg_j_1 - slope_p_neg_j_1 *Theta_max_neg_j_1;
	MmaxProject_neg_j_1 = Mmax_neg_j_1 + slope_pc_neg_j_1*Theta_max_neg_j_1;

	//initially I zero everything   
	U = cU = 0;

	Ri = 0;
	Mi = 0;
	Di = 0;
	Ri_1 = 0;
	Mi_1 = 0;
	Di_1 = 0;
	TangentK = Ke;
	cTangentK = Ke;

	cRi = 0;
	cMi = 0;
	cDi = 0;
	cRi_1 = 0;
	cMi_1 = 0;
	cDi_1 = 0;

	Rreversal = cRreversal = 0;
	Mreversal = cMreversal = 0;

	beta_Sx = 0.0;
	slope_p_x = slope_p_pos_j_1;
	Mpe_x = Mpe_pos_j_1;
	Theta_y_x = Mpe_x / Ke;
	MpeProject_x = Mpe_x - slope_p_x * Theta_y_x;

	beta_S_j_1 = cbeta_S_j_1 = 0;
	beta_C_j_1 = cbeta_C_j_1 = 0;
	beta_K_j_1 = cbeta_K_j_1 = 0;
	beta_F_j_1 = cbeta_F_j_1 = 0;

	Excursion_Flag = cExcursion_Flag = 0;
	Reversal_Flag = cReversal_Flag = 0;
	Yield_Flag = cYield_Flag = 0;
	Fail_FlagPos = cFail_FlagPos = 0;
	Fail_FlagNeg = cFail_FlagNeg = 0;
	Mrpos_Flag = cMrpos_Flag = 0;
	Mrneg_Flag = cMrneg_Flag = 0;
	Energy_Flag = cEnergy_Flag = 0;

	Energy_Excrsni_1 = cEnergy_Excrsni_1 = 0;
	Energy_Excrsn = cEnergy_Excrsn = 0;
	Energy_Rev = cEnergy_Rev = 0;
	Energy_total = cEnergy_total = 0;

	cK_j_1 = Ke;
	cTheta_y_pos_j_1 = Mpe_pos0 / Ke;
	cTheta_max_pos_j_1 = Theta_y_pos0 + Theta_p_pos0;
	cslope_p_pos_j_1 = (Mmax_pos0 - Mpe_pos0) / (Theta_p_pos0);
	cslope_pc_pos_j_1 = Mmax_pos0 / (Theta_pc_pos0);
	cMpe_pos_j_1 = Mpe_pos0;
	cMmax_pos_j_1 = MmaxMpe_pos0*Mpe_pos0;
	cMpeProject_pos_j_1 = Mmax_pos_j_1 - slope_p_pos_j_1 *Theta_max_pos_j_1;
	cMmaxProject_pos_j_1 = Mmax_pos_j_1 + slope_pc_pos_j_1*Theta_max_pos_j_1;

	cTheta_y_neg_j_1 = Mpe_neg0 / Ke;
	cTheta_max_neg_j_1 = Theta_y_neg0 + Theta_p_neg0;
	cslope_p_neg_j_1 = (Mmax_neg0 - Mpe_neg0) / (Theta_p_neg0);
	cslope_pc_neg_j_1 = Mmax_neg0 / (Theta_pc_neg0);
	cMpe_neg_j_1 = Mpe_neg0;
	cMmax_neg_j_1 = MmaxMpe_neg0*Mpe_neg0;
	cMpeProject_neg_j_1 = Mmax_neg_j_1 - slope_p_neg_j_1 *Theta_max_neg_j_1;
	cMmaxProject_neg_j_1 = Mmax_neg_j_1 + slope_pc_neg_j_1*Theta_max_neg_j_1;

	//cout << " revertToStart:" << endln; //<< " U=" << U << " Ri=" << Ri << " TanK=" << TangentK << endln;

	return 0;
}

UniaxialMaterial *
IMKBilin::getCopy(void)
{
	IMKBilin *theCopy = new IMKBilin(this->getTag(), Ke,
		Theta_p_pos0, Theta_pc_pos0, Theta_u_pos0, Mpe_pos0, MmaxMpe_pos0, ResM_pos0,
		Theta_p_neg0, Theta_pc_neg0, Theta_u_neg0, Mpe_neg0, MmaxMpe_neg0, ResM_neg0,
		LAMBDA_S, LAMBDA_C, LAMBDA_K, c_S, c_C, c_K, D_pos, D_neg);

	//cout << " getCopy" << endln;

	theCopy->Mr_pos0 = Mr_pos0;
	theCopy->Mr_neg0 = Mr_neg0;

	theCopy->U = U;
	theCopy->cU = cU;

	theCopy->Ri = Ri;
	theCopy->Mi = Mi;
	theCopy->Di = Di;
	theCopy->Ri_1 = Ri_1;
	theCopy->Mi_1 = Mi_1;
	theCopy->Di_1 = Di_1;

	theCopy->Rreversal = Rreversal;
	theCopy->Mreversal = Mreversal;
	theCopy->TangentK = TangentK;

	theCopy->K_j_1 = K_j_1;
	theCopy->Theta_y_pos_j_1 = Theta_y_pos_j_1;
	theCopy->Theta_max_pos_j_1 = Theta_max_pos_j_1;
	theCopy->slope_p_pos_j_1 = slope_p_pos_j_1;
	theCopy->slope_pc_pos_j_1 = slope_pc_pos_j_1;
	theCopy->Mpe_pos_j_1 = Mpe_pos_j_1;
	theCopy->MpeProject_pos_j_1 = MpeProject_pos_j_1;
	theCopy->Mmax_pos_j_1 = Mmax_pos_j_1;
	theCopy->MmaxProject_pos_j_1 = MmaxProject_pos_j_1;

	theCopy->Theta_y_neg_j_1 = Theta_y_neg_j_1;
	theCopy->Theta_max_neg_j_1 = Theta_max_neg_j_1;
	theCopy->slope_p_neg_j_1 = slope_p_neg_j_1;
	theCopy->slope_pc_neg_j_1 = slope_pc_neg_j_1;
	theCopy->Mpe_neg_j_1 = Mpe_neg_j_1;
	theCopy->MpeProject_neg_j_1 = MpeProject_neg_j_1;
	theCopy->Mmax_neg_j_1 = Mmax_neg_j_1;
	theCopy->MmaxProject_neg_j_1 = MmaxProject_neg_j_1;

	theCopy->beta_S_j_1 = beta_S_j_1;
	theCopy->beta_C_j_1 = beta_C_j_1;
	theCopy->beta_K_j_1 = beta_K_j_1;
	theCopy->beta_F_j_1 = beta_F_j_1;

	theCopy->Excursion_Flag = Excursion_Flag;
	theCopy->Reversal_Flag = Reversal_Flag;
	theCopy->Yield_Flag = Yield_Flag;
	theCopy->Fail_FlagPos = Fail_FlagPos;
	theCopy->Fail_FlagNeg = Fail_FlagNeg;
	theCopy->Mrpos_Flag = Mrpos_Flag;
	theCopy->Mrneg_Flag = Mrneg_Flag;
	theCopy->Energy_Flag = Energy_Flag;

	theCopy->Energy_Excrsni_1 = Energy_Excrsni_1;
	theCopy->Energy_Excrsn = Energy_Excrsn;
	theCopy->Energy_Rev = Energy_Rev;
	theCopy->Energy_total = Energy_total;


	theCopy->cMr_pos0 = cMr_pos0;
	theCopy->cMr_neg0 = cMr_neg0;

	theCopy->cRi = cRi;
	theCopy->cMi = cMi;
	theCopy->cDi = cDi;
	theCopy->cRi_1 = cRi_1;
	theCopy->cMi_1 = cMi_1;
	theCopy->cDi_1 = cDi_1;

	theCopy->cRreversal = cRreversal;
	theCopy->cMreversal = cMreversal;
	theCopy->cTangentK = cTangentK;

	theCopy->cK_j_1 = cK_j_1;
	theCopy->cTheta_y_pos_j_1 = cTheta_y_pos_j_1;
	theCopy->cTheta_max_pos_j_1 = cTheta_max_pos_j_1;
	theCopy->cslope_p_pos_j_1 = cslope_p_pos_j_1;
	theCopy->cslope_pc_pos_j_1 = cslope_pc_pos_j_1;
	theCopy->cMpe_pos_j_1 = cMpe_pos_j_1;
	theCopy->cMpeProject_pos_j_1 = cMpeProject_pos_j_1;
	theCopy->cMmax_pos_j_1 = cMmax_pos_j_1;
	theCopy->cMmaxProject_pos_j_1 = cMmaxProject_pos_j_1;

	theCopy->cTheta_y_neg_j_1 = cTheta_y_neg_j_1;
	theCopy->cTheta_max_neg_j_1 = cTheta_max_neg_j_1;
	theCopy->cslope_p_neg_j_1 = cslope_p_neg_j_1;
	theCopy->cslope_pc_neg_j_1 = cslope_pc_neg_j_1;
	theCopy->cMpe_neg_j_1 = cMpe_neg_j_1;
	theCopy->cMpeProject_neg_j_1 = cMpeProject_neg_j_1;
	theCopy->cMmax_neg_j_1 = cMmax_neg_j_1;
	theCopy->cMmaxProject_neg_j_1 = cMmaxProject_neg_j_1;

	theCopy->cbeta_S_j_1 = cbeta_S_j_1;
	theCopy->cbeta_C_j_1 = cbeta_C_j_1;
	theCopy->cbeta_K_j_1 = cbeta_K_j_1;
	theCopy->cbeta_F_j_1 = cbeta_F_j_1;

	theCopy->cExcursion_Flag = cExcursion_Flag;
	theCopy->cReversal_Flag = cReversal_Flag;
	theCopy->cYield_Flag = cYield_Flag;
	theCopy->cFail_FlagPos = cFail_FlagPos;
	theCopy->cFail_FlagNeg = cFail_FlagNeg;
	theCopy->cMrpos_Flag = cMrpos_Flag;
	theCopy->cMrneg_Flag = cMrneg_Flag;
	theCopy->cEnergy_Flag = cEnergy_Flag;

	theCopy->cEnergy_Excrsni_1 = cEnergy_Excrsni_1;
	theCopy->cEnergy_Excrsn = cEnergy_Excrsn;
	theCopy->cEnergy_Rev = cEnergy_Rev;
	theCopy->cEnergy_total = cEnergy_total;

	return theCopy;
}

int IMKBilin::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	cout << " sendSelf" << endln;

	static Vector data(120);
	data(0) = this->getTag();
	data(1) = Ke;
	data(2) = Theta_p_pos0;
	data(3) = Theta_pc_pos0;
	data(4) = Theta_u_pos0;
	data(5) = Mpe_pos0;
	data(6) = MmaxMpe_pos0;
	data(7) = ResM_pos0;
	data(8) = Theta_p_neg0;
	data(9) = Theta_pc_neg0;
	data(10) = Theta_u_neg0;
	data(11) = Mpe_neg0;
	data(12) = MmaxMpe_neg0;
	data(13) = ResM_neg0;
	data(14) = LAMBDA_S;
	data(15) = LAMBDA_C;
	data(16) = LAMBDA_K;
	data(17) = LAMBDA_F;
	data(18) = c_S;
	data(19) = c_C;
	data(20) = c_K;
	data(21) = c_F;
	data(22) = D_pos;
	data(23) = D_neg;

	data(24) = Ri;
	data(25) = Mi;
	data(26) = Di;
	data(27) = Ri_1;
	data(28) = Mi_1;
	data(29) = Di_1;
	data(30) = Rreversal;
	data(31) = Mreversal;
	data(32) = TangentK;

	data(33) = K_j_1;
	data(34) = Theta_y_pos_j_1;
	data(35) = Theta_max_pos_j_1;
	data(36) = slope_p_pos_j_1;
	data(37) = slope_pc_pos_j_1;
	data(38) = Mpe_pos_j_1;
	data(39) = MpeProject_pos_j_1;
	data(40) = Mmax_pos_j_1;
	data(41) = MmaxProject_pos_j_1;

	data(42) = Theta_y_neg_j_1;
	data(43) = Theta_max_neg_j_1;
	data(44) = slope_p_neg_j_1;
	data(45) = slope_pc_neg_j_1;
	data(46) = Mpe_neg_j_1;
	data(47) = MpeProject_neg_j_1;
	data(48) = Mmax_neg_j_1;
	data(49) = MmaxProject_neg_j_1;

	data(50) = beta_S_j_1;
	data(51) = beta_C_j_1;
	data(52) = beta_K_j_1;
	data(53) = beta_F_j_1;

	data(54) = Ref_Energy_S;
	data(55) = Ref_Energy_C;
	data(56) = Ref_Energy_K;
	data(57) = Ref_Energy_F;

	data(58) = Excursion_Flag;
	data(59) = Reversal_Flag;
	data(60) = Yield_Flag;
	data(61) = Fail_FlagPos;
	data(62) = Fail_FlagNeg;
	data(63) = Mrpos_Flag;
	data(64) = Mrneg_Flag;
	data(65) = Energy_Flag;

	data(66) = Energy_Excrsni_1;
	data(67) = Energy_Excrsn;
	data(68) = Energy_Rev;
	data(69) = Energy_total;

	data(70) = cRi;
	data(71) = cMi;
	data(72) = cDi;
	data(73) = cRi_1;
	data(74) = cMi_1;
	data(75) = cDi_1;
	data(76) = cRreversal;
	data(77) = cMreversal;
	data(78) = cTangentK;

	data(79) = cK_j_1;
	data(80) = cTheta_y_pos_j_1;
	data(81) = cTheta_max_pos_j_1;
	data(82) = cslope_p_pos_j_1;
	data(83) = cslope_pc_pos_j_1;
	data(84) = cMpe_pos_j_1;
	data(85) = cMpeProject_pos_j_1;
	data(86) = cMmax_pos_j_1;
	data(87) = cMmaxProject_pos_j_1;

	data(88) = cTheta_y_neg_j_1;
	data(89) = cTheta_max_neg_j_1;
	data(90) = cslope_p_neg_j_1;
	data(91) = cslope_pc_neg_j_1;
	data(92) = cMpe_neg_j_1;
	data(93) = cMpeProject_neg_j_1;
	data(94) = cMmax_neg_j_1;
	data(95) = cMmaxProject_neg_j_1;

	data(96) = cbeta_S_j_1;
	data(97) = cbeta_C_j_1;
	data(98) = cbeta_K_j_1;
	data(99) = cbeta_F_j_1;

	data(100) = cExcursion_Flag;
	data(101) = cReversal_Flag;
	data(102) = cYield_Flag;
	data(103) = cFail_FlagPos;
	data(104) = cFail_FlagNeg;
	data(105) = cMrpos_Flag;
	data(106) = cMrneg_Flag;
	data(107) = cEnergy_Flag;

	data(108) = cEnergy_Excrsni_1;
	data(109) = cEnergy_Excrsn;
	data(110) = cEnergy_Rev;
	data(111) = cEnergy_total;

	data(112) = Mr_pos0;
	data(113) = Mr_neg0;
	data(114) = cMr_pos0;
	data(115) = cMr_neg0;

	data(116) = U;
	data(117) = cU;

	data(118) = n;
	data(119) = Roffset;

	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "IMKBilin::sendSelf() - failed to send data\n";

	return res;
}

int IMKBilin::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector data(120);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);

	if (res < 0) {
		opserr << "IMKBilin::recvSelf() - failed to receive data\n";
		this->setTag(0);
	}
	else {
		cout << " recvSelf" << endln;
		this->setTag((int)data(0));
		Ke = data(1);
		Theta_p_pos0 = data(2);
		Theta_pc_pos0 = data(3);
		Theta_u_pos0 = data(4);
		Mpe_pos0 = data(5);
		MmaxMpe_pos0 = data(6);
		ResM_pos0 = data(7);
		Theta_p_neg0 = data(8);
		Theta_pc_neg0 = data(9);
		Theta_u_neg0 = data(10);
		Mpe_neg0 = data(11);
		MmaxMpe_neg0 = data(12);
		ResM_neg0 = data(13);
		LAMBDA_S = data(14);
		LAMBDA_C = data(15);
		LAMBDA_K = data(16);
		LAMBDA_F = data(17);
		c_S = data(18);
		c_C = data(19);
		c_K = data(20);
		c_F = data(21);
		D_pos = data(22);
		D_neg = data(23);

		Ri = data(24);
		Mi = data(25);
		Di = data(26);
		Ri_1 = data(27);
		Mi_1 = data(28);
		Di_1 = data(29);
		Rreversal = data(30);
		Mreversal = data(31);
		TangentK = data(32);

		K_j_1 = data(33);
		Theta_y_pos_j_1 = data(34);
		Theta_max_pos_j_1 = data(35);
		slope_p_pos_j_1 = data(36);
		slope_pc_pos_j_1 = data(37);
		Mpe_pos_j_1 = data(38);
		MpeProject_pos_j_1 = data(39);
		Mmax_pos_j_1 = data(40);
		MmaxProject_pos_j_1 = data(41);

		Theta_y_neg_j_1 = data(42);
		Theta_max_neg_j_1 = data(43);
		slope_p_neg_j_1 = data(44);
		slope_pc_neg_j_1 = data(45);
		Mpe_neg_j_1 = data(46);
		MpeProject_neg_j_1 = data(47);
		Mmax_neg_j_1 = data(48);
		MmaxProject_neg_j_1 = data(49);

		beta_S_j_1 = data(50);
		beta_C_j_1 = data(51);
		beta_K_j_1 = data(52);
		beta_F_j_1 = data(53);

		Ref_Energy_S = data(54);
		Ref_Energy_C = data(55);
		Ref_Energy_K = data(56);
		Ref_Energy_F = data(57);

		Excursion_Flag = data(58);
		Reversal_Flag = data(59);
		Yield_Flag = data(60);
		Fail_FlagPos = data(61);
		Fail_FlagNeg = data(62);
		Mrpos_Flag = data(63);
		Mrneg_Flag = data(64);
		Energy_Flag = data(65);

		Energy_Excrsni_1 = data(66);
		Energy_Excrsn = data(67);
		Energy_Rev = data(68);
		Energy_total = data(69);

		cRi = data(70);
		cMi = data(71);
		cDi = data(72);
		cRi_1 = data(73);
		cMi_1 = data(74);
		cDi_1 = data(75);

		cRreversal = data(76);
		cMreversal = data(77);
		cTangentK = data(78);

		cK_j_1 = data(79);
		cTheta_y_pos_j_1 = data(80);
		cTheta_max_pos_j_1 = data(81);
		cslope_p_pos_j_1 = data(82);
		cslope_pc_pos_j_1 = data(83);
		cMpe_pos_j_1 = data(84);
		cMpeProject_pos_j_1 = data(85);
		cMmax_pos_j_1 = data(86);
		cMmaxProject_pos_j_1 = data(87);

		cTheta_y_neg_j_1 = data(88);
		cTheta_max_neg_j_1 = data(89);
		cslope_p_neg_j_1 = data(90);
		cslope_pc_neg_j_1 = data(91);
		cMpe_neg_j_1 = data(92);
		cMpeProject_neg_j_1 = data(93);
		cMmax_neg_j_1 = data(94);
		cMmaxProject_neg_j_1 = data(95);

		cbeta_S_j_1 = data(96);
		cbeta_C_j_1 = data(97);
		cbeta_K_j_1 = data(98);
		cbeta_F_j_1 = data(99);

		cExcursion_Flag = data(100);
		cReversal_Flag = data(101);
		cYield_Flag = data(102);
		cFail_FlagPos = data(103);
		cFail_FlagNeg = data(104);
		cMrpos_Flag = data(105);
		cMrneg_Flag = data(106);
		cEnergy_Flag = data(107);

		cEnergy_Excrsni_1 = data(108);
		cEnergy_Excrsn = data(109);
		cEnergy_Rev = data(110);
		cEnergy_total = data(111);

		Mr_pos0 = data(112);
		Mr_neg0 = data(113);
		cMr_pos0 = data(114);
		cMr_neg0 = data(115);

		U = data(116);
		cU = data(117);

		n = data(118);
		Roffset = data(119);
	}

	return res;
}

void IMKBilin::Print(OPS_Stream &s, int flag)
{
	cout << "IMKBilin tag: " << this->getTag() << endln;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%% SPLINE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/* 
// ***************************************************************************
double IMKBilin::r8_max(double x, double y)
// ****************************************************************************

{
	if (y < x)
	{
		return x;
	}
	else
	{
		return y;
	}
}

// ***************************************************************************
double IMKBilin::r8_min(double x, double y)
// ****************************************************************************

{
	if (y < x)
	{
		return y;
	}
	else
	{
		return x;
	}
}

// ***************************************************************************
int IMKBilin::i4_max(int i1, int i2)
// ****************************************************************************

{
	int value;

	if (i2 < i1)
	{
		value = i1;
	}
	else
	{
		value = i2;
	}
	return value;
}

// ***************************************************************************
double IMKBilin::pchst(double arg1, double arg2)
// ****************************************************************************

{
	double value;

	if (arg1 == 0.0)
	{
		value = 0.0;
	}
	else if (arg1 < 0.0)
	{
		if (arg2 < 0.0)
		{
			value = 1.0;
		}
		else if (arg2 == 0.0)
		{
			value = 0.0;
		}
		else if (0.0 < arg2)
		{
			value = -1.0;
		}
	}
	else if (0.0 < arg1)
	{
		if (arg2 < 0.0)
		{
			value = -1.0;
		}
		else if (arg2 == 0.0)
		{
			value = 0.0;
		}
		else if (0.0 < arg2)
		{
			value = 1.0;
		}
	}

	return value;
}

// **************************************************************************
int IMKBilin::chfev(double x1, double x2, double f1, double f2, double d1, double d2,
	int ne, double xe[], double fe[], int next[])
// ****************************************************************************

{
	double c2;
	double c3;
	double del1;
	double del2;
	double delta;
	double h;
	int i;
	int ierr;
	double x;
	double xma;
	double xmi;

	if (ne < 1)
	{
		ierr = -1;
		cerr << "\n";
		cerr << "CHFEV - Fatal error!\n";
		cerr << "  Number of evaluation points is less than 1.\n";
		cerr << "  NE = " << ne << "\n";
		return ierr;
	}

	h = x2 - x1;

	if (h == 0.0)
	{
		ierr = -2;
		cerr << "\n";
		cerr << "CHFEV - Fatal error!\n";
		cerr << "  The interval [X1,X2] is of zero length.\n";
		return ierr;
	}
	//
	//  Initialize.
	//
	ierr = 0;
	next[0] = 0;
	next[1] = 0;
	xmi = r8_min(0.0, h);
	xma = r8_max(0.0, h);
	//
	//  Compute cubic coefficients expanded about X1.
	//
	delta = (f2 - f1) / h;
	del1 = (d1 - delta) / h;
	del2 = (d2 - delta) / h;
	c2 = -(del1 + del1 + del2);
	c3 = (del1 + del2) / h;
	//
	//  Evaluation loop.
	//
	for (i = 0; i < ne; i++)
	{
		x = xe[i] - x1;
		fe[i] = f1 + x * (d1 + x * (c2 + x * c3));
		//
		//  Count the extrapolation points.
		//
		if (x < xmi)
		{
			next[0] = next[0] + 1;
		}

		if (xma < x)
		{
			next[1] = next[1] + 1;
		}

	}

	return 0;
}

// ***************************************************************************
void IMKBilin::spline_pchip_set(int n, double x[], double f[], double d[])
// ****************************************************************************80

{
	double del1;
	double del2;
	double dmax;
	double dmin;
	double drat1;
	double drat2;
	double dsave;
	double h1;
	double h2;
	double hsum;
	double hsumt3;
	int i;
	int ierr;
	int nless1;
	double temp;
	double w1;
	double w2;
	//
	//  Check the arguments.
	//
	if (n < 2)
	{
		ierr = -1;
		cerr << "\n";
		cerr << "SPLINE_PCHIP_SET - Fatal error!\n";
		cerr << "  Number of data points less than 2.\n";
		exit(ierr);
	}

	for (i = 1; i < n; i++)
	{
		if (x[i] <= x[i - 1])
		{
			ierr = -3;
			cerr << "\n";
			cerr << "SPLINE_PCHIP_SET - Fatal error!\n";
			cerr << "  X array not strictly increasing.\n";
			exit(ierr);
		}
	}

	ierr = 0;
	nless1 = n - 1;
	h1 = x[1] - x[0];
	del1 = (f[1] - f[0]) / h1;
	dsave = del1;
	//
	//  Special case N=2, use linear interpolation.
	//
	if (n == 2)
	{
		d[0] = del1;
		d[n - 1] = del1;
		return;
	}
	//
	//  Normal case, 3 <= N.
	//
	h2 = x[2] - x[1];
	del2 = (f[2] - f[1]) / h2;
	//
	//  Set D(1) via non-centered three point formula, adjusted to be
	//  shape preserving.
	//
	hsum = h1 + h2;
	w1 = (h1 + hsum) / hsum;
	w2 = -h1 / hsum;
	d[0] = w1 * del1 + w2 * del2;

	if (pchst(d[0], del1) <= 0.0)
	{
		d[0] = 0.0;
	}
	//
	//  Need do this check only if monotonicity switches.
	//
	else if (pchst(del1, del2) < 0.0)
	{
		dmax = 3.0 * del1;

		if (fabs(dmax) < fabs(d[0]))
		{
			d[0] = dmax;
		}

	}
	//
	//  Loop through interior points.
	//
	for (i = 2; i <= nless1; i++)
	{
		if (2 < i)
		{
			h1 = h2;
			h2 = x[i] - x[i - 1];
			hsum = h1 + h2;
			del1 = del2;
			del2 = (f[i] - f[i - 1]) / h2;
		}
		//
		//  Set D(I)=0 unless data are strictly monotonic.
		//
		d[i - 1] = 0.0;

		temp = pchst(del1, del2);

		if (temp < 0.0)
		{
			ierr = ierr + 1;
			dsave = del2;
		}
		//
		//  Count number of changes in direction of monotonicity.
		//
		else if (temp == 0.0)
		{
			if (del2 != 0.0)
			{
				if (pchst(dsave, del2) < 0.0)
				{
					ierr = ierr + 1;
				}
				dsave = del2;
			}
		}
		//
		//  Use Brodlie modification of Butland formula.
		//
		else
		{
			hsumt3 = 3.0 * hsum;
			w1 = (hsum + h1) / hsumt3;
			w2 = (hsum + h2) / hsumt3;
			dmax = r8_max(fabs(del1), fabs(del2));
			dmin = r8_min(fabs(del1), fabs(del2));
			drat1 = del1 / dmax;
			drat2 = del2 / dmax;
			d[i - 1] = dmin / (w1 * drat1 + w2 * drat2);
		}
	}
	//
	//  Set D(N) via non-centered three point formula, adjusted to be
	//  shape preserving.
	//
	w1 = -h2 / hsum;
	w2 = (h2 + hsum) / hsum;
	d[n - 1] = w1 * del1 + w2 * del2;

	if (pchst(d[n - 1], del2) <= 0.0)
	{
		d[n - 1] = 0.0;
	}
	else if (pchst(del1, del2) < 0.0)
	{
		//
		//  Need do this check only if monotonicity switches.
		//
		dmax = 3.0 * del2;

		if (fabs(dmax) < fabs(d[n - 1]))
		{
			d[n - 1] = dmax;
		}

	}
	return;
}

// ***************************************************************************
void IMKBilin::spline_pchip_val(int n, double x[], double f[], double d[],
	int ne, double xe[], double fe[])
// ****************************************************************************

{
	int i;
	int ierc;
	int ierr;
	int ir;
	int j;
	int j_first;
	int j_new;
	int j_save;
	int next[2];
	int nj;
	//
	//  Check arguments.
	//
	if (n < 2)
	{
		ierr = -1;
		cerr << "\n";
		cerr << "SPLINE_PCHIP_VAL - Fatal error!\n";
		cerr << "  Number of data points less than 2.\n";
		exit(ierr);
	}

	for (i = 1; i < n; i++)
	{
		if (x[i] <= x[i - 1])
		{
			ierr = -3;
			cerr << "\n";
			cerr << "SPLINE_PCHIP_VAL - Fatal error!\n";
			cerr << "  X array not strictly increasing.\n";
			exit(ierr);
		}
	}

	if (ne < 1)
	{
		ierr = -4;
		cerr << "\n";
		cerr << "SPLINE_PCHIP_VAL - Fatal error!\n";
		cerr << "  Number of evaluation points less than 1.\n";
		return;
	}

	ierr = 0;
	//
	//  Loop over intervals.
	//  The interval index is IL = IR-1.
	//  The interval is X(IL) <= X < X(IR).
	//
	j_first = 1;
	ir = 2;

	for (; ; )
	{
		//
		//  Skip out of the loop if have processed all evaluation points.
		//
		if (ne < j_first)
		{
			break;
		}
		//
		//  Locate all points in the interval.
		//
		j_save = ne + 1;

		for (j = j_first; j <= ne; j++)
		{
			if (x[ir - 1] <= xe[j - 1])
			{
				j_save = j;
				if (ir == n)
				{
					j_save = ne + 1;
				}
				break;
			}
		}
		//
		//  Have located first point beyond interval.
		//
		j = j_save;

		nj = j - j_first;
		//
		//  Skip evaluation if no points in interval.
		//
		if (nj != 0)
		{
			//
			//  Evaluate cubic at XE(J_FIRST:J-1).
			//
			ierc = chfev(x[ir - 2], x[ir - 1], f[ir - 2], f[ir - 1], d[ir - 2], d[ir - 1],
				nj, xe + j_first - 1, fe + j_first - 1, next);

			if (ierc < 0)
			{
				ierr = -5;
				cerr << "\n";
				cerr << "SPLINE_PCHIP_VAL - Fatal error!\n";
				cerr << "  Error return from CHFEV.\n";
				exit(ierr);
			}
			//
			//  In the current set of XE points, there are NEXT(2) to the right of X(IR).
			//
			if (next[1] != 0)
			{
				if (ir < n)
				{
					ierr = -5;
					cerr << "\n";
					cerr << "SPLINE_PCHIP_VAL - Fatal error!\n";
					cerr << "  IR < N.\n";
					exit(ierr);
				}
				//
				//  These are actually extrapolation points.
				//
				ierr = ierr + next[1];

			}
			//
			//  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
			//
			if (next[0] != 0)
			{
				//
				//  These are actually extrapolation points.
				//
				if (ir <= 2)
				{
					ierr = ierr + next[0];
				}
				else
				{
					j_new = -1;

					for (i = j_first; i <= j - 1; i++)
					{
						if (xe[i - 1] < x[ir - 2])
						{
							j_new = i;
							break;
						}
					}

					if (j_new == -1)
					{
						ierr = -5;
						cerr << "\n";
						cerr << "SPLINE_PCHIP_VAL - Fatal error!\n";
						cerr << "  Could not bracket the data point.\n";
						exit(ierr);
					}
					//
					//  Reset J.  This will be the new J_FIRST.
					//
					j = j_new;
					//
					//  Now find out how far to back up in the X array.
					//
					for (i = 1; i <= ir - 1; i++)
					{
						if (xe[j - 1] < x[i - 1])
						{
							break;
						}
					}
					//
					//  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .
					//
					//  Reset IR, recognizing that it will be incremented before cycling.
					//
					ir = i4_max(1, i - 1);
				}
			}

			j_first = j;
		}

		ir = ir + 1;

		if (n < ir)
		{
			break;
		}

	}

	return;
}
// *************************************************************************** */
