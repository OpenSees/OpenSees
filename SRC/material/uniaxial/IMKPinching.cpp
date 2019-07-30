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
#include <IMKPinching.h>
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

static int numIMKPinchingMaterials = 0;

void *
OPS_IMKPinching(void)
{
	if (numIMKPinchingMaterials == 0) {
		numIMKPinchingMaterials++;
		OPS_Error("\nIMK Model with Pinched Response - Code by H. ELJISR & A. ELKADY (July-2018)\n", 1);
	}

	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial = 0;

	int    iData[1];
	double dData[25];
	int numData = 1;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial IMKPinching tag" << endln;
		return 0;
	}

	numData = 25;

	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid Args want: uniaxialMaterial IMKPinching tag? Ke? ";
		opserr << "Up_pos? Upc_pos? Uu_pos? Fy_pos? FmaxFy_pos? ResF_pos? ";
		opserr << "Up_neg? Upc_neg? Uu_neg? Fy_neg? FmaxFy_neg? ResF_neg? ";
		opserr << "LamdaS? LamdaC? LamdaA? LamdaK? Cs? Cc? Ca? Ck? D_pos? D_neg? ";
		opserr << "KappaF? KappaD?";
		return 0;
	}



	// Parsing was successful, allocate the material
	theMaterial = new IMKPinching(iData[0],
		dData[0],
		dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
		dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
		dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20],
		dData[21], dData[22], dData[23], dData[24]);

	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type IMKPinching Material\n";
		return 0;
	}

	return theMaterial;
}

IMKPinching::IMKPinching(int tag, double p_Ke,
	double p_Up_pos, double p_Upc_pos, double p_Uu_pos, double p_Fy_pos, double p_FmaxFy_pos, double p_ResF_pos,
	double p_Up_neg, double p_Upc_neg, double p_Uu_neg, double p_Fy_neg, double p_FmaxFy_neg, double p_ResF_neg,
	double p_LAMBDA_S, double p_LAMBDA_C, double p_LAMBDA_A, double p_LAMBDA_K, double p_c_S, double p_c_C, double p_c_A, double p_c_K, double p_D_pos, double p_D_neg, double p_kappaF, double p_kappaD)
	: UniaxialMaterial(tag, 0), Ke(p_Ke),
	Up_pos(p_Up_pos), Upc_pos(p_Upc_pos), Uu_pos(p_Uu_pos), Fy_pos(p_Fy_pos), FmaxFy_pos(p_FmaxFy_pos), ResF_pos(p_ResF_pos),
	Up_neg(p_Up_neg), Upc_neg(p_Upc_neg), Uu_neg(p_Uu_neg), Fy_neg(p_Fy_neg), FmaxFy_neg(p_FmaxFy_neg), ResF_neg(p_ResF_neg),
	LAMBDA_S(p_LAMBDA_S), LAMBDA_C(p_LAMBDA_C), LAMBDA_A(p_LAMBDA_A), LAMBDA_K(p_LAMBDA_K), c_S(p_c_S), c_C(p_c_C), c_A(p_c_A), c_K(p_c_K), D_pos(p_D_pos), D_neg(p_D_neg), kappaF(p_kappaF), kappaD(p_kappaD)
{
	this->revertToStart();
}

IMKPinching::IMKPinching()
	:UniaxialMaterial(0, 0), Ke(0),
	Up_pos(0), Upc_pos(0), Uu_pos(0), Fy_pos(0), FmaxFy_pos(0), ResF_pos(0),
	Up_neg(0), Upc_neg(0), Uu_neg(0), Fy_neg(0), FmaxFy_neg(0), ResF_neg(0),
	LAMBDA_S(0), LAMBDA_C(0), LAMBDA_A(0), LAMBDA_K(0), c_S(0), c_C(0), c_A(0), c_K(0), D_pos(0), D_neg(0), kappaF(0), kappaD(0)
{
	this->revertToStart();
}

IMKPinching::~IMKPinching()
{
	// does nothing
}

int IMKPinching::setTrialStrain(double strain, double strainRate)
{
	//all variables to the last commit
	this->revertToLastCommit();

	//state determination algorithm: defines the current force and tangent stiffness
	U = strain; //set trial displacement
	ui_1 = ui;
	fi_1 = fi;
	ui = U;

	//cout << "***********************" << endln;
	//cout << "  Excurion=" << Excursion_Flag << " Failure=" << Failure_Flag << endln;
	//cout << "  STEP: ui_1=" << ui_1 << " ui=" << ui << " fi_1=" << fi_1 << " fi=" << fi << endln;
	//cout << "  +VE: Uy= " << Uy_pos_j_1 << " Umax= " << Umax_pos_j_1 << " Upeak= " << Upeak_pos_j_1  << " Ubp=" << Ubp_pos_j_1 << " Fpeak= " << Fpeak_pos_j_1 << " Fbp=" << Fbp_pos_j_1 << " KrelA=" << KrelA_pos_j_1 << " KrelB=" << KrelB_pos_j_1  << endln;
	//cout << "  -VE: Uy=" << Uy_neg_j_1 << " Umax=" << Umax_neg_j_1 << " Upeak=" << Upeak_neg_j_1  << " Ubp=" << Ubp_neg_j_1 << " Fpeak=" << Fpeak_neg_j_1 << " Fbp=" << Fbp_neg_j_1 << " KrelA=" << KrelA_neg_j_1 << " KrelB=" << KrelB_neg_j_1  << endln;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	// Incremental deformation at current step
	du = ui - ui_1;
	// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (Failure_Flag != 1) {
		//// Positive Loading
		if (fi_1 >= 0) {
			// Early reloading before new excursion occurs
			if ((du > 0) && (du_i_1 < 0) && (fi_1 > 0)) {
				//cout << "  OK6**" << endln;

				// Reloading stiffness KrelA if reloading is to the left of the break point
				if (ui <= Ubp_pos_j_1) {
					KrelA_pos_j_1 = (Fbp_pos_j_1 - fi_1) / (Ubp_pos_j_1 - ui_1);
					// Reloading stiffness KrelA if reloading is to the right of the break point
				}
				else {
					KrelB_pos_j_1 = (Fpeak_pos_j_1 - fi_1) / (Upeak_pos_j_1 - ui_1);
				}
			}
			// Deformation in first reloading branch KrelA
			if ((ui <= Ubp_pos_j_1) && (du > 0)) {
				//cout << "  OK5**" << endln;
				df = KrelA_pos_j_1*du;
				Umaxp = ui;
				// Deformation in second reloading branch KrelB
			}
			else if ((ui <= Upeak_pos_j_1) && (du > 0)) {
				//cout << "  OK4**" << endln;
				df = (KrelB_pos_j_1*ui + (Fpeak_pos_j_1 - KrelB_pos_j_1*Upeak_pos_j_1)) - fi_1;
				Umaxp = ui;
				// Deformation in post-yield branch of the backbone
			}
			else if ((ui <= Umax_pos_j_1) && (du > 0)) {
				//cout << "  OK3**" << endln;
				df = (Fy_pos_j_1 + Kp_pos_j_1*(ui - Uy_pos_j_1)) - fi_1;
				Umaxp = ui;
				// Deformation in the post-capping branch of the backbone
			}
			else if ((ui >= Umax_pos_j_1) && (du > 0)) {
				//cout << "  OK2**" << endln;

				// Deformation in residual branch of backbone
				if (ui > Ures_pos_j_1) {
					df = Fres_pos_j_1 - fi_1;
					if (Fres_pos_j_1 == 0) {
						//cout << "Res Fail" << endln;
						Failure_Flag = 1;
					}
					// Deformation in softening branch of the backbone
				}
				else {
					df = (Fmax_pos_j_1 + Kpc_pos_j_1*(ui - Umax_pos_j_1)) - fi_1;
				}
				Umaxp = ui;
				// Deformation in the unloading branch
			}
			else {
				df = Kul_j_1*du;
				//cout << "  OK1**" << endln;

			}
			// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
				//// Negative Loading
		}
		else {
			// Early reloading before new excursion occurs
			if ((du < 0) && (du_i_1 > 0) && (fi_1 < 0)) {
				// Reloading stiffness KrelA if reloading is to the right of the break point
				if (ui >= Ubp_neg_j_1) {
					KrelA_neg_j_1 = (Fbp_neg_j_1 - fi_1) / (Ubp_neg_j_1 - ui_1);
					// Reloading stiffness KrelA if reloading is to the left of the break point
				}
				else {
					KrelB_neg_j_1 = (Fpeak_neg_j_1 - fi_1) / (Upeak_neg_j_1 - ui_1);
				}
			}
			// Deformation in first reloading branch KrelA
			if ((ui >= Ubp_neg_j_1) && (du < 0)) {
				df = KrelA_neg_j_1*du;
				Umaxn = ui;
				// Deformation in second reloading branch KrelB
			}
			else if ((ui >= Upeak_neg_j_1) && (du < 0)) {
				df = (KrelB_neg_j_1*ui + (Fpeak_neg_j_1 - KrelB_neg_j_1*Upeak_neg_j_1)) - fi_1;
				Umaxn = ui;
				// Deformation in post-yield branch of the backbone
			}
			else if ((ui >= Umax_neg_j_1) && (du < 0)) {
				df = (Fy_neg_j_1 + Kp_neg_j_1*(ui - Uy_neg_j_1)) - fi_1;
				Umaxn = ui;
				// Deformationin the post-capping branch of the backbone
			}
			else if ((ui <= Umax_neg_j_1) && (du < 0)) {
				// Deformation in residual branch of backbone
				if (ui < Ures_neg_j_1) {
					df = Fres_neg_j_1 - fi_1;
					if (Fres_neg_j_1 == 0) {
						//cout << "  Res Fail" << endln;
						Failure_Flag = 1;
					}
					// Deformation in the softening branch of the backbone
				}
				else {
					df = Fmax_neg_j_1 + Kpc_neg_j_1*(ui - Umax_neg_j_1) - fi_1;
				}
				Umaxn = ui;
			}
			else {
				// Deformation in the unloading branch
				df = Kul_j_1 *du;
			}
		}
		// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
			//// Deterioration Parameters
			// Internal energy increment
		dEi = 0.5*(df + 2 * fi_1)*(du);
		//cout << "  ENERGY: dEi=" << dEi << " Kul=" << Kul_j_1 << " du=" << du << " df=" << df << endln;

		// Positive excursion flag
		if ((fi_1 + df >= 0) && (fi_1 < 0)) {
			Excursion_Flag = 1;
			// Negative excursion flag
		}
		else if ((fi_1 + df <= 0) && (fi_1 > 0)) {
			Excursion_Flag = 1;
		}
		else {
			Excursion_Flag = 0;
		}
		// Update beta parameters at new excursion
		if (Excursion_Flag == 1) {
			// Total energy dissipated in all previous excursions
			Epj = Energy_Acc + dEi;
			// Energy dissipated in previous excursion
			Ei = Epj - Energy_Diss;
			betaS = pow((Ei / (EtS - Epj)), c_S);
			betaC = pow((Ei / (EtC - Epj)), c_C);
			betaA = pow((Ei / (EtA - Epj)), c_A);
		}
		else {
			// Total energy dissipated in all previous excursions
			Epj = Energy_Diss;
			betaS = 0;
			betaC = 0;
			betaA = 0;
		}
		// Onset of unloading
		Unloading_Flag = du*du_i_1 < 0 && (ui_1 >= Upeak_pos_j_1 || ui_1 <= Upeak_neg_j_1);
		if (Unloading_Flag == 1) {
			// Total energy dissipated until point of unloading
			EpjK = dEi + Energy_Acc - 0.5*(pow((fi_1 + df), 2)) / Kul_j_1;
			// Energy dissipated in current excursion until point of unloading
			EiK = EpjK - Energy_Diss;
			betaK = pow((EiK / (EtK - EpjK)), c_K);
		}
		else {
			betaK = 0;
		}
		// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
			//// Target Peak Deformation
		New_Peak_Pos_Flag = 0;
		New_Peak_Neg_Flag = 0;
		// Update target peak deformation for positive loading
		if (Umaxp >= Upeak_pos_j_1) {
			New_Peak_Pos_Flag = 1;
			Upeak_pos_j_1 = Umaxp;
			Fpeak_pos_j_1 = fi_1 + df;
			// Plastic offset for positive loading
			Plastic_Offset_pos_j_1 = Upeak_pos_j_1 - Fpeak_pos_j_1 / Kul_j_1;
		}
		// Update target peak deformation for negative loading
		if (Umaxn <= Upeak_neg_j_1) {
			New_Peak_Neg_Flag = 1;
			Upeak_neg_j_1 = Umaxn;
			Fpeak_neg_j_1 = fi_1 + df;
			// Plastic offset for negative loading
			Plastic_Offset_neg_j_1 = Upeak_neg_j_1 - Fpeak_neg_j_1 / Kul_j_1;
		}
		// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
			//// Update Positive Backbone and Target Peak Point
		if (Excursion_Flag == 1) {
			// Positive loading backbone
			if (fi_1 < 0) {
				// Basic strength deterioration: Yield point
				Uy_pos_j_1 = std::max(Uy_pos_j_1 - Fy_pos_j_1 *betaS* D_pos / Ke, Fres_pos_j_1 / Ke);
				Fy_pos_j_1 = std::max(Fy_pos_j_1 *(1 - betaS* D_pos), Fres_pos_j_1);
				// Basic strength deterioration: Post-yield Stiffness
				if (Fy_pos_j_1 != Fres_pos_j_1) {
					Kp_pos_j_1 = Kp_pos_j_1 *(1 - betaS* D_pos);
				}
				else {
					Kp_pos_j_1 = 0;
				}
				// Basic strength deterioration: Capping Point
				sPCsp = (Fy_pos_j_1 - Uy_pos_j_1 *Kp_pos_j_1 - Fmax_pos_j_1 + Kpc_pos_j_1*Umax_pos_j_1) / (Kpc_pos_j_1 - Kp_pos_j_1);
				Fmax_pos_j_1 = Fmax_pos_j_1 + (sPCsp - Umax_pos_j_1)*Kpc_pos_j_1;
				Umax_pos_j_1 = sPCsp;
				// Post-capping strength deterioration: Capping point
				sPCpcp = max(Umax_pos_j_1 + betaC* D_pos*(Fmax_pos_j_1 - Kpc_pos_j_1*Umax_pos_j_1) / (Kpc_pos_j_1 - Kp_pos_j_1), Uy_pos_j_1);
				Fmax_pos_j_1 = Fmax_pos_j_1 + (sPCpcp - Umax_pos_j_1)*Kp_pos_j_1;
				Umax_pos_j_1 = sPCpcp;
				// Accelerated reloading stiffness deterioration: Target peak deformation point
				Upeak_pos_j_1 = (1 + betaA* D_pos)*Upeak_pos_j_1;
				// Target peak deformation in reloading branch of the updated backbone
				if (Upeak_pos_j_1 <= Uy_pos_j_1) {
					Fpeak_pos_j_1 = Ke*Upeak_pos_j_1;
					// Target peak deformation in post-yield branch of the updated backbone
				}
				else if (Upeak_pos_j_1 <= Umax_pos_j_1) {
					Fpeak_pos_j_1 = Kp_pos_j_1 *(Upeak_pos_j_1 - Uy_pos_j_1) + Fy_pos_j_1;
					// Target peak deformation in post-capping branch of the updated backbone
				}
				else {
					Fpeak_pos_j_1 = max(Kpc_pos_j_1*(Upeak_pos_j_1 - Umax_pos_j_1) + Fmax_pos_j_1, Fres_pos_j_1);
				}
				// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
						//// Update Negative Backbone and Target Peak Point
			}
			else {
				// Basic strength deterioration: Yield point
				Uy_neg_j_1 = min(Uy_neg_j_1 - Fy_neg_j_1 *betaS* D_neg / Ke, Fres_neg_j_1 / Ke);
				Fy_neg_j_1 = min(Fy_neg_j_1 *(1 - betaS* D_neg), Fres_neg_j_1);
				// Basic strength deterioration: Post-yield stiffness
				if (Fy_neg_j_1 != Fres_neg_j_1) {
					Kp_neg_j_1 = Kp_neg_j_1 *(1 - betaS* D_neg);
				}
				else {
					Kp_neg_j_1 = 0;
				}
				// Basic strength deterioration: Capping point
				sPCsn = (Fy_neg_j_1 - Uy_neg_j_1 *Kp_neg_j_1 - Fmax_neg_j_1 + Kpc_neg_j_1*Umax_neg_j_1) / (Kpc_neg_j_1 - Kp_neg_j_1);
				Fmax_neg_j_1 = Fmax_neg_j_1 + (sPCsn - Umax_neg_j_1)*Kpc_neg_j_1;
				Umax_neg_j_1 = sPCsn;
				// Post-capping strength deterioration: Capping point
				sPCpcn = min(Umax_neg_j_1 + betaC* D_neg*(Fmax_neg_j_1 - Kpc_neg_j_1*Umax_neg_j_1) / (Kpc_neg_j_1 - Kp_neg_j_1), Uy_neg_j_1);
				Fmax_neg_j_1 = Fmax_neg_j_1 + (sPCpcn - Umax_neg_j_1)*Kp_neg_j_1;
				Umax_neg_j_1 = sPCpcn;
				// Accelerated reloading stiffness deterioration: Target peak deformation point
				Upeak_neg_j_1 = (1 + betaA* D_neg)*Upeak_neg_j_1;
				// Target peak deformation in reloading branch of the updated backbone
				if (Upeak_neg_j_1 >= Uy_neg_j_1) {
					Fpeak_neg_j_1 = Ke*Upeak_neg_j_1;
					// Target peak deformation in post-yield branch of the updated backbone
				}
				else if (Upeak_neg_j_1 >= Umax_neg_j_1) {
					Fpeak_neg_j_1 = Kp_neg_j_1 *(Upeak_neg_j_1 - Uy_neg_j_1) + Fy_neg_j_1;
					// Target peak deformation in post-capping branch of the updated backbone
				}
				else {
					Fpeak_neg_j_1 = min(Kpc_neg_j_1*(Upeak_neg_j_1 - Umax_neg_j_1) + Fmax_neg_j_1, Fres_neg_j_1);
				}
			}
		}
		// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
			//// Update Krel Based on New Peak Targets at New Positive Excursion
		if ((fi_1 + df >= 0) && (fi_1 < 0)) {
			Fp = kappaF*Fpeak_pos_j_1;
			// Deformation at reloading
			u0 = ui_1 - (fi_1) / Kul_j_1;
			if (u0 < 0) {
				// Deformation at break point
				Ubp_pos_j_1 = (1 - kappaD)*Plastic_Offset_pos_j_1;
				// Force at break point
				Fbp_pos_j_1 = Fp*(Ubp_pos_j_1 - u0) / (Upeak_pos_j_1 - u0);
			}
			// Reloading is to the left of the break point
			if (u0 < Ubp_pos_j_1) {
				// Reloading stiffness KrelA after new excursion
				KrelA_pos_j_1 = (Fbp_pos_j_1) / (Ubp_pos_j_1 - u0);
				// Reloading stiffness KrelB after new excursion
				KrelB_pos_j_1 = (Fpeak_pos_j_1 - Fbp_pos_j_1) / (Upeak_pos_j_1 - Ubp_pos_j_1);
				df = ((ui - u0)*KrelA_pos_j_1) - fi_1;
				// Reloading is to the right of the break point
			}
			else {
				// Reloading stiffness after new excursion
				KrelB_pos_j_1 = (Fpeak_pos_j_1) / (Upeak_pos_j_1 - u0);
				df = ((ui - u0)*KrelB_pos_j_1) - fi_1;
			}
			// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
				//// Update Krel Based on New Peak Targets at New Negative Excursion
		}
		else if ((fi_1 + df < 0) && (fi_1 > 0)) {
			Fp = kappaF*Fpeak_neg_j_1;
			// Deformation at reloading
			u0 = ui_1 - (fi_1) / Kul_j_1;
			if (u0 > 0) {
				// Deformation at break point
				Ubp_neg_j_1 = (1 - kappaD)*Plastic_Offset_neg_j_1;
				// Force at break point
				Fbp_neg_j_1 = Fp*(Ubp_neg_j_1 - u0) / (Upeak_neg_j_1 - u0);
			}
			// Reloading is to the right of the break point
			if (u0 > Ubp_neg_j_1) {
				// Reloading stiffness KrelA after new excursion
				KrelA_neg_j_1 = (Fbp_neg_j_1) / (Ubp_neg_j_1 - u0);
				// Reloading stiffness KrelB after new excursion
				KrelB_neg_j_1 = (Fpeak_neg_j_1 - Fbp_neg_j_1) / (Upeak_neg_j_1 - Ubp_neg_j_1);
				df = ((ui - u0)*KrelA_neg_j_1) - fi_1;
				// Reloading is to the left of the break point
			}
			else {
				// Reloading stiffness after new excursion
				KrelB_neg_j_1 = (Fpeak_neg_j_1) / (Upeak_neg_j_1 - u0);
				df = ((ui - u0)*KrelB_neg_j_1) - fi_1;
			}
		}
		// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
			//// Update Unloading Stiffness
		if (Unloading_Flag == 1) {
			Kul_j_1 = (1 - betaK)*Kul_j_1;
			df = Kul_j_1*du;
		}
		// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
			//// Update Deformation at Residual Points
			// Deformation at residual onset for positive backbone
		Ures_pos_j_1 = (Fres_pos_j_1 - Fmax_pos_j_1 + Kpc_pos_j_1 * Umax_pos_j_1) / Kpc_pos_j_1;
		// Deformation at residual onset for negative backbone
		Ures_neg_j_1 = (Fres_neg_j_1 - Fmax_neg_j_1 + Kpc_neg_j_1 * Umax_neg_j_1) / Kpc_neg_j_1;
		// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
			//// Force
		fi = fi_1 + df;

		// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
		   //// Refine Peaks
		if (New_Peak_Pos_Flag == 1) {
			Fpeak_pos_j_1 = fi;
		}
		else if (New_Peak_Neg_Flag == 1) {
			Fpeak_neg_j_1 = fi;
		}
		// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
			//// Failure
			// Failure criteria (Tolerance = 1//)
		FailS = ((betaS < -0.01) || (betaS > 1.01));
		FailC = ((betaC < -0.01) || (betaC > 1.01));
		FailA = ((betaA < -0.01) || (betaA > 1.01));
		FailK = ((betaK < -0.01) || (betaK > 1.01));
		//cout << "  ENERGY: EtS=" << EtS << " EtC=" << EtC << " EtA=" << EtA << " EtK=" << EtK << endln;
		//cout << "  ENERGY: dEi=" << dEi << " Ei=" << Ei << " Epj=" << Epj << " Energy_Diss=" << Energy_Diss << " Energy_Acc=" << Energy_Acc << endln;
		//cout << "  ENERGY: betaS=" << betaS << " betaC=" << betaC << " betaA=" << betaA << " betaK=" << betaK << endln;
		//cout << "  FAIL:   FailS=" << FailS << " FailC=" << FailC << " FailA=" << FailA << " FailK=" << FailK << endln;

		if (FailS || FailC || FailA || FailK) {
			//cout << "  Energy Fail" << endln;
			Failure_Flag = 1;
		}
		if ((ui >= 0.0) && (ui >= Uu_pos)) {
			fi = 0;
			//cout << "  Rotation Fail" << endln;
			Failure_Flag = 1;
		}
		else if ((ui < 0.0) && (ui <= -Uu_neg)) {
			fi = 0;
			//cout << "  Rotation Fail" << endln;
			Failure_Flag = 1;
		}
		if ((Fpeak_pos_j_1 == 0) || (Fpeak_neg_j_1 == 0)) {
			fi = 0;
			//cout << "  Strength Fail" << endln;
			Failure_Flag = 1;
		}
	}
	else {
		fi = 0;
		dEi = 0;
		Epj = Energy_Acc + dEi;
	}


	//// Energy
	Energy_Acc = Energy_Acc + dEi; 	// Total internal energy accumulated until current increment
	Energy_Diss = Epj; 				// Total energy dissipated in all previous excursions

	//// Update Variables
	du_i_1 = du;
	Umaxp = 0;
	Umaxn = 0;

	// Tangent Stiffeness Calculation
	if (fi == fi_1) {
		TangentK = pow(10., -6);
		ki = pow(10., -6);
	}	
	
	if ((ui == ui_1)) {
		ki = Ke;
		fi = fi_1;
		TangentK = Ke;
	}
	else {
		ki = (fi - fi_1) / (du);
		TangentK = (fi - fi_1) / (du);
	}



	//cout << "  fi=" << fi << endln;
	//cout << "***********************" << endln;

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%% END OF MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	return 0;
}

double IMKPinching::getStress(void)
{
	//cout << " getStress" << endln;
	return (fi);
}

double IMKPinching::getTangent(void)
{
	//cout << " getTangent" << endln;
	return (TangentK);
}

double IMKPinching::getInitialTangent(void)
{
	//cout << " getInitialTangent" << endln;
	return (Ke);
}

double IMKPinching::getStrain(void)
{
	//cout << " getStrain" << endln;
	return (U);
}

int IMKPinching::commitState(void)
{
	//cout << " commitState" << endln;

	//commit trial  variables

	cU = U;

	cui = ui;
	cfi = fi;
	cui_1 = ui_1;
	cfi_1 = fi_1;

	cTangentK = TangentK;

	cdu_i_1 = du_i_1;

	cUy_pos_j_1 = Uy_pos_j_1;
	cUmax_pos_j_1 = Umax_pos_j_1;
	cUu_pos_j_1 = Uu_pos_j_1;
	cFy_pos_j_1 = Fy_pos_j_1;
	cFmax_pos_j_1 = Fmax_pos_j_1;
	cUpeak_pos_j_1 = Upeak_pos_j_1;
	cFpeak_pos_j_1 = Fpeak_pos_j_1;
	cUbp_pos_j_1 = Ubp_pos_j_1;
	cFbp_pos_j_1 = Fbp_pos_j_1;
	cUres_pos_j_1 = Ures_pos_j_1;
	cFres_pos_j_1 = Fres_pos_j_1;
	cKp_pos_j_1 = Kp_pos_j_1;
	cKpc_pos_j_1 = Kpc_pos_j_1;
	cKrelA_pos_j_1 = KrelA_pos_j_1;
	cKrelB_pos_j_1 = KrelB_pos_j_1;
	cPlastic_Offset_pos_j_1 = Plastic_Offset_pos_j_1;

	cUy_neg_j_1 = Uy_neg_j_1;
	cUmax_neg_j_1 = Umax_neg_j_1;
	cUu_neg_j_1 = Uu_neg_j_1;
	cFy_neg_j_1 = Fy_neg_j_1;
	cFmax_neg_j_1 = Fmax_neg_j_1;
	cUpeak_neg_j_1 = Upeak_neg_j_1;
	cFpeak_neg_j_1 = Fpeak_neg_j_1;
	cUbp_neg_j_1 = Ubp_neg_j_1;
	cFbp_neg_j_1 = Fbp_neg_j_1;
	cUres_neg_j_1 = Ures_neg_j_1;
	cFres_neg_j_1 = Fres_neg_j_1;
	cKp_neg_j_1 = Kp_neg_j_1;
	cKpc_neg_j_1 = Kpc_neg_j_1;
	cKrelA_neg_j_1 = KrelA_neg_j_1;
	cKrelB_neg_j_1 = KrelB_neg_j_1;
	cPlastic_Offset_neg_j_1 = Plastic_Offset_neg_j_1;

	cKul_j_1 = Kul_j_1;

	cFailure_Flag = Failure_Flag;
	cExcursion_Flag = Excursion_Flag;

	cEnergy_Acc = Energy_Acc;
	cEnergy_Diss = Energy_Diss;

	cUmaxp = Umaxp;
	cUmaxn = Umaxn;

	return 0;
}

int IMKPinching::revertToLastCommit(void)
{
	//cout << " revertToLastCommit" << endln;

	//the opposite of commit trial history variables
	U = cU;

	ui = cui;
	fi = cfi;
	ui_1 = cui_1;
	fi_1 = cfi_1;

	TangentK = cTangentK;

	du_i_1 = cdu_i_1;

	Uy_pos_j_1 = cUy_pos_j_1;
	Umax_pos_j_1 = cUmax_pos_j_1;
	Uu_pos_j_1 = cUu_pos_j_1;
	Fy_pos_j_1 = cFy_pos_j_1;
	Fmax_pos_j_1 = cFmax_pos_j_1;
	Upeak_pos_j_1 = cUpeak_pos_j_1;
	Fpeak_pos_j_1 = cFpeak_pos_j_1;
	Ubp_pos_j_1 = cUbp_pos_j_1;
	Fbp_pos_j_1 = cFbp_pos_j_1;
	Ures_pos_j_1 = cUres_pos_j_1;
	Fres_pos_j_1 = cFres_pos_j_1;
	Kp_pos_j_1 = cKp_pos_j_1;
	Kpc_pos_j_1 = cKpc_pos_j_1;
	KrelA_pos_j_1 = cKrelA_pos_j_1;
	KrelB_pos_j_1 = cKrelB_pos_j_1;
	Plastic_Offset_pos_j_1 = cPlastic_Offset_pos_j_1;

	Uy_neg_j_1 = cUy_neg_j_1;
	Umax_neg_j_1 = cUmax_neg_j_1;
	Uu_neg_j_1 = cUu_neg_j_1;
	Fy_neg_j_1 = cFy_neg_j_1;
	Fmax_neg_j_1 = cFmax_neg_j_1;
	Upeak_neg_j_1 = cUpeak_neg_j_1;
	Fpeak_neg_j_1 = cFpeak_neg_j_1;
	Ubp_neg_j_1 = cUbp_neg_j_1;
	Fbp_neg_j_1 = cFbp_neg_j_1;
	Ures_neg_j_1 = cUres_neg_j_1;
	Fres_neg_j_1 = cFres_neg_j_1;
	Kp_neg_j_1 = cKp_neg_j_1;
	Kpc_neg_j_1 = cKpc_neg_j_1;
	KrelA_neg_j_1 = cKrelA_neg_j_1;
	KrelB_neg_j_1 = cKrelB_neg_j_1;
	Plastic_Offset_neg_j_1 = cPlastic_Offset_neg_j_1;

	Kul_j_1 = cKul_j_1;

	Failure_Flag = cFailure_Flag;
	Excursion_Flag = cExcursion_Flag;

	Energy_Acc = cEnergy_Acc;
	Energy_Diss = cEnergy_Diss;

	Umaxp = cUmaxp;
	Umaxn = cUmaxn;

	return 0;
}

int IMKPinching::revertToStart(void)
{
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ONE TIME CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

	Failure_Flag = 0;
	Excursion_Flag = 0;
	Unloading_Flag = 0;
	New_Peak_Pos_Flag = 0;
	New_Peak_Neg_Flag = 0;

	betaS = 0;
	betaC = 0;
	betaK = 0;
	betaA = 0;

	Uy_pos = Fy_pos / Ke;
	Umax_pos = Uy_pos + Up_pos;
	Fmax_pos = FmaxFy_pos*Fy_pos;
	Kp_pos = (Fmax_pos - Fy_pos) / Up_pos;
	Kpc_pos = Fmax_pos / Upc_pos;

	Uy_neg = Fy_neg / Ke;
	Umax_neg = Uy_neg + Up_neg;
	Fmax_neg = FmaxFy_neg*Fy_neg;
	Kp_neg = (Fmax_neg - Fy_neg) / Up_neg;
	Kpc_neg = Fmax_neg / Upc_neg;

	Upeak_pos_j_1 = Uy_pos;
	Fpeak_pos_j_1 = Fy_pos;
	Upeak_neg_j_1 = -Uy_neg;
	Fpeak_neg_j_1 = -Fy_neg;

	Uy_pos_j_1 = Uy_pos;
	Fy_pos_j_1 = Fy_pos;
	Kp_pos_j_1 = Kp_pos;
	Kpc_pos_j_1 = -Kpc_pos;
	Uy_neg_j_1 = -Uy_neg;
	Fy_neg_j_1 = -Fy_neg;
	Kp_neg_j_1 = Kp_neg;
	Kpc_neg_j_1 = -Kpc_neg;

	Umax_pos_j_1 = Umax_pos;
	Fmax_pos_j_1 = Fmax_pos;
	Fres_pos_j_1 = Fy_pos*ResF_pos;
	Umax_neg_j_1 = -Umax_neg;
	Fmax_neg_j_1 = -Fmax_neg;
	Fres_neg_j_1 = -Fy_neg*ResF_neg;

	KrelA_pos_j_1 = Ke;
	KrelA_neg_j_1 = Ke;
	KrelB_pos_j_1 = Ke;
	KrelB_neg_j_1 = Ke;

	Kul_j_1 = Ke;

	Ures_pos_j_1 = (Fres_pos_j_1 - Fmax_pos_j_1) / Kpc_pos_j_1 + Umax_pos_j_1;
	Ures_neg_j_1 = (Fres_neg_j_1 - Fmax_neg_j_1) / Kpc_neg_j_1 + Umax_neg_j_1;

	Ubp_pos_j_1 = 0;
	Ubp_neg_j_1 = 0;

	Fbp_pos_j_1 = 0;
	Fbp_neg_j_1 = 0;

	Energy_Acc = 0.0;
	Energy_Diss = 0.0;

	Umaxp = 0.0;
	Umaxn = 0.0;

	EtS = LAMBDA_S *Fy_pos;
	EtC = LAMBDA_C *Fy_pos;
	EtA = LAMBDA_A *Fy_pos;
	EtK = LAMBDA_K *Fy_pos;

	cPlastic_Offset_pos_j_1 = 0;
	cPlastic_Offset_neg_j_1 = 0;

	cdu_i_1 = 0;

	cUpeak_pos_j_1 = Uy_pos;
	cFpeak_pos_j_1 = Fy_pos;
	cUpeak_neg_j_1 = -Uy_neg;
	cFpeak_neg_j_1 = -Fy_neg;

	cUy_pos_j_1 = Uy_pos;
	cFy_pos_j_1 = Fy_pos;
	cKp_pos_j_1 = Kp_pos;
	cKpc_pos_j_1 = -Kpc_pos;
	cUy_neg_j_1 = -Uy_neg;
	cFy_neg_j_1 = -Fy_neg;
	cKp_neg_j_1 = Kp_neg;
	cKpc_neg_j_1 = -Kpc_neg;

	cUmax_pos_j_1 = Umax_pos;
	cFmax_pos_j_1 = Fmax_pos;
	cFres_pos_j_1 = Fy_pos*ResF_pos;
	cUmax_neg_j_1 = -Umax_neg;
	cFmax_neg_j_1 = -Fmax_neg;
	cFres_neg_j_1 = -Fy_neg*ResF_neg;

	cKrelA_pos_j_1 = Ke;
	cKrelA_neg_j_1 = Ke;
	cKrelB_pos_j_1 = Ke;
	cKrelB_neg_j_1 = Ke;

	cKul_j_1 = Ke;

	cUres_pos_j_1 = (Fres_pos_j_1 - Fmax_pos_j_1) / Kpc_pos_j_1 + Umax_pos_j_1;
	cUres_neg_j_1 = (Fres_neg_j_1 - Fmax_neg_j_1) / Kpc_neg_j_1 + Umax_neg_j_1;

	cUbp_pos_j_1 = 0;
	cUbp_neg_j_1 = 0;

	cFbp_pos_j_1 = 0;
	cFbp_neg_j_1 = 0;

	//initially I zero everything   
	U = cU = 0;
	ui = 0;
	fi = 0;
	ui_1 = 0;
	fi_1 = 0;
	cui = 0;
	cfi = 0;
	cui_1 = 0;
	cfi_1 = 0;

	TangentK = Ke;
	cTangentK = Ke;
	//cout << " revertToStart:" << endln; //<< " U=" << U << " Ri=" << Ri << " TanK=" << TangentK << endln;

	return 0;
}

UniaxialMaterial *
IMKPinching::getCopy(void)
{
	IMKPinching *theCopy = new IMKPinching(this->getTag(), Ke,
		Uy_pos, Umax_pos, Uu_pos, Fy_pos, FmaxFy_pos, ResF_pos,
		Uy_neg, Umax_neg, Uu_neg, Fy_neg, FmaxFy_neg, ResF_neg,
		LAMBDA_S, LAMBDA_C, LAMBDA_A, LAMBDA_K, c_S, c_C, c_A, c_K, D_pos, D_neg, kappaF, kappaD);

	//cout << " getCopy" << endln;


	theCopy->U = U;
	theCopy->cU = cU;

	theCopy->TangentK = TangentK;

	theCopy->ui = ui;
	theCopy->fi = fi;
	theCopy->ui_1 = ui_1;
	theCopy->fi_1 = fi_1;
	theCopy->du_i_1 = du_i_1;

	theCopy->Uy_pos_j_1 = Uy_pos_j_1;
	theCopy->Umax_pos_j_1 = Umax_pos_j_1;
	theCopy->Uu_pos_j_1 = Uu_pos_j_1;
	theCopy->Fy_pos_j_1 = Fy_pos_j_1;
	theCopy->Fmax_pos_j_1 = Fmax_pos_j_1;
	theCopy->Upeak_pos_j_1 = Upeak_pos_j_1;
	theCopy->Fpeak_pos_j_1 = Fpeak_pos_j_1;
	theCopy->Ubp_pos_j_1 = Ubp_pos_j_1;
	theCopy->Fbp_pos_j_1 = Fbp_pos_j_1;
	theCopy->Ures_pos_j_1 = Ures_pos_j_1;
	theCopy->Fres_pos_j_1 = Fres_pos_j_1;
	theCopy->Kp_pos_j_1 = Kp_pos_j_1;
	theCopy->Kpc_pos_j_1 = Kpc_pos_j_1;
	theCopy->KrelA_pos_j_1 = KrelA_pos_j_1;
	theCopy->KrelB_pos_j_1 = KrelB_pos_j_1;
	theCopy->Plastic_Offset_pos_j_1 = Plastic_Offset_pos_j_1;

	theCopy->Uy_neg_j_1 = Uy_neg_j_1;
	theCopy->Umax_neg_j_1 = Umax_neg_j_1;
	theCopy->Uu_neg_j_1 = Uu_neg_j_1;
	theCopy->Fy_neg_j_1 = Fy_neg_j_1;
	theCopy->Fmax_neg_j_1 = Fmax_neg_j_1;
	theCopy->Upeak_neg_j_1 = Upeak_neg_j_1;
	theCopy->Fpeak_neg_j_1 = Fpeak_neg_j_1;
	theCopy->Ubp_neg_j_1 = Ubp_neg_j_1;
	theCopy->Fbp_neg_j_1 = Fbp_neg_j_1;
	theCopy->Ures_neg_j_1 = Ures_neg_j_1;
	theCopy->Fres_neg_j_1 = Fres_neg_j_1;
	theCopy->Kp_neg_j_1 = Kp_neg_j_1;
	theCopy->Kpc_neg_j_1 = Kpc_neg_j_1;
	theCopy->KrelA_neg_j_1 = KrelA_neg_j_1;
	theCopy->KrelB_neg_j_1 = KrelB_neg_j_1;
	theCopy->Plastic_Offset_neg_j_1 = Plastic_Offset_neg_j_1;

	theCopy->Kul_j_1 = Kul_j_1;

	theCopy->Failure_Flag = Failure_Flag;
	theCopy->Excursion_Flag = Excursion_Flag;

	theCopy->Energy_Acc = Energy_Acc;
	theCopy->Energy_Diss = Energy_Diss;

	theCopy->Umaxp = Umaxp;
	theCopy->Umaxn = Umaxn;



	theCopy->cTangentK = cTangentK;

	theCopy->cui = cui;
	theCopy->cfi = cfi;
	theCopy->cui_1 = cui_1;
	theCopy->cfi_1 = cfi_1;
	theCopy->cdu_i_1 = cdu_i_1;

	theCopy->cUy_pos_j_1 = cUy_pos_j_1;
	theCopy->cUmax_pos_j_1 = cUmax_pos_j_1;
	theCopy->cUu_pos_j_1 = cUu_pos_j_1;
	theCopy->cFy_pos_j_1 = cFy_pos_j_1;
	theCopy->cFmax_pos_j_1 = cFmax_pos_j_1;
	theCopy->cUpeak_pos_j_1 = cUpeak_pos_j_1;
	theCopy->cFpeak_pos_j_1 = cFpeak_pos_j_1;
	theCopy->cUbp_pos_j_1 = cUbp_pos_j_1;
	theCopy->cFbp_pos_j_1 = cFbp_pos_j_1;
	theCopy->cUres_pos_j_1 = cUres_pos_j_1;
	theCopy->cFres_pos_j_1 = cFres_pos_j_1;
	theCopy->cKp_pos_j_1 = cKp_pos_j_1;
	theCopy->cKpc_pos_j_1 = cKpc_pos_j_1;
	theCopy->cKrelA_pos_j_1 = cKrelA_pos_j_1;
	theCopy->cKrelB_pos_j_1 = cKrelB_pos_j_1;
	theCopy->cPlastic_Offset_pos_j_1 = cPlastic_Offset_pos_j_1;

	theCopy->cUy_neg_j_1 = cUy_neg_j_1;
	theCopy->cUmax_neg_j_1 = cUmax_neg_j_1;
	theCopy->cUu_neg_j_1 = cUu_neg_j_1;
	theCopy->cFy_neg_j_1 = cFy_neg_j_1;
	theCopy->cFmax_neg_j_1 = cFmax_neg_j_1;
	theCopy->cUpeak_neg_j_1 = cUpeak_neg_j_1;
	theCopy->cFpeak_neg_j_1 = cFpeak_neg_j_1;
	theCopy->cUbp_neg_j_1 = cUbp_neg_j_1;
	theCopy->cFbp_neg_j_1 = cFbp_neg_j_1;
	theCopy->cUres_neg_j_1 = cUres_neg_j_1;
	theCopy->cFres_neg_j_1 = cFres_neg_j_1;
	theCopy->cKp_neg_j_1 = cKp_neg_j_1;
	theCopy->cKpc_neg_j_1 = cKpc_neg_j_1;
	theCopy->cKrelA_neg_j_1 = cKrelA_neg_j_1;
	theCopy->cKrelB_neg_j_1 = cKrelB_neg_j_1;
	theCopy->cPlastic_Offset_neg_j_1 = cPlastic_Offset_neg_j_1;

	theCopy->cKul_j_1 = cKul_j_1;

	theCopy->cFailure_Flag = cFailure_Flag;
	theCopy->cExcursion_Flag = cExcursion_Flag;

	theCopy->cEnergy_Acc = cEnergy_Acc;
	theCopy->cEnergy_Diss = cEnergy_Diss;

	theCopy->cUmaxp = cUmaxp;
	theCopy->cUmaxn = cUmaxn;

	return theCopy;
}

int IMKPinching::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	cout << " sendSelf" << endln;

	static Vector data(149);
	data(0) = this->getTag();
	data(1) = Ke;
	data(2) = Uy_pos;
	data(3) = Umax_pos;
	data(4) = Uu_pos;
	data(5) = Fy_pos;
	data(6) = FmaxFy_pos;
	data(7) = ResF_pos;
	data(8) = Uy_neg;
	data(9) = Umax_neg;
	data(10) = Uu_neg;
	data(11) = Fy_neg;
	data(12) = FmaxFy_neg;
	data(13) = ResF_neg;
	data(14) = LAMBDA_S;
	data(15) = LAMBDA_C;
	data(16) = LAMBDA_A;
	data(17) = LAMBDA_K;
	data(18) = c_S;
	data(19) = c_C;
	data(20) = c_A;
	data(21) = c_K;
	data(22) = D_pos;
	data(23) = D_neg;
	data(24) = kappaF;
	data(25) = kappaD;

	data(26) = ui;
	data(27) = fi;
	data(28) = ui_1;
	data(29) = fi_1;
	data(30) = du_i_1;

	data(31) = Uy_pos_j_1;
	data(32) = Umax_pos_j_1;
	data(33) = Uu_pos_j_1;
	data(34) = Fy_pos_j_1;
	data(35) = Fmax_pos_j_1;
	data(36) = Upeak_pos_j_1;
	data(37) = Fpeak_pos_j_1;
	data(38) = Ubp_pos_j_1;
	data(39) = Fbp_pos_j_1;
	data(40) = Ures_pos_j_1;
	data(41) = Fres_pos_j_1;
	data(42) = Kp_pos_j_1;
	data(43) = Kpc_pos_j_1;
	data(44) = KrelA_pos_j_1;
	data(45) = KrelB_pos_j_1;
	data(46) = Plastic_Offset_pos_j_1;

	data(47) = Uy_neg_j_1;
	data(48) = Umax_neg_j_1;
	data(49) = Uu_neg_j_1;
	data(50) = Fy_neg_j_1;
	data(51) = Fmax_neg_j_1;
	data(52) = Upeak_neg_j_1;
	data(53) = Fpeak_neg_j_1;
	data(54) = Ubp_neg_j_1;
	data(55) = Fbp_neg_j_1;
	data(56) = Ures_neg_j_1;
	data(57) = Fres_neg_j_1;
	data(58) = Kp_neg_j_1;
	data(59) = Kpc_neg_j_1;
	data(60) = KrelA_neg_j_1;
	data(61) = KrelB_neg_j_1;
	data(62) = Plastic_Offset_neg_j_1;

	data(63) = Kul_j_1;
	data(64) = Failure_Flag;
	data(65) = Excursion_Flag;

	data(66) = Energy_Acc;
	data(67) = Energy_Diss;
	data(68) = Umaxp;
	data(69) = Umaxn;

	data(70) = u0;
	data(71) = du;
	data(72) = df;

	data(73) = Unloading_Flag;
	data(74) = New_Peak_Pos_Flag;
	data(75) = New_Peak_Neg_Flag;

	data(76) = Fp;

	data(77) = FailS;
	data(78) = FailC;
	data(79) = FailA;
	data(80) = FailK;

	data(81) = Ei;
	data(82) = dEi;
	data(83) = Epj;
	data(84) = EpjK;
	data(85) = EiK;

	data(86) = c_S;
	data(87) = c_C;
	data(88) = c_A;
	data(89) = c_K;
	data(90) = EtS;
	data(91) = EtC;
	data(92) = EtA;
	data(93) = EtK;
	data(94) = betaS;
	data(95) = betaC;
	data(96) = betaA;
	data(97) = betaK;
	data(98) = sPCsp;
	data(99) = sPCpcp;

	data(100) = TangentK;

	data(101) = Uy_pos;
	data(102) = Umax_pos;
	data(103) = Fmax_pos;
	data(104) = Kp_pos;
	data(105) = Kpc_pos;

	data(106) = Uy_neg;
	data(107) = Umax_neg;
	data(108) = Fmax_neg;
	data(109) = Kp_neg;
	data(110) = Kpc_neg;

	data(111) = cui;
	data(112) = cfi;
	data(113) = cui_1;
	data(114) = cfi_1;
	data(115) = cdu_i_1;

	data(116) = cUy_pos_j_1;
	data(117) = cUmax_pos_j_1;
	data(118) = cUu_pos_j_1;
	data(119) = cFy_pos_j_1;
	data(120) = cFmax_pos_j_1;
	data(121) = cUpeak_pos_j_1;
	data(122) = cFpeak_pos_j_1;
	data(123) = cUbp_pos_j_1;
	data(124) = cFbp_pos_j_1;
	data(125) = cUres_pos_j_1;
	data(126) = cFres_pos_j_1;
	data(127) = cKp_pos_j_1;
	data(128) = cKpc_pos_j_1;
	data(129) = cKrelA_pos_j_1;
	data(130) = cKrelB_pos_j_1;
	data(131) = cPlastic_Offset_pos_j_1;

	data(132) = cUy_neg_j_1;
	data(133) = cUmax_neg_j_1;
	data(134) = cUu_neg_j_1;
	data(135) = cFy_neg_j_1;
	data(136) = cFmax_neg_j_1;
	data(137) = cUpeak_neg_j_1;
	data(138) = cFpeak_neg_j_1;
	data(139) = cUbp_neg_j_1;
	data(140) = cFbp_neg_j_1;
	data(141) = cUres_neg_j_1;
	data(142) = cFres_neg_j_1;
	data(143) = cKp_neg_j_1;
	data(144) = cKpc_neg_j_1;
	data(145) = cKrelA_neg_j_1;
	data(146) = cKrelB_neg_j_1;
	data(147) = cPlastic_Offset_neg_j_1;

	data(148) = cKul_j_1;


	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "IMKPinching::sendSelf() - failed to send data\n";

	return res;
}

int IMKPinching::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector data(149);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);

	if (res < 0) {
		opserr << "IMKPinching::recvSelf() - failed to receive data\n";
		this->setTag(0);
	}
	else {
		cout << " recvSelf" << endln;
		this->setTag((int)data(0));
		Ke = data(1);
		Up_pos = data(2);
		Upc_pos = data(3);
		Uu_pos = data(4);
		Fy_pos = data(5);
		FmaxFy_pos = data(6);
		ResF_pos = data(7);
		Up_neg = data(8);
		Upc_neg = data(9);
		Uu_neg = data(10);
		Fy_neg = data(11);
		FmaxFy_neg = data(12);
		ResF_neg = data(13);
		LAMBDA_S = data(14);
		LAMBDA_C = data(15);
		LAMBDA_A = data(16);
		LAMBDA_K = data(17);
		c_S = data(18);
		c_C = data(19);
		c_A = data(20);
		c_K = data(21);
		D_pos = data(22);
		D_neg = data(23);
		kappaF = data(24);
		kappaD = data(25);

		ui = data(26);
		fi = data(27);
		ui_1 = data(28);
		fi_1 = data(29);
		du_i_1 = data(30);

		Uy_pos_j_1 = data(31);
		Umax_pos_j_1 = data(32);
		Uu_pos_j_1 = data(33);
		Fy_pos_j_1 = data(34);
		Fmax_pos_j_1 = data(35);
		Upeak_pos_j_1 = data(36);
		Fpeak_pos_j_1 = data(37);
		Ubp_pos_j_1 = data(38);
		Fbp_pos_j_1 = data(39);
		Ures_pos_j_1 = data(40);
		Fres_pos_j_1 = data(41);
		Kp_pos_j_1 = data(42);
		Kpc_pos_j_1 = data(43);
		KrelA_pos_j_1 = data(44);
		KrelB_pos_j_1 = data(45);
		Plastic_Offset_pos_j_1 = data(46);

		Uy_neg_j_1 = data(47);
		Umax_neg_j_1 = data(48);
		Uu_neg_j_1 = data(49);
		Fy_neg_j_1 = data(50);
		Fmax_neg_j_1 = data(51);
		Upeak_neg_j_1 = data(52);
		Fpeak_neg_j_1 = data(53);
		Ubp_neg_j_1 = data(54);
		Fbp_neg_j_1 = data(55);
		Ures_neg_j_1 = data(56);
		Fres_neg_j_1 = data(57);
		Kp_neg_j_1 = data(58);
		Kpc_neg_j_1 = data(59);
		KrelA_neg_j_1 = data(60);
		KrelB_neg_j_1 = data(61);
		Plastic_Offset_neg_j_1 = data(62);

		Kul_j_1 = data(63);

		Failure_Flag = data(64);
		Excursion_Flag = data(65);

		Energy_Acc = data(66);
		Energy_Diss = data(67);
		Umaxp = data(68);
		Umaxn = data(69);

		u0 = data(70);
		du = data(71);
		df = data(72);

		Unloading_Flag = data(73);
		New_Peak_Pos_Flag = data(74);
		New_Peak_Neg_Flag = data(75);

		Fp = data(76);

		FailS = data(77);
		FailC = data(78);
		FailA = data(79);
		FailK = data(80);

		Ei = data(81);
		dEi = data(82);
		Epj = data(83);
		EpjK = data(84);
		EiK = data(85);

		c_S = data(86);
		c_C = data(87);
		c_A = data(88);
		c_K = data(89);
		EtS = data(90);
		EtC = data(91);
		EtA = data(92);
		EtK = data(93);
		betaS = data(94);
		betaC = data(95);
		betaA = data(96);
		betaK = data(97);
		sPCsp = data(98);
		sPCpcp = data(99);

		TangentK = data(100);

		Uy_pos = data(101);
		Umax_pos = data(102);
		Fmax_pos = data(103);
		Kp_pos = data(104);
		Kpc_pos = data(105);

		Uy_neg = data(106);
		Umax_neg = data(107);
		Fmax_neg = data(108);
		Kp_neg = data(109);
		Kpc_neg = data(110);

		cui			= data(111);
		cfi			= data(112);
		cui_1		= data(113);
		cfi_1		= data(114);
		cdu_i_1		= data(115);

		cUy_pos_j_1		= data(116);
		cUmax_pos_j_1	= data(117);
		cUu_pos_j_1		= data(118);
		cFy_pos_j_1		= data(119);
		cFmax_pos_j_1	= data(120);
		cUpeak_pos_j_1	= data(121);
		cFpeak_pos_j_1	= data(122);
		cUbp_pos_j_1	= data(123);
		cFbp_pos_j_1	= data(124);
		cUres_pos_j_1	= data(125);
		cFres_pos_j_1	= data(126);
		cKp_pos_j_1		= data(127);
		cKpc_pos_j_1	= data(128);
		cKrelA_pos_j_1	= data(129);
		cKrelB_pos_j_1	= data(130);
		cPlastic_Offset_pos_j_1 = data(131);

		cUy_neg_j_1		= data(132);
		cUmax_neg_j_1	= data(133);
		cUu_neg_j_1		= data(134);
		cFy_neg_j_1		= data(135);
		cFmax_neg_j_1	= data(136);
		cUpeak_neg_j_1	= data(137);
		cFpeak_neg_j_1	= data(138);
		cUbp_neg_j_1	= data(139);
		cFbp_neg_j_1	= data(140);
		cUres_neg_j_1	= data(141);
		cFres_neg_j_1	= data(142);
		cKp_neg_j_1		= data(143);
		cKpc_neg_j_1	= data(144);
		cKrelA_neg_j_1	= data(145);
		cKrelB_neg_j_1	= data(146);
		cPlastic_Offset_neg_j_1 = data(147);

		cKul_j_1 = data(148);
	}

	return res;
}

void IMKPinching::Print(OPS_Stream &s, int flag)
{
	cout << "IMKPinching tag: " << this->getTag() << endln;
}

