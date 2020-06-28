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
#include <IMKPeakOriented.h>
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

static int numIMKPeakOrientedMaterials = 0;

void *
OPS_IMKPeakOriented()
{
	if (numIMKPeakOrientedMaterials == 0) {
		numIMKPeakOrientedMaterials++;
		OPS_Error("IMK Model with Peak-Oriented Response - Code by A. ELKADY & H. ELJISR (July 2020)\n", 1);
	}

	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial = 0;

	int    iData[1];
	double dData[23];
	int numData = 1;

	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial IMKPeakOriented tag" << endln;
		return 0;
	}

	numData = 23;


	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid Args want: uniaxialMaterial IMKPeakOriented tag? Ke? ";
		opserr << "Up_pos? Upc_pos? Uu_pos? Fy_pos? FmaxFy_pos? ResF_pos? ";
		opserr << "Up_neg? Upc_neg? Uu_neg? Fy_neg? FmaxFy_neg? ResF_neg? ";
		opserr << "LamdaS? LamdaC? LamdaA? LamdaK? Cs? Cc? Ca? Ck? D_pos? D_neg? ";
		return 0;
	}



	// Parsing was successful, allocate the material
	theMaterial = new IMKPeakOriented(iData[0],
		dData[0],
		dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
		dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
		dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20],
		dData[21], dData[22]);

	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type IMKPeakOriented Material\n";
		return 0;
	}

	return theMaterial;
}

IMKPeakOriented::IMKPeakOriented(int tag, double p_Ke,
	double p_Up_pos, double p_Upc_pos, double p_Uu_pos, double p_Fy_pos, double p_FmaxFy_pos, double p_ResF_pos,
	double p_Up_neg, double p_Upc_neg, double p_Uu_neg, double p_Fy_neg, double p_FmaxFy_neg, double p_ResF_neg,
	double p_LAMBDA_S, double p_LAMBDA_C, double p_LAMBDA_A, double p_LAMBDA_K, double p_c_S, double p_c_C, double p_c_A, double p_c_K, double p_D_pos, double p_D_neg)
	: UniaxialMaterial(tag, 0), Ke(p_Ke),
	Up_pos(p_Up_pos), Upc_pos(p_Upc_pos), Uu_pos(p_Uu_pos), Fy_pos(p_Fy_pos), FmaxFy_pos(p_FmaxFy_pos), ResF_pos(p_ResF_pos),
	Up_neg(p_Up_neg), Upc_neg(p_Upc_neg), Uu_neg(p_Uu_neg), Fy_neg(p_Fy_neg), FmaxFy_neg(p_FmaxFy_neg), ResF_neg(p_ResF_neg),
	LAMBDA_S(p_LAMBDA_S), LAMBDA_C(p_LAMBDA_C), LAMBDA_A(p_LAMBDA_A), LAMBDA_K(p_LAMBDA_K), c_S(p_c_S), c_C(p_c_C), c_A(p_c_A), c_K(p_c_K), D_pos(p_D_pos), D_neg(p_D_neg)
{
	this->revertToStart();
}

IMKPeakOriented::IMKPeakOriented()
	:UniaxialMaterial(0, 0), Ke(0),
	Up_pos(0), Upc_pos(0), Uu_pos(0), Fy_pos(0), FmaxFy_pos(0), ResF_pos(0),
	Up_neg(0), Upc_neg(0), Uu_neg(0), Fy_neg(0), FmaxFy_neg(0), ResF_neg(0),
	LAMBDA_S(0), LAMBDA_C(0), LAMBDA_A(0), LAMBDA_K(0), c_S(0), c_C(0), c_A(0), c_K(0), D_pos(0), D_neg(0)
{
	this->revertToStart();
}

IMKPeakOriented::~IMKPeakOriented()
{
	// does nothing
}

int IMKPeakOriented::setTrialStrain(double strain, double strainRate)
{
	//all variables to the last commit
	this->revertToLastCommit();

	//state determination algorithm: defines the current force and tangent stiffness
	U = strain; //set trial displacement
	ui_1 = ui;
	fi_1 = fi;
	ui = U;

	//cout << "***********************" << endln;
	//cout << "  +VE: Uy= " << Uy_pos_j_1 << " Umax= " << Umax_pos_j_1 << " Upeak= " << Upeak_pos_j_1 << " Fpeak= " << Fpeak_pos_j_1 << " Krel=" << Krel_j_1 << endln;
	//cout << "  -VE: Uy=" << Uy_neg_j_1 << " Umax=" << Umax_neg_j_1 << " Upeak=" << Upeak_neg_j_1  << " Fpeak=" << Fpeak_neg_j_1 << " Krel=" << Krel_j_1  << endln;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////  MAIN CODE //////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// Incremental deformation at current step
	du = ui - ui_1;


	if (Failure_Flag != 1) {
		
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		////////////////// INITIAL FLAGS CHECKS AND MAIN POINTS COORDINATES ///////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////


		// CHECK FOR UNLOADING
		if ((fi_1 > 0) && (du <= 0) && (du*du_i_1 <= 0)) {
			Unloading_Flag = 1;
			Reversal_Flag = 1;
			Reloading_Flag = 0;
			K_check=(FLastPeak_pos_j_1-fi_1)/(ULastPeak_pos_j_1-ui_1);
			if ((K_check >=1.05*Kul_j_1) || (K_check <=0.95*Kul_j_1)) { // a tailored criteria to avoid registering last peak points during small unload/reload excursions on the unloading branch 
				FLastPeak_pos_j_1 = fi_1;
				ULastPeak_pos_j_1 = ui_1;
			}
		}
		else if ((fi_1 < 0) && (du > 0) && (du*du_i_1 <= 0)) {
			Unloading_Flag = 1;
			Reversal_Flag = 1;
			Reloading_Flag = 0;
			K_check=(FLastPeak_neg_j_1-fi_1)/(ULastPeak_neg_j_1-ui_1);
			if ((K_check >=1.01*Kul_j_1) || (K_check <=0.99*Kul_j_1)) {
				FLastPeak_neg_j_1 = fi_1;
				ULastPeak_neg_j_1 = ui_1;
			}
		}
		else {
			Reversal_Flag = 0;
		}

		// CHECK FOR RELOADING
		if      ((fi_1 > 0) && (du > 0) && (du_i_1 < 0)) {
			Reloading_Flag = 1;
			Unloading_Flag = 0;
		}
		else if ((fi_1 < 0) && (du < 0) && (du_i_1 > 0)) {
			Reloading_Flag = 1;
			Unloading_Flag = 0;
		}


		// CHECK FOR NEW EXCURSION
		if		((fi_1 < 0) && (fi_1 + du * Kul_j_1 >= 0)) {
			Excursion_Flag = 1;
			Reloading_Flag = 0;
			Unloading_Flag = 0;
			u0 = ui_1 - (fi_1 / Kul_j_1);
		}
		else if ((fi_1 > 0) && (fi_1 + du * Kul_j_1 <= 0)) {
			Excursion_Flag = 1;
			Reloading_Flag = 0;
			Unloading_Flag = 0;
			u0 = ui_1 - (fi_1 / Kul_j_1);
		}
		else {
			Excursion_Flag = 0;
		}

		// UPDATE GLOBAL PEAK POINTS
		if ((fi_1 >= 0) && (ui_1 >= Upeak_pos_j_1)) {
			Upeak_pos_j_1 = ui_1;
			Fpeak_pos_j_1 = fi_1;
		} else if ((fi_1 < 0) && (ui_1 <= Upeak_neg_j_1)) {
			Upeak_neg_j_1 = ui_1;
			Fpeak_neg_j_1 = fi_1;
		}
		
		// CHECK FOR YIELDING
		if     ((Upeak_pos_j_1 > Uy_pos_j_1) || (Upeak_neg_j_1 < Uy_neg_j_1)) {
			Yield_Flag = 1;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		/////////////////// UPDATE DETERIORATION PARAMETERS AND BACKBONE CURVE ////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////

		// UPDATE DETERIORATION PARAMETERS AT EACH NEW EXCURSION	
		//cout << "  ENERGY: dEi=" << dEi << " Kul=" << Kul_j_1 << " du=" << du << " df=" << df << endln;

		if (Excursion_Flag == 1) {
			//Epj = Energy_Acc + dEi;
			Ei = fmax(0, Energy_Acc - Energy_Diss);
			betaS = pow((Ei / (EtS - Energy_Acc)), c_S);
			betaC = pow((Ei / (EtC - Energy_Acc)), c_C);
			betaA = pow((Ei / (EtA - Energy_Acc)), c_A);
			Energy_Diss = Energy_Acc;
		}
		else {
			//Epj = Energy_Diss;
			betaS = 0;
			betaC = 0;
			betaA = 0;
		}

		if (Reversal_Flag == 1) {
			EpjK = Energy_Acc - 0.5*(fi_1 / Kul_j_1)*fi_1;
			EiK = Energy_Acc - Energy_Diss + 0.5*(fi_1 / Kul_j_1)*fi_1;
			betaK = pow((EiK / (EtK - EpjK)), c_K);
			Kul_j_1 = Kul_j_1 * (1 - betaK);
		}
		else {
			betaK = 0;
		}

		// Update Positive Backbone and Target Peak Point
		if (Excursion_Flag == 1) {
			// Positive loading backbone
			if ((fi_1 < 0) && (Yield_Flag==1)) {
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
				if (Upeak_pos_j_1 <= Uy_pos_j_1) {	
					Fpeak_pos_j_1 = Ke*Upeak_pos_j_1;
					// Target peak deformation in post-yield branch of the updated backbone
				} else if (Upeak_pos_j_1 <= Umax_pos_j_1) {
					Fpeak_pos_j_1 = Kp_pos_j_1 *(Upeak_pos_j_1 - Uy_pos_j_1) + Fy_pos_j_1;
					// Target peak deformation in post-capping branch of the updated backbone
				} else {
					Fpeak_pos_j_1 = max(Kpc_pos_j_1*(Upeak_pos_j_1 - Umax_pos_j_1) + Fmax_pos_j_1, Fres_pos_j_1);
				}
			} 
			else if ((fi_1 >= 0) && (Yield_Flag==1)) {
				// Update Negative Backbone and Target Peak Point
				// Basic strength deterioration: Yield point
				Uy_neg_j_1 = min(Uy_neg_j_1 - Fy_neg_j_1 *betaS* D_neg / Ke, Fres_neg_j_1 / Ke);
				Fy_neg_j_1 = min(Fy_neg_j_1 *(1 - betaS* D_neg), Fres_neg_j_1);
				// Basic strength deterioration: Post-yield stiffness
				if (Fy_neg_j_1 != Fres_neg_j_1) {
					Kp_neg_j_1 = Kp_neg_j_1 *(1 - betaS* D_neg);
				} else {
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

		// Update Deformation at Residual Points
		Ures_pos_j_1 = (Fres_pos_j_1 - Fmax_pos_j_1 + Kpc_pos_j_1 * Umax_pos_j_1) / Kpc_pos_j_1;
		Ures_neg_j_1 = (Fres_neg_j_1 - Fmax_neg_j_1 + Kpc_neg_j_1 * Umax_neg_j_1) / Kpc_neg_j_1;
		
		// CHECK TARGET POINT: LAST CYCLE PEAK or GLOBAL PEAK
		if (Excursion_Flag == 1) {
			if (du >= 0) {
				Krel_LastPeak   = FLastPeak_pos_j_1 / (ULastPeak_pos_j_1 - u0);
				Krel_GlobalPeak = Fpeak_pos_j_1 	/ (Upeak_pos_j_1 	 - u0);
			}
			else {
				Krel_LastPeak   = FLastPeak_neg_j_1 / (ULastPeak_neg_j_1 - u0);
				Krel_GlobalPeak = Fpeak_neg_j_1 	/ (Upeak_neg_j_1 	 - u0);
			}
			
			if ((du>=0) && (FLastPeak_pos_j_1 >= Fpeak_pos_j_1)) {
				TargetPeak_Flag=0;
			} else if ((du<=0) && (FLastPeak_neg_j_1 <= Fpeak_neg_j_1)) {
				TargetPeak_Flag=0;            
			} else if (abs(Krel_LastPeak) <= abs(Krel_GlobalPeak)) {
				TargetPeak_Flag = 0;
			}
			else {
				TargetPeak_Flag = 1;
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////// COMPUTE FORCE INCREMENT /////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////

		// Positive Force
		if (fi_1 + du * Kul_j_1 >= 0) {
		   
			// CASE 0: At THE ELASTIC SLOPE
			if ((ui>=0) && (Upeak_pos_j_1 <= Uy_pos_j_1) && (Yield_Flag==0)) {
				if (ui >= Uy_pos_j_1) {
					df = Ke*(Uy_pos_j_1 - ui_1) + Kp_pos_j_1*(ui - Uy_pos_j_1);
				} else {
					df = du * Ke;
				}
				//cout << "  Case = 0+" << endln;

			// CASE 1: EACH NEW EXCURSION
			} else if (Excursion_Flag==1) {
				if (TargetPeak_Flag==0) {
					Krel_j_1  = Fpeak_pos_j_1 	  / (Upeak_pos_j_1 	   - u0);
				} else {
					Krel_j_1  = FLastPeak_pos_j_1 / (ULastPeak_pos_j_1 - u0);
				}
				df = Kul_j_1*(u0 - ui_1) + Krel_j_1*(ui - u0);
				//cout << "  Case = 1+" << endln;

			// CASE 2: WHEN RELOADING
			} else if ((Reloading_Flag==1) && (ui <= ULastPeak_pos_j_1)) {
				df = du * Kul_j_1;
				//cout << "  Case = 2+" << endln;

			// CASE 2: WHEN UNLOADING
			} else if (Unloading_Flag==1) {
				df = du * Kul_j_1;
				//cout << "  Case = 2+" << endln;

			// CASE 3: WHEN RELOADING BUT BETWEEN LAST CYCLE PEAK POINT AND GLOBAL PEAK POINT
			} else if ((Reloading_Flag==1) && (ui >= ULastPeak_pos_j_1) && (ui <= Upeak_pos_j_1)) {
				Krel_j_1  = (Fpeak_pos_j_1 - FLastPeak_pos_j_1) / (Upeak_pos_j_1 - ULastPeak_pos_j_1);
				if (ui_1 <= ULastPeak_pos_j_1) {
					df = Kul_j_1*(ULastPeak_pos_j_1 - ui_1) + Krel_j_1*(ui - ULastPeak_pos_j_1);
				} else {
					df = du * Krel_j_1;
				}
				//cout << "  Case = 3+" << endln;

			// CASE 4: WHEN LOADING IN GENERAL TOWARDS THE TARGET PEAK
			} else if ((du >= 0) && (((TargetPeak_Flag==0) && (ui <= Upeak_pos_j_1)) || ((TargetPeak_Flag==1) && (ui <= ULastPeak_pos_j_1)))) {
				if (TargetPeak_Flag==0) {
					Krel_j_1 = (Fpeak_pos_j_1 	  - fi_1) / (Upeak_pos_j_1 	   - ui_1);
				} else {
					Krel_j_1 = (FLastPeak_pos_j_1 - fi_1) / (ULastPeak_pos_j_1 - ui_1);                
				}
				df = du * Krel_j_1;
				//cout << "  Case = 4+" << endln;

			// CASE 5: WHEN LOADING IN GENERAL TOWARDS THE LAST CYCLE PEAK POINT BUT BEYOND IT
			} else if ((du >= 0) && (TargetPeak_Flag==1) && (ui >= ULastPeak_pos_j_1) && (ui <= Upeak_pos_j_1)) {
				Krel_j_1 = (Fpeak_pos_j_1 - FLastPeak_pos_j_1) / (Upeak_pos_j_1 - ULastPeak_pos_j_1);
				if (ui_1 <= ULastPeak_pos_j_1) {
					df = (FLastPeak_pos_j_1 - fi_1) + Krel_j_1 * (ui - ULastPeak_pos_j_1);
				} else {
					df = du * Krel_j_1;
				}
				//cout << "  Case = 5+" << endln;

			// CASE 6: WHEN LOADING BEYOND THE TARGET PEAK BUT BEFORE THE CAPPING POINT
			} else if ((du >= 0) && (ui <= Umax_pos_j_1))  {
				df = du * Kp_pos_j_1;
				//cout << "  Case = 6+" << endln;

			// CASE 7: WHEN LOADING AND BETWEEN THE CAPPING POINT AND THE RESIDUAL POINT
			} else if ((du > 0) && (ui >= Umax_pos_j_1) && (ui <= Ures_pos_j_1)) {
				if ((ui_1<= Umax_pos_j_1) && (ui >= Umax_pos_j_1)) {
					df = Kp_pos_j_1 * (Umax_pos_j_1 - ui_1) + Kpc_pos_j_1 * (ui - Umax_pos_j_1);
				} else {
					df = du * Kpc_pos_j_1;
				}
				//cout << "  Case = 7+" << endln;

			// CASE 8: WHEN LOADING AND BEYOND THE RESIDUAL POINT
			} else if ((du > 0) && (ui >= Ures_pos_j_1)) {
				df = 0.0;
				if (Fres_pos_j_1 == 0) {
					Failure_Flag = 1;
				}
				//cout << "  Case = 8+" << endln;
			}
		}
		
		// Negative Force
		if (fi_1 + du * Kul_j_1 <= 0) {
			
			// CASE 0: At THE ELASTIC SLOPE
			if ((ui<=0) && (Upeak_neg_j_1 >= Uy_neg_j_1) && (Yield_Flag==0)) {
				if (ui <= Uy_neg_j_1) {
					df = Ke*(Uy_neg_j_1 - ui_1) + Kp_neg_j_1 * (ui - Uy_neg_j_1);
			} else {
					df = du * Ke;
			}
				//cout << "  Case = 0-" << endln;

			// CASE 1: EACH NEW EXCURSION
			} else if (Excursion_Flag==1) {
				if (TargetPeak_Flag==0) {
					Krel_j_1 = Fpeak_neg_j_1 	 / (Upeak_neg_j_1 	  - u0);
				} else {
					Krel_j_1 = FLastPeak_neg_j_1 / (ULastPeak_neg_j_1 - u0);
				}
				df = Kul_j_1 * (u0 - ui_1) + Krel_j_1 * (ui - u0);
				//cout << "  Case = 1-" << endln;

			// CASE 2: WHEN RELOADING
			} else if ((Reloading_Flag==1)  && (ui >= ULastPeak_neg_j_1)) {
				df = du * Kul_j_1;
				//cout << "  Case = 2-" << endln;

			// CASE 2: WHEN UNLOADING
			} else if (Unloading_Flag==1) {
				df = du * Kul_j_1;
				//cout << "  Case = 3-" << endln;

			// CASE 3: WHEN RELOADING BUT BETWEEN LAST CYCLE PEAK POINT AND GLOBAL PEAK POINT
			} else if ((Reloading_Flag==1) && (ui <= ULastPeak_neg_j_1) && (ui >= Upeak_neg_j_1)) {
				Krel_j_1 = (Fpeak_neg_j_1 - FLastPeak_neg_j_1) / (Upeak_neg_j_1 - ULastPeak_neg_j_1);
				if (ui_1 >= ULastPeak_neg_j_1) {
					df = Kul_j_1 * (ULastPeak_neg_j_1 - ui_1) + Krel_j_1 * (ui - ULastPeak_neg_j_1);
				} else {
					df = du * Krel_j_1;
				}
				//cout << "  Case = 4-" << endln;

			// CASE 4: WHEN LOADING IN GENERAL TOWARDS THE TARGET PEAK
			} else if ((du <= 0) && (((TargetPeak_Flag==0) && (ui >= Upeak_neg_j_1)) || ((TargetPeak_Flag==1) && (ui >= ULastPeak_neg_j_1)))) {
				df = du * Krel_j_1;
				//cout << "  Case = 4-" << endln;

			// CASE 5: WHEN LOADING IN GENERAL TOWARDS THE LAST CYCLE PEAK POINT BUT BEYOND IT
			} else if ((du <= 0) && (TargetPeak_Flag==1) && (ui <= ULastPeak_neg_j_1) && (ui >= Upeak_neg_j_1)) {
				Krel_j_1 = (Fpeak_neg_j_1 - FLastPeak_neg_j_1) / (Upeak_neg_j_1 - ULastPeak_neg_j_1);
				if (ui_1 >= ULastPeak_neg_j_1) {
					df = (FLastPeak_neg_j_1 - fi_1) + Krel_j_1 * (ui - ULastPeak_neg_j_1);
				} else {
					df = du * Krel_j_1;
				}
				//cout << "  Case = 5-" << endln;

			// CASE 6: WHEN LOADING BEYOND THE TARGET PEAK BUT BEFORE THE CAPPING POINT
			} else if ((du <= 0) && (ui >= Umax_neg_j_1)) {
				df = du * Kp_neg_j_1;
				//cout << "  Case = 6-" << endln;

			// CASE 7: WHEN LOADING AND BETWEEN THE CAPPING POINT AND THE RESIDUAL POINT
			} else if ((du < 0) && (ui <= Umax_neg_j_1) && (ui >= Ures_neg_j_1)) {
				if ((ui_1>=Umax_neg_j_1) && (ui<=Umax_neg_j_1)) {
					df = Kp_neg_j_1 * (Umax_neg_j_1 - ui_1) + Kpc_neg_j_1 * (ui - Umax_neg_j_1);
				} else {
					df = du * Kpc_neg_j_1;
				}
				//cout << "  Case = 7-" << endln;

			// CASE 8: WHEN LOADING AND BEYOND THE RESIDUAL POINT
			} else if ((du < 0) &&  (ui <= Ures_neg_j_1)) {
				df = 0.0;
				if (Fres_neg_j_1 == 0) {
					Failure_Flag = 1;
				}
				//cout << "  Case = 8-" << endln;

			}
		}
	
	
		// Force
		fi = fi_1 + df;
		//cout << "  Excurion=" << Excursion_Flag << " Failure=" << Failure_Flag << "  Reload=" << Reloading_Flag << " Unload=" << Unloading_Flag << " Yield=" << Yield_Flag << endln;
		//cout << "  STEP: ui_1=" << ui_1 << " ui=" << ui << " fi_1=" << fi_1 << " fi=" << fi << endln;

		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////
		// CHECK FOR FAILURE
		///////////////////////////////////////////////////////////////////////////////////////////		
		///////////////////////////////////////////////////////////////////////////////////////////		
		///////////////////////////////////////////////////////////////////////////////////////////		

		// Failure criteria (Tolerance = 1//)
		FailS = ((betaS < -0.01) || (betaS > 1.01));
		FailC = ((betaC < -0.01) || (betaC > 1.01));
		FailA = ((betaA < -0.01) || (betaA > 1.01));
		FailK = ((betaK < -0.01) || (betaK > 1.01));
		//cout << "  ENERGY: EtS=" << EtS << " EtC=" << EtC << " EtA=" << EtA << " EtK=" << EtK << endln;
		//cout << "  ENERGY: dEi=" << dEi << " Ei=" << Ei << " Energy_Diss=" << Energy_Diss << " Energy_Acc=" << Energy_Acc << endln;
		//cout << "  ENERGY: betaS=" << betaS << " betaC=" << betaC << " betaA=" << betaA << " betaK=" << betaK << endln;
		//cout << "  FAIL:   FailS=" << FailS << " FailC=" << FailC << " FailA=" << FailA << " FailK=" << FailK << endln;

		if (FailS || FailC || FailA || FailK) {
			fi = 0;
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

		dEi = 0.5*(fi + fi_1)*du; // Internal energy increment

	}
	else {
		fi = 0;
		dEi = 0;
		//cout << "  FAILURE OCCURED" << endln;
	}


	//// Energy
	Energy_Acc = Energy_Acc + dEi; 	

	//// Update Variables
	du_i_1 = du;

	// Tangent Stiffeness Calculation
	if (fi == fi_1) {
		TangentK = pow(10., -6);
		ki		 = pow(10., -6);
	}	
	
	if (ui == ui_1) {
		ki		 = Ke;
		fi		 = fi_1;
		TangentK = Ke;
	}
	else {
		ki		 = (fi - fi_1) / (du);
		TangentK = (fi - fi_1) / (du);
	}

	//cout << "  fi=" << fi << endln;
	//cout << "***********************" << endln;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////// END OF MAIN CODE ///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return 0;
}

double IMKPeakOriented::getStress(void)
{
	//cout << " getStress" << endln;
	return (fi);
}

double IMKPeakOriented::getTangent(void)
{
	//cout << " getTangent" << endln;
	return (TangentK);
}

double IMKPeakOriented::getInitialTangent(void)
{
	//cout << " getInitialTangent" << endln;
	return (Ke);
}

double IMKPeakOriented::getStrain(void)
{
	//cout << " getStrain" << endln;
	return (U);
}

int IMKPeakOriented::commitState(void)
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
	cFy_pos_j_1 = Fy_pos_j_1;
	cFmax_pos_j_1 = Fmax_pos_j_1;
	cUpeak_pos_j_1 = Upeak_pos_j_1;
	cFpeak_pos_j_1 = Fpeak_pos_j_1;

	cUres_pos_j_1 = Ures_pos_j_1;
	cFres_pos_j_1 = Fres_pos_j_1;
	cKp_pos_j_1 = Kp_pos_j_1;
	cKpc_pos_j_1 = Kpc_pos_j_1;

	cUy_neg_j_1 = Uy_neg_j_1;
	cUmax_neg_j_1 = Umax_neg_j_1;
	cFy_neg_j_1 = Fy_neg_j_1;
	cFmax_neg_j_1 = Fmax_neg_j_1;
	cUpeak_neg_j_1 = Upeak_neg_j_1;
	cFpeak_neg_j_1 = Fpeak_neg_j_1;

	cUres_neg_j_1 = Ures_neg_j_1;
	cFres_neg_j_1 = Fres_neg_j_1;
	cKp_neg_j_1 = Kp_neg_j_1;
	cKpc_neg_j_1 = Kpc_neg_j_1;

	cKul_j_1 = Kul_j_1;

	cEnergy_Acc = Energy_Acc;
	cEnergy_Diss = Energy_Diss;

	cu0 = u0;
	
	cULastPeak_pos_j_1 = ULastPeak_pos_j_1;
	cFLastPeak_pos_j_1 = FLastPeak_pos_j_1;
	cULastPeak_neg_j_1 = ULastPeak_neg_j_1;
	cFLastPeak_neg_j_1 = FLastPeak_neg_j_1;

	cFailure_Flag		= Failure_Flag;
	cExcursion_Flag		= Excursion_Flag;
	cReloading_Flag		= Reloading_Flag;
	cUnloading_Flag		= Unloading_Flag;
	cTargetPeak_Flag	= TargetPeak_Flag;
	cYield_Flag			= Yield_Flag;
	cReversal_Flag = Reversal_Flag;

	cKrel_j_1 = Krel_j_1;

	return 0;
}

int IMKPeakOriented::revertToLastCommit(void)
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
	Fy_pos_j_1 = cFy_pos_j_1;
	Fmax_pos_j_1 = cFmax_pos_j_1;
	Upeak_pos_j_1 = cUpeak_pos_j_1;
	Fpeak_pos_j_1 = cFpeak_pos_j_1;

	Ures_pos_j_1 = cUres_pos_j_1;
	Fres_pos_j_1 = cFres_pos_j_1;
	Kp_pos_j_1 = cKp_pos_j_1;
	Kpc_pos_j_1 = cKpc_pos_j_1;


	Uy_neg_j_1 = cUy_neg_j_1;
	Umax_neg_j_1 = cUmax_neg_j_1;
	Fy_neg_j_1 = cFy_neg_j_1;
	Fmax_neg_j_1 = cFmax_neg_j_1;
	Upeak_neg_j_1 = cUpeak_neg_j_1;
	Fpeak_neg_j_1 = cFpeak_neg_j_1;

	Ures_neg_j_1 = cUres_neg_j_1;
	Fres_neg_j_1 = cFres_neg_j_1;
	Kp_neg_j_1 = cKp_neg_j_1;
	Kpc_neg_j_1 = cKpc_neg_j_1;


	Kul_j_1 = cKul_j_1;



	Energy_Acc = cEnergy_Acc;
	Energy_Diss = cEnergy_Diss;

	ULastPeak_pos_j_1 = cULastPeak_pos_j_1;
	FLastPeak_pos_j_1 = cFLastPeak_pos_j_1;
	ULastPeak_neg_j_1 = cULastPeak_neg_j_1;
	FLastPeak_neg_j_1 = cFLastPeak_neg_j_1;

	Failure_Flag = cFailure_Flag;
	Excursion_Flag = cExcursion_Flag;
	Reloading_Flag    = cReloading_Flag;
	Unloading_Flag    = cUnloading_Flag;
	TargetPeak_Flag   = cTargetPeak_Flag;
	Yield_Flag   = cYield_Flag;
	Reversal_Flag = cReversal_Flag;

	u0 = cu0;

	Krel_j_1 = cKrel_j_1;
	
	return 0;
}

int IMKPeakOriented::revertToStart(void)
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\\
	//////////////////////////////////////////////////////////////////// ONE TIME CALCULATIONS ////////////////////////////////////////////////////////////////////\\
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/



	betaS = 0;
	betaC = 0;
	betaK = 0;
	betaA = 0;

	Uy_pos   = Fy_pos / Ke;
	Umax_pos = Uy_pos + Up_pos;
	Fmax_pos = FmaxFy_pos*Fy_pos;
	Kp_pos 	 = (Fmax_pos - Fy_pos) / Up_pos;
	Kpc_pos  = Fmax_pos / Upc_pos;

	Uy_neg 	 = Fy_neg / Ke;
	Umax_neg = Uy_neg + Up_neg;
	Fmax_neg = FmaxFy_neg*Fy_neg;
	Kp_neg 	 = (Fmax_neg - Fy_neg) / Up_neg;
	Kpc_neg  = Fmax_neg / Upc_neg;

	Upeak_pos_j_1 = Uy_pos;
	Fpeak_pos_j_1 = Fy_pos;
	Upeak_neg_j_1 = -Uy_neg;
	Fpeak_neg_j_1 = -Fy_neg;
	
	Uy_pos_j_1 	=  Uy_pos;
	Fy_pos_j_1  =  Fy_pos;
	Kp_pos_j_1  =  Kp_pos;
	Kpc_pos_j_1 = -Kpc_pos;
	Uy_neg_j_1 	= -Uy_neg;
	Fy_neg_j_1 	= -Fy_neg;
	Kp_neg_j_1 	=  Kp_neg;
	Kpc_neg_j_1 = -Kpc_neg;

	Umax_pos_j_1 = Umax_pos;
	Fmax_pos_j_1 = Fmax_pos;
	Fres_pos_j_1 = Fy_pos*ResF_pos;
	Umax_neg_j_1 = -Umax_neg;
	Fmax_neg_j_1 = -Fmax_neg;
	Fres_neg_j_1 = -Fy_neg*ResF_neg;

	Kul_j_1 = Ke;

	Ures_pos_j_1 = (Fres_pos_j_1 - Fmax_pos_j_1) / Kpc_pos_j_1 + Umax_pos_j_1;
	Ures_neg_j_1 = (Fres_neg_j_1 - Fmax_neg_j_1) / Kpc_neg_j_1 + Umax_neg_j_1;

	Energy_Acc  = cEnergy_Acc = 0.0;
	Energy_Diss = cEnergy_Diss = 0.0;

	u0 = 0.0;

	EtS = LAMBDA_S *Fy_pos;
	EtC = LAMBDA_C *Fy_pos;
	EtA = LAMBDA_A *Fy_pos;
	EtK = LAMBDA_K *Fy_pos;

	Failure_Flag 	= cFailure_Flag	   = 0;
	Excursion_Flag 	= cExcursion_Flag  = 0;
	Unloading_Flag 	= cUnloading_Flag  = 0;
	Reloading_Flag	= cReloading_Flag  = 0;
	TargetPeak_Flag = cTargetPeak_Flag = 0;
	Yield_Flag		= cYield_Flag	   = 0;
	Reversal_Flag	= cReversal_Flag   = 0;

	ULastPeak_pos_j_1 =  Uy_pos;
	FLastPeak_pos_j_1 =  Fy_pos;
	ULastPeak_neg_j_1 = -Uy_neg;
	FLastPeak_neg_j_1 = -Fy_neg;
	
	Krel_j_1 	  = Ke;

	cdu_i_1 = 0;

	cUpeak_pos_j_1 =  Uy_pos;
	cFpeak_pos_j_1 =  Fy_pos;
	cUpeak_neg_j_1 = -Uy_neg;
	cFpeak_neg_j_1 = -Fy_neg;

	cUy_pos_j_1 =	Uy_pos;
	cFy_pos_j_1 =	Fy_pos;
	cKp_pos_j_1 =	Kp_pos;
	cKpc_pos_j_1 = -Kpc_pos;
	cUy_neg_j_1 =  -Uy_neg;
	cFy_neg_j_1 =  -Fy_neg;
	cKp_neg_j_1 =	Kp_neg;
	cKpc_neg_j_1 = -Kpc_neg;

	cUmax_pos_j_1 =  Umax_pos;
	cFmax_pos_j_1 =  Fmax_pos;
	cFres_pos_j_1 =  Fy_pos*ResF_pos;
	cUmax_neg_j_1 = -Umax_neg;
	cFmax_neg_j_1 = -Fmax_neg;
	cFres_neg_j_1 = -Fy_neg*ResF_neg;

	cKul_j_1 = Ke;

	cUres_pos_j_1 = (Fres_pos_j_1 - Fmax_pos_j_1) / Kpc_pos_j_1 + Umax_pos_j_1;
	cUres_neg_j_1 = (Fres_neg_j_1 - Fmax_neg_j_1) / Kpc_neg_j_1 + Umax_neg_j_1;

	cULastPeak_pos_j_1 =  Uy_pos;
	cFLastPeak_pos_j_1 =  Fy_pos;
	cULastPeak_neg_j_1 = -Uy_neg;
	cFLastPeak_neg_j_1 = -Fy_neg;
	
	cKrel_j_1  	   = Ke;
	
	//initially I zero everything   
	U = cU = 0;
	ui	  = 0;
	fi	  = 0;
	ui_1  = 0;
	fi_1  = 0;
	cui   = 0;
	cfi	  = 0;
	cui_1 = 0;
	cfi_1 = 0;

	 TangentK = Ke;
	cTangentK = Ke;
	//cout << " revertToStart:" << endln; //<< " U=" << U << " Ri=" << Ri << " TanK=" << TangentK << endln;

	return 0;
}

UniaxialMaterial *
IMKPeakOriented::getCopy(void)
{
	IMKPeakOriented *theCopy = new IMKPeakOriented(this->getTag(), Ke,
		Uy_pos, Umax_pos, Uu_pos, Fy_pos, FmaxFy_pos, ResF_pos,
		Uy_neg, Umax_neg, Uu_neg, Fy_neg, FmaxFy_neg, ResF_neg,
		LAMBDA_S, LAMBDA_C, LAMBDA_A, LAMBDA_K, c_S, c_C, c_A, c_K, D_pos, D_neg);

	//cout << " getCopy" << endln;

	theCopy->U	= U;
	theCopy->cU = cU;

	theCopy->TangentK = TangentK;

	theCopy->ui		= ui;
	theCopy->fi		= fi;
	theCopy->ui_1	= ui_1;
	theCopy->fi_1	= fi_1;
	theCopy->du_i_1 = du_i_1;

	theCopy->Uy_pos_j_1		= Uy_pos_j_1;
	theCopy->Umax_pos_j_1	= Umax_pos_j_1;
	theCopy->Fy_pos_j_1		= Fy_pos_j_1;
	theCopy->Fmax_pos_j_1	= Fmax_pos_j_1;
	theCopy->Upeak_pos_j_1	= Upeak_pos_j_1;
	theCopy->Fpeak_pos_j_1	= Fpeak_pos_j_1;
	theCopy->Ures_pos_j_1	= Ures_pos_j_1;
	theCopy->Fres_pos_j_1	= Fres_pos_j_1;
	theCopy->Kp_pos_j_1		= Kp_pos_j_1;
	theCopy->Kpc_pos_j_1	= Kpc_pos_j_1;

	theCopy->Uy_neg_j_1		= Uy_neg_j_1;
	theCopy->Umax_neg_j_1	= Umax_neg_j_1;
	theCopy->Fy_neg_j_1		= Fy_neg_j_1;
	theCopy->Fmax_neg_j_1	= Fmax_neg_j_1;
	theCopy->Upeak_neg_j_1	= Upeak_neg_j_1;
	theCopy->Fpeak_neg_j_1	= Fpeak_neg_j_1;
	theCopy->Ures_neg_j_1	= Ures_neg_j_1;
	theCopy->Fres_neg_j_1	= Fres_neg_j_1;
	theCopy->Kp_neg_j_1		= Kp_neg_j_1;
	theCopy->Kpc_neg_j_1	= Kpc_neg_j_1;

	theCopy->Kul_j_1 = Kul_j_1;

	theCopy->Energy_Acc	 = Energy_Acc;
	theCopy->Energy_Diss = Energy_Diss;

	theCopy->u0 = u0;

	theCopy->ULastPeak_pos_j_1 = ULastPeak_pos_j_1;
	theCopy->FLastPeak_pos_j_1 = FLastPeak_pos_j_1;
	theCopy->ULastPeak_neg_j_1 = ULastPeak_neg_j_1;
	theCopy->FLastPeak_neg_j_1 = FLastPeak_neg_j_1;

	theCopy->Failure_Flag 	 = Failure_Flag;
	theCopy->Excursion_Flag  = Excursion_Flag;
	theCopy->Reloading_Flag  = Reloading_Flag;
	theCopy->TargetPeak_Flag = TargetPeak_Flag;
	theCopy->Yield_Flag 	 = Yield_Flag;
	theCopy->Reversal_Flag	 = Reversal_Flag;

	theCopy->Krel_j_1 = Krel_j_1;

	theCopy->cTangentK = cTangentK;

	theCopy->cui = cui;
	theCopy->cfi = cfi;
	theCopy->cui_1 = cui_1;
	theCopy->cfi_1 = cfi_1;
	theCopy->cdu_i_1 = cdu_i_1;

	theCopy->cUy_pos_j_1	= cUy_pos_j_1;
	theCopy->cUmax_pos_j_1	= cUmax_pos_j_1;
	theCopy->cFy_pos_j_1	= cFy_pos_j_1;
	theCopy->cFmax_pos_j_1	= cFmax_pos_j_1;
	theCopy->cUpeak_pos_j_1 = cUpeak_pos_j_1;
	theCopy->cFpeak_pos_j_1 = cFpeak_pos_j_1;
	theCopy->cUres_pos_j_1	= cUres_pos_j_1;
	theCopy->cFres_pos_j_1	= cFres_pos_j_1;
	theCopy->cKp_pos_j_1	= cKp_pos_j_1;
	theCopy->cKpc_pos_j_1	= cKpc_pos_j_1;

	theCopy->cUy_neg_j_1	= cUy_neg_j_1;
	theCopy->cUmax_neg_j_1	= cUmax_neg_j_1;
	theCopy->cFy_neg_j_1	= cFy_neg_j_1;
	theCopy->cFmax_neg_j_1	= cFmax_neg_j_1;
	theCopy->cUpeak_neg_j_1 = cUpeak_neg_j_1;
	theCopy->cFpeak_neg_j_1 = cFpeak_neg_j_1;
	theCopy->cUres_neg_j_1	= cUres_neg_j_1;
	theCopy->cFres_neg_j_1	= cFres_neg_j_1;
	theCopy->cKp_neg_j_1	= cKp_neg_j_1;
	theCopy->cKpc_neg_j_1	= cKpc_neg_j_1;

	theCopy->cKul_j_1 = cKul_j_1;

	theCopy->cEnergy_Acc  = cEnergy_Acc;
	theCopy->cEnergy_Diss = cEnergy_Diss;

	theCopy->cu0 = cu0;
	
	theCopy->cULastPeak_pos_j_1 = cULastPeak_pos_j_1;
	theCopy->cFLastPeak_pos_j_1 = cFLastPeak_pos_j_1;
	theCopy->cULastPeak_neg_j_1 = cULastPeak_neg_j_1;
	theCopy->cFLastPeak_neg_j_1 = cFLastPeak_neg_j_1;	

	theCopy->cFailure_Flag		= cFailure_Flag;
	theCopy->cExcursion_Flag	= cExcursion_Flag;	
	theCopy->cReloading_Flag	= cReloading_Flag;
	theCopy->cTargetPeak_Flag	= cTargetPeak_Flag;
	theCopy->cYield_Flag 		= cYield_Flag;
	theCopy->cReversal_Flag		= cReversal_Flag;

	theCopy->cKrel_j_1 = cKrel_j_1;
	
	return theCopy;
}

int IMKPeakOriented::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	cout << " sendSelf" << endln;

	static Vector data(137);
	data(0) = this->getTag();
	data(1)   = Ke;
	data(2)   = Uy_pos;
	data(3)   = Umax_pos;
	data(4)   = Uu_pos;
	data(5)   = Fy_pos;
	data(6)   = FmaxFy_pos;
	data(7)   = ResF_pos;
	data(8)   = Uy_neg;
	data(9)   = Umax_neg;
	data(10)  = Uu_neg;
	data(11)  = Fy_neg;
	data(12)  = FmaxFy_neg;
	data(13)  = ResF_neg;
	data(14)  = LAMBDA_S;
	data(15)  = LAMBDA_C;
	data(16)  = LAMBDA_A;
	data(17)  = LAMBDA_K;
	data(18)  = c_S;
	data(19)  = c_C;
	data(20)  = c_A;
	data(21)  = c_K;
	data(22)  = D_pos;
	data(23)  = D_neg;
	data(24)  = ui;
	data(25)  = fi;
	data(26)  = ui_1;
	data(27)  = fi_1;
	data(28)  = du_i_1;
			  
	data(29)  = Uy_pos_j_1;
	data(30)  = Umax_pos_j_1;
	data(31)  = Fy_pos_j_1;
	data(32)  = Fmax_pos_j_1;
	data(33)  = Upeak_pos_j_1;
	data(34)  = Fpeak_pos_j_1;
	data(35)  = Ures_pos_j_1;
	data(36)  = Fres_pos_j_1;
	data(37)  = Kp_pos_j_1;
	data(38)  = Kpc_pos_j_1;

	data(39)  = Uy_neg_j_1;
	data(40)  = Umax_neg_j_1;
	data(41)  = Fy_neg_j_1;
	data(42)  = Fmax_neg_j_1;
	data(43)  = Upeak_neg_j_1;
	data(44)  = Fpeak_neg_j_1;
	data(45)  = Ures_neg_j_1;
	data(46)  = Fres_neg_j_1;
	data(47)  = Kp_neg_j_1;
	data(48)  = Kpc_neg_j_1;
			  
	data(49)  = Kul_j_1;

	data(50)  = Failure_Flag;
	data(51)  = Excursion_Flag;
	data(52)  = Unloading_Flag;
	data(53)  = Reloading_Flag;
	data(54)  = TargetPeak_Flag;
	data(55)  = Yield_Flag;

	data(56)  = Energy_Acc;
	data(57)  = Energy_Diss;

	data(58)  = u0;
	data(59)  = du;
	data(60)  = df;

	data(61)  = FailS;
	data(62)  = FailC;
	data(63)  = FailA;
	data(64)  = FailK;

	data(65) = Ei;
	data(66) = dEi;
	data(67) = Epj;
	data(68) = EpjK;
	data(69) = EiK;
	data(70) = c_S;
	data(71) = c_C;
	data(72) = c_A;
	data(73) = c_K;
	data(74) = EtS;
	data(75) = EtC;
	data(76) = EtA;
	data(77) = EtK;
	data(78) = betaS;
	data(79) = betaC;
	data(80) = betaA;
	data(81) = betaK;
	data(82) = sPCsp;
	data(83) = sPCpcp;

	data(84) = TangentK;

	data(85) = Uy_pos;
	data(86) = Umax_pos;
	data(87) = Fmax_pos;
	data(88) = Kp_pos;
	data(89) = Kpc_pos;

	data(90) = Uy_neg;
	data(91) = Umax_neg;
	data(92) = Fmax_neg;
	data(93) = Kp_neg;
	data(94) = Kpc_neg;

	data(95) = cui;
	data(96) = cfi;
	data(97) = cui_1;
	data(98) = cfi_1;
	data(99) = cdu_i_1;

	data(100) = cUy_pos_j_1;
	data(101) = cUmax_pos_j_1;
	data(102) = cFy_pos_j_1;
	data(103) = cFmax_pos_j_1;
	data(104) = cUpeak_pos_j_1;
	data(105) = cFpeak_pos_j_1;
	data(106) = cUres_pos_j_1;
	data(107) = cFres_pos_j_1;
	data(108) = cKp_pos_j_1;
	data(109) = cKpc_pos_j_1;

	data(110) = cUy_neg_j_1;
	data(111) = cUmax_neg_j_1;
	data(112) = cFy_neg_j_1;
	data(113) = cFmax_neg_j_1;
	data(114) = cUpeak_neg_j_1;
	data(115) = cFpeak_neg_j_1;
	data(116) = cUres_neg_j_1;
	data(117) = cFres_neg_j_1;
	data(118) = cKp_neg_j_1;
	data(119) = cKpc_neg_j_1;

	data(120) = cKul_j_1;

	data(121) = cULastPeak_pos_j_1;
	data(122) = cFLastPeak_pos_j_1;
	data(123) = cULastPeak_neg_j_1;
	data(124) = cFLastPeak_neg_j_1;

	data(125) = cFailure_Flag;
	data(126) = cExcursion_Flag;
	data(127) = cReloading_Flag;
	data(128) = cUnloading_Flag;
	data(129) = cTargetPeak_Flag;
	data(130) = cYield_Flag;

	data(131) = cKrel_j_1;

	data(132) = Krel_LastPeak;
	data(133) = Krel_GlobalPeak;
	data(134) = K_check;

	data(135) = cReversal_Flag;
	data(136) = Reversal_Flag;

	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "IMKPeakOriented::sendSelf() - failed to send data\n";

	return res;
}

int IMKPeakOriented::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector data(137);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);

	if (res < 0) {
		opserr << "IMKPeakOriented::recvSelf() - failed to receive data\n";
		this->setTag(0);
	}
	else {
		cout << " recvSelf" << endln;
		this->setTag((int)data(0));
		Ke					= data(1);
		Up_pos				= data(2);
		Upc_pos				= data(3);
		Uu_pos				= data(4);
		Fy_pos				= data(5);
		FmaxFy_pos			= data(6);
		ResF_pos			= data(7);
		Up_neg				= data(8);
		Upc_neg				= data(9);
		Uu_neg				= data(10);
		Fy_neg				= data(11);
		FmaxFy_neg			= data(12);
		ResF_neg			= data(13);
		LAMBDA_S			= data(14);
		LAMBDA_C			= data(15);
		LAMBDA_A			= data(16);
		LAMBDA_K			= data(17);
		c_S					= data(18);
		c_C					= data(19);
		c_A					= data(20);
		c_K					= data(21);
		D_pos				= data(22);
		D_neg				= data(23);
		ui					= data(24);
		fi					= data(25);
		ui_1				= data(26);
		fi_1				= data(27);
		du_i_1				= data(28);
		Uy_pos_j_1			= data(29);
		Umax_pos_j_1		= data(30);
		Fy_pos_j_1			= data(31);
		Fmax_pos_j_1		= data(32);
		Upeak_pos_j_1		= data(33);
		Fpeak_pos_j_1		= data(34);
		Ures_pos_j_1		= data(35);
		Fres_pos_j_1		= data(36);
		Kp_pos_j_1			= data(37);
		Kpc_pos_j_1			= data(38);
		Uy_neg_j_1			= data(39);
		Umax_neg_j_1		= data(40);
		Fy_neg_j_1			= data(41);
		Fmax_neg_j_1		= data(42);
		Upeak_neg_j_1		= data(43);
		Fpeak_neg_j_1		= data(44);
		Ures_neg_j_1		= data(45);
		Fres_neg_j_1		= data(46);
		Kp_neg_j_1			= data(47);
		Kpc_neg_j_1			= data(48);
		Failure_Flag		= data(49);
		Excursion_Flag		= data(50);
		Reloading_Flag		= data(51);
		Unloading_Flag		= data(52);
		TargetPeak_Flag		= data(53);
		Yield_Flag			= data(54);
		Kul_j_1				= data(55);
		Energy_Acc			= data(56);
		Energy_Diss			= data(57);
		u0					= data(58);
		du					= data(59);
		df					= data(60);
		FailS				= data(61);
		FailC				= data(62);
		FailA				= data(63);
		FailK				= data(64);
		Ei					= data(65);
		dEi					= data(66);
		Epj					= data(67);
		EpjK				= data(68);
		EiK					= data(79);
		c_S					= data(70);
		c_C					= data(71);
		c_A					= data(72);
		c_K					= data(73);
		EtS					= data(74);
		EtC					= data(75);
		EtA					= data(76);
		EtK					= data(77);
		betaS				= data(78);
		betaC				= data(79);
		betaA				= data(80);
		betaK				= data(81);
		sPCsp				= data(82);
		sPCpcp				= data(83);
		TangentK			= data(84);
		Uy_pos				= data(85);
		Umax_pos			= data(86);
		Fmax_pos			= data(87);
		Kp_pos				= data(88);
		Kpc_pos				= data(89);
		Uy_neg				= data(90);
		Umax_neg			= data(91);
		Fmax_neg			= data(92);
		Kp_neg				= data(93);
		Kpc_neg				= data(94);
		cui					= data(95);
		cfi					= data(96);
		cui_1				= data(97);
		cfi_1				= data(98);
		cdu_i_1				= data(99);
		cUy_pos_j_1			= data(100);
		cUmax_pos_j_1		= data(101);
		cFy_pos_j_1			= data(102);
		cFmax_pos_j_1		= data(103);
		cUpeak_pos_j_1		= data(104);
		cFpeak_pos_j_1		= data(105);
		cUres_pos_j_1		= data(106);
		cFres_pos_j_1		= data(107);
		cKp_pos_j_1			= data(108);
		cKpc_pos_j_1		= data(109);
		cUy_neg_j_1			= data(110);
		cUmax_neg_j_1		= data(111);
		cFy_neg_j_1			= data(112);
		cFmax_neg_j_1		= data(113);
		cUpeak_neg_j_1		= data(114);
		cFpeak_neg_j_1		= data(115);
		cUres_neg_j_1		= data(116);
		cFres_neg_j_1		= data(117);
		cKp_neg_j_1			= data(118);
		cKpc_neg_j_1		= data(119);
		cKul_j_1			= data(120);
		cULastPeak_pos_j_1	= data(121);
		cFLastPeak_pos_j_1	= data(122);
		cULastPeak_neg_j_1	= data(123);
		cFLastPeak_neg_j_1	= data(124);
		cFailure_Flag		= data(125);
		cExcursion_Flag		= data(126);
		cReloading_Flag		= data(127);
		cUnloading_Flag		= data(128);
		cTargetPeak_Flag	= data(129);
		cYield_Flag   		= data(130);
		cKrel_j_1      		= data(131);
		 Krel_LastPeak		= data(132);
		 Krel_GlobalPeak	= data(133);
		 K_check			= data(134);
		 cReversal_Flag = data(135);
		 Reversal_Flag = data(136);
	}

	return res;
}

void IMKPeakOriented::Print(OPS_Stream &s, int flag)
{
	cout << "IMKPeakOriented tag: " << this->getTag() << endln;
}