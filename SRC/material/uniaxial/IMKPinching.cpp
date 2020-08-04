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
OPS_IMKPinching()
{
	if (numIMKPinchingMaterials == 0) {
		numIMKPinchingMaterials++;
		OPS_Error("IMK Model with Pinched Response - Code by A. ELKADY & H. ELJISR (July 2020)\n", 1);
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
		opserr << "LamdaS? LamdaC? LamdaA? LamdaK? Cs? Cc? Ca? Ck? D_pos? D_neg? kappaF? kappaD? ";
		return 0;
	}



	// Parsing was successful, allocate the material
	theMaterial = new IMKPinching(iData[0],
		dData[0],
		dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
		dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],
		dData[13], dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20],
		dData[21], dData[22],dData[23], dData[24]);

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
			K_check = (FLastPeak_pos_j_1 - fi_1) / (ULastPeak_pos_j_1 - ui_1);
			if ((K_check >= 1.05*Kul_j_1) || (K_check <= 0.95*Kul_j_1)) { // a tailored criteria to avoid registering last peak points during small unload/reload excursions on the unloading branch 
				FLastPeak_pos_j_1 = fi_1;
				ULastPeak_pos_j_1 = ui_1;
			}
		}
		else if ((fi_1 < 0) && (du > 0) && (du*du_i_1 <= 0)) {
			Unloading_Flag = 1;
			Reversal_Flag = 1;
			Reloading_Flag = 0;
			K_check = (FLastPeak_neg_j_1 - fi_1) / (ULastPeak_neg_j_1 - ui_1);
			if ((K_check >= 1.01*Kul_j_1) || (K_check <= 0.99*Kul_j_1)) {
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
		} 
		else if ((fi_1 < 0) && (ui_1 <= Upeak_neg_j_1)) {
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

		if (Excursion_Flag == 1) {
			Ei = fmax(0,Energy_Acc - Energy_Diss);
			betaS = pow((Ei / (EtS - Energy_Acc)), c_S);
			betaC = pow((Ei / (EtC - Energy_Acc)), c_C);
			betaA = pow((Ei / (EtA - Energy_Acc)), c_A);
			Energy_Diss = Energy_Acc;
		}
		else {
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
		
		// CHECK TARGET POINT: LAST CYCLE PEAK or GLOBAL PEAK (i.e., Modified Clough rule, see Mahin & Bertero 1975)
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
			}
			else if (abs(Krel_LastPeak) <= abs(Krel_GlobalPeak)) {
				TargetPeak_Flag = 0;
			}
			else if ((du >= 0) && (abs((ULastPeak_pos_j_1 - Upeak_pos_j_1) / Upeak_pos_j_1) < 0.05) && (abs(Krel_LastPeak) <= 1.05*abs(Krel_GlobalPeak))) {
				TargetPeak_Flag = 0;
			}
			else if ((du <= 0) && (abs((ULastPeak_neg_j_1 - Upeak_neg_j_1) / Upeak_neg_j_1) < 0.05) && (abs(Krel_LastPeak) <= 1.05*abs(Krel_GlobalPeak))) {
				TargetPeak_Flag = 0;
			}
			else {
				TargetPeak_Flag = 1;
			}
		}
		
		// COMPUTE PINCHING POINT COORDINATES AT EACH NEW EXCURSION
		if (Excursion_Flag==1) {
			if (du>0) {
				Upl = Upeak_pos_j_1 - (Fpeak_pos_j_1/Kul_j_1);
				Ubp = (1-kappaD) * Upl ;
				Fbp = kappaF * Fpeak_pos_j_1 * abs((Ubp -u0)/(Upeak_pos_j_1 -u0));
			} 
			else if (du<0) {
				Upl = Upeak_neg_j_1 - (Fpeak_neg_j_1/Kul_j_1);
				Ubp = (1-kappaD) * Upl ;
				Fbp = kappaF * Fpeak_neg_j_1 * abs((Ubp -u0)/(Upeak_neg_j_1 -u0));
			}
		}
		//cout << "  Upl ="  << Upl << "  Ubp =" << Ubp << "  Fbp =" << Fbp << "  kF =" << kappaF << "  kD =" << kappaD << endln;

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
				if ((TargetPeak_Flag==0) && (u0>=Ubp)) {
					Krel_j_1 =  Fpeak_pos_j_1 / (Upeak_pos_j_1 - u0);
				}
				else if ((TargetPeak_Flag==0) && (u0<=Ubp)) {
					Krel_j_1 =  Fbp / (Ubp - u0);
				}
				else if (TargetPeak_Flag==1) {
					Krel_j_1 = FLastPeak_pos_j_1 / (ULastPeak_pos_j_1 - u0);
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
			} else if ((Reloading_Flag==1) && (ui >= ULastPeak_pos_j_1) && (ui <= Upeak_pos_j_1) && (FLastPeak_pos_j_1 <= Fpeak_pos_j_1)) {
				if (TargetPeak_Flag==1) {
					if ((FLastPeak_pos_j_1 <= Fbp) && (ui <= Ubp)) {
						Krel_j_1  = (Fbp-FLastPeak_pos_j_1)/(Ubp-ULastPeak_pos_j_1); 
					} else if ((FLastPeak_pos_j_1 <= Fbp) && (ui >= Ubp)) {
						Krel_j_1 = (Fpeak_pos_j_1-Fbp)/(Upeak_pos_j_1-Ubp);
					} else {
						Krel_j_1 = (Fpeak_pos_j_1-FLastPeak_pos_j_1)/(Upeak_pos_j_1-ULastPeak_pos_j_1);
					}
					df = du * Krel_j_1;
				}
				else if (TargetPeak_Flag==0) {
					Krel_j_1 =  Fbp / (Ubp - u0);
					if (ui_1 <= ULastPeak_pos_j_1) {
						df = Kul_j_1*(ULastPeak_pos_j_1 -ui_1) + Krel_j_1 *(ui -ULastPeak_pos_j_1);
					} else if (ui<= Ubp) {
						df = du * Krel_j_1;
					} else if (ui>= Ubp) {
						Krel_j_1 =  (Fpeak_pos_j_1 - Fbp) / (Upeak_pos_j_1 - Ubp);
						df = du * Krel_j_1;
					}
				}
				//cout << "  Case = 3+" << endln;

			// CASE 4: WHEN LOADING IN GENERAL TOWARDS THE TARGET PEAK
			} else if ((du >= 0) && (ui <= Upeak_pos_j_1)) {
				if ((TargetPeak_Flag==0) && (ui>=Ubp)) {
					Krel_j_1 =  (Fpeak_pos_j_1 - fi_1) / (Upeak_pos_j_1 - ui_1);
				}
				else if ((TargetPeak_Flag==0) && (ui<=Ubp)) { 
					Krel_j_1 = (Fbp) / (Ubp - u0);
				}
				else if ((TargetPeak_Flag==1) && (ui<=ULastPeak_pos_j_1)) { 
					Krel_j_1 = (FLastPeak_pos_j_1) / (ULastPeak_pos_j_1 - u0);
				}
				else if ((TargetPeak_Flag==1) && (ui>=ULastPeak_pos_j_1)) { 
					if ((FLastPeak_pos_j_1 <= Fbp) && (ui <= Ubp)) {
						Krel_j_1 = (Fbp-FLastPeak_pos_j_1)/(Ubp-ULastPeak_pos_j_1);
					} else if ((FLastPeak_pos_j_1 <= Fbp) && (ui >= Ubp)) {
						Krel_j_1 = (Fpeak_pos_j_1-Fbp)/(Upeak_pos_j_1-Ubp);
					} else {
						Krel_j_1 = (Fpeak_pos_j_1-FLastPeak_pos_j_1)/(Upeak_pos_j_1-ULastPeak_pos_j_1);
					}
				}
				df = du * Krel_j_1;
				//cout << "  Case = 4+" << endln;

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
				if ((TargetPeak_Flag==0) && (u0<=Ubp)) {
					Krel_j_1 =  Fpeak_neg_j_1 / (Upeak_neg_j_1 - u0);
				}
				else if ((TargetPeak_Flag==0) && (u0>=Ubp)) {
					Krel_j_1 =  Fbp / (Ubp - u0);
				}
				else if (TargetPeak_Flag==1) {
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
				//cout << "  Case = 2-" << endln;

			// CASE 3: WHEN RELOADING BUT BETWEEN LAST CYCLE PEAK POINT AND GLOBAL PEAK POINT
			} else if ((Reloading_Flag==1) && (ui <= ULastPeak_neg_j_1) && (ui >= Upeak_neg_j_1)  && (FLastPeak_neg_j_1 >= Fpeak_neg_j_1)) {
				if (TargetPeak_Flag==1)  {
					if ((FLastPeak_neg_j_1 >= Fbp) && (ui <= Ubp)) {
						Krel_j_1 = (Fbp-FLastPeak_neg_j_1)/(Ubp-ULastPeak_neg_j_1);
					} else if ((FLastPeak_neg_j_1 >= Fbp) && (ui <= Ubp)) {
						Krel_j_1 = (Fpeak_neg_j_1-Fbp)/(Upeak_neg_j_1-Ubp);
					} else {
						Krel_j_1 = (Fpeak_neg_j_1-FLastPeak_neg_j_1)/(Upeak_neg_j_1-ULastPeak_neg_j_1);
					}
					df = du * Krel_j_1;
				}
				else if (TargetPeak_Flag==0) {
					Krel_j_1 =  Fbp / (Ubp - u0);
					if (ui_1 >= ULastPeak_neg_j_1) {
						df = Kul_j_1*(ULastPeak_neg_j_1 -ui_1) + Krel_j_1 *(ui -ULastPeak_neg_j_1);
					} else if (ui>= Ubp) {
						df = du * Krel_j_1;
					} else if (ui<= Ubp) {
						Krel_j_1 =  (Fpeak_neg_j_1 - Fbp) / (Upeak_neg_j_1 - Ubp);
						df = du * Krel_j_1;
					}
				}
				//cout << "  Case = 3-" << endln;

			// CASE 4: WHEN LOADING IN GENERAL TOWARDS THE TARGET PEAK
			} else if ((du <= 0) && (ui >= Upeak_neg_j_1)) {
				if ((TargetPeak_Flag==0) && (ui<=Ubp)) {
					Krel_j_1 =  (Fpeak_neg_j_1 - fi_1) / (Upeak_neg_j_1 - ui_1);
				}
				else if ((TargetPeak_Flag==0) && (ui>=Ubp)) {
					Krel_j_1 =  (Fbp) / (Ubp - u0);
				}
				else if ((TargetPeak_Flag==1) && (ui>=ULastPeak_neg_j_1)) {
					Krel_j_1 = (FLastPeak_neg_j_1) / (ULastPeak_neg_j_1 - u0);
				}
				else if ((TargetPeak_Flag==1) && (ui<=ULastPeak_neg_j_1)) {
					if ((FLastPeak_neg_j_1 >= Fbp) && (ui <= Ubp)) {
						Krel_j_1 = (Fbp-FLastPeak_neg_j_1)/(Ubp-ULastPeak_neg_j_1);
					} else if ((FLastPeak_neg_j_1 >= Fbp) && (ui <= Ubp)) {
						Krel_j_1 = (Fpeak_neg_j_1-Fbp)/(Upeak_neg_j_1-Ubp);
					} else {
						Krel_j_1 = (Fpeak_neg_j_1-FLastPeak_neg_j_1)/(Upeak_neg_j_1-ULastPeak_neg_j_1);
					}
				}
				df = du * Krel_j_1;
				//cout << "  Case = 4-" << endln;

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
		//cout << "  ENERGY: dEi=" << dEi << " Ei=" << Ei  << " Energy_Diss=" << Energy_Diss << " Energy_Acc=" << Energy_Acc << endln;
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
	cReversal_Flag		= Reversal_Flag;

	cKrel_j_1 = Krel_j_1;

	cUbp = Ubp;
	cFbp = Fbp;

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
	
	Upl = cUpl;	
	Ubp = cUbp;
	Fbp = cFbp;
	
	return 0;
}

int IMKPinching::revertToStart(void)
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
	Yield_Flag 		= cYield_Flag	   = 0;
	Reversal_Flag	= cReversal_Flag   = 0;

	ULastPeak_pos_j_1 =  Uy_pos;
	FLastPeak_pos_j_1 =  Fy_pos;
	ULastPeak_neg_j_1 = -Uy_neg;
	FLastPeak_neg_j_1 = -Fy_neg;
	
	Krel_j_1 	  = Ke;

	Upl = 0.0;
	Ubp = 0.0;
	Fbp = 0.0;
	
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
	
	cUpl = 0.0;
	cUbp = 0.0;
	cFbp = 0.0;
	
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
IMKPinching::getCopy(void)
{
	IMKPinching *theCopy = new IMKPinching(this->getTag(), Ke,
		Uy_pos, Umax_pos, Uu_pos, Fy_pos, FmaxFy_pos, ResF_pos,
		Uy_neg, Umax_neg, Uu_neg, Fy_neg, FmaxFy_neg, ResF_neg,
		LAMBDA_S, LAMBDA_C, LAMBDA_A, LAMBDA_K, c_S, c_C, c_A, c_K, D_pos, D_neg, kappaF, kappaD);

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
	theCopy->Yield_Flag		 = Yield_Flag;
	theCopy->Reversal_Flag   = Reversal_Flag;

	theCopy->Krel_j_1 = Krel_j_1;
	
	theCopy->Upl = Upl;
	theCopy->Ubp = Ubp;
	theCopy->Fbp = Fbp;

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
	
	theCopy->cUpl = cUpl;
	theCopy->cUbp = cUbp;
	theCopy->cFbp = cFbp;
	
	return theCopy;
}

int IMKPinching::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	cout << " sendSelf" << endln;

	static Vector data(144);
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
	data(24) = kappaF;
	data(25) = kappaD;
	data(26)  = ui;
	data(27)  = fi;
	data(28)  = ui_1;
	data(29)  = fi_1;
	data(30)  = du_i_1;
	data(31)  = Uy_pos_j_1;
	data(32)  = Umax_pos_j_1;
	data(33)  = Fy_pos_j_1;
	data(34)  = Fmax_pos_j_1;
	data(35)  = Upeak_pos_j_1;
	data(36)  = Fpeak_pos_j_1;
	data(37)  = Ures_pos_j_1;
	data(38)  = Fres_pos_j_1;
	data(39)  = Kp_pos_j_1;
	data(40)  = Kpc_pos_j_1;
	data(41)  = Uy_neg_j_1;
	data(42)  = Umax_neg_j_1;
	data(43)  = Fy_neg_j_1;
	data(44)  = Fmax_neg_j_1;
	data(45)  = Upeak_neg_j_1;
	data(46)  = Fpeak_neg_j_1;
	data(47)  = Ures_neg_j_1;
	data(48)  = Fres_neg_j_1;
	data(49)  = Kp_neg_j_1;
	data(50)  = Kpc_neg_j_1;
	data(51)  = Kul_j_1;
	data(52)  = Failure_Flag;
	data(53)  = Excursion_Flag;
	data(54)  = Unloading_Flag;
	data(55)  = Reloading_Flag;
	data(56)  = TargetPeak_Flag;
	data(57)  = Yield_Flag;
	data(58)  = Energy_Acc;
	data(59)  = Energy_Diss;
	data(60)  = u0;
	data(61)  = du;
	data(62)  = df;
	data(63)  = FailS;
	data(64)  = FailC;
	data(65)  = FailA;
	data(66)  = FailK;
	data(67)  = Ei;
	data(68)  = dEi;
	data(69)  = Epj;
	data(70)  = EpjK;
	data(71)  = EiK;
	data(72)  = c_S;
	data(73)  = c_C;
	data(74)  = c_A;
	data(75)  = c_K;
	data(76)  = EtS;
	data(77)  = EtC;
	data(78)  = EtA;
	data(79)  = EtK;
	data(80)  = betaS;
	data(81)  = betaC;
	data(82)  = betaA;
	data(83)  = betaK;
	data(84)  = sPCsp;
	data(85)  = sPCpcp;
	data(86)  = TangentK;
	data(87)  = Uy_pos;
	data(88)  = Umax_pos;
	data(89)  = Fmax_pos;
	data(90)  = Kp_pos;
	data(91)  = Kpc_pos;
	data(92)  = Uy_neg;
	data(93)  = Umax_neg;
	data(94)  = Fmax_neg;
	data(95)  = Kp_neg;
	data(96)  = Kpc_neg;
	data(97)  = cui;
	data(98)  = cfi;
	data(99)  = cui_1;
	data(100) = cfi_1;
	data(101) = cdu_i_1;
	data(102) = cUy_pos_j_1;
	data(103) = cUmax_pos_j_1;
	data(104) = cFy_pos_j_1;
	data(105) = cFmax_pos_j_1;
	data(106) = cUpeak_pos_j_1;
	data(107) = cFpeak_pos_j_1;
	data(108) = cUres_pos_j_1;
	data(109) = cFres_pos_j_1;
	data(110) = cKp_pos_j_1;
	data(111) = cKpc_pos_j_1;
	data(112) = cUy_neg_j_1;
	data(113) = cUmax_neg_j_1;
	data(114) = cFy_neg_j_1;
	data(115) = cFmax_neg_j_1;
	data(116) = cUpeak_neg_j_1;
	data(117) = cFpeak_neg_j_1;
	data(118) = cUres_neg_j_1;
	data(119) = cFres_neg_j_1;
	data(120) = cKp_neg_j_1;
	data(121) = cKpc_neg_j_1;
	data(122) = cKul_j_1;
	data(123) = cULastPeak_pos_j_1;
	data(124) = cFLastPeak_pos_j_1;
	data(125) = cULastPeak_neg_j_1;
	data(126) = cFLastPeak_neg_j_1;
	data(127) = cFailure_Flag;
	data(128) = cExcursion_Flag;
	data(129) = cReloading_Flag;
	data(130) = cUnloading_Flag;
	data(131) = cTargetPeak_Flag;
	data(132) = cYield_Flag;
	data(133) = cKrel_j_1;
	data(134) = Krel_LastPeak;
	data(135) = Krel_GlobalPeak;
	data(136) = K_check;
	data(137) = Upl;
	data(138) = Ubp;
	data(139) = Fbp;
	data(140) = cUbp;
	data(141) = cFbp;
	data(142) = Reversal_Flag;
	data(143) = cReversal_Flag;

	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "IMKPinching::sendSelf() - failed to send data\n";

	return res;
}

int IMKPinching::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector data(144);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);

	if (res < 0) {
		opserr << "IMKPinching::recvSelf() - failed to receive data\n";
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
		kappaF				= data(24);
		kappaD	 			= data(25);
		ui					= data(26);
		fi					= data(27);
		ui_1				= data(28);
		fi_1				= data(29);
		du_i_1				= data(30);
		Uy_pos_j_1			= data(31);
		Umax_pos_j_1		= data(32);
		Fy_pos_j_1			= data(33);
		Fmax_pos_j_1		= data(34);
		Upeak_pos_j_1		= data(35);
		Fpeak_pos_j_1		= data(36);
		Ures_pos_j_1		= data(37);
		Fres_pos_j_1		= data(38);
		Kp_pos_j_1			= data(39);
		Kpc_pos_j_1			= data(40);
		Uy_neg_j_1			= data(41);
		Umax_neg_j_1		= data(42);
		Fy_neg_j_1			= data(43);
		Fmax_neg_j_1		= data(44);
		Upeak_neg_j_1		= data(45);
		Fpeak_neg_j_1		= data(46);
		Ures_neg_j_1		= data(47);
		Fres_neg_j_1		= data(48);
		Kp_neg_j_1			= data(49);
		Kpc_neg_j_1			= data(50);
		Failure_Flag		= data(51);
		Excursion_Flag		= data(52);
		Reloading_Flag		= data(53);
		Unloading_Flag		= data(54);
		TargetPeak_Flag		= data(55);
		Yield_Flag			= data(56);
		Kul_j_1				= data(57);
		Energy_Acc			= data(58);
		Energy_Diss			= data(59);
		u0					= data(60);
		du					= data(61);
		df					= data(62);
		FailS				= data(63);
		FailC				= data(64);
		FailA				= data(65);
		FailK				= data(66);
		Ei					= data(67);
		dEi					= data(68);
		Epj					= data(79);
		EpjK				= data(70);
		EiK					= data(71);
		c_S					= data(72);
		c_C					= data(73);
		c_A					= data(74);
		c_K					= data(75);
		EtS					= data(76);
		EtC					= data(77);
		EtA					= data(78);
		EtK					= data(79);
		betaS				= data(80);
		betaC				= data(81);
		betaA				= data(82);
		betaK				= data(83);
		sPCsp				= data(84);
		sPCpcp				= data(85);
		TangentK			= data(86);
		Uy_pos				= data(87);
		Umax_pos			= data(88);
		Fmax_pos			= data(89);
		Kp_pos				= data(90);
		Kpc_pos				= data(91);
		Uy_neg				= data(92);
		Umax_neg			= data(93);
		Fmax_neg			= data(94);
		Kp_neg				= data(95);
		Kpc_neg				= data(96);
		cui					= data(97);
		cfi					= data(98);
		cui_1				= data(99);
		cfi_1				= data(100);
		cdu_i_1				= data(101);
		cUy_pos_j_1			= data(102);
		cUmax_pos_j_1		= data(103);
		cFy_pos_j_1			= data(104);
		cFmax_pos_j_1		= data(105);
		cUpeak_pos_j_1		= data(106);
		cFpeak_pos_j_1		= data(107);
		cUres_pos_j_1		= data(108);
		cFres_pos_j_1		= data(109);
		cKp_pos_j_1			= data(110);
		cKpc_pos_j_1		= data(111);
		cUy_neg_j_1			= data(112);
		cUmax_neg_j_1		= data(113);
		cFy_neg_j_1			= data(114);
		cFmax_neg_j_1		= data(115);
		cUpeak_neg_j_1		= data(116);
		cFpeak_neg_j_1		= data(117);
		cUres_neg_j_1		= data(118);
		cFres_neg_j_1		= data(119);
		cKp_neg_j_1			= data(120);
		cKpc_neg_j_1		= data(121);
		cKul_j_1			= data(122);
		cULastPeak_pos_j_1	= data(123);
		cFLastPeak_pos_j_1	= data(124);
		cULastPeak_neg_j_1	= data(125);
		cFLastPeak_neg_j_1	= data(126);
		cFailure_Flag		= data(127);
		cExcursion_Flag		= data(128);
		cReloading_Flag		= data(129);
		cUnloading_Flag		= data(130);
		cTargetPeak_Flag	= data(131);
		cYield_Flag   		= data(132);
		cKrel_j_1      		= data(133);
		Krel_LastPeak		= data(134);
		Krel_GlobalPeak		= data(135);
		K_check				= data(136);
		Upl  				= data(137);
		Ubp  				= data(138);
		Fbp  				= data(139);
		cUbp 				= data(140);
		cFbp				= data(141);		
		cReversal_Flag		= data(142);
		Reversal_Flag		= data(143);
	}

	return res;
}

void IMKPinching::Print(OPS_Stream &s, int flag)
{
	cout << "IMKPinching tag: " << this->getTag() << endln;
}