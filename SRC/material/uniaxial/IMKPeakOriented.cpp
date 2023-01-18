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
		OPS_Error("IMK with Peak-Oriented Response - Code by AE-KI (Oct22)\n", 1);
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
        opserr << "posUp_0? posUpc_0? posUu_0? posFy_0? posFcapFy_0? posFresFy_0? ";
        opserr << "negUp_0? negUpc_0? negUu_0? negFy_0? negFcapFy_0? negFresFy_0? ";
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
    double p_posUp_0, double p_posUpc_0, double p_posUu_0, double p_posFy_0, double p_posFcapFy_0, double p_posFresFy_0,
    double p_negUp_0, double p_negUpc_0, double p_negUu_0, double p_negFy_0, double p_negFcapFy_0, double p_negFresFy_0,
    double p_LAMBDA_S, double p_LAMBDA_C, double p_LAMBDA_A, double p_LAMBDA_K, double p_c_S, double p_c_C, double p_c_A, double p_c_K, double p_D_pos, double p_D_neg)
    : UniaxialMaterial(tag, 0), Ke(p_Ke),
    posUp_0(p_posUp_0), posUpc_0(p_posUpc_0), posUu_0(p_posUu_0), posFy_0(p_posFy_0), posFcapFy_0(p_posFcapFy_0), posFresFy_0(p_posFresFy_0),
    negUp_0(p_negUp_0), negUpc_0(p_negUpc_0), negUu_0(p_negUu_0), negFy_0(p_negFy_0), negFcapFy_0(p_negFcapFy_0), negFresFy_0(p_negFresFy_0),
    LAMBDA_S(p_LAMBDA_S), LAMBDA_C(p_LAMBDA_C), LAMBDA_A(p_LAMBDA_A), LAMBDA_K(p_LAMBDA_K), c_S(p_c_S), c_C(p_c_C), c_A(p_c_A), c_K(p_c_K), D_pos(p_D_pos), D_neg(p_D_neg)
{
	this->revertToStart();
}

IMKPeakOriented::IMKPeakOriented()
    :UniaxialMaterial(0, 0), Ke(0),
    posUp_0(0), posUpc_0(0), posUu_0(0), posFy_0(0), posFcapFy_0(0), posFresFy_0(0),
    negUp_0(0), negUpc_0(0), negUu_0(0), negFy_0(0), negFcapFy_0(0), negFresFy_0(0),
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
    U		= strain; //set trial displacement
    double  Ui_1	= Ui;
    double  Fi_1	= Fi;
    Ui		= U;
    double  dU	    = Ui - Ui_1;    // Incremental deformation at current step
    double  dEi     = 0;
    ki      = Ktangent;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////  MAIN CODE //////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (Failure_Flag) {     // When a failure has already occured
        Fi 	= 0;
        dEi	= 0;
    } else if (dU == 0) {   // When deformation doesn't change from the last
        Fi 	= Fi_1;
        dEi	= 0;
    } else {
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// FLAG RAISE ////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    //      Excursion_Flag: When crossing X-axis.  Evokes re-considering of the deteriorations and which peak to go for.
    //      Reversal_Flag:  When unloading starts. Evokes re-condiersing of the stiffness deterioration and peak point registration.
        double  dF,betaS=0,betaC=0,betaA=0,betaK=0;
        int     exBranch       	= Branch;
        bool    Excursion_Flag 	= false;
        bool    Reversal_Flag  	= false;
        if (Branch == 0) {                          // Elastic Branch
            if (Ui > posUy) {                           // Yield in Positive
                posYield_Flag   = true;
                Branch 	= 5;
            } else if (Ui < negUy) {                    // Yield in Negative
                negYield_Flag   = true;
                Branch 	= 15;
            }
        } else if (Branch == 1) {                   // Unloading Branch
            if (Fi_1*(Fi_1+dU*Kunload) <= 0) {          // Crossing X-axis and Reloading Towards Opposite Direction
                Excursion_Flag 	= true;
            } else if (Fi_1 > 0 && Ui > posUlocal) {    // Back to Reloading (Positive)
                Branch  = 4;                            // Towards Global Peak
            } else if (Fi_1 < 0 && Ui < negUlocal) {    // Back to Reloading (Negative)
                Branch  = 14;                           // Towards Global Peak
            }
        } else if (Fi_1*dU < 0) {                   // Reversal from Reloading or Backbone Section
            Reversal_Flag  	= true;
            Branch 	= 1;
        }
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////// UPDATE PEAK POINTS ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        if (Reversal_Flag) {
            if ( Fi_1 > 0 ){
                posUlocal	= Ui_1;           // UPDATE LOCAL
                posFlocal	= Fi_1;
                if ( Ui_1 > posUglobal ) {    // UPDATE GLOBAL
                    posUglobal   	= Ui_1;
                    posFglobal   	= Fi_1;
                }
            } else {
                negUlocal	= Ui_1;           // UPDATE LOCAL
                negFlocal	= Fi_1;
                if ( Ui_1 < negUglobal ) {    // UPDATE GLOBAL
                    negUglobal   	= Ui_1;
                    negFglobal   	= Fi_1;
                }
            }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////// UPDATE UNLOADING STIFFNESS ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
            double  EpjK    = engAcml            - 0.5*(Fi_1 / Kunload)*Fi_1;
            double  EiK     = engAcml - engDspt  - 0.5*(Fi_1 / Kunload)*Fi_1;
            betaK           = pow( (EiK / (engRefK - EpjK)), c_K );
            Kunload        *= (1 - betaK);
            Ktangent        = Kunload;
        // Detect unloading completed in a step.
            if (Fi_1*(Fi_1+dU*Kunload) <= 0) {
                exBranch        = 1;
                Excursion_Flag  = true;
                Reversal_Flag   = false;
            }
        }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////// UPDATE BACKBONE CURVE /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        // UPDATE DETERIORATION PARAMETERS AT EACH NEW EXCURSION
        if ( Excursion_Flag ) {
            double  Ei      = fmax(0, engAcml - engDspt);
            betaS   = pow((Ei / (engRefS - engAcml)), c_S);
            betaC   = pow((Ei / (engRefC - engAcml)), c_C);
            betaA   = pow((Ei / (engRefA - engAcml)), c_A);
            engDspt = engAcml;
        // Positive
            if (dU > 0 && posYield_Flag) {
            // Yield Point
                posFy   *= (1 - betaS * D_pos);
                posKp   *= (1 - betaS * D_pos); // Post-Yield Stiffness
                if (posFy < posFres) {
                    posFy   = posFres;
                    posKp   = 0;
                }
                posUy       = posFy / Ke;
            // Capping Point
                posUcap     = ((1 - betaC * D_pos)*(posFcap-posKpc*posUcap) - (posFy-posKp*posUy))/(posKp-posKpc);
                posFcap     = posFy + posKp * (posUcap - posUy);
            // Accelerated reloading stiffness deterioration: Target peak deformation point
                posUglobal  *= (1 + betaA * D_pos);
            // Global Peak on the Updated Backbone
                if (posUglobal < posUy) {           // Elastic Branch
                    posFglobal = Ke * posUglobal;
                }
                else if (posUglobal < posUcap) {    // Post-Yield Branch
                    posFglobal = posFy   + posKp  * (posUglobal - posUy);
                }
                else {                              // Post-Capping Branch
                    posFglobal = posFcap + posKpc * (posUglobal - posUcap);
                    if (posFglobal < posFres) {     // Residual Branch
                        posFglobal  = posFres;
                    }
                }
                posUres = (posFres - posFcap + posKpc * posUcap) / posKpc;
            }
        // Negative
            else if (dU < 0 && negYield_Flag){
            // Yield Point
                negFy	*= (1 - betaS * D_neg);
                negKp	*= (1 - betaS * D_neg); // Post-Yield Stiffness
                if (negFy > negFres) {
                    negFy	= negFres;
                    negKp	= 0;
                }
                negUy		= negFy / Ke;
            // Capping Point
                negUcap     = ((1 - betaC * D_neg)*(negFcap-negKpc*negUcap) - (negFy-negKp*negUy))/(negKp-negKpc);
                negFcap     = negFy + negKp * (negUcap - negUy);
            // Accelerated reloading stiffness deterioration: Target peak deformation point
                negUglobal	*= (1 + betaA * D_neg);
            // Global Peak on the Updated Backbone
                if (negUy < negUglobal) {           // Elastic Branch
                    negFglobal	= Ke * negUglobal; 
                }
                else if (negUcap < negUglobal) {    // Post-Yield Branch
                    negFglobal	= negFy   + negKp  * (negUglobal - negUy);
                }
                else {                              // Post-Capping Branch
                    negFglobal  = negFcap + negKpc * (negUglobal - negUcap);
                    if (negFres < negFglobal) {     // Residual Branch
                        negFglobal  = negFres;
                    }
                }
                negUres  	= (negFres - negFcap + negKpc * negUcap) / negKpc;
            }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// RELOADING TARGET DETERMINATION /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
            double  u0 	= Ui_1 - (Fi_1 / Kunload);
            double  Klocal, Kglobal;
            if (dU > 0) {
                Kglobal   	= posFglobal    / (posUglobal - u0);
                Klocal    	= posFlocal     / (posUlocal  - u0);
                if ( u0 < posUlocal && posFlocal < posFglobal && Klocal > Kglobal) {
                    Branch      = 3;
                    Ktangent     = Klocal;
                }
                else {
                    Branch      = 4;
                    Ktangent     = Kglobal;
                }
            }
            else {
                Kglobal     = negFglobal    / (negUglobal - u0);
                Klocal    	= negFlocal     / (negUlocal  - u0);
                if ( u0 > negUlocal && negFlocal > negFglobal && Klocal > Kglobal) {
                    Branch      = 13;
                    Ktangent     = Klocal;
                }
                else {
                    Branch      = 14;
                    Ktangent     = Kglobal;
                }
            }
        }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ////////////////// BRANCH SHIFT CHECK /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    //  Branch
    //      0:  Elastic
    //      1:  Unloading Branch
    //      3:  Towards Local Peak      +
    //      4:  Towards Global Peak     +
    //      5:  Towards Capping Point   +
    //      6:  Towards Residual Point  +
    //      7:  Residual Branch         +
    //      13: Towards Local Peak      -
    //      14: Towards Global Peak     -
    //      15: Towards Capping Point   -
    //      16: Towards Residual Point  -
    //      17: Residual Branch         -
    // Branch shifting from 3 -> 4 -> 5 -> 6 -> 7
        if (Branch == 3 && Ui > posUlocal) {
            exBranch    = 3;
            Branch      = 4;
        }
        if (Branch == 4 && Ui > posUglobal) {
            posYield_Flag   = true;
            exBranch    = 4;
            Branch 	    = 5;
        }
        if (Branch == 5 && Ui > posUcap) {
            exBranch    = 5;
            Branch 	    = 6;
        }
        if (Branch == 6 && Ui > posUres) {
            exBranch    = 6;
            Branch 	    = 7;
        }

        if (Branch == 13 && Ui < negUlocal) {
            exBranch    = 13;
            Branch      = 14;
        }
        if (Branch == 14 && Ui < negUglobal) {
            negYield_Flag   = true;
            exBranch    = 14;
            Branch 	    = 15;
        }
        if (Branch == 15 && Ui < negUcap) {
            exBranch    = 15;
            Branch 	    = 16;
        }
        if (Branch == 16 && Ui < negUres) {
            exBranch    = 16;
            Branch 	    = 17;
        }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// COMPUTE FORCE INCREMENT /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Without Branch Change
        if (Branch == exBranch || Branch == 1) {
            dF	= dU*Ktangent;
        }
    // With Branch Change
    // Positive Force
        else if (exBranch == 1 && Excursion_Flag) {
            double  u0 	= Ui_1 - (Fi_1 / Kunload);
            dF 	= 0             - Fi_1 + Ktangent*  (Ui - u0);
        }
        // CASE 4: WHEN RELOADING BUT BETWEEN LAST CYCLE PEAK POINT AND GLOBAL PEAK POINT
        // CASE 5: WHEN LOADING IN GENERAL TOWARDS THE TARGET PEAK
        // CASE 6: WHEN LOADING IN GENERAL TOWARDS THE LAST CYCLE PEAK POINT BUT BEYOND IT
        else if (Branch == 4) {
            Ktangent   	= (posFglobal   - posFlocal)    / (posUglobal - posUlocal);
            dF      	= posFlocal     - Fi_1          + Ktangent*(Ui - posUlocal);
        }
        // CASE 7: WHEN LOADING BEYOND THE TARGET PEAK BUT BEFORE THE CAPPING POINT
        else if (Branch == 5 && exBranch == 0) {
            dF 	= posFy       - Fi_1 + posKp*   (Ui - posUy);
            Ktangent	= posKp;
        }
        else if (Branch == 5) {
            dF 	= posFglobal  - Fi_1 + posKp*   (Ui - posUglobal);
            Ktangent	= posKp;
        }
        // CASE 8: WHEN LOADING AND BETWEEN THE CAPPING POINT AND THE RESIDUAL POINT
        else if (Branch == 6 && exBranch == 5) {
            dF 	= posFcap     - Fi_1 + posKpc*  (Ui - posUcap);
            Ktangent	= posKpc;
        }
        else if (Branch == 6) {
            dF 	= posFglobal  - Fi_1 + posKpc*  (Ui - posUglobal);
            Ktangent	= posKpc;
        }
        // CASE 9: WHEN LOADING AND BEYOND THE RESIDUAL POINT
        else if (Branch == 7) {
            dF 	= posFres     - Fi_1;
            Ktangent	= 0;
    // Negative Force
        }
        else if (Branch == 14) {
        // CASE 4: WHEN RELOADING BUT BETWEEN LAST CYCLE PEAK POINT AND GLOBAL PEAK POINT
        // CASE 5: WHEN LOADING IN GENERAL TOWARDS THE TARGET PEAK
        // CASE 6: WHEN LOADING IN GENERAL TOWARDS THE LAST CYCLE PEAK POINT BUT BEYOND IT
            Ktangent   	= (negFglobal - negFlocal) / (negUglobal - negUlocal);
            dF         	= negFlocal - Fi_1 + Ktangent*(Ui - negUlocal);
        }
        // CASE 7: WHEN LOADING BEYOND THE TARGET PEAK BUT BEFORE THE CAPPING POINT
        else if (Branch == 15 && exBranch == 0) {
            dF 	= negFy - Fi_1 + negKp*(Ui - negUy);
            Ktangent	= negKp;
        }
        else if (Branch == 15) {
            dF 	= negFglobal - Fi_1 + negKp*(Ui - negUglobal);
            Ktangent	= negKp;
        }
        // CASE 8: WHEN LOADING AND BETWEEN THE CAPPING POINT AND THE RESIDUAL POINT
        else if (Branch == 16 && exBranch == 15) {
            dF 	= negFcap - Fi_1 + negKpc*(Ui - negUcap);
            Ktangent	= negKpc;
        }
        else if (Branch == 16) {
            dF 	= negFglobal - Fi_1 + negKpc*(Ui - negUglobal);
            Ktangent	= negKpc;
        }
        // CASE 9: WHEN LOADING AND BEYOND THE RESIDUAL POINT
        else if (Branch == 17) {
            dF 	= negFres - Fi_1;
            Ktangent	= 0;
        }
    // Branch Change check
        // if (Branch!=exBranch) {
        //  std::cout << exBranch << " -> " << Branch << "\n";
        // }
// Force
        Fi	= Fi_1 + dF;
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // CHECK FOR FAILURE
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        // Failure criteria (Tolerance	= 1//)
    // I have no idea about why it can' t be 0 nor 1.
        bool    FailS	= ( betaS < -0.01 || betaS > 1.01 );
        bool    FailC	= ( betaC < -0.01 || betaC > 1.01 );
        bool	FailA	= ( betaA < -0.01 || betaA > 1.01 );
        bool	FailK	= ( betaK < -0.01 || betaK > 1.01 );
        bool	FailPp 	= ( posFglobal == 0               );
        bool	FailPn 	= ( negFglobal == 0               );
        bool	FailDp 	= ( Ui >  posUu_0                 );
        bool	FailDn 	= ( Ui < -negUu_0                 );
        bool	FailRp 	= ( Branch ==  7 && posFres == 0  );
        bool	FailRn 	= ( Branch == 17 && negFres == 0  );
        if (FailS||FailC||FailA||FailK||FailPp||FailPn||FailRp||FailRn||FailDp||FailDn) {
            Failure_Flag    = true;
        }
        if (Failure_Flag) {
            Fi 	= 0;
        }
        dEi	= 0.5*(Fi + Fi_1)*dU;   // Internal energy increment
        if (ki!=Ktangent) {
            ki  = (Fi - Fi_1) / dU;
        }
    }
    engAcml	+= dEi; 	            // Energy
    if (ki==0) {
        ki  = 1e-6;
    }
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
	return (Fi);
}

double IMKPeakOriented::getTangent(void)
{
    //cout << " getTangent" << endln;
    return (ki);
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
// 3 State
    cU		    = U;
    cUi	    	= Ui;
    cFi	        = Fi;
// 2 Stiffness
    cKtangent	= Ktangent;
    cKunload	= Kunload;
// 12 Pos U and F
    cPosUy	    = posUy;
    cPosFy	    = posFy;
    cPosUcap	= posUcap;
    cPosFcap	= posFcap;
    cPosUlocal	= posUlocal;
    cPosFlocal	= posFlocal;
    cPosUglobal	= posUglobal;
    cPosFglobal	= posFglobal;
    cPosUres	= posUres;
    cPosFres	= posFres;
    cPosKp	    = posKp;
    cPosKpc	    = posKpc;
// 12 Neg U and F
    cNegUy	    = negUy;
    cNegFy	    = negFy;
    cNegUcap	= negUcap;
    cNegFcap	= negFcap;
    cNegUlocal	= negUlocal;
    cNegFlocal	= negFlocal;
    cNegUglobal	= negUglobal;
    cNegFglobal	= negFglobal;
    cNegUres	= negUres;
    cNegFres	= negFres;
    cNegKp	    = negKp;
    cNegKpc	    = negKpc;
// 2 Energy
    cEngAcml	= engAcml;
    cEngDspt	= engDspt;
// 4 Flag
    cFailure_Flag		= Failure_Flag;
    cPosYield_Flag      = posYield_Flag;
    cNegYield_Flag      = negYield_Flag;
    cBranch				= Branch;
    return 0;
}

int IMKPeakOriented::revertToLastCommit(void)
{
    //cout << " revertToLastCommit" << endln;
    //the opposite of commit trial history variables
// 3 State Variables
    U	            = cU;
    Ui	            = cUi;
    Fi	            = cFi;
// 2 Stiffness
    Ktangent	    = cKtangent;
    Kunload	        = cKunload;
// 12 Positive U and F
    posUy	        = cPosUy;
    posFy	        = cPosFy;
    posUcap	        = cPosUcap;
    posFcap	        = cPosFcap;
    posUlocal	    = cPosUlocal;
    posFlocal	    = cPosFlocal;
    posUglobal	    = cPosUglobal;
    posFglobal	    = cPosFglobal;
    posUres	        = cPosUres;
    posFres	        = cPosFres;
    posKp	        = cPosKp;
    posKpc	        = cPosKpc;
// 12 Negative U and F
    negUy	        = cNegUy;
    negFy	        = cNegFy;
    negUcap	        = cNegUcap;
    negFcap	        = cNegFcap;
    negUlocal	    = cNegUlocal;
    negFlocal	    = cNegFlocal;
    negUglobal	    = cNegUglobal;
    negFglobal	    = cNegFglobal;
    negUres	        = cNegUres;
    negFres	        = cNegFres;
    negKp       	= cNegKp;
    negKpc	        = cNegKpc;
// 2 Energy
    engAcml	        = cEngAcml;
    engDspt	        = cEngDspt;
// 4 Flag
    Failure_Flag	= cFailure_Flag;
    posYield_Flag   = cPosYield_Flag;
    negYield_Flag   = cNegYield_Flag;
    Branch		    = cBranch;
    return 0;
}

int IMKPeakOriented::revertToStart(void)
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\\
    //////////////////////////////////////////////////////////////////// ONE TIME CALCULATIONS ////////////////////////////////////////////////////////////////////\\
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
// 14 Initial Values
    posUy_0  	= posFy_0   / Ke;
    posUcap_0	= posUy_0   + posUp_0;
    posFcap_0	= posFcapFy_0*posFy_0;
    posKp_0 	= (posFcap_0 - posFy_0) / posUp_0;
    posKpc_0 	= posFcap_0 / posUpc_0;
    negUy_0 	= negFy_0   / Ke;
    negUcap_0	= negUy_0   + negUp_0;
    negFcap_0	= negFcapFy_0*negFy_0;
    negKp_0 	= (negFcap_0 - negFy_0) / negUp_0;
    negKpc_0 	= negFcap_0 / negUpc_0;
    engRefS	    = LAMBDA_S  * posFy_0;
    engRefC	    = LAMBDA_C  * posFy_0;
    engRefA	    = LAMBDA_A  * posFy_0;
    engRefK	    = LAMBDA_K  * posFy_0;
// 12 Positive U and F
    posUy		= cPosUy	    = posUy_0;
    posFy    	= cPosFy 	    = posFy_0;
    posUcap  	= cPosUcap	    = posUcap_0;
    posFcap  	= cPosFcap	    = posFcap_0;
    posUlocal	= cPosUlocal	= posUy_0;
    posFlocal	= cPosFlocal	= posFy_0;
    posUglobal	= cPosUglobal	= posUy_0;
    posFglobal	= cPosFglobal	= posFy_0;
    posFres  	= cPosFres	    = posFy_0*posFresFy_0;

    posKp    	= cPosKp 	    =  posKp_0;
    posKpc   	= cPosKpc   	= -posKpc_0;
    posUres	    = cPosUres	    = (posFres - posFcap) / posKpc + posUcap;
// 12 Negative U and F
    negUy    	= cNegUy	    = -negUy_0;
    negFy    	= cNegFy	    = -negFy_0;
    negUcap  	= cNegUcap	    = -negUcap_0;
    negFcap  	= cNegFcap	    = -negFcap_0;
    negUlocal	= cNegUlocal	= -negUy_0;
    negFlocal	= cNegFlocal	= -negFy_0;
    negUglobal	= cNegUglobal	= -negUy_0;
    negFglobal	= cNegFglobal	= -negFy_0;
    negFres  	= cNegFres	    = -negFy_0*negFresFy_0;

    negKp    	= cNegKp	    =  negKp_0;
    negKpc   	= cNegKpc	    = -negKpc_0;
    negUres	    = cNegUres	    = (negFres - negFcap) / negKpc + negUcap;
// 3 State Values
    U	        = cU	        = 0;
    Ui      	= cUi 	        = 0;
    Fi 	        = cFi 	        = 0;
// 2 Stiffness
    Ktangent	= cKtangent	    = Ke;
    Kunload	    = cKunload	    = Ke;
// 4 Flag
    Failure_Flag 	= cFailure_Flag	  	= false;
    posYield_Flag   = cPosYield_Flag    = false;
    negYield_Flag   = cNegYield_Flag    = false;
    Branch      	= cBranch         	= 0;
// 2 Energy
    engAcml 	= cEngAcml	    = 0.0;
    engDspt	    = cEngDspt	    = 0.0;
    //cout << " revertToStart:" << endln; //<< " U=" << U << " Ui=" << Ui << " TanK=" << Ktangent << endln;
    return 0;
}

UniaxialMaterial *
IMKPeakOriented::getCopy(void)
{
    IMKPeakOriented *theCopy = new IMKPeakOriented(this->getTag(), Ke,
        posUy_0, posUcap_0, posUu_0, posFy_0, posFcapFy_0, posFresFy_0,
        negUy_0, negUcap_0, negUu_0, negFy_0, negFcapFy_0, negFresFy_0,
        LAMBDA_S, LAMBDA_C, LAMBDA_A, LAMBDA_K, c_S, c_C, c_A, c_K, D_pos, D_neg);

    //cout << " getCopy" << endln;
// 3 State Values
    theCopy->U          = U;
    theCopy->Ui         = Ui;
    theCopy->Fi         = Fi;
// 2 Stiffness
    theCopy->Ktangent  = Ktangent;
    theCopy->Kunload   = Kunload;
// 12 Positive U and F
    theCopy->posUy      = posUy;
    theCopy->posFy      = posFy;
    theCopy->posUcap    = posUcap;
    theCopy->posFcap    = posFcap;
    theCopy->posUlocal  = posUlocal;
    theCopy->posFlocal  = posFlocal;
    theCopy->posUglobal = posUglobal;
    theCopy->posFglobal = posFglobal;
    theCopy->posUres    = posUres;
    theCopy->posFres    = posFres;
    theCopy->posKp      = posKp;
    theCopy->posKpc     = posKpc;
// 12 Negative U and F
    theCopy->negUy      = negUy;
    theCopy->negFy      = negFy;
    theCopy->negUcap    = negUcap;
    theCopy->negFcap    = negFcap;
    theCopy->negUlocal  = negUlocal;
    theCopy->negFlocal  = negFlocal;
    theCopy->negUglobal = negUglobal;
    theCopy->negFglobal = negFglobal;
    theCopy->negUres    = negUres;
    theCopy->negFres    = negFres;
    theCopy->negKp      = negKp;
    theCopy->negKpc     = negKpc;
// 2 Energy
    theCopy->engAcml    = engAcml;
    theCopy->engDspt    = engDspt;
// 4 Flag
    theCopy->Failure_Flag   = Failure_Flag;
    theCopy->posYield_Flag  = posYield_Flag;
    theCopy->negYield_Flag  = negYield_Flag;
    theCopy->Branch     = Branch;
// 3 State
    theCopy->cU         = cU;
    theCopy->cUi        = cUi;
    theCopy->cFi        = cFi;
// 2 Stiffness
    theCopy->cKtangent = cKtangent;
    theCopy->cKunload  = cKunload;
// 12 Positive U and F
    theCopy->cPosUy     = cPosUy;
    theCopy->cPosFy     = cPosFy;
    theCopy->cPosUcap   = cPosUcap;
    theCopy->cPosFcap   = cPosFcap;
    theCopy->cPosUlocal = cPosUlocal;
    theCopy->cPosFlocal = cPosFlocal;
    theCopy->cPosUglobal= cPosUglobal;
    theCopy->cPosFglobal= cPosFglobal;
    theCopy->cPosUres   = cPosUres;
    theCopy->cPosFres   = cPosFres;
    theCopy->cPosKp     = cPosKp;
    theCopy->cPosKpc    = cPosKpc;
// 12 Negative U and F
    theCopy->cNegUy     = cNegUy;
    theCopy->cNegFy     = cNegFy;
    theCopy->cNegUcap   = cNegUcap;
    theCopy->cNegFcap   = cNegFcap;
    theCopy->cNegUglobal= cNegUglobal;
    theCopy->cNegFglobal= cNegFglobal;
    theCopy->cNegUlocal = cNegUlocal;
    theCopy->cNegFlocal = cNegFlocal;
    theCopy->cNegUres   = cNegUres;
    theCopy->cNegFres   = cNegFres;
    theCopy->cNegKp     = cNegKp;
    theCopy->cNegKpc    = cNegKpc;
// 2 Energy
    theCopy->cEngAcml   = cEngAcml;
    theCopy->cEngDspt   = cEngDspt;
// 4 Flag
    theCopy->cFailure_Flag  = cFailure_Flag;
    theCopy->cPosYield_Flag = cPosYield_Flag;
    theCopy->cNegYield_Flag = cNegYield_Flag;
    theCopy->cBranch    = Branch;
    return theCopy;
}

int IMKPeakOriented::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	cout << " sendSelf" << endln;

    static Vector data(137);
    data(0) = this->getTag();
// 23 Fixed Input Material Parameters 1-23
    data(1)  	= Ke;
    data(2)  	= posUp_0;
    data(3)  	= posUpc_0;
    data(4)  	= posUu_0;
    data(5)  	= posFy_0;
    data(6)  	= posFcapFy_0;
    data(7)  	= posFresFy_0;
    data(8)  	= negUp_0;
    data(9)  	= negUpc_0;
    data(10) 	= negUu_0;
    data(11) 	= negFy_0;
    data(12) 	= negFcapFy_0;
    data(13) 	= negFresFy_0;
    data(14) 	= LAMBDA_S;
    data(15) 	= LAMBDA_C;
    data(16) 	= LAMBDA_A;
    data(17) 	= LAMBDA_K;
    data(18) 	= c_S;
    data(19) 	= c_C;
    data(20) 	= c_A;
    data(21) 	= c_K;
    data(22) 	= D_pos;
    data(23) 	= D_neg;
// 14 Initial Values 25-38
    data(25)	= posUy_0;
    data(26)	= posUcap_0;
    data(27)	= posFcap_0;
    data(28)	= posKp_0;
    data(29)	= posKpc_0;
    data(30)	= negUy_0;
    data(31)	= negUcap_0;
    data(32)	= negFcap_0;
    data(33)	= negKp_0;
    data(34)	= negKpc_0;
    data(35)	= engRefS;
    data(36)	= engRefC;
    data(37)	= engRefA;
    data(38)	= engRefK;
// 3 State Variables 41-43
    data(41)    = U;
    data(42) 	= Ui;
    data(43) 	= Fi;
// 2 Stiffness 45-47
    data(45)	= Ktangent;
    data(46) 	= Kunload;
// 12 Positive U and F 51-62
    data(51) 	= posUy;
    data(52) 	= posFy;
    data(53) 	= posUcap;
    data(54) 	= posFcap;
    data(55) 	= posUlocal;
    data(56) 	= posFlocal;
    data(57) 	= posUglobal;
    data(58) 	= posFglobal;
    data(59) 	= posUres;
    data(60) 	= posFres;
    data(51) 	= posKp;
    data(62) 	= posKpc;
// 12 Negative U and F 65-76
    data(65) 	= negUy;
    data(66) 	= negFy;
    data(67) 	= negUcap;
    data(68) 	= negFcap;
    data(69) 	= negUlocal;
    data(70) 	= negFlocal;
    data(71) 	= negUglobal;
    data(72) 	= negFglobal;
    data(73) 	= negUres;
    data(74) 	= negFres;
    data(75) 	= negKp;
    data(76) 	= negKpc;
// 4 Flag 81-84
    data(81) 	= Failure_Flag;
    data(82)    = posYield_Flag;
    data(83)    = negYield_Flag;
    data(84) 	= Branch;
// 2 Energy 85-86
    data(85) 	= engAcml;
    data(86) 	= engDspt;
// 3 State Variables 91-93
    data(91)    = cU;
    data(92)	= cUi;
    data(93)	= cFi;
// 2 Stiffness 95-97
    data(95)    = cKtangent;
    data(96)	= cKunload;
// 12 Positive U and F 101-112
    data(101)	= cPosUy;
    data(102)	= cPosFy;
    data(103)	= cPosUcap;
    data(104)	= cPosFcap;
    data(105)	= cPosUlocal;
    data(106)	= cPosFlocal;
    data(107)	= cPosUglobal;
    data(108)	= cPosFglobal;
    data(109)	= cPosUres;
    data(110)	= cPosFres;
    data(111)	= cPosKp;
    data(112)	= cPosKpc;
// 12 Negative U and F 115-126
    data(115)	= cNegUy;
    data(116)	= cNegFy;
    data(117)	= cNegUcap;
    data(118)	= cNegFcap;
    data(119)	= cNegUlocal;
    data(120)	= cNegFlocal;
    data(121)	= cNegUglobal;
    data(122)	= cNegFglobal;
    data(123)	= cNegUres;
    data(124)	= cNegFres;
    data(125)	= cNegKp;
    data(126)	= cNegKpc;
// 4 Flag 131-134
    data(131)	= cFailure_Flag;
    data(132)   = cPosYield_Flag;
    data(133)   = cNegYield_Flag;
    data(134)	= cBranch;
// 2 Energy 135-136
    data(135)   = cEngAcml;
    data(136)   = cEngDspt;
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
    // 23 Fixed Input Material Parameters
        Ke				= data(1);
        posUp_0			= data(2);
        posUpc_0		= data(3);
        posUu_0			= data(4);
        posFy_0			= data(5);
        posFcapFy_0		= data(6);
        posFresFy_0		= data(7);
        negUp_0			= data(8);
        negUpc_0		= data(9);
        negUu_0			= data(10);
        negFy_0			= data(11);
        negFcapFy_0		= data(12);
        negFresFy_0		= data(13);
        LAMBDA_S		= data(14);
        LAMBDA_C		= data(15);
        LAMBDA_A		= data(16);
        LAMBDA_K		= data(17);
        c_S				= data(18);
        c_C				= data(19);
        c_A				= data(20);
        c_K				= data(21);
        D_pos			= data(22);
        D_neg			= data(23);
    // 14 Initial Values
        posUy_0			= data(25);
        posUcap_0		= data(26);
        posFcap_0		= data(27);
        posKp_0			= data(28);
        posKpc_0		= data(29);
        negUy_0			= data(30);
        negUcap_0		= data(31);
        negFcap_0		= data(32);
        negKp_0			= data(33);
        negKpc_0		= data(34);
        engRefS			= data(35);
        engRefC			= data(36);
        engRefA			= data(37);
        engRefK			= data(38);
    // 3 State Variables
        U               = data(41);
        Ui				= data(42);
        Fi				= data(43);
    // 2 Stiffness
        Ktangent		= data(45);
        Kunload	    	= data(46);
    // 12 Positive U and F
        posUy			= data(51);
        posFy			= data(52);
        posUcap		    = data(53);
        posFcap		   	= data(54);
        posUlocal	    = data(55);
        posFlocal	    = data(56);
        posUglobal	    = data(57);
        posFglobal	    = data(58);
        posUres	    	= data(59);
        posFres		   	= data(60);
        posKp			= data(61);
        posKpc		    = data(62);
    // 12 Negative U and F
        negUy			= data(65);
        negFy			= data(65);
        negUcap		    = data(67);
        negFcap		   	= data(68);
        negUlocal	    = data(69);
        negFlocal  	    = data(70);
        negUglobal	    = data(71);
        negFglobal  	= data(72);
        negUres		   	= data(73);
        negFres		    = data(74);
        negKp			= data(75);
        negKpc		    = data(76);
    // 4 Flag
        Failure_Flag	= data(81);
        posYield_Flag   = data(82);
        negYield_Flag   = data(83);
        Branch			= data(84);
    // 2 Energy
        engAcml			= data(85);
        engDspt		    = data(86);
    // 3 State Variables
        cU              = data(91);
        cUi				= data(92);
        cFi				= data(93);
    // 2 Stiffness
        cKtangent       = data(95);
        cKunload		= data(96);
    // 12 Positive U and F
        cPosUy			= data(101);
        cPosFy			= data(102);
        cPosUcap		= data(103);
        cPosFcap		= data(104);
        cPosUlocal		= data(105);
        cPosFlocal		= data(106);
        cPosUglobal		= data(107);
        cPosFglobal		= data(108);
        cPosUres		= data(109);
        cPosFres		= data(110);
        cPosKp			= data(111);
        cPosKpc			= data(112);
    // 12 Negative U and F
        cNegUy			= data(115);
        cNegFy			= data(116);
        cNegUcap		= data(117);
        cNegFcap		= data(118);
        cNegUlocal		= data(119);
        cNegFlocal		= data(120);
        cNegUglobal		= data(121);
        cNegFglobal		= data(122);
        cNegUres		= data(123);
        cNegFres		= data(124);
        cNegKp			= data(125);
        cNegKpc			= data(126);
    // 4 Flag
        cFailure_Flag	= data(131);
        cPosYield_Flag  = data(132);
        cNegYield_Flag  = data(133);
        cBranch			= data(134);
    // 2 Energy
        cEngAcml        = data(135);
        cEngDspt        = data(136);
    }

	return res;
}

void IMKPeakOriented::Print(OPS_Stream &s, int flag)
{
	cout << "IMKPeakOriented tag: " << this->getTag() << endln;
}