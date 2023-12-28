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
        opserr << "IMK with Peak-Oriented Response - Code by AE_KI (Sep23)\n";
    }

    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial = 0;

    int    iData[1];
    double dData[23];
    int numInt = 1;

    if (OPS_GetIntInput(&numInt, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial IMKPeakOriented tag" << endln;
        return 0;
    }

    int numDouble = 23;


    if (OPS_GetDoubleInput(&numDouble, dData) != 0) {
        opserr << "Invalid Args want: uniaxialMaterial IMKPeakOriented tag? Ke? ";
        opserr << "posUp_0? posUpc_0? posUu_0? posFy_0? posFcapFy_0? posFresFy_0? ";
        opserr << "negUp_0? negUpc_0? negUu_0? negFy_0? negFcapFy_0? negFresFy_0? ";
        opserr << "LamdaS? LamdaC? LamdaA? LamdaK? Cs? Cc? Ca? Ck? D_pos? D_neg? ";
        return 0;
    }


    // Parsing was successful, allocate the material
    theMaterial = new IMKPeakOriented(iData[0], dData[0],
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
    // Make sure these are all positive 
    if (negUp_0 < 0)
       negUp_0 = -negUp_0;
    if (negUpc_0 < 0)
       negUpc_0 = -negUpc_0;
    if (negUu_0 < 0)
       negUu_0 = -negUu_0;
    if (negFy_0 < 0)
       negFy_0 = -negFy_0;
    
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
    const double Ui_1 = Ui;
    const double Fi_1 = Fi;
    Ui = strain; //set trial displacement
    const double dU = Ui - Ui_1;    // Incremental deformation at current step
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////  MAIN CODE //////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (Failure_Flag) {     // When a failure has already occured
        Fi = 0;
    } else if (dU == 0) {   // When deformation doesn't change from the last
        Fi = Fi_1;
    } else {
        double betaS = 0, betaC = 0, betaK = 0, betaA = 0;
        bool FailS = false, FailC = false, FailK = false, FailA = false;
        const bool onBackbone = (Branch > 1);
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////// WHEN REVERSAL /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        if ( (onBackbone && Fi_1 * dU < 0) || (onBackbone && Fi_1 == 0 && Ui_1 * dU <= 0) ) {
            Branch = 1;
    /////////////////////////// UPDATE PEAK POINTS ////////////////////////////////////////////
            if ( Fi_1 > 0 ){
                posUlocal = Ui_1;           // UPDATE LOCAL
                posFlocal = Fi_1;
                if ( Ui_1 > posUglobal ) {    // UPDATE GLOBAL
                    posUglobal = Ui_1;
                    posFglobal = Fi_1;
                }
            } else {
                negUlocal = Ui_1;           // UPDATE LOCAL
                negFlocal = Fi_1;
                if ( Ui_1 < negUglobal ) {    // UPDATE GLOBAL
                    negUglobal = Ui_1;
                    negFglobal = Fi_1;
                }
            }
    /////////////////// UPDATE UNLOADING STIFFNESS ////////////////////////////////////////////
            const double  EpjK = engAcml - 0.5 * (Fi_1 / Kunload) * Fi_1;
            const double  EiK = engAcml - engDspt - 0.5 * (Fi_1 / Kunload) * Fi_1;
            betaK = pow( (EiK / (engRefK - EpjK)), c_K );
            FailK = (betaK > 1);
            betaK = betaK < 0 ? 0 : (betaK > 1 ? 1 : betaK);
            Kunload *= (1 - betaK);
        }
        Fi = Fi_1 + Kunload * dU;
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////// WHEN NEW EXCURSION /////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        if (Branch == 1 && Fi_1 * Fi <= 0.0) {
    /////////////////// UPDATE BACKBONE CURVE /////////////////////////////////////////////////
            const double Ei = max(0.0, engAcml - engDspt);
            betaS = pow((Ei / (engRefS - engAcml)), c_S);
            betaC = pow((Ei / (engRefC - engAcml)), c_C);
            betaA = pow((Ei / (engRefA - engAcml)), c_A);
            FailS = (betaS > 1);
            FailC = (betaC > 1);
            FailA = (betaA > 1);
            betaS = betaS < 0 ? 0 : (betaS > 1 ? 1 : betaS);
            betaC = betaC < 0 ? 0 : (betaC > 1 ? 1 : betaC);
            betaA = betaA < 0 ? 0 : (betaA > 1 ? 1 : betaA);
            engDspt = engAcml;
        // Positive
            if (dU > 0) {
                double FcapProj = posFcap - posKpc * posUcap;
            // Yield Point
                posFy *= (1 - betaS * D_pos);
                posKp *= (1 - betaS * D_pos); // Post-Yield Stiffness
                FcapProj *= (1 - betaC * D_pos);
                posUglobal *= (1 + betaA * D_pos); // Accelerated Reloading Stiffness
                posUy = posFy / Ke;
            // Capping Point
                const double FyProj = posFy - posKp*posUy;
                posUcap = posKp <= posKpc ? 0 : (FcapProj - FyProj) / (posKp - posKpc);
                posFcap = FyProj + posKp*posUcap;
            // When a part of backbone is beneath the residual strength
            // Global Peak on the Updated Backbone
                if (posUglobal < posUy) {           // Elastic Branch
                    posFglobal = Ke * posUglobal;
                }
                else if (posUglobal < posUcap) {    // Post-Yield Branch
                    posFglobal = posFy + posKp * (posUglobal - posUy);
                }
                else {                              // Post-Capping Branch
                    posFglobal = posFcap + posKpc * (posUglobal - posUcap);
                }
                if (posFglobal < posFres) {     // Residual Branch
                    posFglobal = posFres;
                }
                posUres = (posFres - posFcap + posKpc * posUcap) / posKpc;
            }
        // Negative
            else {
                double FcapProj = negFcap - negKpc * negUcap;
            // Yield Point
                negFy *= (1 - betaS * D_neg);
                negKp *= (1 - betaS * D_neg); // Post-Yield Stiffness
                FcapProj *= (1 - betaC * D_neg);
                negUglobal *= (1 + betaA * D_neg); // Accelerated Reloading Stiffness
                negUy = negFy / Ke;
            // Capping Point
                const double FyProj = negFy - negKp * negUy;
                negUcap = negKp <= negKpc ? 0 : (FcapProj - FyProj) / (negKp - negKpc);
                negFcap = FyProj + negKp * negUcap;
            // When a part of backbone is beneath the residual strength
            // Global Peak on the Updated Backbone
                if (negUy < negUglobal) {           // Elastic Branch
                    negFglobal = Ke * negUglobal;
                }
                else if (negUcap < negUglobal) {    // Post-Yield Branch
                    negFglobal = negFy + negKp * (negUglobal - negUy);
                }
                else {                              // Post-Capping Branch
                    negFglobal = negFcap + negKpc * (negUglobal - negUcap);
                }
                if (negFres < negFglobal) {     // Residual Branch
                    negFglobal = negFres;
                }
                negUres = (negFres - negFcap + negKpc * negUcap) / negKpc;
            }
    ////////////////////////// RELOADING TARGET DETERMINATION /////////////////////////////////
            if (dU > 0) {
                const double u0 = Ui_1 - (Fi_1 / Kunload);
                const double Kglobal = posFglobal / (posUglobal - u0);
                const double Klocal = posFlocal / (posUlocal - u0);
                if ( u0 < posUlocal && posFlocal < posFglobal && Klocal > Kglobal) {
                    Branch = 3;
                    Kreload = Klocal;
                } else {
                    Branch = 4;
                    Kreload = Kglobal;
                }
            }
            else {
                const double u0 = Ui_1 - (Fi_1 / Kunload);
                const double Kglobal = negFglobal / (negUglobal - u0);
                const double Klocal = negFlocal / (negUlocal - u0);
                if ( u0 > negUlocal && negFlocal > negFglobal && Klocal > Kglobal) {
                    Branch = 13;
                    Kreload = Klocal;
                } else {
                    Branch = 14;
                    Kreload = Kglobal;
                }
            }
        }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ////////////////// BRANCH SHIFT CHECK AND TANGENT STIFFNESS UPDATE ////////////////////////
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
        // int exBranch = Branch;
        if (Branch == 0 && Ui > posUy) {            // Yield in Positive
            Branch = 5;
        } else if (Branch == 0 && Ui < negUy) {     // Yield in Negative
            Branch = 15;
        } else if (Branch == 1 && Fi_1 > 0 && Ui > posUlocal) {
            const double Kglobal = (posFglobal - posFlocal) / (posUglobal - posUlocal);
            Kreload = Kglobal;
            Branch = 4;                        // Towards Global Peak
        } else if (Branch == 1 && Fi_1 < 0 && Ui < negUlocal) {    // Back to Reloading (Negative)
            const double Kglobal = (negFglobal - negFlocal) / (negUglobal - negUlocal);
            Kreload = Kglobal;
            Branch = 14;                   // Towards Global Peak
        }
    // Positive
        if (Branch == 3 && Ui > posUlocal) {
            Kreload = (posFglobal - posFlocal) / (posUglobal - posUlocal);
            Branch = 4;
        }
        if (Branch == 4 && Ui > posUglobal) {
            Branch = 5;
        }
        if (Branch == 5 && Ui > posUcap) {
            Branch = 6;
        }
        if (Branch == 6 && Ui > posUres) {
            Branch = 7;
        }
    // Negative
        if (Branch == 13 && Ui < negUlocal) {
            Kreload = (negFglobal - negFlocal) / (negUglobal - negUlocal);
            Branch = 14;
        }
        if (Branch == 14 && Ui < negUglobal) {
            Branch = 15;
        }
        if (Branch == 15 && Ui < negUcap) {
            Branch = 16;
        }
        if (Branch == 16 && Ui < negUres) {
            Branch = 17;
        }
    // Branch Change check
        // if (Branch!=exBranch) {
        //     std::cout << exBranch << " -> " << Branch << "\n";
        // }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////// COMPUTE FORCE BASED ON BRANCH /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        if (Branch == 0) {
            Fi = Ke * Ui;
        } else if (Branch == 1) {
            Fi = Fi_1 + Kunload * dU;
    // Positive
        } else if (Branch == 3) {
            Fi = posFlocal + Kreload * (Ui - posUlocal);
        } else if (Branch == 4) {
            Fi = posFglobal + Kreload * (Ui - posUglobal);
        } else if (Branch == 5) {
            Fi = posFcap + posKp * (Ui - posUcap);
        } else if (Branch == 6) {
            Fi = posFcap + posKpc * (Ui - posUcap);
        } else if (Branch == 7) {
            Fi = posFres;
    // Negative
        } else if (Branch == 13) {
            Fi = negFlocal + Kreload * (Ui - negUlocal);
        } else if (Branch == 14) {
            Fi = negFglobal + Kreload * (Ui - negUglobal);
        } else if (Branch == 15) {
            Fi = negFcap + negKp * (Ui - negUcap);
        } else if (Branch == 16) {
            Fi = negFcap + negKpc * (Ui - negUcap);
        } else if (Branch == 17) {
            Fi = negFres;
        }
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    // CHECK FOR FAILURE
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
        const bool FailPp = ( posFglobal == 0 );
        const bool FailPn = ( negFglobal == 0 );
        const bool FailDp = ( dU > 0 && Ui >=  posUu_0 );
        const bool FailDn = ( dU < 0 && Ui <= -negUu_0 );
        const bool FailRp = ( Branch == 7 && Fi <= 0);
        const bool FailRn = ( Branch == 17 && Fi >= 0);
        if (FailS || FailC || FailA || FailK || FailPp || FailPn || FailRp || FailRn || FailDp || FailDn) {
            Fi = 0;
            Failure_Flag = true;
        }

        engAcml += 0.5 * (Fi + Fi_1) * dU;   // Internal energy increment

        KgetTangent = (Fi - Fi_1) / dU;
    }
    if (KgetTangent == 0) {
        KgetTangent = 1e-6;
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
    return (KgetTangent);
}

double IMKPeakOriented::getInitialTangent(void)
{
    //cout << " getInitialTangent" << endln;
    return (Ke);
}

double IMKPeakOriented::getStrain(void)
{
    //cout << " getStrain" << endln;
    return (Ui);
}

int IMKPeakOriented::commitState(void)
{
    //cout << " commitState" << endln;

    //commit trial  variables
// 12 Pos U and F
    cPosUy = posUy;
    cPosFy = posFy;
    cPosUcap = posUcap;
    cPosFcap = posFcap;
    cPosUlocal = posUlocal;
    cPosFlocal = posFlocal;
    cPosUglobal = posUglobal;
    cPosFglobal = posFglobal;
    cPosUres = posUres;
    cPosFres = posFres;
    cPosKp = posKp;
    cPosKpc = posKpc;
// 12 Neg U and F
    cNegUy = negUy;
    cNegFy = negFy;
    cNegUcap = negUcap;
    cNegFcap = negFcap;
    cNegUlocal = negUlocal;
    cNegFlocal = negFlocal;
    cNegUglobal = negUglobal;
    cNegFglobal = negFglobal;
    cNegUres = negUres;
    cNegFres = negFres;
    cNegKp = negKp;
    cNegKpc = negKpc;
// 3 State
    cUi = Ui;
    cFi = Fi;
// 2 Stiffness
    cKreload = Kreload;
    cKunload = Kunload;
// 2 Energy
    cEngAcml = engAcml;
    cEngDspt = engDspt;
// 2 Flag
    cFailure_Flag = Failure_Flag;
    cBranch = Branch;
    return 0;
}

int IMKPeakOriented::revertToLastCommit(void)
{
    //cout << " revertToLastCommit" << endln;
    //the opposite of commit trial history variables
// 12 Positive U and F
    posUy = cPosUy;
    posFy = cPosFy;
    posUcap = cPosUcap;
    posFcap = cPosFcap;
    posUlocal = cPosUlocal;
    posFlocal = cPosFlocal;
    posUglobal = cPosUglobal;
    posFglobal = cPosFglobal;
    posUres = cPosUres;
    posFres = cPosFres;
    posKp = cPosKp;
    posKpc = cPosKpc;
// 12 Negative U and F
    negUy = cNegUy;
    negFy = cNegFy;
    negUcap = cNegUcap;
    negFcap = cNegFcap;
    negUlocal = cNegUlocal;
    negFlocal = cNegFlocal;
    negUglobal = cNegUglobal;
    negFglobal = cNegFglobal;
    negUres = cNegUres;
    negFres = cNegFres;
    negKp = cNegKp;
    negKpc = cNegKpc;
// 3 State Variables
    Ui = cUi;
    Fi = cFi;
// 2 Stiffness
    Kreload = cKreload;
    Kunload = cKunload;
// 2 Energy
    engAcml = cEngAcml;
    engDspt = cEngDspt;
// 2 Flag
    Failure_Flag = cFailure_Flag;
    Branch = cBranch;
    return 0;
}

int IMKPeakOriented::revertToStart(void)
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\\
    //////////////////////////////////////////////////////////////////// ONE TIME CALCULATIONS ////////////////////////////////////////////////////////////////////\\
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
// 14 Initial Values
    posUy_0 = posFy_0 / Ke;
    posUcap_0 = posUy_0 + posUp_0;
    posFcap_0 = posFcapFy_0*posFy_0;
    posKp_0 = (posFcap_0 - posFy_0) / posUp_0;
    posKpc_0 = posFcap_0 / posUpc_0;
    negUy_0 = negFy_0 / Ke;
    negUcap_0 = negUy_0 + negUp_0;
    negFcap_0 = negFcapFy_0*negFy_0;
    negKp_0 = (negFcap_0 - negFy_0) / negUp_0;
    negKpc_0 = negFcap_0 / negUpc_0;
    engRefS = LAMBDA_S * posFy_0;
    engRefC = LAMBDA_C * posFy_0;
    engRefA = LAMBDA_A * posFy_0;
    engRefK = LAMBDA_K * posFy_0;
// 12 Positive U and F
    posUy = cPosUy = posUy_0;
    posFy = cPosFy = posFy_0;
    posUcap = cPosUcap = posUcap_0;
    posFcap = cPosFcap = posFcap_0;
    posUlocal = cPosUlocal = posUy_0;
    posFlocal = cPosFlocal = posFy_0;
    posUglobal = cPosUglobal = posUy_0;
    posFglobal = cPosFglobal = posFy_0;
    posFres = cPosFres = posFy_0*posFresFy_0;
    posKp = cPosKp =  posKp_0;
    posKpc = cPosKpc = -posKpc_0;
    posUres = cPosUres = (posFres - posFcap) / posKpc + posUcap;
// 12 Negative U and F
    negUy = cNegUy = -negUy_0;
    negFy = cNegFy = -negFy_0;
    negUcap = cNegUcap = -negUcap_0;
    negFcap = cNegFcap = -negFcap_0;
    negUlocal = cNegUlocal = -negUy_0;
    negFlocal = cNegFlocal = -negFy_0;
    negUglobal = cNegUglobal = -negUy_0;
    negFglobal = cNegFglobal = -negFy_0;
    negFres = cNegFres = -negFy_0*negFresFy_0;
    negKp = cNegKp =  negKp_0;
    negKpc = cNegKpc = -negKpc_0;
    negUres = cNegUres = (negFres - negFcap) / negKpc + negUcap;
// 3 State Values
    Ui = cUi = 0;
    Fi = cFi = 0;
// 2 Stiffness
    Kreload = cKreload = Ke;
    Kunload = cKunload = Ke;
    KgetTangent = Ke;
// 2 Energy
    engAcml = cEngAcml = 0.0;
    engDspt = cEngDspt = 0.0;
// 2 Flag
    Failure_Flag = cFailure_Flag = false;
    Branch = cBranch = 0;
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
// 12 Positive U and F
    theCopy->posUy = posUy;
    theCopy->posFy = posFy;
    theCopy->posUcap = posUcap;
    theCopy->posFcap = posFcap;
    theCopy->posUlocal = posUlocal;
    theCopy->posFlocal = posFlocal;
    theCopy->posUglobal = posUglobal;
    theCopy->posFglobal = posFglobal;
    theCopy->posUres = posUres;
    theCopy->posFres = posFres;
    theCopy->posKp = posKp;
    theCopy->posKpc = posKpc;
// 12 Negative U and F
    theCopy->negUy = negUy;
    theCopy->negFy = negFy;
    theCopy->negUcap = negUcap;
    theCopy->negFcap = negFcap;
    theCopy->negUlocal = negUlocal;
    theCopy->negFlocal = negFlocal;
    theCopy->negUglobal = negUglobal;
    theCopy->negFglobal = negFglobal;
    theCopy->negUres = negUres;
    theCopy->negFres = negFres;
    theCopy->negKp = negKp;
    theCopy->negKpc = negKpc;
// 3 State Values
    theCopy->Ui = Ui;
    theCopy->Fi = Fi;
// 2 Stiffness
    theCopy->Kreload = Kreload;
    theCopy->Kunload = Kunload;
// 2 Energy
    theCopy->engAcml = engAcml;
    theCopy->engDspt = engDspt;
// 2 Flag
    theCopy->Failure_Flag = Failure_Flag;
    theCopy->Branch = Branch;
// 12 Positive U and F
    theCopy->cPosUy = cPosUy;
    theCopy->cPosFy = cPosFy;
    theCopy->cPosUcap = cPosUcap;
    theCopy->cPosFcap = cPosFcap;
    theCopy->cPosUlocal = cPosUlocal;
    theCopy->cPosFlocal = cPosFlocal;
    theCopy->cPosUglobal= cPosUglobal;
    theCopy->cPosFglobal= cPosFglobal;
    theCopy->cPosUres = cPosUres;
    theCopy->cPosFres = cPosFres;
    theCopy->cPosKp = cPosKp;
    theCopy->cPosKpc = cPosKpc;
// 12 Negative U and F
    theCopy->cNegUy = cNegUy;
    theCopy->cNegFy = cNegFy;
    theCopy->cNegUcap = cNegUcap;
    theCopy->cNegFcap = cNegFcap;
    theCopy->cNegUglobal= cNegUglobal;
    theCopy->cNegFglobal= cNegFglobal;
    theCopy->cNegUlocal = cNegUlocal;
    theCopy->cNegFlocal = cNegFlocal;
    theCopy->cNegUres = cNegUres;
    theCopy->cNegFres = cNegFres;
    theCopy->cNegKp = cNegKp;
    theCopy->cNegKpc = cNegKpc;
// 3 State
    theCopy->cUi = cUi;
    theCopy->cFi = cFi;
// 2 Stiffness
    theCopy->cKreload = cKreload;
    theCopy->cKunload = cKunload;
// 2 Energy
    theCopy->cEngAcml = cEngAcml;
    theCopy->cEngDspt = cEngDspt;
// 2 Flag
    theCopy->cFailure_Flag = cFailure_Flag;
    theCopy->cBranch = cBranch;
    return theCopy;
}

int IMKPeakOriented::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    cout << " sendSelf" << endln;

    static Vector data(137);
    data(0) = this->getTag();
// 23 Fixed Input Material Parameters 1-25
    data(1) = Ke;
    data(2) = posUp_0;
    data(3) = posUpc_0;
    data(4) = posUu_0;
    data(5) = posFy_0;
    data(6) = posFcapFy_0;
    data(7) = posFresFy_0;
    data(8) = negUp_0;
    data(9) = negUpc_0;
    data(10) = negUu_0;
    data(11) = negFy_0;
    data(12) = negFcapFy_0;
    data(13) = negFresFy_0;
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
// 14 Initial Values 31-44
    data(31) = posUy_0;
    data(32) = posUcap_0;
    data(33) = posFcap_0;
    data(34) = posKp_0;
    data(35) = posKpc_0;
    data(36) = negUy_0;
    data(37) = negUcap_0;
    data(38) = negFcap_0;
    data(39) = negKp_0;
    data(40) = negKpc_0;
    data(41) = engRefS;
    data(42) = engRefC;
    data(43) = engRefA;
    data(44) = engRefK;
// 12 Positive U and F 51-62
    data(51) = posUy;
    data(52) = posFy;
    data(53) = posUcap;
    data(54) = posFcap;
    data(55) = posUlocal;
    data(56) = posFlocal;
    data(57) = posUglobal;
    data(58) = posFglobal;
    data(59) = posUres;
    data(60) = posFres;
    data(51) = posKp;
    data(62) = posKpc;
// 3 State Variables 63-65
    data(64) = Ui;
    data(65) = Fi;
// 2 Stiffness 66 67
    data(66) = Kreload;
    data(67) = Kunload;
// 2 Energy 68 69
    data(68) = engAcml;
    data(69) = engDspt;
// 12 Negative U and F 71-82
    data(71) = negUy;
    data(72) = negFy;
    data(73) = negUcap;
    data(74) = negFcap;
    data(75) = negUlocal;
    data(76) = negFlocal;
    data(77) = negUglobal;
    data(78) = negFglobal;
    data(79) = negUres;
    data(80) = negFres;
    data(81) = negKp;
    data(82) = negKpc;
// 2 Flag 85 86
    data(85) = Failure_Flag;
    data(86) = Branch;
// 12 Positive U and F 101-112
    data(101) = cPosUy;
    data(102) = cPosFy;
    data(103) = cPosUcap;
    data(104) = cPosFcap;
    data(105) = cPosUlocal;
    data(106) = cPosFlocal;
    data(107) = cPosUglobal;
    data(108) = cPosFglobal;
    data(109) = cPosUres;
    data(110) = cPosFres;
    data(111) = cPosKp;
    data(112) = cPosKpc;
// 3 State Variables 113-115
    data(114) = cUi;
    data(115) = cFi;
// 2 Stiffness 116 117
    data(116) = cKreload;
    data(117) = cKunload;
// 2 Energy 118 119
    data(118) = cEngAcml;
    data(119) = cEngDspt;
// 12 Negative U and F 121-132
    data(121) = cNegUy;
    data(122) = cNegFy;
    data(123) = cNegUcap;
    data(124) = cNegFcap;
    data(125) = cNegUlocal;
    data(126) = cNegFlocal;
    data(127) = cNegUglobal;
    data(128) = cNegFglobal;
    data(129) = cNegUres;
    data(130) = cNegFres;
    data(131) = cNegKp;
    data(132) = cNegKpc;
// 2 Flag 135 136
    data(135) = cFailure_Flag;
    data(136) = cBranch;
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
        Ke = data(1);
        posUp_0 = data(2);
        posUpc_0 = data(3);
        posUu_0 = data(4);
        posFy_0 = data(5);
        posFcapFy_0 = data(6);
        posFresFy_0 = data(7);
        negUp_0 = data(8);
        negUpc_0 = data(9);
        negUu_0 = data(10);
        negFy_0 = data(11);
        negFcapFy_0 = data(12);
        negFresFy_0 = data(13);
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
    // 14 Initial Values
        posUy_0 = data(31);
        posUcap_0 = data(32);
        posFcap_0 = data(33);
        posKp_0 = data(34);
        posKpc_0 = data(35);
        negUy_0 = data(36);
        negUcap_0 = data(37);
        negFcap_0 = data(38);
        negKp_0 = data(39);
        negKpc_0 = data(40);
        engRefS = data(41);
        engRefC = data(42);
        engRefA = data(43);
        engRefK = data(44);
    // 12 Positive U and F
        posUy = data(51);
        posFy = data(52);
        posUcap = data(53);
        posFcap = data(54);
        posUlocal = data(55);
        posFlocal = data(56);
        posUglobal = data(57);
        posFglobal = data(58);
        posUres = data(59);
        posFres = data(60);
        posKp = data(61);
        posKpc = data(62);
    // 3 State Variables
        Ui = data(64);
        Fi = data(65);
    // 2 Stiffness
        Kreload = data(66);
        Kunload = data(67);
    // 2 Energy
        engAcml = data(68);
        engDspt = data(69);
    // 12 Negative U and F
        negUy = data(71);
        negFy = data(72);
        negUcap = data(73);
        negFcap = data(74);
        negUlocal = data(75);
        negFlocal = data(76);
        negUglobal = data(77);
        negFglobal = data(78);
        negUres = data(79);
        negFres = data(80);
        negKp = data(81);
        negKpc = data(82);
    // 2 Flag
        Failure_Flag = data(85);
        Branch = data(86);
    // 12 Positive U and F
        cPosUy = data(101);
        cPosFy = data(102);
        cPosUcap = data(103);
        cPosFcap = data(104);
        cPosUlocal = data(105);
        cPosFlocal = data(106);
        cPosUglobal = data(107);
        cPosFglobal = data(108);
        cPosUres = data(109);
        cPosFres = data(110);
        cPosKp = data(111);
        cPosKpc = data(112);
    // 3 State Variables
        cUi = data(114);
        cFi = data(115);
    // 2 Stiffness
        cKreload = data(116);
        cKunload = data(117);
    // 2 Energy
        cEngAcml = data(118);
        cEngDspt = data(119);
    // 12 Negative U and F
        cNegUy = data(121);
        cNegFy = data(122);
        cNegUcap = data(123);
        cNegFcap = data(124);
        cNegUlocal = data(125);
        cNegFlocal = data(126);
        cNegUglobal = data(127);
        cNegFglobal = data(128);
        cNegUres = data(129);
        cNegFres = data(130);
        cNegKp = data(131);
        cNegKpc = data(132);
    // 2 Flag
        cFailure_Flag = data(135);
        cBranch = data(136);
    }

    return res;
}

void IMKPeakOriented::Print(OPS_Stream &s, int flag)
{
    cout << "IMKPeakOriented tag: " << this->getTag() << endln;
}
