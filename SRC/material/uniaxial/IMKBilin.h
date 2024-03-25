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

//Modified Ibarra-Medina-Krawinkler with Bilin Hysteretic Response

//**********************************************************************                                                                     
// Code Developed by: Ahmed Elkady
// Lecturer,	University of Southampton
// Last Updated: October 24th 2022
//**********************************************************************

#ifndef IMKBilin_h
#define IMKBilin_h

#include <UniaxialMaterial.h>

class IMKBilin : public UniaxialMaterial
{
public:
    IMKBilin(int tag, double Ke,
        double posUp_0,  double posUpc_0, double posUu_0,  double posFy_0, double posFcapFy_0, double posFresFy_0,
        double negUp_0,  double negUpc_0, double negUu_0,  double negFy_0, double negFcapFy_0, double negFresFy_0,
        double LAMBDA_S, double LAMBDA_C, double LAMBDA_K, double c_S, double c_C, double c_K, double D_pos, double D_neg);
    IMKBilin();
    ~IMKBilin();
    const char *getClassType(void) const { return "IMKBilin"; };
    int setTrialStrain(double strain, double strainRate = 0.0);
    double  getStrain(void);
    double  getStress(void);
    double  getTangent(void);
    double  getInitialTangent(void);
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    UniaxialMaterial *getCopy(void);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);


protected:

private:
// 21 Fixed input material parameters
    double  Ke;
    double  posUp_0;
    double  posUpc_0;
    double  posUu_0;
    double  posFy_0;
    double  posFcapFy_0;
    double  posFresFy_0;
    double  negUp_0;
    double  negUpc_0;
    double  negUu_0;
    double  negFy_0;
    double  negFcapFy_0;
    double  negFresFy_0;
    double  LAMBDA_S;
    double  LAMBDA_C;
    double  LAMBDA_K;
    double  c_S;
    double  c_C;
    double  c_K;
    double  D_pos;
    double  D_neg;
// 13 Initial Variables
    double  posUy_0,        negUy_0;
    double  posUcap_0,      negUcap_0;
    double  posFcap_0,      negFcap_0;
    double  posKp_0,        negKp_0;
    double  posKpc_0,       negKpc_0;
    double  engRefS;
    double  engRefC;
    double  engRefK;
// History Variables 
// 7 Positive U and F
    double  posUy,          cPosUy;
    double  posFy,          cPosFy;
    double  posUcap,        cPosUcap;
    double  posFcap,        cPosFcap;
    double  posFres,        cPosFres;
    double  posKp,          cPosKp;
    double  posKpc,         cPosKpc;
// 7 Negative U and F
    double  negUy,          cNegUy;
    double  negFy,          cNegFy;
    double  negUcap,        cNegUcap;
    double  negFcap,        cNegFcap;
    double  negFres,        cNegFres;
    double  negKp,          cNegKp;
    double  negKpc,         cNegKpc;
// 3 State Variables 
    double  U,              cU;
    double  Ui,             cUi;
    double  Fi,             cFi;
// 2 Stiffness
    double  KgetTangent;
    double  Kunload,        cKunload;
// 2 Energy
    double  engAcml,        cEngAcml;
    double  engDspt,        cEngDspt;
// 2 Flag
    int     Failure_State,  cFailure_State;
    bool    onBackbone,     cOnBackbone;
};

#endif

//////////////////////////// Variables Definitions /////////////////////////////
/*
Ke 					Initial elastic stiffness
posUp_0 		Initial pre-capping plastic rotation in the +ve loading direction
posUpc_0   	Initial post-capping plastic rotation in the +ve loading direction
posUu_0 		Ultimate rotation in the +ve loading direction
posFy_0 			Initial effective plastic moment in the +ve loading direction
posFcapFy_0 		Initial maximum-to-effective plastic moment ratio in the +ve loading direction
posFresFy_0 			Residual moment in the +ve loading direction
negUp_0 		Initial pre-capping plastic rotation in the -ve loading direction
negUpc_0   	Initial post-capping plastic rotation in the -ve loading direction
negUu_0    	Ultimate rotation in the -ve loading direction
negFy_0        	Initial effective plastic moment in the -ve loading direction
negFcapFy_0    	Initial maximum-to-effective plastic moment ratio in the -ve loading direction
negFresFy_0       	Residual moment in the -ve loading direction
LAMBDA_S 			Cyclic deterioration parameter for strength deterioration 
LAMBDA_C 			Cyclic deterioration parameter for post-capping strength deterioration
LAMBDA_K 			Cyclic deterioration parameter for unloading stiffness deterioration 
c_S 				Rate of strength deterioration.
c_C 				Rate of post-capping strength deterioration.
c_K 				Rate of unloading stiffness deterioration.
D_pos 				Rate of cyclic deterioration in the +ve loading direction
D_neg 				Rate of cyclic deterioration in the -ve loading direction
n					Parameter identifying the offset rotation on the unloading side
Roffset				Offset rotation identifying the Smooth Transition region
LAMBDA_F			Cyclic deterioration parameter for Smooth Transition deterioration 
c_F					Cyclic deterioration parameter for Smooth Transition deterioration
Ui 					Rotation at current step
Fi 					Moment at current step 
Di 					Rotation Direction at current step 
Ui_1  				Rotation at previous step
Fi_1            	Moment at previous step 
Di_1            	Rotation Direction at previous step
Ulocal 			Rotation at direction reversal points
Flocal 			Moment at direction reversal points
Ktangent 			Tangent stiffness
Kunload 				Unloading stiffness in the previous excursion
posUy 	Yielding rotation in the previous +ve excursion
posUcap 	Capping point rotation in the previous +ve excursion
posKp 	Pre-capping slope in the previous +ve excursion
posKpc 	Post-capping slope in the previous +ve excursion
posFy 		Effective plastic moment in the previous +ve excursion
posFyProj  Projed effective plastic moment in the previous +ve excursion
posFcap        Maximum moment in the previous +ve excursion
posFcapProj Projed maximum  moment in the previous +ve excursion
negUy     Yielding rotation in the previous -ve excursion
negUcap   Capping point rotation in the previous -ve excursion
negKp     Pre-capping slope in the previous -ve excursion
negKpc    Post-capping slope in the previous -ve excursion
negFy         Effective plastic moment in the previous -ve excursion
negFyProj  Projed effective plastic moment in the previous -ve 
negFcap        Maximum moment in the previous -ve excursion
negFcapProj Projed maximum  moment in the previous -ve excursion
betaS    
betaC    
betaK    
betaF 	  
engRefS    	Refernence energy for strength deterioration
engRefC    	Refernence energy for post-capping strength deterioration
engRefK    	Refernence energy for unloading stiffness deterioration
engRefF    	Refernence energy for Smooth Transition deterioration
Excursion_Flag 		Flag for Excursion occurrence (i.e.,	crossing the x-axis)
Reversal_Flag 		Flag for Loading direction reversal occurrence
Yield_Flag 			Flag for Yielding occurrence
posFail_State 		Flag for reaching the ultimate rotation in the +ve loading direction
negFail_State 		Flag for reaching the ultimate rotation in the -ve loading direction
Mrpos_Flag 			Flag for reaching the residual moment in the +ve loading direction
Mrneg_Flag 			Flag for reaching the residual moment in the -ve loading direction
Failure_State 		Flag for reaching the reference energy
engDspt	Dissipated energy in previous excursion
Ei 		Dissipated energy in current excursion
engRvrs 			Total dissipated energy till previous load reversal point
engAcml 		Total dissipated energy till current step
*/
