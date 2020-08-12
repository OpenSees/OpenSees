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

//Modified Ibarra-Medina-Krawinkler with Peak-Oriented Hysteretic Response

//**********************************************************************                                                                     
// Code Developed by: Ahmed Elkady and Hammad ElJisr
// Last Updated: July 2020
//**********************************************************************

#ifndef IMKPinching_h
#define IMKPinching_h

#include <UniaxialMaterial.h>

class IMKPinching : public UniaxialMaterial
{
public:
	IMKPinching(int tag, double Ke,
		double Uy_pos, double Umax_pos, double Uu_pos, double Fy_pos, double FmaxFy_pos, double ResF_pos,
		double Uy_neg, double Umax_neg, double Uu_neg, double Fy_neg, double FmaxFy_neg, double ResF_neg,
		double LAMBDA_S, double LAMBDA_C, double LAMBDA_A, double LAMBDA_K, double c_S, double c_C, double c_A, double c_K, double D_pos, double D_neg, double kappaF, double kappaD);
	IMKPinching();
	~IMKPinching();
	const char *getClassType(void) const { return "IMKPinching"; };
	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);
	UniaxialMaterial *getCopy(void);
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	void Print(OPS_Stream &s, int flag = 0);


protected:

private:
	//my functions

	//Fixed input material parameters 	
	double Ke;
	double Up_pos;
	double Upc_pos;
	double Uu_pos;
	double Fy_pos;
	double FmaxFy_pos;
	double ResF_pos;
	double Up_neg;
	double Upc_neg;
	double Uu_neg;
	double Fy_neg;
	double FmaxFy_neg;
	double ResF_neg;
	double LAMBDA_S;
	double LAMBDA_C;
	double LAMBDA_A;
	double LAMBDA_K;
	double c_S;
	double c_C;
	double c_A;
	double c_K;
	double D_pos;
	double D_neg;
	double kappaF;
	double kappaD;
	
	//State variables 
	double U, cU;
	
	//History variables 

	double ui, 		cui;
	double fi, 		cfi;
	double ui_1,	cui_1;
	double fi_1, 	cfi_1;
	double du_i_1,	cdu_i_1;	

	double Uy_pos_j_1,				cUy_pos_j_1;
	double Umax_pos_j_1,			cUmax_pos_j_1;
	double Fy_pos_j_1,				cFy_pos_j_1;
	double Fmax_pos_j_1,			cFmax_pos_j_1;
	double Upeak_pos_j_1,			cUpeak_pos_j_1;
	double Fpeak_pos_j_1,			cFpeak_pos_j_1;
	double Ures_pos_j_1,			cUres_pos_j_1;
	double Fres_pos_j_1,			cFres_pos_j_1;
	double Kp_pos_j_1,				cKp_pos_j_1;
	double Kpc_pos_j_1,				cKpc_pos_j_1;


	double Uy_neg_j_1,               cUy_neg_j_1;
	double Umax_neg_j_1,             cUmax_neg_j_1;
	double Fy_neg_j_1,               cFy_neg_j_1;
	double Fmax_neg_j_1,             cFmax_neg_j_1;
	double Upeak_neg_j_1,            cUpeak_neg_j_1;
	double Fpeak_neg_j_1,            cFpeak_neg_j_1;

	double Ures_neg_j_1,             cUres_neg_j_1;
	double Fres_neg_j_1,             cFres_neg_j_1;
	double Kp_neg_j_1,               cKp_neg_j_1;
	double Kpc_neg_j_1,              cKpc_neg_j_1;

	double Kul_j_1, cKul_j_1;

	double Energy_Acc,  cEnergy_Acc;
	double Energy_Diss, cEnergy_Diss;
	
	double u0, cu0;

	double du;
	double df;
	
	double FailS;
	double FailC;
	double FailA;
	double FailK;

	double Ei, cEi;
	double dEi;
	double Epj;
	double EpjK;
	double EiK;

	double c_cS;
	double c_cC;
	double c_cA;
	double c_cK;

	double EtS;
	double EtC;
	double EtA;
	double EtK;

	double betaS;
	double betaC;
	double betaA;
	double betaK;

	double sPCsp,  sPCsn;
	double sPCpcp, sPCpcn;

	double TangentK, cTangentK, ki;

	double Uy_pos,		Uy_neg;
	double Umax_pos,	Umax_neg;
	double Fmax_pos,	Fmax_neg;
	double Kpc_pos,		Kpc_neg;
	double Kp_pos,		Kp_neg;
	
	double ULastPeak_pos_j_1,		cULastPeak_pos_j_1;
	double FLastPeak_pos_j_1,		cFLastPeak_pos_j_1;
	double ULastPeak_neg_j_1,		cULastPeak_neg_j_1;
	double FLastPeak_neg_j_1,		cFLastPeak_neg_j_1;
	
	double Failure_Flag,	cFailure_Flag;
	double Excursion_Flag,	cExcursion_Flag;
	double Reloading_Flag,	cReloading_Flag;
	double TargetPeak_Flag,	cTargetPeak_Flag;
	double Unloading_Flag,	cUnloading_Flag;
	double Yield_Flag,		cYield_Flag;
	double Reversal_Flag,   cReversal_Flag;

	double Krel_j_1,	    cKrel_j_1;	

	double Krel_LastPeak;
	double Krel_GlobalPeak;
	double K_check;

	double Upl,		cUpl;
	double Ubp,		cUbp;
	double Fbp,		cFbp;
	
};

#endif