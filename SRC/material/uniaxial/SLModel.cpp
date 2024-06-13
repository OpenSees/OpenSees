#include <math.h>
#include <elementAPI.h>
#include <Vector.h>
#include <Channel.h>
#include <float.h>
#include <iostream>


#include "SLModel.h"

#include <OPS_Globals.h>


static int numSLModel = 0;

void *
OPS_SLModel()
{
  // print out some KUDO's
  if (numSLModel == 0) {
    numSLModel++;
    //opserr << "SLModel version 2019.2\n";
	opserr << "SLModel version 2023.03\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  /*double dData[3];*/
  double dData[16];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  SLModel tag" << endln;
    return 0;
  }

  /*numData = 3;*/
  numData = 16;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    //opserr << "Invalid Args want: uniaxialMaterial SLModel tag? Dt?, sgm_ini?, OP_Material?";
	  opserr << "Invalid Args want: uniaxialMaterial SLModel tag? Dt?, E?, sigmaY0?, C?, gamma?, Qinf?, b?, sigmaC?, epsiC?, Ed1?, Ed2?,sigmaDM, aSigma?, aE?, lambda1Degrad?,cDegrad?";
    return 0;	
  }

  // create a new material
  /*theMaterial = new SLModel(iData[0], dData[0], dData[1], dData[2]);*/   
  theMaterial = new SLModel(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], dData[15]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type SLModel\n";
    return 0;
  }

  // return the material
  return theMaterial;
}


//MAT_TAG_SLModel or 0
//SLModel::SLModel(int tag, double Dt_temp, double sgm_ini_temp, double OP_Material_temp)
//:UniaxialMaterial(tag, MAT_TAG_SLModel), Dt(Dt_temp), sgm_ini(sgm_ini_temp), OP_Material(OP_Material_temp)
SLModel::SLModel(int tag, double Dt, double E, double sgm_ini, double c, double gamma, double q, double beta, double sigmaC, double epsiC, double Ed1, double Ed2, double sigmaDM,
	double aSigma, double aE, double lambda1Degrad, double cDegrad)
	:UniaxialMaterial(tag, MAT_TAG_SLModel), Dt(Dt), E(E), sgm_ini(sgm_ini), c(c), gamma(gamma), q(q), beta(beta), sigmaC(sigmaC), epsiC(epsiC), Ed1(Ed1),
	Ed2(Ed2), sigmaDM(sigmaDM), aSigma(aSigma), aE(aE), lambda1Degrad(lambda1Degrad), cDegrad(cDegrad)
{
	this->revertToStart();
}

SLModel::SLModel()
:UniaxialMaterial(0, MAT_TAG_SLModel), Dt(0.0), E(0.0), sgm_ini(0.0), c(0.0), gamma(0.0), q(0.0), beta(0.0), sigmaC(0.0), epsiC(0.0), Ed1(0.0),Ed2(0.0), sigmaDM(0.0),
aSigma(0.0), aE(0.0), lambda1Degrad(0.0), cDegrad(0.0)
{
	 this->revertToStart();
}

SLModel::~SLModel()
{
	// does nothing
}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
int 
SLModel::setTrialStrain(double strain, double strainRate)/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////TODO
{
	
	//all variables to the last commit
	this->revertToLastCommit(); // a = C_a 
	
	// local variables
	double True_tEpsp2, TrueEpsPrev, p_teps_d, Deltap_teps;
	double TrueSgmPinch;
	
	// set trial strain
	 neps = strain; //set trial displacement
	 teps = log(1.0+neps);
	 
	// ignore very small strain increment 
	double deltaD; 
	deltaD = strain - C_neps;
	if (fabs(deltaD) < 1.0e-18 && strain != 0.0) {
		return 0;
	}
	
	
	///////////////////////////////////////////////////////// not used
	if (iInitial <= 0) {
		iInitial++;
	}
	/////////////////////////////////////////////////////////
	
	
	// pre-buckling process
	if (status == 1) {
		if (teps <= teps_prev) {
			if (teps >= yteps_n) {
				status = 1;
				StrainHardeningFunc();
				
				//*** check capping stress ********************************************************************************************
				if (nsgm < cSgmc) {
					//deterioration slope modification
					cEpsc = neps_prev+(neps-neps_prev)/(nsgm-nsgm_prev)*(cSgmc-nsgm_prev);
					cEpsd1 = cEpsc-(cSgmc-cSgmd1)/cEd1;
					cSgmd2 = cSgmd1-cEd2*cEpsd1;
					cEpsd2 = -cSgmd2/cEd2;
					
					//peak point for tensile excursion
					tEpsp = cEpsc;
					tSgmp = sgm_ini;
					
					//reference point for back bone curve offset
					refEps = cEpsc-cIniEpsc; 
					cEpsy = refEps-sgm_ini/E;
					
					//status
					status = 4;
				}
				//**********************************************************************************************************************
			} else if (teps < yteps_n) {
				status = 3;
				StrainHardeningFunc();
				YieldPointFunc();
				
				//*** check capping stress ********************************************************************************************
				if (nsgm < cSgmc) {
					//deterioration slope modification
					cEpsc = neps_prev+(neps-neps_prev)/(nsgm-nsgm_prev)*(cSgmc-nsgm_prev);
					cEpsd1 = cEpsc-(cSgmc-cSgmd1)/cEd1;
					cSgmd2 = cSgmd1-cEd2*cEpsd1;
					cEpsd2 = -cSgmd2/cEd2;
					
					//peak point for tensile excursion
					tEpsp = cEpsc;
					tSgmp = sgm_ini;
					
					//reference point for back bone curve offset
					refEps = cEpsc-cIniEpsc; 
					cEpsy = refEps-sgm_ini/E;
					
					//status
					status = 4;
				}
				//**********************************************************************************************************************
			}
		} else { 
			if (teps <= yteps_p) {
				status = 1;
				StrainHardeningFunc();
			} else if (teps > yteps_p) {
				status = 2;
				StrainHardeningFunc();
				YieldPointFunc();
			}
		}
	} else if (status == 2) {
		if (teps >= teps_prev) {
			status = 2;
			StrainHardeningFunc();
			YieldPointFunc();
		} else if (teps < teps_prev && teps >= yteps_n) { 
			status = 1;
			StrainHardeningFunc();
		} else if (teps < teps_prev && teps < yteps_n) {
			status = 3;
			StrainHardeningFunc();
			YieldPointFunc();
		}  
	} else if (status == 3) {
		if (teps <= teps_prev) {
			status = 3;
			StrainHardeningFunc();
			YieldPointFunc();
				//*** check capping stress ********************************************************************************************
				if (nsgm < cSgmc) {
					//deterioration slope modification
					cEpsc = neps_prev+(neps-neps_prev)/(nsgm-nsgm_prev)*(cSgmc-nsgm_prev);
					cEpsd1 = cEpsc-(cSgmc-cSgmd1)/cEd1;
					cSgmd2 = cSgmd1-cEd2*cEpsd1;
					cEpsd2 = -cSgmd2/cEd2;
					
					//peak point for tensile excursion
					tEpsp = cEpsc;
					tSgmp = sgm_ini;
					
					//reference point for back bone curve offset
					refEps = cEpsc-cIniEpsc; 
					cEpsy = refEps-sgm_ini/E;
					
					//status
					status = 4;
				}
				//**********************************************************************************************************************
		} else if (teps > teps_prev && teps <= yteps_p) {
			status = 1;
			StrainHardeningFunc();
		} else if (teps > teps_prev && teps > yteps_p) {
			status = 2;
			StrainHardeningFunc();
			YieldPointFunc();
		} 
	}
	
	// post-buckling process
	if (status >= 4 && status <= 999) {
		
		if (neps < neps_prev) {
			if (cEpsy < neps) {
				//unloading
				status = 9;
				nsgm = cSgmy-cEu*(cEpsy-neps);
				//Tangent = cEu;
				Tangent = (nsgm-nsgm_prev)/(neps-neps_prev);
			} else if (cEpsc < neps && neps <= cEpsy) {
				//post-yield
				status = 10;
				nsgm = cSgmy+cEsth*(neps-cEpsy);
				//Tangent = cEsth;
				Tangent = (nsgm-nsgm_prev)/(neps-neps_prev);
			} else if (cEpsd1 < neps && neps <= cEpsc) {
				//1st deterioration slope
				status = 4;
				nsgm = cSgmc+cEd1*(neps-cEpsc);
				//Tangent = cEd1;
				Tangent = (nsgm-nsgm_prev)/(neps-neps_prev);
			} else if (cEpsd2 < neps && neps <= cEpsd1) {
				//2nd deterioration slope
				status = 5;
				nsgm = cEd2*neps+cSgmd2;
				//Tangent = cEd2;
				Tangent = (nsgm-nsgm_prev)/(neps-neps_prev);
			} else if (neps <= cEpsd2) {
				//zero stress    
				status = 1000;
				nsgm = -0.0001;
				Tangent = 1.0e-10;
				//std::cout << "*********************** 1111 **************************\n";
			}
			
			//plastic stain
			p_neps = neps-nsgm/tEu;
			
		} else if (neps > neps_prev) {
			if (neps < tEpsy) {
				//unloading
				status = 6;
				nsgm = tSgmy-tEu*(tEpsy-neps);
				//Tangent = tEu;
				Tangent = (nsgm-nsgm_prev)/(neps-neps_prev);
			} else if (tEpsy <= neps && neps < tEpsp2) {
				//reloading
				status = 7;
				nsgm = tSgmy+tEr2*(neps-tEpsy);
				//Tangent = tEr2;
				Tangent = (nsgm-nsgm_prev)/(neps-neps_prev);
			} else if (tEpsp2 <= neps) {
				//strain hardening
				status = 8;
				
				//true stress/strain at pinching point
				True_tEpsp2 = log(1.0+tEpsp2);
				
				//true stress/strain @ previous step
				TrueEpsPrev = log(1.0+neps_prev);
				
				//true strain @ current step
				teps = log(1.0+neps);
				
				//cumulative plastic strain
				if (TrueEpsPrev < True_tEpsp2) {
					 p_teps_d = teps-True_tEpsp2;
					tsgm = nsgm_prev*(1.0+tEpsp2);
				} else {
					p_teps_d = teps-TrueEpsPrev;
					tsgm = nsgm_prev*(1.0+neps_prev);
				}
				
				//increment of plastic strain 
				Deltap_teps = p_teps_d/5.0;
				
				for (int i=1;i<=5;i++) {
					//kinematic hardening
					alf_d = c/sgm_0*(tsgm-alf)*(Deltap_teps)-gamma*alf*(Deltap_teps);
					alf = alf+alf_d; //plus
					
					//isotropic hardening
					cum_p_teps = cum_p_teps+fabs(Deltap_teps);
					sgm_0 = sgm_ini+q*(1-exp(-1*beta*cum_p_teps));
					
					//true stress
					tsgm = alf+sgm_0; //plus
				}
				
				//engineering stress
				nsgm = tsgm/exp(teps);
				
				//Tangent
				Tangent = (nsgm-nsgm_prev)/(neps-neps_prev);
				
			}
			
			// plastic stain
			p_neps = neps-nsgm/cEu;
		}
		
	} else if (status == 1000) {
		status = 1000;
		nsgm = -0.00001;
		Tangent = 1.0e-10;
		//std::cout << "*********************** 2222 **************************\n";
	}
	
	
	//Cumulative plastic strain, Energy Dissipation
	if (status == 1 || status == 2 || status == 3) {
		DeltaE = 0.0;
	} else if (status == 1000) {
		DeltaE = 0.0;
	} else {
		DeltaE = fabs(p_neps-p_neps_prev)*fabs(nsgm+nsgm_prev)/2.0;
	}
	
	
	if (Et1-TotalE < 0.0) {
		Beta1 = 0.0;
		Alpha1 = Alpha1*(1-Beta1);
	} else if (DeltaE > Et1-TotalE) {
		Beta1 = 0.0;
		Alpha1 = Alpha1*(1-Beta1);
	} else {
		Beta1 = pow(DeltaE/(Et1-TotalE),c1);
		Alpha1 = Alpha1*(1-Beta1);
	}
	
	
	if (Et2-TotalE < 0.0) {
		Beta2 = 0.0;
		Alpha2 = Alpha2*(1-Beta2);
	} else if (DeltaE > Et2-TotalE) {
		Beta2 = 0.0;
		Alpha2 = Alpha2*(1-Beta2);
	} else {
		Beta2 = pow(DeltaE/(Et2-TotalE),c2);
		Alpha2 = Alpha2*(1-Beta2);
	}
	
	if (Et3-TotalE < 0.0) {
		Beta3 = 0.0;
		Alpha3 = Alpha3*(1-Beta3);
	} else if (DeltaE > Et3-TotalE) {
		Beta3 = 0.0;
		Alpha3 = Alpha3*(1-Beta3);
	} else {
		Beta3 = pow(DeltaE/(Et3-TotalE),c3);
		Alpha3 = Alpha3*(1-Beta3);
	}
	
	TotalE = TotalE+DeltaE;
	
	
	
	//Update back bone curve
	if (status == 4 || status == 5 || status == 10) {
		BackBoneTenFunc();
		
		//true stress @ pinching point
		TrueSgmPinch = tSgmp*(1.0+tEpsp2);
		//update kinematic hardening
		alf = TrueSgmPinch-sgm_0;
	} else if (status == 7 || status == 8) {
		BackBoneCompFunc();
	} else if (status == 9) {
		BackBoneTen2Func();
	} else if (status == 6) {
		BackBoneComp2Func();
	}
	
	//store stress and strain
	teps_prev = teps;
	neps_prev = neps;
	tsgm_prev = tsgm;
	nsgm_prev = nsgm;
	
	p_neps_prev = p_neps;
	p_teps_prev = p_teps;
	
	
	
	return 0;
}



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



double 
SLModel::getStrain(void)
{
  return neps;
}

double 
SLModel::getStress(void)
{
  return  nsgm;
}

double 
SLModel::getTangent(void)
{
	//std::cout << "*********" << Tangent << "****************************\n";
	return Tangent;
}

int 
SLModel::commitState(void)
{
	
	C_status = status;
	
	C_E = E;
	C_Dteq = Dteq;
	
	C_q = q;
	C_beta = beta;
	C_c = c;
	C_gamma = gamma;
	
	C_CapYieldStressM = CapYieldStressM;
	C_CapYieldStrainM = CapYieldStrainM;
	C_Ed1EM = Ed1EM;
	C_Ed2EM = Ed2EM;
	C_DetCapStressM = DetCapStressM;
	
	C_p_teps = p_teps;
	C_p_neps = p_neps;
	
	C_p_neps_prev = p_neps_prev;
	C_p_teps_prev = p_teps_prev;
	
	C_cum_p_teps = cum_p_teps;
	
	C_sgm_0 = sgm_0;
	C_alf_d = alf_d;
	C_alf = alf;
	
	C_ytsgm_p = ytsgm_p;
	C_ytsgm_n = ytsgm_n;
	C_yteps_p = yteps_p;
	C_yteps_n = yteps_n;
	
	C_teps = teps;
	C_neps = neps;
	C_tsgm = tsgm;
	C_nsgm = nsgm;
	
	C_teps_prev = teps_prev;
	C_neps_prev = neps_prev;
	C_tsgm_prev = tsgm_prev;
	C_nsgm_prev = nsgm_prev;
	
	C_cEu = cEu;
	C_cSgmy = cSgmy;
	C_cEpsy = cEpsy;
	C_cSgmc = cSgmc;
	C_cEpsc = cEpsc;
	C_cSgmd1 = cSgmd1;
	C_cEpsd1 = cEpsd1;
	C_cSgmd2 = cSgmd2;
	C_cEpsd2 = cEpsd2;
	C_cSgmb = cSgmb;
	C_cSgmd = cSgmd;
	C_cEsth = cEsth;
	C_cEd1 = cEd1;
	C_cEd2 = cEd2;
	
	C_cIniSgmy = cIniSgmy;
	C_cIniEpsy = cIniEpsy;
	C_cIniSgmc = cIniSgmc;
	C_cIniEpsc = cIniEpsc;
	C_cIniEsth = cIniEsth;
	C_cIniEd1 = cIniEd1;
	C_cIniEd2 = cIniEd2;
	C_cIniSgmd1 = cIniSgmd1;
	C_cIniEpsd1 = cIniEpsd1;
	C_cIniSgmb = cIniSgmb;
	C_cIniSgmd = cIniSgmd;
	C_cIniSgmd2 = cIniSgmd2;
	C_cIniEpsd2 = cIniEpsd2;
	
	C_tEu = tEu;
	C_tSgmy = tSgmy;
	C_tEpsy = tEpsy;
	C_tSgmp = tSgmp;
	C_tEpsp = tEpsp;
	C_tEpsp2 = tEpsp2;
	C_tEr = tEr;
	C_tEr2 = tEr2;
	C_refEps = refEps;
	
	C_ay = ay;
	C_au = au;
	
	C_Lambda1 = Lambda1;
	C_c1 = c1;
	C_Lambda2 = Lambda2;
	C_c2 = c2;
	C_Lambda3 = Lambda3;
	C_c3 = c3;
	
	C_Et1 = Et1;
	C_Et2 = Et2;
	C_Et3 = Et3;
	
	C_Beta1 = Beta1;
	C_Beta2 = Beta2;
	C_Beta3 = Beta3;
	
	C_Alpha1 = Alpha1;
	C_Alpha2 = Alpha2;
	C_Alpha3 = Alpha3;
	
	C_TotalE = TotalE;
	C_DeltaE = DeltaE;
	
	C_Tangent = Tangent;
	
	C_iInitial = iInitial;
	
	return 0;
}	


int 
SLModel::revertToLastCommit(void)
{
	//the opposite of commit trial history variables
	status = C_status;
	
	E = C_E;
	Dteq = C_Dteq;
	
	q = C_q;
	beta = C_beta;
	c = C_c;
	gamma = C_gamma;
	
	CapYieldStressM = C_CapYieldStressM;
	CapYieldStrainM = C_CapYieldStrainM;
	Ed1EM = C_Ed1EM;
	Ed2EM = C_Ed2EM;
	DetCapStressM = C_DetCapStressM;
	
	p_teps = C_p_teps;
	p_neps = C_p_neps;
	
	p_neps_prev = C_p_neps_prev;
	p_teps_prev = C_p_teps_prev;
	
	cum_p_teps = C_cum_p_teps;
	
	sgm_0 = C_sgm_0;
	alf_d = C_alf_d;
	alf = C_alf;
	
	ytsgm_p = C_ytsgm_p;
	ytsgm_n = C_ytsgm_n;
	yteps_p = C_yteps_p;
	yteps_n = C_yteps_n;
	
	teps = C_teps;
	neps = C_neps;
	tsgm = C_tsgm;
	nsgm = C_nsgm;
	
	teps_prev = C_teps_prev;
	neps_prev = C_neps_prev;
	tsgm_prev = C_tsgm_prev;
	nsgm_prev = C_nsgm_prev;
	
	cEu = C_cEu;
	cSgmy = C_cSgmy;
	cEpsy = C_cEpsy;
	cSgmc = C_cSgmc;
	cEpsc = C_cEpsc;
	cSgmd1 = C_cSgmd1;
	cEpsd1 = C_cEpsd1;
	cSgmd2 = C_cSgmd2;
	cEpsd2 = C_cEpsd2;
	cSgmb = C_cSgmb;
	cSgmd = C_cSgmd;
	cEsth = C_cEsth;
	cEd1 = C_cEd1;
	cEd2 = C_cEd2;
	
	cIniSgmy = C_cIniSgmy;
	cIniEpsy = C_cIniEpsy;
	cIniSgmc = C_cIniSgmc;
	cIniEpsc = C_cIniEpsc;
	cIniEsth = C_cIniEsth;
	cIniEd1 = C_cIniEd1;
	cIniEd2 = C_cIniEd2;
	cIniSgmd1 = C_cIniSgmd1;
	cIniEpsd1 = C_cIniEpsd1;
	cIniSgmb = C_cIniSgmb;
	cIniSgmd = C_cIniSgmd;
	cIniSgmd2 = C_cIniSgmd2;
	cIniEpsd2 = C_cIniEpsd2;
	
	tEu = C_tEu;
	tSgmy = C_tSgmy;
	tEpsy = C_tEpsy;
	tSgmp = C_tSgmp;
	tEpsp = C_tEpsp;
	tEpsp2 = C_tEpsp2;
	tEr = C_tEr;
	tEr2 = C_tEr2;
	refEps = C_refEps;
	
	ay = C_ay;
	au = C_au;
	
	Lambda1 = C_Lambda1;
	c1 = C_c1;
	Lambda2 = C_Lambda2;
	c2 = C_c2;
	Lambda3 = C_Lambda3;
	c3 = C_c3;
	
	Et1 = C_Et1;
	Et2 = C_Et2;
	Et3 = C_Et3;
	
	Beta1 = C_Beta1;
	Beta2 = C_Beta2;
	Beta3 = C_Beta3;
	
	Alpha1 = C_Alpha1;
	Alpha2 = C_Alpha2;
	Alpha3 = C_Alpha3;
	
	TotalE = C_TotalE;
	DeltaE = C_DeltaE;
	
	Tangent = C_Tangent;
	
	iInitial = C_iInitial;

  return 0;
}


int 
SLModel::revertToStart(void)
{
	// Initialize state variables
	status = C_status = 1;
	
	//E = C_E = 200000.0;
	C_E = E;
	
	Dteq = C_Dteq = (Dt)*sqrt(sgm_ini/E);
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//if (OP_Material == 10) {
	//	//*** BCP ***
	//	q = C_q = 68.0956;
	//	beta = C_beta = 14.072;
	//	c = C_c = 4200.0;
	//	gamma = C_gamma = 19.0;
	//	
	//	CapYieldStressM = C_CapYieldStressM = 1.2345*pow(Dteq,-0.564);
	//	CapYieldStrainM = C_CapYieldStrainM = 8.3096*pow(Dteq,-3.156);
	//	Ed1EM = C_Ed1EM = -0.0267*pow(Dteq,1.5213)*0.8;//**********************************************************************************************************************************************************
	//	Ed2EM = C_Ed2EM = -0.0034*pow(Dteq,0.8134)*0.8;//**********************************************************************************************************************************************************
	//	DetCapStressM = C_DetCapStressM = 0.6225*pow(Dteq,0.0888);
	//	
	//	ay = C_ay = 1.6465*pow(Dteq/sqrt(sgm_ini/E),-1.136);
	//	au = C_au = 5.5813*pow(Dteq/sqrt(sgm_ini/E),-1.732);
	//	
	//	Lambda1 = C_Lambda1 = 31.237*pow(Dt,-1.127);
	//	c1 = C_c1 = 1.0;
	//	
	//	Lambda2 = C_Lambda2 = Lambda1*3.0; //**********************************************************************************************************************************************************
	//	c2 = C_c2 = 1.0;
	//	
	//	Lambda3 = C_Lambda3 = Lambda1*3.0; //**********************************************************************************************************************************************************
	//	//Lambda3 = C_Lambda3 = Lambda1*1000.0; //**********************************************************************************************************************************************************
	//	c3 = C_c3 = 1.0;
	//
	//} else if (OP_Material == 11) {
	//	//*** HYP ***
	//	q = C_q = 68.1;
	//	beta = C_beta = 14.1;
	//	c = C_c = 2100.0;
	//	gamma = C_gamma = 11.0;
	//	
	//	CapYieldStressM = C_CapYieldStressM = 1.157*pow(Dteq,-0.36);
	//	CapYieldStrainM = C_CapYieldStrainM = 8.3096*pow(Dteq,-3.156);
	//	Ed1EM = C_Ed1EM = -0.027*pow(Dteq,2.0)*0.8;//**********************************************************************************************************************************************************
	//	Ed2EM = C_Ed2EM = -0.0034*pow(Dteq,0.8134)*0.8;//**********************************************************************************************************************************************************
	//	DetCapStressM = C_DetCapStressM = 0.575*pow(Dteq,-0.03);
	//	
	//	ay = C_ay = 2.355*pow(Dteq/sqrt(sgm_ini/E),-1.204)*1.0;
	//	au = C_au = 5.2672*pow(Dteq/sqrt(sgm_ini/E),-1.683);
	//	
	//	Lambda1 = C_Lambda1 = 31.237*pow(Dt,-1.127);
	//	c1 = C_c1 = 1.0;
	//	
	//	Lambda2 = C_Lambda2 = Lambda1*3.0; //**********************************************************************************************************************************************************
	//	c2 = C_c2 = 1.0;
	//	
	//	Lambda3 = C_Lambda3 = Lambda1*3.0; //**********************************************************************************************************************************************************
	//	//Lambda3 = C_Lambda3 = Lambda1*1000.0; //**********************************************************************************************************************************************************
	//	c3 = C_c3 = 1.0;
	//
	//
	//} else if (OP_Material == 12) {
	//	//*** BCR ***
	//	q = C_q = 22.4;
	//	beta = C_beta = 7.2;
	//	c = C_c = 2500.0;
	//	gamma = C_gamma = 19.0;
	//	
	//	//CapYieldStressM = C_CapYieldStressM = 1.114*pow(Dteq,-0.314);//**********************************************************************************************************************************************************
	//	CapYieldStressM = C_CapYieldStressM = 1.125*pow(Dteq,-0.27);//**********************************************************************************************************************************************************
	//	CapYieldStrainM = C_CapYieldStrainM = 6.720*pow(Dteq,-3.09);
	//	//Ed1EM = C_Ed1EM = -0.015*pow(Dteq,1.5);//**********************************************************************************************************************************************************
	//	Ed1EM = C_Ed1EM = -0.0133*pow(Dteq,2.0)*1.0;//**********************************************************************************************************************************************************
	//	Ed2EM = C_Ed2EM = -0.003*pow(Dteq,1.221)*0.8;//**********************************************************************************************************************************************************
	//	DetCapStressM = C_DetCapStressM = 0.614*pow(Dteq,-0.041);
	//	
	//	ay = C_ay = 24.0*pow(Dteq/sqrt(sgm_ini/E),-1.8);
	//	au = C_au = 5.4522*pow(Dteq/sqrt(sgm_ini/E),-1.653)*1.0;
	//	
	//	Lambda1 = C_Lambda1 = 17.72*pow(Dt,-0.989)*1.0;
	//	c1 = C_c1 = 1.0;
	//	
	//	Lambda2 = C_Lambda2 = Lambda1*3.0; //**********************************************************************************************************************************************************
	//	c2 = C_c2 = 1.0;
	//	
	//	Lambda3 = C_Lambda3 = Lambda1*3.0; //**********************************************************************************************************************************************************
	//	//Lambda3 = C_Lambda3 = Lambda1*1000.0; //**********************************************************************************************************************************************************
	//	c3 = C_c3 = 1.0;
	//	
	//} else {
	//	std::cout << "*********************** Material ID " << OP_Material << " is not prepared. **************************\n";
	//}
	
	C_q = q;
	C_beta = beta;
	C_c = c;
	C_gamma = gamma;

	sigmaCDivSigmaY = sigmaC / sgm_ini;
	epsiCDivEpsiY = epsiC / (sgm_ini / E);
	Ed1DivE = Ed1 / E;
	Ed2DivE = Ed2 / E;
	sigmaDMDivSigmaC = sigmaDM / sigmaC;

	CapYieldStressM = sigmaCDivSigmaY;
	CapYieldStrainM = epsiCDivEpsiY;
	Ed1EM = C_Ed1EM = Ed1DivE;//**********************************************************************************************************************************************************
	Ed2EM = C_Ed2EM = Ed2DivE;//**********************************************************************************************************************************************************
	DetCapStressM = C_DetCapStressM = sigmaDMDivSigmaC;

	ay = C_ay = aSigma;
	au = C_au = aE;

	Lambda1 = C_Lambda1 = lambda1Degrad;
	c1 = C_c1 = cDegrad;

	Lambda2 = C_Lambda2 = Lambda1 * 3.0; //**********************************************************************************************************************************************************
	c2 = C_c2 = cDegrad;

	Lambda3 = C_Lambda3 = Lambda1 * 3.0; //**********************************************************************************************************************************************************
	//Lambda3 = C_Lambda3 = Lambda1*1000.0; //**********************************************************************************************************************************************************
	c3 = C_c3 = cDegrad;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	p_teps = C_p_teps = 0.0;
	p_neps = C_p_neps = 0.0;
	
	p_neps_prev = C_p_neps_prev = 0.0;
	p_teps_prev = C_p_teps_prev = 0.0;
	
	cum_p_teps = C_cum_p_teps = 0.0;
	
	sgm_0 = C_sgm_0 = sgm_ini;
	alf_d = C_alf_d = 0.0;
	alf = C_alf = 0.0;
	
	ytsgm_p = C_ytsgm_p = sgm_ini;
	ytsgm_n = C_ytsgm_n = -1*sgm_ini;
	yteps_p = C_yteps_p = ytsgm_p/E;;
	yteps_n = C_yteps_n = ytsgm_n/E;;
	
	teps = C_teps = 0.0;
	neps = C_neps = 0.0;
	tsgm = C_tsgm = 0.0;
	nsgm = C_nsgm = 0.0;
	
	teps_prev = C_teps_prev = 0.0;
	neps_prev = C_neps_prev = 0.0;
	tsgm_prev = C_tsgm_prev = 0.0;
	nsgm_prev = C_nsgm_prev = 0.0;
	
	cIniSgmy = C_cIniSgmy = (-1)*sgm_ini;
	cIniEpsy = C_cIniEpsy = cIniSgmy/E;
	cIniSgmc = C_cIniSgmc = (-1)*CapYieldStressM*sgm_ini;
	
	
	
	cIniEpsc = C_cIniEpsc = (-1)*CapYieldStrainM*(sgm_ini/E);
	cIniEsth = C_cIniEsth = (cIniSgmc-cIniSgmy)/(cIniEpsc-cIniEpsy);
	cIniEd1 = C_cIniEd1 = Ed1EM*E;
	cIniEd2 = C_cIniEd2 = Ed2EM*E;
	cIniSgmd1 = C_cIniSgmd1 = DetCapStressM*cIniSgmc;
	cIniEpsd1 = C_cIniEpsd1 = cIniEpsc-(cIniSgmc-cIniSgmd1)/cIniEd1;
	cIniSgmb = C_cIniSgmb = cIniSgmy-cIniEsth*cIniEpsy;
	cIniSgmd = C_cIniSgmd = cIniSgmc-cIniEd1*cIniEpsc;
	cIniSgmd2 = C_cIniSgmd2 = cIniSgmd1-cIniEd2*cIniEpsd1;
	cIniEpsd2 = C_cIniEpsd2 = -cIniSgmd2/cIniEd2;
	
	cSgmb = C_cSgmb = cIniSgmb; 
	cSgmd = C_cSgmd = cIniSgmd;
	
	cEu = C_cEu = E;
	cEsth = C_cEsth = cIniEsth;
	cEd1 = C_cEd1 = cIniEd1;
	cEd2 = C_cEd2 = cIniEd2;
	
	
	cEpsc = C_cEpsc = -(cSgmd-cSgmb)/(cIniEd1-cIniEsth);
	cSgmc = C_cSgmc = cIniEd1*cEpsc+cSgmd;
	
	cSgmd1 = C_cSgmd1 = cIniSgmd1;
	cEpsd1 = C_cEpsd1 = cEpsc-(cSgmc-cSgmd1)/cEd1;
	cSgmd2 = C_cSgmd2 = cSgmd1-cEd2*cEpsd1;
	cEpsd2 = C_cEpsd2 = -cSgmd2/cEd2;
	
	cSgmd = C_cSgmd = cSgmc-cEd1*cEpsc;
	cSgmb = C_cSgmb = cSgmc-cEsth*cEpsc;
	
	cSgmy = C_cSgmy = cEu*cSgmb/(cEu-cEsth);
	cEpsy = C_cEpsy = cSgmb/(cEu-cEsth);

	
	tEu = C_tEu = E;
	tSgmy = C_tSgmy = sgm_ini;
	tEpsy = C_tEpsy = sgm_ini/E;
	tSgmp = C_tSgmp = sgm_ini;
	tEpsp = C_tEpsp = cEpsc;
	tEpsp2 = C_tEpsp2 = cEpsc;
	
	tEr = C_tEr = E;
	tEr2 = C_tEr2 = E;
	
	refEps = C_refEps = 0.0;
	
	
	Et1 = C_Et1 = Lambda1*sgm_ini;
	Et2 = C_Et2 = Lambda2*sgm_ini; 
	Et3 = C_Et3 = Lambda3*sgm_ini;
	
	Beta1 = C_Beta1 = 0.0;
	Beta2 = C_Beta2 = 0.0;
	Beta3 = C_Beta3 = 0.0;
	
	Alpha1 = C_Alpha1 = 1.0;
	Alpha2 = C_Alpha2 = 1.0;
	Alpha3 = C_Alpha3 = 1.0;
	
	TotalE = C_TotalE = 0.0;
	DeltaE = C_DeltaE = 0.0;
	
	Tangent = C_Tangent = 0.0;
	
	iInitial = C_iInitial = 0.0;

  return 0;
}


UniaxialMaterial *
SLModel::getCopy(void)
{
	SLModel *theCopy = new SLModel(this->getTag(), Dt, E, sgm_ini, c, gamma, q, beta, sigmaC, epsiC, Ed1, Ed2, sigmaDM, aSigma, aE, lambda1Degrad, cDegrad);

	//Fixed Model parameters: need to change according to material properties
	theCopy -> status = status;
	
	theCopy -> E = E;
	theCopy -> Dteq = Dteq;
	theCopy -> q = q;
	theCopy -> beta = beta;
	theCopy -> c = c;
	theCopy -> gamma = gamma;
	theCopy -> CapYieldStressM = CapYieldStressM;
	theCopy -> CapYieldStrainM = CapYieldStrainM;
	theCopy -> Ed1EM = Ed1EM;
	theCopy -> Ed2EM = Ed2EM;
	theCopy -> DetCapStressM = DetCapStressM;
	theCopy -> p_teps = p_teps;
	theCopy -> p_neps = p_neps;
	theCopy -> p_neps_prev = p_neps_prev;
	theCopy -> p_teps_prev = p_teps_prev;
	theCopy -> cum_p_teps = cum_p_teps;
	theCopy -> sgm_0 = sgm_0;
	theCopy -> alf_d = alf_d;
	theCopy -> alf = alf;
	theCopy -> ytsgm_p = ytsgm_p;
	theCopy -> ytsgm_n = ytsgm_n;
	theCopy -> yteps_p = yteps_p;
	theCopy -> yteps_n = yteps_n;
	theCopy -> teps = teps;
	theCopy -> neps = neps;
	theCopy -> tsgm = tsgm;
	theCopy -> nsgm = nsgm;
	theCopy -> teps_prev = teps_prev;
	theCopy -> neps_prev = neps_prev;
	theCopy -> tsgm_prev = tsgm_prev;
	theCopy -> nsgm_prev = nsgm_prev;
	theCopy -> cEu = cEu;
	theCopy -> cSgmy = cSgmy;
	theCopy -> cEpsy = cEpsy;
	theCopy -> cSgmc = cSgmc;
	theCopy -> cEpsc = cEpsc;
	theCopy -> cSgmd1 = cSgmd1;
	theCopy -> cEpsd1 = cEpsd1;
	theCopy -> cSgmd2 = cSgmd2;
	theCopy -> cEpsd2 = cEpsd2;
	theCopy -> cSgmb = cSgmb;
	theCopy -> cSgmd = cSgmd;
	theCopy -> cEsth = cEsth;
	theCopy -> cEd1 = cEd1;
	theCopy -> cEd2 = cEd2;
	theCopy -> cIniSgmy = cIniSgmy;
	theCopy -> cIniEpsy = cIniEpsy;
	theCopy -> cIniSgmc = cIniSgmc;
	theCopy -> cIniEpsc = cIniEpsc;
	theCopy -> cIniEsth = cIniEsth;
	theCopy -> cIniEd1 = cIniEd1;
	theCopy -> cIniEd2 = cIniEd2;
	theCopy -> cIniSgmd1 = cIniSgmd1;
	theCopy -> cIniEpsd1 = cIniEpsd1;
	theCopy -> cIniSgmb = cIniSgmb;
	theCopy -> cIniSgmd = cIniSgmd;
	theCopy -> cIniSgmd2 = cIniSgmd2;
	theCopy -> cIniEpsd2 = cIniEpsd2;
	theCopy -> tEu = tEu;
	theCopy -> tSgmy = tSgmy;
	theCopy -> tEpsy = tEpsy;
	theCopy -> tSgmp = tSgmp;
	theCopy -> tEpsp = tEpsp;
	theCopy -> tEpsp2 = tEpsp2;
	theCopy -> tEr = tEr;
	theCopy -> tEr2 = tEr2;
	theCopy -> refEps = refEps;
	theCopy -> ay = ay;
	theCopy -> au = au;
	theCopy -> Lambda1 = Lambda1;
	theCopy -> c1 = c1;
	theCopy -> Lambda2 = Lambda2;
	theCopy -> c2 = c2;
	theCopy -> Lambda3 = Lambda3;
	theCopy -> c3 = c3;
	theCopy -> Et1 = Et1;
	theCopy -> Et2 = Et2;
	theCopy -> Et3 = Et3;
	theCopy -> Beta1 = Beta1;
	theCopy -> Beta2 = Beta2;
	theCopy -> Beta3 = Beta3;
	theCopy -> Alpha1 = Alpha1;
	theCopy -> Alpha2 = Alpha2;
	theCopy -> Alpha3 = Alpha3;
	theCopy -> TotalE = TotalE;
	theCopy -> DeltaE = DeltaE;
	theCopy -> Tangent = Tangent;
	theCopy -> iInitial = iInitial;
	
	theCopy -> C_E = C_E;
	theCopy -> C_Dteq = C_Dteq;
	theCopy -> C_q = C_q;
	theCopy -> C_beta = C_beta;
	theCopy -> C_c = C_c;
	theCopy -> C_gamma = C_gamma;
	theCopy -> C_CapYieldStressM = C_CapYieldStressM;
	theCopy -> C_CapYieldStrainM = C_CapYieldStrainM;
	theCopy -> C_Ed1EM = C_Ed1EM;
	theCopy -> C_Ed2EM = C_Ed2EM;
	theCopy -> C_DetCapStressM = C_DetCapStressM;
	theCopy -> C_p_teps = C_p_teps;
	theCopy -> C_p_neps = C_p_neps;
	theCopy -> C_p_neps_prev = C_p_neps_prev;
	theCopy -> C_p_teps_prev = C_p_teps_prev;
	theCopy -> C_cum_p_teps = C_cum_p_teps;
	theCopy -> C_sgm_0 = C_sgm_0;
	theCopy -> C_alf_d = C_alf_d;
	theCopy -> C_alf = C_alf;
	theCopy -> C_ytsgm_p = C_ytsgm_p;
	theCopy -> C_ytsgm_n = C_ytsgm_n;
	theCopy -> C_yteps_p = C_yteps_p;
	theCopy -> C_yteps_n = C_yteps_n;
	theCopy -> C_teps = C_teps;
	theCopy -> C_neps = C_neps;
	theCopy -> C_tsgm = C_tsgm;
	theCopy -> C_nsgm = C_nsgm;
	theCopy -> C_teps_prev = C_teps_prev;
	theCopy -> C_neps_prev = C_neps_prev;
	theCopy -> C_tsgm_prev = C_tsgm_prev;
	theCopy -> C_nsgm_prev = C_nsgm_prev;
	theCopy -> C_cEu = C_cEu;
	theCopy -> C_cSgmy = C_cSgmy;
	theCopy -> C_cEpsy = C_cEpsy;
	theCopy -> C_cSgmc = C_cSgmc;
	theCopy -> C_cEpsc = C_cEpsc;
	theCopy -> C_cSgmd1 = C_cSgmd1;
	theCopy -> C_cEpsd1 = C_cEpsd1;
	theCopy -> C_cSgmd2 = C_cSgmd2;
	theCopy -> C_cEpsd2 = C_cEpsd2;
	theCopy -> C_cSgmb = C_cSgmb;
	theCopy -> C_cSgmd = C_cSgmd;
	theCopy -> C_cEsth = C_cEsth;
	theCopy -> C_cEd1 = C_cEd1;
	theCopy -> C_cEd2 = C_cEd2;
	theCopy -> C_cIniSgmy = C_cIniSgmy;
	theCopy -> C_cIniEpsy = C_cIniEpsy;
	theCopy -> C_cIniSgmc = C_cIniSgmc;
	theCopy -> C_cIniEpsc = C_cIniEpsc;
	theCopy -> C_cIniEsth = C_cIniEsth;
	theCopy -> C_cIniEd1 = C_cIniEd1;
	theCopy -> C_cIniEd2 = C_cIniEd2;
	theCopy -> C_cIniSgmd1 = C_cIniSgmd1;
	theCopy -> C_cIniEpsd1 = C_cIniEpsd1;
	theCopy -> C_cIniSgmb = C_cIniSgmb;
	theCopy -> C_cIniSgmd = C_cIniSgmd;
	theCopy -> C_cIniSgmd2 = C_cIniSgmd2;
	theCopy -> C_cIniEpsd2 = C_cIniEpsd2;
	theCopy -> C_tEu = C_tEu;
	theCopy -> C_tSgmy = C_tSgmy;
	theCopy -> C_tEpsy = C_tEpsy;
	theCopy -> C_tSgmp = C_tSgmp;
	theCopy -> C_tEpsp = C_tEpsp;
	theCopy -> C_tEpsp2 = C_tEpsp2;
	theCopy -> C_tEr = C_tEr;
	theCopy -> C_tEr2 = C_tEr2;
	theCopy -> C_refEps = C_refEps;
	theCopy -> C_ay = C_ay;
	theCopy -> C_au = C_au;
	theCopy -> C_Lambda1 = C_Lambda1;
	theCopy -> C_c1 = C_c1;
	theCopy -> C_Lambda2 = C_Lambda2;
	theCopy -> C_c2 = C_c2;
	theCopy -> C_Lambda3 = C_Lambda3;
	theCopy -> C_c3 = C_c3;
	theCopy -> C_Et1 = C_Et1;
	theCopy -> C_Et2 = C_Et2;
	theCopy -> C_Et3 = C_Et3;
	theCopy -> C_Beta1 = C_Beta1;
	theCopy -> C_Beta2 = C_Beta2;
	theCopy -> C_Beta3 = C_Beta3;
	theCopy -> C_Alpha1 = C_Alpha1;
	theCopy -> C_Alpha2 = C_Alpha2;
	theCopy -> C_Alpha3 = C_Alpha3;
	theCopy -> C_TotalE = C_TotalE;
	theCopy -> C_DeltaE = C_DeltaE;
	theCopy -> C_Tangent = C_Tangent;
	theCopy -> C_iInitial = C_iInitial;
	
	return theCopy;
}


int 
SLModel::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	static Vector data(182);
	data(0) = this->getTag();
	
	// constaint variables
	data(1) = Dt;
	data(2) = sgm_ini;
	//data(3) = OP_Material;
	data(3) = 0.;
	
	// constaint variables
	data(4) = E;
	data(5) = Dteq;
	data(6) = q;
	data(7) = beta;
	data(8) = c;
	data(9) = gamma;
	data(10) = CapYieldStressM;
	data(11) = CapYieldStrainM;
	data(12) = Ed1EM;
	data(13) = Ed2EM;
	data(14) = DetCapStressM;
	data(15) = ay;
	data(16) = au;
	data(17) = Lambda1;
	data(18) = c1;
	data(19) = Lambda2;
	data(20) = c2;
	data(21) = Lambda3;
	data(22) = c3;
	data(23) = Et1;
	data(24) = Et2;
	data(25) = Et3;
	
	// State variables from last converged state
	data(26) = status;
	data(27) = p_teps;
	data(28) = p_neps;
	data(29) = p_neps_prev;
	data(30) = p_teps_prev;
	data(31) = cum_p_teps;
	data(31) = sgm_0;
	data(33) = alf_d;
	data(34) = alf;
	data(35) = ytsgm_p;
	data(36) = ytsgm_n;
	data(37) = yteps_p;
	data(38) = yteps_n;
	data(39) = teps;
	data(40) = neps;
	data(41) = tsgm;
	data(42) = nsgm;
	data(43) = teps_prev;
	data(44) = neps_prev;
	data(45) = tsgm_prev;
	data(46) = nsgm_prev;
	data(47) = cEu;
	data(48) = cSgmy;
	data(49) = cEpsy;
	data(50) = cSgmc;
	data(51) = cEpsc;
	data(52) = cSgmd1;
	data(53) = cEpsd1;
	data(54) = cSgmd2;
	data(55) = cEpsd2;
	data(56) = cSgmb;
	data(57) = cSgmd;
	data(58) = cEsth;
	data(59) = cEd1;
	data(60) = cEd2;
	data(61) = cIniSgmy;
	data(62) = cIniEpsy;
	data(63) = cIniSgmc;
	data(64) = cIniEpsc;
	data(65) = cIniEsth;
	data(66) = cIniEd1;
	data(67) = cIniEd2;
	data(68) = cIniSgmd1;
	data(69) = cIniEpsd1;
	data(70) = cIniSgmb;
	data(71) = cIniSgmd;
	data(72) = cIniSgmd2;
	data(73) = cIniEpsd2;
	data(74) = tEu;
	data(75) = tSgmy;
	data(76) = tEpsy;
	data(77) = tSgmp;
	data(78) = tEpsp;
	data(79) = tEpsp2;
	data(80) = tEr;
	data(81) = tEr2;
	data(82) = refEps;
	data(83) = Beta1;
	data(84) = Beta2;
	data(85) = Beta3;
	data(86) = Alpha1;
	data(87) = Alpha2;
	data(88) = Alpha3;
	data(89) = TotalE;
	data(90) = DeltaE;
	data(91) = Tangent;
	
	// constaint variables
	data(92) = C_E;
	data(93) = C_Dteq;
	data(94) = C_q;
	data(95) = C_beta;
	data(96) = C_c;
	data(97) = C_gamma;
	data(98) = C_CapYieldStressM;
	data(99) = C_CapYieldStrainM;
	data(100) = C_Ed1EM;
	data(101) = C_Ed2EM;
	data(102) = C_DetCapStressM;
	data(103) = C_ay;
	data(104) = C_au;
	data(105) = C_Lambda1;
	data(106) = C_c1;
	data(107) = C_Lambda2;
	data(108) = C_c2;
	data(109) = C_Lambda3;
	data(110) = C_c3;
	data(111) = C_Et1;
	data(112) = C_Et2;
	data(113) = C_Et3;
	
	// State variables from last converged state
	data(114) = C_status;
	data(115) = C_p_teps;
	data(116) = C_p_neps;
	data(117) = C_p_neps_prev;
	data(118) = C_p_teps_prev;
	data(119) = C_cum_p_teps;
	data(120) = C_sgm_0;
	data(121) = C_alf_d;
	data(122) = C_alf;
	data(123) = C_ytsgm_p;
	data(124) = C_ytsgm_n;
	data(125) = C_yteps_p;
	data(126) = C_yteps_n;
	data(127) = C_teps;
	data(128) = C_neps;
	data(129) = C_tsgm;
	data(130) = C_nsgm;
	data(131) = C_teps_prev;
	data(132) = C_neps_prev;
	data(133) = C_tsgm_prev;
	data(134) = C_nsgm_prev;
	data(135) = C_cEu;
	data(136) = C_cSgmy;
	data(137) = C_cEpsy;
	data(138) = C_cSgmc;
	data(139) = C_cEpsc;
	data(140) = C_cSgmd1;
	data(141) = C_cEpsd1;
	data(142) = C_cSgmd2;
	data(143) = C_cEpsd2;
	data(144) = C_cSgmb;
	data(145) = C_cSgmd;
	data(146) = C_cEsth;
	data(147) = C_cEd1;
	data(148) = C_cEd2;
	data(149) = C_cIniSgmy;
	data(150) = C_cIniEpsy;
	data(151) = C_cIniSgmc;
	data(152) = C_cIniEpsc;
	data(153) = C_cIniEsth;
	data(154) = C_cIniEd1;
	data(155) = C_cIniEd2;
	data(156) = C_cIniSgmd1;
	data(157) = C_cIniEpsd1;
	data(158) = C_cIniSgmb;
	data(159) = C_cIniSgmd;
	data(160) = C_cIniSgmd2;
	data(161) = C_cIniEpsd2;
	data(162) = C_tEu;
	data(163) = C_tSgmy;
	data(164) = C_tEpsy;
	data(165) = C_tSgmp;
	data(166) = C_tEpsp;
	data(167) = C_tEpsp2;
	data(168) = C_tEr;
	data(169) = C_tEr2;
	data(170) = C_refEps;
	data(171) = C_Beta1;
	data(172) = C_Beta2;
	data(173) = C_Beta3;
	data(174) = C_Alpha1;
	data(175) = C_Alpha2;
	data(176) = C_Alpha3;
	data(177) = C_TotalE;
	data(178) = C_DeltaE;
	data(179) = C_Tangent;
	
	data(180) = iInitial;
	data(181) = C_iInitial;
	
	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0) 
		opserr << "SLModel::sendSelf() - failed to send data\n";

	return res;
}

int 
SLModel::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(182);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "SLModel::recvSelf() - failed to recv data\n";
  else {
  	
  	this->setTag(data(0));
	
	// constaint variables
	 Dt= data(1);;
	 sgm_ini = data(2);;
	 //OP_Material = data(3);
	
	// constaint variables
	 E = data(4);
	 Dteq = data(5);
	 q = data(6);
	 beta = data(7);
	 c = data(8);
	 gamma = data(9);
	 CapYieldStressM = data(10);
	 CapYieldStrainM = data(11);
	 Ed1EM = data(12);
	 Ed2EM = data(13);
	 DetCapStressM = data(14);
	 ay = data(15);
	 au = data(16);
	 Lambda1 = data(17);
	 c1 = data(18);
	 Lambda2 = data(19);
	 c2 = data(20);
	 Lambda3 = data(21);
	 c3 = data(22);
	 Et1 = data(23);
	 Et2 = data(24);
	 Et3 = data(25);
	
	// State variables from last converged state
	 status = data(26);
	 p_teps = data(27);
	 p_neps = data(28);
	 p_neps_prev = data(29);
	 p_teps_prev = data(30);
	 cum_p_teps = data(31);
	 sgm_0 = data(31);
	 alf_d = data(33);
	 alf = data(34);
	 ytsgm_p = data(35);
	 ytsgm_n = data(36);
	 yteps_p = data(37);
	 yteps_n = data(38);
	 teps = data(39);
	 neps = data(40);
	 tsgm = data(41);
	 nsgm = data(42);
	 teps_prev = data(43);
	 neps_prev = data(44);
	 tsgm_prev = data(45);
	 nsgm_prev = data(46);
	 cEu = data(47);
	 cSgmy = data(48);
	 cEpsy = data(49);
	 cSgmc = data(50);
	 cEpsc = data(51);
	 cSgmd1 = data(52);
	 cEpsd1 = data(53);
	 cSgmd2 = data(54);
	 cEpsd2 = data(55);
	 cSgmb = data(56);
	 cSgmd = data(57);
	 cEsth = data(58);
	 cEd1 = data(59);
	 cEd2 = data(60);
	 cIniSgmy = data(61);
	 cIniEpsy = data(62);
	 cIniSgmc = data(63);
	 cIniEpsc = data(64);
	 cIniEsth = data(65);
	 cIniEd1 = data(66);
	 cIniEd2 = data(67);
	 cIniSgmd1 = data(68);
	 cIniEpsd1 = data(69);
	 cIniSgmb = data(70);
	 cIniSgmd = data(71);
	 cIniSgmd2 = data(72);
	 cIniEpsd2= data(73);
	 tEu = data(74);
	 tSgmy = data(75);
	 tEpsy = data(76);
	 tSgmp = data(77);
	 tEpsp = data(78);
	 tEpsp2 = data(79);
	 tEr = data(80);
	 tEr2 = data(81);
	 refEps = data(82);
	 Beta1 = data(83);
	 Beta2 = data(84);
	 Beta3 = data(85);
	 Alpha1 = data(86);
	 Alpha2 = data(87);
	 Alpha3 = data(88);
	 TotalE = data(89);
	 DeltaE = data(90);
	 Tangent = data(91);
	
	// constaint variables
	 C_E = data(92);
	 C_Dteq = data(93);
	 C_q = data(94);
	 C_beta = data(95);
	 C_c = data(96);
	 C_gamma = data(97);
	 C_CapYieldStressM = data(98);
	 C_CapYieldStrainM = data(99);
	 C_Ed1EM = data(100);
	 C_Ed2EM = data(101);
	 C_DetCapStressM = data(102);
	 C_ay = data(103);
	 C_au = data(104);
	 C_Lambda1 = data(105);
	 C_c1 = data(106);
	 C_Lambda2 = data(107);
	 C_c2 = data(108);
	 C_Lambda3 = data(109);
	 C_c3 = data(110);
	 C_Et1 = data(111);
	 C_Et2 = data(112);
	 C_Et3 = data(113);
	
	// State variables from last converged state
	 C_status = data(114);
	 C_p_teps = data(115);
	 C_p_neps = data(116);
	 C_p_neps_prev = data(117);
	 C_p_teps_prev = data(118);
	 C_cum_p_teps = data(119);
	 C_sgm_0 = data(120);
	 C_alf_d = data(121);
	 C_alf = data(122);
	 C_ytsgm_p = data(123);
	 C_ytsgm_n = data(124);
	 C_yteps_p = data(125);
	 C_yteps_n = data(126);
	 C_teps = data(127);
	 C_neps = data(128);
	 C_tsgm = data(129);
	 C_nsgm = data(130);
	 C_teps_prev = data(131);
	 C_neps_prev = data(132);
	 C_tsgm_prev = data(133);
	 C_nsgm_prev = data(134);
	 C_cEu = data(135);
	 C_cSgmy = data(136);
	 C_cEpsy = data(137);
	 C_cSgmc = data(138);
	 C_cEpsc = data(139);
	 C_cSgmd1 = data(140);
	 C_cEpsd1 = data(141);
	 C_cSgmd2 = data(142);
	 C_cEpsd2 = data(143);
	 C_cSgmb = data(144);
	 C_cSgmd = data(145);
	 C_cEsth = data(146);
	 C_cEd1 = data(147);
	 C_cEd2 = data(148);
	 C_cIniSgmy = data(149);
	 C_cIniEpsy = data(150);
	 C_cIniSgmc = data(151);
	 C_cIniEpsc = data(152);
	 C_cIniEsth = data(153);
	 C_cIniEd1 = data(154);
	 C_cIniEd2 = data(155);
	 C_cIniSgmd1 = data(156);
	 C_cIniEpsd1 = data(157);
	 C_cIniSgmb = data(158);
	 C_cIniSgmd = data(159);
	 C_cIniSgmd2 = data(160);
	 C_cIniEpsd2 = data(161);
	 C_tEu = data(162);
	 C_tSgmy = data(163);
	 C_tEpsy = data(164);
	 C_tSgmp = data(165);
	 C_tEpsp = data(166);
	 C_tEpsp2 = data(167);
	 C_tEr = data(168);
	 C_tEr2 = data(169);
	 C_refEps = data(170);
	 C_Beta1 = data(171);
	 C_Beta2 = data(172);
	 C_Beta3 = data(173);
	 C_Alpha1 = data(174);
	 C_Alpha2 = data(175);
	 C_Alpha3 = data(176);
	 C_TotalE = data(177);
	 C_DeltaE = data(178);
	 C_Tangent = data(179);
  	

	iInitial = data(180);
	C_iInitial = data(181);
  	
  	}

  return res;
}


void 
SLModel::Print(OPS_Stream &s, int flag)
{
	s << "SLModel tag: " << this->getTag() << endln;
	s << "  Dt: " << Dt << endln; 
	s << "  sgm_ini: " << sgm_ini << endln;
	//s << "  OP_Material: " << OP_Material << endln;
}



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/********************************************************************************************
 *** StrainHardeningFunc
 ********************************************************************************************/
void SLModel::StrainHardeningFunc(void)
{
	//define local variables
	double Deltap_teps;
	
	if (status == 1) {
		
	//plastic strain
	p_teps = p_teps_prev;
	cum_p_teps  = cum_p_teps;
	
	//true stress
	tsgm = (teps-p_teps)*E;
	sgm_0 = sgm_0;
	alf_d = 0.0;
	alf = alf+alf_d;
	
	//nominal stress
	nsgm = tsgm/exp(teps);
	
	//nominal plastic strain
	p_neps = p_neps_prev;
	
	//Tangent
	Tangent = E;
	//Tangent = (nsgm-nsgm_prev)/(neps-neps_prev);
		
	} else if (status == 2) {
		
		//true plastic strain
		p_teps = teps-tsgm/E;
		
		//increment of plastic strain 
		Deltap_teps = (p_teps-p_teps_prev)/5.0;
		
		for (int jj=1;jj<=5;jj++) {
			
			//kinematic hardening component 
			alf_d = c/sgm_0*(tsgm-alf)*(Deltap_teps)-gamma*alf*(Deltap_teps);
			alf = alf+alf_d; //plus
			
			//isotropic hardening component 
			cum_p_teps = cum_p_teps+fabs(Deltap_teps);
			sgm_0 = sgm_ini+q*(1-exp(-1*beta*cum_p_teps));
			
			//total true stress
			tsgm = alf+sgm_0; //plus
			
		}
		
		//nominal stress
		 nsgm = tsgm/exp(teps);
		
		//nominal plastic strain
		p_neps = neps-nsgm/E;
		
		//Tangent
		Tangent = (nsgm-nsgm_prev)/(neps-neps_prev);
		
	} else if (status == 3) {
		
		 //true plastic strain
		p_teps = teps-tsgm/E;
		
		//increment of plastic strain 
		Deltap_teps = (p_teps-p_teps_prev)/5.0;
		
		for (int jj=1;jj<=5;jj++) {
			//kinematic hardening component 
			alf_d = c/sgm_0*(tsgm-alf)*(Deltap_teps)-gamma*alf*(Deltap_teps);
			alf = alf-alf_d; //minus
			
			//isotropic hardening component 
			cum_p_teps = cum_p_teps+fabs(Deltap_teps);
			sgm_0 = sgm_ini+q*(1-exp(-1*beta*cum_p_teps));
			
			//total true stress
			tsgm = alf-sgm_0; //minus
		}
		
		//nominal stress
		nsgm = tsgm/exp(teps);
		
		//nominal plastic strain
		p_neps = neps-nsgm/E;
		
		//Tangent
		Tangent = (nsgm-nsgm_prev)/(neps-neps_prev);
		
	}
}



/********************************************************************************************
 *** YieldPointFunc
 ********************************************************************************************/
void SLModel::YieldPointFunc(void)
{
	if (status == 2) {
		ytsgm_p = tsgm;
		yteps_p = teps;
		ytsgm_n = tsgm-2*sgm_0; //negative
		yteps_n = teps-2*sgm_0/E; //negative
	} else if (status == 3) {
		ytsgm_p = tsgm+2*sgm_0; //positive
		yteps_p = teps+2*sgm_0/E;  //positive
		ytsgm_n = tsgm;
		yteps_n = teps;
	}
}


/********************************************************************************************
 *** post-buckling excursion in compression
 ********************************************************************************************/
void SLModel::BackBoneCompFunc(void)
{
	
	// local variables
	double cEpsOffset;
	double cSgme;
	double TempEps1, TempEps2, TempEps3;
	int itemp;
	
	
	// unloading slope
	if (neps < tEpsp) {
		cEu = E*(au/(au+tEpsp-neps));
		if (cEu > E) {
			cEu = E;
		}
	} else {
		cEu = E;
	}
	

	
	// strength deterioration of basic and 1st deterioration strength
	cSgmb = cIniSgmb*Alpha1;
	cSgmd = cIniSgmd*Alpha1;
	
	// intersection of elasitc stiffness and strain hardening slope
	cSgmy = E*cSgmb/(E-cIniEsth);
	cEpsy = cSgmb/(E-cIniEsth);
	
	// stiffness deterioration
	cEsth = cIniEsth*Alpha1;
	
	// modification of basic strength
	cSgmb = cSgmy-cEsth*cEpsy;
	
	// intersection of strain hardening and 1st deterioration slope
	cEpsc = -(cSgmd-cSgmb)/(cIniEd1-cEsth);
	cSgmc = cIniEd1*cEpsc+cSgmd;
	
	// stiffness deterioration
	cEd1 = cIniEd1*Alpha1;
	
	// modification of 1st deterioration strength
	cSgmd = cSgmc-cEd1*cEpsc;
	
	// strength and stiffness deterioraiton of 2nd deterioration slope
	cEd2 = cIniEd2*Alpha1;
	cSgmd2 = cIniSgmd2*Alpha2;
	cEpsd2 = -cSgmd2/cEd2;
	
	// intersection of 1st deterioration slope and 2nd deterioration slope   
	cEpsd1 = -(cSgmd2-cSgmd)/(cEd2-cEd1);
	cSgmd1 = cEd1*cEpsd1+cSgmd;
	
	
	// offset the back bone curve
	if (refEps >= neps-nsgm/cEu) {
		cEpsOffset = refEps;
	} else {
		cEpsOffset = neps-nsgm/cEu;
		refEps = cEpsOffset;
	}
	
	cEpsy = cEpsy+cEpsOffset;
	cEpsc = cEpsc+cEpsOffset;
	cEpsd1 = cEpsd1+cEpsOffset;
	cEpsd2 = cEpsd2+cEpsOffset;
	
	cSgmd2 = cSgmd1-cEd2*cEpsd1;
	cSgmd = cSgmd1-cEd1*cEpsd1;
	cSgmb = cSgmc-cEsth*cEpsc;
	
	cSgme = nsgm-cEu*neps;
	
	// modification of yield stress when 2nd deterioration slope is higher than capping strength
	
	itemp = 0;
	
	if (cEpsc < cEpsd1) {
		cEpsy = -(E*cEpsOffset-cSgmd2)/(E-cEd2);
		cSgmy = E*(cEpsy-cEpsOffset);
		cEpsc = cEpsy;
		cSgmc = cSgmy;
		cEpsd1 = cEpsy;
		cSgmd1 = cSgmy;
		
		itemp = 1;
	}
	
	//
	TempEps1 = neps-(nsgm-cSgmc)/cEu;
	TempEps2 = neps-(nsgm-cSgmd1)/cEu;
	TempEps3 = neps-(nsgm)/cEu;
	
	double ReductionFactor;
	ReductionFactor = 0.0;///////////////////////////////**************************************************************************************************************************
	
	if (cEpsc <= TempEps1) {
		if (itemp == 0) {
			cEpsy = -(cSgme-cSgmb)/(cEu-cEsth);
			cSgmy = cEu*cEpsy+cSgme;
		} else if (itemp == 1) {
			cEpsy = -(cSgme-cSgmd2)/(cEu-cEd2);
			cSgmy = cEu*cEpsy+cSgme;
			
			cEpsc = cEpsy;
			cSgmc = cSgmy;
			cEpsd1 = cEpsy;
			cSgmd1 = cSgmy;
			
			
			cEpsy = cEpsy-cSgmy/cEu*ReductionFactor;
			cSgmy = cEu*cEpsy+cSgme;
			
			cSgmb = cSgmy-cEsth*cEpsy;
			
			cEpsc = -(cSgmb-cSgmd2)/(cEsth-cEd2);
			cSgmc = cEsth*cEpsc+cSgmb;
			
			cEpsd1 = cEpsc;
		}
	} else if (TempEps1 < cEpsc && cEpsd1 <= TempEps2) {
		
		cEpsy = -(cSgme-cSgmd)/(cEu-cEd1);
		cSgmy = cEu*cEpsy+cSgme;
		cEpsy = cEpsy-cSgmy/cEu*ReductionFactor;
		cSgmy = cEu*cEpsy+cSgme;
		
		cSgmb = cSgmy-cEsth*cEpsy;
		
		double TempcEpsc1, TempcEpsc2;
		TempcEpsc1 = -(cSgmb-cSgmd)/(cEsth-cEd1);
		TempcEpsc2 = -(cSgmb-cSgmd2)/(cEsth-cEd2);
		
		if (TempcEpsc1 < TempcEpsc2) {
			cEpsc = TempcEpsc1;
		} else {
			cEpsc = TempcEpsc2;
		}
		cSgmc = cEsth*cEpsc+cSgmb;
		
	} else if (TempEps2 < cEpsd1 && cEpsd2 <= TempEps3) {
		
		cEpsy = -(cSgme-cSgmd2)/(cEu-cEd2);
		cSgmy = cEu*cEpsy+cSgme;
		cEpsy = cEpsy-cSgmy/cEu*ReductionFactor;
		cSgmy = cEu*cEpsy+cSgme;
		
		cSgmb = cSgmy-cEsth*cEpsy;
		
		cEpsc = -(cSgmb-cSgmd2)/(cEsth-cEd2);
		cSgmc = cEsth*cEpsc+cSgmb;
		
		cEpsd1 = cEpsc;
	}
	
	
	if (nsgm > tSgmp) {
		tEpsp = neps;
		tSgmp = nsgm;
	}
	
	
}


/********************************************************************************************
 *** post-buckling excursion in tension
 ********************************************************************************************/
void SLModel::BackBoneTenFunc(void)
{
	// local variables
	double TempSgm, TempEps, TempSgm2, TempEps2, Err; 
	
	//unloading slope
	if (neps < tEpsp) {
		tEu = E*(au/(au+tEpsp-neps));
		if (tEu > E) {
			tEu = E;
		}
	} else {
		tEu = E;
	}
	
	//yield stress
	
	if (neps < tEpsp) {
		TempSgm = sgm_ini*(ay/(ay+tEpsp-neps));
		if (TempSgm > sgm_ini*0.999999) {
			TempSgm = sgm_ini*0.999999;
		}
	} else {
		TempSgm = sgm_ini*0.999999;
	}
	
    TempEps = neps+(TempSgm-nsgm)/tEu;
	
	for (int i=1;i<=20;i++) {
		TempSgm2 = TempSgm;
		TempEps2 = TempEps;
		
		if (TempEps2 < tEpsp) {
			TempSgm = sgm_ini*(ay/(ay+tEpsp-neps));
			if (TempSgm > sgm_ini*0.999999) {
				TempSgm = sgm_ini*0.999999;
			}
		} else {
			TempSgm = sgm_ini*0.999999;
		}
		
		TempEps = TempEps2+(TempSgm-TempSgm2)/tEu;
		
		tSgmy = TempSgm;
		tEpsy = TempEps;
		
		Err = fabs(TempSgm-TempSgm2);
		if (Err < 1.0e-5) {
			break;
		}
	};
	
	
	//reloading slope
	tEr = (tSgmp-tSgmy)/(tEpsp-tEpsy);
	tEr2 = tEr*Alpha3;
	tEpsp2 = tEpsy+(tSgmp-tSgmy)/tEr2;
	
	
}
	
/********************************************************************************************
 *** Unloading in compression
 ********************************************************************************************/
void SLModel::BackBoneComp2Func(void)
{
	
	// local variables
	double cSgme;
	double TempEps1, TempEps2, TempEps3;
	int itemp;
	double cEpsOffset;
	
	cEu = tEu;
	
	
	
	
	// modification of yield stress when 2nd deterioration slope is higher than capping strength
	
	itemp = 0;
	cEpsOffset = refEps;
	
	if (cEpsc < cEpsd1) {
		cEpsy = -(E*cEpsOffset-cSgmd2)/(E-cEd2);
		cSgmy = E*(cEpsy-cEpsOffset);
		cEpsc = cEpsy;
		cSgmc = cSgmy;
		cEpsd1 = cEpsy;
		cSgmd1 = cSgmy;
		
		itemp = 1;
	}
	
	cSgme = nsgm-cEu*neps;
	
	//
	TempEps1 = neps-(nsgm-cSgmc)/cEu;
	TempEps2 = neps-(nsgm-cSgmd1)/cEu;
	TempEps3 = neps-(nsgm)/cEu;
	
	double ReductionFactor;
	ReductionFactor = 0.0;///////////////////////////////**************************************************************************************************************************
	
	if (cEpsc <= TempEps1) {
		if (itemp == 0) {
			cEpsy = -(cSgme-cSgmb)/(cEu-cEsth);
			cSgmy = cEu*cEpsy+cSgme;
		} else if (itemp == 1) {
			cEpsy = -(cSgme-cSgmd2)/(cEu-cEd2);
			cSgmy = cEu*cEpsy+cSgme;
			
			cEpsc = cEpsy;
			cSgmc = cSgmy;
			cEpsd1 = cEpsy;
			cSgmd1 = cSgmy;
			
			
			cEpsy = cEpsy-cSgmy/cEu*ReductionFactor;
			cSgmy = cEu*cEpsy+cSgme;
			
			cSgmb = cSgmy-cEsth*cEpsy;
			
			cEpsc = -(cSgmb-cSgmd2)/(cEsth-cEd2);
			cSgmc = cEsth*cEpsc+cSgmb;
			
			cEpsd1 = cEpsc;
		}
	} else if (TempEps1 < cEpsc && cEpsd1 <= TempEps2) {
		
		cEpsy = -(cSgme-cSgmd)/(cEu-cEd1);
		cSgmy = cEu*cEpsy+cSgme;
		cEpsy = cEpsy-cSgmy/cEu*ReductionFactor;
		cSgmy = cEu*cEpsy+cSgme;
		
		cSgmb = cSgmy-cEsth*cEpsy;
		
		double TempcEpsc1, TempcEpsc2;
		TempcEpsc1 = -(cSgmb-cSgmd)/(cEsth-cEd1);
		TempcEpsc2 = -(cSgmb-cSgmd2)/(cEsth-cEd2);
		
		if (TempcEpsc1 < TempcEpsc2) {
			cEpsc = TempcEpsc1;
		} else {
			cEpsc = TempcEpsc2;
		}
		cSgmc = cEsth*cEpsc+cSgmb;
		
				
	} else if (TempEps2 < cEpsd1 && cEpsd2 <= TempEps3) {
		
		cEpsy = -(cSgme-cSgmd2)/(cEu-cEd2);
		cSgmy = cEu*cEpsy+cSgme;
		cEpsy = cEpsy-cSgmy/cEu*ReductionFactor;
		cSgmy = cEu*cEpsy+cSgme;
		
		cSgmb = cSgmy-cEsth*cEpsy;
		
		cEpsc = -(cSgmb-cSgmd2)/(cEsth-cEd2);
		cSgmc = cEsth*cEpsc+cSgmb;
		
		cEpsd1 = cEpsc;
	}
	
	
	if (nsgm > tSgmp) {
		tEpsp = neps;
		tSgmp = nsgm;
	}
	
	
}
	
	
	
/********************************************************************************************
 *** Unloading in tension
 ********************************************************************************************/
void SLModel::BackBoneTen2Func(void)
{
	// local variables
	double tSgme, TemptSgmp; 
	
	//unloading slope
    tEu = cEu;
    tSgme = nsgm-tEu*neps;
	
	//yield stress
	
	TemptSgmp = tSgmp-tEpsp2*tEr2;
	tEpsy = -(tSgme-TemptSgmp)/(tEu-tEr2);
	tSgmy = tEu*tEpsy+tSgme;
	
	
}



