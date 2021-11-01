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
/* ********************************************************************
** Code developed at UC San Diego
**
** Programmed by: Rodrigo Carreno (rpcarren@ucsd.edu)
** Supervisor: Jose Restrepo (jrestrepo@ucsd.edu)
**
** ********************************************************************/
// Written by: Rodrigo Carreno
// Created	 : February 2017

// Description: Definition of Dodd-Restrepo Steel model including post
// necking behavior and alternative computations of the Bauschinger curve


#include <elementAPI.h>
#include "SteelDRC.h"

#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

/*
#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif
*/


//OPS_Export void *
void *
OPS_SteelDRC()
{
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial *theMaterial = 0;

	// Extract number of input parameters
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 6) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "uniaxialMaterial SteelDRC ";
		opserr << "tag? Es? fy? eu? fu? esh?";
		opserr << "<-P?><-shPoint?><-omegaFac?><-bausch?><-stiffness?>\n";
		return 0;
	}
	int tag;
	numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING invalid tag\n";
		return 0;
	}
	double data[5];
	numdata = 5;
	if (OPS_GetDoubleInput(&numdata, data) < 0) {
		opserr << "WARNING invalid double data\n";
		return 0;
	}
	// Initialize default value for optional variables
	double Psh = 3.0;							// Default exponent for the strain hardening branch of the backbone curve
	double efract = -1;							// fracture strain post uniform strain in tension (efract = -1 means no drop in stress capacity post uniform strain)
	double omegaFac = 1.0;						// Factor to control the area under the curve in the bauschinger effect
	int bauschType = 1;							// Function to use to generate the bauschinger effect (bauschType = 1  for cubic Bezier curve)
	int stiffopt = 0;							// Type of stiffness to be returned by the material.
	double temp[2] = { 0.0,0.0 };				// Temporary array to store the (esh1,fsh1) point if given as input
	double temp1[2] = { 0.0,1.0 };				//Temporary array to store properties of viscous damper if given as input
	double Dfu = 1.0;							// Length of unloading branch in terms of the yield strength fy
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char * type = OPS_GetString();
		if (strcmp(type, "-Psh") == 0 || strcmp(type, "-psh") == 0 || strcmp(type, "-PSh") == 0 || strcmp(type, "-PSH") == 0) {
			numdata = 1;
			if (OPS_GetDoubleInput(&numdata, &Psh) < 0) {
				opserr << "WARNING invalid double data for -Psh\n";
				return 0;
			}
		}
		else if (strcmp(type, "-shPoint") == 0 || strcmp(type, "-SHpoint") == 0 || strcmp(type, "-shpoint") == 0) {
			numdata = 2;
			if (OPS_GetDoubleInput(&numdata, temp) < 0) {
				opserr << "WARNING invalid double input for -shPoint\n";
				return 0;
			}
		}
		else if (strcmp(type, "-omegaFactor") == 0 || strcmp(type,"-omegaFac") == 0 || strcmp(type, "-OmegaFactor") == 0 || strcmp(type, "-omegafactor") == 0 || strcmp(type, "-OmegaFac") == 0 || strcmp(type, "-omegafac") == 0) {
			numdata = 1;
			if (OPS_GetDoubleInput(&numdata, &omegaFac) < 0) {
				opserr << "WARNING invalid double data for -omegaFactor\n";
				return 0;
			}
		}
		else if (strcmp(type, "-fractStrain") == 0 || strcmp(type, "-FractStrain") == 0 || strcmp(type, "-fractstrain") == 0) {
			numdata = 1;
			if (OPS_GetDoubleInput(&numdata, &efract) < 0) {
				opserr << "WARNING invalid double data for -fractStrain\n";
				return 0;
			}
		}
		else if (strcmp(type, "-bausch") == 0 || strcmp(type, "-Bausch") == 0) {
			numdata = 1;
			if (OPS_GetIntInput(&numdata, &bauschType) < 0) {
				opserr << "WARNING invalid int data for -bausch\n";
				return 0;
			}
		}
		else if (strcmp(type, "-stiffOutput")== 0 || strcmp(type, "-StiffOutput") == 0 || strcmp(type, "-stiffoutput") == 0) {
			numdata = 1;
			if (OPS_GetIntInput(&numdata, &stiffopt) < 0) {
				opserr << "WARNING invalid int data for -stiffOutput\n";
				return 0;
			}
		}
		else if (strcmp(type, "-viscousDamper") == 0 || strcmp(type, "-ViscousDamper") == 0 || strcmp(type, "-viscousdamper") == 0) {
			numdata = 2;
			if (OPS_GetDoubleInput(&numdata, temp1) < 0) {
				opserr << "WARNING invalid double data for -ViscousDamper\n";
				return 0;
			}
		}
		else if (strcmp(type, "-Dfu") == 0 || strcmp(type, "-dfu") == 0 || strcmp(type, "-DFu") == 0) {
			numdata = 1;
			if (OPS_GetDoubleInput(&numdata, &Dfu) < 0) {
				opserr << "WARNING invalid double data for -Dfu\n";
				return 0;
			}
		}
		else {
			opserr << "WARNING SteelDRC: invalid material property \n";
			opserr << "Possible Optional Flags : ";
			opserr << "<-Psh?> <-shPoint> <-omegaFac?> <-fractStrain?> <-bausch?> <-stiffOutput?> <-ViscousDamper?><-Dfu?>\n";
			return 0;
		}
	}

	// Use a different constructor if the (esh1, fsh1) is given as input or, alternatively Psh is given as input
	if (temp[0] != 0 && temp[1] != 0) {
		theMaterial = new SteelDRC(tag, data[0], data[1], data[2], data[3], data[4],
			temp[0], temp[1], efract, omegaFac, bauschType, stiffopt,temp1[0],temp1[1],Dfu);
	}
	else {
		theMaterial = new SteelDRC(tag, data[0], data[1], data[2], data[3], data[4],
			Psh, efract, omegaFac, bauschType, stiffopt,temp1[0],temp1[1],Dfu);
	}
	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type SteelDRC\n";
		return 0;
	}

	// return the material
	return theMaterial;
}

// Material Constructors and destructor
// Simplest class constructor, used to generate copy of the object
SteelDRC::SteelDRC(int tag)
	:UniaxialMaterial(tag, MAT_TAG_SteelDRC) {
}
// Constructor for the case when all the model parameters are given and strain hardening is defined from (esh1, fsh1)
SteelDRC::SteelDRC(int tag, double Es, double fy, double eu, double fu, double esh,
	double esh1, double fsh1, double eft, double omegaFac, int bauschType,
	int stiffoption, double C_visc0, double alpha0, double Dfu) 
	: UniaxialMaterial(tag, MAT_TAG_SteelDRC), E(Es), fyEng(fy), eshEng(esh), fuEng(fu), omegaF(omegaFac),
	bauschFlag(bauschType), Etflag(stiffoption), C_visc(C_visc0), alpha(alpha0), Dfu(Dfu)
{
	// Convert yield strain/stress into natural coordinates
	double aux0[3] = { fy / E, fy,0.0 };
	eng2natural(aux0, 2);
	eyN = aux0[0];
	fyN = aux0[1];
	aux0[0] = eu;
	aux0[1] = fu;
	eng2natural(aux0, 2);
	euN = aux0[0];
	fuN = aux0[1];
	aux0[0] = esh;
	eng2natural(aux0);
	eshN = aux0[0];
	aux0[0] = esh1;
	aux0[1] = fsh1;
	eng2natural(aux0, 2);
	double fshN = fy*exp(eshN);// stress at onset of strain hardening
	double c1 = fshN + fuN*(euN - eshN) - fuN;
	Psh = log((aux0[1] + fuN*(euN - aux0[0]) - fuN) / c1) / log((euN - aux0[0]) / (euN - eshN));
	// 2021-04-18 : Added case where the input for efract is -1 to indicate a horizontal line following uniform strain in tension
	if (eft == -1) {
		eftN = -1;
	}
	else {
		aux0[0] = eft;
		eng2natural(aux0);
		eftN = aux0[0];
	}
	
	// Initialize the state variables
	this->revertToStart();
}

// Constructor for the case when all the model parameters are given and strain hardening is defined from input Psh
SteelDRC::SteelDRC(int tag, double Es, double fy, double eu, double fu, double esh,
	double Psh0, double eft, double omegaFac, int bauschType, int stiffoption, double C_visc0, double alpha0, double Dfu)
	: UniaxialMaterial(tag, MAT_TAG_SteelDRC), E(Es), fyEng(fy), eshEng(esh), fuEng(fu), Psh(Psh0) ,omegaF(omegaFac),
	bauschFlag(bauschType), Etflag(stiffoption), C_visc(C_visc0), alpha(alpha0), Dfu(Dfu)
{
	// Convert yield strain/stress into natural coordinates
	double aux0[3] = { fy / E, fy, 0.0 };
	eng2natural(aux0, 2);
	eyN = aux0[0];
	fyN = aux0[1];
	aux0[0] = eu;
	aux0[1] = fu;
	eng2natural(aux0, 2);
	euN = aux0[0];
	fuN = aux0[1];
	aux0[0] = esh;
	eng2natural(aux0);
	eshN = aux0[0];
	aux0[0] = eft;
	if (eft == -1) {
		eftN = -1;
	}
	else {
		eng2natural(aux0);
		eftN = aux0[0];
	}
	// Initialize the state variables
	this->revertToStart();
}

SteelDRC::~SteelDRC()
{
	// does nothing
}

int
SteelDRC::setTrialStrain(double strain, double strainRate)
{
	// Copy state variables of last commited state into trial state.
	revertToLastCommit();

	// If strain change is very small no update is necessary
	if (fabs(trialStrain - strain) < DBL_EPSILON)
		return 0;

	trialStrain = strain;
	trialStrainRate = strainRate; 

	// Transform input strain into natural coordinates
	double aux[3] = { strain, 0.0, 0.0 };
	eng2natural(aux);
	Teps = aux[0];

	// Define the straining direction and index values for two element state arrays
	int S,K,M;
	if ((Teps - Ceps) > 0)
	{
		S = M = 1;
		K = 0;
	}
	else	{
		S = -1;
		K = 1;
		M = 0;
	}
	// Define index Klmr for straining direction of the last major reversal
	int Klmr;
	switch (Tlmr) {
		case 1:	Klmr = 0;
				break;
		case -1:Klmr = 1;
				break; 
		default: Klmr = -1;
				 break;
	}
	// Compute current unloading slope
	double Eun = E*(0.82 + 1 / (5.55 + 1000 * fabs(Te0max)));
	// Verify if previous strain (Ceps) is a reversal point
	if ((Ceps - Ter)*(Teps - Ceps) < 0)
		State_Reversal(S,K,M,Klmr,Eun);
	// Compute the stress and slope corresponding to the new strain

	State_Determination(S, K, M, Klmr, Eun);
	// Transform resulting stress and slope in natural coordinates to engineering values
	aux[1] = Tsig;
	aux[2] = Ttan;
	natural2eng(aux, 3);
	// Add viscous damper contribution to the stress
	double absRate = fabs(trialStrainRate);
	double signRate;
	if (trialStrainRate >= 0)
		signRate = 1;
	else
		signRate = -1;

	if (absRate > 1E-10) 
		trialStress = aux[1] + signRate*C_visc*pow(absRate, alpha);
	else
		trialStress = aux[1];

	trialTangent = aux[2];
	return 0;
}

double SteelDRC::getStrain(void)
{
	return trialStrain;
}

double SteelDRC::getStrainRate(void)
{
	return trialStrainRate;
}

double SteelDRC::getStress(void)
{
	return trialStress;
}
double SteelDRC::getInitialTangent(void) 
{ 
	return E; 
}


double SteelDRC::getTangent(void)
{
	switch (Etflag) {
		case 0: // tangent stifness
			// The function cannot return a slope less than 05% of the initial slope
			return fmin(fmax(0.005*E,trialTangent),E);
			break;
		case 1: // secant slope between last reversal point and current strain/stress.
			if (fabs(trialStrain - Ter) <= DBL_EPSILON)
				// The function cannot return a slope less than 0.5% of the initial slope
				return fmin(fmax(0.005*E,trialTangent),E);
			else 
				return fmin(fmax(fmax(trialTangent, (trialStress - Tsr) / (trialStrain - Ter)), 0.005*E), E);
			break;
		case 2:
			return E;
			break;
		default:
			opserr << "WARNING SteelDRC invalid output stiffness option, tangent stiffness used";
			return trialTangent;
			break;
	}
	//return trialTangent;
}

int
SteelDRC::commitState(void)
{
	// Commit strain, stress and tangent in natural coordinates
	Ceps = Teps;
	Csig = Tsig;
	Ctan = Ttan;

	// Commit remaining state variables
	Clmr = Tlmr;
	Ce0[0] = Te0[0];
	Ce0[1] = Te0[1];
	Ce0max = Te0max;
	Cer = Ter;
	Csr = Tsr;
	CEr = TEr;
	Cea[0] = Tea[0];
	Cea[1] = Tea[1];
	Csa[0] = Tsa[0];
	Csa[0] = Tsa[0];
	Cerejoin[0] = Terejoin[0];
	Cerejoin[1] = Terejoin[1];
	Csrejoin[0] = Tsrejoin[0];
	Csrejoin[1] = Tsrejoin[1];
	CErejoin[0] = TErejoin[0];
	CErejoin[1] = TErejoin[1];

	CerejoinL[0] = TerejoinL[0];
	CerejoinL[1] = TerejoinL[1];
	CsrejoinL[0] = TsrejoinL[0];
	CsrejoinL[1] = TsrejoinL[1];
	CErejoinL[0] = TErejoinL[0];
	CErejoinL[1] = TErejoinL[1];

	Cerm[0] = Term[0];
	Cerm[1] = Term[1];
	Csrm[0] = Tsrm[0];
	Csrm[1] = Tsrm[1];
	CErm[0] = TErm[0];
	CErm[1] = TErm[1];
	Cfrm[0] = Tfrm[0];
	Cfrm[1] = Tfrm[1];
	Ceam[0] = Team[0];
	Ceam[1] = Team[1];
	Csam[0] = Tsam[0];
	Csam[1] = Tsam[1];
	Cfract = Tfract;
	CshOnset = TshOnset;

	// Commit Strain, Stress and Tangent in Engineering Coordinates.
	commitStrain = trialStrain;
	commitStrainRate = trialStrainRate;
	commitTangent = trialTangent;
	commitStress = trialStress;

	return 0;
}


int
SteelDRC::revertToLastCommit(void)
{
	trialStrain = commitStrain;
	trialStrainRate = commitStrainRate;
	trialTangent = commitTangent;
	trialStress = commitStress;

	// revert the internal state variables to the previously commited state
	Teps = Ceps; 
	Tsig = Csig;
	Ttan = Ctan;

	Tlmr = Clmr;
	Te0[0] = Ce0[0];
	Te0[1] = Ce0[1];
	Te0max = Ce0max; 
	Ter = Cer;
	Tsr = Csr;    
	TEr = CEr;
	Tea[0] = Cea[0];
	Tea[1] = Cea[1];
	Tsa[0] = Csa[0];
	Tsa[0] = Csa[0];
	Terejoin[0] = Cerejoin[0]; 
	Terejoin[1] = Cerejoin[1];
	Tsrejoin[0] = Csrejoin[0];
	Tsrejoin[1] = Csrejoin[1];
	TErejoin[0] = CErejoin[0];
	TErejoin[1] = CErejoin[1];

	TerejoinL[0] = CerejoinL[0];
	TerejoinL[1] = CerejoinL[1];
	TsrejoinL[0] = CsrejoinL[0];
	TsrejoinL[1] = CsrejoinL[1];
	TErejoinL[0] = CErejoinL[0];
	TErejoinL[1] = CErejoinL[1];

	Term[0] = Cerm[0];
	Term[1] = Cerm[1];
	Tsrm[0] = Csrm[0];
	Tsrm[1] = Csrm[1];
	TErm[0] = CErm[0];
	TErm[1] = CErm[1];
	Tfrm[0] = Cfrm[0];
	Tfrm[1] = Cfrm[1];
	Team[0] = Ceam[0];
	Team[1] = Ceam[1];
	Tsam[0] = Csam[0];
	Tsam[1] = Csam[1];
	Tfract = Cfract;
	TshOnset = CshOnset;

	return 0;
}

int
SteelDRC::revertToStart(void)
{
	trialStrain = commitStrain = 0.0;
	trialStrainRate = commitStrainRate = 0.0;
	trialTangent = commitTangent = E;
	trialStress = commitStress = 0.0;

	Teps = Ceps=0;
	Tsig = Csig=0;
	Ttan = Ctan=E;
	Tlmr = Clmr=0;
	Te0[0] = Ce0[0]=0;
	Te0[1] = Ce0[1]=0;

	Te0max = Ce0max=0;
	Ter = Cer=0;
	Tsr = Csr=0;
	TEr = CEr=0;
	Tea[0] = Cea[0]=0;
	Tea[1] = Cea[1]=0;
	Tsa[0] = Csa[0]=0;
	Tsa[1] = Csa[1]=0;
	Terejoin[0] = Cerejoin[0]=0;
	Terejoin[1] = Cerejoin[1]=0;
	Tsrejoin[0] = Csrejoin[0]=0;
	Tsrejoin[1] = Csrejoin[1]=0;
	TErejoin[0] = CErejoin[0]=0;
	TErejoin[1] = CErejoin[1]=0;

	TerejoinL[0] = CerejoinL[0] = euN;
	TerejoinL[1] = CerejoinL[1] = -euN;
	TsrejoinL[0] = CsrejoinL[0] = fuN;
	TsrejoinL[1] = CsrejoinL[1] = -fuN;
	TErejoinL[0] = CErejoinL[0] = fuN;
	TErejoinL[1] = CErejoinL[1] = fuN;

	Term[0] = Cerm[0]=0;
	Term[1] = Cerm[1]=0;
	Tsrm[0] = Csrm[0]=0;
	Tsrm[1] = Csrm[1]=0;
	TErm[0] = CErm[0]=0;
	TErm[1] = CErm[1]=0;
	Tfrm[0] = Cfrm[0]=0;
	Tfrm[1] = Cfrm[1]=0;
	Team[0] = Ceam[0]=0;
	Team[1] = Ceam[1]=0;
	Tsam[0] = Csam[0]=0;
	Tsam[1] = Csam[1]=0;
	Tfract = Cfract=0;
	TshOnset = CshOnset = 0;

	return 0;
}


UniaxialMaterial *
SteelDRC::getCopy(void)
{
	SteelDRC *theCopy =  new SteelDRC(this->getTag());

	// Save material properties into the copy of the material object
	theCopy->fyEng = fyEng;
	theCopy->fuEng = fuEng;
	theCopy->E = E;
	theCopy->omegaF = omegaF;
	theCopy->bauschFlag = bauschFlag;
	theCopy->Etflag = Etflag;
	theCopy->eyN = eyN;
	theCopy->fyN = fyN;
	theCopy->euN = euN;
	theCopy->fuN = fuN;
	theCopy->eshN = eshN;
	theCopy->Psh = Psh;
	theCopy->eftN = eftN;
	theCopy->eshEng = eshEng;
	theCopy->C_visc = C_visc;
	theCopy->alpha = alpha;
	theCopy->Dfu = Dfu;

	// Copy state variables in the trial state into the new object
	theCopy->Teps = Teps;
	theCopy->Tsig = Tsig;
	theCopy->Ttan = Ttan;
	theCopy->Tlmr = Tlmr;
	theCopy->Te0[0] = Te0[0];
	theCopy->Te0[1] = Te0[1];
	theCopy->Te0max = Te0max;
	theCopy->Ter = Ter;
	theCopy->Tsr = Tsr;
	theCopy->TEr = TEr;
	theCopy->Tea[0] = Tea[0];
	theCopy->Tea[1] = Tea[1];
	theCopy->Tsa[0] = Tsa[0];
	theCopy->Tsa[1] = Tsa[1];
	theCopy->Terejoin[0] = Terejoin[0];
	theCopy->Terejoin[1] = Terejoin[1];
	theCopy->Tsrejoin[0] = Tsrejoin[0];
	theCopy->Tsrejoin[1] = Tsrejoin[1];
	theCopy->TErejoin[0] = TErejoin[0];
	theCopy->TErejoin[1] = TErejoin[1];

	theCopy->TerejoinL[0] = TerejoinL[0];
	theCopy->TerejoinL[1] = TerejoinL[1];
	theCopy->TsrejoinL[0] = TsrejoinL[0];
	theCopy->TsrejoinL[1] = TsrejoinL[1];
	theCopy->TErejoinL[0] = TErejoinL[0];
	theCopy->TErejoinL[1] = TErejoinL[1];

	theCopy->Term[0] = Term[0];
	theCopy->Term[1] = Term[1];
	theCopy->Tsrm[0] = Tsrm[0];
	theCopy->Tsrm[1] = Tsrm[1];
	theCopy->TErm[0] = TErm[0];
	theCopy->TErm[1] = TErm[1];
	theCopy->Tfrm[0] = Tfrm[0];
	theCopy->Tfrm[1] = Tfrm[1];
	theCopy->Team[0] = Team[0];
	theCopy->Team[1] = Team[1];
	theCopy->Tsam[0] = Tsam[0];
	theCopy->Tsam[1] = Tsam[1];
	theCopy->Tfract = Tfract;
	theCopy->TshOnset = TshOnset;

	theCopy->trialStrain = trialStrain;
	theCopy->trialStrainRate = trialStrainRate;
	theCopy->trialStress = trialStress;
	theCopy->trialTangent = trialTangent;

	// Copy commited state variables into the new object
	theCopy->Ceps = Ceps;
	theCopy->Csig = Csig;
	theCopy->Ctan = Ctan;
	theCopy->Clmr = Clmr;
	theCopy->Ce0[0] = Ce0[0];
	theCopy->Ce0[1] = Ce0[1];
	theCopy->Ce0max = Ce0max;
	theCopy->Cer = Cer;
	theCopy->Csr = Csr;
	theCopy->CEr = CEr;
	theCopy->Cea[0] = Cea[0];
	theCopy->Cea[1] = Cea[1];
	theCopy->Csa[0] = Csa[0];
	theCopy->Csa[1] = Csa[1];
	theCopy->Cerejoin[0] = Cerejoin[0];
	theCopy->Cerejoin[1] = Cerejoin[1];
	theCopy->Csrejoin[0] = Csrejoin[0];
	theCopy->Csrejoin[1] = Csrejoin[1];
	theCopy->CErejoin[0] = CErejoin[0];
	theCopy->CErejoin[1] = CErejoin[1];

	theCopy->CerejoinL[0] = CerejoinL[0];
	theCopy->CerejoinL[1] = CerejoinL[1];
	theCopy->CsrejoinL[0] = CsrejoinL[0];
	theCopy->CsrejoinL[1] = CsrejoinL[1];
	theCopy->CErejoinL[0] = CErejoinL[0];
	theCopy->CErejoinL[1] = CErejoinL[1];

	theCopy->Cerm[0] = Cerm[0];
	theCopy->Cerm[1] = Cerm[1];
	theCopy->Csrm[0] = Csrm[0];
	theCopy->Csrm[1] = Csrm[1];
	theCopy->CErm[0] = CErm[0];
	theCopy->CErm[1] = CErm[1];
	theCopy->Cfrm[0] = Cfrm[0];
	theCopy->Cfrm[1] = Cfrm[1];
	theCopy->Ceam[0] = Ceam[0];
	theCopy->Ceam[1] = Ceam[1];
	theCopy->Csam[0] = Csam[0];
	theCopy->Csam[1] = Csam[1];
	theCopy->Cfract = Cfract;
	theCopy->CshOnset = CshOnset;

	theCopy->commitStrain = commitStrain;
	theCopy->commitStrainRate = commitStrainRate;
	theCopy->commitStress = commitStress;
	theCopy->commitTangent = commitTangent;

	return theCopy;
}


int
SteelDRC::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	int index = 0;
	static Vector data(70);
	
	data(index++) = this->getTag();
	data(index++) = Teps;
	data(index++) = Ceps;
	data(index++) = Tsig;
	data(index++) = Csig;
	data(index++) = Ttan;
	data(index++) = Ctan;
	data(index++) = Tlmr;
	data(index++) = Clmr;
	data(index++) = Te0[0];
	data(index++) = Ce0[0];
	data(index++) = Te0[1];
	data(index++) = Ce0[1];
	data(index++) = Te0max;
	data(index++) = Ce0max;
	data(index++) = Ter;
	data(index++) = Cer;
	data(index++) = Tsr;
	data(index++) = Csr;
	data(index++) = TEr;
	data(index++) = CEr;
	data(index++) = Tea[0];
	data(index++) = Cea[0];
	data(index++) = Tea[1];
	data(index++) = Cea[1];
	data(index++) = Tsa[0];
	data(index++) = Csa[0];
	data(index++) = Tsa[1];
	data(index++) = Csa[1];
	data(index++) = Terejoin[0];
	data(index++) = Cerejoin[0];
	data(index++) = Terejoin[1];
	data(index++) = Cerejoin[1];
	data(index++) = Tsrejoin[0];
	data(index++) = Csrejoin[0];
	data(index++) = Tsrejoin[1];
	data(index++) = Csrejoin[1];
	data(index++) = TErejoin[0];
	data(index++) = CErejoin[0];
	data(index++) = TErejoin[1];
	data(index++) = CErejoin[1];

	data(index++) = TerejoinL[0];
	data(index++) = CerejoinL[0];
	data(index++) = TerejoinL[1];
	data(index++) = CerejoinL[1];
	data(index++) = TsrejoinL[0];
	data(index++) = CsrejoinL[0];
	data(index++) = TsrejoinL[1];
	data(index++) = CsrejoinL[1];
	data(index++) = TErejoinL[0];
	data(index++) = CErejoinL[0];
	data(index++) = TErejoinL[1];
	data(index++) = CErejoinL[1];

	data(index++) = Term[0];
	data(index++) = Cerm[0];
	data(index++) = Term[1];
	data(index++) = Cerm[1];
	data(index++) = Tsrm[0];
	data(index++) = Csrm[0];
	data(index++) = Tsrm[1];
	data(index++) = Csrm[1];
	data(index++) = TErm[0];
	data(index++) = CErm[0];
	data(index++) = TErm[1];
	data(index++) = CErm[1];
	data(index++) = Tfrm[0];
	data(index++) = Cfrm[0];
	data(index++) = Tfrm[1];
	data(index++) = Cfrm[1];
	data(index++) = Team[0];
	data(index++) = Ceam[0];
	data(index++) = Team[1];
	data(index++) = Ceam[1];
	data(index++) = Tsam[0];
	data(index++) = Csam[0];
	data(index++) = Tsam[1];
	data(index++) = Csam[1];
	data(index++) = Tfract;
	data(index++) = Cfract;
	data(index++) = TshOnset;
	data(index++) = CshOnset;
	data(index++) = commitStrain;
	data(index++) = commitStrainRate;
	data(index++) = commitStress;
	data(index++) = commitTangent;

	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "SteelDRC::sendSelf() - failed to send data\n";

	return res;
}

int
SteelDRC::recvSelf(int cTag, Channel &theChannel,
FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector data(70);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "SteelDRC::recvSelf() - failed to recv data\n";
	else {
		int index = 0;
		this->setTag((int)data(index++));

		Teps = data(index++);
		Ceps = data(index++);
		Tsig = data(index++);
		Csig = data(index++);
		Ttan = data(index++);
		Ctan = data(index++);
		Tlmr = (int)data(index++);
		Clmr = (int)data(index++);
		Te0[0] = data(index++);
		Ce0[0] = data(index++);
		Te0[1] = data(index++);
		Ce0[1] = data(index++);
		Te0max = data(index++);
		Ce0max = data(index++);
		Ter = data(index++);
		Cer = data(index++);
		Tsr = data(index++);
		Csr = data(index++);
		TEr = data(index++);
		CEr = data(index++);
		Tea[0] = data(index++);
		Cea[0] = data(index++);
		Tea[1] = data(index++);
		Cea[1] = data(index++);
		Tsa[0] = data(index++);
		Csa[0] = data(index++);
		Tsa[1] = data(index++);
		Csa[1] = data(index++);
		Terejoin[0] = data(index++);
		Cerejoin[0] = data(index++);
		Terejoin[1] = data(index++);
		Cerejoin[1] = data(index++);
		Tsrejoin[0] = data(index++);
		Csrejoin[0] = data(index++);
		Tsrejoin[1] = data(index++);
		Csrejoin[1] = data(index++);
		TErejoin[0] = data(index++);
		CErejoin[0] = data(index++);
		TErejoin[1] = data(index++);
		CErejoin[1] = data(index++);

		TerejoinL[0] = data(index++);
		CerejoinL[0] = data(index++);
		TerejoinL[1] = data(index++);
		CerejoinL[1] = data(index++);
		TsrejoinL[0] = data(index++);
		CsrejoinL[0] = data(index++);
		TsrejoinL[1] = data(index++);
		CsrejoinL[1] = data(index++);
		TErejoinL[0] = data(index++);
		CErejoinL[0] = data(index++);
		TErejoinL[1] = data(index++);
		CErejoinL[1] = data(index++);

		Term[0] = data(index++);
		Cerm[0] = data(index++);
		Term[1] = data(index++);
		Cerm[1] = data(index++);
		Tsrm[0] = data(index++);
		Csrm[0] = data(index++);
		Tsrm[1] = data(index++);
		Csrm[1] = data(index++);
		TErm[0] = data(index++);
		CErm[0] = data(index++);
		TErm[1] = data(index++);
		CErm[1] = data(index++);
		Tfrm[0] = (int)data(index++);
		Cfrm[0] = (int)data(index++);
		Tfrm[1] = (int)data(index++);
		Cfrm[1] = (int)data(index++);
		Team[0] = data(index++);
		Ceam[0] = data(index++);
		Team[1] = data(index++);
		Ceam[1] = data(index++);
		Tsam[0] = data(index++);
		Csam[0] = data(index++);
		Tsam[1] = data(index++);
		Csam[1] = data(index++);
		Tfract = (int)data(index++);
		Cfract = (int)data(index++);
		TshOnset = (int)data(index++);
		CshOnset = (int)data(index++);
		commitStrain = data(index++);
		commitStrainRate = data(index++);
		commitStress = data(index++);
		commitTangent = data(index++);

		trialStrain = commitStrain;
		trialStrainRate = commitStrainRate;
		trialTangent = commitTangent;
		trialStress = commitStress;
	}

	return res;
}

void
SteelDRC::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
 		double aux0[3] = {euN, fuN, fuN};
 		natural2eng(aux0,1);
 		s << "\t\t\t{";
 		s << "\"name\": \"" << this->getTag() << "\", ";
 		s << "\"type\": \"SteelDRC\", ";
 		s << "\"E\": " << E << ", ";
 		s << "\"fy\": " << fyEng << ", ";
 		s << "\"eu\": " << aux0[0] << ", ";
 		s << "\"fu\": " << fuEng << ", ";
 		s << "\"esh\": " << eshEng << ", ";
 		s << "\"P\": " << Psh << ", ";
 		s << "\"omega\": " << omegaF << ", ";
 		s << "\"bausch\": " << bauschFlag << ", ";    
 		s << "\"stiffness\": " << Etflag << ",";
		s << "\"C_visc\": " << C_visc << ",";
		s << "\"alpha\": " << alpha << ",";
		s << "\"Dfu\": " << Dfu << "}";
 	} else {
 		s << "SteelDRC tag: " << this->getTag() << endln;
 		s << "  stress: " << trialStress << " tangent: " << trialTangent << endln;
 	}
}


void SteelDRC::eng2natural(double *ptr1, int n) {
	// ptr1 points to the first element of an array with n elements
	if (ptr1[0] <= -1) // if strain is less than -1 the true strain is an imaginary number
	{
		ptr1[0] = -DBL_MAX;
		return;
	}
	double aux;
	// Transform slope to natural/true coordinates
	if (n==3)
	{
		aux = ptr1[2] * pow(ptr1[0] + 1, 2.0) + ptr1[1] * (ptr1[0] + 1);
		ptr1[2] = aux;
	}
	// Transform stress to natural/true coordinates
	if (n > 1)
	{
		aux = ptr1[1] * (1 + ptr1[0]);
		ptr1[1] = aux;
	}
	// Transform strain to natural/true coordinates
	aux = log(1+ptr1[0]);
	ptr1[0] = aux;
}

void SteelDRC::natural2eng(double *ptr1, int n){
	double aux;
	// Compute slope in engineering coordinates
	if (n == 3)
	{
		aux = (ptr1[2] - ptr1[1])*exp(-2 * ptr1[0]);
		ptr1[2] = aux;
	}
	// Compute stress in engineering coordinates
	if (n > 1)
	{
		aux = ptr1[1] / exp(ptr1[0]);
		ptr1[1] = aux;
	}
	// Compute strain in engineering coordinates
	aux = exp(ptr1[0]) - 1;
	ptr1[0] = aux;
}

void SteelDRC::skeleton(double eps_N,double &sig_N, double &tan_N){
	// compute straining direction
	double sign0;
	if (eps_N == 0)
		sign0 = 0;
	else
		sign0 = eps_N / fabs(eps_N);
	
	// Elastic branch
	if (sign0*eps_N <= eyN) {
		sig_N = E*exp(eps_N)*(exp(eps_N) - 1);
		tan_N = E*exp(eps_N)*(2 * exp(eps_N) - 1);
	}
	// Yield Plateau
	else if (sign0*eps_N < eshN) {
		sig_N = sign0*fyEng*exp(sign0*eps_N);
		tan_N = fyEng*exp(sign0*eps_N);
	}
	// Strain hardening branch
	else if (sign0*eps_N <= euN) {
		TshOnset = 1;
		double c1 = fyEng*exp(eshN) + fuN*(euN - eshN) - fuN;
		sig_N = sign0*c1*pow((euN - sign0*eps_N) / (euN - eshN), Psh) - sign0*fuN*(euN - sign0*eps_N) + sign0*fuN;
		tan_N = -Psh*c1 / (euN - eshN)*pow((euN - sign0*eps_N) / (euN - eshN), Psh - 1) + fuN;
	}
	// Past uniform strain in tension (localization point)
	else if (eps_N > euN) {
		TshOnset = 1;
		Tfract = 1;
		// if Fracture strain has been set up as -1 by the user, material keeps constant stress capacity after necking
		if (eftN == -1) {
			sig_N = fuEng*exp(eps_N);
			tan_N = fuEng*exp(eps_N);
		}
		else {
			// quadratic curve factor for post-necking response
			double a = -fuN*(1 + eftN - euN) / pow(eftN - euN, 2.0);
			// Stress
			//sig_N = fmax(a*pow(eps_N, 2.0) + fuN*eps_N + sign0*fuN, 0);
			sig_N = fmax(a * (eps_N - eftN) * (eps_N + eftN - 2 * euN) - fuN * (eftN - euN), 0);
			// Tangent modulus
			//tan_N = 2 * a*eps_N + fuN;
			tan_N = 2 * a * (eps_N - euN) + fuN;
			if (sig_N == 0) {
				Tfract = -1; // Flag indicating the material has fractured in tension
				tan_N = 0;
			}
		}
	}
	// Past uniform strain in compression
	else {
		TshOnset = 1;
		Tfract = 2;
		sig_N = fuN*(eps_N + euN - 1.0);
		tan_N = fuN;
	}
}

double SteelDRC::omegaFun(double eam, double sam, double e0, int S, int K) {
	// compute ultimate strain in each direction including shifts by plastic strains in Te0.
	double eush[2] = { euN + Te0[0], -euN + Te0[1] };
	// Compute internal variables sigP and sigT (see Restrepo et. al. 1994)
	double sigP = fuN*(S - eush[K] + eam) - sam;
	double sigT = fuN*(2 - eush[0] + eush[1]);
	double epsp = fabs(0.2*S + e0 - eam) / 0.2;
	double sigN = fabs(sigP / sigT);

	// Compute percent area term Omega
	double omega = (0.001 + 0.00108 / (1.043 - epsp)) / 0.18*(sigN - 0.69) + 0.085;
	// Limit omega to the range : [ 0.05-0.085 ]
	omega = fmin(fmax(omega, 0.05), 0.085);

	//Amplify omega by the factor omegaF
	omega *= omegaF;
	return omega;
}

double SteelDRC::PowerP(double eam,double sam,double e0,int S, int K) {
	double omega = omegaFun(eam, sam, e0, S, K);
	//Compute exponent for the Bauschinger curve
	double Power = 56.689 * pow(omega - 0.077,2) - 4.921*(omega - 0.077) + 0.1;
	return Power;
}
double SteelDRC::PowerP(double omega) {
	//Compute exponent for the Bauschinger curve
	double Power = 56.689 * pow(omega - 0.077, 2) - 4.921*(omega - 0.077) + 0.1;
	return Power;
}


void SteelDRC::bauschMinor(int flag, double *ptA, double *ptB, double omega,double eps_N,double &sig_N,double &tan_N) {
	double Pwr;
	double wcubic;
	double xi[2] = { 0.9,0.9 };
	double wcubicarray[4] = { 1.0, 1.3, 0.7, 0.0 };
	switch (flag) {
	case 0:
		Pwr = PowerP(omega);
		bausch1(eps_N, sig_N, tan_N, ptA, ptB, Pwr);
		break;
	case 1:
		wcubic = bezierWeightCubic(omega);
		wcubicarray[3] = wcubic;
		bauschBezierCubic(eps_N, sig_N, tan_N, ptA, ptB, xi, wcubicarray);
		break;
	case 2:
		bauschNURBS(eps_N, sig_N, tan_N, ptA, ptB, 0);
		break;
	}
}
void SteelDRC::bauschMajor(int flag, double *ptA, double *ptB, double eam, double sam, double e0, int S, int K,double eps_N,double &sig_N,double &tan_N) {
	double Pwr;
	double wcubic;
	double omega;
	double xi[2] = { 0.9, 0.9 };
	double wcubicarray[4] = { 1.0, 1.3, 0.7, 0.0 };
	double b;
	switch (flag) {
	case 0:
		Pwr = PowerP(eam, sam, e0, S, K); // Compute exponent of normalized Bauschinger curve
		bausch1(eps_N, sig_N, tan_N, ptA, ptB, Pwr);
		break;
	case 1:
		omega = omegaFun(eam, sam, e0, S, K);
		wcubic = bezierWeightCubic(omega);
		wcubicarray[3] = wcubic;
		bauschBezierCubic(eps_N, sig_N, tan_N, ptA, ptB, xi, wcubicarray);
		break;
	case 2:
		omega = omegaFun(eam, sam, e0, S, K);	
		b = factorb(omega);
		bauschNURBS(eps_N, sig_N, tan_N, ptA, ptB, b);
	}
}

void SteelDRC::bausch1(double eps_N, double &sig_N,double &tan_N, double *pointA, double *pointB, double Pwr) {
	// Define properties of Newton Raphson algorithm run over strains
	double itNR1 = 20;
	double tolNR1 = 10 * DBL_EPSILON;
	// Define properties of Newton Raphson algorithm run over stresses
	double itNR2 = 20;
	double tolNR2 = 10 * DBL_EPSILON;

	// Define normalized strain limit to iterate using strains or stresses
	double epsNlim = 0.2;

	// Internal constant C1
	double C1 = (pointB[1] - pointA[1] - pointA[2] * (pointB[0] - pointA[0])) / 
				(pointB[1] - pointA[1] - pointB[2] * (pointB[0] - pointA[0]));
	// If input strain is outside the limit strains in pointA and pointB, return NAN for stress and slope
	if ((eps_N > fmax(pointA[0], pointB[0])) || (eps_N < fmin(pointA[0], pointB[0])))
	{
		sig_N = NAN;
		tan_N = NAN;
	}
	// If input strain matches the limit points assign stress and slope accordingly
	else if (fabs(eps_N - pointA[0]) < DBL_EPSILON){
		sig_N = pointA[1];
		tan_N = pointA[2];
	}
	else if (fabs(eps_N - pointB[0]) < DBL_EPSILON) {
		sig_N = pointB[1];
		tan_N = pointB[2];
	}
	else {
		// Initialize internal variable C2
		double C2 = ((pointA[2] - pointB[2])*(eps_N - pointA[0])) / 
			(pointB[1] - pointA[1] - pointB[2] * (pointB[0] - pointA[0]));
		// Initialize normalized strain variables to be used by the Newton-Raphson method
		double epsNold, epsNnew, epsNorm;
		epsNold = epsNnew = epsNorm = (eps_N - pointA[0]) / (pointB[0] - pointA[0]);
		// Start Newton Raphson iteration counter
		int j = 1;
		// Compute convergence indicator FepsN for the initial esimate epsNold
		double FepsN = pow(1 - pow(1 - epsNorm, 2), Pwr) - C1*epsNorm - C2;

		while (fabs(FepsN) > tolNR1 && j <= itNR1){
			// New estimate of normalized strain using NR algorithm
			epsNnew = epsNold - (pow(1 - pow(1 - epsNold, 2), Pwr) - C1*epsNold - C2) /
				(2 * Pwr*pow(1 - pow(1 - epsNold, 2), Pwr - 1)*(1 - epsNold) - C1);
			// If normalized strain is too small to use NR algorithm over the strains break the loop
			if (epsNnew < epsNlim)
				break;
			// IF normalized strain is greater than the upper limit (1.0), reset epsNnew 
			if (epsNnew > 1.0)
				epsNnew = epsNlim;
			// Update parameters for next iterations
			epsNold = epsNorm = epsNnew;
			// Update convergence estimate
			FepsN = pow(1 - pow(1 - epsNorm, 2), Pwr) - C1*epsNorm - C2;
			j++;
		}
		
		// Check if previous loop was terminated due to a small normalized strain
		if (epsNnew < epsNlim) {
			//Initialize normalized stress variables to be used by the Newton-Raphson method.
			double sigNold, sigNorm, sigNnew, FsigN;
			j = 1;
			sigNold = 0.5*(C2 / (1 - C1) + 1);
			sigNorm = sigNold;

			// Compute auxiliary variable
			double auxroot = sqrt(1 - pow(sigNorm, 1 / Pwr));
			// Initial convergence measure for the estimate of sigN
			FsigN = sigNorm - C1*(1 - auxroot) - C2;

			// Loop over value of normalized stress until convergence
			while (fabs(FsigN) > tolNR2 && j <= itNR2){
				sigNnew = sigNold - 2 * Pwr*auxroot*(sigNold - C1*(1 - auxroot) - C2) /
					(2 * Pwr*auxroot - C1*pow(sigNold, 1 / Pwr - 1));

				//Update parameters for next iteration (make sure sigNnew is not larger than 1 to avoid complex numbers)
				sigNorm = sigNold = fmin(sigNnew, 1 - 1E-5);
				
				// auxiliary variable 
				auxroot = sqrt(1 - pow(sigNorm, 1 / Pwr));

				// Convergence measure
				FsigN = sigNorm - C1*(1 - auxroot) - C2;

				//Increase counter
				j++;
			}
			//Compute converged normalized strain based on converged stress
			epsNorm = fmax(1 / C1 * (sigNorm - C2), 0);
		}
		// Compute true stress associated to eps_N
		sig_N = epsNorm*((pointB[1] - pointA[1]) - pointA[2] * (pointB[0] - pointA[0])) + 
				pointA[2] * (eps_N - pointA[0]) + pointA[1];

		//Compute slope in the normalized system. 
		double tanNorm = 2 * Pwr*pow(1 - pow(1 - epsNorm, 2), Pwr - 1)*(1 - epsNorm);
		
		// If the resulting slope is too large, the slope in natural coordinate matches the initial slope of the curve
		if (tanNorm > DBL_MAX)
			tan_N = pointA[2];
		else {
			// Compute slope in transformed coordinates
			double tanT = tanNorm*(((pointB[1] - pointA[1]) - pointB[2] * (pointB[0] - pointA[0]))*(pointA[2] - pointB[2])) /
				(pointA[2] * (pointB[0] - pointA[0]) - (pointB[1] - pointA[1]));
			//Compute tangent slope in natural coordinates
			tan_N = tanT*(pointA[2] - pointB[2]) / (pointA[2] - pointB[2] + tanT) + pointB[2];
		}
	}
}
// Using formulation in Kim & Koutromanos (2016)
void SteelDRC::bauschNURBS(double eps_N, double &sig_N, double &tan_N, double *pointA, double *pointB, double b) {
	double ea = pointA[0];
	double fa = pointA[1];
	double Ea = pointA[2];
	double eb = pointB[0];
	double fb = pointB[1];
	double Eb = pointB[2];
	// In case the slope of the asymptotes at both ends are the same, the bauschinger curve is a straight line
	if (fabs(Ea - Eb) < DBL_EPSILON) {
		sig_N = Ea * (eps_N - ea) + fa;
		tan_N = Ea;
		return;
	}
	// In case b = 0 (as for minor and simple reversals), the simple nurb between both ea and eb is sufficient for the Bauschinger effect
	if (b == 0) {
		nurbs(eps_N, sig_N, tan_N, pointA, pointB);
		return;
	}
	else {
		// Compute the strain of the characteristic point
		//double echar = (Eb*(1 - eb + ea) - fa) / ((Ea - Eb)*(eb - ea))*(eb - ea) + ea;
		double echar = (fb - fa + Ea * ea - Eb * eb) / (Ea - Eb);
		// Compute strains at point 1 and 2 in normalized coordinates
		//double e1pp = b * (2.0 - b)*(Eb*(1.0 - eb + ea) - fa) / ((Ea - Eb)*(eb - ea));
		double e1pp = b * (2.0 - b)*(echar - ea) / (eb - ea);
		double e2pp = fmin(1.0 - e1pp, 0.5);
		// Transform both strains into natural coordinates
		double e1 = e1pp * (eb - ea) + ea;
		double e2 = e2pp * (eb - ea) + ea;
		// Compute the corresponding stresses for 1 and 2
		double f1 = fa + Ea * (e1 - ea);
		double f2 = fb + Eb * (e2 - eb);
		// Save points 1 and 2 into three element arrays, containing the strain, stress and tangent modulus
		double point1[3] = { e1,f1,Ea };
		double point2[3] = { e2,f2,Eb };
		// Compute stress and tangent modulus at characteristic point
		double fchar;
		double Echar;
		nurbs(echar, fchar, Echar, point1, point2);

		// Save characteristic point as three element array, containing the strain, stress and tangent modulus
		double pointChar[3] = { echar,fchar,Echar };

		// Finally, compute the stress and tangent modulus corresponding to eps_N
		double signLoad = (eb >= ea) ? 1.0 : -1.0;

		if (signLoad*eps_N <= signLoad * echar) {
			nurbs(eps_N, sig_N, tan_N, pointA, pointChar);
		}
		else {
			nurbs(eps_N, sig_N, tan_N, pointChar, pointB);
		}
	}
}

void SteelDRC::bauschBezierCubic(double eps_N, double &sig_N, double &tan_N, double *pointA, double *pointB, double *xi, double *wcubic){
	// Extract location of end points of the curve
	double ea = pointA[0];
	double fa = pointA[1];
	double Ea = pointA[2];
	double eb = pointB[0];
	double fb = pointB[1];
	double Eb = pointB[2];

	// Define factor xi, indicating the location of the intermediate reference points of the curve (0<xi<1)
	//double xi[2] = { 0.9, 0.9 };

	// Compute location of the intermediate points in the curve
	double e1 = ea + xi[0] * (fb - fa - Eb*(eb - ea)) / (Ea - Eb);
	double f1 = fa + xi[0] * (Ea*Eb*(ea - eb) + Ea*(fb - fa)) / (Ea - Eb);
	double e2 = eb + xi[1] * (fb - fa - Ea*(eb - ea)) / (Ea - Eb);
	double f2 = fb + xi[1] * (Ea*Eb*(ea - eb) + Eb*(fb - fa)) / (Ea - Eb);

	// Consider the case where the given strain is too close to one of the end points
	// In case the asymptotes at both ends are the same, the bauschinger curve is a straight line
	if (fabs(Ea - Eb) < DBL_EPSILON) {
		sig_N = Ea * (eps_N - ea) + fa;
		tan_N = Ea;
		return;
	}
	// Consider the case eps_N is almost equal to one of the limit strains
	if (fabs(eps_N - ea) < DBL_EPSILON){
		sig_N = fa;
		tan_N = Ea;
		return;
	}
	else if (fabs(eps_N - eb) < DBL_EPSILON){
		sig_N = fb;
		tan_N = Eb;
		return;
	}
	// To compute the parametric variable "t" corresponding to the input strain eps_N, solve a cubic equation
	// with the following terms
	double a3 = -wcubic[0] * (ea - eps_N) + 3 * wcubic[1] * (e1 - eps_N) - 3 * wcubic[2] * (e2 - eps_N) + wcubic[3] * (eb - eps_N);
	double a2 = 3 * wcubic[0] * (ea - eps_N) - 6 * wcubic[1] * (e1 - eps_N) + 3 * wcubic[2] * (e2 - eps_N);
	double a1 = -3 * wcubic[0] * (ea - eps_N) + 3 * wcubic[1] * (e1 - eps_N);
	double a0 = wcubic[0] * (ea - eps_N);


	// Initialize parametric variable
	double t;
	// Consider the case the term a3 is very small, thus a quadratic equation needs to be solved instead
	if (fabs(a3) < fabs(a2)*1E-6) {
		double t1 = (-a1 + sqrt(pow(a1, 2.0) - 4.0 * a2*a0)) / (2.0 * a2);
		double t2 = (-a1 - sqrt(pow(a1, 2.0) - 4.0 * a2*a0)) / (2.0 * a2);
		if (0 <= t1 && t1 <= 1.0) 	t = t1;
		else if (0 <= t2 && t2 <= 1.0) 	t = t2;
		else t = NAN;
	}
	else {
		// Compute auxiliary variables to solve cubic equation, using Cardano's formulation
		double p = -a2 / (3 * a3);
		double q = a1 / (3 * a3) - pow(p, 2.0);
		double r = pow(p, 3.0) + (a1*a2 - 3 * a0*a3) / (6 * pow(a3, 2.0));
		// Compute discriminant under the square root in Cardano's formulation
		double disc = pow(q, 3.0) + pow(r, 2.0);

		
		double sign_r;
		// Extract sign of r, since cubic root of negative number can result in complex output
		if (fabs(r) < DBL_EPSILON) 	sign_r = 0.0;
		else if (r > 0) sign_r = 1.0;
		else sign_r = -1.0;

		// Compute value of parameter t based on the sign of the discriminant for the cubic equation	
		if (fabs(disc) < DBL_EPSILON) {
			double t1 = p + 2 * sign_r*pow(sign_r*r, 1.0 / 3.0);
			double t2 = p - sign_r * pow(sign_r*r, 1.0 / 3.0);
			// Select value of t in the [0,1] interval
			if (0.0 <= t1 && t1 <= 1.0) t = t1;
			else if (0.0 <= t2 && t2 <= 1.0) t = t2;
			else t = NAN;
		}
		else if (disc > 0) {
			double sign1 = (r + sqrt(disc) >= 0.0) ? 1.0 : -1.0;
			double sign2 = (r - sqrt(disc) >= 0.0) ? 1.0 : -1.0;
			t = p + sign1 * pow(sign1*(r + sqrt(disc)), 1.0 / 3.0) + sign2 * pow(sign2*(r - sqrt(disc)), 1.0 / 3.0);
		}
		else {
			double theta = acos(r / sqrt(-pow(q, 3.0)));
			// Compute the three real solutions of the cubic equation
			double t1 = 2 * sqrt(-q)*cos(theta / 3.0) + p;
			double t2 = 2 * sqrt(-q)*cos((theta + 2 * atan(1) * 4) / 3.0) + p;
			double t3 = 2 * sqrt(-q)*cos((theta + 4 * atan(1) * 4) / 3.0) + p;
			// Select value of t in the [0,1] interval
			if (0.0 <= t1 && t1 <= 1.0) t = t1;
			else if (0.0 <= t2 && t2 <= 1.0) t = t2;
			else if (0.0 <= t3 && t3 <= 1.0) t = t3;
			else t = NAN;
		}
	}


	// Define auxiliary vectors to compute the value of the stress and tangent modulus
	double b0[4] = { pow(1 - t,3.0),3 * t*pow(1 - t,2.0),3 * pow(t,2.0)*(1 - t),pow(t,3.0) };
	double b1[4] = { -3 * pow(1 - t,2.0),-6 * t*(1 - t) + 3 * pow(1 - t,2.0),-3 * pow(t,2.0) + 6 * t*(1 - t),3 * pow(t,2.0) };
	double ep[4] = { ea,e1,e2,eb };
	double fp[4] = { fa,f1,f2,fb };

	// Compute the stress 
	double sig_N_num = 0;
	double sig_N_den = 0;

	for (int n = 0; n < 4; n++) {
		sig_N_num += b0[n] * fp[n] * wcubic[n];
		sig_N_den += b0[n] * wcubic[n];
	}
	sig_N = sig_N_num / sig_N_den;

	// Compute tangent modulus
	// - Numerator
	double tan_N_num1 = 0;
	double tan_N_num2 = 0;

	for (int n = 0; n < 4; n++) {
		tan_N_num1 += b1[n] * fp[n] * wcubic[n];
		tan_N_num2 += b1[n] * wcubic[n];
	}
	double tan_N_num = tan_N_num1 * sig_N_den - sig_N_num * tan_N_num2;

	// - Denominator
	double tan_N_den1 = 0;
	double tan_N_den2 = 0;

	for (int n = 0; n < 4; n++) {
		tan_N_den1 += b1[n] * ep[n] * wcubic[n];
		tan_N_den2 += b0[n] * ep[n] * wcubic[n];
	}
	double tan_N_den = tan_N_den1 * sig_N_den - tan_N_den2 * tan_N_num2;

	tan_N = tan_N_num / tan_N_den;
}
void SteelDRC::nurbs(double eps_N, double &sig_N, double &tan_N, double *pointA, double *pointB) {
	// Function to compute NURBS interpolation
	double ea = pointA[0];
	double fa = pointA[1];
	double Ea = pointA[2];
	double eb = pointB[0];
	double fb = pointB[1];
	double Eb = pointB[2];
	// Compute strain in normalized coordinates
	double epp = (eps_N - ea) / (eb - ea);
	// compute normalized strain of interior point
	//double eipp = (Eb*(1 - eb + ea) - fa) / ((Ea - Eb)*(eb - ea));
	double eipp = (fb - fa - Eb * (eb - ea)) / ((Ea - Eb)*(eb - ea));
	// compute parameter u corresponding to epp
	double u = (-eipp + sqrt(pow(eipp, 2.0) + (1 - 2 * eipp)*epp)) / (1 - 2 * eipp);
	// Compute normalized stress of interior point
	double fipp = eipp * Ea*(eb - ea) / (fb - fa);
	// compute stress in normalized coordinates
	double fpp = 2 * u*(1 - u)*fipp + pow(u, 2.0);
	// Stress in natural coordinates
	sig_N = fpp*(fb- fa) + fa;
	// Tangent modulus in normalized coordinates
	double dfpp = ((1 - 2 * u)*fipp + u) / sqrt(pow(eipp,2.0)  + (1 - 2 * eipp)*epp);
	// Transform tangent modulus into natural coordinates
	tan_N = (fb - fa) / (eb - ea)*dfpp;
	return;
}

double SteelDRC::bezierWeightCubic(double omega){
	double wcubic;
	/*if (omega <= 0.06) {
		wcubic = fmin(145.3*exp(-70.63*omega), 5e-7*pow(omega, -5.35));
	}
	else {
		wcubic = 1e-14 / pow(omega, 10.91)*(4 + 2.5e7*pow(omega, 5.56));
	}*/
	wcubic = 643.4*exp(-100 * omega);
	return wcubic;
}
double SteelDRC::factorb(double omega) {
	// Function that computes the factor b required to define the shape of the NURBS curve representing the bauschinger effect (from Kim & Koutromanos (2016))
	double eyEng = fyEng / E; // Yield strain in engineering coordinates
	double fshN = fyEng*exp(eshN);// stress at onset of strain hardening in true coordinates

	// Compute characteristic values of NURBS used to define factor b
	double omega0 = 0.04;
	double omega1 = 0.069*(fuEng - fyEng) / fuEng + 0.0555;
	double omega3 = 0.0753*(fuEng - fyEng) / fuEng + 0.0691;
	double md = 17034 * eyEng - 85.66;
	double omegaj = (0.404- omega3*md) / (10.1 - md);
	double omega2 = (omega3 - omegaj)*(Psh - 1) / (84. + Psh) *omega1 / omegaj + omegaj;
	double b1 = 10.1*(omega1 - omega0);
	double b2 = 10.1*(omegaj - omega0);
	double b3 = md*(omega2 - omega3);

	// The  value of the return factor b depends on the interval where omega is
	if (omega <= omega0 || omega >= omega3) 
		return 0.0;
	if (omega > omega0 && omega <= omega1)
			return (10.1*(omega - omega0));
	if (omega > omega1 && omega <= omega2) {
		double u = (2 * (omega1 - omegaj) + sqrt(pow(2*(omegaj - omega1), 2.0) - 4.0*(omega2 - 2 * omegaj + omega1)*(omega1 - omega))) /
					(2 * (omega2 - 2 * omegaj + omega1));
		return (pow(1 - u, 2.0)*b1 + 2 * u*(1 - u)*b2 + pow(u, 2.0)*b3);
	}
	return (md*(omega - omega3));
}

void SteelDRC::State_Determination(int S, int K, int M, int Klmr, double Eun)
{
	// Initialize three element vectors (strain, stress and tangent modulus) to be used as the boundary conditions of the 
	// Bauschinger curves
	double ptA[3];
	double ptB[3];
	// Case 0: Material has fractured in tension
	if (Tfract == -1) {
		Tsig = 0;
		Ttan = 0;
		return;
	}
	// Case 1 : Material is in virgin state
	if (0 == Te0[0] && 0 == Te0[1]) {
		skeleton(Teps, Tsig, Ttan);
		return;
	}
	// Case 2 : Material is in the linear unloading branch
	if (Tlmr * Ter < Tlmr * Teps && Tlmr * Teps < Tlmr * Tea[Klmr])
	{
		Tsig = Eun * (Teps - Ter) + Tsr;
		Ttan = Eun;
		return;
	}
	// Case 3 : Stress strain curve rejoins shifter skeleton curve 
	if ((!isnan(Terejoin[K]) && (S * Teps > S * Terejoin[K])) || (!isnan(TerejoinL[K]) && (S * Teps > S * TerejoinL[K]))) {
		skeleton(Teps - Te0[K], Tsig, Ttan);
		return;
	}
	
	// ============================================================================================================
	// Log 2018/09/07:Case 3.2 :
	// ============================================================================================================
	// Case S*Tsa[K] > S*Tsrm[K] is considered,
	// it occurs when reversal following first major reversal does
	// not go far enough in the plastic region, thus the following 
	// linear unloading ends at a higher stress than the major 
	// bauschinger curve defining the boundary conditions of curves 
	// in the direction of the first major reversal, instead of the 
	// backbone curve, see State_Reversal code, case 5.

	if ((((S * (Tsa[K] + TErm[K]*(Term[K]-Tea[K]))) > S * Tsrm[K]) ||  ((S * (Tsa[K] + Eun * (Term[K] - Tea[K]))) < S * Tsrm[K])) && (TshOnset == 1)) {
		// Use Bauschinger curve between Tea[K] and the shifted ultimate strain
		ptA[0] = Tea[K];
		ptA[1] = Tsa[K];
		ptA[2] = Eun;

		ptB[0] = TerejoinL[K];
		ptB[1] = TsrejoinL[K];
		ptB[2] = TErejoinL[K];

		bauschMajor(bauschFlag, ptA, ptB, Team[K], Tsam[K], Te0[K], S, K,Teps,Tsig,Ttan);

		return;
	}
	// ============================================================================================================
	// End of 2018/09/07 log
	// ============================================================================================================

	// Case 4 : Strain is within the Bauschinger curve and
	// previous reversal was a major reversal or the strain Teps 
	// has surpassed the major/minor reversal strain in the 
	// current straining direction.
	
	if ((S*Teps > S*Term[K] || Tfrm[K] == 1) && isnan(Terejoin[K]))
	{
		// Define vectors containing the initial and final points/slopes of the
		// Bauschinger curve
		ptA[0] = Team[K];
		ptA[1] = Tsam[K];
		ptA[2] = Eun;
		ptB[0] = TerejoinL[K];
		ptB[1] = TsrejoinL[K];
		ptB[2] = TErejoinL[K];
		bauschMajor(bauschFlag, ptA, ptB, Team[K], Tsam[K], Te0[K], S, K, Teps, Tsig, Ttan);
		return;
	}
	// Case 5 : Bauschinger curve from Tea[K] to Terejoin[K].
	// Last reversal occured in yield plateau and the current 
	// state has not rejoined the shifted backbone curve
	if (!isnan(Terejoin[K])) {
		ptA[0] = Tea[K];
		ptA[1] = Tsa[K];
		ptA[2] = Eun;
		ptB[0] = Terejoin[K];
		ptB[1] = Tsrejoin[K];
		ptB[2] = TErejoin[K];
		bauschMinor(bauschFlag, ptA, ptB,0.04,Teps,Tsig,Ttan);
		return;
	}
	// Case 6 : Bauschinger curve from Tea[K] to Term[K]. Last reversal
	// is a minor or simple reversal and the current strain has not reached
	// the last major (or minor) reversal strain in the current straining 
	// direction (Term[K])
	ptA[0] = Tea[K];
	ptA[1] = Tsa[K];
	ptA[2] = Eun;
	ptB[0] = Term[K];
	ptB[1] = Tsrm[K];
	ptB[2] = TErm[K];
	bauschMinor(bauschFlag, ptA, ptB,0.04,Teps,Tsig,Ttan);
	return;
}

void SteelDRC::State_Reversal(int S,int K, int M, int &Klmr, double &Eun) {
	// Case 1: Reversal occured in linear-elastic zone of the skeleton curve
	// No parameter requires updating
	if ((Te0[0] == 0 && Te0[1] == 0) && (Ceps <= eyN && Ceps >= -eyN)) 
		return;
	// Case 2: Reversal occurred in the unloading branch. No parameter requires updating
	if ((Tlmr*Ter < Tlmr*Ceps) && (Tlmr*Ceps < Tlmr*Tea[Klmr]))
		return;

	// REVERSAL FROM BAUSCHINGER OR SHIFTED SKELETON CURVE
	// ------------------------------------------------------------------------------------------------
	// All following cases require to update the following state parameters
	Ter = Ceps;
	Tsr = Csig;
	TEr = Ctan;
	Tlmr = S;
	Klmr = K;
	Tsa[K] = Tsr + S * Dfu * fyEng;
	
	// If Te0max requires updating, so does Eun
	if (S*(Tsr / Eun - Ter)>Te0max) {
		Te0max = S*(Tsr / Eun - Ter);
		Eun = E*(0.82 + 1 / (5.55 + 1000 * Te0max));
	}
	// With Eun updated, compute the new Tea[K]
	Tea[K] = Ter + S*Dfu * fyEng / Eun;
	
	// Case 3: No strain hardening has occurred yet (but a reversal from yield plateau has 
	// ocurred already). Reversal occurred in Bauschinger curve
	if (TshOnset == 0 && !isnan(Terejoin[M]) && S*Ter > S*Terejoin[M])
		return;
	
	// Case 4: No strain hardening yet. Reversal from the yield Plateau, it will
	// rejoin the backbone curve in the opposite direction
	if (!isnan(Terejoin[M]) && S * Ter <= S * Terejoin[M] && TshOnset == 0) {
		Terejoin[M] = Ter;
		Tsrejoin[M] = Tsr;
		TErejoin[M] = TEr;
		// Update strain shift
		Te0[K] = Ter - Tsr / Eun;
		// Update rejoin strain in the opposite direction
		Terejoin[K] = Te0[K] + Te0[M] - Ter;

		// Update uniform strain point in the opposite direction
		TerejoinL[K] = S * euN + Te0[K];
		skeleton(TerejoinL[K] - Te0[K], TsrejoinL[K], TErejoinL[K]);
		// Reset onset of hardening flag
		TshOnset = 0;

		// Stress and slope at Terejoin[K] are calculated from the backbone curve
		skeleton(Terejoin[K] - Te0[K], Tsrejoin[K], TErejoin[K]);
		return;
	}

	// Case 5 : First major reversal from backbone curve.
	// (from strain hardening or post-localization branch in backbone curve)
	if (!isnan(Terejoin[M]) && TshOnset == 1) {
	//if (!isnan(Terejoin[M]) && !isnan(Terejoin[K]) && TshOnset == 1) { 
		// Save copy of Tfract flag for cases when skeleton curve is called
		int Tfract0;
		// Save reversal strain in case there is a reversal in the linear unloading branch and need to return through same path 
		Terejoin[M] = Ter;
		// Curve will not rejoin backbone curve in the opposite direction
		Terejoin[K] = NAN;

		Te0[K] = Ter - Tsr / Eun;

		// If localization point was exceeded pre-reversal, update localization point  
		if (-S * (Ter-Te0[M]) > euN) {
			TerejoinL[M] = Ter;
			TsrejoinL[M] = Tsr;
			TErejoinL[M] = fmax(TEr,0.0);
		}

		// Update the new target strain in direction after reversal
		double e_amp = -S * (TerejoinL[M]-Te0[M]);
		TerejoinL[K] = Te0[K] + S * e_amp;
		Tfract0 = Tfract;
		skeleton(TerejoinL[K] - Te0[K], TsrejoinL[K], TErejoinL[K]);
		Tfract = Tfract0;
		// slope of TErejoinL is limited to 0 to avoid bauschinger curves with intermediate peaks before reaching the ends of the curve
		if (TErejoinL[K] < 0)
			TErejoinL[K] = 0;

		// This being the first major reversal, I need to define a unloading branch followed by the bauschinger curve coinciding in stress with the most recent reversal point
		// Coordinates of end point of bauschinger curve are known (strain,stress and slope), but only the slope and line defining strain/stress of the initial point are known
		// The point is found iteratively using bisection algorithm (note this iterative process is only performed a limited number of times for the material model, most often only once, 
		// thus it should not negatively impact the processing speed of the model.)

		double eam_min = Te0[M];
		double eam_max = Te0[M] - S * fuN / Eun;
		double eam0 = 0.5 * (eam_min + eam_max);
		double sam0 = Eun * (eam0 - Te0[M]);
		double aux[3] = { Ter,Tsr,Eun };
		
		// Set initial point of bauschinger curve as estimate of end of unloading curve: (eam0,sam0,Eun)
		double ptA[3] = { eam0, sam0, Eun };
		// Final point if bauschinger curve matches the end of the shifted skeleton curve
		double ptB[3] = { TerejoinL[M], TsrejoinL[M], TErejoinL[M] };

		// Compute stress and slope of bauschinger curve at the reversal point
		bauschMajor(bauschFlag, ptA, ptB, eam0, sam0, Te0[M], -S, M, Ter, aux[1], aux[2]);

		// Compute stress difference between aux[1] and reversal point
		double diff = -S * (Tsr - aux[1]);
		double delta = fabs((Tsr - aux[1])/Tsr);
		// Define tolerance for algorithm
		double tol0 = 1E-8;

		// Run bisection method to find optimal ptA to match stress Tsr at Ter

		if (fabs(Tsr) > DBL_EPSILON) {
			while (delta > tol0 && fabs(eam_max-eam_min) > tol0) {
				if (diff > 0){
					eam_min = eam0;
					eam0 = 0.5 * (eam0 + eam_max);
				}
				else {
					eam_max = eam0;
					eam0 = 0.5 * (eam_min + eam0);
				}
					
				// Compute the stress corresponding to eam0
				sam0 = Eun * (eam0 - Te0[M]);
				// Redefine initial point of the bauschinger curve
				ptA[0] = eam0;
				ptA[1] = sam0;
				// Compute stress and slope of bauschinger curve at the reversal point
				bauschMajor(bauschFlag, ptA, ptB, eam0, sam0, Te0[M], -S, M, Ter, aux[1], aux[2]);
				diff = -S * (Tsr - aux[1]);
				delta = fabs((Tsr - aux[1]) / Tsr);
			}
		}

		// Define Tam and Tsam based on results of iterative process
		Team[M] = eam0;
		Tsam[M] = sam0;

		// Redefine vector with coordinates and slope of initial point in bauschinger curve
		ptA[0] = Team[M];
		ptA[1] = Tsam[M];
 		//Team[M] = Te0[M] - S * eyN;
		//Tsam[M] = -S * fyN;

		// Major reversal parameters are updated
		// - Update strain Term[M]
		Term[M] = Ter;
		// - Update stress and slope Tsrm[M] & TErm[M] respectively
		// Compute stress and slope from the bauschinger curve
		bauschMajor(bauschFlag, ptA, ptB, Team[M], Tsam[M], Te0[M], -S, M, Ter, Tsrm[M], TErm[M]);


		// Define coordinates of end of unloading branch after the reversal
		Team[K] = Tea[K];
		Tsam[K] = Tsa[K];
		// Activate minor reversal flag in the opposite direction
		Tfrm[K] = 1;
		// Set flag indicating current reversal is a major reversal
		Tfrm[M] = -1;
		return;
	}

	// Case 6 : Major reversal. Occurs if strain exceeds maximum 
	// historic strain in the current direction, or if stress 
	// difference between last major reversal and current reversal
	// is greater than 2*fyEng or if strain exceeds strain at the 
	// last major reversal in the current straining direction.
	if (S * (Ter - Tsr / Eun) < S * Te0[K] || S * (Tsrm[K] - Tsr) > 2 * Dfu* fyEng || (S * Ter < S * Term[M] && Tfrm[M] == -1)) {
		// Save copy of Tfract for instances where the skeleton function is called
		int Tfract0;
		
		// No more rejoing the skeleton curve at the yield plateau
		Terejoin[K] = NAN;
		Terejoin[M] = NAN; // 2021-10-23 : CHECK IF THIS LINE WORKS

		// Update localization point in direction pre-reversal if reversal strain exceeds previous shifted uniform strain
		if (-S*Ter > -S*TerejoinL[M]) {
			TerejoinL[M] = Ter;
			TsrejoinL[M] = Tsr;
			TErejoinL[M] = fmax(TEr,0.0);
		}
		// Update strain shift if maximum strain in the current direction has been exceeded (1st condition of Case 6)
		if (S * (Ter - Tsr / Eun) < S * Te0[K])
			Te0[K] = Ter - Tsr / Eun;
		// Update major reversal parameters
		Term[M] = Ter;
		Tsrm[M] = Tsr;
		TErm[M] = TEr;
		Team[K] = Tea[K];
		Tsam[K] = Tsa[K];

		if (-S * (TerejoinL[M] - Te0[M]) > (S * (TerejoinL[K] - Te0[K]))) {
			TerejoinL[K] = Te0[K] + Te0[M] - TerejoinL[M];
			Tfract0 = Tfract;
			skeleton(TerejoinL[K] - Te0[K], TsrejoinL[K], TErejoinL[K]);
			Tfract = Tfract0;
			if (TErejoinL[K] < 0.0)
				TErejoinL[K] = 0;
		}
		// Activate minor reversal flag in the opposite direction
		Tfrm[K] = 1;
		// Set flag indicating current reversal is a major reversal
		Tfrm[M] = -1;
		return;
	}
	// Case 7: Minor Reversal
	if (Tfrm[M] == 1 || S*Ter < S*Term[M]){
		// No more rejoining backbone curve at the yield plateau
		Terejoin[K] = NAN;
		Terejoin[M] = NAN; // 2021-10-23 : CHECK IF THIS LINE WORKS
		// Update reversal point opposite to last major reversal
		Term[M] = Ter;
		Tsrm[M] = Tsr;
		TErm[M] = TEr;
		// Deactivate minor reversal flag in the current direction
		Tfrm[M] = 0;
		return;
	}
}
