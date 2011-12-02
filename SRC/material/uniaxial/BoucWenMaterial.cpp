/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2003-03-06 18:34:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BoucWenMaterial.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <BoucWenMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Matrix.h>
#include <Information.h>

BoucWenMaterial::BoucWenMaterial(int tag, 
					double p_alpha,
					double p_ko,
					double p_n,
					double p_gamma,
					double p_beta,
					double p_Ao,
					double p_deltaA,
					double p_deltaNu,
					double p_deltaEta,
					double ptolerance,
					int pMaxNumIter)
:UniaxialMaterial(tag,MAT_TAG_BoucWen),
alpha(p_alpha), ko(p_ko), n(p_n), gamma(p_gamma), beta(p_beta), Ao(p_Ao), 
deltaA(p_deltaA), deltaNu(p_deltaNu), deltaEta(p_deltaEta), tolerance(ptolerance),
maxNumIter(pMaxNumIter)
{
	parameterID = 0;
	SHVs = 0;

	// Initialize variables
    this->revertToStart();
}


BoucWenMaterial::~BoucWenMaterial()
{
	if (SHVs != 0) 
		delete SHVs;
}




double 
BoucWenMaterial::signum(double value)
{
	if (value > 0.0) {
		return 1.0;
	}
	else {
		return -1.0;
	}
}




int 
BoucWenMaterial::setTrialStrain (double strain, double strainRate)
{
	// Set trial strain and compute strain increment
	Tstrain = strain;
	double dStrain = Tstrain - Cstrain;



	// Initial declarations (make sure not to declare class variables here!)
	double TA, Tnu, Teta, Psi, Phi, f, Te_, TA_, Tnu_;
	double Teta_, Phi_, f_, Tznew, Tzold, sign;


	// Newton-Raphson scheme to solve for z_{i+1} := z1
	int count = 0;
	double startPoint = 0.01;
	Tz = startPoint;
	Tzold = startPoint;
	Tznew = 1.0;  
	while ( ( fabs(Tzold-Tznew) > tolerance ) && count<maxNumIter) {

		Te = Ce + (1-alpha)*ko*dStrain*Tz;
		TA = Ao - deltaA*Te;
		Tnu = 1.0 + deltaNu*Te;
		Teta = 1.0 + deltaEta*Te;
		sign = signum(dStrain*Tz);
		Psi = gamma + beta*sign;
		Phi = TA - pow(fabs(Tz),n)*Psi*Tnu;
		f = Tz - Cz - Phi/Teta*dStrain;


		// Evaluate function derivative f' (underscore:=prime)
		Te_ = (1.0-alpha)*ko*dStrain;
		TA_ = -deltaA*Te_;
		Tnu_ = deltaNu*Te_;
		Teta_ = deltaEta*Te_;
		sign = signum(Tz);
		double pow1;
		double pow2;
		if (Tz == 0.0) {
			pow1 = 0.0;
			pow2 = 0.0;
		}
		else {
			pow1 = pow(fabs(Tz),(n-1));
			pow2 = pow(fabs(Tz),n);
		}
		Phi_ = TA_ - n*pow1*sign*Psi*Tnu - pow2*Psi*Tnu_;
		f_ = 1.0 - (Phi_*Teta-Phi*Teta_)/pow(Teta,2.0)*dStrain;


		// Issue warning if derivative is zero
		if ( fabs(f_)<1.0e-10 ) {
			opserr << "WARNING: BoucWenMaterial::setTrialStrain() -- zero derivative " << endln
				<< " in Newton-Raphson scheme" << endln;
		}

		// Take a Newton step
		Tznew = Tz - f/f_;


		// Update the root (but the keep the old for convergence check)
		Tzold = Tz; 
		Tz = Tznew;

		// Update counter
		count++;


		// Issue warning if we didn't converge
		if (count == maxNumIter) {

			opserr << "WARNING: BoucWenMaterial::setTrialStrain() -- did not" << endln
				<< " find the root z_{i+1}, after " << maxNumIter << " iterations" << endln
				<< " and norm: " << fabs(Tzold-Tznew) << endln;
		}

		// Compute stress
		Tstress = alpha*ko*Tstrain + (1-alpha)*ko*Tz;


		// Compute deterioration parameters
		Te = Ce + (1-alpha)*ko*dStrain*Tz;
		TA = Ao - deltaA*Te;
		Tnu = 1.0 + deltaNu*Te;
		Teta = 1.0 + deltaEta*Te;


		// Compute tangent
		if (Tz != 0.0) {
			Psi = gamma + beta*signum(dStrain*Tz);
			Phi = TA - pow(fabs(Tz),n)*Psi*Tnu;
			double b1, b2, b3, b4, b5;
			b1  = (1-alpha)*ko*Tz;
			b2  = (1-alpha)*ko*dStrain;
			b3  = dStrain/Teta;
			b4  = -b3*deltaA*b1 - b3*pow(fabs(Tz),n)*Psi*deltaNu*b1 
				- Phi/(Teta*Teta)*dStrain*deltaEta*b1 + Phi/Teta;
			b5  = 1.0 + b3*deltaA*b2 + b3*n*pow(fabs(Tz),(n-1))*signum(Tz)*Psi*Tnu
				+ b3*pow(fabs(Tz),n)*Psi*deltaNu*b2
				+ Phi/(Teta*Teta)*dStrain*deltaEta*b2;
			double DzDeps = b4/b5;
			Ttangent = alpha*ko + (1-alpha)*ko*DzDeps;
		}
		else {
			Ttangent = alpha*ko + (1-alpha)*ko;
		}

	}

    return 0;
}

double 
BoucWenMaterial::getStress(void)
{
    return Tstress;
}

double 
BoucWenMaterial::getInitialTangent(void)
{
    return ( alpha*ko + (1-alpha)*ko*Ao );
}

double 
BoucWenMaterial::getTangent(void)
{
    return Ttangent;
}

double 
BoucWenMaterial::getStrain(void)
{
    return Tstrain;
}

int 
BoucWenMaterial::commitState(void)
{
    // Commit trial history variables
    Cstrain = Tstrain;
	Cz = Tz;
	Ce = Te;

    return 0;
}

int 
BoucWenMaterial::revertToLastCommit(void)
{
	// Nothing to do here
    return 0;
}

int 
BoucWenMaterial::revertToStart(void)
{
    Tstrain = 0.0;
	Cstrain = 0.0;
	Tz = 0.0;
	Cz = 0.0;
	Te = 0.0;
	Ce = 0.0;
	Tstress = 0.0;
	Ttangent = alpha*ko + (1-alpha)*ko*Ao;

	if (SHVs != 0) 
		SHVs->Zero();

    return 0;
}

UniaxialMaterial *
BoucWenMaterial::getCopy(void)
{
    BoucWenMaterial *theCopy =
	new BoucWenMaterial(this->getTag(), alpha, ko, n, gamma,
						beta, Ao, deltaA, deltaNu, deltaEta,tolerance,maxNumIter);
    	
    theCopy->Tstrain = Tstrain;
    theCopy->Cstrain = Cstrain;
    theCopy->Tz = Tz;
    theCopy->Cz = Cz;
    theCopy->Te = Te;
    theCopy->Ce = Ce;
    theCopy->Tstress = Tstress;
    theCopy->Ttangent = Ttangent;

    return theCopy;
}

int 
BoucWenMaterial::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int 
BoucWenMaterial::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

void 
BoucWenMaterial::Print(OPS_Stream &s, int flag)
{
    s << "BoucWenMaterial, tag: " << this->getTag() << endln;
    s << "  alpha: " << alpha << endln;
    s << "  ko: " << ko << endln;
    s << "  n: " << n << endln;
    s << "  gamma: " << gamma << endln;
    s << "  beta: " << beta << endln;
    s << "  Ao: " << Ao << endln;
    s << "  deltaA: " << deltaA << endln;
    s << "  deltaNu: " << deltaNu << endln;
    s << "  deltaEta: " << deltaEta << endln;
}


int
BoucWenMaterial::setParameter(const char **argv, int argc, Information &info)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0],"alpha") == 0) {
		info.theType = DoubleType;
		return 1;
	}
	if (strcmp(argv[0],"ko") == 0) {
		info.theType = DoubleType;
		return 2;
	}
	if (strcmp(argv[0],"n") == 0) {
		info.theType = DoubleType;
		return 3;
	}
	if (strcmp(argv[0],"gamma") == 0) {
		info.theType = DoubleType;
		return 4;
	}
	if (strcmp(argv[0],"beta") == 0) {
		info.theType = DoubleType;
		return 5;
	}
	if (strcmp(argv[0],"Ao") == 0) {
		info.theType = DoubleType;
		return 6;
	}
	if (strcmp(argv[0],"deltaA") == 0) {
		info.theType = DoubleType;
		return 7;
	}
	if (strcmp(argv[0],"deltaNu") == 0) {
		info.theType = DoubleType;
		return 8;
	}
	if (strcmp(argv[0],"deltaEta") == 0) {
		info.theType = DoubleType;
		return 9;
	}
	else
		opserr << "WARNING: Could not set parameter in BoucWenMaterial. " << endln;
                
	return -1;
}

int
BoucWenMaterial::updateParameter(int parameterID, Information &info)
{

	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->alpha = info.theDouble;
		break;
	case 2:
		this->ko = info.theDouble;
		break;
	case 3:
		this->n = info.theDouble;
		break;
	case 4:
		this->gamma = info.theDouble;
		break;
	case 5:
		this->beta = info.theDouble;
		break;
	case 6:
		this->Ao = info.theDouble;
		break;
	case 7:
		this->deltaA = info.theDouble;
		break;
	case 8:
		this->deltaNu = info.theDouble;
		break;
	case 9:
		this->deltaEta = info.theDouble;
		break;
	default:
		return -1;
	}

	return 0;
}



int
BoucWenMaterial::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}




double
BoucWenMaterial::getStressSensitivity(int gradNumber, bool conditional)
{

	// Declare output variable
	double sensitivity = 0.0;


	// Issue warning if response is zero (then an error will occur)
	if (Tz == 0.0)  {
		opserr << "ERROR: BoucWenMaterial::getStressSensitivity() is called " << endln
			<< " is called with zero hysteretic deformation Tz." << endln;
	}

	// First set values depending on what is random
	double Dalpha = 0.0;
	double Dko = 0.0;
	double Dn = 0.0;
	double Dgamma = 0.0;
	double Dbeta = 0.0;
	double DAo = 0.0;
	double DdeltaA = 0.0;
	double DdeltaNu = 0.0;
	double DdeltaEta = 0.0;

	if (parameterID == 0) { }
	else if (parameterID == 1) {Dalpha=1.0;}
	else if (parameterID == 2) {Dko=1.0;}
	else if (parameterID == 3) {Dn=1.0;}
	else if (parameterID == 4) {Dgamma=1.0;}
	else if (parameterID == 5) {Dbeta=1.0;}
	else if (parameterID == 6) {DAo=1.0;}
	else if (parameterID == 7) {DdeltaA=1.0;}
	else if (parameterID == 8) {DdeltaNu=1.0;}
	else if (parameterID == 9) {DdeltaEta=1.0;}


	// Pick up sensitivity history variables for this gradient number
	double DCz = 0.0;
	double DCe = 0.0;
	double DCstrain = 0.0;
	if (SHVs != 0) {
		DCz		 = (*SHVs)(0,(gradNumber-1));
		DCe		 = (*SHVs)(1,(gradNumber-1));
		DCstrain = (*SHVs)(2,(gradNumber-1));
	}

	
	// Compute sensitivity of z_{i+1} 
	// (use same equations as for the unconditional 
	// sensitivities, just set DTstrain=0.0)
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11;
	double DTstrain = 0.0; 
	double dStrain = Tstrain-Cstrain;
	double TA, Tnu, Teta, DTz, Psi, Phi, DPsi;

	c1  = DCe 
		- Dalpha*ko*dStrain*Tz
		+ (1-alpha)*Dko*dStrain*Tz
		+ (1-alpha)*ko*(DTstrain-DCstrain)*Tz;
	c2  = (1-alpha)*ko*dStrain;
	c3  = DAo - DdeltaA*Te - deltaA*c1;
	c4  = -deltaA*c2;
	c5  = DdeltaNu*Te + deltaNu*c1;
	c6  = deltaNu*c2;
	c7  = DdeltaEta*Te + deltaEta*c1;
	c8  = deltaEta*c2;
	TA = Ao - deltaA*Te;
	Tnu = 1.0 + deltaNu*Te;
	Teta = 1.0 + deltaEta*Te;
	Psi = gamma + beta*signum(dStrain*Tz);
	DPsi= Dgamma + Dbeta*signum(dStrain*Tz);
	Phi = TA - pow(fabs(Tz),n)*Psi*Tnu;
	c9  = dStrain/Teta;
	c10 = DCz + c9*c3 - c9*pow(fabs(Tz),n)*Dn*log(fabs(Tz))*Psi*Tnu
		- c9*pow(fabs(Tz),n)*DPsi*Tnu - c9*pow(fabs(Tz),n)*Psi*c5
		- Phi/(Teta*Teta)*c7*dStrain + Phi/Teta*(DTstrain-DCstrain);
	c11 = 1.0 - c9*c4 + c9*pow(fabs(Tz),n)*Psi*c6
		+ c9*pow(fabs(Tz),n)*n/fabs(Tz)*signum(Tz)*Psi*Tnu
		+ Phi/(Teta*Teta)*c8*dStrain;

	DTz = c10/c11;

	sensitivity = Dalpha*ko*Tstrain
		        + alpha*Dko*Tstrain
				- Dalpha*ko*Tz
				+ (1-alpha)*Dko*Tz
				+ (1-alpha)*ko*DTz;

	return sensitivity;
}



double
BoucWenMaterial::getTangentSensitivity(int gradNumber)
{
	return 0.0;
}

double
BoucWenMaterial::getDampTangentSensitivity(int gradNumber)
{
	return 0.0;
}

double
BoucWenMaterial::getStrainSensitivity(int gradNumber)
{
	return 0.0;
}

double
BoucWenMaterial::getRhoSensitivity(int gradNumber)
{
	return 0.0;
}


int
BoucWenMaterial::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(3,numGrads);
	}

	// First set values depending on what is random
	double Dalpha = 0.0;
	double Dko = 0.0;
	double Dn = 0.0;
	double Dgamma = 0.0;
	double Dbeta = 0.0;
	double DAo = 0.0;
	double DdeltaA = 0.0;
	double DdeltaNu = 0.0;
	double DdeltaEta = 0.0;

	if (parameterID == 1) {Dalpha=1.0;}
	else if (parameterID == 2) {Dko=1.0;}
	else if (parameterID == 3) {Dn=1.0;}
	else if (parameterID == 4) {Dgamma=1.0;}
	else if (parameterID == 5) {Dbeta=1.0;}
	else if (parameterID == 6) {DAo=1.0;}
	else if (parameterID == 7) {DdeltaA=1.0;}
	else if (parameterID == 8) {DdeltaNu=1.0;}
	else if (parameterID == 9) {DdeltaEta=1.0;}


	// Pick up sensitivity history variables for this gradient number
	double DCz = 0.0;
	double DCe = 0.0;
	double DCstrain = 0.0;
	if (SHVs != 0) {
		DCz		 = (*SHVs)(0,(gradNumber-1));
		DCe		 = (*SHVs)(1,(gradNumber-1));
		DCstrain = (*SHVs)(2,(gradNumber-1));
	}

	
	// Compute sensitivity of z_{i+1} 
	// (use same equations as for the unconditional 
	// sensitivities, just set DTstrain=0.0)
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11;
	double DTstrain = TstrainSensitivity; 
	double dStrain = Tstrain-Cstrain;
	double TA, Tnu, Teta, DTz, Psi, Phi, DPsi, DTe;

	c1  = DCe 
		- Dalpha*ko*dStrain*Tz
		+ (1-alpha)*Dko*dStrain*Tz
		+ (1-alpha)*ko*(DTstrain-DCstrain)*Tz;
	c2  = (1-alpha)*ko*dStrain;
	c3  = DAo - DdeltaA*Te - deltaA*c1;
	c4  = -deltaA*c2;
	c5  = DdeltaNu*Te + deltaNu*c1;
	c6  = deltaNu*c2;
	c7  = DdeltaEta*Te + deltaEta*c1;
	c8  = deltaEta*c2;
	TA = Ao - deltaA*Te;
	Tnu = 1.0 + deltaNu*Te;
	Teta = 1.0 + deltaEta*Te;
	Psi = gamma + beta*signum(dStrain*Tz);
	DPsi= Dgamma + Dbeta*signum(dStrain*Tz);
	Phi = TA - pow(fabs(Tz),n)*Psi*Tnu;
	c9  = dStrain/Teta;
	c10 = DCz + c9*c3 - c9*pow(fabs(Tz),n)*Dn*log(fabs(Tz))*Psi*Tnu
		- c9*pow(fabs(Tz),n)*DPsi*Tnu - c9*pow(fabs(Tz),n)*Psi*c5
		- Phi/(Teta*Teta)*c7*dStrain + Phi/Teta*(DTstrain-DCstrain);
	c11 = 1.0 - c9*c4 + c9*pow(fabs(Tz),n)*Psi*c6
		+ c9*pow(fabs(Tz),n)*n/fabs(Tz)*signum(Tz)*Psi*Tnu
		+ Phi/(Teta*Teta)*c8*dStrain;

	DTz = c10/c11;

	DTe = DCe - Dalpha*ko*dStrain*Tz
		+ (1-alpha)*Dko*dStrain*Tz
		+ (1-alpha)*ko*(DTstrain-DCstrain)*Tz
		+ (1-alpha)*ko*dStrain*DTz;


	// Save sensitivity history variables
	(*SHVs)(0,(gradNumber-1)) = DTz;
	(*SHVs)(1,(gradNumber-1)) = DTe;
	(*SHVs)(2,(gradNumber-1)) = DTstrain;

	return 0;
}


