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
                                                                        
// $Revision: 1.6 $
// $Date: 2008-08-26 16:22:17 $
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
#include <Parameter.h>
#include <string.h>
#include <elementAPI.h>

void* OPS_BoucWenMaterial()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 10) {
	opserr << "WARNING: Insufficient arguements\n";
	opserr << "Want: uniaxialMaterial BoucWen tag? alpha? ko? n? gamma?" << endln 
	       << " beta? Ao? deltaA? deltaNu? deltaEta?" << endln;
	return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING invalid tag\n";
	return 0;
    }

    double data[10] = {0,0,0,0,0,0,0,0,0,1.0e-8};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 10) {
	numdata = 10;
    }
    if (OPS_GetDoubleInput(&numdata,data)) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }

    int maxNumIter = 20;
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 0) {
	numdata = 1;
	if (OPS_GetIntInput(&numdata,&maxNumIter) < 0) {
	    opserr << "WARNING invalid int inputs\n";
	    return 0;
	}
    }

    UniaxialMaterial* mat = new BoucWenMaterial(tag,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],maxNumIter);
    if (mat == 0) {
	opserr << "WARNING: failed to create Boucwenmaterial material\n";
	return 0;
    }

    return mat;
}

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

//SAJalali
BoucWenMaterial::BoucWenMaterial()
	:UniaxialMaterial(0, MAT_TAG_BoucWen)
{
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
	// SAJalali
	static Vector data(21);
	data(0) = alpha;
	data(1) = ko;
	data(2) = n;
	data(3) = gamma;
	data(4) = beta;
	data(5) = Ao;
	data(6) = deltaA;
	data(7) = deltaNu;
	data(8) = deltaEta;
	data(9) = Tstrain;
	data(10) = Cstrain;
	data(11) = Tz;
	data(12) = Cz;
	data(13) = Te;
	data(14) = Ce;
	data(15) = Tstress;
	data(16) = Ttangent;
	data(17) = tolerance;
	data(18) = maxNumIter;
	data(19) = this->getTag();
	data(20) = parameterID;

	if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
		opserr << "BoucWenMaterial::sendSelf() - failed to send Vector\n";
		return -1;
	}

	return 0;
}

int 
BoucWenMaterial::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	// SAJalali
	static Vector data(21);

	if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
		opserr << "BoucWenMaterial::recvSelf() - failed to recvSelf\n";
		return -1;
	}

	alpha = data(0);
	ko = data(1);
	n = data(2);
	gamma = data(3);
	beta = data(4);
	Ao = data(5);
	deltaA = data(6);
	deltaNu = data(7);
	deltaEta = data(8);
	Tstrain = data(9);
	Cstrain = data(10);
	Tz = data(11);
	Cz = data(12);
	Te = data(13);
	Ce = data(14);
	Tstress = data(15);
	Ttangent = data(16);
	tolerance = data(17);
	maxNumIter = data(18);
	this->setTag((int)data(19));
	parameterID = data(20);
	
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
BoucWenMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"alpha") == 0)
    return param.addObject(1, this);
  
  if (strcmp(argv[0],"ko") == 0)
    return param.addObject(2, this);
  
  if (strcmp(argv[0],"n") == 0)
    return param.addObject(3, this);
  
  if (strcmp(argv[0],"gamma") == 0)
    return param.addObject(4, this);
    
  if (strcmp(argv[0],"beta") == 0)
    return param.addObject(5, this);
  
  if (strcmp(argv[0],"Ao") == 0)
    return param.addObject(6, this);
  
  if (strcmp(argv[0],"deltaA") == 0)
    return param.addObject(7, this);
    
  if (strcmp(argv[0],"deltaNu") == 0)
    return param.addObject(8, this);
  
  if (strcmp(argv[0],"deltaEta") == 0)
    return param.addObject(9, this);

  return -1;
}

int
BoucWenMaterial::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case 1:
		this->alpha = info.theDouble;
		return 0;
	case 2:
		this->ko = info.theDouble;
		return 0;
	case 3:
		this->n = info.theDouble;
		return 0;
	case 4:
		this->gamma = info.theDouble;
		return 0;
	case 5:
		this->beta = info.theDouble;
		return 0;
	case 6:
		this->Ao = info.theDouble;
		return 0;
	case 7:
		this->deltaA = info.theDouble;
		return 0;
	case 8:
		this->deltaNu = info.theDouble;
		return 0;
	case 9:
		this->deltaEta = info.theDouble;
		return 0;
	default:
		return -1;
	}
}



int
BoucWenMaterial::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}




double
BoucWenMaterial::getStressSensitivity(int gradIndex, bool conditional)
{

	// Declare output variable
	double sensitivity = 0.0;


	// Issue warning if response is zero (then an error will occur)
	// changed by Quan & Michele 02-07-2006
	if (Tz == 0.0)  {
		if (Tstrain == 0.0) {
			sensitivity = 0.0;
			return sensitivity;
		}
		else {
		opserr << "ERROR: BoucWenMaterial::getStressSensitivity() is called " << endln
			<< " is called with zero hysteretic deformation Tz." << endln;
		}
	} // end changes by  Quan & Michele 02-07-2006

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
		DCz		 = (*SHVs)(0,gradIndex);
		DCe		 = (*SHVs)(1,gradIndex);
		DCstrain = (*SHVs)(2,gradIndex);
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
BoucWenMaterial::getTangentSensitivity(int gradIndex)
{
	return 0.0;
}

double
BoucWenMaterial::getDampTangentSensitivity(int gradIndex)
{
	return 0.0;
}

double
BoucWenMaterial::getStrainSensitivity(int gradIndex)
{
	return 0.0;
}

double
BoucWenMaterial::getRhoSensitivity(int gradIndex)
{
	return 0.0;
}


int
BoucWenMaterial::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
//  Quan & Michele Apr. 2006
	if (Tz == 0){return 0;}
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
		DCz		 = (*SHVs)(0,gradIndex);
		DCe		 = (*SHVs)(1,gradIndex);
		DCstrain = (*SHVs)(2,gradIndex);
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
	(*SHVs)(0,gradIndex) = DTz;
	(*SHVs)(1,gradIndex) = DTe;
	(*SHVs)(2,gradIndex) = DTstrain;

	return 0;
}


double
BoucWenMaterial::getInitialTangentSensitivity(int gradIndex)
{

	double dAlphadh =0.0;
	double dKodh =0.0;
	double dAodh=0.0;

	if (parameterID == 1) {
		dAlphadh = 1.0;
	}
	if (parameterID == 2) {
		dKodh =1.0;
	}
	
	if (parameterID == 6) {
		dAodh=1.0;
	}

    //return ( alpha*ko + (1-alpha)*ko*Ao );
	
	return dAlphadh*ko +alpha*dKodh -dAlphadh*ko*Ao+ (1-alpha)*dKodh*Ao+ (1-alpha)*ko*dAodh;
	
}


