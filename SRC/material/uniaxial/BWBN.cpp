#include <elementAPI.h>
#include "BWBN.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>


void *
OPS_BWBN(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData1[1];
  double dData[13];
  int    iData2[1];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData1) != 0) {
    opserr << "WARNING invalid uniaxialMaterial BWBN tag" << endln;
    return 0;
  }

  numData = 13;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Double Values\n";
    return 0;	
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, iData2) != 0) {
    opserr << "WARNING invalid maxNumIter" << endln;
    return 0;
  }

  theMaterial = new BWBN(iData1[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
			 dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],iData2[0]);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type BWBN\n";
    return 0;
  }

  return theMaterial;
}



BWBN::BWBN(int tag, 
	   double p_alpha,
	   double p_ko,
	   double p_n,
	   double p_gamma,
	   double p_beta,
	   double p_Ao,
	   double p_q,
	   double p_zetas,
	   double p_p,
	   double p_Shi,
	   double p_deltaShi,
	   double p_lamda,
	   double ptolerance,
	   int pMaxNumIter)
 :UniaxialMaterial(tag,MAT_TAG_BWBN),
  alpha(p_alpha), ko(p_ko), n(p_n), gamma(p_gamma), beta(p_beta), Ao(p_Ao), q(p_q), 
  zetas(p_zetas), p(p_p), Shi(p_Shi), deltaShi(p_deltaShi), lamda(p_lamda), tolerance(ptolerance),
  maxNumIter(pMaxNumIter)
{
  // Initialize variables
  this->revertToStart();
}


BWBN::~BWBN()
{
  // does nothing
}

double 
BWBN::signum(double value)
{
  if (value > 0.0) {
    return 1.0;
  }
  else {
    return -1.0;
  }
}


int 
BWBN::setTrialStrain (double strain, double strainRate)
{
	// Set trial strain and compute strain increment
	Tstrain = strain;
	double dStrain = Tstrain - Cstrain;



	// Initial declarations (make sure not to declare class variables here!)
	double Psi, Phi, f, Te_;
	double Phi_, f_, Tznew, Tzold, sign;
	double Tzeta1, Tzeta2, Tzeta1_, Tzeta2_, h, h_,zu;


	// Newton-Raphson scheme to solve for z_{i+1} := z1
	int count = 0;
	double startPoint = 0.01;
	Tz = startPoint;
	Tzold = startPoint;
	Tznew = 1.0;  
	while ( ( fabs(Tzold-Tznew) > tolerance ) && count<maxNumIter) {

		Te = Ce + (1.0-alpha)*ko*dStrain*Tz;
		sign = signum(dStrain*Tz);
		Tzeta1 = zetas*(1-exp(-p*Te));
		Tzeta2 = (Shi+deltaShi*Te)*(lamda+Tzeta1);
		zu = pow(1/(beta+gamma),1/n);
		h = 1.0-Tzeta1*exp(-pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2));
		Psi = gamma + beta*sign;
		Phi = Ao - pow(fabs(Tz),n)*Psi;
		f = Tz - Cz - Phi*h*dStrain;


		// Evaluate function derivative f' (underscore:=prime)
		Te_ = (1.0-alpha)*ko*dStrain;
		Tzeta1_ = zetas*p*exp(-p*Te)*Te_;
		Tzeta2_ = Shi*Tzeta1_+lamda*deltaShi*Te_+deltaShi*Te*Tzeta1_+deltaShi*Te_*Tzeta1;
		//h_ = -exp(-pow(Tz*signum(dStrain)-q*zu,2)/(Tzeta2*Tzeta2))*(Tzeta1_-Tzeta1*2*(Tz*signum(dStrain)-q*zu)/(Tzeta2*Tzeta2)+Tzeta2_*Tzeta1*2*pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2*Tzeta2));
		h_ = -exp(-pow(Tz*signum(dStrain)-q*zu,2)/(Tzeta2*Tzeta2))*(Tzeta1_-Tzeta1*2*(Tz*signum(dStrain)-q*zu)*signum(dStrain)/(Tzeta2*Tzeta2)+Tzeta2_*Tzeta1*2*pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2*Tzeta2));
		sign = signum(Tz);
		double pow1;
		if (Tz == 0.0) {
			pow1 = 0.0;
		}
		else {
			pow1 = pow(fabs(Tz),(n-1));
		}
		Phi_ = - n*pow1*sign*Psi;
		f_ = 1.0 - (Phi_*h+Phi*h_)*dStrain;


		// Issue warning if derivative is zero
		if ( fabs(f_)<1.0e-10 ) {
			opserr << "WARNING: BWBN::setTrialStrain() -- zero derivative " << endln
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

			opserr << "WARNING: BWBN::setTrialStrain() -- did not" << endln
				<< " find the root z_{i+1}, after " << maxNumIter << " iterations" << endln
				<< " and norm: " << fabs(Tzold-Tznew) << endln;
		}

		// Compute stress
		Tstress = alpha*ko*Tstrain + (1-alpha)*ko*Tz;


		// Compute deterioration parameters
		Te = Ce + (1-alpha)*ko*dStrain*Tz;
		Tzeta1 = zetas*(1-exp(-p*Te));
		Tzeta2 = (Shi+deltaShi*Te)*(lamda+Tzeta1);
		

		// Compute tangent
		if (Tz != 0.0) {
			Psi = gamma + beta*signum(dStrain*Tz);
			Phi = Ao - pow(fabs(Tz),n)*Psi;
			double b1, b2, b3, b4, b5, b6, b7, b8, b9;
			b1 = (1-alpha)*ko*Tz;
			b2 = zetas*p*exp(-p*Te)*b1;
			b3 = Shi*b2+lamda*deltaShi*b1+deltaShi*Te*b2+deltaShi*b1*Tzeta1;
			b4 = -exp(-pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2))*(b2+b3*Tzeta1*2*pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2*Tzeta2)); 
			h = 1.0-Tzeta1*exp(-pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2));
			
			b5 = (1.0-alpha)*ko*dStrain;
			b6 = zetas*p*exp(-p*Te)*b5;
			b7 = Shi*b6+lamda*deltaShi*b5+deltaShi*Te*b6+deltaShi*b5*Tzeta1;
			//b8 = -exp(-pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2))*(b6-Tzeta1*2*(Tz*signum(dStrain)-q*zu)/(Tzeta2*Tzeta2)+b7*Tzeta1*2*pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2*Tzeta2));
			b8 = -exp(-pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2))*(b6-Tzeta1*2*(Tz*signum(dStrain)-q*zu)*signum(dStrain)/(Tzeta2*Tzeta2)+b7*Tzeta1*2*pow((Tz*signum(dStrain)-q*zu),2)/(Tzeta2*Tzeta2*Tzeta2));
			sign = signum(Tz);
			pow1 = pow(fabs(Tz),(n-1));
			b9 = - n*pow1*sign*Psi;
			double DzDeps = (h*Phi-b4*Phi)/(1.0 - (b9*h+Phi*b8)*dStrain);
			Ttangent = alpha*ko + (1-alpha)*ko*DzDeps;
			//Ttangent = Tstress/Tstrain;
		}
		else {
			Ttangent = alpha*ko + (1-alpha)*ko;
		}

	}

    return 0;
}

double 
BWBN::getStress(void)
{
    return Tstress;
}

double 
BWBN::getInitialTangent(void)
{
    return ( alpha*ko + (1-alpha)*ko*Ao );
}


double 
BWBN::getTangent(void)
{
    return Ttangent;
}

double 
BWBN::getStrain(void)
{
    return Tstrain;
}

int 
BWBN::commitState(void)
{
    // Commit trial history variables
    Cstrain = Tstrain;
	Cz = Tz;
	Ce = Te;

    return 0;
}

int 
BWBN::revertToLastCommit(void)
{
	// Nothing to do here
    return 0;
}

int 
BWBN::revertToStart(void)
{
    Tstrain = 0.0;
	Cstrain = 0.0;
	Tz = 0.0;
	Cz = 0.0;
	Te = 0.0;
	Ce = 0.0;
	Tstress = 0.0;
	Ttangent = alpha*ko + (1-alpha)*ko*Ao;

    return 0;
}

UniaxialMaterial *
BWBN::getCopy(void)
{
    BWBN *theCopy =
	new BWBN(this->getTag(), alpha, ko, n, gamma,
						beta, Ao, q, zetas, p, Shi, deltaShi, lamda,tolerance,maxNumIter);
    	
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
BWBN::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int 
BWBN::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

void 
BWBN::Print(OPS_Stream &s, int flag)
{
    s << "BWBN, tag: " << this->getTag() << endln;
    s << "  alpha: " << alpha << endln;
    s << "  ko: " << ko << endln;
    s << "  n: " << n << endln;
    s << "  gamma: " << gamma << endln;
    s << "  beta: " << beta << endln;
    s << "  Ao: " << Ao << endln;
	s << "  q: " << q << endln;
    s << "  deltaA: " << zetas << endln;
    s << "  deltaNu: " << p << endln;
    s << "  deltaEta: " << Shi << endln;
	s << "  deltaNu: " << deltaShi << endln;
    s << "  deltaEta: " << lamda << endln;
}


