// ExponReducing.cpp: implementation of the ExponReducing class.
//
//////////////////////////////////////////////////////////////////////

#include "ExponReducing.h"
#include <stdlib.h>

#define MAT_TAG_EXPON -1
#define DEBG 0
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ExponReducing::ExponReducing(int tag, double kp0, double alfa)
:PlasticHardeningMaterial(tag,MAT_TAG_EXPON),
  Kp0(kp0), alpha(alfa), resFactor(0.0)
{
}

ExponReducing::ExponReducing(int tag, double kp0, double alfa, double min_fact)
:PlasticHardeningMaterial(tag,MAT_TAG_EXPON),
  Kp0(kp0), alpha(alfa), resFactor(min_fact)
{
//	opserr << "ResFact = " <<  res_fact << endln; opserr << "\a";
}


ExponReducing::~ExponReducing()
{

}


double ExponReducing::getTrialPlasticStiffness()
{
	double K ;//= Kp0*exp(-1*val_trial*alpha);

	// if x0 and Kp0 is a const:
	// K = Kp0(1.0  - exp(-alpha*x0 + alpha*val_trial));	
	// K = Kp0*(1.0 - exp(-alpha + alpha*val_trial));
	
	// for pinching type stuff
	K = residual*Kp0*(1 - exp(-1*alpha*val_trial));

	if(sFactor != 1.0)
		K = Kp0*sFactor;
	
	if(K < (Kp0*resFactor))
		K = Kp0*resFactor;

//	opserr << "K = " << K << ", sFactor = " << sFactor << endln;
	
	if(K <0.0)
	{
		opserr << "Ri = " << val_trial << ", Factor = " << K/Kp0 << ", res_fact = " << resFactor << endln;
		opserr << "\a";
	}
	
	return K;
}


void ExponReducing::Print(OPS_Stream &s, int flag)
{
	s << "MultiLinear, Tag = " << getTag() << endln;
	s << "Kp0 = " << Kp0 << endln;
	s << "Alpha = " <<  alpha << endln;
}

PlasticHardeningMaterial *ExponReducing::getCopy(void)
{
 	PlasticHardeningMaterial *theMat = new ExponReducing(getTag(), Kp0, alpha, resFactor);
    return theMat;
}

