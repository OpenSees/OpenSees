//
//
//////////////////////////////////////////////////////////////////////

#include "PlasticHardeningMaterial.h"
#include <math.h>

#define evolDebug 0
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PlasticHardeningMaterial::PlasticHardeningMaterial(int tag, int classTag)
:Material(tag, classTag), val_hist(0), val_trial(0),
 sFactor(1.0), residual(1.0)
{
}

PlasticHardeningMaterial::~PlasticHardeningMaterial()
{

}

int PlasticHardeningMaterial::setTrialValue(double xVal, double factor)
{
	sFactor = factor;
	val_trial = xVal;
	if(val_trial < 0) val_trial = 0;
	return 0;	
}

double PlasticHardeningMaterial::getTrialValue(void)
{
	return val_trial;
}

int PlasticHardeningMaterial::setTrialIncrValue(double dxVal)
{
	sFactor = 1.0;
	val_trial = val_hist + dxVal;
	if(val_trial < 0) val_trial = 0;
	return 0;	
}


void PlasticHardeningMaterial::setResidual(double res)
{
	residual = res;
}

int PlasticHardeningMaterial::commitState ()
{
	val_hist = val_trial;
	sFactor = 1.0;
	//cerr << "------ Ep value = " <<  val_hist << endln;
	return 0;	
}

int PlasticHardeningMaterial::revertToLastCommit (void)
{
	val_trial = val_hist;
	return 0;	
}

int PlasticHardeningMaterial::revertToStart (void)
{
	val_trial = 0;
	val_hist = 0;
	return 0;	
}


Response *PlasticHardeningMaterial::setResponse (char **argv, int argc, Information &matInformation)
{
	return 0;
}

int PlasticHardeningMaterial::getResponse (int responseID, Information &matInformation)
{
	return -1;
}

void PlasticHardeningMaterial::Print(OPS_Stream &s, int flag)
{
	s << "+Plastic Material: " << getTag() << endln;
	s << "|  sFactor =  " << sFactor << endln;
	s << "|  val_hist = " << val_hist << ", val_trial = " << val_trial << endln;
	return;
}


