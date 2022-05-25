#include <DelayMaterial.h>
#include <ID.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

void*
OPS_DelayMaterial(void)
{
	// pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	// get Input Values
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs != 4) {
		opserr << "WARNING! Invalid number of args in uniaxialMaterial Delay" << endln;
		opserr << "want: uniaxialMaterial Delay tag? tag0? ts? dts?" << endln;
		return 0;
	}

	int iData[2];
	double dData[2];

	int numData = 2;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING! Invalid tags for uniaxialMaterial Delay" << endln;
		return 0;
	}

	numData = 2;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING! Invalid time parameters for uniaxialMaterial Delay " << iData[0] << endln;
		return 0;
	}

	// create pointer to reference material model
	UniaxialMaterial* theMat = OPS_getUniaxialMaterial(iData[1]);

	if (theMat == 0) {
		opserr << "WARNING! Reference material model " << iData[1] << " does not exist" << endln;
		opserr << "uniaxialMaterial Delay: " << iData[0] << endln;
		
		delete theMat;
	}

	// allocate material model
	theMaterial = new DelayMaterial(iData[0], theMat, dData[0], dData[1]);

	if (theMaterial == 0) {
		opserr << "WARNING! Could not create uniaxialMaterial of type Delay\n";
		return 0;
	}

	return theMaterial;
}


DelayMaterial::DelayMaterial(int tag, UniaxialMaterial *theMaterialModel, double ts, double dts)
	:UniaxialMaterial(tag, MAT_TAG_Delay),
	trialStrain(0.0), trialStrainRate(0.0), theModel(0),
	t_s(ts), dt_s(dts), strain_s(0.0), active(0)
{
	theModel = theMaterialModel->getCopy();
}


DelayMaterial::DelayMaterial()
	:UniaxialMaterial(0, MAT_TAG_Delay),
	trialStrain(0.0), trialStrainRate(0.0), theModel(0),
	t_s(0.0), dt_s(0.0), strain_s(0.0), active(0)
{

}


DelayMaterial::~DelayMaterial()
{

}


int
DelayMaterial::setTrialStrain(double strain, double strainRate)
{
	int err = 0;

	trialStrain = strain;
	trialStrainRate = strainRate;

	// determine strain_s for material model
	if (active == 0 && ops_t >= t_s) {
		strain_s = trialStrain;
		err += theModel->setTrialStrain(trialStrain - strain_s, trialStrainRate);
	}

	// set trial strains in the material model if active
	if (active == 1)	// material replacement
		err += theModel->setTrialStrain(trialStrain - strain_s, trialStrainRate);

	return err;
}

double
DelayMaterial::getStrain(void)
{
	return trialStrain;
}

double
DelayMaterial::getStrainRate(void)
{
	return trialStrainRate;
}

double
DelayMaterial::getStress(void)
{
	// determine material model contribution factor
	double R = 0.0;

	if (active || (active == 0 && ops_t >= t_s)) {
		if (ops_t >= t_s + dt_s)
			R = 1.0;
		else if (ops_t >= t_s)
			R = (ops_t - t_s) / dt_s;
	}

	// get the stress
	return (R * theModel->getStress());
}



double
DelayMaterial::getTangent(void)
{
	// determine material model contribution factor
	double R = 0.0;

	if (active) {
		if (ops_t >= t_s + dt_s)
			R = 1.0;
		else if (ops_t >= t_s)
			R = (ops_t - t_s) / dt_s;
	}

	// get the tangent
	return (R * theModel->getTangent());
}

double
DelayMaterial::getInitialTangent(void)
{
	// material is not active initially
	return 0.0;
}

double
DelayMaterial::getDampTangent(void)
{
	// determine material model contribution factor
	double R = 0.0;

	if (active) {
		if (ops_t >= t_s + dt_s)
			R = 1.0;
		else if (ops_t >= t_s)
			R = (ops_t - t_s) / dt_s;
	}

	// get the dampTan
	return (R * theModel->getDampTangent());
}

int
DelayMaterial::commitState(void)
{
	//opserr << "committed\n" << endln;

	int err = 0;

	if (active == 0 && ops_t >= t_s) {
		//opserr << "NOTE DelayMaterial::commitState() - tag: " << this->getTag() << endln;
		//opserr << "Wrapped material model activated!" << endln;

		strain_s = trialStrain;
		active = true;
	}

	// invoke commitState() on active material
	if (active)
		err = theModel->commitState();

	if (err != 0) {
		opserr << "WARNING DelayMaterial::commitState() - tag: " << this->getTag() << endln;
		opserr << "- Wrapped material model failed to commitState():";
		theModel->Print(opserr);
	}

	return err;
}

int
DelayMaterial::revertToLastCommit(void)
{
	int err = 0;

	// make switched false if switch time has not been exceeded
	if (ops_t < t_s)
		active == false;

	// invoke revertToLastCommit() on active material
	if (active)
		err = theModel->revertToLastCommit();

	if (err != 0) {
		opserr << "WARNING DelayMaterial::revertToLastCommit() - tag: " << this->getTag() << endln;
		opserr << "- Wrapped material model failed to revertToLastCommit():";
		theModel->Print(opserr);
	}

	return err;
}


int
DelayMaterial::revertToStart(void)
{
	int err = 0;

	trialStrain = 0.0;
	trialStrainRate = 0.0;

	active = false;

	// invoke revertToStart() on MaterialModel object
	err = theModel->revertToStart();

	if (err != 0) {
		opserr << "WARNING DelayMaterial::revertToStart() - tag: " << this->getTag() << endln;
		opserr << "- Wrapped material model failed to revertToStart():";
		theModel->Print(opserr);
	}

	return err;
}


UniaxialMaterial *
DelayMaterial::getCopy(void)
{
	DelayMaterial *theCopy = new DelayMaterial(this->getTag(), theModel, t_s, dt_s);

	theCopy->trialStrain = trialStrain;
	theCopy->trialStrainRate = trialStrainRate;
	theCopy->strain_s = strain_s;
	theCopy->active = active;

	return theCopy;
}


int
DelayMaterial::sendSelf(int cTag, Channel &theChannel)
{
	opserr << "DelayMaterial::sendSelf() - not implemented!\n";
	return 0;
}

int
DelayMaterial::recvSelf(int cTag, Channel &theChannel,
FEM_ObjectBroker &theBroker)
{
	opserr << "DelayMaterial::recvSelf() - not implemented!\n";
	return 0;
}

void
DelayMaterial::Print(OPS_Stream &s, int flag)
{
	s << "Delay material tag: " << this->getTag() << endln;
	for (int i = 0; i < 2; i++) {
		s << " ";
		theModel->Print(s, flag);
	}
}
