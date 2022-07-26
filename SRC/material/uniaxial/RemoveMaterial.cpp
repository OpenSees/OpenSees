/* Written by: Mohammad Salehi (mohammad.salehi@rice.edu)
** Created: 2020
** Description: The source code for Remove material model used to simulate repaired RC members
**
**
** Reference:
**
** Mohammad Salehi, Petros Sideris, and Reginald DesRoches (2022)
** “Numerical modeling of repaired reinforced concrete bridge columns”
** Engineering Structures, 253: 113801
*/

#include <RemoveMaterial.h>
#include <ID.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

void*
OPS_RemoveMaterial(void)
{
	// pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	// get Input Values
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs != 4) {
		opserr << "WARNING! Invalid number of args in uniaxialMaterial Remove" << endln;
		opserr << "want: uniaxialMaterial Delay tag? tag0? ts? dts?" << endln;
		return 0;
	}

	int iData[2];
	double dData[2];

	int numData = 2;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING! Invalid tags for uniaxialMaterial Remove" << endln;
		return 0;
	}

	numData = 2;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING! Invalid time parameters for uniaxialMaterial Remove " << iData[0] << endln;
		return 0;
	}

	// create pointer to reference material model
	UniaxialMaterial* theMat = OPS_getUniaxialMaterial(iData[1]);

	if (theMat == 0) {
		opserr << "WARNING! Reference material model " << iData[1] << " does not exist" << endln;
		opserr << "uniaxialMaterial Remove: " << iData[0] << endln;

		delete theMat;
	}

	// allocate material model
	theMaterial = new RemoveMaterial(iData[0], theMat, dData[0], dData[1]);

	if (theMaterial == 0) {
		opserr << "WARNING! Could not create uniaxialMaterial of type Remove\n";
		return 0;
	}

	return theMaterial;
}


RemoveMaterial::RemoveMaterial(int tag, UniaxialMaterial *theMaterialModel, double ts, double dts)
	:UniaxialMaterial(tag, MAT_TAG_Remove),
	trialStrain(0.0), trialStrainRate(0.0), theModel(0),
	t_s(ts), dt_s(dts)
{
	
	theModel = theMaterialModel->getCopy();
}


RemoveMaterial::RemoveMaterial()
	:UniaxialMaterial(0, MAT_TAG_Remove),
	trialStrain(0.0), trialStrainRate(0.0), theModel(0),
	t_s(0.0), dt_s(0.0)
{

}


RemoveMaterial::~RemoveMaterial()
{

}


int
RemoveMaterial::setTrialStrain(double strain, double strainRate)
{
	int err = 0;

	trialStrain = strain;
	trialStrainRate = strainRate;

	// set trial strains in the material model if active
	if (ops_t < t_s + dt_s)
		err += theModel->setTrialStrain(trialStrain, trialStrainRate);

	return err;
}

double
RemoveMaterial::getStrain(void)
{
	return trialStrain;
}

double
RemoveMaterial::getStrainRate(void)
{
	return trialStrainRate;
}

double
RemoveMaterial::getStress(void)
{
	// determine material model contribution factor
	double R = 1.0;

	if (ops_t >= t_s + dt_s)
		R = 0.0;
	else if (ops_t >= t_s)
		R = 1.0 - (ops_t - t_s) / dt_s;

	// get the stress
	return (R * theModel->getStress());
}



double
RemoveMaterial::getTangent(void)
{
	// determine material model contribution factor
	double R = 1.0;

	if (ops_t >= t_s + dt_s)
		R = 0.0;
	else if (ops_t >= t_s)
		R = 1.0 - (ops_t - t_s) / dt_s;

	// get the tangent
	return (R * theModel->getTangent());
}

double
RemoveMaterial::getInitialTangent(void)
{
	// material is active initially
	return theModel->getInitialTangent();
}

double
RemoveMaterial::getDampTangent(void)
{
	// determine material model contribution factor
	double R = 1.0;

	if (ops_t >= t_s + dt_s)
		R = 0.0;
	else if (ops_t >= t_s)
		R = 1.0 - (ops_t - t_s) / dt_s;

	// get the dampTan
	return (R * theModel->getDampTangent());
}

int
RemoveMaterial::commitState(void)
{
	//opserr << "committed\n" << endln;

	int err = 0;

	// invoke commitState() on active material
	if (ops_t < t_s + dt_s)
		err = theModel->commitState();

	if (err != 0) {
		opserr << "WARNING RemoveMaterial::commitState() ";
		opserr << "- MaterialModel failed to commitState():";
		theModel->Print(opserr);
	}

	return err;
}

int
RemoveMaterial::revertToLastCommit(void)
{
	int err = 0;

	// invoke revertToLastCommit() on active material
	if (ops_t < t_s + dt_s)
		err = theModel->revertToLastCommit();

	if (err != 0) {
		opserr << "WARNING RemoveMaterial::revertToLastCommit() ";
		opserr << "- MaterialModel failed to revertToLastCommit():";
		theModel->Print(opserr);
	}

	return err;
}


int
RemoveMaterial::revertToStart(void)
{
	int err = 0;

	trialStrain = 0.0;
	trialStrainRate = 0.0;

	// invoke revertToStart() on MaterialModel object
	err = theModel->revertToStart();

	if (err != 0) {
		opserr << "WARNING RemoveMaterial::revertToStart() ";
		opserr << "- MaterialModel failed to revertToStart():";
		theModel->Print(opserr);
	}

	return err;
}


UniaxialMaterial *
RemoveMaterial::getCopy(void)
{
	RemoveMaterial *theCopy = new RemoveMaterial(this->getTag(), theModel, t_s, dt_s);

	theCopy->trialStrain = trialStrain;
	theCopy->trialStrainRate = trialStrainRate;

	return theCopy;
}


int
RemoveMaterial::sendSelf(int cTag, Channel &theChannel)
{
	opserr << "RemoveMaterial::sendSelf() - not implemented!\n";
	return 0;
}

int
RemoveMaterial::recvSelf(int cTag, Channel &theChannel,
FEM_ObjectBroker &theBroker)
{
	opserr << "RemoveMaterial::recvSelf() - not implemented!\n";
	return 0;
}

void
RemoveMaterial::Print(OPS_Stream &s, int flag)
{
	s << "Remove material tag: " << this->getTag() << endln;
	for (int i = 0; i < 2; i++) {
		s << " ";
		theModel->Print(s, flag);
	}
}
