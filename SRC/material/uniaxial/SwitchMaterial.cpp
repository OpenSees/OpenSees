#include <SwitchMaterial.h>
#include <ID.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>
#include <MaterialResponse.h>
#include <elementAPI.h>

void*
OPS_SwitchMaterial(void)
{
	// pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	// get Input Values
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs != 6) {
		opserr << "WARNING! Invalid number of args in uniaxialMaterial Switch" << endln;
		opserr << "want: uniaxialMaterial Switch tag? tag1? tag2? ts? dts? type?" << endln;
		return 0;
	}

	int iData[4];
	double dData[3];

	int numData = 3;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING! Invalid tags for uniaxialMaterial Switch" << endln;
		return 0;
	}

	numData = 2;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING! Invalid time parameters for uniaxialMaterial Switch " << iData[0] << endln;
		return 0;
	}

	const char* type = OPS_GetString();
	if (strcmp(type, "replace") == 0)
		iData[3] = 1;
	else if (strcmp(type, "modify") == 0)
		iData[3] = 2;

	// create an array of pointers to reference material models
	UniaxialMaterial** theMats = new UniaxialMaterial * [2];

	for (int i = 0; i < 2; i++) {
		UniaxialMaterial* theMat = OPS_getUniaxialMaterial(iData[i + 1]);

		if (theMat == 0) {
			opserr << "WARNING! Material model " << iData[i + 1] << " does not exist" << endln;
			opserr << "uniaxialMaterial Switch: " << iData[0] << endln;

			delete[] theMats;
		}
		else
			theMats[i] = theMat;
	}

	// allocate material model
	theMaterial = new SwitchMaterial(iData[0], theMats, dData[0], dData[1], iData[3]);

	if (theMaterial == 0) {
		opserr << "WARNING! Could not create uniaxialMaterial of type Switch\n";
		return 0;
	}

	delete[] theMats;

	return theMaterial;
}


SwitchMaterial::SwitchMaterial(int tag, UniaxialMaterial **theMaterialModels, double ts, double dts, int st)
	:UniaxialMaterial(tag, MAT_TAG_Switch),
	trialStrain(0.0), trialStrainRate(0.0), theModels(0),
	t_s(ts), dt_s(dts), sType(st), strain_s(0.0), strain_res(0.0), switched(0)
{
	// create an array (theModels) to store copies of the MaterialModels
	theModels = new UniaxialMaterial *[2];

	if (theModels == 0) {
		opserr << "FATAL SwitchMaterial::SwitchMaterial() ";
		opserr << "- ran out of memory for array of size: " << 2 << "\n";
		exit(-1);
	}

	// into the newly created array store a ponter to a copy
	// of the UniaxialMaterial stored in theMaterialModels
	for (int i = 0; i < 2; i++) {
		theModels[i] = theMaterialModels[i]->getCopy();
	}
}


SwitchMaterial::SwitchMaterial()
	:UniaxialMaterial(0, MAT_TAG_Switch),
	trialStrain(0.0), trialStrainRate(0.0), theModels(0),
	t_s(0.0), dt_s(0.0), sType(1), strain_s(0.0), strain_res(0.0), switched(0)
{

}


SwitchMaterial::~SwitchMaterial()
{
	// invoke the destructor on each MaterialModel object
	for (int i = 0; i < 2; i++)
		delete theModels[i];

	// now we can delete the array
	if (theModels != 0) // just in case blank constructor called and no recvSelf()
		delete[] theModels;
}


int
SwitchMaterial::setTrialStrain(double strain, double strainRate)
{
	int err = 0;

	trialStrain = strain;
	trialStrainRate = strainRate;

	// determine strain_s for new material model
	if (switched == 0) {
		if (sType == 1 && ops_t >= t_s + dt_s) {	// material replacement
				strain_s = trialStrain;
		}
		else if (sType == 2 && ops_t >= t_s) {	// material modification

			if (theModels[0]->getStress() < -1.e-12) {	// if current stress is compressive
				// iterate to find new model's strain for which its stress is equal to the old model's
				err += theModels[0]->setTrialStrain(trialStrain, trialStrainRate);	// set old model's strain
				double stress1 = theModels[0]->getStress();					// get old model's current stress

				double strain2 = 0.0;									// first guess for strain in new model that gives stress1
				err += theModels[1]->setTrialStrain(strain2, trialStrainRate);	// set new model's strain
				double stress2 = theModels[1]->getStress();				// get new model's current stress

				double tan2;							// new model's tangent
				double d_stress = stress2 - stress1;	// stress difference

				for (int i = 0; i < 20; i++) {	// iterate
					if (abs(d_stress) < 1.e-8) {		// check convergence
						strain_s = trialStrain - strain2;
						break;
					}
					else if (i == 19) {	// not converged
						opserr << "WARNING SwitchMaterial::setTrialStrain()\n";
						opserr << "could not find the new material's switch strain after 20 iterations\n\n";

						return -1;
					}

					// correct strain2 using Newton-Raphson
					tan2 = theModels[1]->getTangent();

					//opserr << i << '\t' << strain2 << '\t' << stress1 << '\t' << stress2 << '\t' << tan2 << endln;

					strain2 -= (d_stress / tan2);

					// update d_stress
					err += theModels[1]->setTrialStrain(strain2, trialStrainRate);	// set new model's strain
					stress2 = theModels[1]->getStress();				// get new model's current stress
					d_stress = stress2 - stress1;
				}
			}
			else {
				strain_s = strain_res;
			}
		}
	}

	// set trial strains in material models that are currently active
	if (sType == 1) {	// material replacement
		if (ops_t < t_s + dt_s)		// old model
			err += theModels[0]->setTrialStrain(trialStrain, trialStrainRate);
		else						// new model
			err += theModels[1]->setTrialStrain(trialStrain - strain_s, trialStrainRate);
	}
	else if (sType == 2) {	// material modification
		if (ops_t < t_s + dt_s)		// old model
			err += theModels[0]->setTrialStrain(trialStrain, trialStrainRate);
		
		if (ops_t >= t_s)	// new model
			err += theModels[1]->setTrialStrain(trialStrain - strain_s, trialStrainRate);
	}

	return err;
}

double
SwitchMaterial::getStrain(void)
{
	return trialStrain;
}

double
SwitchMaterial::getStrainRate(void)
{
	return trialStrainRate;
}

double
SwitchMaterial::getStress(void)
{
	// determine material model contribution factors
	double R1 = 1.0;
	double R2 = 0.0;

	if (sType == 1) {	// material replacement
		if (ops_t >= t_s + 2.0 * dt_s) {
			R1 = 0.0;
			R2 = 1.0;
		}
		else if (ops_t >= t_s + dt_s) {
			R1 = 0.0;
			R2 = (ops_t - (t_s + dt_s)) / dt_s;
		}
		else if (ops_t >= t_s) {
			R1 = 1.0 - (ops_t - t_s) / dt_s;
		}
	}
	else if (sType == 2) {	// material modification
		if (ops_t >= t_s + dt_s) {
			R1 = 0.0;
			R2 = 1.0;
		}
		else if (ops_t >= t_s) {
			R1 = 1.0 - (ops_t - t_s) / dt_s;
			R2 = 1.0 - R1;
		}
	}

	// get the stress = sum of stress in MaterialModel objects
	double stress = R1 * theModels[0]->getStress() + R2 * theModels[1]->getStress();

	return stress;
}



double
SwitchMaterial::getTangent(void)
{
	// determine material model contribution factors
	double R1 = 1.0;
	double R2 = 0.0;

	if (sType == 1) {	// material replacement
		if (ops_t >= t_s + 2.0 * dt_s) {
			R1 = 0.0;
			R2 = 1.0;
		}
		else if (ops_t >= t_s + dt_s) {
			R1 = 0.0;
			R2 = (ops_t - (t_s + dt_s)) / dt_s;
		}
		else if (ops_t >= t_s) {
			R1 = 1.0 - (ops_t - t_s) / dt_s;
		}
	}
	else if (sType == 2) {	// material modification
		if (ops_t >= t_s + dt_s) {
			R1 = 0.0;
			R2 = 1.0;
		}
		else if (ops_t >= t_s) {
			R1 = 1.0 - (ops_t - t_s) / dt_s;
			R2 = 1.0 - R1;
		}
	}

	// get the tangent = sum of tan in MaterialModel objects
	double tan = R1 * theModels[0]->getTangent() + R2 * theModels[1]->getTangent();

	return tan;
}

double
SwitchMaterial::getInitialTangent(void)
{
	// only old model's initial tangent
	return theModels[0]->getInitialTangent();
}

double
SwitchMaterial::getDampTangent(void)
{
	// determine material model contribution factors
	double R1 = 1.0;
	double R2 = 0.0;

	if (sType == 1) {	// material replacement
		if (ops_t >= t_s + 2.0 * dt_s) {
			R1 = 0.0;
			R2 = 1.0;
		}
		else if (ops_t >= t_s + dt_s) {
			R1 = 0.0;
			R2 = (ops_t - (t_s + dt_s)) / dt_s;
		}
		else if (ops_t >= t_s) {
			R1 = 1.0 - (ops_t - t_s) / dt_s;
		}
	}
	else if (sType == 2) {	// material modification
		if (ops_t >= t_s + dt_s) {
			R1 = 0.0;
			R2 = 1.0;
		}
		else if (ops_t >= t_s) {
			R1 = 1.0 - (ops_t - t_s) / dt_s;
			R2 = 1.0 - R1;
		}
	}

	// get the dampTan = sum of dampTan in MaterialModel objects
	double eta = R1 * theModels[0]->getDampTangent() + R2 * theModels[1]->getDampTangent();

	return eta;
}

int
SwitchMaterial::commitState(void)
{
	//opserr << "committed\n" << endln;

	int err0, err = 0;

	if (switched == 0) {
		if (sType == 1 && ops_t >= t_s + dt_s)	// material replacement
			switched = true;
		else if (sType == 2 && ops_t >= t_s)	// material modification
			switched = true;
	}

	// invoke commitState() on each active material
	if (sType == 1) {	// material replacement
		if (ops_t < t_s + dt_s) {		// old model
			err0 = theModels[0]->commitState();

			if (err0 != 0) {
				opserr << "WARNING SwitchMaterial::commitState() ";
				opserr << "- MaterialModel failed to commitState():";
				theModels[0]->Print(opserr);

				err += err0;
			}
		}
		else {						// new model
			err0 = theModels[1]->commitState();

			if (err0 != 0) {
				opserr << "WARNING SwitchMaterial::commitState() ";
				opserr << "- MaterialModel failed to commitState():";
				theModels[1]->Print(opserr);

				err += err0;
			}
		}
	}
	else if (sType == 2) {	// material modification
		if (ops_t < t_s + dt_s) {		// old model
			err0 = theModels[0]->commitState();

			if (err0 != 0) {
				opserr << "WARNING SwitchMaterial::commitState() ";
				opserr << "- MaterialModel failed to commitState():";
				theModels[0]->Print(opserr);

				err += err0;
			}
		}

		if (ops_t >= t_s) {	// new model
			err0 = theModels[1]->commitState();

			if (err0 != 0) {
				opserr << "WARNING SwitchMaterial::commitState() ";
				opserr << "- MaterialModel failed to commitState():";
				theModels[1]->Print(opserr);

				err += err0;
			}
		}

		// record strain_res
		if (ops_t < t_s && theModels[0]->getStress() < 0.0)
			strain_res = trialStrain;
	}

	return err;
}

int
SwitchMaterial::revertToLastCommit(void)
{
	// make switched false if switch time has not been exceeded
	if (sType == 1 && ops_t < t_s + dt_s)	// material replacement
		switched == false;
	else if (sType == 2 && ops_t < t_s)	// material modification
		switched == false;

	// invoke revertToLastCommit() on each active material
	if (sType == 1) {	// material replacement
		if (ops_t < t_s + dt_s)		// old model
			if (theModels[0]->revertToLastCommit() != 0) {
				opserr << "WARNING SwitchMaterial::revertToLastCommit() ";
				opserr << "- MaterialModel failed to revertToLastCommit():";
				theModels[0]->Print(opserr);
			}
			else						// new model
				if (theModels[1]->revertToLastCommit() != 0) {
					opserr << "WARNING SwitchMaterial::revertToLastCommit() ";
					opserr << "- MaterialModel failed to revertToLastCommit():";
					theModels[1]->Print(opserr);
				}
	}
	else if (sType == 2) {	// material modification
		if (ops_t < t_s + dt_s)		// old model
			if (theModels[0]->revertToLastCommit() != 0) {
				opserr << "WARNING SwitchMaterial::revertToLastCommit() ";
				opserr << "- MaterialModel failed to revertToLastCommit():";
				theModels[0]->Print(opserr);
			}

		if (ops_t >= t_s)	// new model
			if (theModels[1]->revertToLastCommit() != 0) {
				opserr << "WARNING SwitchMaterial::revertToLastCommit() ";
				opserr << "- MaterialModel failed to revertToLastCommit():";
				theModels[1]->Print(opserr);
			}
	}

	return 0;
}


int
SwitchMaterial::revertToStart(void)
{
	trialStrain = 0.0;
	trialStrainRate = 0.0;

	switched = false;

	// invoke revertToStart() on each of local MaterialModel objects
	for (int i = 0; i < 2; i++)
		if (theModels[i]->revertToStart() != 0) {
			opserr << "WARNING SwitchMaterial::revertToStart() ";
			opserr << "- MaterialModel failed to revertToStart():";
			theModels[i]->Print(opserr);
		}

	return 0;
}


UniaxialMaterial *
SwitchMaterial::getCopy(void)
{
	SwitchMaterial *theCopy = new SwitchMaterial(this->getTag(), theModels, t_s, dt_s, sType);

	theCopy->trialStrain = trialStrain;
	theCopy->trialStrainRate = trialStrainRate;
	theCopy->strain_s = strain_s;
	theCopy->strain_res = strain_res;
	theCopy->switched = switched;

	return theCopy;
}


int
SwitchMaterial::sendSelf(int cTag, Channel &theChannel)
{
	opserr << "SwitchMaterial::sendSelf() - not implemented!\n";
	return 0;
}

int
SwitchMaterial::recvSelf(int cTag, Channel &theChannel,
FEM_ObjectBroker &theBroker)
{
	opserr << "SwitchMaterial::recvSelf() - not implemented!\n";
	return 0;
}

void
SwitchMaterial::Print(OPS_Stream &s, int flag)
{
	s << "Switch material tag: " << this->getTag() << endln;
	for (int i = 0; i < 2; i++) {
		s << " ";
		theModels[i]->Print(s, flag);
	}
}
