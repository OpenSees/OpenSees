// $Revision: 1.1 $
// $Date: 2000-12-19 03:35:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/FluidSolidPorousMaterial.cpp,v $
                                                                        
// Written: ZHY

//
// FluidSolidPorousMaterial.cpp
// -------------------
//
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include "FluidSolidPorousMaterial.h"
#include <Information.h>
#include <MaterialResponse.h>

int FluidSolidPorousMaterial::loadStage = 0;
double FluidSolidPorousMaterial::AtmoPress = 0.;
Vector * FluidSolidPorousMaterial::workV = 0;
Matrix * FluidSolidPorousMaterial::workM = 0;
const Vector zeroVector(6);

FluidSolidPorousMaterial::FluidSolidPorousMaterial (int tag, int nd, NDMaterial * soilMat,
                                      double combinedBulkModul, double atm)
 : NDMaterial(tag, MAT_TAG_FluidSolidPorousMaterial)
{
	if (combinedBulkModul < 0) {
		cerr << "WARNING:FluidSolidPorousMaterial::FluidSolidPorousMaterial: combinedBulkModulus < 0" << endl;
	  cerr << "Will reset to 0." <<endl;
    combinedBulkModul = 0.;
  }
	ndm = nd;
	loadStage = 0;  //default
  theSoilMaterial = soilMat;
	AtmoPress = atm;
  combinedBulkModulus = combinedBulkModul;
  trialExcessPressure = currentExcessPressure = 0.;
	trialVolumeStrain = currentVolumeStrain = 0.;
	if (ndm==2) {
		workV = new Vector(3);
		workM = new Matrix(3,3);
	} 
	else if (ndm==3) {
		workV = new Vector(6);
		workM = new Matrix(6,6);
	} 
}
   

FluidSolidPorousMaterial::FluidSolidPorousMaterial () 
 : NDMaterial(0,MAT_TAG_FluidSolidPorousMaterial), theSoilMaterial()
{
	ndm = 3; 
	combinedBulkModulus = 0.;
  trialExcessPressure = currentExcessPressure = 0.;
	trialVolumeStrain = currentVolumeStrain = 0.;
}


FluidSolidPorousMaterial::FluidSolidPorousMaterial (const FluidSolidPorousMaterial & a)
 : NDMaterial(a.getTag(),MAT_TAG_FluidSolidPorousMaterial)
{
	ndm = a.ndm;
	combinedBulkModulus = a.combinedBulkModulus;
  theSoilMaterial = a.theSoilMaterial->getCopy();
  trialExcessPressure = a.trialExcessPressure;
	currentExcessPressure = a.currentExcessPressure;
	trialVolumeStrain = a.trialVolumeStrain;
	currentVolumeStrain = a.currentVolumeStrain;
}


FluidSolidPorousMaterial::~FluidSolidPorousMaterial ()
{
  delete theSoilMaterial;
}


int FluidSolidPorousMaterial::setTrialStrain (const Vector &strain)
{
	if (ndm==2 && strain.Size()==3)
		trialVolumeStrain = strain[0]+strain[1];
	else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = strain[0]+strain[1]+strain[2];
	else {
		cerr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal("");
	}

  return theSoilMaterial->setTrialStrain(strain);
}


int FluidSolidPorousMaterial::setTrialStrain (const Vector &strain, const Vector &rate)
{
	if (ndm==2 && strain.Size()==3)
		trialVolumeStrain = strain[0]+strain[1];
	else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = strain[0]+strain[1]+strain[2];
	else {
		cerr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal("");
	}

  return theSoilMaterial->setTrialStrain(strain, rate);
}


int FluidSolidPorousMaterial::setTrialStrainIncr (const Vector &strain)
{
	if (ndm==2 && strain.Size()==3)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1];
	else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1]+strain[2];
	else {
		cerr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal("");
	}

  return theSoilMaterial->setTrialStrainIncr(strain);
}


int FluidSolidPorousMaterial::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
	if (ndm==2 && strain.Size()==3)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1];
	else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1]+strain[2];
	else {
		cerr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal("");
	}

  return theSoilMaterial->setTrialStrainIncr(strain, rate);
}


const Matrix & FluidSolidPorousMaterial::getTangent (void)
{
  *workM = theSoilMaterial->getTangent();

	if (loadStage != 0) { 
	  for (int i=0; i<ndm; i++) 
		  for (int j=0; j<ndm; j++) 
			  (*workM)(i,j) = (*workM)(i,j) + combinedBulkModulus;
  }

	return *workM;
}


const Vector & FluidSolidPorousMaterial::getStress (void)
{
  *workV = theSoilMaterial->getStress();

	if (loadStage != 0) { 
		trialExcessPressure = currentExcessPressure;
    trialExcessPressure += 
			       (trialVolumeStrain - currentVolumeStrain) * combinedBulkModulus;
	  for (int i=0; i<ndm; i++) 
      (*workV)[i] += trialExcessPressure;
  }

  return *workV;
}


int FluidSolidPorousMaterial::updateParameter(int responseID, Information &info)
{
	loadStage = responseID;
	return 0;
}


const Vector & FluidSolidPorousMaterial::getCommittedStress (void)
{

	*workV = theSoilMaterial->getCommittedStress();

	//if (loadStage != 0) { 
	//  for (int i=0; i<ndm; i++) 
  //    stress.theVector[i] += currentExcessPressure;
  //}

	return *workV;
}


const Vector & FluidSolidPorousMaterial::getCommittedStrain (void)
{
	*workV = theSoilMaterial->getCommittedStrain();
  return *workV;
}


double FluidSolidPorousMaterial::getCommittedPressure (void)
{
	return currentExcessPressure;
}

const Vector & FluidSolidPorousMaterial::getStrain (void)
{
  return theSoilMaterial->getStrain();
}


int FluidSolidPorousMaterial::commitState (void)
{
	currentVolumeStrain = trialVolumeStrain;
	if (loadStage != 0) 
		currentExcessPressure = trialExcessPressure;
	else
    currentExcessPressure = 0.;

	return theSoilMaterial->commitState();
}


int FluidSolidPorousMaterial::revertToLastCommit (void)
{
	return 0;
}


NDMaterial * FluidSolidPorousMaterial::getCopy (void)
{
  FluidSolidPorousMaterial * copy = new FluidSolidPorousMaterial(*this);
	return copy;
}


NDMaterial * FluidSolidPorousMaterial::getCopy (const char *code)
{
	if (strcmp(code,"FluidSolidPorous") == 0) {
     FluidSolidPorousMaterial * copy = new FluidSolidPorousMaterial(*this);
	   return copy;
	}

	return 0;
}


const char * FluidSolidPorousMaterial::getType (void) const
{
  return "FluidSolidPorous";
}


int FluidSolidPorousMaterial::getOrder (void) const
{
	return (ndm == 2) ? 3 : 6;
}


int FluidSolidPorousMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	// Need to implement
	return 0;
}


int FluidSolidPorousMaterial::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)    
{
	// Need to implement
	return 0;
}



Response*
FluidSolidPorousMaterial::setResponse (char **argv, int argc, Information &matInfo)
{
    if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getCommittedStress());

    else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getCommittedStrain());
    
	else if (strcmp(argv[0],"tangent") == 0)
		return new MaterialResponse(this, 3, this->getTangent());

    else if (strcmp(argv[0],"pressure") == 0)
		return new MaterialResponse(this, 4, this->getCommittedPressure());

    else
		return 0;
}


int FluidSolidPorousMaterial::getResponse (int responseID, Information &matInfo)
{
	switch (responseID) {
		case 1:
			return matInfo.setVector(this->getCommittedStress());

		case 2:
			return matInfo.setVector(this->getCommittedStrain());

		case 3:
			return matInfo.setMatrix(this->getTangent());
			
		case 4:
			return matInfo.setDouble(this->getCommittedPressure());

		default:
			return -1;
	}
}


void FluidSolidPorousMaterial::Print(ostream &s, int flag )
{
	s << "FluidSolidPorousMaterial" << endl;
}










