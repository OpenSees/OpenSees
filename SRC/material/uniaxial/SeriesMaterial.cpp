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
                                                                        
// $Revision: 1.7 $
// $Date: 2004-11-24 00:48:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SeriesMaterial.cpp,v $

// Written: MHS
// Created: Sept 2000
//
// Description: This file contains the class definition for 
// SeriesModel. SeriesModel is an aggregation
// of UniaxialMaterial objects all considered acting in Series.
// Uses the same state determination as the beam elements.
// b = [1 1 ... 1]^

#include <SeriesMaterial.h>
#include <ID.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>
#include <MaterialResponse.h>

#include <OPS_Globals.h>

SeriesMaterial::SeriesMaterial(int tag, int num,
			       UniaxialMaterial ** theMaterialModels,
			       int maxIter, double tol)
:UniaxialMaterial(tag,MAT_TAG_SeriesMaterial),
 Tstrain(0.0), Cstrain(0.0), Tstress(0.0), Cstress(0.0),
 Ttangent(0.0), Ctangent(0.0), 
 maxIterations(maxIter), tolerance(tol),
 stress(0), flex(0), strain(0), initialFlag(false),
 numMaterials(num), theModels(0)
{
    theModels = new UniaxialMaterial *[numMaterials];

    if (theModels == 0) {
      opserr << "SeriesMaterial::SeriesMaterial -- failed to allocate material array\n";
      exit(-1);
    }


    int i;
    for (i = 0; i < numMaterials; i++) {
      theModels[i] = theMaterialModels[i]->getCopy();
      if (theModels[i] == 0) {
	opserr << "SeriesMaterial::SeriesMaterial -- failed to get copy of material: " << i << endln;
	exit(-1);
      }
    }

    strain = new double [numMaterials];
    if (strain == 0) {
      opserr << "SeriesMaterial::SeriesMaterial -- failed to allocate strain array\n";
      exit(-1);
    }
    
    stress = new double [numMaterials];
    if (stress == 0) {
      opserr << "SeriesMaterial::SeriesMaterial -- failed to allocate stress array\n";
      exit(-1);
    }

    flex = new double [numMaterials];
    if (flex == 0) {
      opserr << "SeriesMaterial::SeriesMaterial -- failed to allocate flex array\n";
      exit(-1);
    }

    for (i = 0; i < numMaterials; i++) {
      strain[i] = 0.0;
      stress[i] = 0.0;
      flex[i] = 0.0;
    }
}

SeriesMaterial::SeriesMaterial()
:UniaxialMaterial(0,MAT_TAG_SeriesMaterial),
 Tstrain(0.0), Cstrain(0.0), Tstress(0.0), Cstress(0.0),
 Ttangent(0.0), Ctangent(0.0), 
 maxIterations(0), tolerance(0.0),
 stress(0), flex(0), strain(0), initialFlag(false),
 numMaterials(0), theModels(0)
{

}

SeriesMaterial::~SeriesMaterial()
{
    for (int i = 0; i < numMaterials; i++)
		delete theModels[i];

    if (theModels != 0)
		delete [] theModels;

	if (strain)
		delete [] strain;

	if (stress)
		delete [] stress;

	if (flex)
		delete [] flex;

}

int 
SeriesMaterial::setTrialStrain(double newStrain, double strainRate)
{
	// Using the incremental iterative strain
	double dv = newStrain-Tstrain;

	Tstrain = newStrain;

	// Stress increment using tangent at last iteration
	double dq = Ttangent*dv;

	// Update stress 
	Tstress += dq;

	for (int j = 0; j < maxIterations; j++) {

		// Set to zero for integration
		double f = 0.0;
		double vr = 0.0;

		for (int i = 0; i < numMaterials; i++) {

			// Stress unbalance in material i
			double ds = Tstress - stress[i];

			// Strain increment
			double de = flex[i]*ds;

			if (initialFlag == true)
				strain[i] += de;

			// Update material i
			theModels[i]->setTrialStrain(strain[i]);

			// Get updated stress from material i
			stress[i] = theModels[i]->getStress();

			// Get updated flexibility from material i
			flex[i] = theModels[i]->getTangent();
			if (fabs(flex[i]) > 1.0e-12)
				flex[i] = 1.0/flex[i];
			else
				flex[i] = (flex[i] < 0.0) ? -1.0e12 : 1.0e12;

			// Stress unbalance in material i
			ds = Tstress - stress[i];

			// Residual strain in material i
			de = flex[i]*ds;

			// Integrate flexibility ...
			f += flex[i];

			// ... and integrate residual strain
			vr += strain[i] + de;
		}

		// Updated series tangent
		if (fabs(f) > 1.0e-12)
			Ttangent = 1.0/f;
		else
			Ttangent = (f < 0.0) ? -1.0e12 : 1.0e12;

		// Residual deformation
		dv = Tstrain - vr;

		// Stress increment
		dq = Ttangent*dv;

		if (fabs(dq*dv) < tolerance)
			break;
	}

	// Updated stress
	Tstress += dq;
	
	initialFlag = true;

	return 0;
}

double 
SeriesMaterial::getStrain(void)
{
    return Tstrain;
}

double 
SeriesMaterial::getStress(void)
{
	return Tstress;
}

double 
SeriesMaterial::getTangent(void)
{
	return Ttangent;
}


double 
SeriesMaterial::getInitialTangent(void)
{
  double kf = 0.0;
  double k  = 0.0;

  if (numMaterials != 0)
    kf = theModels[0]->getInitialTangent();
  
  for (int i=1; i<numMaterials; i++) {
    
    k = theModels[i]->getInitialTangent();
    if ((kf + k) != 0.0)
      kf = kf*k/(kf+k);
    else
      return 0.0;
  }

  return kf;
}


int 
SeriesMaterial::commitState(void)
{
	int err = 0;

	Cstrain = Tstrain;
	Cstress = Tstress;
	Ctangent = Ttangent;

	for (int i = 0; i < numMaterials; i++)
		err += theModels[i]->commitState();

	// Commented out for the same reason it was taken
	// out of the NLBeamColumn commitState() -- MHS
	//initialFlag = false;

	return err;
}

int 
SeriesMaterial::revertToLastCommit(void)
{
	int err = 0;

	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = Ctangent;
    
	for (int i = 0; i < numMaterials; i++) {
		err += theModels[i]->revertToLastCommit();

		strain[i] = theModels[i]->getStrain();
		stress[i] = theModels[i]->getStress();
		flex[i] = theModels[i]->getTangent();
		
		if (fabs(flex[i]) > 1.0e-12)
			flex[i] = 1.0/flex[i];
		else
			flex[i] = (flex[i] < 0.0) ? -1.0e12 : 1.0e12;
	}

	initialFlag = false;

	return err;
}


int 
SeriesMaterial::revertToStart(void)
{
	int err = 0;

	Cstrain = 0.0;
	Cstress = 0.0;
	Ctangent = 0.0;

	for (int i = 0; i < numMaterials; i++) {
		err += theModels[i]->revertToLastCommit();

		strain[i] = 0.0;
		stress[i] = 0.0;
		flex[i] = 0.0;
	}

    return err;    
}



UniaxialMaterial *
SeriesMaterial::getCopy(void)
{
    SeriesMaterial *theCopy = new 
		SeriesMaterial(this->getTag(), numMaterials, theModels,
			maxIterations, tolerance);

    theCopy->Cstrain = Cstrain;
    theCopy->Cstress = Cstress;
    theCopy->Ctangent = Ctangent;
    theCopy->initialFlag = initialFlag;
    
    for (int i = 0; i < numMaterials; i++) {
      theCopy->strain[i] = strain[i];
      theCopy->stress[i] = stress[i];
      theCopy->flex[i] = flex[i];
    }
    
    return theCopy;
}


int 
SeriesMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  int dataTag = this->getDbTag();
  
  static Vector data(5);
  
  data(0) = this->getTag();
  data(1) = numMaterials;
  data(2) = (initialFlag) ? 1.0 : 0.0;
  data(3) = maxIterations;
  data(4) = tolerance;
  
  res = theChannel.sendVector(dataTag, cTag, data);
  if (res < 0) {
    opserr << "SeriesMaterial::sendSelf -- failed to send data Vector\n";
    return res;
  }
  
  ID classTags(2*numMaterials);
  
  int i;
  for (i = 0; i < numMaterials; i++) {
    classTags(i) = theModels[i]->getClassTag();
    
    int dbTag = theModels[i]->getDbTag();
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	theModels[i]->setDbTag(dbTag);
    }
    
    classTags(i+numMaterials) = dbTag;
  }
  
  res = theChannel.sendID(dataTag, cTag, classTags);
  if (res < 0) {
    opserr << "SeriesMaterial::sendSelf -- failed to send classTags ID\n";
    return res;
  }
  
  for (i = 0; i < numMaterials; i++) {
    res = theModels[i]->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "SeriesMaterial::sendSelf -- failed to send UniaxialMaterial: " << i << endln;
      return res;
    }
  }
  
  return res;
}

int 
SeriesMaterial::recvSelf(int cTag, Channel &theChannel, 
				FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  static Vector data(5);

  res = theChannel.recvVector(dataTag, cTag, data);
  if (res < 0) {
    opserr << "SeriesMaterial::recvSelf -- failed to receive data Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  initialFlag = (data(2) == 1.0) ? true : false;
  maxIterations = (int)data(3);
  tolerance = data(4);

  // if number of materials != new number we must alolocate space for data
  int i;
  if (numMaterials != (int)data(1)) {
    
    // free up old memory if allocated
    if (theModels != 0) {
      for (i = 0; i < numMaterials; i++)
	if (theModels[i] != 0)
	  delete theModels[i];
      delete [] theModels;
    }

    if (strain != 0)
      delete [] strain;

    if (stress != 0)
      delete [] stress;

    if (flex != 0)
      delete [] flex;

    // allocate new memory for data
    numMaterials = (int)data(1);
    theModels = new UniaxialMaterial *[numMaterials];
    if (theModels == 0) {
      opserr << "SeriesMaterial::recvSelf -- failed to allocate UniaxialMaterial array\n";
      return -1;
    }

    for (i = 0; i < numMaterials; i++)
      theModels[i] = 0;

    strain = new double [numMaterials];
    if (strain == 0) {
      opserr << "SeriesMaterial::recvSelf -- failed to allocate strain array\n";
      return -1;
    }

    stress = new double [numMaterials];
    if (stress == 0) {
      opserr << "SeriesMaterial::recvSelf -- failed to allocate stress array\n";
      return -1;
    }

    flex = new double [numMaterials];
    if (flex== 0) {
      opserr << "SeriesMaterial::recvSelf -- failed to allocate flex array\n";
      return -1;
    }
  }

  ID classTags(2*numMaterials);
  res = theChannel.recvID(dataTag, cTag, classTags);
  if (res < 0) {
    opserr << "SeriesMaterial::recvSelf -- failed to receive classTags ID\n";
    return res;
  }

  for (i = 0; i < numMaterials; i++) {
    int matClassTag = classTags(i);

    if (theModels[i] == 0)
      theModels[i] = theBroker.getNewUniaxialMaterial(matClassTag);
    
    else if (theModels[i]->getClassTag() != matClassTag) {
      delete theModels[i];
      theModels[i] = theBroker.getNewUniaxialMaterial(matClassTag);
    }
    
    if (theModels[i] == 0) {
      opserr << "SeriesMaterial::recvSelf -- failed to get a newUniaxialMaterial\n";
      return -1;
    }
    
    theModels[i]->setDbTag(classTags(i+numMaterials));
    res = theModels[i]->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "SeriesMaterial::recvSelf -- failed to receive UniaxialMaterial: " << i << endln;
      return res;
    }
  }
  
  return res;
}

void 
SeriesMaterial::Print(OPS_Stream &s, int flag)
{
    s << "\nSeriesMaterial, tag: " << this->getTag() << endln;
    s << "\tUniaxial Componenets" << endln;
    for (int i = 0; i < numMaterials; i++)
		s << "\t\tUniaxial Material, tag: " << theModels[i]->getTag() << endln;
}

Response*
SeriesMaterial::setResponse(const char **argv, int argc,
			    Information &info)
{
  // See if the response is one of the defaults
  Response *res = UniaxialMaterial::setResponse(argv, argc, info);
  if (res != 0)
    return res;

  if (strcmp(argv[0],"strains") == 0)
    return new MaterialResponse(this, 100, Vector(numMaterials));

  else if (strcmp(argv[0],"material") == 0 ||
	   strcmp(argv[0],"component") == 0) {
    if (argc > 1) {
      int matNum = atoi(argv[1]) - 1;
      if (matNum >= 0 && matNum < numMaterials)
	return theModels[matNum]->setResponse(&argv[2], argc-2, info);
      else
	return 0;
    }
    else
      return 0;
  }
  
  else
    return this->UniaxialMaterial::setResponse(argv, argc, info);
}

int
SeriesMaterial::getResponse(int responseID, Information &info)
{
  Vector strains(strain, numMaterials);

  switch (responseID) {
  case 100:
    return info.setVector(strains);

  default:
    return this->UniaxialMaterial::getResponse(responseID, info);
  }
}
