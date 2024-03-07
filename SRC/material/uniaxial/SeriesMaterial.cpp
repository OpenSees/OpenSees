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
                                                                        
// $Revision: 1.9 $
// $Date: 2007-02-02 01:19:30 $
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
#include <math.h>
#include <float.h>

#include <string.h>
#include <elementAPI.h>

#include <OPS_Globals.h>

void *
OPS_SeriesMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int maxIter = 1;
  double tol = 1e-10;
  
  // Read optional args first
  int numOptionalArgs = 0;  
  int numArgs = OPS_GetNumRemainingInputArgs();
  while (OPS_GetNumRemainingInputArgs() > 0) {
    std::string type = OPS_GetString();
    if (type == "-iter") {
      numOptionalArgs++;
      int numData = 1;
      if (OPS_GetNumRemainingInputArgs() > 1) {
	if (OPS_GetIntInput(&numData, &maxIter) < 0) {
	  opserr << "WARNING: failed to get maxIter" << endln;
	  return 0;
	}
	numOptionalArgs++;
	if (OPS_GetDoubleInput(&numData, &tol) < 0) {
	  opserr << "WARNING: failed to get tol" << endln;
	  return 0;
	}
	numOptionalArgs++;	  
      }
    }
  }

  if (numArgs > 0) {
    OPS_ResetCurrentInputArg(-numArgs);
  }
  numArgs = numArgs - numOptionalArgs;
    
  if (numArgs < 2) {
    opserr << "Invalid #args,  want: uniaxialMaterial Series $tag $tag1 $tag2 ... $tagN <-iter maxIter tol>" << endln;
    return 0;
  }

  int *iData = new int[numArgs];
  UniaxialMaterial **theMats = new UniaxialMaterial *[numArgs-1];
    
  if (OPS_GetIntInput(&numArgs, iData) != 0) {
    opserr << "WARNING invalid data for uniaxialMaterial Series" << endln;
    return 0;
  }

  for (int i=1; i<numArgs; i++) {
    UniaxialMaterial *theMat = OPS_GetUniaxialMaterial(iData[i]);
    if (theMat == 0) {
      opserr << "WARNING no existing material with tag " << iData[i] 
	     << " for uniaxialMaterial Series " << iData[0] << endln;
      delete [] iData;
      delete [] theMats;
      return 0;
    }
    theMats[i-1] = theMat;
  }

  // Parsing was successful, allocate the material
  theMaterial = new SeriesMaterial(iData[0], numArgs-1, theMats, maxIter, tol);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Series" << endln;
    return 0;
  }
  
  delete [] iData;
  delete [] theMats;
  
  return theMaterial;
}

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
      opserr << "SeriesMaterial::SeriesMaterial -- failed to allocate material array" << endln;
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
      opserr << "SeriesMaterial::SeriesMaterial -- failed to allocate strain array" << endln;
      exit(-1);
    }
    
    stress = new double [numMaterials];
    if (stress == 0) {
      opserr << "SeriesMaterial::SeriesMaterial -- failed to allocate stress array" << endln;
      exit(-1);
    }

    flex = new double [numMaterials];
    if (flex == 0) {
      opserr << "SeriesMaterial::SeriesMaterial -- failed to allocate flex array" << endln;
      exit(-1);
    }

    for (i = 0; i < numMaterials; i++) {
      strain[i] = 0.0;
      stress[i] = 0.0;
      flex[i] = 0.0;
    }

    Ttangent = this->getInitialTangent();
    Ctangent = Ttangent;
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


	if (fabs(dv) < DBL_EPSILON)
	  return 0;

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

    theCopy->Tstrain = Tstrain;
    theCopy->Cstrain = Cstrain;
    theCopy->Tstress = Tstress;
    theCopy->Cstress = Cstress;    
    theCopy->Ttangent = Ttangent;
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
  int dataTag = this->getDbTag();
  
  static Vector data(5);
  
  data(0) = this->getTag();
  data(1) = numMaterials;
  data(2) = (initialFlag) ? 1.0 : 0.0;
  data(3) = maxIterations;
  data(4) = tolerance;
  
  if (theChannel.sendVector(dataTag, cTag, data) < 0) {
    opserr << "SeriesMaterial::sendSelf -- failed to send data Vector" << endln;
    return -1;
  }
  
  ID classTags(2*numMaterials);
  
  for (int i = 0; i < numMaterials; i++) {
    classTags(i) = theModels[i]->getClassTag();
    
    int dbTag = theModels[i]->getDbTag();
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	theModels[i]->setDbTag(dbTag);
    }
    
    classTags(i+numMaterials) = dbTag;
  }
  
  if (theChannel.sendID(dataTag, cTag, classTags) < 0) {
    opserr << "SeriesMaterial::sendSelf -- failed to send classTags ID" << endln;
    return -2;
  }

  // Note: will not clash with Vector of length 5 above
  Vector stateData(3*numMaterials);
  for (int i = 0; i < numMaterials; i++) {
    stateData(i               ) = strain[i];
    stateData(i+  numMaterials) = stress[i];
    stateData(i+2*numMaterials) = flex[i];    
  }
  
  if (theChannel.sendVector(dataTag, cTag, stateData) < 0) {
    opserr << "SeriesMaterial::sendSelf -- failed to send stateData Vector" << endln;
    return -3;
  }
  
  for (int i = 0; i < numMaterials; i++) {
    if (theModels[i]->sendSelf(cTag, theChannel) < 0) {
      opserr << "SeriesMaterial::sendSelf -- failed to send UniaxialMaterial: " << i << endln;
      return -4;
    }
  }
  
  return 0;
}

int 
SeriesMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dataTag = this->getDbTag();

  static Vector data(5);

  if (theChannel.recvVector(dataTag, cTag, data) < 0) {
    opserr << "SeriesMaterial::recvSelf -- failed to receive data Vector" << endln;
    return -1;
  }

  this->setTag((int)data(0));
  initialFlag = (data(2) == 1.0) ? true : false;
  maxIterations = (int)data(3);
  tolerance = data(4);

  // if number of materials != new number we must alolocate space for data
  if (numMaterials != (int)data(1)) {
    
    // free up old memory if allocated
    if (theModels != 0) {
      for (int i = 0; i < numMaterials; i++)
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
      opserr << "SeriesMaterial::recvSelf -- failed to allocate UniaxialMaterial array" << endln;
      return -2;
    }

    for (int i = 0; i < numMaterials; i++)
      theModels[i] = 0;

    strain = new double [numMaterials];
    if (strain == 0) {
      opserr << "SeriesMaterial::recvSelf -- failed to allocate strain array" << endln;
      return -3;
    }

    stress = new double [numMaterials];
    if (stress == 0) {
      opserr << "SeriesMaterial::recvSelf -- failed to allocate stress array" << endln;
      return -3;
    }

    flex = new double [numMaterials];
    if (flex == 0) {
      opserr << "SeriesMaterial::recvSelf -- failed to allocate flex array" << endln;
      return -3;
    }
  }

  ID classTags(2*numMaterials);
  if (theChannel.recvID(dataTag, cTag, classTags) < 0) {
    opserr << "SeriesMaterial::recvSelf -- failed to receive classTags ID" << endln;
    return -4;
  }

  Vector stateData(3*numMaterials);
  if (theChannel.recvVector(dataTag, cTag, stateData) < 0) {
    opserr << "SeriesMaterial::recvSelf -- failed to receive stateData Vector" << endln;
    return -5;
  }
  for (int i = 0; i < numMaterials; i++) {
    strain[i] = stateData(i               );
    stress[i] = stateData(i+  numMaterials);
    flex[i]   = stateData(i+2*numMaterials);
  }
  
  for (int i = 0; i < numMaterials; i++) {
    int matClassTag = classTags(i);

    if (theModels[i] == 0)
      theModels[i] = theBroker.getNewUniaxialMaterial(matClassTag);
    
    else if (theModels[i]->getClassTag() != matClassTag) {
      delete theModels[i];
      theModels[i] = theBroker.getNewUniaxialMaterial(matClassTag);
    }
    
    if (theModels[i] == 0) {
      opserr << "SeriesMaterial::recvSelf -- failed to get a newUniaxialMaterial" << endln;
      return -6;
    }
    
    theModels[i]->setDbTag(classTags(i+numMaterials));
    if (theModels[i]->recvSelf(cTag, theChannel, theBroker) < 0) {
      opserr << "SeriesMaterial::recvSelf -- failed to receive UniaxialMaterial: " << i << endln;
      return -7;
    }
  }
  
  return 0;
}

void 
SeriesMaterial::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "\nSeriesMaterial, tag: " << this->getTag() << endln;
        s << "\tUniaxial Components" << endln;
        for (int i = 0; i < numMaterials; i++)
            s << "\t\tUniaxial Material, tag: " << theModels[i]->getTag() << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"SeriesMaterial\", ";
        s << "\"materials\": [";
        for (int i = 0; i < numMaterials - 1; i++)
            s << "\"" << theModels[i]->getTag() << "\", ";
        s << "\"" << theModels[numMaterials - 1]->getTag() << "\"]}";
    }
}

Response*
SeriesMaterial::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{

  Response *theResponse = 0;

  if (strcmp(argv[0],"strains") == 0) {
    for (int i=0; i<numMaterials; i++) {
      theOutput.tag("UniaxialMaterialOutput");
      theOutput.attr("matType", this->getClassType());
      theOutput.attr("matTag", this->getTag());
      theOutput.tag("ResponseType", "eps11");
      theOutput.endTag();
    }

    theResponse =  new MaterialResponse(this, 100, Vector(numMaterials));
  }
  else if (strcmp(argv[0],"material") == 0 ||
	   strcmp(argv[0],"component") == 0) {
    if (argc > 1) {
      int matNum = atoi(argv[1]) - 1;
      if (matNum >= 0 && matNum < numMaterials)
	theResponse =  theModels[matNum]->setResponse(&argv[2], argc-2, theOutput);
    }
  }

  if (theResponse != 0)
    return theResponse;
  else
    return UniaxialMaterial::setResponse(argv, argc, theOutput);
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
