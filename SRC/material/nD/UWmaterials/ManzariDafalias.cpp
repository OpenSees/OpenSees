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
                                                                       
// Created: Pedro Arduino, UW, 11.2011
//
// Description: This file contains the implementation of the ManzariDafalias class.

#include <ManzariDafalias.h>
#include <ManzariDafalias3D.h>
#include <ManzariDafaliasPlaneStrain.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <elementAPI.h>

static int numManzariDafaliasMaterials = 0;

void *
OPS_NewManzariDafaliasMaterial(void)
{
  if (numManzariDafaliasMaterials == 0) {
    numManzariDafaliasMaterials++;
    opserr << "ManzariDafalias nDmaterial - Written: P.Arduino, C.McGann, U.Washington\n";
  }

  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 10) {
    opserr << "Want: nDMaterial ManzariDafalias tag? " << endln;
    return 0;	
  }
  
  int tag;
  double dData[9];

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial ManzariDafalias material  tag" << endln;
    return 0;
  }
  numData = 9;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial ManzariDafalias material  with tag: " << tag << endln;
    return 0;
  }

  theMaterial = new ManzariDafalias(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], 
                                            dData[6], dData[7], dData[8]);
  
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial ManzariDafalias material  with tag: " << tag << endln;
  }

  return theMaterial;
}

// full constructor
ManzariDafalias::ManzariDafalias(int tag, int classTag, double mDen)
  : NDMaterial(tag,ND_TAG_ManzariDafalias)
{
}

   
// null constructor
ManzariDafalias ::ManzariDafalias() 
  : NDMaterial()
{
}

// destructor
ManzariDafalias::~ManzariDafalias()
{
}

NDMaterial*
ManzariDafalias::getCopy(const char *type)
{
	if (strcmp(type,"PlanStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0) {
		ManzariDafaliasPlaneStrain *clone;
		clone = new ManzariDafaliasPlaneStrain(this->getTag(), massDen);
		return clone;
	} else if (strcmp(type,"ThreeDimensional")==0 || strcmp(type,"3D") ==0) {
		ManzariDafalias3D *clone;
     	clone = new ManzariDafalias3D(this->getTag(), massDen);
	 	return clone;
  	} else {
	  	opserr << "ManzariDafalias::getCopy failed to get copy: " << type << endln;
	  	return 0;
  	}
}

int 
ManzariDafalias::commitState(void)
{
	return 0;
}
 
int ManzariDafalias::revertToLastCommit (void)
{
    return 0;
}

int ManzariDafalias::revertToStart(void)
{
	// added: C.McGann, U.Washington for InitialStateAnalysis
	if (ops_InitialStateAnalysis) {
		// do nothing, keep state variables from last step
	} else {
		// normal call for revertToStart (not initialStateAnalysis)
    	this->initialize();
	}

    return 0;
}

NDMaterial*
ManzariDafalias::getCopy (void)
{
	opserr << "ManzariDafalias::getCopy -- subclass responsibility\n"; 
  	exit(-1);
  	return 0;
}

const char*
ManzariDafalias::getType (void) const
{
    opserr << "ManzariDafalias::getType -- subclass responsibility\n";
    exit(-1);
    return 0;
}

int
ManzariDafalias::getOrder (void) const
{
    opserr << "ManzariDafalias::getOrder -- subclass responsibility\n";
    exit(-1);
    return 0;
}

Response*
ManzariDafalias::setResponse (const char **argv, int argc, OPS_Stream &output)
{
	if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->GetState());
	else
		return 0;
}

int
ManzariDafalias::getResponse(int responseID, Information &matInfo)
{
	switch (responseID) {
		case -1:
			return -1;
		case 1:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getStress();
			return 0;
		case 2:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getStrain();
			return 0;
		case 3:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = GetState();
			return 0;
		default:
			return -1;
	}
}

int
ManzariDafalias::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state into a vector object
  static Vector data(8);
  int cnt = 0;

  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "ManzariDafalias::sendSelf - failed to send vector to channel\n";
    return -1;
  }
  

  return 0;
 
}

int 
ManzariDafalias::recvSelf(int commitTag, Channel &theChannel, 
                                         FEM_ObjectBroker &theBroker)    
{
  // recv the vector object from the channel which defines material param and state
  static Vector data(7);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "ManzariDafalias::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag(data(cnt++));

  return 0;

}

void ManzariDafalias::Print(OPS_Stream &s, int flag )
{
  s << "ManzariDafalias" << endln;
}

int
ManzariDafalias::setParameter(const char **argv, int argc, Parameter &param)
{
  	if (argc < 2)
    	return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0],"updateMaterialStage") == 0) {
			return param.addObject(1, this);
		}
	}

    return -1;
}

int
ManzariDafalias::updateParameter(int responseID, Information &info)
{
	// called updateMaterialStage in tcl file
	if (responseID == 1) {
		mElastFlag = info.theInt;
	}
	// called materialState in tcl file
	if (responseID == 5) {
		mElastFlag = info.theDouble;
	}
	
	return 0;
}
