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
                                                                        
// Written: Chris McGann
//          February 2011
                                                                      
// Description: This file contains the implementation for the InitialStateAnalysisWrapper class.
//              This wrapper can be used with any nDMaterial, and enables the use of the 
//              InitialStateAnalysis command for setting initial conditions.

#include <InitialStateAnalysisWrapper.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <ID.h>

#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <elementAPI.h>

static int numInitialStateAnalysisWrapperMaterials = 0;

void * OPS_ADD_RUNTIME_VPV(OPS_InitialStateAnalysisWrapperMaterial)
{
	if (numInitialStateAnalysisWrapperMaterials == 0) {
		numInitialStateAnalysisWrapperMaterials++;
		opserr << "InitialStateAnalysisWrapper nDmaterial - Written: C.McGann, P.Arduino, P.Mackenzie-Helnwein, U.Washington\n";
  	}

  	NDMaterial *theMaterial = 0;

  	int numArgs = OPS_GetNumRemainingInputArgs();

  	if (numArgs < 2) {
    	opserr << "Want: nDMaterial InitialStateAnalysisWrapper tag? nDMatTag? numDim?" << endln;
    	return 0;	
  	}
  
  	int iData[3];

  	int numData = 3;
  	if (OPS_GetInt(&numData, iData) != 0) {
    	opserr << "WARNING invalid integer data: nDMaterial InitialStateAnalysisWrapper with tag: " << iData[0] << endln;
    	return 0;
  	}

	int matID = iData[1];
	NDMaterial *theMainMaterial = OPS_getNDMaterial(matID);
	if (theMainMaterial == 0) {
		opserr << "WARNING: For InitialStateAnalysisWrapper " << iData[0] << endln;
		opserr << "Material: " << matID << "not found\n";
		return 0;
	}

	// parsing was successful, allocate the material
  	theMaterial = new InitialStateAnalysisWrapper(iData[0], *theMainMaterial, iData[2]);
  
  	if (theMaterial == 0) {
    	opserr << "WARNING ran out of memory for nDMaterial InitialStateAnalysisWrapper with tag: " << iData[0] << endln;
  	}

  	return theMaterial;
}

// full constructor
InitialStateAnalysisWrapper::InitialStateAnalysisWrapper(int tag, NDMaterial &mainMat, int ndim)
  : NDMaterial(tag, ND_TAG_InitialStateAnalysisWrapper),
    theMainMaterial(0),
    mEpsilon_o(3*ndim-3),
    mStrain(3*ndim-3)
{
  mDIM = ndim;
  mEpsilon_o.Zero();
  mStrain.Zero();
  
  // get copy of the main material
  if (mDIM == 2) {
    theMainMaterial = mainMat.getCopy("PlaneStrain");
  } else if (mDIM == 3) {
    theMainMaterial = mainMat.getCopy("ThreeDimensional");
  } else {
    opserr << "Incompatible number of dimensions for InitialStateAnalysisWrapper - want 2 or 3" << endln;
  }
}

// null constructor
InitialStateAnalysisWrapper::InitialStateAnalysisWrapper()
  : NDMaterial(0, ND_TAG_InitialStateAnalysisWrapper),
    theMainMaterial(0),
    mEpsilon_o(3),
    mStrain(3)
{
  mEpsilon_o.Zero();
  mStrain.Zero();
  mDIM = 0;
}

// destructor
InitialStateAnalysisWrapper::~InitialStateAnalysisWrapper()
{
	if (theMainMaterial != 0) {
		delete theMainMaterial;
	}
}

// clone material
NDMaterial*
InitialStateAnalysisWrapper::getCopy(const char *type)
{
  return this->getCopy();
}

NDMaterial*
InitialStateAnalysisWrapper::getCopy(void)
{
  // new instance of class
  InitialStateAnalysisWrapper *clone;
  // make copy
  clone = new InitialStateAnalysisWrapper(this->getTag(), *theMainMaterial, mDIM);
  
  return clone;
}

const char*
InitialStateAnalysisWrapper::getType(void) const
{
  return theMainMaterial->getType();
}

int
InitialStateAnalysisWrapper::getOrder(void) const
{
  return theMainMaterial->getOrder();
}

int
InitialStateAnalysisWrapper::commitState(void)
{
	return theMainMaterial->commitState();
}

int
InitialStateAnalysisWrapper::revertToLastCommit(void)
{
	return theMainMaterial->revertToLastCommit();
}

int
InitialStateAnalysisWrapper::revertToStart(void)
{
  // update epsilon_o when InitialStateAnalysis off is called
  if (ops_InitialStateAnalysis) {
    mEpsilon_o += mStrain;
  }
  return theMainMaterial->revertToStart();
}

int
InitialStateAnalysisWrapper::setTrialStrain(const Vector &strain_from_element)
// this function receives the strain from the element and sends strain to material
{
	// add epsilon_o to the element strain
	mStrain = strain_from_element + mEpsilon_o;
	
	// send the sum to the main material
	theMainMaterial->setTrialStrain(mStrain);

	return 0;
}

double
InitialStateAnalysisWrapper::getRho(void)
// this function gets the mass density from the main material
{
	return theMainMaterial->getRho();
}

const Vector&
InitialStateAnalysisWrapper::getStrain()
// this function sends the strain back to the element
{
  return theMainMaterial->getStrain();
}

const Vector&
InitialStateAnalysisWrapper::getStress()
// this function sends the stress back to the element
{
	return theMainMaterial->getStress();
}

const Matrix&
InitialStateAnalysisWrapper::getTangent()
// this function sends the tangent back to the element
{
	return theMainMaterial->getTangent();
}

const Matrix&
InitialStateAnalysisWrapper::getInitialTangent()
// this function sends the initial tangent back to the element
{
	return theMainMaterial->getInitialTangent();
}

int
InitialStateAnalysisWrapper::getMainClassTag()
// this function sends the class tag of the main material
{
	return theMainMaterial->getClassTag();
}

int
InitialStateAnalysisWrapper::sendSelf(int commitTag, Channel &theChannel)
{
  int res;
  int dataTag = this->getDbTag();
  
  static ID data(4);
  data(0) = this->getTag();
  data(1) = theMainMaterial->getClassTag();
  
  int matDbTag = theMainMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0) {
      theMainMaterial->setDbTag(matDbTag);
    }
  }
  data(2) = matDbTag;
  data(3) = mDIM;
  
  res = theChannel.sendID(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING InitialStateAnalysisWrapper::sendSelf() - " << this->getTag() << " failed to send data\n";
    return -1;
  }


  int dim = 3*mDIM-3;
  Vector oData(2*dim);
  for (int i=0; i<dim; i++) {			
    oData(i) = mStrain(i);
    oData(i+dim) = mEpsilon_o(i);
  }
  
  res = theChannel.sendVector(dataTag, commitTag, oData);
  if (res < 0) {
    opserr << "WARNING InitialStateAnalysisWrapper::sendSelf() - " << this->getTag() << " failed to send Initial State\n";
    return -1;
  }
  
  res = theMainMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "WARNING: InitialStateAnalysisWrapper - " << this->getTag() << " - failed to send vector data to channel" << endln;
    return res;
  }
  
  return res;
}

int
InitialStateAnalysisWrapper::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{

  int res;
  int dataTag = this->getDbTag();
  
  
  static ID data(4);
  res = theChannel.recvID(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING InitialStateAnalysisWrapper::recvSelf() - failed to receive Vector\n";
    return -1;
  }
  
  this->setTag(data(0));
  int matClassTag = data(1);
  int matDbTag = data(2);
  
  mDIM = data(3);
  mEpsilon_o.resize(3*mDIM-3);
  mStrain.resize(3*mDIM-3);

  int dim = 3*mDIM-3;
  Vector oData(2*dim);
  res = theChannel.recvVector(dataTag, commitTag, oData);
  if (res < 0) {
    opserr << "WARNING InitialStateAnalysisWrapper::recvSelf() - failed to receive Vector\n";
    return -1;
  }

  for (int i=0; i<dim; i++) {			
    mStrain(i) = oData(i);
    mEpsilon_o(i) = oData(i+dim);
  }
  
  // check if material object exists and that it is the right type
  if ((theMainMaterial == 0) || (theMainMaterial->getClassTag() != matClassTag)) {
    
    // if old, delete
    if (theMainMaterial != 0) {
      delete theMainMaterial;
    }
    
    // create new material object
    theMainMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMainMaterial == 0) {
      opserr << "InitialStateAnalysisWrapper::recvSelf() - " <<
	"Broker could not create nDMaterial of classType: " << matClassTag << endln;
      exit(-1);
    }
  }
  
  // set material dBtag and receive the material
  theMainMaterial->setDbTag(matDbTag);
  res = theMainMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "WARNING InitialStateAnalysisWrapper::recvSelf() - " << this->getTag() << " failed to receive its Material\n";
    return -3;
  }
  
  return res;
}

void
InitialStateAnalysisWrapper::Print(OPS_Stream &s, int flag)
{
	s << "InitialStateAnalysisWrapper Material Tag: " << this-getTag() << endln;
	s << "wrapping the material: \n";
	theMainMaterial->Print(s, flag);
	
	return;
}

int
InitialStateAnalysisWrapper::setParameter(const char **argv, int argc, Parameter &param)
{
	// this allows for the use of updateMaterialStage command for the main material
	if (strcmp(argv[0], "updateMaterialStage") == 0) {
		if (argc < 2) {
			return -1;
		}
		int matTag = atoi(argv[1]);
		if (this->getTag() == matTag) {
			return param.addObject(1,this);
		} else {
			return -1;
		}
	} else if (strcmp(argv[0], "shearModulus") == 0) {
		if (argc < 2) {
			return -1;
		}
		int matTag = atoi(argv[1]);
		if (this->getTag() == matTag) {
			return param.addObject(10,this);
		} else {
			return -1;
		}
	} else if (strcmp(argv[0], "bulkModulus") == 0) {
		if (argc < 2) {
			return -1;
		}
		int matTag = atoi(argv[1]);
		if (this->getTag() == matTag) {
			return param.addObject(11,this);
		} else {
			return -1;
		}
	} else if (strcmp(argv[0], "frictionAngle") == 0) {
		if (argc < 2) {
			return -1;
		}
		int matTag = atoi(argv[1]);
		if (this->getTag() == matTag) {
			return param.addObject(12,this);
		} else {
			return -1;
		}
	} else if (strcmp(argv[0], "cohesion") == 0) {
		if (argc < 2) {
			return -1;
		}
		int matTag = atoi(argv[1]);
		if (this->getTag() == matTag) {
			return param.addObject(13,this);
		} else {
			return -1;
		}
	}

	return -1;
}

int
InitialStateAnalysisWrapper::updateParameter(int responseID, Information &info)
{
	// routes the updateParameter call to the main material
	theMainMaterial->updateParameter(responseID, info);
	
	return 0;
}

Response*
InitialStateAnalysisWrapper::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	return theMainMaterial->setResponse(argv, argc, output);
}

int
InitialStateAnalysisWrapper::getResponse(int responseID, Information &matInfo)
{
	return theMainMaterial->getResponse(responseID, matInfo);
}
