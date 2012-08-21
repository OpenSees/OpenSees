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
#define OPS_Export 

static int numInitialStateAnalysisWrapperMaterials = 0;

OPS_Export void *
OPS_NewInitialStateAnalysisWrapperMaterial(void)
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
	NDMaterial *theMainMaterial = OPS_GetNDMaterial(matID);
	if (theMainMaterial == 0) {
		opserr << "WARNING: For InitialStateAnalysisWrapper " << iData[0] << endln;
		opserr << "Material: " << matID << "not found\n";
		return 0;
	}

	// parsing was sucessfull, allocate the material
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
	// send int data
	static ID idData(4);
	idData(0) = this->getTag();
	idData(1) = theMainMaterial->getClassTag();
	idData(2) = mDIM;
	int matDbTag = theMainMaterial->getDbTag();
	// check if main material has a database tag before sending to database channel
	if (matDbTag == 0) {
		matDbTag = theChannel.getDbTag();
		if (matDbTag != 0) {
			theMainMaterial->setDbTag(matDbTag);
		}
	}
	idData(3) = matDbTag;

	int res = 0;
	res += theChannel.sendID(this->getDbTag(), commitTag, idData);
	if (res < 0) {
		opserr << "WARNING: InitialStateAnalysisWrapper - " << this->getTag() << " - failed to send ID data to channel" << endln;
		return res;
    }

	// send double data
	int vecSize = 3*mDIM - 3;
	Vector data(2*vecSize);

	int i;
	int cnt = 0;
	for (i = 0; i < vecSize; i++) data(cnt+i) = mEpsilon_o(i);
	cnt = cnt+i+1;
	for (i = 0; i < vecSize; i++) data(cnt+i) = mStrain(i);

	res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: InitialStateAnalysisWrapper - " << this->getTag() << " - failed to send vector data to channel" << endln;
		return res;
    }

	// instruct the main material to send itself
	res += theMainMaterial->sendSelf(commitTag, theChannel);
	if (res < 0) {
		opserr << "WARNING: InitialStateAnalysisWrapper - " << this->getTag() << " - failed to send to main material" << endln;
		return res;
    }

	return res;
}

int
InitialStateAnalysisWrapper::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	
    // receive int data
	static ID idData(2);
	res += theChannel.recvID(this->getDbTag(), commitTag, idData);
	if (res < 0) {
		opserr << "WARNING InitialStateAnalysisWrapper - " << this->getTag() << " - could not receive ID data" << endln;
		return res;
	}

	// set int member variables
	this->setTag((int)idData(0));
	mDIM = idData(2);

	// check if material object exists and that it is the right type
	int matClassTag = idData(1);
	int matDbTag = idData(3);
	if ((theMainMaterial == 0) || (theMainMaterial->getClassTag() != matClassTag)) {
		
		// if old, delete
		if (theMainMaterial != 0) {
			delete theMainMaterial;
		}

		// create new material object
		theMainMaterial = theBroker.getNewNDMaterial(matClassTag);
		if (theMainMaterial == 0) {
			opserr << "InitialStateAnalysisWrapper - " << this->getTag() << " - Broker could not create nDMaterial - " << matClassTag << endln;
			exit(-1);
		}
	}
	
	// set material dBtag and receive the material
	theMainMaterial->setDbTag(matDbTag);
	res += theMainMaterial->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		opserr << "WARNING InitialStateAnalysisWrapper::recvSelf() - " << this->getTag() << " failed to receive its Material\n";
		return res;
	}

	// receive double data
	int vecSize = 3*mDIM - 3;
	Vector data(2*vecSize);
	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: InitialStateAnalysisWrapper - " << this->getTag() << " - could not receive vector data" << endln;
		return res;
	}

	// set vector member variables
	int i;
	int cnt = 0;
	for (i = 0; i < vecSize; i++) mEpsilon_o(i) = data(cnt+i);
	cnt = cnt+i+1;
	for (i = 0; i < vecSize; i++) mStrain(i) = data(cnt+i);

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
