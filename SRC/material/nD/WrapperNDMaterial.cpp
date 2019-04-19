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
                                                                        
// $Revision: 1.2 $
// $Date: 2008-12-19 17:28:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/WrapperNDMaterial.cpp,v $

// Written: fmk                                                                         

#include <WrapperNDMaterial.h>
#include <string.h>
#include <stdlib.h>

extern modelState theModelState;

WrapperNDMaterial::WrapperNDMaterial(const char *name, matObject *theMat_, int mType)
  :NDMaterial(theMat_->tag,ND_TAG_WrapperNDMaterial),
   funcName(0),
   theMat(theMat_), matType(mType), dataSize(0), data(0),
   strain(0), stress(0), tangent(0), initTangent(0)
{
  funcName = new char[strlen(name)+1];
  if (funcName != 0)
    strcpy(funcName, name);

  if (matType == OPS_PLANESTRESS_TYPE || matType == OPS_PLANESTRAIN_TYPE) {
    dataSize = 3;
    data = new double[24];
    strain = new Vector(&data[0],3);
    stress = new Vector(&data[3],3);
    tangent = new Matrix(&data[6],3,3);
    initTangent = new Matrix(&data[15],3,3);
  } else if (matType == OPS_THREEDIMENSIONAL_TYPE) {
    dataSize = 6;
    data = new double[84];
    strain = new Vector(&data[0],6);
    stress = new Vector(&data[6],6);
    tangent = new Matrix(&data[12],6,6);
    initTangent = new Matrix(&data[48],6,6); // GR: corrected from "tangent" to "initTangent"
  } else {
    opserr << "FATAL:WrapperNDMaterial::WrapperNDMaterial - unknown material type: " << matType << endln;
    exit(-1);
  }

  int isw = ISW_FORM_TANG_AND_RESID;
  int error = 0;
  for (int i=0; i<dataSize; i++) {
    data[i]=0.0;
  }

  double *strainData = &data[0];
  double *stressData = &data[dataSize];
  double *initData = &data[2*dataSize + dataSize*dataSize];

  theMat->matFunctPtr(theMat, &theModelState, strainData, initData, stressData, &isw, &error);

  for (int i=0; i<dataSize*dataSize; i++) {
    data[2*dataSize+i]=data[2*dataSize+dataSize*dataSize+i]=0.0;
  }
}

WrapperNDMaterial::WrapperNDMaterial()
  :NDMaterial(0, ND_TAG_WrapperNDMaterial),
   funcName(0),
   theMat(0), matType(0), dataSize(0), data(0),
   strain(0), stress(0), tangent(0), initTangent(0)
{
  
}

// destructor
WrapperNDMaterial::~WrapperNDMaterial()
{

  /*opserr << "WrapperNDMaterial::~WrapperNDMaterial()\n";*/

  if (funcName != 0)
    delete [] funcName;

  if (theMat->theParam != 0)
    delete [] theMat->theParam;

  if (theMat->cState != 0)
    delete [] theMat->cState;

  if (theMat->tState != 0)
    delete [] theMat->tState;

  delete theMat;

  if (data != 0)
    delete [] data;

  if (strain != 0)
    delete strain;

  if (stress != 0)
    delete stress;

  if (tangent != 0)
    delete tangent;

  if (initTangent != 0)
    delete initTangent;
}

int 
WrapperNDMaterial::setTrialStrain (const Vector & theStrain)
{
  int isw = ISW_FORM_TANG_AND_RESID;
  int error = 0;
  *strain = theStrain;

  double *strainData = &data[0];
  double *stressData = &data[dataSize];
  double *tangData = &data[2*dataSize];

  theMat->matFunctPtr(theMat, &theModelState, strainData, tangData, stressData, &isw, &error);

  return error;
}

const Vector & 
WrapperNDMaterial::getStrain (void)
{
  return *strain;
}
const Vector & 
WrapperNDMaterial::getStress (void)
{
  return *stress;
}

const Matrix & 
WrapperNDMaterial::getTangent (void)
{
  return *tangent;
}

const Matrix  & 
WrapperNDMaterial::getInitialTangent(void)
{
  return *initTangent;
}

int 
WrapperNDMaterial::commitState (void)
{
  int isw = ISW_COMMIT;
  int error = 0;

  double *strainData = &data[0];
  double *stressData = &data[dataSize];
  double *tangData = &data[2*dataSize];
  theMat->matFunctPtr(theMat, &theModelState, strainData, tangData, stressData, &isw, &error);

  return error;
}

int 
WrapperNDMaterial::revertToLastCommit (void)
{
  int isw = ISW_REVERT;
  int error = 0;

  double *strainData = &data[0];
  double *stressData = &data[dataSize];
  double *tangData = &data[2*dataSize];
  theMat->matFunctPtr(theMat, &theModelState, strainData, tangData, stressData, &isw, &error);
  return error;
}

int 
WrapperNDMaterial::revertToStart (void)
{
  int isw = ISW_REVERT_TO_START;
  int error = 0;

  double *strainData = &data[0];
  double *stressData = &data[dataSize];
  double *tangData = &data[2*dataSize];
  theMat->matFunctPtr(theMat, &theModelState, strainData, tangData, stressData, &isw, &error);

  return error;
}


const char*
WrapperNDMaterial::getType (void) const
{
  if (matType == OPS_PLANESTRESS_TYPE) 
    return "PlaneStress";
  else if (matType == OPS_PLANESTRAIN_TYPE) 
    return "PlaneStrain";
  else if (matType == OPS_THREEDIMENSIONAL_TYPE) 
    return "ThreeDimensional";

  return "UNKNOWN";
}


NDMaterial *
WrapperNDMaterial::getCopy (void) 
{
    matObject *theMatObject = new matObject;
    theMatObject->tag = theMat->tag;
    theMatObject->nParam = theMat->nParam;
    theMatObject->nState = theMat->nState;

    OPS_AllocateMaterial(theMatObject);

    for (int i=0; i<theMat->nParam; i++)
      theMatObject->theParam[i] = theMat->theParam[i];

    for (int i=0; i<theMat->nState; i++) {
      theMatObject->cState[i] = theMat->cState[i];
      theMatObject->tState[i] = theMat->tState[i];
    }

    theMatObject->matFunctPtr = theMat->matFunctPtr;

    WrapperNDMaterial *theResult = new WrapperNDMaterial(funcName, theMatObject, matType);
    return theResult;
}

NDMaterial *
WrapperNDMaterial::getCopy (const char *code) 
{
  int ok = -1;

  if (matType == OPS_PLANESTRESS_TYPE) {
    if (strcmp(code, "PlaneStress") == 0)
      ok = 0;
  } else if (matType == OPS_PLANESTRAIN_TYPE) {
    if (strcmp(code, "PlaneStrain") == 0)
      ok = 0;
  } else if (matType == OPS_THREEDIMENSIONAL_TYPE) {
    if (strcmp(code, "PlaneStrain") == 0)
      ok = 0;
	if (strcmp(code, "ThreeDimensional") == 0) // GR added 3D material support
      ok = 0;
  }

  if (ok != 0) {
    opserr << "WrapperNDMaterial::getCopy - unknown code type: " << code << endln;
    return 0;
  }


  matObject *theMatObject = new matObject;
  theMatObject->tag = theMat->tag;
  theMatObject->nParam = theMat->nParam;
  theMatObject->nState = theMat->nState;
  theMatObject->matType = theMat->matType; // GR added 3D material support
  
  OPS_AllocateMaterial(theMatObject);
  
  for (int i=0; i<theMat->nParam; i++)
    theMatObject->theParam[i] = theMat->theParam[i];
  
  for (int i=0; i<theMat->nState; i++) {
    theMatObject->cState[i] = theMat->cState[i];
    theMatObject->tState[i] = theMat->tState[i];
  }
  
  theMatObject->matFunctPtr = theMat->matFunctPtr;
  
  WrapperNDMaterial *theResult = new WrapperNDMaterial(funcName, theMatObject, matType);
  return theResult;
}


int 
WrapperNDMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}
int 
WrapperNDMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}
void 
WrapperNDMaterial::Print(OPS_Stream &s, int flag)
{
  s << "WrapperNDMaterial - wrapping function  matTag: " << theMat->tag << endln;
}



