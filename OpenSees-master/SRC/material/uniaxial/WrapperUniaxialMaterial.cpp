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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/WrapperUniaxialMaterial.cpp,v $

// Written: fmk                                                                         

#include <WrapperUniaxialMaterial.h>
#include <string.h>

extern modelState theModelState;

WrapperUniaxialMaterial::WrapperUniaxialMaterial(const char *name, matObject *theMat_)
  :UniaxialMaterial(theMat_->tag,MAT_TAG_WrapperUniaxialMaterial),
  funcName(0),
  theMat(theMat_),
  strain(0.0), stress(0.0), tangent(0.0), initTangent(0.0)
{
  /*opserr << "WrapperMaterial::WrapperMaterial() " << theMat->tag << endln; */

  funcName = new char[strlen(name)+1];
  if (funcName != 0)
    strcpy(funcName, name);

  int isw = ISW_FORM_TANG_AND_RESID;
  int error = 0;
  theMat->matFunctPtr(theMat, &theModelState, &strain, &initTangent, &stress, &isw, &error);

  tangent=initTangent;
}

WrapperUniaxialMaterial::WrapperUniaxialMaterial()
  :UniaxialMaterial(0, MAT_TAG_WrapperUniaxialMaterial),
  funcName(0),
  theMat(0),
  strain(0.0), stress(0.0), tangent(0.0), initTangent(0.0)
{
  
}

// destructor
WrapperUniaxialMaterial::~WrapperUniaxialMaterial()
{

  /*opserr << "WrapperUniaxialMaterial::~WrapperUniaxialMaterial()\n";*/

  if (funcName != 0)
    delete [] funcName;

  if (theMat->theParam != 0)
    delete [] theMat->theParam;

  if (theMat->cState != 0)
    delete [] theMat->cState;

  if (theMat->tState != 0)
    delete [] theMat->tState;

  delete theMat;
}

int 
WrapperUniaxialMaterial::setTrialStrain (double theStrain, double strainRate)
{
  int isw = ISW_FORM_TANG_AND_RESID;
  int error = 0;
  strain = theStrain;
  theMat->matFunctPtr(theMat, &theModelState, &strain, &tangent, &stress, &isw, &error);

  return error;
}

int 
WrapperUniaxialMaterial::setTrial (double theStrain, double &theStress, double &theTangent, double strainRate)
{
  int isw = ISW_FORM_TANG_AND_RESID;
  int error = 0;
  strain = theStrain;
  theMat->matFunctPtr(theMat, &theModelState, &strain, &tangent , &stress, &isw, &error);
  theTangent = tangent;
  theStress = stress;

  return error;
}

double 
WrapperUniaxialMaterial::getStrain (void)
{
  return strain;
}
double 
WrapperUniaxialMaterial::getStress (void)
{
  return stress;
}

double 
WrapperUniaxialMaterial::getTangent (void)
{
  return tangent;
}

double 
WrapperUniaxialMaterial::getInitialTangent(void)
{
  return initTangent;
}

double 
WrapperUniaxialMaterial::getDampTangent (void)
{
  return 0;
}

int 
WrapperUniaxialMaterial::commitState (void)
{
  int isw = ISW_COMMIT;
  int error = 0;
  theMat->matFunctPtr(theMat, &theModelState, &strain, &tangent, &stress, &isw, &error);

  return error;
}

int 
WrapperUniaxialMaterial::revertToLastCommit (void)
{
  int isw = ISW_REVERT;
  int error = 0;
  theMat->matFunctPtr(theMat, &theModelState, &strain, &tangent , &stress, &isw, &error);
  return error;
}

int 
WrapperUniaxialMaterial::revertToStart (void)
{
  int isw = ISW_REVERT_TO_START;
  int error = 0;
  theMat->matFunctPtr(theMat, &theModelState, &strain, &tangent , &stress, &isw, &error);
  return error;
}

UniaxialMaterial *
WrapperUniaxialMaterial::getCopy (void) 
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

    WrapperUniaxialMaterial *theResult = new WrapperUniaxialMaterial(funcName, theMatObject);
    return theResult;
}

int 
WrapperUniaxialMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}
int 
WrapperUniaxialMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}
void 
WrapperUniaxialMaterial::Print(OPS_Stream &s, int flag)
{
  s << "WrapperUniaxialMaterial - wrapping function  matTag: " << theMat->tag << endln;
}



