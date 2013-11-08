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
                                                                        
// $Revision: 1.14 $
// $Date: 2008-08-26 16:45:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/SectionForceDeformation.cpp,v $
                                                                        
                                                                        
// Written: MHS 
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for SectionForceDeformation.
//
// What: "@(#) SectionForceDeformation.C, revA"

#include <SectionForceDeformation.h>
#include <Information.h>
#include <Matrix.h>
#include <Vector.h>
#include <MaterialResponse.h>

#include <string.h>

double invert2by2Matrix(const Matrix &a, Matrix &b);
double invert3by3Matrix(const Matrix &a, Matrix &b);
void invertMatrix(int n, const Matrix &a, Matrix &b);


#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>

static MapOfTaggedObjects theSectionForceDeformationObjects;

bool OPS_addSectionForceDeformation(SectionForceDeformation *newComponent) {
  return theSectionForceDeformationObjects.addComponent(newComponent);
}

SectionForceDeformation *OPS_getSectionForceDeformation(int tag) {

  TaggedObject *theResult = theSectionForceDeformationObjects.getComponentPtr(tag);
  if (theResult == 0) {
    opserr << "SectionForceDeformation *getSectionForceDeformation(int tag) - none found with tag: " << tag << endln;
    return 0;
  }
  SectionForceDeformation *theMat = (SectionForceDeformation *)theResult;

  return theMat;
}

void OPS_clearAllSectionForceDeformation(void) {
  theSectionForceDeformationObjects.clearAll();
}


SectionForceDeformation::SectionForceDeformation(int tag, int classTag)
  :Material(tag,classTag), fDefault(0), sDefault(0)
{

}

SectionForceDeformation::~SectionForceDeformation()
{
  if (fDefault != 0)
    delete fDefault;
  if (sDefault != 0)
    delete sDefault;
}

const Matrix&
SectionForceDeformation::getSectionFlexibility ()
{
  int order = this->getOrder();
  
  if (fDefault == 0) {		
    fDefault = new Matrix(order,order);
    if (fDefault == 0) {
      opserr << "SectionForceDeformation::getSectionFlexibility -- failed to allocate flexibility matrix\n";
      exit(-1);
    }
  }

  const Matrix &k = this->getSectionTangent();
  
  switch(order) {
  case 1:
    if (k(0,0) != 0.0)
      (*fDefault)(0,0) = 1.0/k(0,0);
    break;
  case 2:
    invert2by2Matrix(k,*fDefault);
    break;
  case 3:
    invert3by3Matrix(k,*fDefault);
    break;
  default:
    invertMatrix(order,k,*fDefault);
    break;
  }

  return *fDefault;
}

const Matrix&
SectionForceDeformation::getInitialFlexibility ()
{
  int order = this->getOrder();
  
  if (fDefault == 0) {		
    fDefault = new Matrix(order,order);
    if (fDefault == 0) {
      opserr << "SectionForceDeformation::getInitialFlexibility -- failed to allocate flexibility matrix\n";
      exit(-1);
    }
  }
  
  const Matrix &k = this->getInitialTangent();
  
  switch(order) {
  case 1:
    if (k(0,0) != 0.0)
      (*fDefault)(0,0) = 1.0/k(0,0);
    break;
  case 2:
    invert2by2Matrix(k,*fDefault);
    break;
  case 3:
    invert3by3Matrix(k,*fDefault);
    break;
  default:
    invertMatrix(order,k,*fDefault);
    break;
  }
  
  return *fDefault;
}

double 
SectionForceDeformation::getRho(void) 
{
  return 0.0 ;
}

Response*
SectionForceDeformation::setResponse(const char **argv, int argc,
				     OPS_Stream &output)
{
  const ID &type = this->getType();
  int typeSize = this->getOrder();
  
  Response *theResponse =0;

  output.tag("SectionOutput");
  output.attr("secType", this->getClassType());
  output.attr("secTag", this->getTag());

  // deformations
  if (strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"deformation") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","kappaY");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","gammaZ");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","theta");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 1, this->getSectionDeformation());
  
  // forces
  } else if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","My");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","Vz");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","T");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 2, this->getStressResultant());
  
  // force and deformation
  } else if (strcmp(argv[0],"forceAndDeformation") == 0) { 
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","kappaY");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","gammaZ");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","theta");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    for (int j=0; j<typeSize; j++) {
      int code = type(j);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","My");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","Vz");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","T");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    
    theResponse =  new MaterialResponse(this, 4, Vector(2*this->getOrder()));
  }
  
  else if (strcmp(argv[0],"stiffness") == 0) {
    theResponse =  new MaterialResponse(this, 12, this->getSectionTangent());
  }

  else if (strcmp(argv[0],"flexibility") == 0) {
    theResponse =  new MaterialResponse(this, 13, this->getSectionFlexibility());
  }
  

  output.endTag(); // SectionOutput
  return theResponse;
}

int 
SectionForceDeformation::getResponse(int responseID, Information &secInfo)
{
  switch (responseID) {
  case 1:
    return secInfo.setVector(this->getSectionDeformation());
    
  case 2:
    return secInfo.setVector(this->getStressResultant());
    
  case 4: {
    Vector &theVec = *(secInfo.theVector);
    const Vector &e = this->getSectionDeformation();
    const Vector &s = this->getStressResultant();
    int order = this->getOrder();
    for (int i = 0; i < order; i++) {
      theVec(i) = e(i);
      theVec(i+order) = s(i);
    }
    
    return secInfo.setVector(theVec);
  }

  case 12:
    return secInfo.setMatrix(this->getSectionTangent());

  case 13:
    return secInfo.setMatrix(this->getSectionFlexibility());

  default:
    return -1;
  }
}

int 
SectionForceDeformation::getResponseSensitivity(int responseID, int gradIndex,
						Information &secInfo)
{
  Vector &theVec = *(secInfo.theVector);

  switch (responseID) {
  case 1:
    theVec = this->getSectionDeformationSensitivity(gradIndex);
    return secInfo.setVector(theVec);
    
  case 2: {
    const Matrix &ks = this->getSectionTangent();
    const Vector &dedh = this->getSectionDeformationSensitivity(gradIndex);
    const Vector &dsdh = this->getStressResultantSensitivity(gradIndex, true);
    theVec.addMatrixVector(0.0, ks, dedh, 1.0);
    theVec.addVector(1.0, dsdh, 1.0);
    return secInfo.setVector(theVec);
  }

  default:
    return -1;
  }
}

// AddingSensitivity:BEGIN ////////////////////////////////////////
const Vector &
SectionForceDeformation::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  if (sDefault == 0)
    sDefault = new Vector (this->getOrder());
  return *sDefault;
}

const Vector &
SectionForceDeformation::getSectionDeformationSensitivity(int gradIndex)
{
  if (sDefault == 0)
    sDefault = new Vector (this->getOrder());
  return *sDefault;
}

const Matrix &
SectionForceDeformation::getSectionTangentSensitivity(int gradIndex)
{
  int order = this->getOrder();
  
  if (fDefault == 0) {		
    fDefault = new Matrix(order,order);
    if (fDefault == 0) {
      opserr << "SectionForceDeformation::getSectionTangentSensitivity -- failed to allocate matrix\n";
      exit(-1);
    }
  }

  fDefault->Zero();

  return *fDefault;
}

const Matrix &
SectionForceDeformation::getInitialTangentSensitivity(int gradIndex)
{
  int order = this->getOrder();
  
  if (fDefault == 0) {		
    fDefault = new Matrix(order,order);
    if (fDefault == 0) {
      opserr << "SectionForceDeformation::getInitialTangentSensitivity -- failed to allocate matrix\n";
      exit(-1);
    }
  }

  fDefault->Zero();

  return *fDefault;
}

const Matrix&
SectionForceDeformation::getSectionFlexibilitySensitivity(int gradIndex)
{
  int order = this->getOrder();
  
  if (fDefault == 0) {		
    fDefault = new Matrix(order,order);
    if (fDefault == 0) {
      opserr << "SectionForceDeformation::getSectionFlexibilitySensitivity -- failed to allocate matrix\n";
      exit(-1);
    }
  }

  const Matrix &dksdh = this->getSectionTangentSensitivity(gradIndex);
  
  const Matrix &fs = this->getSectionFlexibility();

  *fDefault = (fs * dksdh * fs) * -1;

  return *fDefault;
}

const Matrix&
SectionForceDeformation::getInitialFlexibilitySensitivity(int gradIndex)
{
  int order = this->getOrder();
  
  if (fDefault == 0) {		
    fDefault = new Matrix(order,order);
    if (fDefault == 0) {
      opserr << "SectionForceDeformation::getInitialFlexibilitySensitivity -- failed to allocate matrix\n";
      exit(-1);
    }
  }
  
  const Matrix &dksdh = this->getInitialTangentSensitivity(gradIndex);
  
  const Matrix &fs = this->getInitialFlexibility();

  *fDefault = (fs * dksdh * fs) * -1;
  
  return *fDefault;
}

double
SectionForceDeformation::getRhoSensitivity(int gradIndex)
{
  return 0.0;
}

int
SectionForceDeformation::commitSensitivity(const Vector& defSens,
					   int gradIndex, int numGrads)
{
  return -1;
}
// AddingSensitivity:END ///////////////////////////////////////////

//--- Adding Thermal Functions:[BEGIN]   by UoE OpenSees Group ----//
int
SectionForceDeformation::setTrialSectionDeformation (const Vector&) //JZ
{
  opserr << "SectionForceDeformation::setTrialSectionDeformation(strain) - should not be called\n";
  return -1;
}

int
SectionForceDeformation::setTrialSectionDeformation(const Vector& nouse, const Vector &data) //JZ
{
  opserr << "SectionForceDeformation::setTrialSectionDeformationTemperature (strain, tData) - should not be called\n";
  return -1;
}

//static Vector errRes(3);

const Vector &
SectionForceDeformation::getTemperatureStress(const Vector &tData) //PK
{
  opserr << "SectionForceDeformation::getTemperatureStress(double *dataMixed) - should not be called\n";
  return this->getStressResultant();
}
//--- Adding Thermal Functions:[END]   by UoE OpenSees Group ----//
