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

#include <elementAPI.h>
#include <DummyStream.h>
#include <Element.h>
#include <Domain.h>

#include <string.h>

#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>
#include <MapOfTaggedObjectsIter.h>

static MapOfTaggedObjects theSectionForceDeformationObjects;

bool OPS_addSectionForceDeformation(SectionForceDeformation *newComponent) {
  return theSectionForceDeformationObjects.addComponent(newComponent);
}

bool OPS_removeSectionForceDeformation(int tag)
{
    TaggedObject* obj = theSectionForceDeformationObjects.removeComponent(tag);
    if (obj != 0) {
	delete obj;
	return true;
    }
    return false;
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

void OPS_printSectionForceDeformation(OPS_Stream &s, int flag) {

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\"sections\": [\n";    
    MapOfTaggedObjectsIter theObjects = theSectionForceDeformationObjects.getIter();
    theObjects.reset();
    TaggedObject *theObject;
    int count = 0;
    int numComponents = theSectionForceDeformationObjects.getNumComponents();
    while ((theObject = theObjects()) != 0) {
      SectionForceDeformation *theSection = (SectionForceDeformation *)theObject;
      theSection->Print(s, flag);
      if (count < numComponents-1)
	s << ",\n";
      count++;
    }
    s << "\n\t\t]";
  }
}

int OPS_sectionLocation()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - sectionLocation eleTag? <secNum?> \n";
	return -1;
    }

    //opserr << "sectionLocation: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 1;
    int tag;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING sectionLocation eleTag? <secNum?> - could not read int input? \n";
	return -1;
    }

    int secNum = 0;
    if (OPS_GetNumRemainingInputArgs() > 0) {
      if (OPS_GetIntInput(&numdata, &secNum) < 0) {
	opserr << "WARNING sectionLocation eleTag? <secNum?> - could not read int input? \n";
	return -1;
      }
    }
    
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionLocation element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "integrationPoints";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Vector &theVec = *(info.theVector);
    int Np = theVec.Size();

    if (secNum > 0 && secNum <= Np) { // One IP
      double value = theVec(secNum-1);
      numdata = 1;
      if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }
    } else { // All IPs in a list
      std::vector<double> data(Np);
      for (int i = 0; i < Np; i++)
	data[i] = theVec(i);
      numdata = Np;
      if (OPS_SetDoubleOutput(&numdata, &data[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }      
    }        

    delete theResponse;

    return 0;
}

int OPS_sectionWeight()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - sectionWeight eleTag? <secNum?> \n";
	return -1;
    }

    //opserr << "sectionWeight: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 1;
    int tag;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING sectionWeight eleTag? <secNum?> - could not read int input? \n";
	return -1;
    }

    int secNum = 0;
    if (OPS_GetNumRemainingInputArgs() > 0) {
      if (OPS_GetIntInput(&numdata, &secNum) < 0) {
	opserr << "WARNING sectionWeight eleTag? <secNum?> - could not read int input? \n";
	return -1;
      }
    }
    
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionWeight element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "integrationWeights";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Vector &theVec = *(info.theVector);
    int Np = theVec.Size();

    if (secNum > 0 && secNum <= Np) { // One IP
      double value = theVec(secNum-1);
      numdata = 1;
      if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }
    } else { // All IPs in a list
      std::vector<double> data(Np);
      for (int i = 0; i < Np; i++)
	data[i] = theVec(i);
      numdata = Np;
      if (OPS_SetDoubleOutput(&numdata, &data[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }      
    }    

    delete theResponse;

    return 0;
}

int OPS_sectionTag()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - sectionTag eleTag? <secNum?> \n";
	return -1;
    }

    //opserr << "sectionLocation: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 1;
    int tag;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING sectionTag eleTag? <secNum?> - could not read int input? \n";
	return -1;
    }

    int secNum = 0;
    if (OPS_GetNumRemainingInputArgs() > 0) {
      if (OPS_GetIntInput(&numdata, &secNum) < 0) {
	opserr << "WARNING sectionTag eleTag? <secNum?> - could not read int input? \n";
	return -1;
      }
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionTag - element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "sectionTags";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const ID &theID = *(info.theID);
    int Np = theID.Size();

    if (secNum > 0 && secNum <= Np) { // One IP
      int value = theID(secNum-1);
      numdata = 1;
      if (OPS_SetIntOutput(&numdata, &value, true) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }
    } else { // All IPs in a list
      std::vector<int> data(Np);
      for (int i = 0; i < Np; i++)
	data[i] = theID(i);
      numdata = Np;
      if (OPS_SetIntOutput(&numdata, &data[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }      
    }

    delete theResponse;

    return 0;
}

int OPS_sectionDisplacement()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - sectionDisplacement eleTag? secNum? \n";
	return -1;
    }

    //opserr << "sectionWeight: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 2;
    int data[2];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING sectionDisplacement eleTag? secNum? <-local>- could not read int input? \n";
	return -1;
    }

    int tag = data[0];
    int secNum = data[1];
    bool local = false;
    
    if (OPS_GetNumRemainingInputArgs() > 0) {
      const char* localGlobal = OPS_GetString();
      if (strstr(localGlobal,"local") != 0)
	local = true;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionDisplacement element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 2;
    char a[80] = "sectionDisplacements";
    const char *argvv[2];
    argvv[0] = a;
    if (local)
      argvv[1] = "local";
    else
      argvv[1] = "global";

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Matrix &theMatrix = *(info.theMatrix);
    if (secNum <= 0 || secNum > theMatrix.noRows()) {
	opserr << "WARNING invalid secNum\n";
	delete theResponse;
	return -1;
    }

    double value[3];
    value[0] = theMatrix(secNum-1,0);
    value[1] = theMatrix(secNum-1,1);
    value[2] = theMatrix(secNum-1,2);        

    numdata = 3;
    if (OPS_SetDoubleOutput(&numdata, &value[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

SectionForceDeformation::SectionForceDeformation(int tag, int classTag)
  :Material(tag,classTag), fDefault(0), sDefault(0)
{

}

SectionForceDeformation::SectionForceDeformation()
    : Material(0, 0), fDefault(0), sDefault(0)
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
  default:
    k.Invert(*fDefault);
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
  default:
    k.Invert(*fDefault);
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
      case SECTION_RESPONSE_FXX:
          output.tag("ResponseType", "epsXX");
          break;
      case SECTION_RESPONSE_FYY:
          output.tag("ResponseType", "epsYY");
          break;
      case SECTION_RESPONSE_FXY:
          output.tag("ResponseType", "epsXY");
          break;
      case SECTION_RESPONSE_MXX:
          output.tag("ResponseType", "kappaXX");
          break;
      case SECTION_RESPONSE_MYY:
          output.tag("ResponseType", "kappaYY");
          break;
      case SECTION_RESPONSE_MXY:
          output.tag("ResponseType", "kappaXY");
          break;
      case SECTION_RESPONSE_VXZ:
          output.tag("ResponseType", "gammaXZ");
          break;
      case SECTION_RESPONSE_VYZ:
          output.tag("ResponseType", "gammaYZ");
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
      case SECTION_RESPONSE_FXX:
          output.tag("ResponseType", "Fxx");
          break;
      case SECTION_RESPONSE_FYY:
          output.tag("ResponseType", "Fyy");
          break;
      case SECTION_RESPONSE_FXY:
          output.tag("ResponseType", "Fxy");
          break;
      case SECTION_RESPONSE_MXX:
          output.tag("ResponseType", "Mxx");
          break;
      case SECTION_RESPONSE_MYY:
          output.tag("ResponseType", "Myy");
          break;
      case SECTION_RESPONSE_MXY:
          output.tag("ResponseType", "Mxy");
          break;
      case SECTION_RESPONSE_VXZ:
          output.tag("ResponseType", "Vxz");
          break;
      case SECTION_RESPONSE_VYZ:
          output.tag("ResponseType", "Vyz");
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
      case SECTION_RESPONSE_FXX:
          output.tag("ResponseType", "epsXX");
          break;
      case SECTION_RESPONSE_FYY:
          output.tag("ResponseType", "epsYY");
          break;
      case SECTION_RESPONSE_FXY:
          output.tag("ResponseType", "epsXY");
          break;
      case SECTION_RESPONSE_MXX:
          output.tag("ResponseType", "kappaXX");
          break;
      case SECTION_RESPONSE_MYY:
          output.tag("ResponseType", "kappaYY");
          break;
      case SECTION_RESPONSE_MXY:
          output.tag("ResponseType", "kappaXY");
          break;
      case SECTION_RESPONSE_VXZ:
          output.tag("ResponseType", "gammaXZ");
          break;
      case SECTION_RESPONSE_VYZ:
          output.tag("ResponseType", "gammaYZ");
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
      case SECTION_RESPONSE_FXX:
          output.tag("ResponseType", "Fxx");
          break;
      case SECTION_RESPONSE_FYY:
          output.tag("ResponseType", "Fyy");
          break;
      case SECTION_RESPONSE_FXY:
          output.tag("ResponseType", "Fxy");
          break;
      case SECTION_RESPONSE_MXX:
          output.tag("ResponseType", "Mxx");
          break;
      case SECTION_RESPONSE_MYY:
          output.tag("ResponseType", "Myy");
          break;
      case SECTION_RESPONSE_MXY:
          output.tag("ResponseType", "Mxy");
          break;
      case SECTION_RESPONSE_VXZ:
          output.tag("ResponseType", "Vxz");
          break;
      case SECTION_RESPONSE_VYZ:
          output.tag("ResponseType", "Vyz");
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
SectionForceDeformation::setTrialSectionDeformation(const Vector& nouse, const Vector &data) //JZ
{
  opserr << "SectionForceDeformation::setTrialSectionDeformationTemperature (strain, tData) - should not be called\n";
  return -1;
}

static Vector errRes(3);

const Vector &
SectionForceDeformation::getTemperatureStress(const Vector &tData) //PK
{
  opserr << "SectionForceDeformation::getTemperatureStress(double *dataMixed) - should not be called\n";
  errRes.resize(this->getStressResultant().Size());
  return errRes;
  //  return this->getStressResultant();
}
//--- Adding Thermal Functions:[END]   by UoE OpenSees Group ----//

const Vector& SectionForceDeformation::getThermalElong(void)
{
    opserr << "SectionForceDeformation::getThermalElong() - should not be called\n";
    errRes.resize(this->getStressResultant().Size());
    return errRes;
}
