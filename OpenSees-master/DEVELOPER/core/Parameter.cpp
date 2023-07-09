/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
// $Date: 2008-08-26 15:43:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/Parameter.cpp,v $

#include <classTags.h>
#include <Parameter.h>
#include <DomainComponent.h>

Parameter::Parameter(int passedTag,
		     DomainComponent *parentObject,
		     const char **argv, int argc)
  :TaggedObject(passedTag), MovableObject(PARAMETER_TAG_Parameter),
   parameterID(0), theObjects(0), numObjects(0), maxNumObjects(0),
   theComponents(0), numComponents(0), maxNumComponents(0),
   gradIndex(-1)
{
  theInfo.theDouble = 1.0;
  int ok = -1;

  maxNumObjects = initialSize;
  maxNumComponents = initialSize;

  theComponents = new DomainComponent *[maxNumComponents];

  theObjects = new MovableObject *[maxNumObjects];

  parameterID = new int[maxNumObjects];

  for (int i =0; i < maxNumObjects; i++) {
    theObjects[i] = 0;
    parameterID[i] = 0;
  }

  if (parentObject != 0) {
    ok = parentObject->setParameter(argv, argc, *this);
    theComponents[0] = parentObject;
    numComponents = 1;
  }
  else {
    ok = 0; // Created a blank parameter but don't want a warning below
  }

  if (ok < 0) {
    opserr << "Parameter::Parameter "<< this->getTag() <<" -- error encountered while attempting to identify parameter" << endln;
    for (int i = 0; i < argc; i++)
      opserr << argv[i] << ' ';
    opserr << endln;
  }
}

Parameter::Parameter(const Parameter &param):
  TaggedObject(param.getTag()), MovableObject(PARAMETER_TAG_Parameter),
   theComponents(0), numComponents(0), maxNumComponents(0)
{
  theInfo = param.theInfo;
  numComponents = param.numComponents;
  maxNumComponents = param.maxNumComponents;
  numObjects = param.numObjects;
  maxNumObjects = param.maxNumObjects;
  gradIndex = param.gradIndex;

  theComponents = new DomainComponent *[maxNumComponents];
  int i;
  for (i = 0; i < numComponents; i++)
    theComponents[i] = param.theComponents[i];

  theObjects = new MovableObject *[maxNumObjects];
  parameterID = new int[maxNumObjects];

  for (i = 0; i < numObjects; i++) {
    theObjects[i] = param.theObjects[i];
    parameterID[i] = param.parameterID[i];
  }
  for ( ; i < maxNumObjects; i++) {
    theObjects[i] = 0;
    parameterID[i] = 0;
  }
}

Parameter::Parameter(int tag, int classTag)
  :TaggedObject(tag), MovableObject(classTag),
   parameterID(0), theObjects(0), numObjects(0), maxNumObjects(0),
   theComponents(0), numComponents(0), maxNumComponents(0),
   gradIndex(-1)
{

}


Parameter::Parameter()
  :TaggedObject(0), MovableObject(PARAMETER_TAG_Parameter), 
   theObjects(0), 
   theComponents(0), numComponents(0), maxNumComponents(0),
   numObjects(0), maxNumObjects(0), parameterID(0), gradIndex(-1)
{

}


Parameter::~Parameter()
{
  if (theComponents != 0)
    delete [] theComponents;

  if (theObjects != 0)
    delete [] theObjects;

  if (parameterID != 0)
    delete [] parameterID;	
}

int
Parameter::addComponent(int tag, const char **argv, int argc)
{
  opserr << "Parameter::addComponent() - failed\n";
  return -1;
}

int
Parameter::addComponent(DomainComponent *parentObject,
			const char **argv, int argc)
{
  if (numComponents == maxNumComponents) {
    maxNumComponents += expandSize;
    DomainComponent **newComponents = new DomainComponent *[maxNumComponents];

    for (int i = 0; i < numComponents; i++) 
      newComponents[i] = theComponents[i];

    if (theComponents != 0) 
      delete [] theComponents;

    theComponents = newComponents;
  }

  theComponents[numComponents] = parentObject;
  numComponents++;

  //return (parentObject != 0) ?
  //  parentObject->setParameter(argv, argc, *this) : -1;

  int ok = -1;

  int oldNumObjects = numObjects;
  if (parentObject != 0)
    ok = parentObject->setParameter(argv, argc, *this);

  if (numObjects == oldNumObjects || ok < 0) {
    opserr << "Parameter::addComponent "<< this->getTag() <<" -- no objects were able to identify parameter" << endln;
    for (int i = 0; i < argc; i++)
      opserr << argv[i] << ' ';
    opserr << endln;
    return -1;
  }
  return 0;
}

int
Parameter::update(int newValue)
{
  theInfo.theInt = newValue;

  int ok = 0;

  for (int i = 0; i < numObjects; i++)
    ok += theObjects[i]->updateParameter(parameterID[i], theInfo);

  return ok;
}

int
Parameter::update(double newValue)
{
  theInfo.theDouble = newValue;

  int ok = 0;

  for (int i = 0; i < numObjects; i++)
    ok += theObjects[i]->updateParameter(parameterID[i], theInfo);
  
  return ok;
}

int
Parameter::activate(bool active)
{
  int ok = 0;
  for (int i = 0; i < numObjects; i++) {
    ok += theObjects[i]->activateParameter(active ? parameterID[i] : 0);
  }

  return ok;
}

double
Parameter::getSensitivity(int index)
{
  //return 1.0;
  return (index == gradIndex) ? 1.0 : 0.0;
}

void
Parameter::Print(OPS_Stream &s, int flag)  
{
  s << "Parameter, tag = " << this->getTag() << " " << this->getValue() << endln;
  //s << "\tparameterID = " << parameterID << endln;
}


int
Parameter::addObject(int paramID, MovableObject *object)
{
  if (numObjects == maxNumObjects) {
    maxNumObjects += expandSize;
    MovableObject **newObjects = new MovableObject*[maxNumObjects];
    int *newParameterID = new int[maxNumObjects];

    for (int i = 0; i < numObjects; i++) {
      newObjects[i] = theObjects[i];
      newParameterID[i] = parameterID[i];
    }

    if (theObjects != 0) delete [] theObjects;
    if (parameterID != 0) delete [] parameterID;

    theObjects = newObjects;
    parameterID = newParameterID;
  }

  parameterID[numObjects] = paramID;
  theObjects[numObjects] = object;
  numObjects++;

  return 0;
}


void
Parameter::setDomain(Domain *theDomain)
{
  return;
}

int 
Parameter::sendSelf(int commitTag, Channel &theChannel)
{
  opserr << "Parameter::sendSelf - not yet implemented\n";
  //return -1;
  return 0;
}

int 
Parameter::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

int 
Parameter::clean(void)
{
  for (int i = 0; i < numObjects; i++) {
    theObjects[i] = 0;
  }
  for (int i = 0; i < numComponents; i++) {
    theComponents[i] = 0;
  }

  numObjects = 0;
  numComponents = 0;
  currentValue = 0.0;

  return 0;
}
