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

// $Revision: 1.5 $
// $Date: 2009-08-19 17:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CrdTransf.cpp,v $

// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Description: This file contains the implementation for the CrdTransf class.
// CrdTransf provides the abstraction of a frame 
// coordinate transformation. It is an abstract base class and 
// thus no objects of its type can be instatiated. 

#include <CrdTransf.h>
#include <ID.h>
#include <Vector.h>

#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>
#include <CrdTransfResponse.h>

static MapOfTaggedObjects theCrdTransfObjects;

bool 
OPS_addCrdTransf(CrdTransf *newComponent) {
  return theCrdTransfObjects.addComponent(newComponent);
}

bool OPS_removeCrdTransf(int tag)
{
    TaggedObject* obj = theCrdTransfObjects.removeComponent(tag);
    if (obj != 0) {
	delete obj;
	return true;
    }
    return false;
}

CrdTransf *
OPS_getCrdTransf(int tag) {

  TaggedObject *theResult = theCrdTransfObjects.getComponentPtr(tag);
  if (theResult == 0) {
    opserr << "CrdTransf *getCrdTransf(int tag) - none found with tag: " << tag << endln;
    return 0;
  }
  CrdTransf *theSeries = (CrdTransf *)theResult;

  return theSeries;
}

void 
OPS_clearAllCrdTransf(void) {
  theCrdTransfObjects.clearAll();
}


void OPS_printCrdTransf(OPS_Stream &s, int flag) {
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\"crdTransformations\": [\n";        
    MapOfTaggedObjectsIter theObjects = theCrdTransfObjects.getIter();
    theObjects.reset();
    TaggedObject *theObject;
    int count = 0;
    int numComponents = theCrdTransfObjects.getNumComponents();    
    while ((theObject = theObjects()) != 0) {
      CrdTransf *theTransf = (CrdTransf *)theObject;
      theTransf->Print(s, flag);
      if (count < numComponents-1)
	s << ",\n";
      count++;      
    }
    s << "\n\t\t]";
  }
}

ID OPS_getAllCrdTransfTags() {

    ID allCrdTransfTags(0);
      
    MapOfTaggedObjectsIter theObjects = theCrdTransfObjects.getIter();
    theObjects.reset();
    TaggedObject *theObject;

    while ((theObject = theObjects()) != 0) {
      CrdTransf *theTransf = (CrdTransf *)theObject;    
      allCrdTransfTags.insert(theTransf->getTag());
    }

    return allCrdTransfTags;
}



// constructor:
CrdTransf::CrdTransf(int tag, int classTag):TaggedObject(tag), MovableObject(classTag)
{
}

// destructor:
CrdTransf::~CrdTransf()
{
}

int
CrdTransf::getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis)
{
  xAxis.Zero();
  yAxis.Zero();
  zAxis.Zero();
  
  return 0;
}

int
CrdTransf::getRigidOffsets(Vector &offsets)
{
  offsets.Zero();
  
  return 0;
}

Response*
CrdTransf::setResponse(const char **argv, int argc, OPS_Stream &theHandler)
{
  if (argc < 1)
    return 0;

  Response *theResponse = 0;
  
  if (strcmp(argv[0],"xaxis") == 0 || strcmp(argv[0],"xlocal") == 0)
    theResponse = new CrdTransfResponse(this, 201, Vector(3));
  
  if (strcmp(argv[0],"yaxis") == 0 || strcmp(argv[0],"ylocal") == 0)
    theResponse = new CrdTransfResponse(this, 202, Vector(3));
  
  if (strcmp(argv[0],"zaxis") == 0 || strcmp(argv[0],"zlocal") == 0)
    theResponse = new CrdTransfResponse(this, 203, Vector(3));

  if (strcmp(argv[0],"offsets") == 0 || strcmp(argv[0],"rigidOffsets") == 0)
    theResponse = new CrdTransfResponse(this, 204, Vector(6));
  
  return theResponse;
}

int
CrdTransf::getResponse(int responseID, Information &eleInfo)
{
    if (responseID >= 201 && responseID <= 203) {
        static Vector xlocal(3);
        static Vector ylocal(3);
        static Vector zlocal(3);
        
        this->getLocalAxes(xlocal, ylocal, zlocal);
        
        if (responseID == 201)
            return eleInfo.setVector(xlocal);
        else if (responseID == 202)
            return eleInfo.setVector(ylocal);
        else if (responseID == 203)
            return eleInfo.setVector(zlocal);
        else
            return -1;
    }
    if (responseID == 204) {
      static Vector offsets(6);

      this->getRigidOffsets(offsets);

      return eleInfo.setVector(offsets);
    }
    else
        return -1;
}

const Vector &
CrdTransf::getBasicDisplSensitivity(int gradNumber)
{
    opserr << "WARNING CrdTransf::getBasicDisplSensitivity() - this method "
        << " should not be called." << endln;
    
    static Vector dummy(1);
    return dummy;
}

const Vector &
CrdTransf::getGlobalResistingForceShapeSensitivity(const Vector &pb,
						   const Vector &p0,
						   int gradNumber)
{
    opserr << "ERROR CrdTransf::getGlobalResistingForceSensitivity() - has not been"
        << " implemented yet for the chosen transformation." << endln;
    
    static Vector dummy(1);
    return dummy;
}

const Vector &
CrdTransf::getGlobalResistingForceShapeSensitivity(const Vector &pb,
						   const Vector &p0)
{
    opserr << "ERROR CrdTransf::getGlobalResistingForceSensitivity() - has not been"
        << " implemented yet for the chosen transformation." << endln;
    
    static Vector dummy(1);
    return dummy;
}

const Vector &
CrdTransf::getBasicTrialDispShapeSensitivity(void)
{
    opserr << "ERROR CrdTransf::getBasicTrialDispShapeSensitivity() - has not been"
        << " implemented yet for the chosen transformation." << endln;
    
    static Vector dummy(1);
    return dummy;
}



// --Quan
const Vector &
CrdTransf::getBasicDisplSensitivity(int gradNumber, int)
{
    opserr << "WARNING CrdTransf::getBasicDisplSensitivity() - this method "
        << " should not be called." << endln;
    
    static Vector dummy(1);
    return dummy;
}
