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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-04-14 21:34:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSection.cpp,v $
                                                                        
// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the class definition for 
// FiberSection.h. FiberSection provides the abstraction of a 
// section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

#include <stdlib.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSection.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>

// constructors:
FiberSection::FiberSection(int tag, int num, Fiber **fibers): 
  SectionForceDeformation(tag, SEC_TAG_Fiber),
  numFibers(num), theFibers(0), sizeFibers(num),
  e(0), eCommit(0), s(0), ks(0), order(0), code(0), otherDbTag(0)
{
    theFibers = new Fiber *[numFibers];

    if (theFibers == 0) {
      opserr << "FiberSection::FiberSection -- failed to allocate Fiber pointers";
      exit(-1);
    }

    for (int i = 0; i < numFibers; i++) {
      theFibers[i] = fibers[i]->getCopy();

      if (theFibers[i] == 0) {
	opserr << "FiberSection::FiberSection -- failed to get copy of Fiber";
	exit(-1);
      }
			      
    }

	order = theFibers[0]->getOrder();

	e = new Vector(order);
	eCommit = new Vector(order);
	s = new Vector(order);
	ks = new Matrix(order,order);

	code = new ID(order);
	*code = theFibers[0]->getType();
}

// constructor for blank object that recvSelf needs to be invoked upon
FiberSection::FiberSection():
  SectionForceDeformation(0, SEC_TAG_Fiber),
  numFibers(0), theFibers(0), sizeFibers(0),
  e(0), eCommit(0), s(0), ks(0), order(0), code(0), otherDbTag(0)
{

}

FiberSection::FiberSection(int tag, int num):
  SectionForceDeformation(tag, SEC_TAG_Fiber),
  numFibers(0), theFibers(0), sizeFibers(num),
  e(0), eCommit(0), s(0), ks(0), order(0), code(0), otherDbTag(0)
{
	// check for zero -- if zero make initial size 2
	if (sizeFibers == 0)
		sizeFibers = 2;

	// create the array of fiber pointers
	theFibers = new Fiber *[sizeFibers];

	if (theFibers == 0) {
	  opserr << "FiberSection::FiberSection -- failed to allocate Fiber pointers";
	  exit(-1);
	}

	// zero the pointers
	for (int i = 0; i < sizeFibers; i++)
		theFibers[i] =0;
}

int
FiberSection::addFiber(Fiber &newFiber)
{
	if (order == 0) {
		order = newFiber.getOrder();

		e = new Vector(order);
		eCommit = new Vector(order);
		s = new Vector(order);
		ks = new Matrix(order,order);

		code = new ID(order);
		*code = newFiber.getType();
	}

	if (numFibers < sizeFibers) {
		// space available in array .. set new pointer and increment number
		theFibers[numFibers] = &newFiber;
		numFibers++;
	}
	else {
		// need to create a larger array
		int newSize = 2*numFibers;
		if (newSize == 0) 
			newSize = 2; // in case failed in constructor

		Fiber **newArray = new Fiber *[newSize]; 

	    if (newArray == 0) {
	      opserr << "FiberSection::addFiber -- failed to allocate Fiber pointers";
	      exit(-1);

			return -1;
		}

		// set the new size of the array
		sizeFibers = newSize;

		// copy the old pointers
		for (int i = 0; i < numFibers; i++)
			newArray[i] = theFibers[i];

		// add the new pointer
		newArray[numFibers] = &newFiber;
		numFibers++;

		// zero the last elements of the array
		for (int j = numFibers; j < newSize; j++) 
			newArray[j] = 0;
    
		delete [] theFibers;
		
		theFibers = newArray;
	}
  
	return 0;
}

// destructor:
FiberSection::~FiberSection()
{
	if (theFibers) {
		for (int i = 0; i < numFibers; i++)
			if (theFibers[i])
				delete theFibers[i];
      
		delete [] theFibers;
	}

	if (e)
		delete e;
	if (eCommit)
		delete eCommit;
	if (s)
		delete s;
	if (ks)
		delete ks;
	if (code)
		delete code;
}

int FiberSection::setTrialSectionDeformation (const Vector &deforms)
{
	int res = 0;

	*e = deforms;

	for (int i = 0; i < numFibers; i++)
		res += theFibers[i]->setTrialFiberStrain(deforms);

	return res;
}

const Vector&
FiberSection::getSectionDeformation(void)
{
	return *e;
}

const Matrix&
FiberSection::getSectionTangent(void)
{
	ks->Zero();
 
	for (int i = 0; i < numFibers; i++)
		ks->addMatrix(1.0, theFibers[i]->getFiberTangentStiffContr(), 1.0);

	return *ks;
}

const Matrix&
FiberSection::getInitialTangent(void)
{
	ks->Zero();
 
	for (int i = 0; i < numFibers; i++)
	  // Should be initial stiffness contribution
	  ks->addMatrix(1.0, theFibers[i]->getFiberTangentStiffContr(), 1.0);

	return *ks;
}

const Vector&
FiberSection::getStressResultant(void)
{
	s->Zero();
 
	for (int i = 0; i < numFibers; i++)
		s->addVector(1.0, theFibers[i]->getFiberStressResultants(), 1.0);

	return *s;
}

SectionForceDeformation*
FiberSection::getCopy(void)
{
	FiberSection *theCopy = new FiberSection (this->getTag(), numFibers, theFibers);
  
	*(theCopy->eCommit) = *eCommit;

	return theCopy;
}

const ID&
FiberSection::getType ()
{
    return *code;
}

int
FiberSection::getOrder () const
{
    return order;
}

int
FiberSection::commitState(void)
{
	// commit the fibers state
	int err = 0;

	for (int i = 0; i < numFibers; i++)
		err += theFibers[i]->commitState();

	*eCommit = *e;

	return err;
}

int
FiberSection::revertToLastCommit(void)
{
	int err = 0;

	// Last committed section deformations
	*e = *eCommit;

	s->Zero();
	ks->Zero();

	for (int i = 0; i < numFibers; i++) {
		// Revert the fiber ...
		err += theFibers[i]->revertToLastCommit();
		
		// ... now recompute the section stress resultant and tangent stiffness
		theFibers[i]->setTrialFiberStrain(*e);
		s->addVector(1.0, theFibers[i]->getFiberStressResultants(), 1.0);
		ks->addMatrix(1.0, theFibers[i]->getFiberTangentStiffContr(), 1.0);
	}

	return err;
}

int
FiberSection::revertToStart(void)
{
	// revert the fibers to start    
	int err = 0;
   
	for (int i = 0; i < numFibers; i++)
		err += theFibers[i]->revertToStart();
    
	eCommit->Zero();

	return err;
}

int
FiberSection::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	// Need two dbTags since data and dbTag ID vectors could be the same size
	if (otherDbTag == 0)
		otherDbTag = theChannel.getDbTag();

	// Create an ID with tag, number and size of fibers, and section order
	static ID data(5);

	data(0) = this->getTag();
	data(1) = numFibers;
	data(2) = sizeFibers;
	data(3) = order;
	data(4) = otherDbTag;

	// Send the data ID
	res += theChannel.sendID(this->getDbTag(), commitTag, data);
	if (res < 0) {
	  opserr << "FiberSection::sendSelf -- failed to send data ID";
	  return res;
	}

	if (order > 0) {
		// Send the committed section deformations
		res += theChannel.sendVector(this->getDbTag(), commitTag, *eCommit);
		if (res < 0) {
		  opserr << "FiberSection::sendSelf -- failed to send section deformations";
		  return res;
		}
	}

	// If there are any fibers to send ...
	if (numFibers > 0) {
		// Create an ID for fiber dbTags
		ID dbTags(numFibers+1);
		
		int i;
		int dbTag;
		for (i = 0; i < numFibers; i++) {
			dbTag = theFibers[i]->getDbTag();

			if (dbTag == 0) {
				dbTag = theChannel.getDbTag();
				if (dbTag != 0)
					theFibers[i]->setDbTag(dbTag);
			}

			dbTags(i) = dbTag;
		}
	
		// Get class tag from first fiber
		dbTags(numFibers) = theFibers[0]->getClassTag();

		// Send the dbTags ID
		res += theChannel.sendID(otherDbTag, commitTag, dbTags);
		if (res < 0) {
		  opserr << "FiberSection::sendSelf -- failed to send dbTags ID";
		  return res;
		}

		// Ask the fibers to send themselves
		for (i = 0; i < numFibers; i++) {
			res += theFibers[i]->sendSelf(commitTag, theChannel);
			if (res < 0) {
			  opserr << "FiberSection::sendSelf -- failed to send Fiber %d";
			  return res;
			}
		}
	}

	return res;
}

int
FiberSection::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;

	// Create an ID for tag, number and size of fibers, and section order
	static ID data(5);

	// Receive the data ID
	res += theChannel.recvID(this->getDbTag(), commitTag, data);
	if (res < 0) {
	  opserr << "FiberSection::recvSelf -- failed to receive data ID";
	  return res;
	}

	this->setTag(data(0));
	numFibers = data(1);
	sizeFibers = data(2);
	order = data(3);
	otherDbTag = data(4);

	// Check that section variables are of the right size, allocate if needed
	if (order > 0) {
		if (e == 0)	
			e = new Vector(order);
		if (eCommit == 0)	
			eCommit = new Vector(order);
		if (s == 0) 
			s = new Vector(order);
		if (ks == 0) 
			ks = new Matrix(order,order);
		if (code == 0) 
			code = new ID(order);

		if (e->Size() != order) {
			delete e;
			e = new Vector(order);
		}
		if (eCommit->Size() != order) {
			delete eCommit;
			eCommit = new Vector(order);
		}
		if (s->Size() != order) {
			delete s;
			s = new Vector(order);
		}
		if (ks->noRows() != order) {
			delete ks;
			ks = new Matrix(order,order);
		}
		if (code->Size() != order) {
			delete code;
			code = new ID(order);
		}

		// Receive the committed section deformation
		res += theChannel.recvVector(this->getDbTag(), commitTag, *eCommit);
		if (res < 0) {
		  opserr << "FiberSection::recvSelf -- failed to receive section deformations";
		  return res;
		}

		*e = *eCommit;
	}

	if (numFibers > 0) {
		// Create an ID for fiber dbTags
		ID dbTags(numFibers+1);

		// Receive the dbTags ID
		res += theChannel.recvID(otherDbTag, commitTag, dbTags);
		if (res < 0) {
		  opserr << "FiberSection::recvSelf -- failed to receive dbTags ID";
		  return res;
		}

		int i;

		if (theFibers == 0) {
			theFibers = new Fiber *[sizeFibers];
			if (theFibers == 0) {
			  opserr << "FiberSection::recvSelf -- failed to allocate Fiber pointers";
			  return -1;
			}
			for (i = 0; i < sizeFibers; i++)
				theFibers[i] = 0;
		}

		int fiberClassTag = dbTags(numFibers);
		for (i = 0; i < numFibers; i++) {
			// Check for null fiber, allocate if so
			if (theFibers[i] == 0)
				theFibers[i] = theBroker.getNewFiber(fiberClassTag);
		
			// Check that the Fiber is of the right type; if not, delete
			// the current one and get a new one of the right type
			else if (theFibers[i]->getClassTag() != fiberClassTag) {
				delete theFibers[i];
				theFibers[i] = theBroker.getNewFiber(fiberClassTag);
			}

			// Check if either allocation failed
			if (theFibers[i] == 0) {
			  opserr << "FiberSection::recvSelf -- could not get Fiber %d";
			  return -1;
			}

			// Now, receive the Fiber
			theFibers[i]->setDbTag(dbTags(i));
			res += theFibers[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
			  opserr << "FiberSection::recvSelf -- could not receive Fiber %d";
			  return res;
			}
		}

		// Reset the section code
		*code = theFibers[0]->getType();
	}

	//this->revertToLastCommit();
	
	return res;
}

void
FiberSection::Print(OPS_Stream &s, int flag)
{
	s << "\nFiberSection, tag: " << this->getTag() << endln;
	s << "\tSection code: " << *code;
	s << "\tNumber of Fibers: " << numFibers << endln;
	
	if (flag == 1)
	  for (int i = 0; i < numFibers; i++)
	    theFibers[i]->Print(s, flag);


}

Response*
FiberSection::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("SectionOutput");
  output.attr("secType", this->getClassType());
  output.attr("secTag", this->getTag());

  // Check if fiber response is requested
  if (strcmp(argv[0],"fiber") == 0) {
    int key = 0;
    int passarg = 2;
    
    if (argc <= 2)          // not enough data input
      return 0;
    else if (argc <= 3)		// fiber number was input directly
      key = atoi(argv[1]);
    else {                  // fiber near-to coordinate specified
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double ySearch, zSearch;
      theFibers[0]->getFiberLocation(ySearch,zSearch);
      double closestDist = sqrt( pow(ySearch-yCoord,2) +
				 pow(zSearch-zCoord,2) );
      double distance;
      for (int j = 1; j < numFibers; j++) {
	theFibers[j]->getFiberLocation(ySearch,zSearch);
	distance = sqrt( pow(ySearch-yCoord,2) +
			 pow(zSearch-zCoord,2) );
	if (distance < closestDist) {
	  closestDist = distance;
	  key = j;
	}
      }
      theFibers[key]->getFiberLocation(ySearch,zSearch);
      passarg = 3;
    }
	
    if (key < numFibers && key >= 0) {
	output.tag("FiberOutput");
	double yLoc, zLoc;
	theFibers[key]->getFiberLocation(yLoc, zLoc);
	output.attr("yLoc", yLoc);
	output.attr("zLoc",0.0);
	output.attr("area", zLoc);
	theResponse = theFibers[key]->setResponse(&argv[passarg],argc-passarg, output);
	output.endTag();
    }    
  }

  output.endTag();
  return SectionForceDeformation::setResponse(argv,argc,output);
}

int 
FiberSection::getResponse(int responseID, Information &sectInfo)
{
	// Just call the base class method ... don't need to define
	// this function, but keeping it here just for clarity
	return SectionForceDeformation::getResponse(responseID, sectInfo);
}
