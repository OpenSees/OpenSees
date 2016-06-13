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
                                                                        
// $Revision: 1.1 $
// $Date: 2010-05-13 00:15:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/CompositeResponse.cpp,v $
                                                                        
// Written: fmk
// Created: 05/10
//
// Description: This file contains the CompositeResponse class interface.
// A CompositeResponse is a container holding a number of response objects.

#include <CompositeResponse.h>
#include <Vector.h>
#include <ID.h>

CompositeResponse::CompositeResponse(void)
  :Response(), theResponses(0), numResponses(0)
{

}

CompositeResponse::~CompositeResponse()
{
  for (int i=0; i<numResponses; i++)
    delete theResponses[i];

  delete [] theResponses;
}


int 
CompositeResponse::addResponse(Response *nextResponse)
{

  if (nextResponse == 0)
    return 0;

  //
  // check valid type double, Vector, int, ID only
  // resize theType accordiingly (vector or ID)
  //

  Information&otherType = nextResponse->getInformation();

  if (otherType.theType == DoubleType || 
      otherType.theType == VectorType) {

    if (myInfo.theType == UnknownType) {
      myInfo.theType = VectorType;
      myInfo.theVector = new Vector();
    }
    
    if (myInfo.theType != VectorType) {
      opserr << "WARNING: CompositeResponse::addResponse() - mismatching type, no responses will be addeed\n";
      return -1;
    }
    
    int curSize = myInfo.theVector->Size();

    if (otherType.theType == DoubleType)
      curSize++;
    else
      curSize += otherType.theVector->Size();

    myInfo.theVector->resize(curSize);

  } else if (otherType.theType == IntType || 
	     otherType.theType == IdType) {

    if (myInfo.theType == UnknownType) {
      myInfo.theID = new ID();
      myInfo.theType = IdType;    
    }
    
    if (myInfo.theType != IdType) {
      opserr << "WARNING: CompositeResponse::addResponse() - mismatching type, no responses will be addeed\n";
      return -1;
    }
    
    int curSize = myInfo.theID->Size();
    if (otherType.theType == IntType)
      curSize++;
    else
      curSize += otherType.theID->Size();
    
    myInfo.theID->resize(curSize);
  }


  //
  // now add the response to the array
  //

  Response **theNextResponses = new Response *[numResponses+1];
  if (theNextResponses == 0) {
    opserr << "WARNING: CompositeResponse::addResponse() - out of memory, no responses will be added\n";
    return -1;
  }
   
  for (int i=0; i<numResponses; i++)
    theNextResponses[i] = theResponses[i];
  if (theResponses != 0)
    delete [] theResponses;
  theResponses = theNextResponses;
  theResponses[numResponses] = nextResponse;
  numResponses++;

  return numResponses;
}


int 
CompositeResponse::getResponse(void)
{
  int res = 0;

  //
  // invoke getResponse on all responses & add the data to myInfo
  //

  int currentLoc = 0;
  for (int i=0; i<numResponses; i++) {
    Response *theResponse = theResponses[i];
    res += theResponse->getResponse();

    Information&otherType = theResponse->getInformation();

    if (otherType.theType == DoubleType || otherType.theType == VectorType) {
      if (otherType.theType == DoubleType)
	(*myInfo.theVector)(currentLoc++) = otherType.theDouble;
      else {
	int otherSize = otherType.theVector->Size();
	for (int i=0; i<otherSize; i++, currentLoc++) 
	  (*myInfo.theVector)(currentLoc) = (*otherType.theVector)(i);
      }
    } else if (otherType.theType == IntType || otherType.theType == IdType) {
      if (otherType.theType == IntType) {
	(*myInfo.theID)(currentLoc++) = otherType.theInt;
      }
      else {
	int otherSize = otherType.theID->Size();
	for (int i=0; i<otherSize; i++, currentLoc++) 
	  (*myInfo.theID)(currentLoc) = (*otherType.theID)(i);
      }    
    }    
  }
  
  return res;
}


