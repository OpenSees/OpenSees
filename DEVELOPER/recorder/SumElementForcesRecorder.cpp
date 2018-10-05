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
                                                                        
// $Revision$
// $Date$
// $Source$
                                                                        
// Written: fmk 
// Created: 05/10
//
// Description: This file contains the class implementatation of SumElementForcesRecorder.
//
// What: "@(#) SumElementForcesRecorder.C, revA"

#include "SumElementForcesRecorder.h"
#include <elementAPI.h>

#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <string.h>
#include <Response.h>
#include <Message.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <StandardStream.h>
#include <BinaryFileStream.h>
#include <DataFileStream.h>

#include <elementAPI.h>

SumElementForcesRecorder::SumElementForcesRecorder()
:Recorder(-1),
 numEle(0), theElements(0), eleID(0), 
 theDomain(0), theOutput(0),
 echoTimeFlag(true), data(0)
{

}

SumElementForcesRecorder::SumElementForcesRecorder(const ID ele,
						   bool echoTime,
						   OPS_Stream *theoutput)
 :Recorder(-1),
 numEle(0), theElements(0), eleID(ele), 
 theDomain(0), theOutput(theoutput),
 echoTimeFlag(echoTime), data(0)
{

}


SumElementForcesRecorder::~SumElementForcesRecorder()
{
  if (theElements != 0)
    delete [] theElements;

  if (data != 0)
    delete data;

  if (theOutput != 0)
    delete theOutput;
}

int 
SumElementForcesRecorder::record(int commitTag, double timeStamp)
{
  // check for initialization
  if (data == 0) {
    opserr << "SumElementForcesRecorder::record() - setDomain() has not been called\n";
    return -1;
  }

  // zero the data vector
  data->Zero();

  int forceSize = data->Size();
  int startLoc = 0;


  // write the time if echTimeFlag set
  if (echoTimeFlag == true) {
    (*data)(0) = timeStamp;
    forceSize -= 1;
    startLoc = 1;
  }
    
  //
  // for each element that has been added to theElements add force contribution
  //

  for (int i=0; i< numEle; i++) {
    if (theElements[i] != 0) {
      int loc = startLoc;
      const Vector &force = theElements[i]->getResistingForce();
      int forceSize = force.Size();
      for (int j=0; j<forceSize; j++, loc++)
	  (*data)(loc) += force(j);
    }
  } 

  //
  // send the response vector to the output handler for o/p
  //

  theOutput->write(*data);

  // successful completion - return 0
  return 0;
}

int
SumElementForcesRecorder::restart(void)
{
  return 0;
}


int
SumElementForcesRecorder::domainChanged(void)
{
  if (theDomain != 0)
    this->setDomain(*theDomain);
}


int 
SumElementForcesRecorder::setDomain(Domain &theDom)
{
  theDomain = &theDom;

  // delete old data
  if (theElements != 0)
    delete [] theElements;
  if (data != 0)
    delete data;

  // set numEle
  numEle = eleID.Size();
  if (numEle == 0) {
    opserr << "WARNING SumElementForcesRecorder::initialize() - no elements tags passed in input!\n";
    return 0; 
  }

  // create theElements, an array of pointers to elements
  theElements = new Element *[numEle];
  if (theElements == 0) {
    opserr << "WARNING SumElementForcesRecorder::initialize() - out of memory\n";
    numEle = 0; // set numEle = 0, in case record() still called
    return -1;
  }

  // zero theElements array
  for (int k=0; k<numEle; k++)
    theElements[k] = 0;

  //
  // loop over the list of elements, 
  //    if element exists add it's pointer o the array
  //    get its resisting force, check size to determine compatible with others
  // 

  int sizeArray = -1;

  for (int i=0; i<numEle; i++) {
    int eleTag = eleID(i);

    Element *theEle = theDomain->getElement(eleTag);
    
    if (theEle != 0) {

      const Vector &force = theEle->getResistingForce();
      int forceSize = force.Size();
      if (sizeArray == -1) {
	sizeArray = forceSize;
	theElements[i] = theEle;
      } else if (sizeArray != forceSize) {
	opserr << "WARNING: forces mismatch - element: " << eleTag << " will not be included\n";
      } else {
	theElements[i] = theEle;	
      }
    }
  }

  // if echTimeFlag is set, add room for the time to be output
  if (echoTimeFlag == true)
    sizeArray++;

  // create the vector to hold the data
  data = new Vector(sizeArray);
  
  if (data == 0 || data->Size() != sizeArray) {
    opserr << "SumElementForcesRecorder::initialize() - out of memory\n";
    delete [] theElements; 
    theElements = 0; 
    numEle = 0;
  }

  return 0;
}

static char myClassType[] = {"SumElementForcesRecorder"};

const char *
SumElementForcesRecorder::getClassType(void) const
{
  return myClassType;
}

int
SumElementForcesRecorder::sendSelf(int commitTag, Channel &theChannel)
{
  // send in an ID (integar array) to the receiving object the following:
  //   recorder tag
  //   size of eleID
  //   class tag of handler
  //   echoTimeFlag

  static ID idData(5);
  idData(0) = this->getTag();;
  idData(1) = eleID.Size();
  idData(2) = theOutput->getClassTag();
  if (echoTimeFlag == true)
    idData(3) = 1;
  else
    idData(3) = 0;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "SumElementForcesRecorder::recvSelf() - failed to recv idData\n";
    return -1;
  }

  // send eleID  to receiving object
  if (theChannel.sendID(0, commitTag, eleID) < 0) {
    opserr << "SumElementForcesRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }

  // send theOutput to receiving object
  if (theOutput->sendSelf(commitTag, theChannel) < 0) {
    opserr << "SumElementForcesRecorder::sendSelf() - failed to send theOutput\n";
    return -1;
  }

  return 0;
}

int 
SumElementForcesRecorder::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  // receive from the sending object the ID
  static ID idData(5);
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "SumElementForcesRecorder::recvSelf() - failed to recv idData\n";
    return -1;
  }

  // with the data received
  //   setTag
  //   resize the eleID array
  //   set echoTimeFlag
  //   get an outputHandler

  this->setTag(idData(0));
  eleID.resize(idData(1));
  idData(2) = theOutput->getClassTag();
  if (idData(3)  == 0)
    echoTimeFlag = true;
  else
    echoTimeFlag = false;

  if (theOutput != 0 && theOutput->getClassTag() != idData(4)) 
    delete theOutput;

  theOutput = theBroker.getPtrNewStream(idData(4));
  if (theOutput == 0) {
    opserr << "SumElementForcesRecorder::recvSelf() - failed to get Output of correct type\n";
    return -1;
  }

  // receive eleID
  if (theChannel.recvID(0, commitTag, eleID) < 0) {
    opserr << "SumElementForcesRecorder::recvSelf() - failed to recv eleID\n";
    return -1;
  }

  // get theOutput to receive data
  if (theOutput->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "SumElementForcesRecorder::sendSelf() - failed to send theOutput\n";
    return -1;
  }  


  return 0;
}


#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numSumElementForcesREcorder = 0;

OPS_Export void *
OPS_SumElementForcesRecorder()
{
  Recorder *theRecorder = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  // check for quick return, possibly parallel case
  if (numRemainingArgs == 0) {
    Recorder *theRecorder = new SumElementForcesRecorder();
  }

  //
  // parse args
  //

  int numEle = 0, eleTag;
  ID eleID(0);

  OPS_Stream *theOutputStream = 0;
  int outMode = 0;    // standard stream
  bool echoTime = false;

  bool doneParsingArgs = false;
  char data[100];
  char outputName[200];
  char **eleArgs = 0;
  int numEleArgs = 0;


  while (numRemainingArgs > 0) {
    if (OPS_GetString(data,100) < 0) 
      return 0;

    // output to standard file
    if (strcmp(data,"-file") == 0)  { 
      outMode = 1;
      if (OPS_GetString(outputName,200) < 0)
	return 0;
      numRemainingArgs -= 2;
    }

    // output to binary file
    else if (strcmp(data,"-binary") == 0) { 
      outMode = 2;
      if (OPS_GetString(outputName,200) < 0)
	return 0;
      numRemainingArgs -= 2;
    }

    // echo domain time stamp in output
    else if (strcmp(data,"-time") == 0) {
      echoTime = true;
      numRemainingArgs -= 1;
    }

    // read the list of elements & place in an ID
    else if ((strcmp(data,"-ele") == 0) ||
	     (strcmp(data,"-eles") == 0) ||
	     (strcmp(data,"-element") == 0)) {

      numRemainingArgs --;      
      int one = 1;
      while (numRemainingArgs > 0 && OPS_GetIntInput(&one, &eleTag) == 0) {
	eleID[numEle] = eleTag;
	numEle++;
	numRemainingArgs--;
      }
      doneParsingArgs = true;
    }
    
    else
      doneParsingArgs = true;
  }

  //
  // create the output handler
  //

  if (outMode == 0) 
    theOutputStream = new StandardStream();
  if (outMode == 1)
    theOutputStream = new DataFileStream(outputName);
  else if (outMode == 2) 
    theOutputStream = new BinaryFileStream(outputName);

  //
  // create the recorder
  //

  theRecorder = new SumElementForcesRecorder(eleID,
					     echoTime, 
					     theOutputStream);
  
  // return it
  return theRecorder;
}
