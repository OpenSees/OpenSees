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
// $Date: 2009-03-20 22:44:17 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewRecorder/ElementSumRecorder.cpp,v $
                                                                        
// Written: fmk 
// Created: 09/99
//
// Description: This file contains the class implementatation of ElementSumRecorder.
//
// What: "@(#) ElementSumRecorder.C, revA"

#include "ElementSumRecorder.h"
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <string.h>
#include <Response.h>
#include <FE_Datastore.h>
#include <Message.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <StandardStream.h>
#include <DataFileStream.h>
#include <XmlFileStream.h>
#include <BinaryFileStream.h>

#include <elementAPI.h>

#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

// init function called once on first loading of package

OPS_Export void
localInit() 
{
  OPS_Error("ElementSumRecorder - Written by fmk UC Berkeley Copyright 2009\n", 1);
}

OPS_Export void *
OPS_ElementSumRecorder()
{
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
  
  while (doneParsingArgs  == false) {

    if (OPS_GetCharInput(data,100) < 0) 
      return 0;

    // output to standard file
    if (strcmp(data,"-file") == 0)  { 
      outMode = 1;
      if (OPS_GetCharInput(outputName,200) < 0)
	return 0;
    }

    // output to xml file
    else if (strcmp(data,"-xml") == 0) { 
      outMode = 2;
      if (OPS_GetCharInput(outputName,200) < 0)
	return 0;
    }

    // output to binary file
    else if (strcmp(data,"-binary") == 0) { 
      outMode = 2;
      if (OPS_GetCharInput(outputName,200) < 0)
	return 0;
    }

    // echo domain time stamp in output
    else if (strcmp(data,"-time") == 0) {
      echoTime = true;
    }


    // read the list of elements & place in an ID
    else if ((strcmp(data,"-ele") == 0) ||
	     (strcmp(data,"-eles") == 0) ||
	     (strcmp(data,"-element") == 0)) {
      
      int one = 1;
      while (OPS_GetIntInput(&one, &eleTag) == 0) {
	eleID[numEle] = eleTag;
	numEle++;
      }
      doneParsingArgs = true;
    }
    
    else
      doneParsingArgs = true;
  }

  doneParsingArgs = false;

  // 
  // now parsing end of line & creating list that will be 
  // provided to setResponse in element.
  //

  while (doneParsingArgs  == false) {
    if (OPS_GetCharInput(data,100) < 0) 
      doneParsingArgs = true;
    else {

      char **nextEleArgs = new char *[numEleArgs+1];
      for (int i=0; i<numEleArgs; i++)
	nextEleArgs[i] = eleArgs[i];

      char *nextArg = new char[strlen(data)+1];
      strcpy(nextArg, data);

      nextEleArgs[numEleArgs] = nextArg;

      if (eleArgs != 0)
	delete [] eleArgs;

      eleArgs = nextEleArgs;
      numEleArgs++;
    }
  }


  //
  // create the output handler
  //
  
  if (outMode == 0) 
    theOutputStream = new StandardStream();
  if (outMode == 1)
    theOutputStream = new DataFileStream(outputName);
  else if (outMode == 2) 
    theOutputStream = new XmlFileStream(outputName);
  else if (outMode == 3) 
    theOutputStream = new BinaryFileStream(outputName);

  //
  // create the recorder
  //

  Recorder *theRecorder = new ElementSumRecorder(eleID,
						 eleArgs,
						 numEleArgs,
						 echoTime, 
						 theOutputStream);

  // return it

  return theRecorder;
}



ElementSumRecorder::ElementSumRecorder()
:Recorder(-1),
 numEle(0), eleID(0), theResponses(0), 
 theDomain(0), theOutputHandler(0),
 echoTimeFlag(true), data(0), 
 responseArgs(0), numArgs(0)
{

}

ElementSumRecorder::ElementSumRecorder(const ID ele,
				       char **argv, 
				       int argc,
				       bool echoTime, 
				       OPS_Stream *theOutputHandler)
:Recorder(-1),
 numEle(0), eleID(0), theResponses(0), 
 theDomain(0), theOutputHandler(theOutputHandler),
 echoTimeFlag(echoTime), data(0),
 responseArgs(0), numArgs(0)
{
  if (ele != 0) {
    numEle = ele.Size();
    eleID = new ID(ele);
    if (eleID == 0 || eleID->Size() != numEle)
      opserr << "ElementSumRecorder::ElementSumRecorder() - out of memory\n";
  } 

  //
  // create a copy of the response request
  //

  responseArgs = new char *[argc];
  if (responseArgs == 0) {
    opserr << "ElementSumRecorder::ElementSumRecorder() - out of memory\n";
    numEle = 0;
  }
  
  for (int i=0; i<argc; i++) {
    responseArgs[i] = new char[strlen(argv[i])+1];
    if (responseArgs[i] == 0) {
      delete [] responseArgs;
      opserr << "ElementSumRecorder::ElementSumRecorder() - out of memory\n";
      numEle = 0;
    }
    strcpy(responseArgs[i], argv[i]);
  }
  
  numArgs = argc;
}


ElementSumRecorder::~ElementSumRecorder()
{
  theOutputHandler->endTag(); // Data
  theOutputHandler->endTag(); // OpenSeesOutput

  if (theOutputHandler != 0)
    delete theOutputHandler;

  //
  // invoke the destructor on the response objects
  //

  if (eleID != 0)
    delete eleID;

  if (theResponses != 0) {
    for (int i = 0; i < numEle; i++)
      delete theResponses[i];
    delete [] theResponses;
  }

  if (data != 0)
    delete data;
  
  // 
  // invoke destructor on response args
  //
  
  for (int i=0; i<numArgs; i++)
    delete [] responseArgs[i];
  delete [] responseArgs;

}


int 
ElementSumRecorder::record(int commitTag, double timeStamp)
{
  // 
  // check that initialization has been done
  //

  int result = 0;
  data->Zero();

  if (echoTimeFlag == true) 
    (*data)(0) = timeStamp;
    
  //
  // for each element if responses exist, put them in response vector
  //

  for (int i=0; i< numEle; i++) {
    if (theResponses[i] != 0) {
      // ask the element for the reponse
      int res;
      if (( res = theResponses[i]->getResponse()) < 0)
	result += res;
      else {
	Information &eleInfo = theResponses[i]->getInformation();
	const Vector &eleData = eleInfo.getData();
	int loc;
	if (echoTimeFlag == true) 
	  loc = 1;
	else 
	  loc = 0;
	for (int j=0; j<eleData.Size(); j++)
	  (*data)(loc++) += eleData(j);
      }
    } 
  }

  //
  // send the response vector to the output handler for o/p
  //
  
  theOutputHandler->write(*data);

  // succesfull completion - return 0
  return result;
}

int
ElementSumRecorder::restart(void)
{
  if (data != 0)
    data->Zero();
  return 0;
}


int 
ElementSumRecorder::setDomain(Domain &theDom)
{
  theDomain = &theDom;
  
  if (theDomain == 0)
    return -1;

  if (theResponses != 0) {
    delete theResponses;
    theResponses = 0;
  }

  theOutputHandler->tag("OpenSeesOutput");

  int numDbColumns = 0;
  if (echoTimeFlag == true) {
    theOutputHandler->tag("TimeOutput");
    theOutputHandler->tag("ResponseType", "time");
    theOutputHandler->endTag(); // TimeOutput
  }

  // Set the response objects:
  //   1. create an array of pointers for them
  //   2. iterate over the elements invoking setResponse() to get the new objects & determine size of data
  //      all elements must of course have same response size, if not error
  //

  // allocate memory for Reponses & set to 0
  theResponses = new Response *[numEle];
  if (theResponses == 0) {
    opserr << "ElementSumRecorder::initialize() - out of memory\n";
    return -1;
  }
  
  for (int k=0; k<numEle; k++)
    theResponses[k] = 0;
  
  // loop over ele & set Reponses
  for (int i=0; i<numEle; i++) {
    Element *theEle = theDomain->getElement((*eleID)(i));
    if (theEle == 0) {
      theResponses[i] = 0;
    } else {
      theResponses[i] = theEle->setResponse((const char **)responseArgs, numArgs, *theOutputHandler);
      if (theResponses[i] != 0) {
	// from the response type determine no of cols for each
	Information &eleInfo = theResponses[i]->getInformation();
	const Vector &eleData = eleInfo.getData();
	if (i == 0)
	  numDbColumns = eleData.Size();
	else {
	  if (numDbColumns != eleData.Size()) {
	    opserr << numDbColumns << " " << eleData.Size() << " " << (*eleID)(i) << endln;
	    opserr << "ElementSumRecorder::initialize() - elements do not all have same response size\n";
	    return -1;
	  }
	}
      }
    }
  }

  if (echoTimeFlag == true)
    numDbColumns ++;


  // create the vector to hold the data
  data = new Vector(numDbColumns);
  
  if (data == 0) {
    opserr << "ElementSumRecorder::initialize() - out of memory\n";
    return -1;
  }
  
  theOutputHandler->tag("Data");
  
  return 0;
}

int
ElementSumRecorder::sendSelf(int commitTag, Channel &theChannel)
{
  if (theChannel.isDatastore() == 1) {
    opserr << "ElementSumRecorder::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  //
  // into an ID, place & send (*eleID) size, numArgs and length of all responseArgs
  //

  static ID idData(6);
  if (eleID != 0)
    idData(0) = eleID->Size();
  else
    idData(0) = 0;

  idData(1) = numArgs;

  int msgLength = 0;
  for (int i=0; i<numArgs; i++) 
    msgLength += strlen(responseArgs[i])+1;

  idData(2) = msgLength;

  if (theOutputHandler != 0) {
    idData(3) = theOutputHandler->getClassTag();
  } else 
    idData(3) = 0;

  if (echoTimeFlag == true)
    idData(4) = 1;
  else
    idData(4) = 0;


  idData(5) = this->getTag();

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "ElementSumRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }

  //
  // send the eleID
  //

  if (eleID != 0)
    if (theChannel.sendID(0, commitTag, *eleID) < 0) {
      opserr << "ElementSumRecorder::sendSelf() - failed to send idData\n";
      return -1;
    }

  //
  // create a single char array holding all strings
  //    will use string terminating character to differentiate strings on other side
  //

  if (msgLength ==  0) {
    opserr << "ElementSumRecorder::sendSelf() - no data to send!!\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "ElementSumRecorder::sendSelf() - out of memory\n";
    return -1;
  }

  char *currentLoc = allResponseArgs;
  for (int j=0; j<numArgs; j++) {
    strcpy(currentLoc, responseArgs[j]);
    currentLoc += strlen(responseArgs[j]);
    currentLoc++;
  }

  //
  // send this single char array
  //


  Message theMessage(allResponseArgs, msgLength);
  if (theChannel.sendMsg(0, commitTag, theMessage) < 0) {
    opserr << "ElementSumRecorder::sendSelf() - failed to send message\n";
    return -1;
  }

  //
  // invoke sendSelf() on the output handler
  //

  if (theOutputHandler == 0 || theOutputHandler->sendSelf(commitTag, theChannel) < 0) {
    opserr << "ElementSumRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  //
  // clean up & return success
  //

  delete [] allResponseArgs;
  return 0;
}

int 
ElementSumRecorder::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  if (theChannel.isDatastore() == 1) {
    opserr << "ElementSumRecorder::recvSelf() - does not recv data to a datastore\n";
    return -1;
  }

  if (responseArgs != 0) {
    for (int i=0; i<numArgs; i++)
      delete [] responseArgs[i];
  
    delete [] responseArgs;
  }

  //
  // into an ID of size 2 recv eleID size and length of all responseArgs
  //

  static ID idData(6);
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "ElementSumRecorder::recvSelf() - failed to recv idData\n";
    return -1;
  }

  int eleSize = idData(0);
  numArgs = idData(1);
  int msgLength = idData(2);

  this->setTag(idData(5));

  if (idData(4) == 1)
    echoTimeFlag = true;
  else
    echoTimeFlag = false;    

  numEle = eleSize;

  //
  // resize & recv the eleID
  //

  if (eleSize != 0) {
    eleID = new ID(eleSize);
    if (eleID == 0) {
      opserr << "ElementSumRecorder::recvSelf() - failed to recv idData\n";
      return -1;
    }
    if (theChannel.recvID(0, commitTag, *eleID) < 0) {
      opserr << "ElementSumRecorder::recvSelf() - failed to recv idData\n";
      return -1;
    }
  }

  //
  // recv the single char array of element response args
  //

  if (msgLength == 0) {
    opserr << "ElementSumRecorder::recvSelf() - 0 sized string for responses\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "ElementSumRecorder::recvSelf() - out of memory\n";
    return -1;
  }

  Message theMessage(allResponseArgs, msgLength);
  if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
    opserr << "ElementSumRecorder::recvSelf() - failed to recv message\n";
    return -1;
  }

  //
  // now break this single array into many
  // 

  responseArgs = new char *[numArgs];
  if (responseArgs == 0) {
    opserr << "ElementSumRecorder::recvSelf() - out of memory\n";
    return -1;
  }

  char *currentLoc = allResponseArgs;
  for (int j=0; j<numArgs; j++) {

    int argLength = strlen(currentLoc)+1;

    responseArgs[j] = new char[argLength];
    if (responseArgs[j] == 0) {
      opserr << "ElementSumRecorder::recvSelf() - out of memory\n";
      return -1;
    }

    strcpy(responseArgs[j], currentLoc);
    currentLoc += argLength;
  }

  //
  // create a new handler object and invoke recvSelf() on it
  //

  if (theOutputHandler != 0)
    delete theOutputHandler;

  theOutputHandler = theBroker.getPtrNewStream(idData(3));
  if (theOutputHandler == 0) {
    opserr << "NodeRecorder::sendSelf() - failed to get a data output handler\n";
    return -1;
  }

  if (theOutputHandler->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "NodeRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  //
  // clean up & return success
  //

  delete [] allResponseArgs;
  return 0;
}


