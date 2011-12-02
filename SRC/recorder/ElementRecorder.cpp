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
                                                                        
// $Revision: 1.24 $
// $Date: 2005-12-06 00:21:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/ElementRecorder.cpp,v $
                                                                        
// Written: fmk 
// Created: 09/99
//
// Description: This file contains the class implementatation of ElementRecorder.
//
// What: "@(#) ElementRecorder.C, revA"

#include <ElementRecorder.h>
#include <Domain.h>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <string.h>
#include <Response.h>
#include <FE_Datastore.h>
#include <DataOutputHandler.h>
#include <OPS_Globals.h>
#include <Message.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

ElementRecorder::ElementRecorder()
:Recorder(RECORDER_TAGS_ElementRecorder),
 numEle(0), eleID(0), theResponses(0), 
 theDomain(0), theHandler(0),
 echoTimeFlag(true), deltaT(0), nextTimeStampToRecord(0.0), data(0), 
 initializationDone(false), responseArgs(0), numArgs(0)
{

}

ElementRecorder::ElementRecorder(const ID &ele,
				 const char **argv, 
				 int argc,
				 bool echoTime, 
				 Domain &theDom, 
				 DataOutputHandler &theOutputHandler,
				 double dT)
:Recorder(RECORDER_TAGS_ElementRecorder),
 numEle(ele.Size()), eleID(ele), theResponses(0), 
 theDomain(&theDom), theHandler(&theOutputHandler),
 echoTimeFlag(echoTime), deltaT(dT), nextTimeStampToRecord(0.0), data(0),
 initializationDone(false), responseArgs(0), numArgs(0)
{
  //
  // create a copy of the response request
  //

  responseArgs = new char *[argc];
  if (responseArgs == 0) {
    opserr << "ElementRecorder::ElementRecorder() - out of memory\n";
    numEle = 0;
  }
  
  for (int i=0; i<argc; i++) {
    responseArgs[i] = new char[strlen(argv[i])+1];
    if (responseArgs[i] == 0) {
      delete [] responseArgs;
      opserr << "ElementRecorder::ElementRecorder() - out of memory\n";
      numEle = 0;
    }
    strcpy(responseArgs[i], argv[i]);
  }
  
  numArgs = argc;
}


ElementRecorder::~ElementRecorder()
{
  //
  // invoke the destructor on the response objects
  //

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


  if (theHandler != 0)
    delete theHandler;
}


int 
ElementRecorder::record(int commitTag, double timeStamp)
{
  // 
  // check that initialization has been done
  //

  if (initializationDone == false) {
    if (this->initialize() != 0) {
      opserr << "ElementRecorder::record() - failed to initialize\n";
      return -1;
    }
  }

  int result = 0;
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {

    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;

    int loc = 0;
    if (echoTimeFlag == true) 
      (*data)(loc++) = timeStamp;
    
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
	  for (int j=0; j<eleData.Size(); j++)
	    (*data)(loc++) = eleData(j);
	}
      } 
    }

    //
    // send the response vector to the output handler for o/p
    //

    theHandler->write(*data);
  }
  
  // succesfull completion - return 0
  return result;
}

int
ElementRecorder::restart(void)
{
  if (data != 0)
    data->Zero();
  return 0;
}


int 
ElementRecorder::setDomain(Domain &theDom)
{
  theDomain = &theDom;
  return 0;
}

int
ElementRecorder::sendSelf(int commitTag, Channel &theChannel)
{
  if (theChannel.isDatastore() == 0) {
    opserr << "ElementRecorder::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  //
  // into an ID, place & send eleID size, numArgs and length of all responseArgs
  //

  static ID idData(4);
  idData(0) = eleID.Size();
  idData(1) = numArgs;

  int msgLength = 0;
  for (int i=0; i<numArgs; i++) 
    msgLength += strlen(responseArgs[i])+1;
  idData(2) = msgLength;

  if (theHandler != 0) {
    idData(3) = theHandler->getClassTag();
  } else 
    idData(3) = 0;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "ElementRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }

  //
  // send the eleID
  //

  if (eleID.Size() == 0 || theChannel.sendID(0, commitTag, eleID) < 0) {
    opserr << "ElementRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }

  //
  // create a single char array holding all strings
  //    will use string terminating character to differentiate strings on other side
  //

  if (msgLength ==  0) {
    opserr << "ElementRecorder::sendSelf() - no data to send!!\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "ElementRecorder::sendSelf() - out of memory\n";
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
    opserr << "ElementRecorder::sendSelf() - failed to send message\n";
    return -1;
  }

  //
  // invoke sendSelf() on the output handler
  //

  if (theHandler == 0 || theHandler->sendSelf(commitTag, theChannel) < 0) {
    opserr << "ElementRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  //
  // clean up & return success
  //

  delete [] allResponseArgs;
  return 0;
}

int 
ElementRecorder::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  if (theChannel.isDatastore() == 0) {
    opserr << "ElementRecorder::recvSelf() - does not recv data to a datastore\n";
    return -1;
  }

  //
  // into an ID of size 2 recv eleID size and length of all responseArgs
  //

  static ID idData(4);
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "ElementRecorder::recvSelf() - failed to recv idData\n";
    return -1;
  }

  int eleSize = idData(0);
  int numArgs = idData(1);
  int msgLength = idData(2);

  numEle = eleSize;

  //
  // resize & recv the eleID
  //

  if (eleSize == 0) {
    opserr << "ElementRecorder::recvSelf() - 0 sized eleID\n";
    return -1;
  }

  int *eleData = new int [eleSize];
  eleID.setData(eleData, eleSize, true);
  if (theChannel.recvID(0, commitTag, eleID) < 0) {
    opserr << "ElementRecorder::recvSelf() - failed to recv idData\n";
    return -1;
  }

  //
  // recv the single char array of element response args
  //

  if (msgLength == 0) {
    opserr << "ElementRecorder::recvSelf() - 0 sized string for responses\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "ElementRecorder::recvSelf() - out of memory\n";
    return -1;
  }

  Message theMessage(allResponseArgs, msgLength);
  if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
    opserr << "ElementRecorder::recvSelf() - failed to recv message\n";
    return -1;
  }

  //
  // now break this single array into many
  // 

  responseArgs = new char *[numArgs];
  if (responseArgs == 0) {
    opserr << "ElementRecorder::recvSelf() - out of memory\n";
    return -1;
  }

  char *currentLoc = allResponseArgs;
  for (int j=0; j<numArgs; j++) {
    int argLength = strlen(currentLoc)+1;
    responseArgs[j] = new char[argLength];
    if (responseArgs[j] == 0) {
      opserr << "ElementRecorder::recvSelf() - out of memory\n";
      return -1;
    }
    opserr << argLength << responseArgs[j] << endln;
    strcpy(responseArgs[j], currentLoc);
    currentLoc += argLength;
  }

  //
  // create a new handler object and invoke recvSelf() on it
  //

  if (theHandler != 0)
    delete theHandler;

  theHandler = theBroker.getPtrNewDataOutputHandler(idData(3));
  if (theHandler == 0) {
    opserr << "NodeRecorder::sendSelf() - failed to get a data output handler\n";
    return -1;
  }

  if (theHandler->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "NodeRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  //
  // clean up & return success
  //

  delete [] allResponseArgs;
  return 0;
}

int 
ElementRecorder::initialize(void)
{
  if (numEle == 0 || theDomain == 0)
    return 0;

  // Set the response objects:
  //   1. create an array of pointers for them
  //   2. iterate over the elements invoking setResponse() to get the new objects & determine size of data
  //

  int numDbColumns = 0;
  if (echoTimeFlag == true) 
    numDbColumns = 1;  // 1 for the pseudo-time

  theResponses = new Response *[numEle];
  for (int k=0; k<numEle; k++)
    theResponses[k] = 0;

  Information eleInfo(1.0);
  int i;
  for (i=0; i<numEle; i++) {
    Element *theEle = theDomain->getElement(eleID(i));
    if (theEle == 0) {
      theResponses[i] = 0;
    } else {
      theResponses[i] = theEle->setResponse((const char **)responseArgs, numArgs, eleInfo);
      if (theResponses[i] != 0) {
	// from the response type determine no of cols for each
	Information &eleInfo = theResponses[i]->getInformation();
	const Vector &eleData = eleInfo.getData();
	numDbColumns += eleData.Size();
      }
    }
  }

  //
  // now create the columns strings for the data description
  // for each element do a getResponse() 
  //
  char **dbColumns = new char *[numDbColumns];
  static char aColumn[1012]; // assumes a column name will not be longer than 256 characters
  char *newColumn = new char[5];

  int counter = 0;
  if (echoTimeFlag == true) {
    sprintf(newColumn, "%s","time");  
    dbColumns[0] = newColumn;
    counter = 1;
  }
  
  int lengthString = 0;
  for (i=0; i<numArgs; i++)
    lengthString += strlen(responseArgs[i])+1;
  char *dataToStore = new char[lengthString];
  lengthString = 0;
  for (int j=0; j<numArgs; j++) {
    int argLength = strlen(responseArgs[j]);
    strcpy(&dataToStore[lengthString], responseArgs[j]);
    if (j<(numArgs-1)) {
      lengthString += argLength;
      dataToStore[lengthString] = ' ';
      lengthString ++;
    } else
      lengthString += argLength+1;
  }
  
  for (i=0; i<eleID.Size(); i++) {
    int eleTag = eleID(i);
    int numVariables = 0;
    if (theResponses[i]!= 0) {
      const Information &eleInfo = theResponses[i]->getInformation();
      
      if (eleInfo.theType == IntType || eleInfo.theType == DoubleType) {
	// create column heading for single data item for element
	numVariables = 0;
	sprintf(aColumn, "Element%d_%s", eleTag, dataToStore);
	int lenColumn = strlen(aColumn);
	char *newColumn = new char[lenColumn+1];
	strcpy(newColumn, aColumn);
	dbColumns[counter] = newColumn;
	counter++;
      }
      
      else if (eleInfo.theType == VectorType) 
	numVariables = eleInfo.theVector->Size();
      else if (eleInfo.theType == IdType) 
	numVariables = eleInfo.theID->Size();
      else if (eleInfo.theType == MatrixType) 
	numVariables = eleInfo.theMatrix->noRows()* eleInfo.theMatrix->noCols();

      // create the column headings for multiple data for the element
      for (int j=1; j<=numVariables; j++) {
	sprintf(aColumn, "Element%d_%s_%d",eleTag, dataToStore, j);
	int lenColumn = strlen(aColumn);
	char *newColumn = new char[lenColumn+1];
	strcpy(newColumn, aColumn);
	dbColumns[counter] = newColumn;
	counter++;
      }
    }
  }

  // replace spaces with undescore for tables
  for (i=0; i<numDbColumns; i++) {
    char *data = dbColumns[i];
    int length = strlen(data);
    for (int j=0; j<length; j++)
      if (data[j] == ' ') data[j]='_';
  }

  //
  // call open in the handler with the data description
  //

  theHandler->open(dbColumns, numDbColumns);

  //
  // clean up the data description
  //

  if (dbColumns != 0) {

    for (int i=0; i<numDbColumns; i++) 
      delete [] dbColumns[i];

      delete [] dbColumns;
  }
  
  delete [] dataToStore;

  // create the vector to hold the data
  data = new Vector(numDbColumns);

  if (data == 0) {
    opserr << "ElementRecorder::initialize() - out of memory\n";
    return -1;
  }

  

  initializationDone = true;
  return 0;
}
