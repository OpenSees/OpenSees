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
                                                                        
// $Revision: 1.15 $
// $Date: 2006-01-18 19:43:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/EnvelopeElementRecorder.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the class implementatation of 
// EnvelopeElementRecorder.
//
// What: "@(#) EnvelopeElementRecorder.C, revA"

#include <EnvelopeElementRecorder.h>
#include <Domain.h>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <string.h>
#include <Response.h>
#include <FE_Datastore.h>
#include <Information.h>
#include <DataOutputHandler.h>

#include <Message.h>
#include <FEM_ObjectBroker.h>

EnvelopeElementRecorder::EnvelopeElementRecorder()
:Recorder(RECORDER_TAGS_EnvelopeElementRecorder),
 numEle(0), eleID(0), theResponses(0), theDomain(0),
 theHandler(0), deltaT(0), nextTimeStampToRecord(0.0), 
 data(0), currentData(0), first(true),
 initializationDone(false), responseArgs(0), numArgs(0), echoTimeFlag(false)
{

}


EnvelopeElementRecorder::EnvelopeElementRecorder(const ID &theEleID, 
						 const char **argv, 
						 int argc,
						 Domain &theDom, 
						 DataOutputHandler &theOutputHandler,
						 double dT, bool echoTime)
:Recorder(RECORDER_TAGS_EnvelopeElementRecorder),
 numEle(theEleID.Size()), eleID(theEleID), theResponses(0), theDomain(&theDom),
 theHandler(&theOutputHandler), deltaT(dT), nextTimeStampToRecord(0.0), 
 data(0), currentData(0), first(true),
 initializationDone(false), responseArgs(0), numArgs(0), echoTimeFlag(echoTime)
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

EnvelopeElementRecorder::~EnvelopeElementRecorder()
{
  //
  // write the data
  //
  
  if (theHandler != 0 && currentData != 0) {
    for (int i=0; i<3; i++) {
      int size = currentData->Size();
      for (int j=0; j<size; j++)
	(*currentData)(j) = (*data)(i,j);
      theHandler->write(*currentData);
    }
  }

  //
  // clean up the memory
  //

  if (theResponses != 0) {
    for (int i = 0; i < numEle; i++) 
      if (theResponses[i] != 0)
	delete theResponses[i];

    delete [] theResponses;
  }

  
  if (data != 0)
    delete data;
  
  if (currentData != 0)
    delete currentData;

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
EnvelopeElementRecorder::record(int commitTag, double timeStamp)
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
    // for each element do a getResponse() & put the result in current data
    for (int i=0; i< numEle; i++) {
      if (theResponses[i] != 0) {
	// ask the element for the reponse
	int res;
	if (( res = theResponses[i]->getResponse()) < 0)
	  result += res;
	else {
	  // from the response determine no of cols for each
	  Information &eleInfo = theResponses[i]->getInformation();
	  const Vector &eleData = eleInfo.getData();
	  for (int j=0; j<eleData.Size(); j++) 
	    (*currentData)(loc++) = eleData(j);
	}
      }
    }

    if (echoTimeFlag == false) {
      // check if max or min
      // check if currentData modifies the saved data
      int size = currentData->Size();
      if (first == true) {
	for (int i=0; i<size; i++) {
	  (*data)(0,i) = (*currentData)(i);
	  (*data)(1,i) = (*currentData)(i);
	  (*data)(2,i) = fabs((*currentData)(i));
	  first = false;
	} 
      } else {
	for (int i=0; i<size; i++) {
	  double value = (*currentData)(i);
	  if ((*data)(0,i) > value) {
	    (*data)(0,i) = value;
	    double absValue = fabs(value);
	    if ((*data)(2,i) < absValue) 
	      (*data)(2,i) = absValue;
	  } else if ((*data)(1,i) < value) {
	    (*data)(1,i) = value;
	    double absValue = fabs(value);
	    if ((*data)(2,i) < absValue) 
	      (*data)(2,i) = absValue;
	  }
	}
      }
    } else {
      // check if max or min
      // check if currentData modifies the saved data
      int size = currentData->Size();
      size /= 2;
      if (first == true) {
	for (int i=0; i<size; i++) {
	  (*data)(0,i*2) = timeStamp;
	  (*data)(1,i*2) = timeStamp;
	  (*data)(2,i*2) = timeStamp;
	  (*data)(0,i*2+1) = (*currentData)(i);
	  (*data)(1,i*2+1) = (*currentData)(i);
	  (*data)(2,i*2+1) = fabs((*currentData)(i));
	  first = false;
	} 
      } else {
	for (int i=0; i<size; i++) {
	  double value = (*currentData)(i);
	  if ((*data)(0,2*i+1) > value) {
	    (*data)(0,2*i) = timeStamp;
	    (*data)(0,2*i+1) = value;
	    double absValue = fabs(value);
	    if ((*data)(2,2*i+1) < absValue) {
	      (*data)(2,2*i) = timeStamp;
	      (*data)(2,2*i+1) = absValue;
	    }
	  } else if ((*data)(1,2*i+1) < value) {
	    (*data)(1,2*i) = timeStamp;
	    (*data)(1,2*i+1) = value;
	    double absValue = fabs(value);
	    if ((*data)(2,2*i+1) < absValue) {
	      (*data)(2,2*i) = timeStamp;
	      (*data)(2,2*i+1) = absValue;
	    }
	  }
	}
      }
    }
  }
  // succesfull completion - return 0
  return result;
}

int
EnvelopeElementRecorder::restart(void)
{
  data->Zero();
  first = true;
  return 0;
}

int 
EnvelopeElementRecorder::setDomain(Domain &theDom)
{
  theDomain = &theDom;
  return 0;
}

int
EnvelopeElementRecorder::sendSelf(int commitTag, Channel &theChannel)
{
  if (theChannel.isDatastore() == 1) {
    opserr << "EnvelopeElementRecorder::sendSelf() - does not send data to a datastore\n";
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
    opserr << "EnvelopeElementRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }

  //
  // send the eleID
  //

  if (eleID.Size() == 0 || theChannel.sendID(0, commitTag, eleID) < 0) {
    opserr << "EnvelopeElementRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }

  //
  // create a single char array holding all strings
  //    will use string terminating character to differentiate strings on other side
  //

  if (msgLength ==  0) {
    opserr << "EnvelopeElementRecorder::sendSelf() - no data to send!!\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "EnvelopeElementRecorder::sendSelf() - out of memory\n";
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
    opserr << "EnvelopeElementRecorder::sendSelf() - failed to send message\n";
    return -1;
  }

  //
  // invoke sendSelf() on the output handler
  //

  if (theHandler == 0 || theHandler->sendSelf(commitTag, theChannel) < 0) {
    opserr << "EnvelopeElementRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  //
  // clean up & return success
  //

  delete [] allResponseArgs;
  return 0;
}


int 
EnvelopeElementRecorder::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  if (theChannel.isDatastore() == 1) {
    opserr << "EnvelopeElementRecorder::recvSelf() - does not recv data to a datastore\n";
    return -1;
  }

  //
  // into an ID of size 2 recv eleID size and length of all responseArgs
  //

  static ID idData(4);
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "EnvelopeElementRecorder::recvSelf() - failed to recv idData\n";
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
    opserr << "EnvelopeElementRecorder::recvSelf() - 0 sized eleID\n";
    return -1;
  }

  int *eleData = new int [eleSize];
  eleID.setData(eleData, eleSize, true);
  if (theChannel.recvID(0, commitTag, eleID) < 0) {
    opserr << "EnvelopeElementRecorder::recvSelf() - failed to recv idData\n";
    return -1;
  }

  //
  // recv the single char array of element response args
  //

  if (msgLength == 0) {
    opserr << "EnvelopeElementRecorder::recvSelf() - 0 sized string for responses\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "EnvelopeElementRecorder::recvSelf() - out of memory\n";
    return -1;
  }

  Message theMessage(allResponseArgs, msgLength);
  if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
    opserr << "EnvelopeElementRecorder::recvSelf() - failed to recv message\n";
    return -1;
  }

  //
  // now break this single array into many
  // 

  responseArgs = new char *[numArgs];
  if (responseArgs == 0) {
    opserr << "EnvelopeElementRecorder::recvSelf() - out of memory\n";
    return -1;
  }

  char *currentLoc = allResponseArgs;
  for (int j=0; j<numArgs; j++) {
    int argLength = strlen(currentLoc)+1;
    responseArgs[j] = new char[argLength];
    if (responseArgs[j] == 0) {
      opserr << "EnvelopeElementRecorder::recvSelf() - out of memory\n";
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
EnvelopeElementRecorder::initialize(void) 
{
  if (numEle == 0 || theDomain == 0)
    return 0;

  // Set the response objects:
  //   1. create an array of pointers for them
  //   2. iterate over the elements invoking setResponse() to get the new objects & determine size of data
  //

  theResponses = new Response *[numEle];
  for (int k=0; k<numEle; k++)
    theResponses[k] = 0;

  Information eleInfo(1.0);
  int numDbColumns = 0;

  for (int ii=0; ii<numEle; ii++) {
    Element *theEle = theDomain->getElement(eleID(ii));
    if (theEle == 0) {
      theResponses[ii] = 0;
    } else {
      theResponses[ii] = theEle->setResponse((const char **)responseArgs, numArgs, eleInfo);
      if (theResponses[ii] != 0) {
	// from the response type determine no of cols for each      
	Information &eleInfo = theResponses[ii]->getInformation();
	const Vector &eleData = eleInfo.getData();
	numDbColumns += eleData.Size();
      }
    }
  }

  //
  // create the matrix & vector that holds the data
  //

  if (echoTimeFlag == true)
    numDbColumns *= 2;

  data = new Matrix(3, numDbColumns);
  currentData = new Vector(numDbColumns);
  if (data == 0 || currentData == 0) {
    opserr << "EnvelopeElementRecorder::EnvelopeElementRecorder() - out of memory\n";
    exit(-1);
  }

  //
  // now create the columns strings for the database
  // for each element do a getResponse() & print the result
  //

  char **dbColumns = new char *[numDbColumns];
  static char aColumn[1012]; // assumes a column name will not be longer than 256 characters
  
  int lengthString = 0;
  for (int l=0; l<numArgs; l++)
    lengthString += strlen(responseArgs[l])+1;
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

  int counter = 0;
  for (int i=0; i<eleID.Size(); i++) {
    int eleTag = eleID(i);
    int numVariables = 0;
    if (theResponses[i]!= 0) {
      const Information &eleInfo = theResponses[i]->getInformation();
      
      if (eleInfo.theType == IntType || eleInfo.theType == DoubleType) {
	// create column heading for single data item for element
	if (echoTimeFlag == true) {
	  int lenColumn = strlen("time");
	  char *newColumn = new char[lenColumn+1];
	  strcpy(newColumn, "time");
	  dbColumns[counter] = newColumn;
	  counter++;
	}

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

      // create the column headings for multiple data for the element
      for (int j=1; j<=numVariables; j++) {
	if (echoTimeFlag == true) {
	  int lenColumn = strlen("time");
	  char *newColumn = new char[lenColumn+1];
	  strcpy(newColumn, "time");
	  dbColumns[counter] = newColumn;
	  counter++;
	}
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
  for (int kk=0; kk<numDbColumns; kk++) {
    char *data = dbColumns[kk];
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

  initializationDone = true;  
  return 0;
}
  
