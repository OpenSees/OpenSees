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
                                                                        
// $Revision: 1.34 $
// $Date: 2009-04-30 23:25:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/NormElementRecorder.cpp,v $
                                                                        
// Written: fmk 
// Created: 09/99
//
// Description: This file contains the class implementatation of NormElementRecorder.
//
// What: "@(#) NormElementRecorder.C, revA"

#include <NormElementRecorder.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <string.h>
#include <Response.h>
#include <FE_Datastore.h>
#include <OPS_Globals.h>
#include <Message.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <math.h>


NormElementRecorder::NormElementRecorder()
:Recorder(RECORDER_TAGS_NormElementRecorder),
 numEle(0), numDOF(0), eleID(0), dof(0), theResponses(0), 
 theDomain(0), theOutputHandler(0),
 echoTimeFlag(true), deltaT(0), nextTimeStampToRecord(0.0), data(0), 
 initializationDone(false), responseArgs(0), numArgs(0), addColumnInfo(0)
{

}

NormElementRecorder::NormElementRecorder(const ID *ele,
					 const char **argv, 
					 int argc,
					 bool echoTime, 
					 Domain &theDom, 
					 OPS_Stream &theOutputHandler,
					 double dT,
					 const ID *indices)
:Recorder(RECORDER_TAGS_NormElementRecorder),
 numEle(0), numDOF(0), eleID(0), dof(0), theResponses(0), 
 theDomain(&theDom), theOutputHandler(&theOutputHandler),
 echoTimeFlag(echoTime), deltaT(dT), nextTimeStampToRecord(0.0), data(0),
 initializationDone(false), responseArgs(0), numArgs(0), addColumnInfo(0)
{

  if (ele != 0) {
    numEle = ele->Size();
    eleID = new ID(*ele);
    if (eleID == 0 || eleID->Size() != numEle)
      opserr << "NormElementRecorder::NormElementRecorder() - out of memory\n";
  } 

  if (indices != 0) {
    dof = new ID(*indices);
    numDOF = dof->Size();
  }

  //
  // create a copy of the response request
  //

  responseArgs = new char *[argc];
  if (responseArgs == 0) {
    opserr << "NormElementRecorder::NormElementRecorder() - out of memory\n";
    numEle = 0;
  }
  
  for (int i=0; i<argc; i++) {
    responseArgs[i] = new char[strlen(argv[i])+1];
    if (responseArgs[i] == 0) {
      delete [] responseArgs;
      opserr << "NormElementRecorder::NormElementRecorder() - out of memory\n";
      numEle = 0;
    }
    strcpy(responseArgs[i], argv[i]);
  }
  
  numArgs = argc;
}


NormElementRecorder::~NormElementRecorder()
{
  if (theOutputHandler != 0) {
    theOutputHandler->endTag(); // Data
    delete theOutputHandler;
  }

  //
  // invoke the destructor on the response objects
  //

  if (eleID != 0)
    delete eleID;

  if (dof != 0)
    delete dof;

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
NormElementRecorder::record(int commitTag, double timeStamp)
{
  // 
  // check that initialization has been done
  //

  if (initializationDone == false) {
    if (this->initialize() != 0) {
      opserr << "NormElementRecorder::record() - failed to initialize\n";
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
	  double normV = 0.0;
	  if (numDOF == 0) {
	    for (int j=0; j<eleData.Size(); j++)
	      normV += eleData(j)*eleData(j);
	  } else {
	    int dataSize = data->Size();
	    for (int j=0; j<numDOF; j++) {
	      int index = (*dof)(j);
	      if (index >= 0 && index < dataSize)
		normV += eleData(index) * eleData(index);		
	    }
	  }
	  (*data)(loc++) = sqrt(normV);
	}
      } 
    }

    //
    // send the response vector to the output handler for o/p
    //
    theOutputHandler->write(*data);
  }
  
  // succesfull completion - return 0
  return result;
}

int
NormElementRecorder::restart(void)
{
  if (data != 0)
    data->Zero();
  return 0;
}


int 
NormElementRecorder::setDomain(Domain &theDom)
{
  theDomain = &theDom;
  return 0;
}

int
NormElementRecorder::sendSelf(int commitTag, Channel &theChannel)
{
  addColumnInfo = 1;

  if (theChannel.isDatastore() == 1) {
    opserr << "NormElementRecorder::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  initializationDone = false;
  //
  // into an ID, place & send (*eleID) size, numArgs and length of all responseArgs
  //

  static ID idData(7);
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
  idData(6) = numDOF;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "NormElementRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }

  static Vector dData(2);
  dData(0) = deltaT;
  dData(1) = nextTimeStampToRecord;
  if (theChannel.sendVector(0, commitTag, dData) < 0) {
    opserr << "NormElementRecorder::sendSelf() - failed to send dData\n";
    return -1;
  }
  
  //
  // send the eleID
  //

  if (eleID != 0)
    if (theChannel.sendID(0, commitTag, *eleID) < 0) {
      opserr << "NormElementRecorder::sendSelf() - failed to send idData\n";
      return -1;
    }

  // send dof
  if (dof != 0)
    if (theChannel.sendID(0, commitTag, *dof) < 0) {
      opserr << "ElementRecorder::sendSelf() - failed to send dof\n";
      return -1;
    }

  //
  // create a single char array holding all strings
  //    will use string terminating character to differentiate strings on other side
  //

  if (msgLength ==  0) {
    opserr << "NormElementRecorder::sendSelf() - no data to send!!\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "NormElementRecorder::sendSelf() - out of memory\n";
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
    opserr << "NormElementRecorder::sendSelf() - failed to send message\n";
    return -1;
  }

  //
  // invoke sendSelf() on the output handler
  //

  if (theOutputHandler == 0 || theOutputHandler->sendSelf(commitTag, theChannel) < 0) {
    opserr << "NormElementRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  //
  // clean up & return success
  //

  delete [] allResponseArgs;
  return 0;
}

int 
NormElementRecorder::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  addColumnInfo = 1;

  if (theChannel.isDatastore() == 1) {
    opserr << "NormElementRecorder::recvSelf() - does not recv data to a datastore\n";
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

  static ID idData(7);
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "NormElementRecorder::recvSelf() - failed to recv idData\n";
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
  numDOF = idData(6);

  static Vector dData(2);
  if (theChannel.recvVector(0, commitTag, dData) < 0) {
    opserr << "NormElementRecorder::sendSelf() - failed to send dData\n";
    return -1;
  }
  deltaT = dData(0);
  nextTimeStampToRecord = dData(1);
  //
  // resize & recv the eleID
  //

  if (eleSize != 0) {
    eleID = new ID(eleSize);
    if (eleID == 0) {
      opserr << "NormElementRecorder::recvSelf() - failed to recv idData\n";
      return -1;
    }
    if (theChannel.recvID(0, commitTag, *eleID) < 0) {
      opserr << "NormElementRecorder::recvSelf() - failed to recv idData\n";
      return -1;
    }
  }


  //
  // resize & recv the dof
  //

  if (numDOF != 0) {
    dof = new ID(numDOF);
    if (dof == 0) {
      opserr << "ElementRecorder::recvSelf() - failed to create dof\n";
      return -1;
    }
    if (theChannel.recvID(0, commitTag, *dof) < 0) {
      opserr << "ElementRecorder::recvSelf() - failed to recv dof\n";
      return -1;
    }
  }

  //
  // recv the single char array of element response args
  //

  if (msgLength == 0) {
    opserr << "NormElementRecorder::recvSelf() - 0 sized string for responses\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "NormElementRecorder::recvSelf() - out of memory\n";
    return -1;
  }

  Message theMessage(allResponseArgs, msgLength);
  if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
    opserr << "NormElementRecorder::recvSelf() - failed to recv message\n";
    return -1;
  }

  //
  // now break this single array into many
  // 

  responseArgs = new char *[numArgs];
  if (responseArgs == 0) {
    opserr << "NormElementRecorder::recvSelf() - out of memory\n";
    return -1;
  }

  char *currentLoc = allResponseArgs;
  for (int j=0; j<numArgs; j++) {

    int argLength = strlen(currentLoc)+1;

    responseArgs[j] = new char[argLength];
    if (responseArgs[j] == 0) {
      opserr << "NormElementRecorder::recvSelf() - out of memory\n";
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

int 
NormElementRecorder::initialize(void)
{
  if (theDomain == 0)
    return 0;

  if (theResponses != 0) {
    for (int i = 0; i < numEle; i++)
      delete theResponses[i];
    delete [] theResponses;
  }

  int numDbColumns = 0;

  // Set the response objects:
  //   1. create an array of pointers for them
  //   2. iterate over the elements invoking setResponse() to get the new objects & determine size of data
  //

  int i =0;
  ID xmlOrder(0,64);
  ID responseOrder(0,64);

  if (eleID != 0) {

    //
    // if we have an eleID we know Reponse size so allocate Response holder & loop over & ask each element
    //

    int eleCount = 0;
    int responseCount = 0;

    if (echoTimeFlag == true && addColumnInfo == 1) {
      xmlOrder[0] = 0;
      responseOrder[0] = 0;
      eleCount = 1;
      responseCount =1;
    }

    // loop over ele & set Reponses
    for (i=0; i<numEle; i++) {
      Element *theEle = theDomain->getElement((*eleID)(i));
      if (theEle != 0) {
	xmlOrder[eleCount] = i+1;
	eleCount++;
      }
    }

    theOutputHandler->setOrder(xmlOrder);

    //
    // do time
    //

    if (echoTimeFlag == true) {
      theOutputHandler->tag("TimeOutput");
      theOutputHandler->tag("ResponseType", "time");
      theOutputHandler->endTag(); // TimeOutput
      numDbColumns += 1;
    }

    //
    // if we have an eleID we know Reponse size so allocate Response holder & loop over & ask each element
    //

    // allocate memory for Reponses & set to 0
    theResponses = new Response *[numEle];
    if (theResponses == 0) {
      opserr << "NormElementRecorder::initialize() - out of memory\n";
      return -1;
    }

    for (int k=0; k<numEle; k++)
      theResponses[k] = 0;

    // loop over ele & set Reponses
    for (i=0; i<numEle; i++) {
      Element *theEle = theDomain->getElement((*eleID)(i));
      if (theEle == 0) {
	theResponses[i] = 0;
      } else {
	theResponses[i] = theEle->setResponse((const char **)responseArgs, numArgs, *theOutputHandler);
	if (theResponses[i] != 0) {
	  // from the response type determine no of cols for each
	  Information &eleInfo = theResponses[i]->getInformation();
	  const Vector &eleData = eleInfo.getData();
	  int dataSize = eleData.Size();
	  if (dataSize > 0) {
	    numDbColumns += 1;
	    if (addColumnInfo == 1) {
	      responseOrder[responseCount++] = i+1;
	    }
	  }
	}
      }
    }

    theOutputHandler->setOrder(responseOrder);

  } else {

    if (echoTimeFlag == true) {
      theOutputHandler->tag("TimeOutput");
      theOutputHandler->tag("ResponseType", "time");
      theOutputHandler->endTag(); // TimeOutput
      numDbColumns += 1;
    }

    //
    // if no eleID we don't know response size so make initial guess & loop over & ask ele
    // if guess to small, we enlarge
    //

    // initial size & allocation
    int numResponse = 0;
    numEle = 12;
    theResponses = new Response *[numEle];

    if (theResponses == 0) {
      opserr << "NormElementRecorder::initialize() - out of memory\n";
      return -1;
    }

    for (int k=0; k<numEle; k++)
      theResponses[k] = 0;

    // loop over ele & set Reponses
    ElementIter &theElements = theDomain->getElements();
    Element *theEle;

    while ((theEle = theElements()) != 0) {
      Response *theResponse = theEle->setResponse((const char **)responseArgs, numArgs, *theOutputHandler);
      if (theResponse != 0) {
	if (numResponse == numEle) {
	  // Why is this created locally and not used? -- MHS
	  Response **theNextResponses = new Response *[numEle*2];
	  if (theNextResponses != 0) {
	    for (int i=0; i<numEle; i++)
	      theNextResponses[i] = theResponses[i];
	    for (int j=numEle; j<2*numEle; j++)
	      theNextResponses[j] = 0;
	  }
	  numEle = 2*numEle;
	  delete [] theNextResponses;
	}
	theResponses[numResponse] = theResponse;

	// from the response type determine no of cols for each
	Information &eleInfo = theResponses[numResponse]->getInformation();
	const Vector &eleData = eleInfo.getData();
	numDbColumns += 1;

	numResponse++;

      }
    }
    numEle = numResponse;
  }

  // create the vector to hold the data
  data = new Vector(numDbColumns);

  if (data == 0) {
    opserr << "NormElementRecorder::initialize() - out of memory\n";
    return -1;
  }
  
  theOutputHandler->tag("Data");
  initializationDone = true;

  return 0;
}
