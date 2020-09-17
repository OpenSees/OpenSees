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
                                                                        
// Written: fmk 
//
// Description: This file contains the class implementatation of 
// EnvelopeElementRecorderRMS.
//
// What: "@(#) EnvelopeElementRecorderRMS.C, revA"

#include <ElementRecorderRMS.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <Response.h>
#include <FE_Datastore.h>
#include <Information.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <MeshRegion.h>

#include <StandardStream.h>
#include <DataFileStream.h>
#include <DataFileStreamAdd.h>
#include <XmlFileStream.h>
#include <BinaryFileStream.h>
#include <DatabaseStream.h>
#include <TCP_Stream.h>

#include <elementAPI.h>

#include <string.h>
#include <stdlib.h>
#include <math.h>

void*
OPS_ElementRecorderRMS()
{
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "WARING: recorder Element ";
        opserr << "-ele <list elements> -file <fileName> -dT <dT> reponse";
        return 0;
    }

    const char** data = 0;
    int nargrem = 0;
    OPS_Stream *theOutputStream = 0;
    const char* filename = 0;

    const int STANDARD_STREAM = 0;
    const int DATA_STREAM = 1;
    const int XML_STREAM = 2;
    const int DATABASE_STREAM = 3;
    const int BINARY_STREAM = 4;
    const int DATA_STREAM_CSV = 5;
    const int TCP_STREAM = 6;
    const int DATA_STREAM_ADD = 7;

    int eMode = STANDARD_STREAM;

    bool echoTimeFlag = false;
    double dT = 0.0;
    bool doScientific = false;

    int precision = 6;

    bool closeOnWrite = false;

    const char *inetAddr = 0;
    int inetPort;

    ID elements(0, 6);
    ID dofs(0, 6);

    while (OPS_GetNumRemainingInputArgs() > 0) {

        const char* option = OPS_GetString();

        if (strcmp(option, "-time") == 0) {
            echoTimeFlag = true;
        }
        else if (strcmp(option, "-load") == 0) {
            echoTimeFlag = true;
        }
        else if (strcmp(option, "-scientific") == 0) {
            doScientific = true;
        }
        else if (strcmp(option, "-file") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                filename = OPS_GetString();
            }
            eMode = DATA_STREAM;
        }
        else if (strcmp(option, "-fileAdd") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                filename = OPS_GetString();
            }
            eMode = DATA_STREAM_ADD;
        }
        else if (strcmp(option, "-closeOnWrite") == 0) {
            closeOnWrite = true;
        }
        else if (strcmp(option, "-csv") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                filename = OPS_GetString();
            }
            eMode = DATA_STREAM_CSV;
        }
        else if (strcmp(option, "-tcp") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                inetAddr = OPS_GetString();
            }
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &inetPort) < 0) {
                    opserr << "WARNING: failed to read inetPort\n";
                    return 0;
                }
            }
            eMode = TCP_STREAM;
        }
        else if (strcmp(option, "-xml") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                filename = OPS_GetString();
            }
            eMode = XML_STREAM;
        }
        else if (strcmp(option, "-binary") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                filename = OPS_GetString();
            }
            eMode = BINARY_STREAM;
        }
        else if (strcmp(option, "-dT") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetDoubleInput(&num, &dT) < 0) {
                    opserr << "WARNING: failed to read dT\n";
                    return 0;
                }
            }
        }
        else if (strcmp(option, "-precision") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &precision) < 0) {
                    opserr << "WARNING: failed to read precision\n";
                    return 0;
                }
            }
        }
        else if (strcmp(option, "-ele") == 0) {
            int numEle = 0;
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                int el;
                if (OPS_GetIntInput(&num, &el) < 0) {
                    break;
                }
                elements[numEle++] = el;
            }
        }
        else if (strcmp(option, "-eleRange") == 0) {
            int start, end;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &start) < 0) {
                    opserr << "WARNING: failed to read start element\n";
                    return 0;
                }
            }
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &end) < 0) {
                    opserr << "WARNING: failed to read end element\n";
                    return 0;
                }
            }
            if (start > end) {
                int swap = end;
                end = start;
                start = swap;
            }
            int numEle = 0;
            for (int i = start; i <= end; i++)
                elements[numEle++] = i;
        }
        else if (strcmp(option, "-region") == 0) {
            int tag;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                if (OPS_GetIntInput(&num, &tag) < 0) {
                    opserr << "WARNING: failed to read region tag\n";
                    return 0;
                }
            }
            Domain *domain = OPS_GetDomain();
            MeshRegion *theRegion = domain->getRegion(tag);
            if (theRegion == 0) {
                opserr << "WARNING: region does not exist\n";
                return 0;
            }
            const ID &eleRegion = theRegion->getElements();
            int numEle = 0;
            for (int i = 0; i < eleRegion.Size(); i++)
                elements[numEle++] = eleRegion(i);
        }
        else if (strcmp(option, "-dof") == 0) {
            int numDOF = 0;
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int num = 1;
                int dof;
                if (OPS_GetIntInput(&num, &dof) < 0) {
                    OPS_ResetCurrentInputArg(-1);
                    break;
                }
                dofs[numDOF++] = dof - 1;
            }
        }
        else {
            // first unknown string then is assumed to start 
            // element response request
            nargrem = 1 + OPS_GetNumRemainingInputArgs();
            data = new const char *[nargrem];
            data[0] = option;
            for (int i = 1; i < nargrem; i++)
                data[i] = OPS_GetString();
        }
    }
    
    // data handler
    if (eMode == DATA_STREAM && filename != 0)
        theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
    else if (eMode == DATA_STREAM_ADD && filename != 0)
        theOutputStream = new DataFileStreamAdd(filename, OVERWRITE, 2, 0, closeOnWrite, precision, doScientific);
    else if (eMode == DATA_STREAM_CSV && filename != 0)
        theOutputStream = new DataFileStream(filename, OVERWRITE, 2, 1, closeOnWrite, precision, doScientific);
    else if (eMode == XML_STREAM && filename != 0)
        theOutputStream = new XmlFileStream(filename);
    //else if (eMode == DATABASE_STREAM && tableName != 0)
    //    theOutputStream = new DatabaseStream(theDatabase, tableName);
    else if (eMode == BINARY_STREAM && filename != 0)
        theOutputStream = new BinaryFileStream(filename);
    else if (eMode == TCP_STREAM && inetAddr != 0)
        theOutputStream = new TCP_Stream(inetPort, inetAddr);
    else
        theOutputStream = new StandardStream();

    theOutputStream->setPrecision(precision);

    Domain* domain = OPS_GetDomain();
    if (domain == 0)
        return 0;
    ElementRecorderRMS* recorder = new ElementRecorderRMS(&elements,
        data, nargrem, *domain, *theOutputStream, dT, &dofs);

    return recorder;
}


ElementRecorderRMS::ElementRecorderRMS()
:Recorder(RECORDER_TAGS_ElementRecorderRMS),
 numEle(0), numDOF(0), eleID(0), dof(0), theResponses(0), theDomain(0),
 theHandler(0), deltaT(0), nextTimeStampToRecord(0.0), 
 runningTotal(0), currentData(0), count(0), 
 initializationDone(false), responseArgs(0), numArgs(0), addColumnInfo(0)
{

}


ElementRecorderRMS::ElementRecorderRMS(const ID *ele, 
				       const char **argv, 
				       int argc,
				       Domain &theDom, 
				       OPS_Stream &theOutputHandler,
				       double dT, 
				       const ID *indexValues)
  :Recorder(RECORDER_TAGS_ElementRecorderRMS),
  numEle(0), eleID(0), numDOF(0), dof(0), theResponses(0), theDomain(&theDom),
  theHandler(&theOutputHandler), deltaT(dT), nextTimeStampToRecord(0.0), 
   runningTotal(0), currentData(0), count(0),
  initializationDone(false), responseArgs(0), numArgs(0), addColumnInfo(0)
{
  opserr << "ElementRMS:: constructor\n";

  if (ele != 0) {
    numEle = ele->Size();
    eleID = new ID(*ele);
    if (eleID == 0 || eleID->Size() != numEle)
      opserr << "ElementRecorder::ElementRecorder() - out of memory\n";
  }
  
  if (indexValues != 0) {
    dof = new ID(*indexValues);
    numDOF = dof->Size();
  }
  
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

ElementRecorderRMS::~ElementRecorderRMS()
{
  //
  // write the data
  //

  if (eleID != 0)
    delete eleID;

  opserr << "ElementRMS DESTRUCTOR\n";

  if (theHandler != 0 && currentData != 0) {

    theHandler->tag("Data"); // Data

    if (runningTotal != 0) {

      opserr << "ElementRMS DESTRUCTOR - runin\n" << runningTotal->Size() << "\n";

      for (int j=0; j<runningTotal->Size(); j++)
	if (count != 0) {
	  double value = (*runningTotal)(j);
	  (*runningTotal)(j) = sqrt(value/(count));
	}
      theHandler->write(*runningTotal);
    }

    theHandler->endTag(); // Data
  }

  if (theHandler != 0)
    delete theHandler;

  if (runningTotal != 0)
    delete runningTotal;
  
  if (currentData != 0)
    delete currentData;

  //
  // clean up the memory
  //

  if (theResponses != 0) {
    for (int i = 0; i < numEle; i++) 
      if (theResponses[i] != 0)
	delete theResponses[i];

    delete [] theResponses;
  }



  // 
  // invoke destructor on response args
  //

  for (int i=0; i<numArgs; i++)
    delete [] responseArgs[i];
  delete [] responseArgs;
}


int 
ElementRecorderRMS::record(int commitTag, double timeStamp)
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
	  //	  for (int j=0; j<eleData.Size(); j++) 
	  //	    (*currentData)(loc++) = eleData(j);
	  if (numDOF == 0) {
	    for (int j=0; j<eleData.Size(); j++)
	      (*currentData)(loc++) = eleData(j);
	  } else {
	    int dataSize = eleData.Size();
	    for (int j=0; j<numDOF; j++) {
	      int index = (*dof)(j);
	      if (index >= 0 && index < dataSize)
		(*currentData)(loc++) = eleData(index);		
	      else
		(*currentData)(loc++) = 0.0;		
	    }
	  }
	}
      }
    }

    // add data contribution to runningTotal
    count++;
    int sizeData = currentData->Size();
    for (int i=0; i<sizeData; i++) {
      double value = (*currentData)(i);
      (*runningTotal)(i) += value*value;
    } 
  }    

  // succesfull completion - return 0
  return result;
}

int
ElementRecorderRMS::restart(void)
{
  runningTotal->Zero();
  count = 0;
  return 0;
}

int 
ElementRecorderRMS::setDomain(Domain &theDom)
{
  theDomain = &theDom;
  return 0;
}

int
ElementRecorderRMS::sendSelf(int commitTag, Channel &theChannel)
{
  addColumnInfo = 1;

  if (theChannel.isDatastore() == 1) {
    opserr << "ElementRecorderRMS::sendSelf() - does not send data to a datastore\n";
    return -1;
  }
  initializationDone = false;
  //
  // into an ID, place & send eleID size, numArgs and length of all responseArgs
  //

  static ID idData(7);
  if (eleID != 0)
    idData(0) = eleID->Size();
  else
    idData(0) = 0;

  idData(1) = numArgs;

  idData(5) = this->getTag();
  idData(6) = numDOF;

  int msgLength = 0;
  for (int i=0; i<numArgs; i++) 
    msgLength += strlen(responseArgs[i])+1;
  idData(2) = msgLength;


  if (theHandler != 0) {
    idData(3) = theHandler->getClassTag();
  } else 
    idData(3) = 0;


  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "ElementRecorderRMS::sendSelf() - failed to send idData\n";
    return -1;
  }

  static Vector dData(1);
  dData(1) = deltaT;
  if (theChannel.sendVector(0, commitTag, dData) < 0) {
    opserr << "ElementRecorderRMS::sendSelf() - failed to send dData\n";
    return -1;
  }


  //
  // send the eleID
  //

  if (eleID != 0)
    if (theChannel.sendID(0, commitTag, *eleID) < 0) {
      opserr << "ElementRecorderRMS::sendSelf() - failed to send idData\n";
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
    opserr << "ElementRecorderRMS::sendSelf() - no data to send!!\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "ElementRecorderRMS::sendSelf() - out of memory\n";
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
    opserr << "ElementRecorderRMS::sendSelf() - failed to send message\n";
    return -1;
  }

  //
  // invoke sendSelf() on the output handler
  //

  if (theHandler == 0 || theHandler->sendSelf(commitTag, theChannel) < 0) {
    opserr << "ElementRecorderRMS::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  //
  // clean up & return success
  //

  delete [] allResponseArgs;
  return 0;
}


int 
ElementRecorderRMS::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  addColumnInfo = 1;

  if (theChannel.isDatastore() == 1) {
    opserr << "ElementRecorderRMS::recvSelf() - does not recv data to a datastore\n";
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
    opserr << "ElementRecorderRMS::recvSelf() - failed to recv idData\n";
    return -1;
  }

  int eleSize = idData(0);
  numArgs = idData(1);
  int msgLength = idData(2);
  numDOF = idData(6);

  this->setTag(idData(5));

  numEle = eleSize;

  static Vector dData(1);
  if (theChannel.recvVector(0, commitTag, dData) < 0) {
    opserr << "ElementRecorderRMS::recvSelf() - failed to recv dData\n";
    return -1;
  }
  deltaT = dData(1);


  //
  // resize & recv the eleID
  //

  if (eleSize != 0) {
    eleID = new ID(eleSize);
    if (eleID == 0) {
      opserr << "ElementRecorder::recvSelf() - failed to recv idData\n";
      return -1;
    }
    if (theChannel.recvID(0, commitTag, *eleID) < 0) {
      opserr << "ElementRecorder::recvSelf() - failed to recv idData\n";
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
    opserr << "ElementRecorderRMS::recvSelf() - 0 sized string for responses\n";
    return -1;
  }

  char *allResponseArgs = new char[msgLength];
  if (allResponseArgs == 0) {
    opserr << "ElementRecorderRMS::recvSelf() - out of memory\n";
    return -1;
  }

  Message theMessage(allResponseArgs, msgLength);
  if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
    opserr << "ElementRecorderRMS::recvSelf() - failed to recv message\n";
    return -1;
  }

  //
  // now break this single array into many
  // 

  responseArgs = new char *[numArgs];
  if (responseArgs == 0) {
    opserr << "ElementRecorderRMS::recvSelf() - out of memory\n";
    return -1;
  }

  char *currentLoc = allResponseArgs;
  for (int j=0; j<numArgs; j++) {

    int argLength = strlen(currentLoc)+1;

    responseArgs[j] = new char[argLength];
    if (responseArgs[j] == 0) {
      opserr << "ElementRecorderRMS::recvSelf() - out of memory\n";
      return -1;
    }

    strcpy(responseArgs[j], currentLoc);
    currentLoc += argLength;
  }

  //
  // create a new handler object and invoke recvSelf() on it
  //

  if (theHandler != 0)
    delete theHandler;

  theHandler = theBroker.getPtrNewStream(idData(3));
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
ElementRecorderRMS::initialize(void) 
{
  if (theDomain == 0)
    return 0;

  if (theResponses != 0) {
    for (int i = 0; i < numEle; i++)
      delete theResponses[i];
    delete [] theResponses;
  }

  int numDbColumns = 0;

  //
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

    // loop over ele & set Reponses
    for (i=0; i<numEle; i++) {
      Element *theEle = theDomain->getElement((*eleID)(i));
      if (theEle != 0) {
	xmlOrder[eleCount] = i+1;
	eleCount++;
      }
    }

    theHandler->setOrder(xmlOrder);

    //
    // if we have an eleID we know Reponse size so allocate Response holder & loop over & ask each element
    //

    // allocate memory for Reponses & set to 0
    theResponses = new Response *[numEle];
    if (theResponses == 0) {
      opserr << "ElementRecorder::initialize() - out of memory\n";
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
	opserr << *theEle;

	for (int a=0; a<numArgs; a++)
	  opserr << responseArgs[i];

	theResponses[i] = theEle->setResponse((const char **)responseArgs, numArgs, *theHandler);
	if (theResponses[i] != 0) {
	  opserr << "HAS RESPONSE\n";
	  // from the response type determine no of cols for each
	  Information &eleInfo = theResponses[i]->getInformation();
	  const Vector &eleData = eleInfo.getData();
	  int dataSize = eleData.Size();
	  if (numDOF == 0)
	    numDbColumns += dataSize;
	  else
	    numDbColumns += numDOF;

	  if (addColumnInfo == 1) {
	    if (numDOF == 0)
	      for (int j=0; j<dataSize; j++)
		responseOrder[responseCount++] = i+1;
	    else
	      for (int j=0; j<numDOF; j++)
		responseOrder[responseCount++] = i+1;
	  }
	} else { opserr << "NO RESPONSE\n";}
      } 
    }

    theHandler->setOrder(responseOrder);

  } else {

    //
    // if no eleID we don't know response size so make initial guess & loop over & ask ele
    // if guess to small, we enlarge
    //

    // initial size & allocation
    int numResponse = 0;
    numEle = 12;
    theResponses = new Response *[numEle];

    if (theResponses == 0) {
      opserr << "ElementRecorder::initialize() - out of memory\n";
      return -1;
    }

    for (int k=0; k<numEle; k++)
      theResponses[k] = 0;

    // loop over ele & set Reponses
    ElementIter &theElements = theDomain->getElements();
    Element *theEle;

    while ((theEle = theElements()) != 0) {
      Response *theResponse = theEle->setResponse((const char **)responseArgs, numArgs, *theHandler);
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
	numDbColumns += eleData.Size();

	numResponse++;

      }
    }
    numEle = numResponse;
  }

  //
  // create the matrix & vector that holds the data
  //

  opserr << numEle << " " << numDbColumns << "\n";

  runningTotal = new Vector(numDbColumns);
  currentData = new Vector(numDbColumns);
  if (runningTotal == 0 || currentData == 0) {
    opserr << "ElementRecorderRMS::ElementRecorderRMS() - out of memory\n";
    exit(-1);
  }
  runningTotal->Zero();
  currentData->Zero();

  initializationDone = true;  
  return 0;
}
  
//added by SAJalali
double ElementRecorderRMS::getRecordedValue(int clmnId, int rowOffset, bool reset)
{
	double res = 0;
	if (!initializationDone)
		return res;
	if (clmnId >= runningTotal->Size())
		return res;
	res = (*runningTotal)(clmnId);
	if (count != 0)
	  res = sqrt(res/count);
	if (reset)
	  count = 0;
	return (res*res)/count;
}
