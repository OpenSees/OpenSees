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
                                                                        


// $Revision: 1.9 $
// $Date: 2004-11-25 00:53:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/EnvelopeNodeRecorder.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the class definition for EnvelopeNodeRecorder.
// A EnvelopeNodeRecorder is used to record the envelop of specified dof responses 
// at a collection of nodes over an analysis. (between commitTag of 0 and
// last commitTag).
//
// What: "@(#) EnvelopeNodeRecorder.C, revA"

#include <math.h>

#include <EnvelopeNodeRecorder.h>
#include <Domain.h>
#include <Node.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <FE_Datastore.h>
#include <FEM_ObjectBroker.h>
#include <DataOutputHandler.h>

#include <string.h>

EnvelopeNodeRecorder::EnvelopeNodeRecorder()
:Recorder(RECORDER_TAGS_EnvelopeNodeRecorder),
 theDofs(0), theNodalTags(0), theNodes(0),
 currentData(0), data(0), 
 theDomain(0), theHandler(0),
 deltaT(0.0), nextTimeStampToRecord(0.0), 
 first(true), initializationDone(false), numValidNodes(0)
{

}

EnvelopeNodeRecorder::EnvelopeNodeRecorder(const ID &dofs, 
					   const ID &nodes, 
					   const char *dataToStore,
					   Domain &theDom,
					   DataOutputHandler &theOutputHandler,
					   double dT)
:Recorder(RECORDER_TAGS_EnvelopeNodeRecorder),
 theDofs(0), theNodalTags(0), theNodes(0),
 currentData(0), data(0), 
 theDomain(&theDom), theHandler(&theOutputHandler),
 deltaT(dT), nextTimeStampToRecord(0.0), 
 first(true), initializationDone(false), numValidNodes(0)
{
  // verify dof are valid 
  int numDOF = dofs.Size();
  theDofs = new ID(0, numDOF);

  int count = 0;
  int i;
  for (i=0; i<numDOF; i++) {
    int dof = dofs(i);
    if (dof >= 0) {
      (*theDofs)[count] = dof;
      count++;
    } else {
      opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - invalid dof  " << dof;
      opserr << " will be ignored\n";
    }
  }


  // 
  // create memory to hold nodal ID's (neeed parallel)
  //

  int numNode = nodes.Size();
  if (numNode != 0) {
    theNodalTags = new ID(nodes);
    
    if (theNodalTags == 0) {
      opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - out of memory\n";
    }
  }

  //
  // set the data flag used as a switch to get the response in a record
  //

  if (dataToStore == 0 || (strcmp(dataToStore, "disp") == 0)) {
    dataFlag = 0;
  } else if ((strcmp(dataToStore, "vel") == 0)) {
    dataFlag = 1;
  } else if ((strcmp(dataToStore, "accel") == 0)) {
    dataFlag = 2;
  } else if ((strcmp(dataToStore, "incrDisp") == 0)) {
    dataFlag = 3;
  } else if ((strcmp(dataToStore, "incrDeltaDisp") == 0)) {
    dataFlag = 4;
  } else if ((strcmp(dataToStore, "unbalance") == 0)) {
    dataFlag = 5;
  } else if ((strncmp(dataToStore, "eigen",5) == 0)) {
    int mode = atoi(&(dataToStore[5]));
    if (mode > 0)
      dataFlag = 10 + mode;
    else
      dataFlag = 6;
  } else {
    dataFlag = 6;
    opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - dataToStore " << dataToStore;
    opserr << "not recognized (disp, vel, accel, incrDisp, incrDeltaDisp)\n";
  }
}


EnvelopeNodeRecorder::~EnvelopeNodeRecorder()
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

  if (theDofs != 0)
    delete theDofs;

  if (theNodalTags != 0)
    delete theNodalTags;
  
  if (theHandler != 0)
    delete theHandler;

  if (currentData != 0)
    delete currentData;

  if (data != 0)
    delete data;

  if (theNodes != 0)
    delete [] theNodes;
}

int 
EnvelopeNodeRecorder::record(int commitTag, double timeStamp)
{
  if (theDomain == 0 || theNodalTags == 0 || theDofs == 0) {
    return 0;
  }

  if (theHandler == 0) {
    opserr << "EnvelopeNodeRecorder::record() - no DataOutputHandler has been set\n";
    return -1;
  }


  if (initializationDone != true) 
    if (this->initialize() != 0) {
      opserr << "EnvelopeNodeRecorder::record() - failed in initialize()\n";
      return -1;
    }

  int numDOF = theDofs->Size();
  int sizeData = currentData->Size();
  
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {
    
    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;
    
    for (int i=0; i<numValidNodes; i++) {
      int cnt = i*numDOF;
      Node *theNode = theNodes[i];
      if (dataFlag == 0) {
	const Vector &response = theNode->getTrialDisp();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);
	  if (response.Size() > dof) {
	    (*currentData)(cnt) = response(dof);
	  }else 
	    (*currentData)(cnt) = 0.0;
	  
	  cnt++;
	}
      } else if (dataFlag == 1) {
	const Vector &response = theNode->getTrialVel();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);
	  if (response.Size() > dof) {
	    (*currentData)(cnt) = response(dof);
	  } else 
	    (*currentData)(cnt) = 0.0;
	  
	  cnt++;
	}
      } else if (dataFlag == 2) {
	const Vector &response = theNode->getTrialAccel();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);
	  if (response.Size() > dof) {
	    (*currentData)(cnt) = response(dof);
	  } else 
	    (*currentData)(cnt) = 0.0;
	  
	  cnt++;
	}
      } else if (dataFlag == 3) {
	const Vector &response = theNode->getIncrDisp();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);
	  if (response.Size() > dof) {
	    (*currentData)(cnt) = response(dof);
	  } else 
	    (*currentData)(cnt) = 0.0;
	  
	  cnt++;
	}
      } else if (dataFlag == 4) {
	const Vector &response = theNode->getIncrDeltaDisp();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);
	  if (response.Size() > dof) {
	    (*currentData)(cnt) = response(dof);
	  } else 
	    (*currentData)(cnt) = 0.0;
	  
	  cnt++;
	}
      } else if (dataFlag == 5) {
	const Vector &theResponse = theNode->getUnbalancedLoad();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);
	  if (theResponse.Size() > dof) {
	    (*currentData)(cnt) = theResponse(dof);
	  } else 
	    (*currentData)(cnt) = 0.0;
	  
	  cnt++;
	}
      } else if (dataFlag > 10) {
	int mode = dataFlag - 10;
	int column = mode - 1;
	const Matrix &theEigenvectors = theNode->getEigenvectors();
	if (theEigenvectors.noCols() > column) {
	  int noRows = theEigenvectors.noRows();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (noRows > dof) {
	      (*currentData)(cnt) = theEigenvectors(dof,column);
	    } else 
	      (*currentData)(cnt) = 0.0;
	    cnt++;		
	  }
	} else {
	  for (int j=0; j<numDOF; j++) {
	    (*currentData)(cnt) = 0.0;
	    cnt++;		
	  }
	}
      }
    }
  }
    
  // check if currentData modifies the saved data
  bool writeIt = false;
  if (first == true) {
    for (int i=0; i<sizeData; i++) {
      (*data)(0,i) = (*currentData)(i);
      (*data)(1,i) = (*currentData)(i);
      (*data)(2,i) = fabs((*currentData)(i));
      first = false;
      writeIt = true;
    } 
  } else {
    for (int i=0; i<sizeData; i++) {
      double value = (*currentData)(i);
      if ((*data)(0,i) > value) {
	(*data)(0,i) = value;
	double absValue = fabs(value);
	if ((*data)(2,i) < absValue) 
	  (*data)(2,i) = absValue;
	writeIt = true;
      } else if ((*data)(1,i) < value) {
	(*data)(1,i) = value;
	double absValue = fabs(value);
	if ((*data)(2,i) < absValue) 
	  (*data)(2,i) = absValue;
	writeIt = true;
      }
    }
  }

  return 0;
}


int
EnvelopeNodeRecorder::restart(void)
{
  data->Zero();
  first = true;
  return 0;
}


int 
EnvelopeNodeRecorder::setDomain(Domain &theDom)
{
  theDomain = &theDom;
  return 0;
}


int
EnvelopeNodeRecorder::sendSelf(int commitTag, Channel &theChannel)
{
  if (theChannel.isDatastore() == 0) {
    opserr << "EnvelopeNodeRecorder::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  static ID idData(4); 
  idData.Zero();
  if (theDofs != 0)
    idData(0) = theDofs->Size();
  if (theNodalTags != 0)
    idData(1) = theNodalTags->Size();
  if (theHandler != 0) {
    idData(2) = theHandler->getClassTag();
  }

  idData(3) = dataFlag;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }

  if (theDofs != 0) 
    if (theChannel.sendID(0, commitTag, *theDofs) < 0) {
      opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send dof id's\n";
      return -1;
    }

  if (theNodalTags != 0)
    if (theChannel.sendID(0, commitTag, *theNodalTags) < 0) {
      opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send nodal tags\n";
      return -1;
    }

  static Vector data(2);
  data(0) = deltaT;
  data(1) = nextTimeStampToRecord;
  if (theChannel.sendVector(0, commitTag, data) < 0) {
    opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send data\n";
    return -1;
  }

  if (theHandler->sendSelf(commitTag, theChannel) < 0) {
    opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  return 0;
}

int 
EnvelopeNodeRecorder::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  if (theChannel.isDatastore() == 0) {
    opserr << "EnvelopeNodeRecorder::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  static ID idData(4); 
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "EnvelopeNodeRecorder::recvSelf() - failed to send idData\n";
    return -1;
  }

  int numDOFs = idData(0);
  int numNodes = idData(1);

  dataFlag = idData(3);

  //
  // get the DOF ID data
  //

  if (theDofs == 0 || theDofs->Size() != numDOFs) {
    if (theDofs != 0)
      delete theDofs;

    if (numDOFs != 0) {
      theDofs = new ID(numDOFs);
      if (theDofs == 0 || theDofs->Size() != numDOFs) {
	opserr << "EnvelopeNodeRecorder::recvSelf() - out of memory\n";
	return -1;
      }	
    }
  }
  if (theDofs != 0)
    if (theChannel.recvID(0, commitTag, *theDofs) < 0) {
      opserr << "EnvelopeNodeRecorder::recvSelf() - failed to recv dof data\n";
      return -1;
    } 

  //
  // get the NODAL tag data
  //

  if (theNodalTags == 0 || theNodalTags->Size() != numNodes) {
    if (theNodalTags != 0)
      delete theNodalTags;

    if (numNodes != 0) {
      theNodalTags = new ID(numNodes);
      if (theNodalTags == 0 || theNodalTags->Size() != numNodes) {
	opserr << "EnvelopeNodeRecorder::recvSelf() - out of memory\n";
	return -1;
      }	
    }
  }
  if (theNodalTags != 0)
    if (theChannel.recvID(0, commitTag, *theNodalTags) < 0) {
      opserr << "EnvelopeNodeRecorder::recvSelf() - failed to recv dof data\n";
      return -1;
    } 


  static Vector data(2);
  data(0) = deltaT;
  data(1) = nextTimeStampToRecord;
  if (theChannel.recvVector(0, commitTag, data) < 0) {
    opserr << "EnvelopeNodeRecorder::sendSelf() - failed to receive data\n";
    return -1;
  }

  if (theHandler != 0)
    delete theHandler;

  theHandler = theBroker.getPtrNewDataOutputHandler(idData(2));
  if (theHandler == 0) {
    opserr << "EnvelopeNodeRecorder::sendSelf() - failed to get a data output handler\n";
    return -1;
  }

  if (theHandler->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  return 0;
}


int
EnvelopeNodeRecorder::initialize(void)
{
  if (theDofs == 0 || theNodalTags == 0 || theDomain == 0) {
    opserr << "EnvelopeNodeRecorder::initialize() - either nodes, dofs or domain has not been set\n";
    return -1;
  }
  
  //
  // create & set nodal array pointer
  //

  if (theNodes != 0) 
    delete [] theNodes;
  
  numValidNodes = 0;
  int i;
  int numNode = theNodalTags->Size();
  for (i=0; i<numNode; i++) {
    int nodeTag = (*theNodalTags)(i);
    Node *theNode = theDomain->getNode(nodeTag);
    if (theNode != 0) {
      numValidNodes++;
    }
  }

  theNodes = new Node *[numValidNodes];
  if (theNodes == 0) {
    opserr << "EnvelopeNodeRecorder::domainChanged - out of memory\n";
    return -1;
  }

  int count = 0;
  for (i=0; i<numNode; i++) {
    int nodeTag = (*theNodalTags)(i);
    Node *theNode = theDomain->getNode(nodeTag);
    if (theNode != 0) {
      theNodes[count] = theNode;
      count++;
    }
  }

  //
  // resize the response vector
  //
  int numValidResponse = numValidNodes*theDofs->Size();

  currentData = new Vector(numValidResponse);
  data = new Matrix(3, numValidResponse);
  data->Zero();

  //
  // need to create the data description, i.e. what each column of data is
  //
  
  char dataToStore[16];
  if (dataFlag == 0)
    strcpy(dataToStore, "disp");
  else if (dataFlag == 1)
    strcpy(dataToStore, "vel");
  else if (dataFlag == 2)
    strcpy(dataToStore, "accel");
  else if (dataFlag == 3)
    strcpy(dataToStore, "deltaDisp");
  else if (dataFlag == 4)
    strcpy(dataToStore, "incrDeltaDisp");
  else if (dataFlag == 5)
    sprintf(dataToStore,"eigen_%d",dataToStore-10);
  else
    strcpy(dataToStore,"invalid");
  
  int numDbColumns = numValidResponse;
  char **dbColumns = new char *[numDbColumns];

  static char aColumn[128]; // assumes a column name will not be longer than 256 characters
  
  int counter = 0;
  for (i=0; i<numValidNodes; i++) {
    int nodeTag = theNodes[i]->getTag();
    for (int j=0; j<theDofs->Size(); j++) {
      int dof = (*theDofs)(j);
      sprintf(aColumn, "Node%d_%s_%d", nodeTag, dataToStore, dof+1);
      int lenColumn = strlen(aColumn);
      char *newColumn = new char[lenColumn+1];
      strcpy(newColumn, aColumn);
      dbColumns[counter] = newColumn;
      counter++;
    }
  }
  
  //
  // call open in the handler with the data description
  //

  if (theHandler != 0)
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
