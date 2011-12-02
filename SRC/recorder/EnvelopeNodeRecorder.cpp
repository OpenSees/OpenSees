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
                                                                        


// $Revision: 1.6 $
// $Date: 2004-01-29 23:30:29 $
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

#include <string.h>
#include <iomanip>
using std::ios;

EnvelopeNodeRecorder::EnvelopeNodeRecorder(const ID &dofs, 
			   const ID &nodes, 
			   Domain &theDom,
			   const char *theFileName,
			   const char *dataToStore,
			   double dT,
			   int startFlag)
:theDofs(0), theNodes(0), 
  currentData(0), data(0), 
  theDomain(&theDom), flag(startFlag), deltaT(dT), nextTimeStampToRecord(0.0), 
  db(0), dbColumns(0), first(true)
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

  // verify the nodes exist 
  count = 0;
  int numNode = nodes.Size();
  theNodes = new ID(1, numNode);
  for (i=0; i<numNode; i++) {
    int nodeTag = nodes(i);
    Node *theNode = theDomain->getNode(nodeTag);
    if (theNode == 0) {
      opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - invalid node  " << nodeTag;
      opserr << " does not exist in domain - will be ignored\n";
    } else {
      (*theNodes)[count++] = nodeTag;
    }
  }

  currentData = new Vector(1 + theNodes->Size()*theDofs->Size());
  data = new Matrix(3, theNodes->Size()*theDofs->Size());
  data->Zero();

  // create char array to store file name
  if (theFileName == 0) {
    opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - no file name passed\n";
    exit(-1);
  }

  int fileNameLength = strlen(theFileName) + 1;
  fileName = new char[fileNameLength];
  if (fileName == 0) {
    opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - out of memory creating string " <<
      fileNameLength << "to long\n";
    exit(-1);
  }

  // copy the strings
  strcpy(fileName, theFileName);    

  if ((strcmp(dataToStore, "disp") == 0)) {
    dataFlag = 0;
  } else if ((strcmp(dataToStore, "vel") == 0)) {
    dataFlag = 1;
  } else if ((strcmp(dataToStore, "accel") == 0)) {
    dataFlag = 2;
  } else if ((strcmp(dataToStore, "incrDisp") == 0)) {
    dataFlag = 3;
  } else if ((strcmp(dataToStore, "incrDeltaDisp") == 0)) {
    dataFlag = 4;
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


EnvelopeNodeRecorder::EnvelopeNodeRecorder(const ID &dofs, 
			   const ID &nodes, 
			   Domain &theDom,
			   FE_Datastore *database,
			   const char *dbTable,
			   const char *dataToStore,
			   double dT,
			   int startFlag)
:theDofs(0), theNodes(0),
  currentData(0), data(0), 
  theDomain(&theDom), flag(startFlag), deltaT(dT), nextTimeStampToRecord(0.0), 
  db(database), dbColumns(0),  first(true)
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

  // verify the nodes exist 
  count = 0;
  int numNode = nodes.Size();
  theNodes = new ID(1, numNode);
  for (i=0; i<numNode; i++) {
    int nodeTag = nodes(i);
    Node *theNode = theDomain->getNode(nodeTag);
    if (theNode == 0) {
      opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - invalid node  " << nodeTag;
      opserr << " does not exist in domain - will be ignored\n";
    } else {
      (*theNodes)[count++] = nodeTag;
    }
  }

  // create char array to store table name
  int fileNameLength = strlen(dbTable) + 1;
  fileName = new char[fileNameLength];
  if (fileName == 0) {
    opserr << "EnvelopeNodeRecorder::EnvelopeNodeRecorder - out of memory creating string " <<
      fileNameLength << "too long\n";
    exit(-1);
  }

  // copy the strings
  strcpy(fileName, dbTable);

  if ((strcmp(dataToStore, "disp") == 0)) {
    dataFlag = 0;
  } else if ((strcmp(dataToStore, "vel") == 0)) {
    dataFlag = 1;
  } else if ((strcmp(dataToStore, "accel") == 0)) {
    dataFlag = 2;
  } else if ((strcmp(dataToStore, "incrDisp") == 0)) {
    dataFlag = 3;
  } else if ((strcmp(dataToStore, "incrDeltaDisp") == 0)) {
    dataFlag = 4;
  } else if ((strncmp(dataToStore, "eigen",5) == 0)) {
    int mode = atoi(&(dataToStore[5]));
    opserr << "MODE: " << mode << endln;
    if (mode > 0)
      dataFlag = 10 + mode;
    else
      dataFlag = 6;
  } else {
    dataFlag = 6;
    opserr << "NodeRecorder::NodeRecorder - dataToStore " << dataToStore;
    opserr << "not recognized (disp, vel, accel, incrDisp, incrDeltaDisp)\n";
  }

  // now create the columns strings for the database
  numDbColumns = 1 + nodes.Size()*dofs.Size();
  dbColumns = new char *[numDbColumns];


  static char aColumn[256]; // assumes a column name will not be longer than 256 characters
  
  char *newColumn = new char[5];
  sprintf(newColumn, "%s","time");  
  dbColumns[0] = newColumn;

  int counter = 1;
  for (i=0; i<theNodes->Size(); i++) {
    int nodeTag = (*theNodes)(i);
    for (int j=0; j<theDofs->Size(); j++) {
      int dof = (*theDofs)(j);
      sprintf(aColumn, "Node%d_%s_%d", nodeTag, dataToStore, dof);
      int lenColumn = strlen(aColumn);
      char *newColumn = new char[lenColumn+1];
      strcpy(newColumn, aColumn);
      dbColumns[counter] = newColumn;
      counter++;
    }
  }
  
  // create the table in the database
  db->createTable(dbTable, numDbColumns, dbColumns);
}


EnvelopeNodeRecorder::~EnvelopeNodeRecorder()
{
  if (theDofs != 0)
    delete theDofs;
  
  if (theNodes != 0)
    delete theNodes;

  if (fileName != 0)
    delete [] fileName;
  
  if (dbColumns != 0) {
    for (int i=0; i<numDbColumns; i++)
      delete [] dbColumns[i];
    
    delete [] dbColumns;
  }
}

int 
EnvelopeNodeRecorder::record(int commitTag, double timeStamp)
{
    // now we go get the displacements from the nodes
    int numDOF = theDofs->Size();
    int numNodes = theNodes->Size();
    
    if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {
      
      if (deltaT != 0.0) 
	nextTimeStampToRecord = timeStamp + deltaT;

      for (int i=0; i<numNodes; i++) {
	int cnt = i*numDOF;
	Node *theNode = theDomain->getNode((*theNodes)(i));
	if (theNode != 0) {
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
	for (int i=0; i<numNodes*numDOF; i++) {
	  (*data)(0,i) = (*currentData)(i);
	  (*data)(1,i) = (*currentData)(i);
	  (*data)(2,i) = fabs((*currentData)(i));
	  first = false;
	  writeIt = true;
	} 
      } else {
	for (int i=0; i<numNodes*numDOF; i++) {
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

      // if any value has changed .. write it out
      if (writeIt == true) {
	if (db == 0) {      
	  
	  // open the file
	  if (fileName != 0)
	  theFile.open(fileName, ios::out);
	  if (theFile.bad()) {
	    opserr << "WARNING - EnvelopeNodeRecorder::EnvelopeNodeRecorder()";
	    opserr << " - could not open file " << fileName << endln;
	  }    
	  
	  for (int i=0; i<3; i++) {
	    for (int j=0; j<numNodes*numDOF; j++)
	      theFile << (*data)(i,j) << " ";
      	    theFile << endln;
	  }
	  theFile.flush();
	  theFile.close(); 	  
	} else {
	  for (int i=0; i<3; i++) {
	    for (int j=1; j<=numNodes*numDOF; j++)
	      (*currentData)(j) = (*data)(i,j);
	    db->insertData(fileName, dbColumns, i, (*currentData));
	  }
	}
      }    
    }

    return 0;
}

int 
EnvelopeNodeRecorder::playback(int commitTag)
{
  opserr << data;
  
  // does nothing
  return 0;
}


void
EnvelopeNodeRecorder::restart(void)
{
  data->Zero();
  first = true;
}





