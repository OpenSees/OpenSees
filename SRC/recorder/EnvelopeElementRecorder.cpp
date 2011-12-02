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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-25 23:34:29 $
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

#include <iomanip>
using std::ios;

EnvelopeElementRecorder::EnvelopeElementRecorder(const ID &eleID, Domain &theDom, 
						 const char **argv, int argc,
						 double dT, const char *theFileName)
:numEle(eleID.Size()), responseID(eleID.Size()), theDomain(&theDom),
 fileName(0), deltaT(dT), nextTimeStampToRecord(0.0), 
 db(0), dbColumns(0), numDbColumns(0), 
 data(0), currentData(0), first(true)
{
  theResponses = new Response *[numEle];
  for (int j=0; j<numEle; j++)
    theResponses[j] = 0;

  Information eleInfo(1.0);

  for (int i=0; i<numEle; i++) {
    Element *theEle = theDom.getElement(eleID(i));
    if (theEle == 0) {
      opserr << "WARNING EnvelopeElementRecorder::EnvelopeElementRecorder() -";
      opserr << " no element with tag: " << eleID(i) << " exists in Domain\n";
      theResponses[i] = 0;
    } else {
      theResponses[i] = theEle->setResponse(argv, argc, eleInfo);
      Information &eleInfo = theResponses[i]->getInformation();
      const Vector &eleData = eleInfo.getData();
      numDbColumns += eleData.Size();
    }
  }

  // create the matrix & vector that holds the data
  data = new Matrix(3, numDbColumns);
  currentData = new Vector(numDbColumns);
  if (data == 0 || currentData == 0) {
    opserr << "EnvelopeElementRecorder::EnvelopeElementRecorder() - out of memory\n";
    exit(-1);
  }
  
  // if file is specified, copy name and open the file
  if (theFileName != 0) {

    // create char array to store file name
    int fileNameLength = strlen(theFileName) + 1;
    fileName = new char[fileNameLength];
    if (fileName == 0) {
      opserr << "EnvelopeElementRecorder::EnvelopeElementRecorder - out of memory creating string " <<
	fileNameLength << endln;
      exit(-1);
    }

    // copy file name string
    strcpy(fileName, theFileName);    
  }
}


EnvelopeElementRecorder::EnvelopeElementRecorder(const ID &eleID, Domain &theDom, 
						 const char **argv, int argc,
						 FE_Datastore *database, 
						 const char *tableName, 
						 double dT)
  :numEle(eleID.Size()), responseID(eleID.Size()), theDomain(&theDom),
   fileName(0), deltaT(dT), nextTimeStampToRecord(0.0), 
   db(database), dbColumns(0), numDbColumns(0), 
   data(0), currentData(0), first(true)
{
  Information eleInfo(1.0);
  int i;
  for (i=0; i<numEle; i++) {
    Element *theEle = theDom.getElement(eleID(i));
    if (theEle == 0) {
      opserr << "WARNING EnvelopeElementRecorder::EnvelopeElementRecorder() -";
      opserr << " no element with tag: " << eleID(i) << " exists in Domain\n";
      theResponses[i] = 0;
    } else {
      theResponses[i] = theEle->setResponse(argv, argc, eleInfo);
      if (theResponses[i] != 0) {
      
	// from the response determine no of cols for each
	Information &eleInfo = theResponses[i]->getInformation();
	const Vector &eleData = eleInfo.getData();
	numDbColumns += eleData.Size();
      }
    }
  }

  // create the matrix & vector that holds the data
  data = new Matrix(3, numDbColumns);
  currentData = new Vector(numDbColumns);
  if (data == 0 || currentData == 0) {
    opserr << "EnvelopeElementRecorder::EnvelopeElementRecorder() - out of memory\n";
    exit(-1);
  }

  // create char array to store table name
  int fileNameLength = strlen(tableName) + 1;
  fileName = new char[fileNameLength];
  if (fileName == 0) {
    opserr << "EnvelopeElementRecorder::EnvelopeElementRecorder - out of memory creating string " << fileNameLength << endln;
    exit(-1);
  }

  // copy the strings
  strcpy(fileName, tableName);    

  // now create the columns strings for the database
  // for each element do a getResponse() & print the result
  dbColumns = new char *[numDbColumns];

  char aColumn[256]; // assumes a column name will not be longer than 256 characters
  
  int counter = 0;
  const char *dataToStore = argv[argc-1];

  for (i=0; i<eleID.Size(); i++) {
    int eleTag = eleID(i);
    int numVariables = 0;
    if (theResponses[i]!= 0) {
      Information &eleInfo = theResponses[i]->getInformation();
      const Vector &eleData = eleInfo.getData();
      numVariables = eleData.Size();

      // create the column headings for multiple data for the element

      for (int j=1; j<=numVariables; j++) {
	sprintf(aColumn, "%s_%d_%d",dataToStore,eleTag,j);
	int lenColumn = strlen(aColumn+1);
	char *newColumn = new char[lenColumn];
	sprintf(newColumn, "%s",aColumn);
	dbColumns[counter] = newColumn;
	counter++;
      }
    }
  }
  
  // create the table in the database
  db->createTable(tableName, numDbColumns, dbColumns);
}
  
EnvelopeElementRecorder::~EnvelopeElementRecorder()
{
    if (theResponses != 0) {
      for (int i = 0; i < numEle; i++)
	delete theResponses[i];
      delete [] theResponses;
    }

    if (fileName != 0)
      delete [] fileName;

    if (data != 0)
      delete data;
    
    if (currentData != 0)
      delete currentData;
}


int 
EnvelopeElementRecorder::record(int commitTag, double timeStamp)
{
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
	  // print results to file or stderr depending on whether
	  // a file was opened
	  
	  // from the response determine no of cols for each
	  Information &eleInfo = theResponses[i]->getInformation();
	  const Vector &eleData = eleInfo.getData();
	  for (int j=0; j<eleData.Size(); j++) 
	    (*currentData)(loc++) = eleData(j);
	}
      }
    }

    // check if max or min
    // check if currentData modifies the saved data
    bool writeIt = false;
    if (first == true) {
      for (int i=0; i<numDbColumns; i++) {
	(*data)(0,i) = (*currentData)(i);
	(*data)(1,i) = (*currentData)(i);
	(*data)(2,i) = fabs((*currentData)(i));
	first = false;
	writeIt = true;
      } 
    } else {
      for (int i=0; i<numDbColumns; i++) {
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
	  for (int j=0; j<numDbColumns; j++) 
	    theFile << (*data)(i,j) << " ";
	  theFile << endln;
	}
	theFile.close(); 	  
      } else {
	for (int i=0; i<3; i++) {
	  for (int j=0; j<numDbColumns; j++)
	    (*currentData)(j) = (*data)(i,j);
	  db->insertData(fileName, dbColumns, i, *currentData);
	}
      }
    }    
  }    

  // succesfull completion - return 0
  return result;
}


int 
EnvelopeElementRecorder::playback(int commitTag)
{
  opserr << data;
  return 0;
}

void 
EnvelopeElementRecorder::restart(void)
{
  data->Zero();
  first = true;
}
