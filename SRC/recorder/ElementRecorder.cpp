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
                                                                        
// $Revision: 1.16 $
// $Date: 2004-05-11 00:07:34 $
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

#include <OPS_Globals.h>

#include <iomanip>
using std::ios;

ElementRecorder::ElementRecorder(const ID &eleID, Domain &theDom, 
				 const char **argv, int argc,
				 bool echoTime, double dT, const char *theFileName)
:numEle(eleID.Size()), responseID(eleID.Size()), theDomain(&theDom),
 echoTimeFlag(echoTime), fileName(0), deltaT(dT), nextTimeStampToRecord(0.0), 
 db(0), dbColumns(0), numDbColumns(0), data(0), destroyDatabase(false)
{
  theResponses = new Response *[numEle];
  for (int j=0; j<numEle; j++)
    theResponses[j] = 0;

  Information eleInfo(1.0);
  for (int i=0; i<numEle; i++) {
    Element *theEle = theDom.getElement(eleID(i));
    if (theEle == 0) {
      opserr << "WARNING ElementRecorder::ElementRecorder() -";
      opserr << " no element with tag: " << eleID(i) << " exists in Domain\n";
      theResponses[i] = 0;
    } else {
      theResponses[i] = theEle->setResponse(argv, argc, eleInfo);
    }
  }

  // if file is specified, copy name and open the file
  if (theFileName != 0) {
    // create char array to store file name
    int fileNameLength = strlen(theFileName) + 1;
    fileName = new char[fileNameLength];
    if (fileName == 0) {
      opserr << "ElementRecorder::ElementRecorder - out of memory creating string " <<
	fileNameLength << " long\n";
      exit(-1);
    }

    // copy file name string
    strcpy(fileName, theFileName);    

    // open the file
    theFile.open(fileName, ios::out);
    if (theFile.bad()) {
      opserr << "WARNING - ElementRecorder::ElementRecorder()";
      opserr << " - could not open file " << fileName << endln;
    }    
  }
  
  // no file .. results will be sent to opserr
  else
    fileName = 0;
}


ElementRecorder::ElementRecorder(const ID &eleID, Domain &theDom, 
				 const char **argv, int argc,
				 bool echoTime, 
				 FE_Datastore *database, 
				 const char *tableName, 
				 double dT,
				 bool invokeDatabaseDestructor)
:numEle(eleID.Size()), responseID(eleID.Size()), theDomain(&theDom),
 echoTimeFlag(echoTime), fileName(0), deltaT(dT), nextTimeStampToRecord(0.0), 
 db(database), dbColumns(0), numDbColumns(0), data(0), 
 destroyDatabase(invokeDatabaseDestructor)
{

  numDbColumns = 0;
  if (echoTimeFlag == true) 
    numDbColumns = 1;  // 1 for the pseudo-time

  theResponses = new Response *[numEle];

  for (int jj=0; jj<numEle; jj++)
    theResponses[jj] = 0;

  Information eleInfo(1.0);
  int i;
  for (i=0; i<numEle; i++) {
    Element *theEle = theDom.getElement(eleID(i));
    if (theEle == 0) {
      opserr << "WARNING ElementRecorder::ElementRecorder() -";
      opserr << " no element with tag: " << eleID(i) << " exists in Domain\n";
      theResponses[i] = 0;
    } else {
      theResponses[i] = theEle->setResponse(argv, argc, eleInfo);
      if (theResponses[i] != 0) {
	// from the response type determine no of cols for each
	Information &eleInfo = theResponses[i]->getInformation();
	const Vector &eleData = eleInfo.getData();
	numDbColumns += eleData.Size();
      }
    }
  }

  // create char array to store table name
  if (db != 0 && tableName != 0) {
    int fileNameLength = strlen(tableName) + 1;
    fileName = new char[fileNameLength];
    if (fileName == 0) {
      opserr << "ElementRecorder::ElementRecorder - out of memory creating string " <<
	fileNameLength << " long\n";
      exit(-1);
    }
    
    // copy the strings
    strcpy(fileName, tableName);    
  } else {
    db = 0;
    fileName = 0;
  }

  // now create the columns strings for the database
  // for each element do a getResponse() & print the result
  dbColumns = new char *[numDbColumns];
  
  static char aColumn[1012]; // assumes a column name will not be longer than 256 characters
  
  char *newColumn = new char[5];
  sprintf(newColumn, "%s","time");  
  dbColumns[0] = newColumn;
  
  int lengthString = 0;
  for (i=0; i<argc; i++)
    lengthString += strlen(argv[i])+1;
  char *dataToStore = new char[lengthString];
  lengthString = 0;
  for (int j=0; j<argc; j++) {
    int argLength = strlen(argv[j]);
    strcpy(&dataToStore[lengthString], argv[j]);
    if (j<(argc-1)) {
      lengthString += argLength;
      dataToStore[lengthString] = ' ';
      lengthString ++;
    } else
      lengthString += argLength+1;
  }
  
  int counter = 1;
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

  // create the table in the database
  if (db != 0 && fileName != 0)
    db->createTable(fileName, numDbColumns, dbColumns);
  else {
    opserr << "ElementRecorder::ElementRecorder - database pointer is NULL\n";
  }
  
  // create the vector to hold the data
  data = new Vector(numDbColumns);
}
  
ElementRecorder::~ElementRecorder()
{
    // close the file
    if (!theFile.bad())
	theFile.close();    

    if (theResponses != 0) {
      for (int i = 0; i < numEle; i++)
	delete theResponses[i];
      delete [] theResponses;
    }

    if (fileName != 0)
      delete [] fileName;

  if (dbColumns != 0) {

    for (int i=0; i<numDbColumns; i++) 
      delete [] dbColumns[i];

      delete [] dbColumns;
  }

  if (destroyDatabase == true && db != 0)
    delete db;
}


int 
ElementRecorder::record(int commitTag, double timeStamp)
{
  int result = 0;
  if (deltaT == 0.0 || timeStamp >= nextTimeStampToRecord) {
      
    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;

    if (db == 0) {
      // print out the pseudo time if requested
      if (echoTimeFlag == true) {
	if (!theFile) 
	  opserr << timeStamp << " ";			
	else 
	  theFile << timeStamp << " ";	
      }
      
      // get the responses and write to file if file or opserr specified
      // for each element do a getResponse() & print the result
      for (int i=0; i< numEle; i++) {
	if (theResponses[i] != 0) {
	  
	// ask the element for the reponse
	  int res;
	  if (( res = theResponses[i]->getResponse()) < 0)
	    result = res;
	  else {
	    // print results to file or stderr depending on whether
	    // a file was opened
	    
	    if (theFile.bad())
	      theResponses[i]->Print(opserr);	    
	    else {
	      theResponses[i]->Print(theFile);
	      theFile << "  ";  // added for OSP
	    }
	  }
	} 
      }
      if (theFile.bad()) 
	opserr << endln;
      else {
	theFile << " \n";
	theFile.flush();
      }
    }

    else {  // send the data to the database

      int loc = 0;
      if (echoTimeFlag == true) 
	(*data)(loc++) = timeStamp;

      // for each element do a getResponse() & print the result
      for (int i=0; i< numEle; i++) {
	if (theResponses[i] != 0) {
	  
	// ask the element for the reponse
	  int res;
	  if (( res = theResponses[i]->getResponse()) < 0)
	    result = res;
	  else {
	    Information &eleInfo = theResponses[i]->getInformation();
	    const Vector &eleData = eleInfo.getData();
	    for (int j=0; j<eleData.Size(); j++)
	      (*data)(loc++) = eleData(j);
	  }
	} 
      }
      
      db->insertData(fileName, dbColumns, commitTag, *data);
    }

  }

  // succesfull completion - return 0
  return result;
}


int 
ElementRecorder::playback(int commitTag)
{
    return 0;
}

void 
ElementRecorder::restart(void)
{
  theFile.close();
  theFile.open(fileName, ios::out);
  if (!theFile) {
    opserr << "WARNING - ElementRecorder::restart() - could not open file ";
    opserr << fileName << endln;
  }    
}
