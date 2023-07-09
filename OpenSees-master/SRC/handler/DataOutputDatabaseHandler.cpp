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
                                                                        
// $Revision: 1.2 $
// $Date: 2004-11-24 22:40:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataOutputDatabaseHandler.cpp,v $
                                                                        
// Written: fmk 
// Date: 10/04
//
// Description: This file contains the class implementation for
// DataOutputDatabaseHandler. 
//
// What: "@(#) DataOutputDatabaseHandler.C, revA"

#include "DataOutputDatabaseHandler.h"
#include <FE_Datastore.h>
#include <Vector.h>

DataOutputDatabaseHandler::DataOutputDatabaseHandler(FE_Datastore *database, const char *tName)
  :DataOutputHandler(DATAHANDLER_TAGS_DataOutputDatabaseHandler),
   theDatabase(database), tableName(0), numColumns(0), columns(0), commitTag(0)
{
  //
  // create memory to store the dataDescription and make a copy of it
  //

  if (tName != 0) {
    tableName = new char [strlen(tName)+1];
    strcpy(tableName, tName);
    if (tableName == 0) {
      opserr << "DataOutputDatabaseHandler::DataOutputDatabaseHandler - out of memory creating copy of tableName: " << tName << endln;
    }
  }
}

DataOutputDatabaseHandler::~DataOutputDatabaseHandler()
{
  if (tableName != 0)
    delete [] tableName;

  if (columns != 0) {
    for (int j=0; j<numColumns; j++)
      delete [] columns[j];
    delete [] columns;
  }
}

int 
DataOutputDatabaseHandler::open(char **dataDescription, int numData)
{
  //
  // check the args are valid & that database has been set
  //

  if (theDatabase == 0) {
    opserr << "DataOutputStreamHandler::open() - database has not been set\n";
    return -1;
  } 

  if (tableName == 0) {
    opserr << "DataOutputDatabaseHandler::open() - no tableName passed or failed to get memory\n";
    return -1;
  }

  if (dataDescription == 0) {
    opserr << "DataOutputDatabaseHandler::open() - no column description data passed\n";
    return -1;
  }

  if (numData < 0) {
    opserr << "DataOutputDatabaseHandler::open() - numColumns (" << numData << ") < 0\n";
    return -1;
  } else
    numColumns = numData;

  //
  // remove old if open has already been called
  //

  if (columns != 0) {
    for (int j=0; j<numColumns; j++)
      delete [] columns[j];
    delete [] columns;
    columns = 0;
  }

  //
  // create memory to store the dataDescription and make a copy of it
  //

  columns = new char *[numColumns];
  if (columns == 0) {
    opserr << "DataOutputDatabaseHandler::open() - out of memory creating array for columns of size (" << numData << ") < 0\n";
    numColumns = 0;
    return -1;
  } 

  // make copy
  for (int i=0; i<numColumns; i++) {
    columns[i] = new char[strlen(dataDescription[i])+1];
    if (columns[i] == 0) {
      opserr << "DataOutputDatabaseHandler::open() - out of memory creating copy of string " << dataDescription[i] << endln;
      for (int j=0; j<i; j++)
	delete [] columns[j];
      delete [] columns;
      numColumns = 0;
      return -1;
    }
    strcpy(columns[i], dataDescription[i]);
  }

  // ask database to create a table
  return theDatabase->createTable(tableName, numColumns, columns);
}

int 
DataOutputDatabaseHandler::write(Vector &data) 
{
  int result = 0;

  if (data.Size() == numColumns)
    if (theDatabase != 0)
      result = theDatabase->insertData(tableName, columns, commitTag, data);
      else {
      opserr << "DataOutputStreamHandler::write() - database has not been set\n";
      return -1;
    } 
  else {
    opserr << "DataOutputStreamHandler::write() - Vector not of correct size or open() has not been called\n" << numColumns << " " << data.Size() << endln;
    return -1;
  }

  commitTag++;

  return result;
}

int 
DataOutputDatabaseHandler::setDatabase(FE_Datastore &database, const char *tName)
{
  if (tName == 0 || strlen(tName) == 0) {
    opserr << "DataOutputDatabaseHandler::DataOutputDatabaseHandler - no tableName passed\n";
    return -1;
  }

  if (tableName != 0)
    delete [] tableName;
  
  tableName = new char [strlen(tName)+1];
  strcpy(tableName, tName);
  if (tableName == 0) {
    opserr << "DataOutputDatabaseHandler::DataOutputDatabaseHandler - out of memory creating copy of tableName: " << tName << endln;
    return -1;
  }

  theDatabase = &database;
  return 0;
}


int 
DataOutputDatabaseHandler::sendSelf(int commitTag, Channel &theChannel)
{
  opserr << "DataOutputDatabaseHandler::sendSelf() - not yet implemented\n";
  return -1;
}
int 
DataOutputDatabaseHandler::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  opserr << "DataOutputDatabaseHandler::sendSelf() - not yet implemented\n";
  return -1;
}
