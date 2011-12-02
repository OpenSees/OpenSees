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
                                                                        
// $Revision: 1.5 $
// $Date: 2006-01-10 18:16:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/NEESData.cpp,v $
                                                                        
// Written: fmk 
// Created: 03/04
//
// Description: This file contains the class definition for NeesData.
// NeesData is a concretes subclass of FE_Datastore. NEESData datastore
// is used to store data from the recorders in a file format that can
// be stored in the NEES data repository.

// NEESData objects currently only implement the createTable() and 
// insertDataIntoTable() methods of the abstract base FE_Datastore
// class. At this time the other methods should never be called.
//
// What: "@(#) NEESData.C, revA"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
//#include <bool.h>
#include <iostream>
#include <iomanip>
using std::ios;
using std::setiosflags;

#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <NEESData.h>

NEESData::NEESData(const char *dataBaseName,
		   Domain &theDomain, 
		   FEM_ObjectBroker &theObjBroker) 
  :FE_Datastore(theDomain, theObjBroker), tables(0), numTables(0)
{
  dataBase = new char [strlen(dataBaseName)+1];
  strcpy(dataBase, dataBaseName);
}

NEESData::~NEESData() 
{
  // clean up all memory allocated
  for (int i=0; i<numTables; i++) {
    NEES_table* nextTable = tables->next;
    int numColumns = tables->numColumns;
    char **columns = tables->columns;
    for (int j=0; j<numColumns; j++)
      delete [] columns[j];
    delete [] columns;
    delete tables;
    tables = nextTable;
  }
}

int
NEESData::getDbTag(void)
{
  opserr << "NEESData::getDbTag(void) - should not be called\n";
  return -1;
}

int 
NEESData::commitState(int commitTag)
{
  opserr << "NEESData::commitState(int commitTag) - should not be called\n";
  return -1;
}

int 
NEESData::sendMsg(int dataTag, int commitTag, 
		       const Message &, 
		       ChannelAddress *theAddress)
{
  opserr << "NEESData::sendMsg(int dataTag, int commitTag,  - should not be called\n";
  return -1;
}		       

int 
NEESData::recvMsg(int dataTag, int commitTag, 
		       Message &, 
		       ChannelAddress *theAddress)
{
  opserr << "NEESData::recvMsg(int dataTag, int commitTag,  - should not be called\n";
  return -1;
}		       


int 
NEESData::sendMatrix(int dataTag, int commitTag, 
		      const Matrix &theMatrix, 
		      ChannelAddress *theAddress)
{
  opserr << "NEESData::sendMatrix(int dataTag, int commitTag,  - should not be called\n";
  return -1;
}		       




int 
NEESData::recvMatrix(int dataTag, int commitTag, 
		      Matrix &theMatrix, 
		      ChannelAddress *theAddress)    
{
  return -1;
}		       



int 
NEESData::sendVector(int dataTag, int commitTag, 
		      const Vector &theVector, 
		      ChannelAddress *theAddress)
{
  opserr << "NEESData::sendVector(int dataTag, int commitTag,  - should not be called\n";
  return -1;
}		       

int 
NEESData::recvVector(int dataTag, int commitTag, 
		      Vector &theVector, 
		      ChannelAddress *theAddress)    
{
  opserr << "NEESData::recvVector(int dataTag, int commitTag,  - should not be called\n";
  return -1;
}		       



int 
NEESData::sendID(int dataTag, int commitTag, 
		      const ID &theID, 
		      ChannelAddress *theAddress)
{
  opserr << "NEESData::sendID(int dataTag, int commitTag,  - should not be called\n";
  return -1;
}

int 
NEESData::recvID(int dataTag, int commitTag, 
		      ID &theID, 
		      ChannelAddress *theAddress)    
{
  opserr << "NEESData::recvID(int dataTag, int commitTag,  - should not be called\n";
  return -1;
}		       


#include <fstream>
using std::ofstream;
using std::cerr;


int 
NEESData::createTable(const char *tableName, int numColumns, char *columns[])
{
  // create a copy of the data file name
  int res = 0;
  int tableNameLength = strlen(tableName);
  char *copyTableName = new char[tableNameLength + 5];
  char *fileName = new char[tableNameLength + 5];

  bool hasOutExtension = false;

  if (fileName == 0 || copyTableName == 0) {
    opserr << "NEESData::insertData - out of memory creating copy of string " << fileName << endln;
    return -1;
  }

  strcpy(copyTableName, tableName);
  if (tableNameLength > 4 && strcmp(".out", &tableName[tableNameLength-4]) == 0) {
    copyTableName[tableNameLength-4] = '\0';
    hasOutExtension = true;
  }

  // first ensure createTable has not already been called with this table name
  NEES_table *t = tables;
  for (int i=0; i<numTables; i++) {
    if (strcmp(t->name, tableName) == 0) {
      opserr << "WARNING: NEESData::createTable - table already exists: " << copyTableName << endln;
      return -1;
    } 
    t = t->next;
  }

  strcpy(fileName, copyTableName);
  strcat(fileName,".out");

  if (numColumns <= 0) {
    opserr << "WARNING: NEESData::createTable - number of data columns < 0 for table name: " << copyTableName << endln;
    delete [] fileName;
    delete [] copyTableName;
    return -1;
  }

  // 
  // test opening the data file
  // 

  strcpy(fileName, copyTableName);    
  strcat(fileName,".out");

  ofstream dataFile;
  dataFile.open(fileName, ios::out | ios::trunc); 

  if (!dataFile.is_open()) {
    opserr << "NEESData::createTable - failed to open file: " << fileName << endln;
    delete [] fileName;
    delete [] copyTableName;
    return -1;
  }

  // 
  // now open the xml file & write the data to it
  // 

  strcpy(fileName, copyTableName);    
  strcat(fileName,".xml");

  ofstream xmlFile;
  xmlFile.open(fileName, ios::out | ios::trunc); 
  
  if (xmlFile.is_open()) {
    // write the xml data
    xmlFile << "<?xml version=\"1.0\"?>\n";
    xmlFile << "<NumericalFileDataDescription>\n";
    xmlFile << "\t<DataFile>\n";
    xmlFile << "\t\t<DataFileName> " << copyTableName << ".out </DataFileName>\n";
    xmlFile << "\t\t<NumberDataColumns> " << numColumns << "</NumberDataColumns>\n";
    xmlFile << "\t</DataFile>\n";
    for (int i=0; i<numColumns; i++) {
      xmlFile << "\t<DataColumnDescription>\n";
      xmlFile << "\t\t<ColumnLocation> " << i+1 << "</ColumnLocation>\n";      
      xmlFile << "\t\t<Description> " << columns[i] << "</Description>\n";
      xmlFile << "\t</DataColumnDescription>\n";
    }
    xmlFile << "</NumericalFileDataDescription>\n";
    xmlFile.close();
    
  } else {
    opserr << "NEESData::createData - failed to open file: " << copyTableName << endln;
    delete [] fileName;
    delete [] copyTableName;
    res = -1;
  }

  // no longer need copyTableName
  delete [] copyTableName;

  // if we can get here we can open the data and have written to the xml file
  // now add a new struct to our list of tables
  NEES_table *nextTable = new NEES_table;
  if (nextTable == 0) {
    opserr << "NEESData::createData - out of memory creating Table structure for table: " << copyTableName << endln;
    res = -1;
  }    
  strcpy(fileName, tableName);

  nextTable->name = fileName;
  nextTable->numColumns = numColumns;
  nextTable->columns = new char *[numColumns];
  for (int ii=0; ii<numColumns; ii++)
    nextTable->columns[ii] = 0;
  nextTable->hasOutExtension = hasOutExtension;

  if (nextTable->columns == 0) {
    opserr << "NEESData::createData - out of memory creating Table structure for table: " << copyTableName << endln;
    delete nextTable;
    res = -1;
  }

  for (int k=0; k<numColumns; k++) {

    nextTable->columns[k] = new char [strlen(columns[k])+1];

    if (nextTable->columns[k] == 0) {
      opserr << "NEESData::createData - out of memory creating Table structure for table: " << copyTableName << endln;
      for (int l=0; l<k-1; l++)
	delete nextTable->columns[k];
      delete nextTable;
      res = -1;    
    }

    strcpy(nextTable->columns[k], columns[k]);
  }    

  nextTable->next = tables;
  tables = nextTable;
  numTables++;

  return 0;
}

int 
NEESData::insertData(const char *tableName, char *columns[], 
		     int commitTag, const Vector &data)
{
  // search the tables for valid table
  NEES_table *t = tables;
  for (int i=0; i<numTables; i++) {
    if (strcmp(t->name, tableName) == 0) 
      i = numTables;
    else 
      t = t->next;
  }

  if ( t == 0) {
    opserr << "NEESData::insertData - table: " << tableName << " has not been created\n";
    return -1;
  }      

  if (t->numColumns != data.Size()) {
    opserr << "NEESData::insertData - incorrect number of columns for table: " << tableName << "\n";
    return -2;
  }

  char *fileName = t->name;
  if (t->hasOutExtension == false)
    strcat(fileName,".out");

  ofstream table;
  table.open(fileName, ios::app); 

  table << setiosflags(ios::scientific);
  table << std::setprecision(16);

  if (table.is_open()) {
    // write the data
    for (int i=0; i<data.Size(); i++) {
      table << data(i) << "\t";
    }
    
    table << "\n";
    table.close();

  } else {
    opserr << "NEESData::insertData - failed to open file: " << fileName << endln;
    return -1;
  }

  strcpy(fileName, tableName);
  return 0;
}


int 
NEESData::getData(const char *tableName, char *columns[], int commitTag, Vector &data)
{
  return 0;
}

