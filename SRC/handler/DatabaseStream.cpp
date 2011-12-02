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
                                                                        
// $Revision: 1.1 $
// $Date: 2006-08-03 23:28:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DatabaseStream.cpp,v $
                                                                        
// Written: fmk 
// Date: 10/04
//
// Description: This file contains the class implementation for
// DatabaseStream. 
//
// What: "@(#) DatabaseStream.C, revA"

#include "DatabaseStream.h"
#include <FE_Datastore.h>
#include <Vector.h>
#include <bool.h>

DatabaseStream::DatabaseStream(FE_Datastore *database, const char *tName)
  :OPS_Stream(OPS_STREAM_TAGS_DatabaseStream),
   theDatabase(database), tableName(0), numColumns(0), columns(0), commitTag(0)
{
  //
  // create memory to store the dataDescription and make a copy of it
  //

  if (tName != 0) {
    tableName = new char [strlen(tName)+1];
    strcpy(tableName, tName);
    if (tableName == 0) {
      opserr << "DatabaseStream::DatabaseStream - out of memory creating copy of tableName: " << tName << endln;
    }
  }
}

DatabaseStream::~DatabaseStream()
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
DatabaseStream::write(Vector &data) 
{
  if (descriptionDone == false) {
    descriptionDone = true;
    this->open();
  }
    
  int result = 0;

  if (data.Size() == numColumns)
    if (theDatabase != 0)
      result = theDatabase->insertData(tableName, columns, commitTag, data);
      else {
      opserr << "DatabaseStream::write() - database has not been set\n";
      return -1;
    } 
  else {
    opserr << "DatabaseStream::write() - Vector not of correct size or open() has not been called\n" << numColumns << " " << data.Size() << endln;
    return -1;
  }

  commitTag++;

  return result;
}


int 
DatabaseStream::setDatabase(FE_Datastore &database, const char *tName)
{
  if (tName == 0 || strlen(tName) == 0) {
    opserr << "DatabaseStream::DatabaseStream - no tableName passed\n";
    return -1;
  }

  if (tableName != 0)
    delete [] tableName;
  
  tableName = new char [strlen(tName)+1];
  strcpy(tableName, tName);
  if (tableName == 0) {
    opserr << "DatabaseStream::DatabaseStream - out of memory creating copy of tableName: " << tName << endln;
    return -1;
  }

  theDatabase = &database;
  return 0;
}


int 
DatabaseStream::tag(const char *tagName)
{
  // output the xml for it to the file

  return 0;
}


int 
DatabaseStream::tag(const char *tagName, const char * value)
{
  // output the xml for it to the file

  return 0;
}


int 
DatabaseStream::endTag()
{

  return 0;
}

int 
DatabaseStream::attr(const char *name, int value)
{

  //  outputDatabase << name << " = " << value << endln;
  
  return 0;
}

int 
DatabaseStream::attr(const char *name, double value)
{
  //  outputDatabase << name << " = " << value << endln;

  return 0;
}

int 
DatabaseStream::attr(const char *name, const char *value)
{

  //  outputDatabase << name << " = " << value << endln;

  return 0;
}


// regular stuff
OPS_Stream& 
DatabaseStream::operator<<(char c)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(unsigned char c)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(signed char c)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(const char *s)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(const unsigned char *s)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(const signed char *s)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(const void *p)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(int n)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(unsigned int n)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(long n)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(unsigned long n)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(short n)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(unsigned short n)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(bool b)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(double n)
{

  return *this;
}

OPS_Stream& 
DatabaseStream::operator<<(float n)
{

  return *this;
}

 

int 
DatabaseStream::sendSelf(int commitTag, Channel &theChannel)
{
  opserr << "DatabaseStream::sendSelf() - not yet implemented\n";
  return -1;
}
int 
DatabaseStream::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  opserr << "DatabaseStream::sendSelf() - not yet implemented\n";
  return -1;
}


int 
DatabaseStream::open()
{
  //
  // check the args are valid & that database has been set
  //

  if (theDatabase == 0) {
    opserr << "DatabaseStream::open() - database has not been set\n";
    return -1;
  } 

  if (tableName == 0) {
    opserr << "DatabaseStream::open() - no tableName passed or failed to get memory\n";
    return -1;
  }

  // ask database to create a table
  return theDatabase->createTable(tableName, numColumns, columns);
}

