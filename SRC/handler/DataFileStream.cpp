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
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataFileStream.cpp,v $


#include <DataFileStream.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <ID.h>
#include <Channel.h>
#include <Message.h>

using std::cerr;
using std::ios;
using std::setiosflags;

DataFileStream::DataFileStream(int indent)
  :OPS_Stream(OPS_STREAM_TAGS_DataFileStream), 
   fileOpen(0), fileName(0), indentSize(indent)
{
  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+1];
  for (int i=0; i<indentSize; i++)
    strcpy(indentString, " ");
}


DataFileStream::DataFileStream(const char *file, openMode mode, int indent)
  :OPS_Stream(OPS_STREAM_TAGS_DataFileStream), 
   fileOpen(0), fileName(0), indentSize(indent)
{
  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+1];
  for (int i=0; i<indentSize; i++)
    strcpy(indentString, " ");

  this->setFile(file, mode);
}


DataFileStream::~DataFileStream()
{
  if (fileOpen == 1)
    theFile.close();
  
  if (fileName != 0)
    delete [] fileName;
}

int 
DataFileStream::setFile(const char *name, openMode mode)
{
  if (name == 0) {
    std::cerr << "DataFileStream::setFile() - no name passed\n";
    return -1;
  }

  // first create a copy of the file name
  if (fileName != 0) {
    if (strcmp(fileName, name) != 0)
      delete [] fileName;
    fileName = 0;
  }
  if (fileName == 0) {
    fileName = new char[strlen(name)+1];
    if (fileName == 0) {
      std::cerr << "DataFileStream::setFile() - out of memory copying name: " << name << std::endl;
      return -1;
    }
    
    // copy the strings
    strcpy(fileName, name);
  }

  // if file already open, close it
  if (fileOpen == 1) {
    theFile.close();
    fileOpen = 0;
  }

  if (mode == OVERWRITE) 
    theFile.open(fileName, ios::out);
  else
    theFile.open(fileName, ios::out| ios::app);

  if (theFile.bad()) {
    std::cerr << "WARNING - DataFileStream::setFile()";
    std::cerr << " - could not open file " << fileName << std::endl;

    return -1;
  } else
    fileOpen = 1;

  if (mode == 0)
    theOpenMode = OVERWRITE;
  else
    theOpenMode = APPEND;

  return 0;
}

int 
DataFileStream::open(void)
{
  // check setFile has been called
  if (fileName == 0) {
    std::cerr << "DataFileStream::open(void) - no file name has been set\n";
    return -1;
  }

  // if file already open, return
  if (fileOpen == 1) {
    return 0;
  }

  // open file
  theFile.open(fileName, ios::out| ios::app);
  if (theFile.bad()) {
    std::cerr << "WARNING - DataFileStream::open()";
    std::cerr << " - could not open file " << fileName << std::endl;

    return -1;
  } else
    fileOpen = 1;

  return 0;
}

int 
DataFileStream::close(void)
{
  if (fileOpen != 0)
    theFile.close();
  fileOpen = 0;
  return 0;
}


int 
DataFileStream::setPrecision(int prec)
{
  if (fileOpen != 0)
    theFile << std::setprecision(prec);

  return 0;
}

int 
DataFileStream::setFloatField(floatField field)
{
  if (field == FIXEDD) {
    if (fileOpen != 0)
      theFile << setiosflags(ios::fixed);
  }
  else if (field == SCIENTIFIC) {
    if (fileOpen != 0)
      theFile << setiosflags(ios::scientific);
  }

  return 0;
}


int 
DataFileStream::tag(const char *tagName)
{
  return 0;
}

int 
DataFileStream::tag(const char *tagName, const char *value)
{
  return 0;
}


int 
DataFileStream::endTag()
{
  return 0;
}

int 
DataFileStream::attr(const char *name, int value)
{
  return 0;
}

int 
DataFileStream::attr(const char *name, double value)
{
  return 0;
}

int 
DataFileStream::attr(const char *name, const char *value)
{
  return 0;
}

int 
DataFileStream::write(Vector &data)
{
  (*this) << data;  

  return 0;
}



OPS_Stream& 
DataFileStream::write(const char *s,int n)
{
  if (fileOpen != 0)
    theFile.write(s, n);

  return *this;
}

OPS_Stream& 
DataFileStream::write(const unsigned char*s,int n)
{
  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
DataFileStream::write(const signed char*s,int n)
{
  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
DataFileStream::write(const void *s, int n)
{
  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(char c)
{
  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(unsigned char c)
{
  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(signed char c)
{
  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(const char *s)
{
  // note that we do the flush so that a "/n" before
  // a crash will cause a flush() - similar to what 
  if (fileOpen != 0) {
    theFile << s;
    theFile.flush();
  }

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(const unsigned char *s)
{
  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(const signed char *s)
{
  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(const void *p)
{
/*
  if (fileOpen != 0)
    theFile << p;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(int n)
{
  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(unsigned int n)
{
  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(long n)
{
/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(unsigned long n)
{
/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(short n)
{
/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(unsigned short n)
{
/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(bool b)
{
/*
  if (fileOpen != 0)
    theFile << b;
*/
  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(double n)
{
  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
DataFileStream::operator<<(float n)
{
  if (fileOpen != 0)
    theFile << n;

  return *this;
}


int 
DataFileStream::sendSelf(int commitTag, Channel &theChannel)
{
  static ID idData(2);
  int fileNameLength = 0;
  if (fileName != 0)
    fileNameLength = strlen(fileName);

  idData(0) = fileNameLength;

  if (theOpenMode == OVERWRITE)
    idData(1) = 0;
  else
    idData(1) = 1;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "DataFileStream::sendSelf() - failed to send id data\n";
    return -1;
  }

  if (fileNameLength != 0) {
    Message theMessage(fileName, fileNameLength);
    if (theChannel.sendMsg(0, commitTag, theMessage) < 0) {
      opserr << "DataFileStream::sendSelf() - failed to send message\n";
      return -1;
    }
  }

  return 0;
}

int 
DataFileStream::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID idData(2);

  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "DataFileStream::recvSelf() - failed to recv id data\n";
    return -1;
  }

  int fileNameLength = idData(0);
  if (idData(1) == 0)
    theOpenMode = OVERWRITE;
  else
    theOpenMode = APPEND;

  if (fileNameLength != 0) {
    if (fileName != 0)
      delete [] fileName;
    fileName = new char[fileNameLength+5];
    if (fileName == 0) {
      opserr << "DataFileStream::recvSelf() - out of memory\n";
      return -1;
    }

    Message theMessage(fileName, fileNameLength);
    if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
      opserr << "DataFileStream::recvSelf() - failed to recv message\n";
      return -1;
    }
    sprintf(&fileName[fileNameLength],".%d",commitTag);

    if (this->setFile(fileName, theOpenMode) < 0) {
      opserr << "DataFileStream::DataFileStream() - setFile() failed\n";
      if (fileName != 0) {
	delete [] fileName;
	fileName = 0;
      }
    }
  }
  
  return 0;
}


void
DataFileStream::indent(void)
{
  if (fileOpen != 0)
    for (int i=0; i<numIndent; i++)
      theFile << indentString;
}
