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
                                                                        
// $Revision: 1.10 $
// $Date: 2007-05-03 21:24:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/FileStream.cpp,v $


#include <FileStream.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <ID.h>
#include <Channel.h>
#include <Message.h>
#include <string.h>


using std::cerr;
using std::ios;
using std::setiosflags;

FileStream::FileStream(int indent)
  :OPS_Stream(OPS_STREAM_TAGS_FileStream), 
   fileOpen(0), fileName(0), filePrecision(6), indentSize(indent), sendSelfCount(0)
{
  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+1];
  for (int i=0; i<indentSize; i++)
    strcpy(indentString, " ");
}


FileStream::FileStream(const char *name, openMode mode, int indent)
  :OPS_Stream(OPS_STREAM_TAGS_FileStream), 
   fileOpen(0), fileName(0), filePrecision(6), indentSize(indent), sendSelfCount(0)
{
  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+1];
  for (int i=0; i<indentSize; i++)
    strcpy(indentString, " ");
  
  this->setFile(name, mode);
}


FileStream::~FileStream()
{
  if (fileOpen == 1)
    theFile.close();
  
  if (fileName != 0)
    delete [] fileName;

  if (indentString != 0)
    delete [] indentString;
}

int 
FileStream::setFile(const char *name, openMode mode)
{
  if (name == 0) {
    std::cerr << "FileStream::setFile() - no name passed\n";
    return -1;
  }

  // first create a copy of the file name
  if (fileName != 0) {
    if (strcmp(fileName, name) != 0)
      delete [] fileName;
    fileName = 0;
  }

  if (fileName == 0) {
    fileName = new char[strlen(name)+5];
    if (fileName == 0) {
      std::cerr << "FileStream::setFile() - out of memory copying name: " << name << std::endl;
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

  /* *************************************************************
  if (mode == OVERWRITE) 
    theFile.open(fileName, ios::out);
  else
    theFile.open(fileName, ios::out| ios::app);

  if (theFile.bad()) {
    std::cerr << "WARNING - FileStream::setFile()";
    std::cerr << " - could not open file " << fileName << std::endl;

    return -1;
  } else
    fileOpen = 1;
  *************************************************************** */

  if (mode == 0)
    theOpenMode = OVERWRITE;
  else
    theOpenMode = APPEND;

  return 0;
}

int 
FileStream::open(void)
{
  // check setFile has been called
  if (fileName == 0) {
    std::cerr << "FileStream::open(void) - no file name has been set\n";
    return -1;
  }

  // if file already open, return
  if (fileOpen == 1) {
    return 0;
  }

  if (sendSelfCount != 0) {
    int fileNameLength = strlen(fileName);
    sprintf(&fileName[fileNameLength],".%d",1);
  }


  if (theOpenMode == OVERWRITE) 
    theFile.open(fileName, ios::out);
  else
    theFile.open(fileName, ios::out| ios::app);

  theOpenMode = APPEND;

  if (theFile.bad()) {
    std::cerr << "WARNING - FileStream::setFile()";
    std::cerr << " - could not open file " << fileName << std::endl;
    fileOpen = 0;
    return -1;
  } else
    fileOpen = 1;


  theFile << std::setprecision(filePrecision);

  return 0;
}

int 
FileStream::close(void)
{
  if (fileOpen != 0)
    theFile.close();

  fileOpen = 0;
  return 0;
}


int 
FileStream::setPrecision(int prec)
{
  //  if (fileOpen == 0)
  //    this->open();

  filePrecision = prec;

  if (fileOpen != 0)
    theFile << std::setprecision(prec);

  return 0;
}

int 
FileStream::setFloatField(floatField field)
{
  if (fileOpen == 0)
    this->open();

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
FileStream::tag(const char *tagName)
{
  if (fileOpen == 0)
    this->open();

  // output the xml for it to the file
  this->indent();
  (*this) << tagName << endln;

  numIndent++;

  return 0;
}

int 
FileStream::tag(const char *tagName, const char *value)
{
  if (fileOpen == 0)
    this->open();

  // output the xml for it to the file
  numIndent++;
  this->indent();
  (*this) << tagName << " = " << value << endln;

  numIndent--;

  return 0;
}


int 
FileStream::endTag()
{
  numIndent--;

  return 0;
}

int 
FileStream::attr(const char *name, int value)
{
  if (fileOpen == 0)
    this->open();

  this->indent();
  (*this) << name << " = " << value << endln;
  
  return 0;
}

int 
FileStream::attr(const char *name, double value)
{
  if (fileOpen == 0)
    this->open();

  this->indent();
  (*this) << name << " = " << value << endln;

  return 0;
}

int 
FileStream::attr(const char *name, const char *value)
{
  if (fileOpen == 0)
    this->open();

  this->indent();
  (*this) << name << " = " << value << endln;

  return 0;
}

int 
FileStream::write(Vector &data)
{
  if (fileOpen == 0)
    this->open();

  this->indent();
  (*this) << data << endln;  

  return 0;
}



OPS_Stream& 
FileStream::write(const char *s,int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write(s, n);

  return *this;
}

OPS_Stream& 
FileStream::write(const unsigned char*s,int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
FileStream::write(const signed char*s,int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
FileStream::write(const void *s, int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
FileStream::operator<<(char c)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(unsigned char c)
{  
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(signed char c)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(const char *s)
{
  if (fileOpen == 0)
    this->open();

  // note that we do the flush so that a "/n" before
  // a crash will cause a flush() - similar to what 
  if (fileOpen != 0) {
    theFile << s;
    theFile.flush();
  }

  return *this;
}
OPS_Stream& 
FileStream::operator<<(const unsigned char *s)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(const signed char *s)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(const void *p)
{
/*
  if (fileOpen != 0)
    theFile << p;
*/
  return *this;
}
OPS_Stream& 
FileStream::operator<<(int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(unsigned int n)
{  
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(long n)
{
/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
FileStream::operator<<(unsigned long n)
{
/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
FileStream::operator<<(short n)
{
/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
FileStream::operator<<(unsigned short n)
{
/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
FileStream::operator<<(bool b)
{
/*
  if (fileOpen != 0)
    theFile << b;
*/
  return *this;
}
OPS_Stream& 
FileStream::operator<<(double n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(float n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << n;

  return *this;
}


int 
FileStream::sendSelf(int commitTag, Channel &theChannel)
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
    opserr << "FileStream::sendSelf() - failed to send id data\n";
    return -1;
  }

  if (fileNameLength != 0) {
    Message theMessage(fileName, fileNameLength);
    if (theChannel.sendMsg(0, commitTag, theMessage) < 0) {
      opserr << "FileStream::sendSelf() - failed to send message\n";
      return -1;
    }
  }

  sendSelfCount++;

  return 0;
}

int 
FileStream::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID idData(2);

  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "FileStream::recvSelf() - failed to recv id data\n";
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
      opserr << "FileStream::recvSelf() - out of memory\n";
      return -1;
    }

    Message theMessage(fileName, fileNameLength);
    if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
      opserr << "FileStream::recvSelf() - failed to recv message\n";
      return -1;
    }
    sprintf(&fileName[fileNameLength],".%d",commitTag);

    if (this->setFile(fileName, theOpenMode) < 0) {
      opserr << "FileStream::FileStream() - setFile() failed\n";
      if (fileName != 0) {
	delete [] fileName;
	fileName = 0;
      }
    }
  }
  
  return 0;
}


void
FileStream::indent(void)
{
  if (fileOpen != 0)
    for (int i=0; i<numIndent; i++)
      theFile << indentString;
}

int FileStream::flush() {
  if (theFile.is_open() && theFile.good()) {
    theFile.flush();
  }
  return 0;
}
