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
// $Date: 2008-01-16 00:19:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/BinaryFileStream.cpp,v $


#include <BinaryFileStream.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <ID.h>
#include <Channel.h>
#include <Message.h>

using std::cerr;
using std::ios;
using std::setiosflags;
using std::ifstream;
using std::string;
using std::getline;

BinaryFileStream::BinaryFileStream()
  :OPS_Stream(OPS_STREAM_TAGS_BinaryFileStream), 
   fileOpen(0), fileName(0), sendSelfCount(0)
{

}


BinaryFileStream::BinaryFileStream(const char *file, openMode mode)
  :OPS_Stream(OPS_STREAM_TAGS_BinaryFileStream), 
   fileOpen(0), fileName(0), sendSelfCount(0)
{
  this->setFile(file, mode);
}


BinaryFileStream::~BinaryFileStream()
{
  if (fileOpen == 1)
    theFile.close();

  if (sendSelfCount != 0) {

    int fileNameLength = strlen(fileName);
    sprintf(&fileName[fileNameLength-2],"");
    
    theFile.open(fileName, ios::out);

    ifstream **theFiles = new ifstream *[sendSelfCount+1];
    
    // open up the files
    for (int i=0; i<=sendSelfCount; i++) {
      theFiles[i] = new ifstream;
      sprintf(&fileName[fileNameLength-2],".%d",i+1);
      theFiles[i]->open(fileName, ios::in);
    }

    // go through each file, reading a line & sending to the output file
    bool done = false;
    string s;

    while (done == false) {
      for (int i=0; i<=sendSelfCount; i++) {
	getline(*(theFiles[i]), s);	
	theFile << s;
	
	if (theFiles[i]->eof()) {
	  done = true;
	  theFiles[i]->close();
	  delete theFiles[i];
	}
      }
      theFile << "\n";
    }

    delete [] theFiles;
    theFile.close();
  }

  if (fileName != 0)
    delete [] fileName;
}

int 
BinaryFileStream::setFile(const char *name, openMode mode)
{
  if (name == 0) {
    std::cerr << "BinaryFileStream::setFile() - no name passed\n";
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
      std::cerr << "BinaryFileStream::setFile() - out of memory copying name: " << name << std::endl;
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

  if (mode == 0)
    theOpenMode = OVERWRITE;
  else
    theOpenMode = APPEND;

  return 0;
}

int 
BinaryFileStream::open(void)
{
  // check setFile has been called
  if (fileName == 0) {
    std::cerr << "BinaryFileStream::open(void) - no file name has been set\n";
    return -1;
  }

  // if file already open, return
  if (fileOpen == 1) {
    return 0;
  }

  if (sendSelfCount != 0) {
    strcat(fileName, ".1");
  }

  if (theOpenMode == OVERWRITE) 
    theFile.open(fileName, ios::out | ios::binary);
  else
    theFile.open(fileName, ios::out | ios::app | ios::binary);

  theOpenMode = APPEND;

  if (theFile.bad()) {
    std::cerr << "WARNING - BinaryFileStream::setFile()";
    std::cerr << " - could not open file " << fileName << std::endl;
    fileOpen = 0;
    return -1;
  } else
    fileOpen = 1;

  return 0;
}

int 
BinaryFileStream::close(void)
{
  if (fileOpen != 0)
    theFile.close();
  fileOpen = 0;

  return 0;
}


int 
BinaryFileStream::setPrecision(int prec)
{
  return 0;
}

int 
BinaryFileStream::setFloatField(floatField field)
{
  return 0;
}


int 
BinaryFileStream::tag(const char *tagName)
{
  return 0;
}

int 
BinaryFileStream::tag(const char *tagName, const char *value)
{
  return 0;
}


int 
BinaryFileStream::endTag()
{
  return 0;
}

int 
BinaryFileStream::attr(const char *name, int value)
{
  return 0;
}

int 
BinaryFileStream::attr(const char *name, double value)
{
  return 0;
}

int 
BinaryFileStream::attr(const char *name, const char *value)
{
  return 0;
}

int 
BinaryFileStream::write(Vector &data)
{
  if (fileOpen == 0)
    this->open();

  (*this) << data;  

  return 0;
}



OPS_Stream& 
BinaryFileStream::write(const char *s,int n)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}

OPS_Stream& 
BinaryFileStream::write(const unsigned char*s,int n)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}

OPS_Stream& 
BinaryFileStream::write(const signed char*s,int n)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}

OPS_Stream& 
BinaryFileStream::write(const void *s, int n)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}

OPS_Stream& 
BinaryFileStream::operator<<(char c)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}

OPS_Stream& 
BinaryFileStream::operator<<(unsigned char c)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}

OPS_Stream& 
BinaryFileStream::operator<<(signed char c)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}

OPS_Stream& 
BinaryFileStream::operator<<(const char *s)
{
  if (fileOpen == 0)
    this->open();

  if (strcmp(s,"\n") == 0) {

    // note that we do the flush so that a "/n" before
    // a crash will cause a flush() - similar to what 
    if (fileOpen != 0) {
      theFile << '\n';
      theFile.flush();
    }
  }

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(const unsigned char *s)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(const signed char *s)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(const void *p)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(int n)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(unsigned int n)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(long n)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(unsigned long n)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(short n)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(unsigned short n)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(bool b)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(double n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write((char *)&n, 8);

  return *this;
}
OPS_Stream& 
BinaryFileStream::operator<<(float n)
{
  if (fileOpen == 0)
    this->open();

  return *this;
}


int 
BinaryFileStream::sendSelf(int commitTag, Channel &theChannel)
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
    opserr << "BinaryFileStream::sendSelf() - failed to send id data\n";
    return -1;
  }

  if (fileNameLength != 0) {
    Message theMessage(fileName, fileNameLength);
    if (theChannel.sendMsg(0, commitTag, theMessage) < 0) {
      opserr << "BinaryFileStream::sendSelf() - failed to send message\n";
      return -1;
    }
  }

  sendSelfCount++;
  
  return 0;
}

int 
BinaryFileStream::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID idData(2);

  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "BinaryFileStream::recvSelf() - failed to recv id data\n";
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
      opserr << "BinaryFileStream::recvSelf() - out of memory\n";
      return -1;
    }

    Message theMessage(fileName, fileNameLength);
    if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
      opserr << "BinaryFileStream::recvSelf() - failed to recv message\n";
      return -1;
    }
    sprintf(&fileName[fileNameLength],".%d",commitTag);

    if (this->setFile(fileName, theOpenMode) < 0) {
      opserr << "BinaryFileStream::BinaryFileStream() - setFile() failed\n";
      if (fileName != 0) {
	delete [] fileName;
	fileName = 0;
      }
    }
  }
  
  return 0;
}


int 
binaryToText(const char *inputFilename, const char *outputFilename)
{
  //
  // open the files
  //

  ifstream input(inputFilename, ios::in | ios::binary);
  ofstream output(outputFilename, ios::out);
  if (input.bad()) {
    std::cerr << "WARNING - BinaryFileStream - binaryToText()";
    std::cerr << " - could not open file " << inputFilename << std::endl;
    output.close();
    return -1;
  }
  if (output.bad()) {
    std::cerr << "WARNING - BinaryFileStream - binaryToText()";
    std::cerr << " - could not open file " << outputFilename << std::endl;
    output.close();
    return -1;
  }

  //
  // until done
  //   read input consisting doubles till \n and write to output file
  //

  double data;
  char *c = (char *)&data;
  int numNumbers = 0;

  while ( !input.eof()) {
    input.read(c, 1);
    if (*c != '\n') {
      input.read(&c[1],7);
      output << data << " ";
      numNumbers++;
    } else {
      if (numNumbers != 0)
	output << "\n";
      numNumbers = 0;
    }
  }
  
  // 
  // close the files
  //

  input.close();
  output.close();

  return 0;
}

int 
textToBinary(const char *inputFilename, const char *outputFilename)
{
  //
  // open the files
  //

  ifstream input(inputFilename, ios::in);
  ofstream output(outputFilename, ios::out | ios::binary);
  if (input.bad()) {
    std::cerr << "WARNING - BinaryFileStream - binaryToText()";
    std::cerr << " - could not open file " << inputFilename << std::endl;
    output.close();
    return -1;
  }
  if (output.bad()) {
    std::cerr << "WARNING - BinaryFileStream - binaryToText()";
    std::cerr << " - could not open file " << outputFilename << std::endl;
    output.close();
    return -1;
  }

  //
  // until done
  //   read input consisting doubles till \n and write to output file
  //

  char data[100];
  char *dataNext;
  double d;

  while ( !input.eof()) {
    string inputLine;
    getline(input, inputLine);
    const char *c = inputLine.data();
    const char *cNext = c;

    int loc = 0;
    int endLoc = inputLine.length();
    int numNumbers = 0;

    while (loc < endLoc) {
      
      int dataCount = 0;


      while ((loc < endLoc) && 
	     (*cNext != ' ') && 
	     (*cNext != '\n')) {
	data[dataCount++] = cNext[0];
	cNext++;
	loc++;
      }
      
      if (dataCount != 0) {
	data[dataCount] = '\n';
	d = strtod(&data[0], &dataNext);
	output.write((char *)&d, 8);
	numNumbers++;
      }
      
      cNext++;
      loc++;
    }
    
    if (numNumbers != 0)
      output << '\n';
  }
  
  // 
  // close the files
  //

  input.close();
  output.close();

  return 0;
}

