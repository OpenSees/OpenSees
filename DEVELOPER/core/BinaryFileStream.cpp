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
// $Date: 2009-04-30 23:23:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/BinaryFileStream.cpp,v $

#include <BinaryFileStream.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <ID.h>
#include <Channel.h>
#include <Message.h>
#include <Matrix.h>

using std::cerr;
using std::ios;
using std::setiosflags;
using std::ifstream;
using std::string;
using std::getline;

BinaryFileStream::BinaryFileStream()
  :OPS_Stream(OPS_STREAM_TAGS_BinaryFileStream), 
   fileOpen(0), fileName(0), sendSelfCount(0),
   theChannels(0), numDataRows(0),
   mapping(0), maxCount(0), sizeColumns(0), theColumns(0), theData(0), theRemoteData(0)
{

}

BinaryFileStream::BinaryFileStream(const char *file, openMode mode)
  :OPS_Stream(OPS_STREAM_TAGS_BinaryFileStream), 
   fileOpen(0), fileName(0), sendSelfCount(0),
   theChannels(0), numDataRows(0),
   mapping(0), maxCount(0), sizeColumns(0), theColumns(0), theData(0), theRemoteData(0)
{
  this->setFile(file, mode);
}

BinaryFileStream::~BinaryFileStream()
{
  if (fileOpen == 1)
    theFile.close();

  if (theChannels != 0) {

    static ID lastMsg(1);
    if (sendSelfCount > 0) {
      for (int i=0; i<sendSelfCount; i++) 
	theChannels[i]->sendID(0, 0, lastMsg);
    } else
	theChannels[0]->recvID(0, 0, lastMsg);
    delete [] theChannels;
  }

  if (fileName != 0)
    delete [] fileName;

  if (sendSelfCount > 0) {

    for (int i=0; i<=sendSelfCount; i++) {
      if (theColumns[i] != 0)
	delete theColumns[i];

      if (theData[i] != 0)
	delete [] theData[i];

      if (theRemoteData[i] != 0)
	delete theRemoteData[i];
    }
    delete [] theData;
    delete [] theRemoteData;
    delete [] theColumns;
    delete sizeColumns;
  }    
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

  //
  // if not parallel, just write the data
  //

  if (sendSelfCount == 0) {
    (*this) << data;  
    return 0;
  }

  //
  // otherwise parallel, send the data if not p0
  //

  if (sendSelfCount < 0) {
    if (data.Size() != 0) {
      return theChannels[0]->sendVector(0, 0, data);
    } else
      return 0;
  }

  //
  // if p0 recv the data & write it out sorted
  //

  // recv data
  for (int i=0; i<=sendSelfCount; i++) {
    int numColumns = (*sizeColumns)(i);
    double *dataI = theData[i];
    if (i == 0) {
      for (int j=0; j<numColumns; j++) {
	dataI[j] = data(j);
      }
    } else { 
      if (numColumns != 0) {
	theChannels[i-1]->recvVector(0, 0, *(theRemoteData[i]));
      }
    }
  }

  Matrix &printMapping = *mapping;

  // write data
  for (int i=0; i<maxCount+1; i++) {
    int fileID = (int)printMapping(0,i);
    int startLoc = (int)printMapping(1,i);
    int numData = (int)printMapping(2,i);
    double *data = theData[fileID];
    //    for (int j=0; j<numData; j++)
    theFile.write((char *)(&data[startLoc]), 8*numData);
  }
  theFile << "\n";

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
BinaryFileStream::write(const double *s, int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0) {
    //    for (int i=0; i<n; i++)
    theFile.write((char *)(&s[0]), 8*n);

    theFile << '\n';
    theFile.flush();
  }
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
  sendSelfCount++;

  Channel **theNextChannels = new Channel *[sendSelfCount];
  for (int i=0; i<sendSelfCount-1; i++)
    theNextChannels[i] = theChannels[i];
  theNextChannels[sendSelfCount-1] = &theChannel;
  if (theChannels != 0)
    delete [] theChannels;
  theChannels = theNextChannels;

  static ID idData(3);
  int fileNameLength = 0;
  if (fileName != 0)
    fileNameLength = strlen(fileName);

  idData(0) = fileNameLength;

  if (theOpenMode == OVERWRITE)
    idData(1) = 0;
  else
    idData(1) = 1;

  idData(2) = sendSelfCount;

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
  
  return 0;
}

int 
BinaryFileStream::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID idData(3);

  sendSelfCount = -1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

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

    int tag = idData(2);

    sprintf(&fileName[fileNameLength],".%d",tag);

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
  /* ORIGINAL
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
  REPLACED WITH: */

  int nLF = 0;
  int intervalLF = 0;
  double aveIntervalLF = 0;
  int numCol;

  int dataLen = 0;
  int numRow;

  while (!input.eof()) {
    input.read(c, 1);
    dataLen++;

    if ((*c == '\n') && (intervalLF%8 == 0)) {    
      aveIntervalLF = (aveIntervalLF * nLF + intervalLF)/(nLF + 1);
      nLF++;
      intervalLF = 0;
    } else {
      intervalLF++;      
    }

  }
  numCol = (int)(aveIntervalLF/8.0 + 0.5);
  numRow = (int)(dataLen/(8.0*numCol+1.0) + 0.5);

  // rewind
  input.clear();
  input.seekg(0,std::ios::beg); 

  // read, output
  for (int ii=0; ii<numRow; ii++) {
    for (int jj=0; jj<numCol; jj++) {

      input.read(&c[0],8);
      output << data ;
      if (jj<(numCol-1)) output << " ";

    }

    input.read(c, 1);//LF
    output << "\n";

  }
  /* END REPLACEMENT BLOCK */

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

int
BinaryFileStream::setOrder(const ID &orderData)
{
  if (sendSelfCount < 0) {
    static ID numColumnID(1);
    int numColumn = orderData.Size();
    numColumnID(0) = numColumn;
    theChannels[0]->sendID(0, 0, numColumnID);
    if (numColumn != 0)
      theChannels[0]->sendID(0, 0, orderData);
  }

  if (sendSelfCount > 0) {      

    sizeColumns = new ID(sendSelfCount+1);
    theColumns = new ID *[sendSelfCount+1];
    theData = new double *[sendSelfCount+1];
    theRemoteData = new Vector *[sendSelfCount+1];
    
    int numColumns = orderData.Size();
    (*sizeColumns)(0) = numColumns;
    if (numColumns != 0) {
      theColumns[0] = new ID(orderData);
      theData[0] = new double [numColumns];
    } else {
      theColumns[0] = 0;
      theData[0] = 0;
    }      
    theRemoteData[0] = 0;

    maxCount = 0;
    if (numColumns != 0)
      maxCount = orderData(numColumns-1);

    // now receive orderData from the other channels
    for (int i=0; i<sendSelfCount; i++) { 
      static ID numColumnID(1);	  
      if (theChannels[i]->recvID(0, 0, numColumnID) < 0) {
	opserr << "BinaryFileStream::setOrder - failed to recv column size for process: " << i+1 << endln;
	return -1;
      }

      int numColumns = numColumnID(0);

      (*sizeColumns)(i+1) = numColumns;
      if (numColumns != 0) {
	theColumns[i+1] = new ID(numColumns);
	if (theChannels[i]->recvID(0, 0, *theColumns[i+1]) < 0) {
	  opserr << "BinaryFileStream::setOrder - failed to recv column data for process: " << i+1 << endln;
	  return -1;
	}
	
	if (numColumns != 0 && (*theColumns[i+1])[numColumns-1] > maxCount)
	  maxCount = (*theColumns[i+1])[numColumns-1];
	
	theData[i+1] = new double [numColumns];
	theRemoteData[i+1] = new Vector(theData[i+1], numColumns);
      } else {
	theColumns[i+1] = 0;
	theData[i+1] = 0;
	theRemoteData[i+1] = 0;
      }
    }

    ID currentLoc(sendSelfCount+1);
    ID currentCount(sendSelfCount+1);
	
    if (mapping != 0)
      delete mapping;

    mapping = new Matrix(3, maxCount+1);

    Matrix &printMapping = *mapping;
	
    for (int i=0; i<=sendSelfCount; i++) {
      currentLoc(i) = 0;
      if (theColumns[i] != 0)
	currentCount(i) = (*theColumns[i])[0];
      else
	currentCount(i) = -1;
    }

    int count =0;
    while (count <= maxCount) {
      for (int i=0; i<=sendSelfCount; i++) {
	if (currentCount(i) == count) {
	  printMapping(0,count) = i;
	  
	  int maxLoc = theColumns[i]->Size();
	  int loc = currentLoc(i);
	  int columnCounter = 0;
	  
	  printMapping(1,count) = loc;
	  
	  while (loc < maxLoc && (*theColumns[i])(loc) == count) {
	    loc++;
	    columnCounter++;
	  }
	  
	  printMapping(2,count) = columnCounter;
	  
	  currentLoc(i) = loc;
	  
	  if (loc < maxLoc)
	    currentCount(i) = (*theColumns[i])(loc);		
	  else
	    currentCount(i) = -1; 		
	}
      }
      count++;
    }
    
    opserr << printMapping;
  }

  return 0;
}
