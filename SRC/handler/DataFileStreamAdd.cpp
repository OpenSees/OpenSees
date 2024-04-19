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
// $Date: 2009-10-13 21:17:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataFileStreamAdd.cpp,v $


#include <DataFileStreamAdd.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <string.h>
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

DataFileStreamAdd::DataFileStreamAdd(int indent)
  :OPS_Stream(OPS_STREAM_TAGS_DataFileStreamAdd), 
   fileOpen(0), fileName(0), indentSize(indent), sendSelfCount(0), theChannels(0), numDataRows(0),
   mapping(0), maxCount(0), sizeColumns(0), theColumns(0), theData(0), theRemoteData(0), doCSV(0)
{
  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+5];
  for (int i=0; i<indentSize; i++)
    strcpy(indentString, " ");
}


DataFileStreamAdd::DataFileStreamAdd(const char *file, openMode mode, int indent, int csv, bool closeWrite, int prec, bool scientific)
  :OPS_Stream(OPS_STREAM_TAGS_DataFileStreamAdd), 
   fileOpen(0), fileName(0), indentSize(indent), sendSelfCount(0), 
   theChannels(0), numDataRows(0),
   mapping(0), maxCount(0), sizeColumns(0), 
   theColumns(0), theData(0), theRemoteData(0), 
   doCSV(csv), closeOnWrite(closeWrite)
{
  thePrecision = prec;
  doScientific = scientific;

  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+1];
  for (int i=0; i<indentSize; i++)
    strcpy(indentString, " ");

  this->setFile(file, mode);
}


DataFileStreamAdd::~DataFileStreamAdd()
{
  if (fileOpen == 1)
    theFile.close();

  if (theChannels != 0) {
    delete [] theChannels;
  }

  if (indentString != 0)
    delete [] indentString;

  if (fileName != 0)
    delete [] fileName;

  if (sendSelfCount > 0) {
    for (int i=0; i<=sendSelfCount; i++) {
      if (theColumns != 0)
	if (theColumns[i] != 0)
	  delete theColumns[i];

      if (theData != 0)
	if (theData[i] != 0)
	  delete [] theData[i];

      if (theRemoteData != 0)
	if (theRemoteData[i] != 0)
	  delete theRemoteData[i];

    }

    if (theData != 0) delete [] theData;
    if (theRemoteData != 0) delete [] theRemoteData;
    if (theColumns != 0) delete [] theColumns;
    if (sizeColumns != 0) delete sizeColumns;
  }    
}

int 
DataFileStreamAdd::setFile(const char *name, openMode mode)
{
  if (name == 0) {
    std::cerr << "DataFileStreamAdd::setFile() - no name passed\n";
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
      std::cerr << "DataFileStreamAdd::setFile() - out of memory copying name: " << name << std::endl;
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
DataFileStreamAdd::open(void)
{
  // check setFile has been called
  if (fileName == 0) {
    std::cerr << "DataFileStreamAdd::open(void) - no file name has been set\n";
    return -1;
  }

  // if file already open, return
  if (fileOpen == 1) {
    return 0;
  }

  if (theOpenMode == OVERWRITE) 
    theFile.open(fileName, ios::out);
  else
    theFile.open(fileName, ios::out| ios::app);

  theOpenMode = APPEND;

  if (theFile.bad()) {
    std::cerr << "WARNING - DataFileStreamAdd::setFile()";
    std::cerr << " - could not open file " << fileName << std::endl;
    fileOpen = 0;
    return -1;
  } else
    fileOpen = 1;

  if (doScientific == true)
    theFile << std::scientific;

  theFile << std::setprecision(thePrecision);

  return 0;
}

int 
DataFileStreamAdd::close(void)
{
  if (fileOpen != 0)
    theFile.close();
  fileOpen = 0;

  return 0;
}


int 
DataFileStreamAdd::setPrecision(int prec)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << std::setprecision(prec);

  return 0;
}

int 
DataFileStreamAdd::setFloatField(floatField field)
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
DataFileStreamAdd::tag(const char *tagName)
{
  return 0;
}

int 
DataFileStreamAdd::tag(const char *tagName, const char *value)
{
  return 0;
}


int 
DataFileStreamAdd::endTag()
{
  return 0;
}

int 
DataFileStreamAdd::attr(const char *name, int value)
{
  return 0;
}

int 
DataFileStreamAdd::attr(const char *name, double value)
{
  return 0;
}

int 
DataFileStreamAdd::attr(const char *name, const char *value)
{
  return 0;
}

int 
DataFileStreamAdd::write(Vector &data)
{
  if (fileOpen == 0)
    this->open();

  //
  // if not parallel, just write the data
  //

  if (sendSelfCount == 0) {
    (*this) << data;  
    if (closeOnWrite == true)
      this->close();
    return 0;
  }

  //
  // otherwise parallel, send the data if not p0
  //

  if (sendSelfCount < 0) {
    if (data.Size() != 0) {
      if ( theChannels[0]->sendVector(0, 0, data) < 0) {
	return -1;
      }
      return 0;
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
	Vector *theV = theRemoteData[i];
	if (theChannels[i-1]->recvVector(0, 0, *theV) < 0) {
	  opserr << "DataFileStreamAdd::write - failed to recv data\n";
	}
      } 
    }
  }

  Matrix &printMapping = *mapping;

  // write data
  ID currentCol(sendSelfCount+1);
  ID currentColLoc(sendSelfCount+1);
  for (int i=0; i<=sendSelfCount; i++) {
    if (theColumns[i] != 0) {
      currentCol(i) = (*theColumns[i])(0);
      currentColLoc(i) = 0;
    } else {
      currentCol(i) = -1;
      currentColLoc(i) = -1;
    }
  }

  //  opserr << "CurrentCOL: " << currentCol;
  //  opserr << "MaxCOUNT: " << maxCount << endln;

  for (int i=0; i<maxCount+1; i++) {
    double val = 0;
    int numData = (int)printMapping(2,i);
    for (int j=0; j<numData; j++) {
      val = 0.;
      for (int k=0; k<=sendSelfCount; k++) {
	if (theColumns[k] != 0) {
	  ID &currentColID = *theColumns[k];
	  if (currentCol(k) == i) {
	    double *theColData = theData[k];
	    val += theColData[currentColLoc(k)+j];
	  }
	}
      }
      theFile << val << " ";
    }
    for (int k=0; k<=sendSelfCount; k++) {
      if (theColumns[k] != 0) {
	if (currentCol(k) == i) {
	  currentCol(k) = currentCol(k)+1;
	  currentColLoc(k) = currentColLoc(k)+numData;
	}
      }
    }
  }
  theFile << "\n";

  /*
  if (doCSV == 0) {
    for (int i=0; i<maxCount+1; i++) {
      int fileID = (int)printMapping(0,i);
      int startLoc = (int)printMapping(1,i);
      int numData = (int)printMapping(2,i);
      double *data = theData[fileID];
      for (int j=0; j<numData; j++)
	theFile << data[startLoc++] << " ";
    }
    theFile << "\n";
  } else {
    for (int i=0; i<maxCount+1; i++) {
      int fileID = (int)printMapping(0,i);
      int startLoc = (int)printMapping(1,i);
      int numData = (int)printMapping(2,i);
      double *data = theData[fileID];
      int nM1 = numData-1;
      for (int j=0; j<numData; j++)
	if ((i ==maxCount) && (j == nM1))
	  theFile << data[startLoc++] << "\n";
	else
	  theFile << data[startLoc++] << ",";
    }
  }
  */

  if (closeOnWrite == true)
    this->close();
  
  return 0;
}



OPS_Stream& 
DataFileStreamAdd::write(const char *s,int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write(s, n);

  return *this;
}

OPS_Stream& 
DataFileStreamAdd::write(const unsigned char*s,int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
DataFileStreamAdd::write(const signed char*s,int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
DataFileStreamAdd::write(const void *s, int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}

OPS_Stream& 
DataFileStreamAdd::write(const double *s, int n)
{
  numDataRows++;

  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0) {
    if (n > 0) {
      if (doCSV == 0) {
	int nm1 = n-1;
	for (int i=0; i<nm1; i++) {
	  theFile << s[i] << " ";
	}
	theFile << s[nm1] << "\n";
      } else { int nm1 = n-1;
	for (int i=0; i<nm1; i++) {
	  theFile << s[i] << ",";
	}
	theFile << s[nm1] << "\n";
      }
    }
  }
  return *this;
}


OPS_Stream& 
DataFileStreamAdd::operator<<(char c)
{  
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(unsigned char c)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(signed char c)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(const char *s)
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
DataFileStreamAdd::operator<<(const unsigned char *s)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(const signed char *s)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(const void *p)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << p;
*/
  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(unsigned int n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(long n)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(unsigned long n)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(short n)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(unsigned short n)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(bool b)
{
  if (fileOpen == 0)
    this->open();

/*
  if (fileOpen != 0)
    theFile << b;
*/
  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(double n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
DataFileStreamAdd::operator<<(float n)
{
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << n;

  return *this;
}


int 
DataFileStreamAdd::sendSelf(int commitTag, Channel &theChannel)
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
    fileNameLength = int(strlen(fileName));

  idData(0) = fileNameLength;

  if (theOpenMode == OVERWRITE)
    idData(1) = 0;
  else
    idData(1) = 1;

  idData(2) = sendSelfCount;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "DataFileStreamAdd::sendSelf() - failed to send id data\n";
    return -1;
  }

  if (fileNameLength != 0) {
    Message theMessage(fileName, fileNameLength);
    if (theChannel.sendMsg(0, commitTag, theMessage) < 0) {
      opserr << "DataFileStreamAdd::sendSelf() - failed to send message\n";
      return -1;
    }
  }
  
  return 0;
}

int 
DataFileStreamAdd::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID idData(3);

  sendSelfCount = -1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "DataFileStreamAdd::recvSelf() - failed to recv id data\n";
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
    fileName = new char[fileNameLength+10];
    if (fileName == 0) {
      opserr << "DataFileStreamAdd::recvSelf() - out of memory\n";
      return -1;
    }

    Message theMessage(fileName, fileNameLength);
    if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
      opserr << "DataFileStreamAdd::recvSelf() - failed to recv message\n";
      return -1;
    }

    int tag = idData(2);

    sprintf(&fileName[fileNameLength],".%d",tag);

    if (this->setFile(fileName, theOpenMode) < 0) {
      opserr << "DataFileStreamAdd::DataFileStreamAdd() - setFile() failed\n";
      if (fileName != 0) {
	delete [] fileName;
	fileName = 0;
      }
    }
  }
  
  return 0;
}


void
DataFileStreamAdd::indent(void)
{
  if (fileOpen != 0)
    for (int i=0; i<numIndent; i++)
      theFile << indentString;
}


int
DataFileStreamAdd::setOrder(const ID &orderData)
{
  if (sendSelfCount == 0)
    return 0;

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

    for (int i=0; i<=sendSelfCount; i++) {
      (*sizeColumns)(i) = 0;
      theColumns[i] = 0;
      theData[i] = 0;
      theRemoteData[i] = 0;
    }
    
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
	opserr << "DataFileStreamAdd::setOrder - failed to recv column size for process: " << i+1 << endln;
	return -1;
      }

      int numColumns = numColumnID(0);

      (*sizeColumns)(i+1) = numColumns;
      if (numColumns != 0) {
	theColumns[i+1] = new ID(numColumns);
	if (theChannels[i]->recvID(0, 0, *theColumns[i+1]) < 0) {
	  opserr << "DataFileStreamAdd::setOrder - failed to recv column data for process: " << i+1 << endln;
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
  }

  return 0;
}

int DataFileStreamAdd::flush() {
  if (theFile.is_open() && theFile.good()) {
    theFile.flush();
  }
  return 0;
}
