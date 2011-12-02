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
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataOutputFileHandler.cpp,v $
                                                                        
// Written: fmk 
// Date: 10/04
//
// Description: This file contains the class implementation for
// DataOutputFileHandler. 
//
// What: "@(#) DataOutputFileHandler.C, revA"

#include "DataOutputFileHandler.h"
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <Message.h>

DataOutputFileHandler::DataOutputFileHandler(const char *theFileName, 
					     echoMode theEMode, 
					     openMode theOMode)
  :DataOutputHandler(DATAHANDLER_TAGS_DataOutputFileHandler),
   fileName(0), theEchoMode(theEMode), theOpenMode(theOMode), numColumns(-1)
{
  if (theFileName != 0) {
    fileName = new char [strlen(theFileName)+1];
    if (fileName == 0) {
    opserr << "DataOutputFileHandler::DataOutputFileHandler() - out of memory\n";      
    } else 
      strcpy(fileName, theFileName);
  }

  if (fileName != 0 && outputFile.setFile(fileName, theOpenMode) < 0) {
    opserr << "DataOutputFileHandler::DataOutputFileHandler() - setFile() failed\n";
    if (fileName != 0) {
      delete [] fileName;
      fileName = 0;
    }
  }
}

DataOutputFileHandler::~DataOutputFileHandler()
{
  if (fileName != 0)
    delete [] fileName;
}

int 
DataOutputFileHandler::open(char **dataDescription, int numData)
{
  if (fileName == 0) {
    opserr << "DataOutputFileHandler::open() - no filename passed in constructor\n";
    return -1;
  }

  if (dataDescription == 0) {
    opserr << "DataOutputFileHandler::open() - no column description data passed\n";
    return -1;
  }

  if (numData < 0) {
    opserr << "DataOutputFileHandler::open() - numColumns (" << numData << ") < 0\n";
    return -1;
  } else
    numColumns = numData;


  if (theEchoMode == DATA_FILE) {
    for (int i=0; i<numData; i++)
      outputFile << dataDescription[i] << " ";
    outputFile << endln;

  } else if (theEchoMode == XML_FILE) {

    // create a copy of the  file name
    int res = 0;
    const char *name = outputFile.getFileName();
    int fileNameLength = strlen(name);
    char *xmlFileName = new char[fileNameLength + 5];
    
    if (xmlFileName == 0) {
      opserr << "DataOutputFileHandler::open - out of memory creating copy of string " << xmlFileName << endln;
      return -1;
  }
    
    strcpy(xmlFileName, name);
    if (fileNameLength > 4 && strcmp(".out", &name[fileNameLength-4]) == 0) {
      xmlFileName[fileNameLength-4] = '\0';
    }
    
    strcat(xmlFileName,".xml");
    
    FileStream xmlFile;
    if (xmlFile.setFile(xmlFileName, OVERWRITE) == 0) {
      
      // write the xml data
      xmlFile << "<?xml version=\"1.0\"?>\n";
      xmlFile << "<NumericalFileDataDescription>\n";
      xmlFile << "\t<DataFile>\n";
      xmlFile << "\t\t<DataFileName> " << name << "</DataFileName>\n";
      xmlFile << "\t\t<NumberDataColumns> " << numData << "</NumberDataColumns>\n";
      xmlFile << "\t</DataFile>\n";
      for (int i=0; i<numData; i++) {
	xmlFile << "\t<DataColumnDescription>\n";
	xmlFile << "\t\t<ColumnLocation> " << i+1 << "</ColumnLocation>\n";      
	xmlFile << "\t\t<Description> " << dataDescription[i] << "</Description>\n";
	xmlFile << "\t</DataColumnDescription>\n";
      }
      xmlFile << "</NumericalFileDataDescription>\n";
      xmlFile.close();
      
    } else {
      opserr << "DataOutputFileHandler::open - failed to open cml file: " << xmlFileName << endln;
      delete [] xmlFileName;
      res = -1;
    }
    
    // no longer need xmlFileName
    delete [] xmlFileName;
  }

  return 0;
}

int 
DataOutputFileHandler::write(Vector &data) 
{
  if (fileName == 0 || numColumns < 0) {
    opserr << "DataOutputFileHandler::write() - no filename or data description has been set\n";
    return -1;
  }

  if (data.Size() == numColumns)
    outputFile << data;
  else {
    opserr << fileName;
    opserr << "DataOutputStreamHandler::write() - Vector not of correct size\n";
    return -1;
  }

  return 0;
}
 
int 
DataOutputFileHandler::sendSelf(int commitTag, Channel &theChannel)
{
  static ID idData(3);
  int fileNameLength = 0;
  if (fileName != 0)
    fileNameLength = strlen(fileName);

  idData(0) = fileNameLength;

  if (theOpenMode == OVERWRITE)
    idData(1) = 0;
  else
    idData(1) = 1;

  if (theEchoMode == NONE)
    idData(2) = 0;
  else if (theEchoMode == DATA_FILE)
    idData(2) = 1;
  else
    idData(2) = 2;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "DataOutputFileHandler::sendSelf() - failed to send id data\n";
    return -1;
  }

  if (fileNameLength != 0) {
    Message theMessage(fileName, fileNameLength);
    if (theChannel.sendMsg(0, commitTag, theMessage) < 0) {
      opserr << "DataOutputFileHandler::sendSelf() - failed to send message\n";
      return -1;
    }
  }

  return 0;
}

int 
DataOutputFileHandler::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID idData(3);

  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "DataOutputFileHandler::recvSelf() - failed to recv id data\n";
    return -1;
  }

  int fileNameLength = idData(0);
  if (idData(1) == 0)
    theOpenMode = OVERWRITE;
  else
    theOpenMode = APPEND;


  if (idData(2) == 0)
    theEchoMode = NONE;
  else if (idData(2) == 1)
    theEchoMode = DATA_FILE;
  else
    theEchoMode = XML_FILE;

  if (fileNameLength != 0) {
    if (fileName != 0)
      delete [] fileName;
    fileName = new char[fileNameLength+5];
    if (fileName == 0) {
      opserr << "DataOutputFileHandler::recvSelf() - out of memory\n";
      return -1;
    }

    Message theMessage(fileName, fileNameLength);
    if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
      opserr << "DataOutputFileHandler::recvSelf() - failed to recv message\n";
      return -1;
    }
    sprintf(&fileName[fileNameLength],".%d",commitTag);

    if (fileName != 0 && outputFile.setFile(fileName, theOpenMode) < 0) {
      opserr << "DataOutputFileHandler::DataOutputFileHandler() - setFile() failed\n";
      if (fileName != 0) {
	delete [] fileName;
	fileName = 0;
      }
    }
  }
  return 0;
}
