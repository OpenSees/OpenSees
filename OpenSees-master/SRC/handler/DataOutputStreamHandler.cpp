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
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataOutputStreamHandler.cpp,v $
                                                                        
// Written: fmk 
// Date: 10/04
//
// Description: This file contains the class implementation for
// DataOutputStreamHandler. 
//
// What: "@(#) DataOutputStreamHandler.C, revA"

#include "DataOutputStreamHandler.h"
#include <Vector.h>

DataOutputStreamHandler::DataOutputStreamHandler(bool echoD)
  :DataOutputHandler(DATAHANDLER_TAGS_DataOutputStreamHandler), 
   echoDescription(echoD), numColumns(0)
{

}

DataOutputStreamHandler::~DataOutputStreamHandler()
{
    // does nothing
}

int 
DataOutputStreamHandler::open(char **dataDescription, int numData)
{
  if (numData < 0) {
    opserr << "DataOutputStreamHandler::open() - numColumns (" << numData << ") < 0\n";
    return -1;
  } else
    numColumns = numData;

  if (echoDescription == true) {
    for (int i=0; i<numData; i++)
      opserr << dataDescription[i] << " ";
    opserr << endln;
  }
  
  return 0;
}

int 
DataOutputStreamHandler::write(Vector &data) 
{
  if (data.Size() == numColumns)
    outputStream << data;
  else {
    opserr << "DataOutputStreamHandler::write() - Vector not of correct size\n";
    return -1;
  }
  
  return 0;
}

int 
DataOutputStreamHandler::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int 
DataOutputStreamHandler::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}
