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
// $Date: 2004-11-24 22:45:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/PatternRecorder.cpp,v $

// Written: MHS 
// Created: 2002
//
// Description: This file contains the class definition for PatternRecorder.
// A PatternRecorder is used to record the load values from a pattern.

#include <stdlib.h>
#include <string.h>

#include <PatternRecorder.h>
#include <Domain.h>
#include <LoadPattern.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <FE_Datastore.h>

#include <fstream>
using std::ifstream;

#include <iomanip>
using std::ios;

PatternRecorder::PatternRecorder(int pattern,
				 Domain &theDom,
				 const char *theFileName,
				 double dT,
				 double rTolDt,
				 int startFlag)
  :Recorder(RECORDER_TAGS_PatternRecorder),
   thePattern(pattern), theDomain(&theDom),
   flag(startFlag), deltaT(dT), relDeltaTTol(rTolDt), nextTimeStampToRecord(0.0)
{
  // create char array to store file name
  int fileNameLength = strlen(theFileName) + 1;
  fileName = new char[fileNameLength];
  if (fileName == 0) {
    opserr << "PatternRecorder::PatternRecorder - out of memory creating string" << endln;
  }
  
  // copy the strings
  strcpy(fileName, theFileName);    
  
  // open the file
  theFile.open(fileName, ios::out);
  if (theFile.bad()) {
    opserr << "WARNING - PatternRecorder::PatternRecorder()";
    opserr << " - could not open file " << fileName << endln;
  }
  //theFile << setiosflags(ios::scientific) << setprecision(15);
}

PatternRecorder::~PatternRecorder()
{
  if (!theFile)
    theFile.close();
  
  if (fileName != 0)
    delete [] fileName;
}

int 
PatternRecorder::record(int commitTag, double timeStamp)
{
  double value = 0.0;
  
  // where relDeltaTTol is the maximum reliable ratio between analysis time step and deltaT
  // and provides tolerance for floating point precision (see floating-point-tolerance-for-recorder-time-step.md)
    if (deltaT == 0.0 || timeStamp - nextTimeStampToRecord >= -deltaT * relDeltaTTol) {

    if (deltaT != 0.0)
      nextTimeStampToRecord = timeStamp + deltaT;
    
    LoadPattern *pattern = theDomain->getLoadPattern(thePattern);
    if (pattern != 0) {
      value = pattern->getLoadFactor();
    }
  }
  
  // write them to the file
  if (flag == 1)
    theFile << timeStamp << " ";
  else if (flag == 2)
    theFile << timeStamp << " ";
  
  theFile << value << " ";
  
  theFile << endln;
  theFile.flush();
  
  return 0;
}

int 
PatternRecorder::playback(int commitTag)
{
  if (theFile.bad())
    return 0;
  
  // close o/p file to ensure all buffered data gets written to file
  theFile.close(); 
  
  // open a stream for reading from the file
  ifstream inputFile;
  inputFile.open(fileName, ios::in);
  if (inputFile.bad()) {
    opserr << "WARNING - PatternRecorder::playback() - could not open file ";
    opserr << fileName << endln;
    return -1;
  }   
  
  double data;
  // read file up until line we want
  for (int i=0; i<(commitTag-1); i++)
    // now read in a line
    if (flag == 1 || flag == 2) {
      inputFile >> data;
      inputFile >> data;
    }
  
  // now read in our line and print out
  if (flag == 1 || flag == 2) {
    inputFile >> data;
    opserr << data << " ";
    inputFile >> data;
    opserr << data << " ";
    opserr << endln;
  }
  inputFile.close();
  
  // open file again for writing
  theFile.open(fileName, ios::app);
  if (theFile.bad()) {
    opserr << "WARNING - PatternRecorder::playback() - could not open file ";
    opserr << fileName << endln;
    return -1;
  }    
  
  // does nothing
  return 0;
}

int
PatternRecorder::restart(void)
{
  theFile.close();
  theFile.open(fileName, ios::out);
  if (theFile.bad()) {
    opserr << "WARNING - PatternRecorder::restart() - could not open file ";
    opserr << fileName << endln;
  }
  return 0;
}
