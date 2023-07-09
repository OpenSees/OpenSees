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
                                                                        
// $Revision: 1.12 $
// $Date: 2010-02-04 00:33:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/PathTimeSeriesThermal.cpp,v $
                                                                        
//This file is based on the pathtimeSeries
// Purpose: This file contains the class definition for PathTimeSeriesThermal.
// PathTimeSeriesThermal is a concrete class. A PathTimeSeriesThermal object provides
// a linear time series. the factor is given by the pseudoTime and 
// a constant factor provided in the constructor.
//
// What: "@(#) PathTimeSeriesThermal.C, revA"
                                                             
 //Modified by Liming Jiang [http://openseesforfire.github.io]

#include <PathTimeSeriesThermal.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <math.h>

#include <fstream>
using std::ifstream;

#include <iomanip>
using std::ios;

PathTimeSeriesThermal::PathTimeSeriesThermal()	
  :TimeSeries(TSERIES_TAG_PathTimeSeriesThermal),CurrentFactors(0),
   thePath(0), time(0), currentTimeLoc(0),numCols(0),numRows(0), 
   cFactor(0.0), dbTag1(0), dbTag2(0), lastSendCommitTag(-1), TempOut(true)
{
  // does nothing
}



PathTimeSeriesThermal::PathTimeSeriesThermal(int tag,
			       const char *fileName, int DataNum, bool tempOut,
			       double theFactor)
  :TimeSeries(tag, TSERIES_TAG_PathTimeSeriesThermal),CurrentFactors(0),
   thePath(0), time(0), currentTimeLoc(0), numCols(DataNum),numRows(0), 
   cFactor(theFactor), dbTag1(0), dbTag2(0), lastChannel(0), TempOut(tempOut)
{


  // determine the number of data points
  int numDataPoints = 0;
  int numRows = 0;
  double dataPoint;

  
  ifstream theFile;
  
  // first open and go through file counting entries
  theFile.open(fileName, ios::in);
  if (theFile.bad() || !theFile.is_open()) {
    opserr << "WARNING - PathTimeSeriesThermal::PathTimeSeriesThermal()";
    opserr << " - could not open file " << fileName << endln;
  }
  else {
    while (theFile >> dataPoint) {
      numDataPoints++;
      //theFile >> dataPoint;	// Read in second value of pair
    }
  }

  if ((numDataPoints % (numCols+1)) != 0) {
    opserr << "WARNING - PathTimeSeriesThermal::PathTimeSeriesThermal()";
    opserr << " - num data entries in file NOT COMPATIBLE! " << fileName << endln;
  }
  numRows= numDataPoints/(numCols+1);


  theFile.close();

  // create a vector and read in the data
  if (numRows != 0) {
    
    // now create the two vector
    thePath = new Matrix(numRows,numCols);
    time = new Vector(numRows);
	CurrentFactors = new Vector(numCols);
    
    // ensure did not run out of memory creating copies
    if (thePath == 0 || thePath->noCols() == 0 || thePath->noRows() == 0 || time == 0 || time->Size() == 0) {
      
      opserr << "WARNING PathTimeSeriesThermal::PathTimeSeriesThermal() - out of memory\n ";
      if (thePath != 0)
	delete thePath;
      if (time != 0)
	delete time;
      thePath = 0;
      time = 0;
    }
    
    // first open the file and read in the data
    ifstream theFile1;
    theFile1.open(fileName, ios::in);
    if (theFile1.bad() || !theFile1.is_open()) {
      
      opserr << "WARNING - PathTimeSeriesThermal::PathTimeSeriesThermal()";
      opserr << " - could not open file " << fileName << endln;
      delete thePath;
      delete time;
      thePath = 0;
      time =0;
    }
    else { // read in the time and then read the value
      int RowCount = 0;
      while (theFile1 >> dataPoint) {
	(*time)(RowCount) = dataPoint;
	for(int j=0;j<numCols;j++){
	theFile1 >> dataPoint;
	if(TempOut)
	(*thePath)(RowCount,j) = dataPoint-20;
	else
	(*thePath)(RowCount,j) = dataPoint;

	}

	RowCount++;
      }
      
      // finally close the file
      theFile1.close();
    } 
  }
}


PathTimeSeriesThermal::PathTimeSeriesThermal(int tag,
                                             int DataNum, bool tempOut,
                                             double theFactor)
:TimeSeries(tag, TSERIES_TAG_PathTimeSeriesThermal),CurrentFactors(0),
thePath(0), time(0), currentTimeLoc(0),
cFactor(theFactor), dbTag1(0), dbTag2(0), lastChannel(0), TempOut(tempOut)
{
  
  numRows = 1;
  numCols = DataNum;

 
  
  
  // create a vector and read in the data
  if (numRows != 0) {
    
    // now create the two vector
    thePath = new Matrix(numRows,numCols);
    time = new Vector(numRows);
	CurrentFactors = new Vector(numCols);
    
    // ensure did not run out of memory creating copies
    if (thePath == 0 || thePath->noCols() == 0 || thePath->noRows() == 0 || time == 0 || time->Size() == 0) {
      
      opserr << "WARNING PathTimeSeriesThermal::PathTimeSeriesThermal() - out of memory\n ";
      if (thePath != 0)
	delete thePath;
      if (time != 0)
	delete time;
      thePath = 0;
      time = 0;
    }
  }

  (*time)(numRows-1) = 0;

  for (int i =0; i<numCols; i++) {
    (*thePath)(numRows-1, i) = 0.0;
  }
  
  //update it later
  
  
}



PathTimeSeriesThermal::~PathTimeSeriesThermal()
{
  if (thePath != 0)
    delete thePath;
  if (time != 0)
    delete time;
}


int
PathTimeSeriesThermal::WriteResults(double currentTime, const Vector& newData)
{
  numRows++;
 
  Vector BuffTime = *time;
  time->resize(numRows);
  for(int i=0;i<numRows;i++){
	  if(i<numRows-1)
		  (*time)(i)=BuffTime(i);
	  else
		  (*time)(i)= currentTime;
  }
#ifdef _DEBUG
 // opserr<<"PathTimeSeries   "<<newData<<endln;
#endif
  if (newData==0||newData.Size()<numCols) {
    opserr<<"WARNING::PathTimeSeriesThermal received incompatible data when attempring to write the results"<<endln;
    return -1;
  }
  
  Matrix BuffPath = (*thePath);
  thePath->resize(numRows,numCols);
  for(int j=0; j<numRows;j++){
	for (int i =0; i<numCols; i++) {
		if(j<numRows-1)
			(*thePath)(j, i) = BuffPath(j,i);
		else 
			(*thePath)(j, i) = newData(i);
	}
  }
  
  return 0;
}

TimeSeries *
PathTimeSeriesThermal::getCopy(void) 
{
opserr<<"PathTimeSeriesThermal::getCopy() should not be called"<<endln;  
return new PathTimeSeriesThermal();
}


double
PathTimeSeriesThermal::getTimeIncr (double pseudoTime)
{
  // NEED TO FILL IN, FOR NOW return 1.0
  return 1.0;
}

const Vector&
PathTimeSeriesThermal::getFactors(double pseudoTime)
{
  // check for a quick return
  if (thePath == 0)
    return 0.0;

  //opserr<<"PathTimeSeries Tag"<<this->getTag()<<endln;
  // determine indexes into the data array whose boundary holds the time
  double time1 = (*time)(currentTimeLoc);
  // check for another quick return
  if (pseudoTime == time1){
	  for(int j=0;j<numCols;j++){
		  (*CurrentFactors)(j)=cFactor * (*thePath)(currentTimeLoc,j);
	  }
	  return *CurrentFactors;
  }

  int size = time->Size();
  int sizem1 = size - 1;
  int sizem2 = size - 2;
  
  // check we are not at the end
  if (pseudoTime > time1 && currentTimeLoc == sizem1){
    for(int j=0;j<numCols;j++){
		  (*CurrentFactors)(j)=0.0;
	  }
	  return *CurrentFactors;
  }

  if (pseudoTime < time1 && currentTimeLoc == 0){
    for(int j=0;j<numCols;j++){
		  (*CurrentFactors)(j)=(*thePath)(currentTimeLoc,j)*pseudoTime/time1;
	  }
	  return *CurrentFactors;
  }

  // otherwise go find the current interval
  double time2 = (*time)(currentTimeLoc+1);
  if (pseudoTime > time2) {
    while ((pseudoTime > time2) && (currentTimeLoc < sizem2)) {
      currentTimeLoc++;
      time1 = time2;
      time2 = (*time)(currentTimeLoc+1);
    }
    // if pseudo time greater than ending time return 0
    if (pseudoTime > time2){
    for(int j=0;j<numCols;j++){
		  (*CurrentFactors)(j)=0.0;
	  }
	  return *CurrentFactors;
    }

  } else if (pseudoTime < time1) {
    while ((pseudoTime < time1) && (currentTimeLoc > 0)) {
      currentTimeLoc--;
      time2 = time1;	
      time1 = (*time)(currentTimeLoc);
    }
    // if starting time less than initial starting time return 0
    if (pseudoTime < time1){
	  for(int j=0;j<numCols;j++){
		//(*CurrentFactors)(j)=0.0;
		(*CurrentFactors)(j)= (*thePath)(currentTimeLoc,j)*pseudoTime/time1;
	  }
	  return *CurrentFactors;
    }
  }
  
  for(int j=0;j<numCols;j++){
   double value1 = (*thePath)(currentTimeLoc,j);
   double value2 = (*thePath)(currentTimeLoc+1,j);
  (*CurrentFactors)(j)= cFactor*(value1 + (value2-value1)*(pseudoTime-time1)/(time2 - time1));
  }
  return *CurrentFactors;
}

double
PathTimeSeriesThermal::getDuration()
{
  if (thePath == 0)
  {
    opserr << "WARNING -- PathTimeSeriesThermal::getDuration() on empty Vector" << endln;
	return 0.0;
  }

  int lastIndex = time->Size(); // index to last entry in time vector
  return ((*time)[lastIndex-1]);
}

double
PathTimeSeriesThermal::getPeakFactor()
{
  if (thePath == 0)
  {
    opserr << "WARNING -- PathTimeSeriesThermal::getPeakFactor() on empty Vector" << endln;
	return 0.0;
  }

  
  return 0.0;
}

int
PathTimeSeriesThermal::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  Vector data(5);
  data(0) = cFactor;
  data(1) = -1;
  
  if (thePath != 0) {
	  int size = thePath->noRows();
    data(1) = size;
    if (dbTag1 == 0) {
      dbTag1 = theChannel.getDbTag();
      dbTag2 = theChannel.getDbTag();
    }
    data(2) = dbTag1;
    data(3) = dbTag2;
  }

  if ((lastSendCommitTag == -1) && (theChannel.isDatastore() == 1)) {
    lastSendCommitTag = commitTag;
  }

  data(4) = lastSendCommitTag;
  
  int result = theChannel.sendVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "PathTimeSeriesThermal::sendSelf() - channel failed to send data\n";
    return result;
  }

  // we only send the vector data if this is the first time it is sent to the database
  // or the channel is for sending the data to a remote process

  if (lastChannel != &theChannel || (lastSendCommitTag == commitTag) || (theChannel.isDatastore() == 0)) {  

    lastChannel = &theChannel;

    if (thePath != 0) {
		opserr<<"Remote access to PathTimeSeriesThermal(sendSelf) is deactivated by UoE Group"<<endln;
		/*     
		result = theChannel.sendVector(dbTag1, commitTag, *thePath);
      if (result < 0) {
	opserr << "PathTimeSeriesThermal::sendSelf() - ";
	opserr << "channel failed to send the Path Vector\n";
	return result; 
	
      }
	  */
    }
    if (time != 0) {
      result = theChannel.sendVector(dbTag2, commitTag, *time);
      if (result < 0) {
	opserr << "PathTimeSeriesThermal::sendSelf() - ";
	opserr << "channel failed to send the Path Vector\n";
	return result;  
      }
    }
    return 0;
  }
  return 0;
}


int 
PathTimeSeriesThermal::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector data(5);
  int result = theChannel.recvVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "PathTimeSeriesThermal::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    return result;
  }
  cFactor = data(0);
  int size = data(1);
  lastSendCommitTag = data(4);
/*
  // get the data cvector, only receive them once as they cannot change
  if (thePath == 0 && size > 0) {
    dbTag1 = data(2);
    dbTag2 = data(3);
    thePath = new Vector(size);
    time = new Vector(size);
    if (thePath == 0 || time == 0 ||
	thePath->Size() == 0 || time->Size() == 0) {

      opserr << "PathTimeSeriesThermal::recvSelf() - ran out of memory";
      opserr << " a Vector of size: " <<  size << endln;  
      if (thePath != 0)
	delete thePath;
      if (time != 0)
	delete time;
      thePath = 0;
      time = 0;
      return -1;
    }
    result = theChannel.recvVector(dbTag1, lastSendCommitTag, *thePath);    
    if (result < 0) {
      opserr << "PathTimeSeriesThermal::recvSelf() - ";
      opserr << "channel failed to receive the Path Vector\n";
      return result;  
    }
    result = theChannel.recvVector(dbTag2, lastSendCommitTag, *time);    
    if (result < 0) {
      opserr << "PathTimeSeriesThermal::recvSelf() - ";
      opserr << "channel failed to receive the time Vector\n";
      return result;  
    }
  }
  */
  return 0;    
}



void
PathTimeSeriesThermal::Print(OPS_Stream &s, int flag)
{
    s << "Path Time Series: constant factor: " << cFactor;
    if (flag == 1 && thePath != 0) {
      s << " specified path: " << *thePath;
      s << " specified time: " << *time;
    }
}
