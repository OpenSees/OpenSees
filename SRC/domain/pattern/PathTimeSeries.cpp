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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/PathTimeSeries.cpp,v $
                                                                        
                                                                        
// File: ~/domain/pattern/PathTimeSeries.C
//
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for PathTimeSeries.
// PathTimeSeries is a concrete class. A PathTimeSeries object provides
// a linear time series. the factor is given by the pseudoTime and 
// a constant factor provided in the constructor.
//
// What: "@(#) PathTimeSeries.C, revA"


#include <PathTimeSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>

#include <fstream.h>

PathTimeSeries::PathTimeSeries()	
  :TimeSeries(TSERIES_TAG_PathTimeSeries),
   thePath(0), time(0), currentTimeLoc(0), 
   cFactor(0.0), dbTag1(0), dbTag2(0)
{
  // does nothing
}

		   
PathTimeSeries::PathTimeSeries(const Vector &theLoadPath, 
			       const Vector &theTimePath, 
			       double theFactor)
  :TimeSeries(TSERIES_TAG_PathTimeSeries),
   thePath(0), time(0), currentTimeLoc(0), 
   cFactor(theFactor), dbTag1(0), dbTag2(0)

{
  // check vectors are of same size
  if (theLoadPath.Size() != theTimePath.Size()) {
    cerr << "WARNING PathTimeSeries::PathTimeSeries() - vector containing data ";
    cerr << "points for path and time are not of the same size\n";
  } else {

    // create copies of the vectors
    thePath = new Vector(theLoadPath);
    time = new Vector(theTimePath);

    // ensure did not run out of memory creating copies
    if (thePath == 0 || thePath->Size() == 0 ||
	time == 0 || time->Size() == 0) {
      
      cerr << "WARNING PathTimeSeries::PathTimeSeries() - out of memory\n ";
      if (thePath != 0)
	delete thePath;
      if (time != 0)
	delete time;
      thePath = 0;
      time = 0;
    }
  }
}


PathTimeSeries::PathTimeSeries(char *filePathName, 
			       char *fileTimeName, 
			       double theFactor)
  :TimeSeries(TSERIES_TAG_PathTimeSeries),
   thePath(0), time(0), currentTimeLoc(0), 
   cFactor(theFactor), dbTag1(0), dbTag2(0)
{

  // determine the number of data points
  int numDataPoints1 =0;
  int numDataPoints2 =0;
  double dataPoint;
  ifstream theFile;
  
  // first open and go through file containg path
  theFile.open(filePathName, ios::in);
  if (!theFile) {
    cerr << "WARNING - PathTimeSeries::PathTimeSeries()";
    cerr << " - could not open file " << filePathName << endl;
  } else {
    while (theFile >> dataPoint)
      numDataPoints1++;
  }   
  theFile.close();

  // now open and go through file containg time
  theFile.open(fileTimeName, ios::in);
  if (!theFile) {
    cerr << "WARNING - PathTimeSeries::PathTimeSeries()";
    cerr << " - could not open file " << fileTimeName << endl;
  } else {
    while (theFile >> dataPoint)
      numDataPoints2++;
  }   
  theFile.close();

  // check number of data entries in both are the same
  if (numDataPoints1 != numDataPoints2) {
    cerr << "WARNING PathTimeSeries::PathTimeSeries() - files containing data ";
    cerr << "points for path and time do not contain same number of points\n";
  } else {


    // create a vector and read in the data
    if (numDataPoints1 != 0) {

      // now create the two vector
      thePath = new Vector(numDataPoints1);
      time = new Vector(numDataPoints1);

      // ensure did not run out of memory creating copies
      if (thePath == 0 || thePath->Size() == 0 ||
	  time == 0 || time->Size() == 0) {
	  
        cerr << "WARNING PathTimeSeries::PathTimeSeries() - out of memory\n ";
	if (thePath != 0)
	  delete thePath;
	if (time != 0)
	  delete time;
	thePath = 0;
	time = 0;
      }
      
      // first open the path file and read in the data
      theFile.open(filePathName, ios::in);
      if (!theFile) {
	cerr << "WARNING - PathTimeSeries::PathTimeSeries()";
	cerr << " - could not open file " << filePathName << endl;
	delete thePath;
	delete time;
	thePath = 0;
	time =0;
      } else { // read in the path data and then do the time
	int count = 0;
	while (theFile >> dataPoint) {
	  (*thePath)(count) = dataPoint;
	  count++;
	}

	// finally close the file
	theFile.close();

	// now open the time file and read in the data
	theFile.open(fileTimeName, ios::in);
	if (!theFile) {
	  cerr << "WARNING - PathTimeSeries::PathTimeSeries()";
	  cerr << " - could not open file " << fileTimeName << endl;
	  delete thePath;
	  delete time;
	  thePath = 0;
	  time =0;
	} else { // read in the data
	  int count = 0;
	  while (theFile >> dataPoint) {
	    (*time)(count) = dataPoint;
	    count++;
	  }
	} // read in the data 
      }   // read in the path data and then do the time
    }
  }
}




PathTimeSeries::~PathTimeSeries()
{
  if (thePath != 0)
    delete thePath;
  if (time != 0)
    delete time;
}

double
PathTimeSeries::getTimeIncr (double pseudoTime)
{
  // NEED TO FILL IN, FOR NOW return 1.0
  return 1.0;
}

double
PathTimeSeries::getFactor(double pseudoTime)
{
  // check for a quick return
  if (thePath == 0)
    return 0.0;

  // determine indexes into the data array whose boundary holds the time
  double time1 = (*time)(currentTimeLoc);

  // check for another quick return
  if (pseudoTime == time1)
    return cFactor * (*thePath)[currentTimeLoc];

  int size = time->Size();
  int sizem1 = size - 1;
  int sizem2 = size - 2;
  
  // check we are not at the end
  if (pseudoTime > time1 && currentTimeLoc == sizem1)
    return 0.0;

  if (pseudoTime < time1 && currentTimeLoc == 0)
    return 0.0;

  // otherwise go find the current interval
  double time2 = (*time)(currentTimeLoc+1);
  if (pseudoTime > time2) {
    while ((pseudoTime > time2) && (currentTimeLoc < sizem2)) {
      currentTimeLoc++;
      time1 = time2;
      time2 = (*time)(currentTimeLoc+1);
    }
    // if pseudo time greater than ending time reurn 0
    if (pseudoTime > time2)
      return 0.0;

  } else if (pseudoTime < time1) {
    while ((pseudoTime < time1) && (currentTimeLoc > 0)) {
      currentTimeLoc--;
      time2 = time1;	
      time1 = (*time)(currentTimeLoc);
    }
    // if starting time less than initial starting time return 0
    if (pseudoTime < time1)
      return 0.0;
  }

  double value1 = (*thePath)[currentTimeLoc];
  double value2 = (*thePath)[currentTimeLoc+1];
  return cFactor*(value1 + (value2-value1)*(pseudoTime-time1)/(time2 - time1));
}

double
PathTimeSeries::getDuration()
{
  if (thePath == 0)
  {
    cerr << "WARNING -- PathTimeSeries::getDuration() on empty Vector" << endl;
	return 0.0;
  }

  int lastIndex = time->Size(); // index to last entry in time vector
  return ((*time)[lastIndex-1]);
}

double
PathTimeSeries::getPeakFactor()
{
  if (thePath == 0)
  {
    cerr << "WARNING -- PathTimeSeries::getPeakFactor() on empty Vector" << endl;
	return 0.0;
  }

  double peak = fabs((*thePath)[0]);
  int num = thePath->Size();
  double temp;

  for (int i = 1; i < num; i++)
  {
	temp = fabs((*thePath)[i]);
	if (temp > peak)
	  peak = temp;
  }
  
  return (peak*cFactor);
}

int
PathTimeSeries::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  Vector data(4);
  data(0) = cFactor;
  data(1) = -1;
  
  if (thePath != 0) {
    int size = thePath->Size();
    data(1) = size;
    if (dbTag1 == 0) {
      dbTag1 = theChannel.getDbTag();
      dbTag2 = theChannel.getDbTag();
    }
    data(2) = dbTag1;
    data(3) = dbTag2;
  }
  
  int result = theChannel.sendVector(dbTag,commitTag, data);
  if (result < 0) {
    cerr << "PathTimeSeries::sendSelf() - channel failed to send data\n";
    return result;
  }
  
  if (thePath != 0) {
    result = theChannel.sendVector(dbTag1, commitTag, *thePath);
    if (result < 0) {
      cerr << "PathTimeSeries::sendSelf() - ";
      cerr << "channel failed to send tha Path Vector\n";
      return result;  
    }
  }
  if (thePath != 0) {
    result = theChannel.sendVector(dbTag2, commitTag, *time);
    if (result < 0) {
      cerr << "PathTimeSeries::sendSelf() - ";
      cerr << "channel failed to send tha Path Vector\n";
      return result;  
    }
  }
  return 0;
}


int 
PathTimeSeries::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector data(4);
  int result = theChannel.recvVector(dbTag,commitTag, data);
  if (result < 0) {
    cerr << "PathTimeSeries::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    return result;
  }
  cFactor = data(0);
  int size = data(1);


  // get the other path vector
  if (size > 0) {
    dbTag1 = data(2);
    dbTag2 = data(3);
    thePath = new Vector(size);
    time = new Vector(size);
    if (thePath == 0 || time == 0 ||
	thePath->Size() == 0 || time->Size() == 0) {

      cerr << "PathTimeSeries::recvSelf() - ran out of memory";
      cerr << " a Vector of size: " <<  size << endl;  
      if (thePath != 0)
	delete thePath;
      if (time != 0)
	delete time;
      thePath = 0;
      time = 0;
      return -1;
    }
    result = theChannel.recvVector(dbTag1,commitTag, *thePath);    
    if (result < 0) {
      cerr << "PathTimeSeries::recvSelf() - ";
      cerr << "channel failed to receive tha Path Vector\n";
      return result;  
    }
    result = theChannel.recvVector(dbTag2,commitTag, *time);    
    if (result < 0) {
      cerr << "PathTimeSeries::recvSelf() - ";
      cerr << "channel failed to receive tha time Vector\n";
      return result;  
    }
  }
  return 0;    
}



void
PathTimeSeries::Print(ostream &s, int flag)
{
    s << "Path Time Series: constant factor: " << cFactor;
    if (flag == 1 && thePath != 0) {
      s << " specified path: " << *thePath;
      s << " specified time: " << *time;
    }
}
