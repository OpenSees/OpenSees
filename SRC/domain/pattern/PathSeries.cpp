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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/PathSeries.cpp,v $
                                                                        
                                                                        
// File: ~/domain/pattern/PathSeries.C
//
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for PathSeries.
// PathSeries is a concrete class. A PathSeries object provides
// a linear time series. the factor is given by the pseudoTime and 
// a constant factor provided in the constructor.
//
// What: "@(#) PathSeries.C, revA"


#include <PathSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>

#include <fstream.h>
#include <iostream.h>

PathSeries::PathSeries()	
  :TimeSeries(TSERIES_TAG_PathSeries),
   thePath(0), pathTimeIncr(0.0), cFactor(0.0), otherDbTag(0)
{
  // does nothing
}

		   
PathSeries::PathSeries(const Vector &theLoadPath, 
		       double theTimeIncr, 
		       double theFactor)
  :TimeSeries(TSERIES_TAG_PathSeries),
   thePath(0), pathTimeIncr(theTimeIncr), cFactor(theFactor), otherDbTag(0)
{
  // create a copy of the vector containg path points
  thePath = new Vector(theLoadPath);

  // ensure we did not run out of memory
  if (thePath == 0 || thePath->Size() == 0) {
    cerr << "PathSeries::PathSeries() - ran out of memory constructing";
    cerr << " a Vector of size: " <<  theLoadPath.Size() << endl;
    if (thePath != 0)
      delete thePath; 
    thePath = 0;
  }
}


PathSeries::PathSeries(char *fileName, 
		       double theTimeIncr, 
		       double theFactor)
  :TimeSeries(TSERIES_TAG_PathSeries),
   thePath(0), pathTimeIncr(theTimeIncr), cFactor(theFactor)
{
  // determine the number of data points
  int numDataPoints =0;
  double dataPoint;
  ifstream theFile;
  theFile.open(fileName);

  if (!theFile) {
    cerr << "WARNING - PathSeries::PathSeries()";
    cerr << " - could not open file " << fileName << endl;
  } else {
    while (theFile >> dataPoint)
      numDataPoints++;
  }   
  theFile.close();


  // create a vector and read in the data
  if (numDataPoints != 0) {

    // first open the file
    theFile.open(fileName, ios::in);
    if (!theFile) {
      cerr << "WARNING - PathSeries::PathSeries()";
      cerr << " - could not open file " << fileName << endl;
    } else {

      // now create the vector
      thePath = new Vector(numDataPoints);

      // ensure we did not run out of memory
      if (thePath == 0 || thePath->Size() == 0) {
	cerr << "PathSeries::PathSeries() - ran out of memory constructing";
	cerr << " a Vector of size: " << numDataPoints << endl;

	if (thePath != 0)
	  delete thePath;
	thePath = 0;
      }

      // read the data from the file
      else {
	int count = 0;
	while (theFile >> dataPoint) {
	  (*thePath)(count) = dataPoint;
	  count++;
	}
      }

      // finally close the file
      theFile.close();
    }
  }
}


PathSeries::~PathSeries()
{
  if (thePath != 0)
    delete thePath;
}

double
PathSeries::getFactor(double pseudoTime)
{
  // check for a quick return
  if (pseudoTime < 0.0 || thePath == 0)
    return 0.0;

  // determine indexes into the data array whose boundary holds the time
  double incr = pseudoTime/pathTimeIncr; 
  int incr1 = floor(incr);
  int incr2 = incr1+1;

  if (incr2 >= thePath->Size())
    return 0.0;
  else {
    double value1 = (*thePath)[incr1];
    double value2 = (*thePath)[incr2];
    return cFactor*(value1 + (value2-value1)*(pseudoTime/pathTimeIncr - incr1));
  }
}


double
PathSeries::getDuration()
{
  if (thePath == 0)
  {
    cerr << "WARNING -- PathSeries::getDuration() on empty Vector" << endl;
	return 0.0;
  }
  return (thePath->Size() * pathTimeIncr);
}

double
PathSeries::getPeakFactor()
{
  if (thePath == 0)
  {
    cerr << "WARNING -- PathSeries::getPeakFactor() on empty Vector" << endl;
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
PathSeries::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  Vector data(4);
  data(0) = cFactor;
  data(1) = pathTimeIncr;
  data(2) = -1;
  
  if (thePath != 0) {
    int size = thePath->Size();
    data(2) = size;
    if (otherDbTag == 0)
      otherDbTag = theChannel.getDbTag();
    data(3) = otherDbTag;
  }
  
  int result = theChannel.sendVector(dbTag,commitTag, data);
  if (result < 0) {
    cerr << "PathSeries::sendSelf() - channel failed to send data\n";
    return result;
  }
  
  if (thePath != 0) {
    result = theChannel.sendVector(otherDbTag, commitTag, *thePath);
    if (result < 0) {
      cerr << "PathSeries::sendSelf() - ";
      cerr << "channel failed to send tha Path Vector\n";
      return result;  
    }
  }
  return 0;
}


int 
PathSeries::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector data(4);
  int result = theChannel.recvVector(dbTag,commitTag, data);
  if (result < 0) {
    cerr << "PathSeries::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    return result;
  }
  cFactor = data(0);
  pathTimeIncr = data(1);
  int size = data(2);
  otherDbTag = data(3);

  // get the other path vector
  if (size > 0) {
    thePath = new Vector(size);
    if (thePath == 0 || thePath->Size() == 0) {
      cerr << "PathSeries::recvSelf() - ran out of memory";
      cerr << " a Vector of size: " <<  size << endl;  
      if (thePath != 0)
	delete thePath;
      thePath = 0;
      return -1;
    }
    result = theChannel.recvVector(otherDbTag,commitTag, *thePath);    
    if (result < 0) {
      cerr << "PathSeries::recvSelf() - ";
      cerr << "channel failed to receive tha Path Vector\n";
      return result;  
    }
  }
  return 0;    
}



void
PathSeries::Print(ostream &s, int flag)
{
    //s << "Path Time Series: constant factor: " << cFactor;
    //s << "  time Incr: " << pathTimeIncr << endl;
    if (flag == 1 && thePath != 0)
      //s << " specified path: " << *thePath;
	  s << *thePath;
}
