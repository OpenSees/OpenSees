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
                                                                        
// $Revision$
// $Date$
// $URL$
                                                                        
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for PathSeries.
// PathSeries is a concrete class. A PathSeries object provides
// a linear time series. the factor is given by the pseudoTime and 
// a constant factor provided in the constructor.

#include <PathSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Parameter.h>

#include <fstream>
using std::ifstream;

#include <iomanip>
using std::ios;

#include <PathTimeSeries.h>
#include <elementAPI.h>
#include <string>

void *OPS_PathSeries() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "insufficient arguments: PathSeries\n";
        return 0;
    }

    // get tag
    int tag = 0;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) < 0) {
        opserr << "tag is not specified\n";
        return 0;
    }

    // get other inputs
    double factor = 1.0, dt = -1.0;
    const char *fileTime = 0, *filePath = 0, *fileName = 0;
    std::vector<double> values, times;
    bool useLast = false, prependZero = false;
    double startTime = 0.0;

    // check inputs
    TimeSeries *theSeries = 0;
    int loc = 2;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        // next arg
        const char *arg = OPS_GetString();
        loc++;

        // check arg
        if (strcmp(arg, "-dt") == 0 ||
            strcmp(arg, "-dT") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: dt is not specified\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &dt) < 0) {
                opserr << "WARNING: failed to get dt\n";
                return 0;
            }
            loc++;

        } else if (strcmp(arg, "-tag") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: tag is not specified\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetIntInput(&numData, &tag) < 0) {
                opserr << "WARNING: failed to get tag\n";
                return 0;
            }
            loc++;

        } else if (strcmp(arg, "-values") == 0) {
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                numData = 1;
                if (OPS_GetDoubleInput(&numData, &val) <
                    0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                values.push_back(val);
                loc++;
            }
        } else if (strcmp(arg, "-factor") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING no factor is given\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &factor) < 0) {
                opserr << "WARNING invalid factor\n";
                return 0;
            }
            loc++;

        } else if (strcmp(arg, "-useLast") == 0) {
            useLast = true;

        } else if (strcmp(arg, "-prependZero") == 0) {
            prependZero = true;

        } else if (strcmp(arg, "-startTime") == 0 ||
                   strcmp(arg, "-tStart") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr
                    << "WARNING no start time is given\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &startTime) <
                0) {
                opserr << "WARNING invalid start time\n";
                return 0;
            }
            loc++;

        } else if (strcmp(arg, "-filePath") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING no file path is given\n";
                return 0;
            }
            filePath = OPS_GetString();
            loc++;

        } else if (strcmp(arg, "-fileTime") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING no file time is given\n";
                return 0;
            }
            fileTime = OPS_GetString();
            loc++;

        } else if (strcmp(arg, "-file") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING no file is given\n";
                return 0;
            }
            fileName = OPS_GetString();
            loc++;

        } else if (strcmp(arg, "-time") == 0) {
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                numData = 1;
                if (OPS_GetDoubleInput(&numData, &val) <
                    0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                times.push_back(val);
                loc++;
            }
        }
    }

    if (dt > 0 && values.empty() == false) {
        Vector thePath(&values[0], (int)values.size());
        return new PathSeries(tag, thePath, dt, factor,
                              useLast, prependZero,
                              startTime);

    } else if (dt > 0 && filePath != 0) {
        return new PathSeries(tag, filePath, dt, factor,
                              useLast, prependZero,
                              startTime);

    } else if (times.empty() == false &&
               values.empty() == false) {
        Vector thePath(&values[0], (int)values.size());
        Vector theTime(&times[0], (int)times.size());
        return new PathTimeSeries(tag, thePath, theTime,
                                  factor);

    } else if (fileTime != 0 && filePath != 0) {
        return new PathTimeSeries(tag, filePath, fileTime,
                                  factor);
    } else if (fileName != 0) {
        return new PathTimeSeries(tag, fileName, factor,
                                  useLast);
    }

    opserr << "WARNING choice of options for path series "
              "is invalid\n";
    return 0;
}

PathSeries::PathSeries()	
  :TimeSeries(TSERIES_TAG_PathSeries),
   thePath(0), pathTimeIncr(0.0), cFactor(0.0), otherDbTag(0), lastSendCommitTag(-1), startTime(0.0), parameterID(0)
{
  // does nothing
}

PathSeries::PathSeries(int tag,
		       const Vector &theLoadPath, 
		       double theTimeIncr, 
		       double theFactor,
		       bool last,
               bool prependZero,
               double tStart)
  :TimeSeries(tag, TSERIES_TAG_PathSeries),
   thePath(0), pathTimeIncr(theTimeIncr), cFactor(theFactor),
   otherDbTag(0), lastSendCommitTag(-1), useLast(last), startTime(tStart), parameterID(0)
{
  // create a copy of the vector containing path points
  if (prependZero == false) {
    thePath = new Vector(theLoadPath);
  } else if (theLoadPath != 0) {
    // prepend a zero value
    thePath = new Vector(1 + theLoadPath.Size());
    thePath->Assemble(theLoadPath, 1);
  }

  // ensure we did not run out of memory
  if (thePath == 0 || thePath->Size() == 0) {
    opserr << "PathSeries::PathSeries() - ran out of memory constructing";
    opserr << " a Vector of size: " <<  theLoadPath.Size() << endln;
    if (thePath != 0)
      delete thePath; 
    thePath = 0;
  }
}

PathSeries::PathSeries(int tag,
		       const char *fileName, 
		       double theTimeIncr, 
		       double theFactor,
		       bool last,
               bool prependZero,
               double tStart)
  :TimeSeries(tag, TSERIES_TAG_PathSeries),
   thePath(0), pathTimeIncr(theTimeIncr), cFactor(theFactor),
   otherDbTag(0), lastSendCommitTag(-1), useLast(last), startTime(tStart), parameterID(0)
{
  // determine the number of data points .. open file and count num entries
  int numDataPoints = 0;
  double dataPoint;

  ifstream theFile;
  theFile.open(fileName);

  if (theFile.bad() || !theFile.is_open()) {
    opserr << "WARNING - PathSeries::PathSeries()";
    opserr << " - could not open file " << fileName << endln;
  } else {
    while (theFile >> dataPoint)
      numDataPoints++;
  }   
  theFile.close();

  // create a vector and read in the data
  if (numDataPoints != 0) {

    // increment size if we need to prepend a zero value
    if (prependZero == true)
      numDataPoints++;
    
    // first open the file
    ifstream theFile1;
    theFile1.open(fileName, ios::in);
    if (theFile1.bad() || !theFile1.is_open()) {
      opserr << "WARNING - PathSeries::PathSeries()";
      opserr << " - could not open file " << fileName << endln;
    } else {

      // now create the vector
      thePath = new Vector(numDataPoints);

      // ensure we did not run out of memory
      if (thePath == 0 || thePath->Size() == 0) {
	opserr << "PathSeries::PathSeries() - ran out of memory constructing";
	opserr << " a Vector of size: " << numDataPoints << endln;

	if (thePath != 0)
	  delete thePath;
	thePath = 0;
      }

      // read the data from the file
      else {
	int count = 0;
    if (prependZero == true)
      count++;
	while (theFile1 >> dataPoint) {
	  (*thePath)(count) = dataPoint;
	  count++;
	}
      }

      // finally close the file
      theFile1.close();
    }
  }
}

PathSeries::~PathSeries()
{
  if (thePath != 0)
    delete thePath;
}

TimeSeries *
PathSeries::getCopy(void) {
	if (thePath != 0)
  return new PathSeries(this->getTag(), *thePath, pathTimeIncr, cFactor,
                        useLast, false, startTime);
	else
		return 0;
}

double
PathSeries::getFactor(double pseudoTime)
{
  // check for a quick return
  if (pseudoTime < startTime || thePath == 0)
    return 0.0;

  // determine indexes into the data array whose boundary holds the time
  double incr = (pseudoTime-startTime)/pathTimeIncr; 
  long long incr1 = floor(incr);
  long long incr2 = incr1+1;
  int size = thePath->Size();

  if (incr2 >= size) {
    if (useLast == false)
      return 0.0;
    else
      return cFactor*(*thePath)[size-1];
  } else {
    double value1 = (*thePath)[incr1];
    double value2 = (*thePath)[incr2];
    return cFactor*(value1 + (value2-value1)*(incr - incr1));
  }
}

double
PathSeries::getDuration()
{
  if (thePath == 0)
  {
    opserr << "WARNING -- PathSeries::getDuration() on empty Vector" << endln;
	return 0.0;
  }
  return (thePath->Size()*pathTimeIncr);
  //return (startTime + thePath->Size()*pathTimeIncr);
}

double
PathSeries::getStartTime() {

	return startTime;
}


double
PathSeries::getPeakFactor()
{
  if (thePath == 0)
  {
    opserr << "WARNING -- PathSeries::getPeakFactor() on empty Vector" << endln;
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

  Vector data(7);
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

  if ((lastSendCommitTag == -1) && (theChannel.isDatastore() == 1)) {
    lastSendCommitTag = commitTag;
  }

  data(4) = lastSendCommitTag;

  if (useLast == true)
    data(5) = 1;
  else
    data(5) = 0;

  data(6) = startTime;

  int result = theChannel.sendVector(dbTag, commitTag, data);
  if (result < 0) {
    opserr << "PathSeries::sendSelf() - channel failed to send data\n";
    return result;
  }

  // we only send the vector data if this is the first time it is sent to the database
  // or the channel is for sending the data to a remote process

  if ((lastSendCommitTag == commitTag) || (theChannel.isDatastore() == 0)) {
    if (thePath != 0) {
      result = theChannel.sendVector(otherDbTag, commitTag, *thePath);
      if (result < 0) {
	opserr << "PathSeries::sendSelf() - ";
	opserr << "channel failed to send the Path Vector\n";
	return result;  
      }
    }
  }

  return 0;
}

int 
PathSeries::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  Vector data(7);
  int result = theChannel.recvVector(dbTag, commitTag, data);
  if (result < 0) {
    opserr << "PathSeries::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    return result;
  }

  cFactor = data(0);
  pathTimeIncr = data(1);
  int size = data(2);
  otherDbTag = data(3);
  lastSendCommitTag = data(4);

  if (data(5) == 1)
    useLast = true;
  else
    useLast = false;

  startTime = data(6);
  
  // get the path vector, only receive it once as it can't change
  if (thePath == 0 && size > 0) {
    thePath = new Vector(size);
    if (thePath == 0 || thePath->Size() == 0) {
      opserr << "PathSeries::recvSelf() - ran out of memory";
      opserr << " a Vector of size: " <<  size << endln;  
      if (thePath != 0)
	delete thePath;
      thePath = 0;
      return -1;
    }

    result = theChannel.recvVector(otherDbTag, lastSendCommitTag, *thePath);    
    if (result < 0) {
      opserr << "PathSeries::recvSelf() - ";
      opserr << "channel failed to receive the Path Vector\n";
      return result;  
    }
  }

  return 0;    
}

void
PathSeries::Print(OPS_Stream &s, int flag)
{
    //s << "Path Time Series: constant factor: " << cFactor;
    //s << "  time Incr: " << pathTimeIncr << endln;
    if (flag == 1 && thePath != 0)
      //s << " specified path: " << *thePath;
	  s << *thePath;
}

    // AddingSensitivity:BEGIN //////////////////////////////////////////
double
PathSeries::getFactorSensitivity(double pseudoTime)
{
  // check for a quick return
  if (pseudoTime < startTime || thePath == 0)
    return 0.0;

  if (parameterID != 1)
    return 0.0;
  
  // determine indexes into the data array whose boundary holds the time
  double incr = (pseudoTime-startTime)/pathTimeIncr; 
  long long incr1 = floor(incr);
  long long incr2 = incr1+1;
  int size = thePath->Size();

  if (incr2 >= size) {
    if (useLast == false)
      return 0.0;
    else
      return (*thePath)[size-1];
  } else {
    double value1 = (*thePath)[incr1];
    double value2 = (*thePath)[incr2];
    return (value1 + (value2-value1)*(incr - incr1));
  }
  }

int 
PathSeries::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strncmp(argv[0],"factor",80) == 0) {
    param.setValue(cFactor);
    return param.addObject(1, this);
  }

  return -1;
}
   
int 
PathSeries::updateParameter(int parameterID, Information &info)
{
  if (parameterID == 1) {
    cFactor = info.theDouble;
    return 0;
  }

  return -1;
}

int
PathSeries::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}
