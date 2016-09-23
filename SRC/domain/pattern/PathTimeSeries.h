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
                                                                        
// $Revision: 1.7 $
// $Date: 2010-02-04 00:34:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/PathTimeSeries.h,v $
                                                                        
                                                                        
#ifndef PathTimeSeries_h
#define PathTimeSeries_h

// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the class definition for PathTimeSeries.
// PathTimeSeries is a TimeSeries class which linear interpolates the
// load factor using user specified control points provided in a vector object.
// the points in the vector are given at time points specified in another vector.
// object. 
//
// What: "@(#) PathTimeSeries.h, revA"

#include <TimeSeries.h>

class Vector;

class PathTimeSeries : public TimeSeries
{
  public:
    // constructors  
  PathTimeSeries(int tag,
		 const Vector &thePath, 
		 const Vector &theTime, 
		 double cfactor = 1.0,
         bool useLast = false);
  
  PathTimeSeries(int tag,
		 const char *fileNamePath, 
		 const char *fileNameTime, 
		 double cfactor = 1.0,
         bool useLast = false);
  
  PathTimeSeries(int tag,
		 const char *fileName,
		 double cfactor = 1.0,
         bool useLast = false);

    PathTimeSeries();    
    
    // destructor    
    ~PathTimeSeries();
    
    TimeSeries *getCopy(void);

    // method to get factor
    double getFactor(double pseudoTime);
    double getDuration ();
    double getPeakFactor ();
    double getTimeIncr (double pseudoTime);

    // methods for output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
    
  private:
    Vector *thePath;      // vector containg the data points
    Vector *time;		  // vector containg the time values of data points
    int currentTimeLoc;   // current location in time
    double cFactor;       // additional factor on the returned load factor
    int dbTag1, dbTag2;   // additional database tags needed for vector objects
    int lastSendCommitTag;
    Channel *lastChannel;
    bool useLast;
};

#endif

