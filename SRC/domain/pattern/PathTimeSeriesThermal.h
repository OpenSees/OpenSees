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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/PathTimeSeriesThermal.h,v $
                                                                        
                                                                        
#ifndef PathTimeSeriesThermal_h
#define PathTimeSeriesThermal_h

// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the class definition for PathTimeSeriesThermal.
// PathTimeSeriesThermal is a TimeSeries class which linear interpolates the
// load factor using user specified control points provided in a vector object.
// the points in the vector are given at time points specified in another vector.
// object. 
//
// What: "@(#) PathTimeSeriesThermal.h, revA"
//Modified by Liming for multi-column path timeseiries input.

#include <TimeSeries.h>

class Vector;

class PathTimeSeriesThermal : public TimeSeries
{
  public:
    // constructors     
  
  PathTimeSeriesThermal(int tag,
		 const char *fileName,int DataNum =9, bool tempOut =true ,
		 double cfactor = 1.0);
  PathTimeSeriesThermal(int tag,
		 int DataNum =9, bool tempOut =true ,
		 double cfactor = 1.0);

    PathTimeSeriesThermal();    
    
    // destructor    
    ~PathTimeSeriesThermal();
    
    TimeSeries *getCopy(void);
  
   int WriteResults(double currentTime, const Vector& newData);

    // method to get factor
    const Vector& getFactors(double pseudoTime);
	  double getFactor(double pseudoTime) {return 0;};
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
	int numCols;
	int numRows;
    Matrix *thePath;      // vector containing the data points
	Vector *CurrentFactors;
    Vector *time;		  // vector containing the time values of data points
    int currentTimeLoc;   // current location in time
    double cFactor;       // additional factor on the returned load factor
    int dbTag1, dbTag2;   // additional database tags needed for vector objects
    int lastSendCommitTag;
	int TempOut;

    Channel *lastChannel;
};

#endif

