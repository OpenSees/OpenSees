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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/PathSeries.h,v $
                                                                        
                                                                        
#ifndef PathSeries_h
#define PathSeries_h

// File: ~/analysis/integrator/PathSeries.h
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the class definition for PathSeries.
// PathSeries is a TimeSeries class which linear interpolates the
// load factor using user specified control points provided in a vector object.
// the points in the vector are given at regular time increments pathTimeIncr
// apart. (could be provided in another vector if different)
//
// What: "@(#) PathSeries.h, revA"

#include <TimeSeries.h>
#include <iostream.h>

class Vector;

class PathSeries : public TimeSeries
{
  public:
    // constructors  
    PathSeries(const Vector &thePath, 
	       double pathTimeIncr = 1.0, 
	       double cfactor = 1.0);
    PathSeries(char *fileName, 
	       double pathTimeIncr = 1.0, 
	       double cfactor = 1.0);    
    PathSeries();    
    
    // destructor    
    ~PathSeries();

    // method to get factor
    double getFactor(double pseudoTime);
	double getDuration ();
	double getPeakFactor ();
	double getTimeIncr (double pseudoTime) {return pathTimeIncr;}

    // methods for output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    void Print(ostream &s, int flag =0);    
    
  protected:
    
  private:
    Vector *thePath;      // vector containg the data points
    double pathTimeIncr;  // specifies the time increment used in load path vector
    double cFactor;       // additional factor on the returned load factor
    int otherDbTag;       // a database tag needed for the vector object
};

#endif

