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
                                                                        
// $Revision: 1.3 $
// $Date: 2010-02-04 00:34:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/RectangularSeries.h,v $
                                                                        
                                                                        
#ifndef RectangularSeries_h
#define RectangularSeries_h

// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for RectangularSeries.
// RectangularSeries is a concrete class. A RectangularSeries object provides
// a rectangular pulse time series. the factor is given by the pseudoTime and 
// a constant factor provided in the constructor, the duration by tStart
// and tFinal;
//
// What: "@(#) RectangularSeries.h, revA"

#include <TimeSeries.h>

class RectangularSeries : public TimeSeries
{
  public:
    // constructors
    RectangularSeries(int tag, double tStart, double tFinal, double cFactor = 1.0);
    RectangularSeries();
    
    // destructor
    ~RectangularSeries();

    TimeSeries *getCopy(void);
    
    // method to get load factor
    double getFactor(double pseudoTime);
    double getDuration () {return (tFinish - tStart);}
    double getPeakFactor () {return cFactor;}
    double getTimeIncr (double pseudoTime) {return (tFinish - tStart);}
    
    // methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
	
  private:
    double tStart;    // start time of the pulse
	double tFinish;   // finish time of the pulse
    double cFactor;   // factor = cFactor for tStart <= t <= tFinish
};

#endif
