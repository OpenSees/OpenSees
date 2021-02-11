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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/LinearSeries.h,v $
                                                                        
                                                                        
#ifndef LinearSeries_h
#define LinearSeries_h

// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for LinearSeries.
// LinearSeries is a concrete class. A LinearSeries object provides
// a linear time series. the factor is given by the pseudoTime and 
// a constant factor provided in the constructor.
//
// What: "@(#) LinearSeries.h, revA"

#include <TimeSeries.h>

class LinearSeries : public TimeSeries
{
  public:
    LinearSeries(int tag =0, double cFactor =1.0);

    ~LinearSeries();

    TimeSeries *getCopy(void);

    // method to get load factor
    double getFactor(double pseudoTime);

    // None of the following functions should be invoked on this type
    // of object
    double getDuration () {return 0.0;} // dummy function
    double getPeakFactor () {return cFactor;} // dummy function
    double getTimeIncr (double pseudoTime) {return 1.0;} // dummy function
    
    // methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
	
  private:
    double cFactor;  // factor = pseudoTime * cFactor
};

#endif
