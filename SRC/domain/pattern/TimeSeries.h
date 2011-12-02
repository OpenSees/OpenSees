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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TimeSeries.h,v $
                                                                        
                                                                        
#ifndef TimeSeries_h
#define TimeSeries_h

// File: ~/domain/pattern/TimeSeries.h
//
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for TimeSeries.
// TimeSeries is an abstract class. A TimeSeries object is used to
// determine the load factor to be applied to the loads in a pettern.
// to the model. 
//
// What: "@(#) TimeSeries.h, revA"

#include <MovableObject.h>
#include <iostream.h>

class TimeSeries : public MovableObject
{
  public:
    TimeSeries(int classTag);
    virtual ~TimeSeries();

    // pure virtual functions
    virtual double getFactor (double pseudoTime) = 0;
    virtual double getDuration () = 0;
    virtual double getPeakFactor () = 0;

    virtual double getTimeIncr (double pseudoTime) = 0;
    // This is defined to be the time increment from the argument
    // 'pseudoTime' to the NEXT point in the time series path
    // THIS MAY CHANGE -- MAY BE BETTER TO GET THE TIME INCREMENT
    // FROM THE PREVIOUS POINT IN THE PATH UP TO 'pseudoTime', WILL
    // DECIDE ONCE GroundMotionIntegrator IS IMPLEMENTED

    virtual void Print(ostream &s, int flag = 0) = 0;        

  protected:
	
  private:
};

#endif
