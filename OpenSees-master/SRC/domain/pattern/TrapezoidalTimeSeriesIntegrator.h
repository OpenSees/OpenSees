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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:01:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TrapezoidalTimeSeriesIntegrator.h,v $
                                                                        
                                                                        
// File: src/domain/TrapezoidalTimeSeriesIntegrator.h
// 
// Written: MHS
// Created: 10/99
// Revision: A
//
// Description: This file contains the class definition for 
// a TrapezoidalTimeSeriesIntegrator, which integrates a
// ground motion TimeSeries using the trapezoidal rule.
//
// What: "@(#) TrapezoidalTimeSeriesIntegrator.h, revA"

#ifndef TrapezoidalTimeSeriesIntegrator_h
#define TrapezoidalTimeSeriesIntegrator_h

#include <TimeSeriesIntegrator.h>

class TrapezoidalTimeSeriesIntegrator : public TimeSeriesIntegrator
{
   public:
      TrapezoidalTimeSeriesIntegrator();

      ~TrapezoidalTimeSeriesIntegrator();

      TimeSeries* integrate(TimeSeries *theTimeSeries, double delta);

      // methods for output    
      int sendSelf(int commitTag, Channel &theChannel);
      int recvSelf(int commitTag, Channel &theChannel, 
		   FEM_ObjectBroker &theBroker);

      void Print(OPS_Stream &s, int flag =0);

   protected:

   private:
};

#endif
