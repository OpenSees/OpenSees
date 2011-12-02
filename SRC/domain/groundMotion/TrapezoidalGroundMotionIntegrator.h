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
// $Date: 2000-09-15 08:23:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/TrapezoidalGroundMotionIntegrator.h,v $
                                                                        
                                                                        
// File: src/domain/TrapezoidalGroundMotionIntegrator.h
// 
// Written: MHS
// Created: 10/99
// Revision: A
//
// Description: This file contains the class definition for 
// a TrapezoidalGroundMotionIntegrator, which integrates a
// ground motion TimeSeries using the trapezoidal rule.
//
// What: "@(#) TrapezoidalGroundMotionIntegrator.h, revA"

#ifndef TrapezoidalGroundMotionIntegrator_h
#define TrapezoidalGroundMotionIntegrator_h

#include <GroundMotionIntegrator.h>

class TrapezoidalGroundMotionIntegrator : public GroundMotionIntegrator
{
   public:
      TrapezoidalGroundMotionIntegrator (double del = 0.01);

      ~TrapezoidalGroundMotionIntegrator ();

      TimeSeries* integrate ();

	  // methods for output    
      int sendSelf(int commitTag, Channel &theChannel);
      int recvSelf(int commitTag, Channel &theChannel, 
		           FEM_ObjectBroker &theBroker);

      void Print(ostream &s, int flag =0);

   protected:

   private:

};

#endif
