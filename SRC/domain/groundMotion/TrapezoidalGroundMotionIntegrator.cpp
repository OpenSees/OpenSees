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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/TrapezoidalGroundMotionIntegrator.cpp,v $
                                                                        
                                                                        
// File: ~/domain/groundMotion/TrapezoidalGroundMotionIntegrator.cpp
// 
// Written: MHS
// Created: 10/99
// Revision: A
//
// Description: This file contains the class definition for 
// a TrapezoidalGroundMotionIntegrator, which integrates a
// ground motion TimeSeries using the trapezoidal rule.
//
// What: "@(#) TrapezoidalGroundMotionIntegrator.cpp, revA"

#include <iostream.h>

#include <TrapezoidalGroundMotionIntegrator.h>
#include <Vector.h>
#include <Channel.h>
#include <PathSeries.h>

TrapezoidalGroundMotionIntegrator::TrapezoidalGroundMotionIntegrator
(double del) :
GroundMotionIntegrator(GROUND_MOTION_INTEGRATOR_TAG_Trapezoidal, del)
{

}

TrapezoidalGroundMotionIntegrator::~TrapezoidalGroundMotionIntegrator ()
{

}

TimeSeries*
TrapezoidalGroundMotionIntegrator::integrate ()
{	
   // Check for zero time step, before dividing to get number of steps
   if (0.0 == delta)
   {
	  cerr << "WARNING -- Attempting to integrate with zero time step\n"
		   << "TrapezoidalGroundMotionIntegrator::integrate()" << endl;
	  return 0;
   }

   // Add one to get ceiling out of type cast
   int numSteps = (int)(theSeries->getDuration()/delta + 1.0);

   Vector *theIntegratedValues = new Vector (numSteps);

   // Check that the Vector was allocated properly
   if (theIntegratedValues == 0 || theIntegratedValues->Size() == 0)
   {
	  cerr << "WARNING -- Ran out of memory allocating Vector\n"
		   << "TrapezoidalGroundMotionIntegrator::integrate()" << endl;

	  if (theIntegratedValues != 0)
	     delete theIntegratedValues;

	  theIntegratedValues = 0;
   }
   else
   {
	  int i;                // Counter for indexing
      double dummyTime;     // Dummy variable for integrating
	  double previousValue; // Temporary storage to avoid accessing same value twice
	                        // through identical method calls
	  double currentValue;
      
      // Set the first point
      // Assuming initial condition is zero, i.e. F(0) = 0
	  (*theIntegratedValues)[0] = theSeries->getFactor(0.0) * delta * 0.5;

	  previousValue = (*theIntegratedValues)[0];

	  dummyTime = delta;

	  for (i = 1; i < numSteps; i++, dummyTime += delta)
	  {
		 currentValue = theSeries->getFactor(dummyTime);
	     
		 // Apply the trapezoidal rule to update the integrated value
		 (*theIntegratedValues)[i] = (*theIntegratedValues)[i-1] +
			 delta*0.5 * (currentValue + previousValue);
		 
		 previousValue = currentValue;
	  }

	  // Set the last point
	  (*theIntegratedValues)[i] = (*theIntegratedValues)[i-1] +
		   delta*0.5 * (theSeries->getFactor(dummyTime));
   }

   // Set the method return value
   PathSeries *returnSeries = new PathSeries (*theIntegratedValues, delta);

   if (returnSeries == 0)
   {
	  cerr << "WARNING -- Ran out of memory creating PathSeries\n"
		   << "TrapezoidalGroundMotionIntegrator::integrate()" << endl;

	  returnSeries = 0;
   }

   return returnSeries;
}

int
TrapezoidalGroundMotionIntegrator::sendSelf(int commitTag, Channel &theChannel)
{
   // Need to implement, return dummy value for now
   return 1;
}

int
TrapezoidalGroundMotionIntegrator::recvSelf(int commitTag, Channel &theChannel, 
		           FEM_ObjectBroker &theBroker)
{
   // Need to implement, return dummy value for now
   return 1;
}

void
TrapezoidalGroundMotionIntegrator::Print(ostream &s, int flag)
{
   // Need to implement, return for now
   return;
}
