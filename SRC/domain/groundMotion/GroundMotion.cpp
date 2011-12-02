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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/GroundMotion.cpp,v $
                                                                        
                                                                        
// File: ~/earthquake/GroundMotion.C
//
// Written: fmk 
// Created: 05/98
// Revision: A

#include <iostream.h>

#include <GroundMotion.h>
#include <TrapezoidalGroundMotionIntegrator.h>

//  GroundMotionModel(Domain &theDomain);
//	constructor
GroundMotion::GroundMotion(int theClassTag)
:MovableObject(theClassTag)
{
   theAccelSeries = 0;
   theVelSeries = 0;
   theDispSeries = 0;

   theIntegrator = 0;
}

GroundMotion::~GroundMotion()
{
    
}

void
GroundMotion::setIntegrator (GroundMotionIntegrator *integrator)
{
   theIntegrator = integrator;
}

TimeSeries*
GroundMotion::integrateAccel (double delta)
{
   if (theIntegrator == 0)
   {
	  cerr << "WARNING -- No GroundMotionIntegrator specified" << endl;
	  return 0;
   }

   theIntegrator->setSeries (theAccelSeries);
   return theIntegrator->integrate();
}

TimeSeries*
GroundMotion::integrateVel (double delta)
{
   if (theIntegrator == 0)
   {
	  cerr << "WARNING -- No GroundMotionIntegrator specified" << endl;
      return 0;
   }

   theIntegrator->setSeries (theVelSeries);
   return theIntegrator->integrate();
}
