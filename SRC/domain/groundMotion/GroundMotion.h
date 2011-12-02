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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/GroundMotion.h,v $
                                                                        
                                                                        
// File: ~/earthquake/GroundMotion.h
// 
// Written: fmk 
// Created: 05/98
// Revision: A
//
// Description: This file contains the class definition for GroundMotion.
// GroundMotion is an abstract base class, i.e. no objects of it's
// type can be created. A GroundMotion describes the motion of the ground
// in a dynamic analysis.
//
// What: "@(#) GroundMotion.h, revA"

#ifndef GroundMotion_h
#define GroundMotion_h

#include <MovableObject.h>
#include <TimeSeries.h>
#include <GroundMotionIntegrator.h>


class GroundMotion : public MovableObject
{
  public:
    GroundMotion(int classTag);
    virtual ~GroundMotion();

    virtual double getDuration (void) = 0;

    virtual double getPeakAccel () = 0;
    virtual double getPeakVel () = 0;
    virtual double getPeakDisp () = 0;

    virtual double getAccel (double time) = 0;
    virtual double getVel (double time) = 0;
    virtual double getDisp (double time) = 0;

    void setIntegrator (GroundMotionIntegrator *integrator);
    TimeSeries *integrateAccel (double delta = 0.01); 
    TimeSeries *integrateVel (double delta = 0.01);
    
  protected:
    TimeSeries *theAccelSeries;  // Ground acceleration
    TimeSeries *theVelSeries;	 // Ground velocity
    TimeSeries *theDispSeries;	 // Ground displacement

    GroundMotionIntegrator *theIntegrator;

  private:

};

#endif

