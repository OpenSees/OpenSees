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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/GroundMotionIntegrator.h,v $
                                                                        
                                                                        
/* ****************************************************************** **
**              G3 - Framework for Seismic Simulation                 **
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
**                                                                    **
** Version: 0                                                         **
** Tag: RED                                                           **
** ****************************************************************** */
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
// File: ~/earthquake/GroundMotionIntegrator.h
// 
// Written: MHS
// Created: 10/99
// Revision: A
//
// Description: This file contains the class definition for GroundMotionIntegrator.
// GroundMotionIntegrator is an abstract base class, i.e. no objects of it's
// type can be created. A GroundMotionIntegrator describes the numerical integration
// of a ground motion record used in a dynamic analysis.
//
// What: "@(#) GroundMotion.h, revA"

#ifndef GroundMotionIntegrator_h
#define GroundMotionIntegrator_h

#include <MovableObject.h>
#include <TimeSeries.h>

class GroundMotion;

class GroundMotionIntegrator : public MovableObject
{
  public:
    GroundMotionIntegrator(int classTag, double del = 0.01);

    virtual ~GroundMotionIntegrator();

	virtual TimeSeries* integrate () = 0;

	void setSeries (TimeSeries *series) {theSeries = series;}

  protected:
    TimeSeries *theSeries;
	double delta;

  private:
};

#endif

