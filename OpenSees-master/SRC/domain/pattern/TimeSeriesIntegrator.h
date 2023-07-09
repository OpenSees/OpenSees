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
// $Date: 2003-03-04 00:48:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TimeSeriesIntegrator.h,v $
                                                                        
// Written: MHS
// Created: 10/99
// Revision: A
//
// Description: This file contains the class definition for TimeSeriesIntegrator.
// TimeSeriesIntegrator is an abstract base class, i.e. no objects of it's
// type can be created. A TimeSeriesIntegrator describes the numerical integration
// of a ground motion record used in a dynamic analysis.
//
// What: "@(#) TimeSeries.h, revA"

#ifndef TimeSeriesIntegrator_h
#define TimeSeriesIntegrator_h

#include <MovableObject.h>
#include <TimeSeries.h>

class TimeSeries;

class TimeSeriesIntegrator : public MovableObject
{
  public:
    TimeSeriesIntegrator(int classTag);

    virtual ~TimeSeriesIntegrator();

    virtual TimeSeries* integrate(TimeSeries *theSeries, double delta) = 0;

  protected:

  private:
};

#endif

