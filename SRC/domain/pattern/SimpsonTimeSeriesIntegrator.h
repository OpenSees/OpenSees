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

// $Revision$
// $Date$
// $URL$

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 05/15
// Revision: A
//
// Description: This file contains the class definition for 
// a SimpsonTimeSeriesIntegrator, which integrates a
// ground motion TimeSeries using the Simpson's rule.

#ifndef SimpsonTimeSeriesIntegrator_h
#define SimpsonTimeSeriesIntegrator_h

#include <TimeSeriesIntegrator.h>

class SimpsonTimeSeriesIntegrator : public TimeSeriesIntegrator
{
public:
    SimpsonTimeSeriesIntegrator();
    
    ~SimpsonTimeSeriesIntegrator();
    
    TimeSeries* integrate(TimeSeries *theTimeSeries, double delta);
    
    // methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
        FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
protected:
    
private:
};

#endif
