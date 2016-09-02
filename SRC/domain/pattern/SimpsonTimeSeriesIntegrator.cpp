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

// Written: Andreas Schellenberg
// Created: 05/15
// Revision: A
//
// Description: This file contains the class definition for 
// a SimpsonTimeSeriesIntegrator, which integrates a
// ground motion TimeSeries using the trapezoidal rule.


#include <SimpsonTimeSeriesIntegrator.h>
#include <Vector.h>
#include <Channel.h>
#include <PathSeries.h>

void* OPS_SimpsonTimeSeriesIntegrator()
{
    return new SimpsonTimeSeriesIntegrator();
}

SimpsonTimeSeriesIntegrator::SimpsonTimeSeriesIntegrator()
    : TimeSeriesIntegrator(TIMESERIES_INTEGRATOR_TAG_Simpson)
{
    
}


SimpsonTimeSeriesIntegrator::~SimpsonTimeSeriesIntegrator()
{
    
}


TimeSeries* SimpsonTimeSeriesIntegrator::integrate(TimeSeries *theSeries, double delta)
{
    // check for zero time step, before dividing to get number of steps
    if (delta <= 0.0)  {
        opserr << "SimpsonTimeSeriesIntegrator::integrate() - attempting to integrate time step "
            << delta << "<= 0.0.\n";
        return 0;
    }
    
    // check a TimeSeries object was passed
    if (theSeries == 0)  {
        opserr << "SimpsonTimeSeriesIntegrator::integrate() - no TimeSeries passed.\n";
        return 0;
    }
    
    // add one to get ceiling out of type cast
    int numSteps = (int)(theSeries->getDuration()/delta + 1.0);
    
    // create new vector for integrated values
    Vector *theInt = new Vector(numSteps);
    
    // check that the Vector was allocated properly
    if (theInt == 0 || theInt->Size() == 0)  {
        opserr << "SimpsonTimeSeriesIntegrator::integrate() - ran out of memory allocating Vector of size "
            << numSteps << endln;
        
        if (theInt != 0)
            delete theInt;
        
        return 0;
    }
    
    double t = 0.0;
    double fi, fj, fk;
    
    // set the first two integrated values (assume that f(0) = 0)
    fi = theSeries->getFactor(0.0);
    fj = theSeries->getFactor(delta);
    fk = theSeries->getFactor(2.0*delta);
    (*theInt)[0] = 0.0;
    (*theInt)[1] = delta/12.0*(5.0*fi + 8.0*fj - fk);
    
    // calculate remaining integrated values
    for (int i=2; i<numSteps-1; i++)  {
        
        (*theInt)[i] = (*theInt)[i-2] + delta/3.0*(fi + 4.0*fj + fk);
        
        // update function values
        fi = fj;
        fj = fk;
        fk = theSeries->getFactor((i+1)*delta);
    }
    
    // calculate the last integrated value
    (*theInt)[numSteps-1] = (*theInt)[numSteps-3] + delta/3.0*(fi + 4.0*fj + fk);
    
    // set the method return value
    PathSeries *returnSeries = new PathSeries(0, *theInt, delta, true);
    
    if (returnSeries == 0)  {
        opserr << "SimpsonTimeSeriesIntegrator::integrate() - ran out of memory creating PathSeries.\n";
        return 0;
    }
    
    return returnSeries;
}


int SimpsonTimeSeriesIntegrator::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}


int SimpsonTimeSeriesIntegrator::recvSelf(int commitTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    return 0;
}


void SimpsonTimeSeriesIntegrator::Print(OPS_Stream &s, int flag)
{
    return;
}
