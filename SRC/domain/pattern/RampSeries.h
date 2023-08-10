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

// Written: Codi McKee (Texas A&M University)
// Created: 07/2023

#ifndef RampSeries_h
#define RampSeries_h

#include <TimeSeries.h>

class RampSeries : public TimeSeries
{
public:
    // constructors
    RampSeries(int tag, double tStart, double tRamp, double offsetFact, double smoothFact, double cFactor);

    RampSeries();

    // destructor
    ~RampSeries();

    TimeSeries *getCopy(); 

    // method to get load factor
    double getFactor(double pseudoTime);
    double getDuration () {return tStart+tRamp;}
    double getPeakFactor () {return cFactor;}
    double getTimeIncr (double pseudoTime) {return tStart + tRamp;}
    double getStartTime() { return tStart; }
    // methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

protected:

private:
    double tStart{};      // start time of time series (sec)
    double tRamp{};       // Duration of ramp (sec)
    double smoothFact{};  // Smoothing factor
    double offsetFact{};  // Offset amount
    double cFactor{};     // amplitude of sig series
};

#endif
