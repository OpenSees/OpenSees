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

// $Revision: 1.3 $
// $Date: 2010-02-04 00:34:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TrigSeries.h,v $


#ifndef TrigSeries_h
#define TrigSeries_h

// Written: MHS 
// Created: July 2000
// Revision: A
//
// Purpose: This file contains the class definition for TrigSeries.
// TrigSeries is a concrete class. A TrigSeries object provides
// a sine time series. the factor is given by the pseudoTime (x),
// pulse period (T), phase shift (phi), i.e. by sin(2*PI*(x-xo)/T + phi),
// and a constant factor provided in the constructor,
// the duration by tStart and tFinal;
//
// What: "@(#) TrigSeries.h, revA"

#include <TimeSeries.h>

class TrigSeries : public TimeSeries
{
public:
    // constructors
    TrigSeries(int tag,
        double tStart, 
        double tFinish,
        double period, 
        double phaseShift, 
        double cFactor = 1.0,
        double zeroShift = 0.0);

    TrigSeries();

    // destructor
    ~TrigSeries();

    TimeSeries *getCopy(); 

    // method to get load factor
    double getFactor(double pseudoTime);
    double getDuration () {return tFinish-tStart;}
    double getPeakFactor () {return cFactor;}
    double getTimeIncr (double pseudoTime) {return tFinish-tStart;}
    double getStartTime() { return tStart; } 
    // methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    double getFactorSensitivity(double pseudoTime);
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int parameterID);
    // AddingSensitivity:BEGIN //////////////////////////////////////////  

protected:

private:
    double tStart;      // start time of time series (sec)
    double tFinish;     // end time of time series (sec)
    double period;      // period of trig series (sec)
    double phaseShift;  // phase shift of trig series (rad)
    double cFactor;     // amplitude of trig series
    double zeroShift;   // zero shift of trig series
  int parameterID;
};

#endif
