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
// $Date: 2010-02-04 00:34:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/PulseSeries.h,v $

#ifndef PulseSeries_h
#define PulseSeries_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/04
// Revision: A
//
// Purpose: This file contains the class definition for PulseSeries.
// PulseSeries is a concrete class. A PulseSeries object provides
// a pulse time series. The factor is given by the pseudoTime (t),
// pulse period (T), pulse width (pw) and phase shift (phi),
// and a constant factor provided in the constructor,
// the duration by tStart and tFinal;
//
// What: "@(#) PulseSeries.h, revA"

#include <TimeSeries.h>

class PulseSeries : public TimeSeries
{
public:
    // constructors
    PulseSeries(int tag,
        double tStart,
        double tFinish,
        double period,
        double pWidth,
        double phaseShift,
        double cFactor = 1.0,
        double zeroShift = 0.0);

    PulseSeries();

    // destructor
    ~PulseSeries();

    TimeSeries *getCopy();

    // method to get load factor
    double getFactor(double pseudoTime);
    double getDuration () {return tFinish-tStart;}
    double getPeakFactor () {return cFactor;}
    double getTimeIncr (double pseudoTime) {return tFinish-tStart;}

    // methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

protected:

private:
    double tStart;      // start time of time series (sec)
    double tFinish;     // end time of time series (sec)
    double period;      // period of pulse series (sec)
    double pWidth;      // width of pulse series [0..1] 
    double phaseShift;  // phase shift of pulse series (sec)
    double cFactor;     // amplitude of pulse series
    double zeroShift;   // zero shift of pulse series
};

#endif
