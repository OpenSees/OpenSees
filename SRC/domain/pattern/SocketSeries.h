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

#ifndef SocketSeries_h
#define SocketSeries_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 03/23
// Revision: A
//
// Description: This file contains the class definition for SocketSeries.
// SocketSeries is a concrete class. A SocketSeries object provides
// a time series by receiving factors through a socket.

#include <TimeSeries.h>
#include <Vector.h>

class Channel;


class SocketSeries : public TimeSeries
{
public:
    // constructors
    SocketSeries(int tag,
        unsigned int port,
        bool udp = false,
        bool checkEndianness = false,
        double socketTimeIncr = 1.0,
        double cFactor = 1.0);
    SocketSeries();
    
    // destructor
    ~SocketSeries();
    
    TimeSeries* getCopy();
    
    // method to get load factor
    double getFactor(double pseudoTime);
    
    // none of the following functions should be invoked on this type of object
    double getDuration() { return 0.0; } // dummy function
    double getPeakFactor() { return cFactor; } // dummy function
    double getTimeIncr(double pseudoTime) { return socketTimeIncr; }
    
    // methods for output    
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel,
        FEM_ObjectBroker& theBroker);
    
    void Print(OPS_Stream& s, int flag = 0);
    
protected:

private:
    double cFactor;
    double socketTimeIncr;

    int recvSize;
    Vector newFactor;
    Channel* theChannel;
    double previousTime;
    double previousFactor;
};

#endif
