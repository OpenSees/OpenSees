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

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 03/23
// Revision: A
//
// Description: This file contains the class definition for SocketSeries.
// SocketSeries is a concrete class. A SocketSeries object provides
// a time series by receiving factors through a socket.

#include <SocketSeries.h>

#include <Vector.h>
#include <Channel.h>
#include <TCP_Socket.h>
#include <UDP_Socket.h>

#include <elementAPI.h>


void* OPS_SocketSeries()
{
    // pointer to the time series that will be returned
    TimeSeries* theSeries = 0;
    
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING invalid number of arguments\n";
        opserr << "Want: timeSeries Socket tag port "
            << "<-udp> <-dt> <-factor f)\n";
        return 0;
    }
    
    // get tag
    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) != 0) {
        opserr << "WARNING invalid timeSeries Socket tag\n";
        return 0;
    }
    
    // get the port
    int ipPort;
    if (OPS_GetIntInput(&numdata, &ipPort) != 0) {
        opserr << "WARNING invalid ipPort \n";
        opserr << "timeSeries Socket series: " << tag << endln;
        return 0;
    }
    
    // optional parameters
    const char* type;
    bool udp = false;
    bool checkEndianness = false;
    double dt = 1.0;
    double cFactor = 1.0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        type = OPS_GetString();
        // udp socket
        if (strcmp(type, "-udp") == 0 ||
            strcmp(type, "-UDP") == 0) {
            udp = true;
        }
        // dt
        else if (strcmp(type, "-dt") == 0 ||
            strcmp(type, "-dT") == 0) {
            if (OPS_GetDoubleInput(&numdata, &dt) < 0) {
                opserr << "WARNING: invalid -dt value\n";
                opserr << "timeSeries Socket series: " << tag << endln;
                return 0;
            }
        }
        // factor
        else if (strcmp(type, "-factor") == 0) {
            if (OPS_GetDoubleInput(&numdata, &cFactor) < 0) {
                opserr << "WARNING: invalid -cFactor value\n";
                opserr << "timeSeries Socket series: " << tag << endln;
                return 0;
            }
        }
    }
    
    // now create the time series
    theSeries = new SocketSeries(tag, (unsigned int)ipPort, udp,
                                 checkEndianness, dt, cFactor);
    if (theSeries == 0) {
        opserr << "WARNING ran out of memory creating time series\n";
        opserr << "timeSeries Socket series: " << tag << endln;
        return 0;
    }
    
    return theSeries;
}


SocketSeries::SocketSeries(int tag,
    unsigned int port, bool udp, bool checkEndianness,
    double theTimeIncr, double theFactor)
    : TimeSeries(tag, TSERIES_TAG_SocketSeries),
    socketTimeIncr(theTimeIncr), cFactor(theFactor),
    recvSize(1), newFactor(1), theChannel(0),
    previousTime(0.0), previousFactor(0.0)
{
    // setup the connection
    if (udp)
        theChannel = new UDP_Socket(port, checkEndianness);
    else
        theChannel = new TCP_Socket(port, checkEndianness);
    
    if (theChannel != 0) {
        opserr << "\nChannel successfully created: "
            << "Waiting for client to send time series values...\n";
    } else {
        opserr << "SocketSeries::SocketSeries() "
            << "- failed to create channel\n";
        exit(-1);
    }
    if (theChannel->setUpConnection() != 0) {
        opserr << "SocketSeries::SocketSeries() "
            << "- failed to setup connection\n";
        delete theChannel;
        theChannel = 0;
    }
}


SocketSeries::SocketSeries()
    : TimeSeries(TSERIES_TAG_SocketSeries),
    socketTimeIncr(0.0), cFactor(1.0),
    recvSize(1), newFactor(1), theChannel(0),
    previousTime(0.0), previousFactor(0.0)
{
    // does nothing
}


SocketSeries::~SocketSeries()
{
    if (theChannel != 0)
        delete theChannel;
}


TimeSeries* SocketSeries::getCopy()
{
    return new SocketSeries(*this);
}


double SocketSeries::getFactor(double pseudoTime)
{
    double factor;
    
    // get new factor if time has advanced by >= dt
    if (pseudoTime >= previousTime + socketTimeIncr) {
        if (theChannel->recvVector(0, 0, newFactor) < 0) {
            opserr << "SocketSeries::getFactor() "
                << "- failed to receive data\n";
            return -1;
        }
        previousTime += socketTimeIncr;
    }
    // in case time is at zero get the first factor
    else if (pseudoTime == 0.0) {
        if (theChannel->recvVector(0, 0, newFactor) < 0) {
            opserr << "SocketSeries::getFactor() "
                << "- failed to receive data\n";
            return -1;
        }
    }
    
    factor = cFactor*newFactor(0);
    
    return factor;
}


int SocketSeries::sendSelf(int commitTag, Channel& theChannel)
{
    return -1;
}


int SocketSeries::recvSelf(int commitTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    return -1;
}


void SocketSeries::Print(OPS_Stream& s, int flag)
{
    s << "Socket Time Series" << endln;
    s << "\tFactor: " << cFactor << endln;
}
