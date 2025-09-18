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
// $Date: 2010-04-06 20:16:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TriangleSeries.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/04
// Revision: A
//
// Purpose: This file contains the class implementation of TriangleSeries.
//
// What: "@(#) TriangleSeries.cpp, revA"


#include <TriangleSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>

#include <math.h>

#include <elementAPI.h>
#define OPS_Export 


OPS_Export void *OPS_TriangleSeries()
{
    // Pointer to a uniaxial material that will be returned
    TimeSeries *theSeries = 0;

    int numRemainingArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingArgs < 3) {
        opserr << "WARNING: invalid num args Triangle <tag?> tStart tFinish period <-phaseShift shift> <-factor cFactor> <-zeroShift shift>\n";
        return 0;
    }

    int tag = 0;      // default tag = 0
    double dData[6];
    dData[3] = 0.0;   // default phaseShift = 0.0
    dData[4] = 1.0;   // default cFactor = 1.0
    dData[5] = 0.0;   // default zeroShift = 0.0
    int numData = 0;

    // get tag if provided
    if (numRemainingArgs == 4 || numRemainingArgs == 6 || numRemainingArgs == 8 || numRemainingArgs == 10) {
        numData = 1;
        if (OPS_GetIntInput(&numData, &tag) != 0) {
            opserr << "WARNING invalid series tag in Triangle tag?" << endln;
            return 0;
        }
        numRemainingArgs -= 1;
    }

    numData = 3;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid double data in Triangle Series with tag: " << tag << endln;
        return 0;
    }
    numRemainingArgs -= 3;

    // parse the optional args
    while (numRemainingArgs > 1) {
      const char *argvS = OPS_GetString();

        if (strcmp(argvS,"-shift") == 0 || strcmp(argvS,"-phaseShift") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[3]) != 0) {
                opserr << "WARNING invalid phase shift in Triangle Series with tag?" << tag << endln;
                return 0;
            }
        } else if (strcmp(argvS,"-factor") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[4]) != 0) {
                opserr << "WARNING invalid factor in Triangle Series with tag?" << tag << endln;
                return 0;
            }
        } else if (strcmp(argvS,"-zeroShift") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[5]) != 0) {
                opserr << "WARNING invalid zero shift in Triangle Series with tag?" << tag << endln;
                return 0;
            }
        } else {
            opserr << "WARNING unknown option: " << argvS << "  in Triangle Series with tag?" << tag << endln;      
            return 0;
        }
        numRemainingArgs -= 2;
    }

    theSeries = new TriangleSeries(tag, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);

    if (theSeries == 0) {
        opserr << "WARNING ran out of memory creating Triangle Series with tag: " << tag << "\n";
        return 0;
    }

    return theSeries;
}


TriangleSeries::TriangleSeries(int tag,
    double startTime, 
    double finishTime,
    double T, 
    double phaseshift, 
    double theFactor,
    double zeroshift)
    : TimeSeries(tag, TSERIES_TAG_TriangleSeries),
    tStart(startTime), tFinish(finishTime),
    period(T), phaseShift(phaseshift),
    cFactor(theFactor), zeroShift(zeroshift)
{
    if (period == 0.0)  {
        opserr << "TriangleSeries::TriangleSeries -- input period is zero, setting period to 1\n";
        period = 1;
    }
}


TriangleSeries::TriangleSeries()
    : TimeSeries(TSERIES_TAG_TriangleSeries),
    tStart(0.0), tFinish(0.0),
    period(1.0), phaseShift(0.0),
    cFactor(1.0), zeroShift(0.0)
{
    // does nothing
}


TriangleSeries::~TriangleSeries()
{
    // does nothing
}


TimeSeries *TriangleSeries::getCopy()
{
    return new TriangleSeries(this->getTag(), tStart, tFinish, period,
        phaseShift, cFactor, zeroShift);
}


double TriangleSeries::getFactor(double pseudoTime)
{
    if (tStart <= pseudoTime && pseudoTime <= tFinish)  {
        double slope = cFactor/(period/4);
        double phi = phaseShift - zeroShift/slope;
        double k = (pseudoTime+phi-tStart)/period - floor((pseudoTime+phi-tStart)/period);
        if (k < 0.25)
            return slope*k*period + zeroShift;
        else if (k < 0.75)
            return cFactor - slope*(k-0.25)*period + zeroShift;
        else if (k < 1.00)
            return -cFactor + slope*(k-0.75)*period + zeroShift;
        else
            return 0.0;
    }
    else
        return 0.0;
}


int TriangleSeries::sendSelf(int commitTag, Channel &theChannel)
{
    int dbTag = this->getDbTag();
    Vector data(6);
    data(0) = cFactor;
    data(1) = tStart;
    data(2) = tFinish;
    data(3) = period;
    data(4) = phaseShift;
    data(5) = zeroShift;

    int result = theChannel.sendVector(dbTag,commitTag, data);
    if (result < 0)  {
        opserr << "TriangleSeries::sendSelf() - channel failed to send data\n";
        return result;
    }

    return 0;
}


int TriangleSeries::recvSelf(int commitTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    int dbTag = this->getDbTag();
    Vector data(6);
    int result = theChannel.recvVector(dbTag,commitTag, data);
    if (result < 0)  {
        opserr << "TriangleSeries::recvSelf() - channel failed to receive data\n";
        cFactor    = 1.0;
        tStart     = 0.0;
        tFinish    = 0.0;
        period     = 1.0;
        phaseShift = 0.0;
        zeroShift  = 0.0;
        return result;
    }
    cFactor    = data(0);
    tStart     = data(1);
    tFinish    = data(2);
    period     = data(3);
    phaseShift = data(4);
    zeroShift  = data(5);

    return 0;
}


void TriangleSeries::Print(OPS_Stream &s, int flag)
{
    s << "Triangle Series" << endln;
    s << "\tFactor: " << cFactor << endln;
    s << "\ttStart: " << tStart << endln;
    s << "\ttFinish: " << tFinish << endln;
    s << "\tPeriod: " << period << endln;
    s << "\tPhase Shift: " << phaseShift << endln;
    s << "\tZero Shift: " << zeroShift << endln;
}
