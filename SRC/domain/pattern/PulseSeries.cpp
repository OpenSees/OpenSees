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
// $Date: 2010-02-04 00:34:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/PulseSeries.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/04
// Revision: A
//
// Purpose: This file contains the class implementation of PulseSeries.
//
// What: "@(#) PulseSeries.cpp, revA"


#include <PulseSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>

#include <math.h>

#include <elementAPI.h>

void *OPS_PulseSeries(void)
{
    // Pointer to a uniaxial material that will be returned
    TimeSeries *theSeries = 0;

    int numRemainingArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingArgs < 3) {
        opserr << " Pulse <tag?> tStart tFinish period <-width pulseWidth> <-phaseShift shift> <-factor cFactor> <-zeroShift shift>\n";
        return 0;
    }

    int tag = 0;      // default tag = 0
    double dData[7];
    dData[3] = 0.5;   // default width = 0.5
    dData[4] = 0.0;   // default phaseShift = 0.0
    dData[5] = 1.0;   // default factor = 1.0
    dData[6] = 0.0;   // default zeroShift = 0.0
    int numData = 0;

    // get tag if provided
    if (numRemainingArgs == 4 || numRemainingArgs == 6 || numRemainingArgs == 8 || numRemainingArgs == 10 || numRemainingArgs == 12) {
        numData = 1;
        if (OPS_GetIntInput(&numData, &tag) != 0) {
            opserr << "WARNING invalid series tag in Pulse tag?" << endln;
            return 0;
        }
        numRemainingArgs -= 1;
    }

    numData = 3;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid double data in Pulse Series with tag: " << tag << endln;
        return 0;
    }
    numRemainingArgs -= 2;

    // parse the optional args
    while (numRemainingArgs > 1) {
      const char *argvS = OPS_GetString();

        if (strcmp(argvS,"-width") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[3]) != 0) {
                opserr << "WARNING invalid width in Pulse Series with tag?" << tag << endln;
                return 0;
            }
        } else if (strcmp(argvS,"-shift") == 0 || strcmp(argvS,"-phaseShift") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[4]) != 0) {
                opserr << "WARNING invalid phase shift in Pulse Series with tag?" << tag << endln;
                return 0;
            }
        } else if (strcmp(argvS,"-factor") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[5]) != 0) {
                opserr << "WARNING invalid factor in Pulse Series with tag?" << tag << endln;
                return 0;
            }
        } else if (strcmp(argvS,"-zeroShift") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[6]) != 0) {
                opserr << "WARNING invalid zero shift in Pulse Series with tag?" << tag << endln;
                return 0;
            }
        } else {
            opserr << "WARNING unknown option: " << argvS << "  in Pulse Series with tag?" << tag << endln;      
            return 0;
        }
        numRemainingArgs -= 2;
    }

    // create the object
    theSeries = new PulseSeries(tag, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6]);

    if (theSeries == 0) {
        opserr << "WARNING ran out of memory creating Pulse Series with tag: " << tag << "\n";
        return 0;
    }

    return theSeries;
}


PulseSeries::PulseSeries(int tag,
    double startTime,
    double finishTime,
    double T,
    double pulseWidth,
    double phaseshift,
    double theFactor,
    double zeroshift)
    : TimeSeries(tag, TSERIES_TAG_PulseSeries),
    tStart(startTime), tFinish(finishTime),
    period(T), pWidth(pulseWidth),
    phaseShift(phaseshift), cFactor(theFactor),
    zeroShift(zeroshift)
{
    if (period == 0.0)  {
        opserr << "PulseSeries::PulseSeries -- input period is zero, setting period to 1\n";
        period = 1;
    }
}


PulseSeries::PulseSeries()
    : TimeSeries(TSERIES_TAG_PulseSeries),
    tStart(0.0), tFinish(0.0),
    period(1.0), pWidth(0.5),
    phaseShift(0.0), cFactor(1.0),
    zeroShift(0.0)
{
    // does nothing
}


PulseSeries::~PulseSeries()
{
    // does nothing
}


TimeSeries *PulseSeries::getCopy()
{
    return new PulseSeries(this->getTag(), tStart, tFinish, period,
        pWidth, phaseShift, cFactor, zeroShift);
}


double PulseSeries::getFactor(double pseudoTime)
{	
    if (tStart <= pseudoTime && pseudoTime <= tFinish)  {
        double k = (pseudoTime+phaseShift-tStart)/period - floor((pseudoTime+phaseShift-tStart)/period);
        if (k < pWidth)
            return cFactor + zeroShift;
        else if (k < 1.00)
            return zeroShift;
        else
            return 0.0;
    }
    else
        return 0.0;
}


int PulseSeries::sendSelf(int commitTag, Channel &theChannel)
{
    int dbTag = this->getDbTag();
    Vector data(7);
    data(0) = cFactor;
    data(1) = tStart;	
    data(2) = tFinish;
    data(3) = period;
    data(4) = pWidth;
    data(5) = phaseShift;
    data(6) = zeroShift;

    int result = theChannel.sendVector(dbTag,commitTag, data);
    if (result < 0)  {
        opserr << "PulseSeries::sendSelf() - channel failed to send data\n";
        return result;
    }
    return 0;
}


int PulseSeries::recvSelf(int commitTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    int dbTag = this->getDbTag();
    Vector data(7);
    int result = theChannel.recvVector(dbTag,commitTag, data);
    if (result < 0)  {
        opserr << "PulseSeries::sendSelf() - channel failed to receive data\n";
        cFactor    = 1.0;
        tStart     = 0.0;
        tFinish    = 0.0;
        period     = 1.0;
        pWidth     = 0.5;
        phaseShift = 0.0;
        zeroShift  = 0.0;
        return result;
    }
    cFactor    = data(0);
    tStart     = data(1);
    tFinish    = data(2);
    period     = data(3);
    pWidth     = data(4);
    phaseShift = data(5);
    zeroShift  = data(6);

    return 0;
}


void PulseSeries::Print(OPS_Stream &s, int flag)
{
    s << "Pulse Series" << endln;
    s << "\tFactor: " << cFactor << endln;
    s << "\ttStart: " << tStart << endln;
    s << "\ttFinish: " << tFinish << endln;
    s << "\tPeriod: " << period << endln;
    s << "\tPulse Width: " << pWidth << endln;
    s << "\tPhase Shift: " << phaseShift << endln;
    s << "\tZero Shift: " << zeroShift << endln;
}
