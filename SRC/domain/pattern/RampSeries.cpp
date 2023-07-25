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

#include <RampSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>

#include <math.h>

#include <elementAPI.h>


void* OPS_RampSeries(void)
{
    // Pointer to a uniaxial material that will be returned
    TimeSeries* theSeries = 0;

    int numRemainingArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingArgs < 2) {
        opserr << "WARNING: invalid num args RampSeries <tag?> $tStart $tRamp <-offset $offset?> <-smooth $smoothness? > <-factor $cFactor> \n";
        return 0;
    }

    int tag = 0;      // default tag = 0
    double dData[5]{};
    dData[4] = 1.0; // cFactor default = 1.0

    int numData = 0;

    if (numRemainingArgs == 9 || numRemainingArgs == 7 || numRemainingArgs == 5 || numRemainingArgs == 3) {
        //Parse Tag
        numData = 1;
        if (OPS_GetIntInput(&numData, &tag) != 0) {
            opserr << "WARNING invalid series tag in RampSeries tag " << endln;
            return 0;
        }
    }

    //Parse tStart
    numData = 1;
    if (OPS_GetDouble(&numData, &dData[0]) != 0) {
        opserr << "WARNING invalid tStart? in RampSeries with tag " << tag << endln;
        return 0;
    }
    //Parse tRamp
    if (OPS_GetDouble(&numData, &dData[1]) != 0) {
        opserr << "WARNING invalid tRamp? in RampSeries with tag " << tag << endln;
        return 0;
    }
    if (dData[1] <= 0.0) {
        opserr << "WARNING invalid tRamp? in RampSeries with tag " << tag << " Duration must be larger than 0" << endln;
        return 0;
    }

    // Parse the optional args
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* argvS = OPS_GetString();

        if (strcmp(argvS, "-offset") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[2]) != 0) {
                opserr << "WARNING invalid offset? in RampSeries with tag " << tag << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-smooth") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[3]) != 0) {
                opserr << "WARNING invalid smoothness? in RampSeries with tag " << tag << " Expected value ranging between 0 and 1" << endln;
                return 0;
            }
            if (dData[3] < 0.0 || dData[3] > 1.0) {
                opserr << "WARNING invalidtRamp? in RampSeries with tag " << tag << " Smoothness must be a value ranging from 0 and 1" << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-factor") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[4]) != 0) {
                opserr << "WARNING invalid factor in RampSeries with tag " << tag << endln;
                return 0;
            }
        }
        else {
            opserr << "WARNING unknown option: " << argvS << "  in RampSeries with tag " << tag << endln;
            return 0;
        }
    }

    theSeries = new RampSeries(tag, dData[0], dData[1], dData[2], dData[3], dData[4]);

    if (theSeries == 0) {
        opserr << "WARNING ran out of memory creating Ramp Series with tag: " << tag << "\n";
        return 0;
    }

    return theSeries;
}


RampSeries::RampSeries(int tag, double startTime, double rampTime, double offset, double smooth, double factor)
    : TimeSeries(tag, TSERIES_TAG_RampSeries),
    tStart(startTime), tRamp(rampTime), offsetFact(offset),smoothFact(smooth),cFactor(factor)
{

}


RampSeries::RampSeries()
    : TimeSeries(TSERIES_TAG_RampSeries),
    tStart(0.0), tRamp(0.0), offsetFact(0.0), smoothFact(0.0), cFactor(1.0)
{
    // does nothing
}


RampSeries::~RampSeries()
{
    // does nothing
}

TimeSeries* RampSeries::getCopy()
{
    return new RampSeries(this->getTag(), tStart, tRamp, offsetFact, smoothFact, cFactor);
}


double RampSeries::getFactor(double pseudoTime)
{

    if (pseudoTime <= tStart ) {
        return 0.0 + offsetFact;
    }
    else if (pseudoTime <= tStart  + (smoothFact * tRamp / 2.0)) {
        return (offsetFact + cFactor * (pow((pseudoTime - tStart ) / tRamp, 2.0) * (2.0 / (smoothFact * (2.0 - smoothFact)))));
    }
    else if (pseudoTime <= tStart + tRamp - (smoothFact * tRamp / 2.0)) {
        return (offsetFact + cFactor * (0.5 + (((pseudoTime - tStart  - (tRamp / 2.0)) / tRamp) * (2.0 / (2.0 - smoothFact)))));
    }
    else if (pseudoTime <= tStart + tRamp) {
        return (offsetFact + cFactor * (1.0 - (pow((pseudoTime - tStart - tRamp) / tRamp, 2.0) * (2.0 / (smoothFact * (2.0 - smoothFact))))));
    }
    else {
        return (offsetFact + cFactor * (1.0));
    }
}


int RampSeries::sendSelf(int commitTag, Channel& theChannel)
{
    int dbTag = this->getDbTag();
    Vector data(5);
    data(0) = tStart;
    data(1) = tRamp;
    data(2) = offsetFact;
    data(3) = smoothFact;
    data(4) = cFactor;



    int result = theChannel.sendVector(dbTag, commitTag, data);
    if (result < 0) {
        opserr << "RampSeries::sendSelf() - channel failed to send data\n";
        return result;
    }

    return 0;
}


int RampSeries::recvSelf(int commitTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int dbTag = this->getDbTag();
    Vector data(5);
    int result = theChannel.recvVector(dbTag, commitTag, data);
    if (result < 0) {
        opserr << "RampSeries::recvSelf() - channel failed to receive data\n";
        cFactor = 1.0;
        tStart = 0.0;
        tRamp = 0.0;
        offsetFact = 0.0;
        smoothFact = 0.0;
        return result;
    }
    tStart = data(0);
    tRamp = data(1);
    offsetFact = data(2);
    smoothFact = data(3);
    cFactor = data(4);

    return 0;
}


void RampSeries::Print(OPS_Stream& s, int flag)
{
    s << "Ramp Series" << endln;
    s << "\tFactor: " << cFactor << endln;
    s << "\ttStart: " << tStart << endln;
    s << "\ttRamp: " << tRamp << endln;
    s << "\toffsetFactor: " << offsetFact << endln;
    s << "\tsmoothFactor: " << smoothFact << endln;
}
