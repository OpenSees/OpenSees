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

// $Date: 2023-01-20 20:16:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/MPAccSeries.cpp,v $

// Written: Tang.S
// Created: 2023-01-20
//
// Purpose: This file contains the class implementation of MPAccSeries.
//
// What: "@(#) MPAccSeries.C, revA"


#include <MPAccSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>

#include <math.h>

#include <elementAPI.h>


void *OPS_MPAccSeries(void)
{
    // Pointer to a uniaxial material that will be returned
    TimeSeries *theSeries = 0;

    int numRemainingArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingArgs < 3) {
        opserr << "WARNING: invalid num args MPAcc <tag?> $tStart $tFinish $period <-AFactor> <-gammaMP> <-nuMP>\n";
        return 0;
    }

    int tag = 0;      // default tag = 0
    double dData[6];
    dData[3] = 1.0;   // default gammaMP = 1.0
    dData[4] = 90.0;   // default nuMP = 90.0
    dData[5] = 1.0;   // default AFactor = 1.0
    int numData = 0;

    // get tag if provided
    if (numRemainingArgs == 4 || numRemainingArgs == 6 || numRemainingArgs == 8 || numRemainingArgs == 10) {
        numData = 1;
        if (OPS_GetIntInput(&numData, &tag) != 0) {
            opserr << "WARNING invalid series tag in MPAcc tag?" << endln;
            return 0;
        }
        numRemainingArgs -= 1;
    }

    numData = 3;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid double data in MPAcc Series with tag: " << tag << endln;
        return 0;
    }
    numRemainingArgs -= 3;

    // parse the optional args
    while (numRemainingArgs > 1) {
      const char *argvS = OPS_GetString();

      if (strcmp(argvS,"-gammaMP") == 0) {
	numData = 1;
	if (OPS_GetDouble(&numData, &dData[3]) != 0) {
	  opserr << "WARNING invalid gamma in MPAcc Series with tag?" << tag << endln;
                return 0;
	}
      } else if (strcmp(argvS,"-nuMP") == 0) {
	numData = 1;
            if (OPS_GetDouble(&numData, &dData[4]) != 0) {
                opserr << "WARNING invalid nu in MPAcc Series with tag?" << tag << endln;
                return 0;
            }
        } else if (strcmp(argvS,"-AFactor") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[5]) != 0) {
                opserr << "WARNING invalid amplitude in MPAcc Series with tag?" << tag << endln;
                return 0;
            }
        } else {
            opserr << "WARNING unknown option: " << argvS << "  in MPAcc Series with tag?" << tag << endln;      
            return 0;
        }
        numRemainingArgs -= 2;
    }
    //double tData = dData[4] * asin(1.0) / 90;
    //double tF = dData[3] * dData[2];
    theSeries = new MPAccSeries(tag, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);

    if (theSeries == 0) {
        opserr << "WARNING ran out of memory creating MPAcc Series with tag: " << tag << "\n";
        return 0;
    }

    return theSeries;
}


MPAccSeries::MPAccSeries(int tag,
    double startTime, 
    double finishTime,
    double T, 
    double gamma, 
    double nu,
    double A)
    : TimeSeries(tag, TSERIES_TAG_MPAccSeries),
    tStart(startTime), tFinish(finishTime),
    period(T), gammaMP(gamma),
    nuMP(nu), AFactor(A)
{
    if (period == 0.0) {
        opserr << "MPAccSeries::MPAccSeries -- input period is zero, setting period to PI\n";
        period = 2*asin(1.0);
    }
}


MPAccSeries::MPAccSeries()
    : TimeSeries(TSERIES_TAG_MPAccSeries),
    tStart(0.0), tFinish(0.0),
    period(1.0), gammaMP(0.0),
    nuMP(0.0), AFactor(1.0)
{
    // does nothing
}


MPAccSeries::~MPAccSeries()
{
    // does nothing
}

TimeSeries *MPAccSeries::getCopy()
{
    return new MPAccSeries(this->getTag(), tStart, tFinish, period,
        gammaMP, nuMP, AFactor);
}


double MPAccSeries::getFactor(double pseudoTime)
{
    static double onepi = 2*asin(1.0);

    if (pseudoTime >= tStart && pseudoTime <= tFinish)  {
        double part2 = gammaMP * sin((2 * onepi * pseudoTime) / period - onepi * gammaMP + (nuMP * onepi)/180) * (1 - cos((2 * onepi * pseudoTime) / (period * gammaMP)));
        return ((AFactor * onepi) / (period * gammaMP)) * (sin((2 * onepi * pseudoTime) / (period * gammaMP)) * (cos((2 * onepi * pseudoTime) / period - onepi * gammaMP + (nuMP * onepi)/180)) - part2);
    }
    else
        return 0.0;
}


int MPAccSeries::sendSelf(int commitTag, Channel &theChannel)
{
    int dbTag = this->getDbTag();
    Vector data(6);
    data(0) = AFactor;
    data(1) = tStart;
    data(2) = tFinish;
    data(3) = period;
    data(4) = gammaMP;
    data(5) = nuMP;

    int result = theChannel.sendVector(dbTag,commitTag, data);
    if (result < 0) {
        opserr << "MPAccSeries::sendSelf() - channel failed to send data\n";
        return result;
    }

    return 0;
}


int MPAccSeries::recvSelf(int commitTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    int dbTag = this->getDbTag();
    Vector data(6);
    int result = theChannel.recvVector(dbTag,commitTag, data);
    if (result < 0) {
        opserr << "MPAccSeries::recvSelf() - channel failed to receive data\n";
        AFactor    = 1.0;
        tStart     = 0.0;
        tFinish    = 0.0;
        period     = 1.0;
        gammaMP = 0.0;
        nuMP  = 0.0;
        return result;
    }
    AFactor    = data(0);
    tStart     = data(1);
    tFinish    = data(2);
    period     = data(3);
    gammaMP = data(4);
    nuMP  = data(5);

    return 0;
}


void MPAccSeries::Print(OPS_Stream &s, int flag)
{
    s << "MPAcc Series" << endln;
    s << "\tAFactor: " << AFactor << endln;
    s << "\ttStart: " << tStart << endln;
    s << "\ttFinish: " << tFinish << endln;
    s << "\tPeriod: " << period << endln;
    s << "\tgammaMP: " << gammaMP << endln;
    s << "\tnuMP: " << nuMP << endln;
}
