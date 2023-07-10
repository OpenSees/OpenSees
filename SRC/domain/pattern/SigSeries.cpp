
// Written: Codi McKee (Texas A&M University)
// Created: 14-Jun-2021


#include <SigSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>

#include <math.h>

#include <elementAPI.h>


void* OPS_SigSeries(void)
{
    // Pointer to a uniaxial material that will be returned
    TimeSeries* theSeries = 0;

    int numRemainingArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingArgs < 3) {
        opserr << "WARNING: invalid num args Sig <tag?> $tStart $tFinish $a $b <-factor cFactor> \n";
        return 0;
    }

    int tag = 0;      // default tag = 0
    double dData[5];
    dData[4] = 1.0;

    int numData = 0;

    // get tag if provided
    if (numRemainingArgs == 5) {
        numData = 1;
        if (OPS_GetIntInput(&numData, &tag) != 0) {
            opserr << "WARNING invalid series tag in SigSeries tag?" << endln;
            return 0;
        }
        numRemainingArgs -= 1;
    }

    numData = 4;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid double data in SigSeries with tag: " << tag << endln;
        return 0;
    }

    // parse the optional args
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* argvS = OPS_GetString();

        if (strcmp(argvS, "-factor") == 0) {

            numData = 1;
            if (OPS_GetDouble(&numData, &dData[4]) != 0) {
                opserr << "WARNING invalid factor in SigSeries with tag?" << tag << endln;
                return 0;
            }
        }
        else {
            opserr << "WARNING unknown option: " << argvS << "  in SigSeries with tag?" << tag << endln;
            return 0;
        }
    }

    theSeries = new SigSeries(tag, dData[0], dData[1], dData[4], dData[2], dData[3]);

    if (theSeries == 0) {
        opserr << "WARNING ran out of memory creating Trig Series with tag: " << tag << "\n";
        return 0;
    }

    return theSeries;
}


SigSeries::SigSeries(int tag,
    double startTime,
    double finishTime, double cFactor_, double c1_, double c2_)
    : TimeSeries(tag, TSERIES_TAG_SigSeries),
    tStart(startTime), tFinish(finishTime), cFactor(cFactor_), c1(c1_), c2(c2_)
{

}


SigSeries::SigSeries()
    : TimeSeries(TSERIES_TAG_SigSeries),
    tStart(0.0), tFinish(0.0), cFactor(0.0), c1(0.0), c2(0.0)
{
    // does nothing
}


SigSeries::~SigSeries()
{
    // does nothing
}

TimeSeries* SigSeries::getCopy()
{
    return new SigSeries(this->getTag(), tStart, tFinish, cFactor, c1, c2);
}


double SigSeries::getFactor(double pseudoTime)
{

    if (pseudoTime <= tStart) {
        return 0.0;
    }
    else {
        return (1.0 - exp(-c1 * pow(pseudoTime - tStart, c2))) * cFactor;
        //  return (1.0 + exp(-c1*(tFinish-tStart) *(pseudoTime - (tFinish - tStart)*c2))) * cFactor;
    }
}


int SigSeries::sendSelf(int commitTag, Channel& theChannel)
{
    int dbTag = this->getDbTag();
    Vector data(6);
    data(0) = cFactor;
    data(1) = tStart;
    data(2) = tFinish;
    data(3) = c1;
    data(4) = c2;


    int result = theChannel.sendVector(dbTag, commitTag, data);
    if (result < 0) {
        opserr << "SigSeries::sendSelf() - channel failed to send data\n";
        return result;
    }

    return 0;
}


int SigSeries::recvSelf(int commitTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int dbTag = this->getDbTag();
    Vector data(6);
    int result = theChannel.recvVector(dbTag, commitTag, data);
    if (result < 0) {
        opserr << "SigSeries::recvSelf() - channel failed to receive data\n";
        cFactor = 1.0;
        tStart = 0.0;
        tFinish = 0.0;
        return result;
    }
    cFactor = data(0);
    tStart = data(1);
    tFinish = data(2);
    c1 = data(3);
    c2 = data(4);

    return 0;
}


void SigSeries::Print(OPS_Stream& s, int flag)
{
    s << "Sigmoid Series" << endln;
    s << "\tFactor: " << cFactor << endln;
    s << "\ttStart: " << tStart << endln;
    s << "\ttFinish: " << tFinish << endln;
    s << "\ttc1: " << c1 << endln;
    s << "\ttc2: " << c2 << endln;
}
