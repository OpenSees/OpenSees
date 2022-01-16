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

// $Revision$
// $Date$
// $URL$

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 04/12
// Revision: A
//
// Description: This file contains the class implementation for the
// VelDepMultiLinear friction model.

#include <VelDepMultiLinear.h>
#include <Channel.h>
#include <Information.h>
#include <elementAPI.h>

#include <math.h>


void * OPS_ADD_RUNTIME_VPV(OPS_VelDepMultiLinear)
{
    // pointer to a friction model that will be returned
    FrictionModel *theFrnMdl = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc < 7) {
        opserr << "WARNING invalid number of arguments\n";
        opserr << "Want: frictionModel VelDepMultiLinear tag ";
        opserr << "-vel velocityPoints -frn frictionPoints  ";
        opserr << "(with at least two friction-velocity points)";
        return 0;
    }
    
    int tag[1];
    double velData[64];
    double frnData[64];
    const char *paraStr;
    int numData = 1;
    if (OPS_GetIntInput(&numData, tag) != 0) {
        opserr << "WARNING invalid tag for frictionModel VelDepMultiLinear\n";
        return 0;
    }
    
    // get velocity data points
    numData = (argc - 3) / 2;
    paraStr = OPS_GetString();
    if (strcmp(paraStr, "-vel") == 0) {
        if (OPS_GetDoubleInput(&numData, velData) != 0) {
            opserr << "WARNING invalid velocityPoints\n";
            opserr << "frictionModel VelDepMultiLinear: " << tag[0] << endln;
            return 0;
        }
    }
    else {
        opserr << "WARNING expecting -vel but got " << paraStr << endln;
        opserr << "frictionModel VelDepMultiLinear: " << tag[0] << endln;
        return 0;
    }
    Vector velPts(velData, numData);
    
    // get friction data points
    paraStr = OPS_GetString();
    if (strcmp(paraStr, "-frn") == 0) {
        if (OPS_GetDoubleInput(&numData, frnData) != 0) {
            opserr << "WARNING invalid frictionPoints\n";
            opserr << "frictionModel VelDepMultiLinear: " << tag[0] << endln;
            return 0;
        }
    }
    else {
        opserr << "WARNING expecting -frn but got " << paraStr << endln;
        opserr << "frictionModel VelDepMultiLinear: " << tag[0] << endln;
        return 0;
    }
    Vector frnPts(frnData, numData);
    
    // parsing was successful, allocate the friction model
    theFrnMdl = new VelDepMultiLinear(tag[0], velPts, frnPts);
    if (theFrnMdl == 0) {
        opserr << "WARNING could not create frictionModel of type VelDepMultiLinear\n";
        return 0;
    }
    
    return theFrnMdl;
}


VelDepMultiLinear::VelDepMultiLinear()
    : FrictionModel(0, FRN_TAG_VelDepMultiLinear),
    velocityPoints(1), frictionPoints(1),
    trialID(0), trialIDmin(0), trialIDmax(0),
    numDataPoints(0), mu(0.0), DmuDvel(0.0)
{
    // does nothing
}


VelDepMultiLinear::VelDepMultiLinear(int tag,
    const Vector &velPts, const Vector &frnPts)
    : FrictionModel(tag, FRN_TAG_VelDepMultiLinear),
    velocityPoints(velPts), frictionPoints(frnPts),
    trialID(0), trialIDmin(0), trialIDmax(0),
    mu(0.0), DmuDvel(0.0)
{
    numDataPoints = velocityPoints.Size();
    if (numDataPoints != frictionPoints.Size())  {
        opserr << "VelDepMultiLinear::VelDepMultiLinear() "
            << "- velocity and friction arrays do not have same length.\n";
        exit(-1);
    }
    trialIDmax = numDataPoints - 2;
    
    // check that velocity and friction points are non-negative
    for (int i=0; i<numDataPoints; i++)  {
        if (velocityPoints(i) < 0.0 || frictionPoints(i) < 0.0)  {
            opserr << "VelDepMultiLinear::VelDepMultiLinear - "
                << "the velocity and friction points have to be positive.\n";
            exit(-1);
        }
    }
    
    // check that velocity points are monotonically increasing
    for (int i=0; i<numDataPoints-1; i++)  {
        if (velocityPoints(i) >= velocityPoints(i+1))  {
            opserr << "VelDepMultiLinear::VelDepMultiLinear - "
                << "the velocity points have to increase monotonically.\n";
            exit(-1);
        }
    }
    
    // initialize variables
    this->revertToStart();
}


VelDepMultiLinear::~VelDepMultiLinear()
{
    // does nothing
}


int VelDepMultiLinear::setTrial(double normalForce, double velocity)
{
    trialN   = normalForce;
    trialVel = velocity;
    double absTrialVel = fabs(trialVel);
    
    // find the current interval
    double vel1 = velocityPoints(trialID);
    double vel2 = velocityPoints(trialID+1);
    if (absTrialVel >= vel2 && trialID < trialIDmax)  {
        while (absTrialVel >= vel2 && trialID < trialIDmax)  {
            trialID++;
            vel1 = vel2;
            vel2 = velocityPoints(trialID+1);
        }
    } else if (absTrialVel < vel1 && trialID > trialIDmin)  {
        while (absTrialVel <= vel1 && trialID > trialIDmin)  {
            trialID--;
            vel2 = vel1;
            vel1 = velocityPoints(trialID);
        }
    }
    double mu1 = frictionPoints(trialID);
    double mu2 = frictionPoints(trialID+1);
    
    // get derivative of COF wrt velocity for the selected interval
    DmuDvel = (mu2-mu1)/(vel2-vel1);
    
    // get the COF for the selected interval
    mu = mu1 + DmuDvel*(absTrialVel-vel1);
    
    return 0;
}


double VelDepMultiLinear::getFrictionForce()
{
    if (trialN > 0.0)
        return mu*trialN;
    else
        return 0.0;
}


double VelDepMultiLinear::getFrictionCoeff()
{
    return mu;
}


double VelDepMultiLinear::getDFFrcDNFrc()
{
    if (trialN >= 0.0)
        return mu;
    else
        return 0.0;
}


double VelDepMultiLinear::getDFFrcDVel()
{
    if (trialN > 0.0)
        return DmuDvel*trialN;
    else
        return 0.0;
}


int VelDepMultiLinear::commitState()
{
    return 0;
}


int VelDepMultiLinear::revertToLastCommit()
{
    return 0;
}


int VelDepMultiLinear::revertToStart()
{
    trialN   = 0.0;
    trialVel = 0.0;
    trialID  = 0;
    
    // find the current interval
    double vel1 = velocityPoints(trialID);
    double vel2 = velocityPoints(trialID+1);
    if (trialVel >= vel2 && trialID < trialIDmax)  {
        while (trialVel >= vel2 && trialID < trialIDmax)  {
            trialID++;
            vel1 = vel2;
            vel2 = velocityPoints(trialID+1);
        }
    } else if (trialVel < vel1 && trialID > trialIDmin)  {
        while (trialVel <= vel1 && trialID > trialIDmin)  {
            trialID--;
            vel2 = vel1;
            vel1 = velocityPoints(trialID);
        }
    }
    double mu1 = frictionPoints(trialID);
    double mu2 = frictionPoints(trialID+1);
    
    // get derivative of COF wrt velocity for the selected interval
    DmuDvel = (mu2-mu1)/(vel2-vel1);
    
    // get the COF for the selected interval
    mu = mu1 + DmuDvel*(trialVel-vel1);
    
    return 0;
}


FrictionModel* VelDepMultiLinear::getCopy()
{
    VelDepMultiLinear *theCopy = new VelDepMultiLinear(this->getTag(),
        velocityPoints, frictionPoints);
    theCopy->trialN   = trialN;
    theCopy->trialVel = trialVel;
    theCopy->mu       = mu;
    theCopy->DmuDvel  = DmuDvel;
    
    return theCopy;
}


int VelDepMultiLinear::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(4);
    data(0) = this->getTag();
    data(1) = trialIDmin;
    data(2) = trialIDmax;
    data(3) = numDataPoints;
    
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    res += theChannel.sendVector(this->getDbTag(), cTag, velocityPoints);
    res += theChannel.sendVector(this->getDbTag(), cTag, frictionPoints);
    if (res < 0) 
        opserr << "VelDepMultiLinear::sendSelf() - failed to send data.\n";
    
    return res;
}


int VelDepMultiLinear::recvSelf(int cTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(4);
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0)  {
        opserr << "VelDepMultiLinear::recvSelf() - failed to receive data.\n";
        this->setTag(0);      
        trialIDmin    = 0;
        trialIDmax    = 0;
        numDataPoints = 0;
    }
    else  {
        this->setTag((int)data(0));
        trialIDmin    = (int)data(1);
        trialIDmax    = (int)data(2);
        numDataPoints = (int)data(3);
        
        // receive the velocity and friction arrays
        velocityPoints.resize(numDataPoints);
        frictionPoints.resize(numDataPoints);
        res += theChannel.recvVector(this->getDbTag(), cTag, velocityPoints);
        res += theChannel.recvVector(this->getDbTag(), cTag, frictionPoints);
        if (res < 0) 
            opserr << "VelDepMultiLinear::recvSelf() - failed to receive arrays.\n";
    }
    
    return res;
}


void VelDepMultiLinear::Print(OPS_Stream &s, int flag)
{
    s << "VelDepMultiLinear tag: " << this->getTag() << endln;
    s << "  velocityPoints: " << velocityPoints << endln;
    s << "  frictionPoints: " << frictionPoints << endln;
}
