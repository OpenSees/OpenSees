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

// Developed: Nhan D. Dao (nhan.unr@gmail.com)
// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 08/14
// Revision: A
//
// Description: This file contains the class implementation for the
// VelNormalFrcDep friction model.

#include <VelNormalFrcDep.h>
#include <Channel.h>
#include <Information.h>
#include <elementAPI.h>

#include <math.h>


void * OPS_ADD_RUNTIME_VPV(OPS_VelNormalFrcDep)
{
    // pointer to a friction model that will be returned
    FrictionModel *theFrnMdl = 0;
    
    if (OPS_GetNumRemainingInputArgs() < 9) {
        opserr << "WARNING invalid number of arguments\n";
        opserr << "Want: frictionModel VelNormalFrcDep tag aSlow nSlow aFast nFast alpha0 alpha1 alpha2 maxMuFact\n";
        return 0;
    }
    
    int tag[1];
    double dData[8];
    int numData = 1;
    if (OPS_GetIntInput(&numData, tag) != 0) {
        opserr << "WARNING invalid tag for frictionModel VelNormalFrcDep\n";
        return 0;
    }
    numData = 8;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid data for frictionModel VelNormalFrcDep " << tag[0] << endln;
        return 0;
    }
    
    // parsing was successful, allocate the friction model
    theFrnMdl = new VelNormalFrcDep(tag[0], dData[0], dData[1], dData[2],
        dData[3], dData[4], dData[5], dData[6], dData[7]);
    if (theFrnMdl == 0) {
        opserr << "WARNING could not create frictionModel of type VelNormalFrcDep\n";
        return 0;
    }
    
    return theFrnMdl;
}


VelNormalFrcDep::VelNormalFrcDep()
    : FrictionModel(0, FRN_TAG_VelNormalFrcDep),
    aSlow(0.0), nSlow(1.0), aFast(0.0), nFast(1.0),
    alpha0(0.0), alpha1(0.0), alpha2(0.0), maxMuFact(2.0),
    mu(0.0), DmuDn(0.0), DmuDvel(0.0)
{
    // does nothing
}


VelNormalFrcDep::VelNormalFrcDep(int tag,
    double aslow, double nslow, double afast, double nfast,
    double _alpha0, double _alpha1, double _alpha2, double maxmufact)
    : FrictionModel(tag, FRN_TAG_VelNormalFrcDep),
    aSlow(aslow), nSlow(nslow), aFast(afast), nFast(nfast),
    alpha0(_alpha0), alpha1(_alpha1), alpha2(_alpha2), maxMuFact(maxmufact),
    mu(0.0), DmuDn(0.0), DmuDvel(0.0)
{
    // check that constants are positive and not zero
    if (aSlow <= 0.0  || aFast <= 0.0)  {
        opserr << "VelNormalFrcDep::VelNormalFrcDep - "
            << "the aSlow & aFast constants have to be positive.\n";
        exit(-1);
    }
    // check that exponents are <= 1
    if (nSlow > 1.0 || nFast > 1.0)  {
        opserr << "VelNormalFrcDep::VelNormalFrcDep - "
            << "the exponents n have to be <= 1.0.\n";
        exit(-1);
    }
    
    // initialize variables
    this->revertToStart();
}


VelNormalFrcDep::~VelNormalFrcDep()
{
    // does nothing
}


int VelNormalFrcDep::setTrial(double normalForce, double velocity)
{
    trialN   = normalForce;
    trialVel = velocity;
    
    // determine slow and fast COF
    double muSlow = aSlow*pow(trialN,nSlow-1.0);
    double muFast = aFast*pow(trialN,nFast-1.0);
    
    // determine rate parameter
    double transRate = alpha0 + alpha1*trialN + alpha2*pow(trialN,2.0);
    
    // get COF for given normal force and trial velocity
    double temp1 = exp(-transRate*fabs(trialVel));
    double temp2 = (muFast - muSlow)*temp1;
    mu = muFast - temp2;
    if (mu > maxMuFact*muFast || trialN <= 0.0)
        mu = maxMuFact*muFast;
    
    // get derivative of COF wrt normal force
    double DmuSlowDn = aSlow*(nSlow-1.0)*pow(trialN,nSlow-2.0);
    double DmuFastDn = aFast*(nFast-1.0)*pow(trialN,nFast-2.0);
    DmuDn = DmuFastDn - (DmuFastDn - DmuSlowDn)*temp1
          + (alpha1 + 2.0*alpha2*trialN)*fabs(trialVel)*temp2;
    
    // get derivative of COF wrt velocity
    if (trialVel != 0.0)
        DmuDvel = transRate*trialVel/fabs(trialVel)*temp2;
    else
        DmuDvel = 0.0;
    
    return 0;
}


double VelNormalFrcDep::getFrictionForce()
{
    if (trialN > 0.0)
        return mu*trialN;
    else
        return 0.0;
}


double VelNormalFrcDep::getFrictionCoeff()
{
    return mu;
}


double VelNormalFrcDep::getDFFrcDNFrc()
{
    if (trialN >= 0.0)
        return mu + DmuDn*trialN;
    else
        return 0.0;
}


double VelNormalFrcDep::getDFFrcDVel()
{
    if (trialN > 0.0)
        return DmuDvel*trialN;
    else
        return 0.0;
}


int VelNormalFrcDep::commitState()
{
    return 0;
}


int VelNormalFrcDep::revertToLastCommit()
{
    return 0;
}


int VelNormalFrcDep::revertToStart()
{
    trialN   = 0.0;
    trialVel = 0.0;
    mu       = aSlow*pow(trialN,(nSlow-1.0));
    DmuDn    = 0.0;
    DmuDvel  = 0.0;
    
    return 0;
}


FrictionModel* VelNormalFrcDep::getCopy()
{
    VelNormalFrcDep *theCopy = new VelNormalFrcDep(this->getTag(),
        aSlow, nSlow, aFast, nFast, alpha0, alpha1, alpha2, maxMuFact);
    theCopy->trialN   = trialN;
    theCopy->trialVel = trialVel;
    theCopy->mu       = mu;
    theCopy->DmuDn    = DmuDn;
    theCopy->DmuDvel  = DmuDvel;
    
    return theCopy;
}


int VelNormalFrcDep::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(9);
    data(0) = this->getTag();
    data(1) = aSlow;
    data(2) = nSlow;
    data(3) = aFast;
    data(4) = nFast;
    data(5) = alpha0;
    data(6) = alpha1;
    data(7) = alpha2;
    data(8) = maxMuFact;
    
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) 
        opserr << "VelNormalFrcDep::sendSelf() - failed to send data.\n";
    
    return res;
}


int VelNormalFrcDep::recvSelf(int cTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(9);
    
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0)  {
        opserr << "VelNormalFrcDep::recvSelf() - failed to receive data.\n";
        this->setTag(0);      
        aSlow     = 0.0;
        nSlow     = 1.0;
        aFast     = 0.0;
        nFast     = 1.0;
        alpha0    = 0.0;
        alpha1    = 0.0;
        alpha2    = 0.0;
        maxMuFact = 2.0;
    }
    else  {
        this->setTag((int)data(0));
        aSlow     = data(1);
        nSlow     = data(2);
        aFast     = data(3);
        nFast     = data(4);
        alpha0    = data(5);
        alpha1    = data(6);
        alpha2    = data(7);
        maxMuFact = data(8);
    }
    
    // initialize variables
    this->revertToStart();
    
    return res;
}


void VelNormalFrcDep::Print(OPS_Stream &s, int flag)
{
    s << "VelNormalFrcDep tag: " << this->getTag() << endln;
    s << "  aSlow: " << aSlow << "  nSlow: " << nSlow;
    s << "  aFast: " << aFast << "  nFast: " << nFast << endln;
    s << "  alpha0: " << alpha0 << "  alpha1: " << alpha1;
    s << "  alpha2: " << alpha2 << endln;
    s << "  maxMuFact: " << maxMuFact << endln;
}
