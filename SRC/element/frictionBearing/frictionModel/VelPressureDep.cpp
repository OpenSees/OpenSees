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
// Created: 02/06
// Revision: A
//
// Description: This file contains the class implementation for the
// VelPressureDep friction model.

#include <VelPressureDep.h>
#include <Channel.h>
#include <Information.h>
#include <elementAPI.h>

#include <math.h>


void * OPS_ADD_RUNTIME_VPV(OPS_VelPressureDep)
{
    // pointer to a friction model that will be returned
    FrictionModel *theFrnMdl = 0;
    
    if (OPS_GetNumRemainingInputArgs() < 7) {
        opserr << "WARNING invalid number of arguments\n";
        opserr << "Want: frictionModel VelPressureDep tag muSlow muFast0 A deltaMu alpha transRate\n";
        return 0;
    }
    
    int tag[1];
    double dData[6];
    int numData = 1;
    if (OPS_GetIntInput(&numData, tag) != 0) {
        opserr << "WARNING invalid tag for frictionModel VelPressureDep\n";
        return 0;
    }
    numData = 6;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid data for frictionModel VelPressureDep " << tag[0] << endln;
        return 0;
    }
    
    // parsing was successful, allocate the friction model
    theFrnMdl = new VelPressureDep(tag[0], dData[0], dData[1], dData[2],
        dData[3], dData[4], dData[5]);
    if (theFrnMdl == 0) {
        opserr << "WARNING could not create frictionModel of type VelPressureDep\n";
        return 0;
    }
    
    return theFrnMdl;
}


VelPressureDep::VelPressureDep()
    : FrictionModel(0, FRN_TAG_VelPressureDep),
    muSlow(0.0), muFast0(0.0), A(0.0),
    deltaMu(0.0), alpha(0.0), transRate(0.0),
    mu(0.0), DmuDn(0.0), DmuDvel(0.0)
{
    // does nothing
}


VelPressureDep::VelPressureDep(int tag, double muslow,
    double mufast0, double a, double deltamu, double _alpha, double transrate)
    : FrictionModel(tag, FRN_TAG_VelPressureDep),
    muSlow(muslow), muFast0(mufast0), A(a), deltaMu(deltamu),
    alpha(_alpha), transRate(transrate),
    mu(0.0), DmuDn(0.0), DmuDvel(0.0)
{
    // check that COF are positive and not zero
    if (muSlow <= 0.0  || muFast0 <= 0.0)  {
        opserr << "VelPressureDep::VelPressureDep - "
            << "the friction coefficients have to be positive.\n";
        exit(-1);
    }
    // check that A is positive and not zero
    if (A <= 0.0)  {
        opserr << "VelPressureDep::VelPressureDep - "
            << "the nominal contact area has to be positive.\n";
        exit(-1);
    }
    // check that transRate is positive
    if (transRate < 0.0)  {
        opserr << "VelPressureDep::VelPressureDep - "
            << "the transition rate has to be positive.\n";
        exit(-1);
    }
    
    // initialize variables
    this->revertToStart();
}


VelPressureDep::~VelPressureDep()
{
    // does nothing
}


int VelPressureDep::setTrial(double normalForce, double velocity)
{
    trialN   = normalForce;
    trialVel = velocity;
    
    // determine COF at high velocities
    double muFast;
    if (trialN > 0.0)
        muFast = muFast0 - deltaMu*tanh(alpha*trialN/A);
    else
        muFast = muFast0;
    
    // get COF for given trial velocity
    double temp1 = exp(-transRate*fabs(trialVel));
    double temp2 = (muFast-muSlow)*temp1;
    mu = muFast - temp2;
    
    // get derivative of COF wrt normal force
    DmuDn = deltaMu*alpha/A/pow(cosh(alpha*trialN/A),2)*(temp1-1.0);
    
    // get derivative of COF wrt velocity
    if (trialVel != 0.0)
        DmuDvel = transRate*trialVel/fabs(trialVel)*temp2;
    else
        DmuDvel = 0.0;
    
    return 0;
}


double VelPressureDep::getFrictionForce()
{
    if (trialN > 0.0)
        return mu*trialN;
    else
        return 0.0;
}


double VelPressureDep::getFrictionCoeff()
{
    return mu;
}


double VelPressureDep::getDFFrcDNFrc()
{
    if (trialN >= 0.0)
        return mu + DmuDn*trialN;
    else
        return 0.0;
}


double VelPressureDep::getDFFrcDVel()
{
    if (trialN > 0.0)
        return DmuDvel*trialN;
    else
        return 0.0;
}


int VelPressureDep::commitState()
{
    return 0;
}


int VelPressureDep::revertToLastCommit()
{
    return 0;
}


int VelPressureDep::revertToStart()
{
    trialN   = 0.0;
    trialVel = 0.0;
    mu       = muSlow;
    DmuDn    = 0.0;
    DmuDvel  = 0.0;
    
    return 0;
}


FrictionModel* VelPressureDep::getCopy()
{
    VelPressureDep *theCopy = new VelPressureDep(this->getTag(),
        muSlow, muFast0, A, deltaMu, alpha, transRate);
    theCopy->trialN   = trialN;
    theCopy->trialVel = trialVel;
    theCopy->mu       = mu;
    theCopy->DmuDn    = DmuDn;
    theCopy->DmuDvel  = DmuDvel;
    
    return theCopy;
}


int VelPressureDep::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(7);
    data(0) = this->getTag();
    data(1) = muSlow;
    data(2) = muFast0;
    data(3) = A;
    data(4) = deltaMu;
    data(5) = alpha;
    data(6) = transRate;
    
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) 
        opserr << "VelPressureDep::sendSelf() - failed to send data.\n";
    
    return res;
}


int VelPressureDep::recvSelf(int cTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(7);
    
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0)  {
        opserr << "VelPressureDep::recvSelf() - failed to receive data.\n";
        this->setTag(0);      
        muSlow    = 0.0;
        muFast0   = 0.0;
        A         = 0.0;
        deltaMu   = 0.0;
        alpha     = 0.0;
        transRate = 0.0;
    }
    else  {
        this->setTag((int)data(0));
        muSlow    = data(1);
        muFast0   = data(2);
        A         = data(3);
        deltaMu   = data(4);
        alpha     = data(5);
        transRate = data(6);
    }
    
    // initialize variables
    this->revertToStart();
    
    return res;
}


void VelPressureDep::Print(OPS_Stream &s, int flag)
{
    s << "VelPressureDep tag: " << this->getTag() << endln;
    s << "  muSlow: " << muSlow << endln;
    s << "  muFast0: " << muFast0 << "  A: " << A << "  deltaMu: " << deltaMu;
    s << "  alpha: " << alpha << endln;
    s << "  transRate: " << transRate << endln;
}
