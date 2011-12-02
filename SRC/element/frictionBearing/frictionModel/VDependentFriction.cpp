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

// $Revision: 1.1 $
// $Date: 2009-04-17 23:02:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/frictionBearing/frictionModel/VDependentFriction.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class implementation for the
// VDependentFriction friction model.

#include <VDependentFriction.h>
#include <Channel.h>
#include <Information.h>

#include <math.h>


VDependentFriction::VDependentFriction()
    : FrictionModel(0, FRN_TAG_VDependentFriction),
    muSlow(0.0), muFast(0.0), transRate(0.0), mu(0.0)
{
    // does nothing
}


VDependentFriction::VDependentFriction (int tag, double muslow,
    double mufast, double transrate)
    : FrictionModel(tag, FRN_TAG_VDependentFriction),
    muSlow(muslow), muFast(mufast), transRate(transrate), mu(0.0)
{
    if (muSlow <= 0.0  || muFast <= 0.0)  {
        opserr << "VDependentFriction::VDependentFriction - "
            << "the friction coefficients have to be positive\n";
        exit(-1);
    }
}


VDependentFriction::~VDependentFriction()
{
    // does nothing
}


int VDependentFriction::setTrial(double normalForce, double velocity)
{	
    trialN   = normalForce;
    trialVel = velocity;
    
    mu = muFast - (muFast-muSlow)*exp(-transRate*fabs(trialVel));
    
    return 0;
}


double VDependentFriction::getFrictionForce(void)
{
    if (trialN > 0.0)
        return mu*trialN;
    else
        return 0.0;
}


double VDependentFriction::getFrictionCoeff(void)
{
    if (trialN > 0.0)
        return mu;
    else
        return 0.0;
}


double VDependentFriction::getDFFrcDNFrc(void)
{
    if (trialN > 0.0)
        return mu;
    else
        return 0.0;
}


int VDependentFriction::commitState(void)
{
    return 0;
}


int VDependentFriction::revertToLastCommit(void)
{
    return 0;
}


int VDependentFriction::revertToStart(void)
{
    trialN   = 0.0;
    trialVel = 0.0;
    mu       = 0.0;
    
    return 0;
}


FrictionModel* VDependentFriction::getCopy(void)
{
    VDependentFriction *theCopy = new VDependentFriction(this->getTag(),
        muSlow, muFast, transRate);
    theCopy->trialN   = trialN;
    theCopy->trialVel = trialVel;
    theCopy->mu = mu;
    
    return theCopy;
}


int VDependentFriction::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(4);
    data(0) = this->getTag();
    data(1) = muSlow;
    data(2) = muFast;
    data(3) = transRate;
    
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) 
        opserr << "VDependentFriction::sendSelf() - failed to send data\n";
    
    return res;
}


int VDependentFriction::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(4);
    
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0)  {
        opserr << "VDependentFriction::recvSelf() - failed to receive data\n";
        this->setTag(0);      
        muSlow    = 0.0;
        muFast    = 0.0;
        transRate = 0.0;
    }
    else  {
        this->setTag((int)data(0));
        muSlow    = data(1);
        muFast    = data(2);
        transRate = data(3);
    }
    
    return res;
}


void VDependentFriction::Print(OPS_Stream &s, int flag)
{
    s << "VDependentFriction tag: " << this->getTag() << endln;
    s << "  muSlow: " << muSlow << endln;
    s << "  muFast: " << muFast << endln;
    s << "  transRate: " << transRate << endln;
}
