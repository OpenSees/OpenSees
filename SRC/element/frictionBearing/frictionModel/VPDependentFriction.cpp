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
// $Source: /usr/local/cvs/OpenSees/SRC/element/frictionBearing/frictionModel/VPDependentFriction.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class implementation for the
// VPDependentFriction friction model.

#include <VPDependentFriction.h>
#include <Channel.h>
#include <Information.h>

#include <math.h>


VPDependentFriction::VPDependentFriction()
    : FrictionModel(0, FRN_TAG_VPDependentFriction),
    muSlow(0.0), muFast0(0.0), A(0.0), deltaMu(0.0),
    alpha(0.0), transRate(0.0), mu(0.0)
{
    // does nothing
}


VPDependentFriction::VPDependentFriction (int tag, double muslow,
    double mufast0, double a, double deltamu, double _alpha, double transrate)
    : FrictionModel(tag, FRN_TAG_VPDependentFriction),
    muSlow(muslow), muFast0(mufast0), A(a), deltaMu(deltamu),
    alpha(_alpha), transRate(transrate), mu(0.0)
{
    if (muSlow <= 0.0  || muFast0 <= 0.0)  {
        opserr << "VPDependentFriction::VPDependentFriction - "
            << "the friction coefficients have to be positive\n";
        exit(-1);
    }
}


VPDependentFriction::~VPDependentFriction()
{
    // does nothing
}


int VPDependentFriction::setTrial(double normalForce, double velocity)
{	
    trialN   = normalForce;
    trialVel = velocity;
    
    double muFast = muFast0 - deltaMu*tanh(alpha*trialN/A);
    mu = muFast - (muFast-muSlow)*exp(-transRate*fabs(trialVel));
    
    return 0;
}


double VPDependentFriction::getFrictionForce(void)
{    
    if (trialN > 0.0)
        return mu*trialN;
    else
        return 0.0;
}


double VPDependentFriction::getFrictionCoeff(void)
{
    if (trialN > 0.0)
        return mu;
    else
        return 0.0;
}


double VPDependentFriction::getDFFrcDNFrc(void)
{
    double dFFdFN = mu + deltaMu*alpha*trialN/A/
        pow(cosh(alpha*trialN/A),2)*
        (exp(-transRate*fabs(trialVel))-1.0);
    
    if (trialN > 0.0) 
        return dFFdFN;
    else
        return 0.0;
}


int VPDependentFriction::commitState(void)
{
    return 0;
}


int VPDependentFriction::revertToLastCommit(void)
{
    return 0;
}


int VPDependentFriction::revertToStart(void)
{
    trialN   = 0.0;
    trialVel = 0.0;
    mu       = 0.0;
    
    return 0;
}


FrictionModel* VPDependentFriction::getCopy(void)
{
    VPDependentFriction *theCopy = new VPDependentFriction(this->getTag(),
        muSlow, muFast0, A, deltaMu, alpha, transRate);
    theCopy->trialN   = trialN;
    theCopy->trialVel = trialVel;
    theCopy->mu = mu;
    
    return theCopy;
}


int VPDependentFriction::sendSelf(int cTag, Channel &theChannel)
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
        opserr << "VPDependentFriction::sendSelf() - failed to send data\n";
    
    return res;
}


int VPDependentFriction::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(7);
    
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0)  {
        opserr << "VPDependentFriction::recvSelf() - failed to receive data\n";
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
    
    return res;
}


void VPDependentFriction::Print(OPS_Stream &s, int flag)
{
    s << "VPDependentFriction tag: " << this->getTag() << endln;
    s << "  muSlow: " << muSlow << endln;
    s << "  muFast0: " << muFast0 << "  A: " << A << "  deltaMu: " << deltaMu;
    s << "  alpha: " << alpha << endln;
    s << "  transRate: " << transRate << endln;
}
