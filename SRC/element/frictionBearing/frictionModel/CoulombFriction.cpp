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
// $Source: /usr/local/cvs/OpenSees/SRC/element/frictionBearing/frictionModel/CoulombFriction.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class implementation for the
// CoulombFriction friction model.

#include <CoulombFriction.h>
#include <Channel.h>
#include <Information.h>

#include <math.h>


CoulombFriction::CoulombFriction()
    : FrictionModel(0, FRN_TAG_CoulombFriction),
    mu(0.0)
{
    // does nothing
}


CoulombFriction::CoulombFriction (int tag, double _mu)
    : FrictionModel(tag, FRN_TAG_CoulombFriction),
    mu(_mu)
{
    if (mu <= 0.0)  {
        opserr << "CoulombFriction::CoulombFriction - "
            << "the friction coefficient has to be positive\n";
        exit(-1);
    }
}


CoulombFriction::~CoulombFriction()
{
    // does nothing
}


int CoulombFriction::setTrial(double normalForce, double velocity)
{	
    trialN   = normalForce;
    trialVel = velocity;
    
    return 0;
}


double CoulombFriction::getFrictionForce(void)
{
    if (trialN > 0.0)
        return mu*trialN;
    else
        return 0.0;
}


double CoulombFriction::getFrictionCoeff(void)
{
    if (trialN > 0.0)
        return mu;
    else
        return 0.0;
}


double CoulombFriction::getDFFrcDNFrc(void)
{
    if (trialN > 0.0)
        return mu;
    else
        return 0.0;
}


int CoulombFriction::commitState(void)
{
    return 0;
}


int CoulombFriction::revertToLastCommit(void)
{
    return 0;
}


int CoulombFriction::revertToStart(void)
{
    trialN   = 0.0;
    trialVel = 0.0;
    
    return 0;
}


FrictionModel* CoulombFriction::getCopy(void)
{
    CoulombFriction *theCopy = new CoulombFriction(this->getTag(), mu);
    theCopy->trialN   = trialN;
    theCopy->trialVel = trialVel;
    
    return theCopy;
}


int CoulombFriction::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(2);
    data(0) = this->getTag();
    data(1) = mu;
    
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) 
        opserr << "CoulombFriction::sendSelf() - failed to send data\n";
    
    return res;
}


int CoulombFriction::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(2);
    
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0)  {
        opserr << "CoulombFriction::recvSelf() - failed to receive data\n";
        this->setTag(0);      
        mu = 0.0;
    }
    else  {
        this->setTag((int)data(0));
        mu = data(1);
    }
    
    return res;
}


void CoulombFriction::Print(OPS_Stream &s, int flag)
{
    s << "CoulombFriction tag: " << this->getTag() << endln;
    s << "  mu: " << mu << endln;
}
