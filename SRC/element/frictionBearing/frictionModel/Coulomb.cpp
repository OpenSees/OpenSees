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
// Coulomb friction model.

#include <Coulomb.h>
#include <Channel.h>
#include <Information.h>
#include <elementAPI.h>

#include <math.h>


void * OPS_ADD_RUNTIME_VPV(OPS_Coulomb)
{
    // pointer to a friction model that will be returned
    FrictionModel *theFrnMdl = 0;
    
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING invalid number of arguments\n";
        opserr << "Want: frictionModel Coulomb tag mu\n";
        return 0;
    }
    
    int tag[1];
    double dData[1];
    int numData = 1;
    if (OPS_GetIntInput(&numData, tag) != 0) {
        opserr << "WARNING invalid tag for frictionModel Coulomb\n";
        return 0;
    }
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid data for frictionModel Coulomb " << tag[0] << endln;
        return 0;
    }
    
    // parsing was successful, allocate the friction model
    theFrnMdl = new Coulomb(tag[0], dData[0]);
    if (theFrnMdl == 0) {
        opserr << "WARNING could not create frictionModel of type Coulomb\n";
        return 0;
    }
    
    return theFrnMdl;
}


Coulomb::Coulomb()
    : FrictionModel(0, FRN_TAG_Coulomb),
    mu(0.0)
{
    // does nothing
}


Coulomb::Coulomb(int tag, double _mu)
    : FrictionModel(tag, FRN_TAG_Coulomb),
    mu(_mu)
{
    // check that COF is positive and not zero
    if (mu <= 0.0)  {
        opserr << "Coulomb::Coulomb - "
            << "the friction coefficient has to be positive.\n";
        exit(-1);
    }
    
    // initialize variables
    this->revertToStart();
}


Coulomb::~Coulomb()
{
    // does nothing
}


int Coulomb::setTrial(double normalForce, double velocity)
{	
    trialN   = normalForce;
    trialVel = velocity;
    
    return 0;
}


double Coulomb::getFrictionForce()
{
    if (trialN > 0.0)
        return mu*trialN;
    else
        return 0.0;
}


double Coulomb::getFrictionCoeff()
{
    return mu;
}


double Coulomb::getDFFrcDNFrc()
{
    if (trialN >= 0.0)
        return mu;
    else
        return 0.0;
}


double Coulomb::getDFFrcDVel()
{
    return 0.0;
}


int Coulomb::commitState()
{
    return 0;
}


int Coulomb::revertToLastCommit()
{
    return 0;
}


int Coulomb::revertToStart()
{
    trialN   = 0.0;
    trialVel = 0.0;
    
    return 0;
}


FrictionModel* Coulomb::getCopy()
{
    Coulomb *theCopy = new Coulomb(this->getTag(), mu);
    theCopy->trialN   = trialN;
    theCopy->trialVel = trialVel;
    
    return theCopy;
}


int Coulomb::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(2);
    data(0) = this->getTag();
    data(1) = mu;
    
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) 
        opserr << "Coulomb::sendSelf() - failed to send data.\n";
    
    return res;
}


int Coulomb::recvSelf(int cTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(2);
    
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0)  {
        opserr << "Coulomb::recvSelf() - failed to receive data.\n";
        this->setTag(0);      
        mu = 0.0;
    }
    else  {
        this->setTag((int)data(0));
        mu = data(1);
    }
    
    // initialize variables
    this->revertToStart();
    
    return res;
}


void Coulomb::Print(OPS_Stream &s, int flag)
{
    s << "Coulomb tag: " << this->getTag() << endln;
    s << "  mu: " << mu << endln;
}
