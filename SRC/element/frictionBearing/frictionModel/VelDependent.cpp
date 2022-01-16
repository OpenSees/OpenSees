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
// VelDependent friction model.

#include <VelDependent.h>
#include <Channel.h>
#include <Information.h>
#include <elementAPI.h>

#include <math.h>


void * OPS_ADD_RUNTIME_VPV(OPS_VelDependent)
{
    // pointer to a friction model that will be returned
    FrictionModel *theFrnMdl = 0;
    
    if (OPS_GetNumRemainingInputArgs() < 4) {
        opserr << "WARNING invalid number of arguments\n";
        opserr << "Want: frictionModel VelDependent tag muSlow muFast transRate\n";
        return 0;
    }
    
    int tag[1];
    double dData[3];
    int numData = 1;
    if (OPS_GetIntInput(&numData, tag) != 0) {
        opserr << "WARNING invalid tag for frictionModel VelDependent\n";
        return 0;
    }
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid data for frictionModel VelDependent " << tag[0] << endln;
        return 0;
    }
    
    // parsing was successful, allocate the friction model
    theFrnMdl = new VelDependent(tag[0], dData[0], dData[1], dData[2]);
    if (theFrnMdl == 0) {
        opserr << "WARNING could not create frictionModel of type VelDependent\n";
        return 0;
    }
    
    return theFrnMdl;
}


VelDependent::VelDependent()
    : FrictionModel(0, FRN_TAG_VelDependent),
    muSlow(0.0), muFast(0.0), transRate(0.0),
    mu(0.0), DmuDvel(0.0)
{
    // does nothing
}


VelDependent::VelDependent(int tag, double muslow,
    double mufast, double transrate)
    : FrictionModel(tag, FRN_TAG_VelDependent),
    muSlow(muslow), muFast(mufast), transRate(transrate),
    mu(0.0), DmuDvel(0.0)
{
    // check that COF are positive and not zero
    if (muSlow <= 0.0  || muFast <= 0.0)  {
        opserr << "VelDependent::VelDependent - "
            << "the friction coefficients have to be positive.\n";
        exit(-1);
    }
    // check that transRate is positive
    if (transRate < 0.0)  {
        opserr << "VelDependent::VelDependent - "
            << "the transition rate has to be positive.\n";
        exit(-1);
    }
    
    // initialize variables
    this->revertToStart();
}


VelDependent::~VelDependent()
{
    // does nothing
}


int VelDependent::setTrial(double normalForce, double velocity)
{
    trialN   = normalForce;
    trialVel = velocity;
    
    // get COF for given trial velocity
    double temp = (muFast-muSlow)*exp(-transRate*fabs(trialVel));
    mu = muFast - temp;
    
    // get derivative of COF wrt velocity
    if (trialVel != 0.0)
        DmuDvel = transRate*trialVel/fabs(trialVel)*temp;
    else
        DmuDvel = 0.0;
    
    return 0;
}


double VelDependent::getFrictionForce()
{
    if (trialN > 0.0)
        return mu*trialN;
    else
        return 0.0;
}


double VelDependent::getFrictionCoeff()
{
    return mu;
}


double VelDependent::getDFFrcDNFrc()
{
    if (trialN >= 0.0)
        return mu;
    else
        return 0.0;
}


double VelDependent::getDFFrcDVel()
{
    if (trialN > 0.0)
        return DmuDvel*trialN;
    else
        return 0.0;
}


int VelDependent::commitState()
{
    return 0;
}


int VelDependent::revertToLastCommit()
{
    return 0;
}


int VelDependent::revertToStart()
{
    trialN   = 0.0;
    trialVel = 0.0;
    mu       = muSlow;
    DmuDvel  = 0.0;
    
    return 0;
}


FrictionModel* VelDependent::getCopy()
{
    VelDependent *theCopy = new VelDependent(this->getTag(),
        muSlow, muFast, transRate);
    theCopy->trialN   = trialN;
    theCopy->trialVel = trialVel;
    theCopy->mu       = mu;
    theCopy->DmuDvel  = DmuDvel;
    
    return theCopy;
}


int VelDependent::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(4);
    data(0) = this->getTag();
    data(1) = muSlow;
    data(2) = muFast;
    data(3) = transRate;
    
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) 
        opserr << "VelDependent::sendSelf() - failed to send data.\n";
    
    return res;
}


int VelDependent::recvSelf(int cTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(4);
    
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0)  {
        opserr << "VelDependent::recvSelf() - failed to receive data.\n";
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
    
    // initialize variables
    this->revertToStart();
    
    return res;
}


void VelDependent::Print(OPS_Stream &s, int flag)
{
    s << "VelDependent tag: " << this->getTag() << endln;
    s << "  muSlow: " << muSlow << endln;
    s << "  muFast: " << muFast << endln;
    s << "  transRate: " << transRate << endln;
}
