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
// $Source: /usr/local/cvs/OpenSees/SRC/element/frictionBearing/frictionModel/VDependentFriction.h,v $

#ifndef VDependentFriction_h
#define VDependentFriction_h         

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class definition for the VDependentFriction
// friction model after Constantinou et al. (1990). In the velocity dependent
// friction model the friction force is given in terms of the friction coefficients
// at low and high velocities and a constant describing the rate of transition.
//
// What: "@(#) VDependentFriction.h, revA"

#include <FrictionModel.h>

class VDependentFriction : public FrictionModel
{
public:
    // constructor
    VDependentFriction();
    VDependentFriction(int tag, double muSlow, double muFast, double transRate);
    
    // destructor
    ~VDependentFriction();
    
    // public methods to set and obtain response
    int setTrial(double normalForce, double velocity = 0.0);
    double getFrictionForce(void);
    double getFrictionCoeff(void);
    double getDFFrcDNFrc(void);
    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    
    FrictionModel *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
protected:

private:
    double muSlow;      // coefficient of friction at low velocity
    double muFast;      // coefficient of friction at high velocity
    double transRate;   // transition rate from low to high velocity
    
    double mu;          // current coefficient of friction
};

#endif
