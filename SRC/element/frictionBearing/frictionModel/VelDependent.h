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

#ifndef VelDependent_h
#define VelDependent_h         

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class definition for the VelDependent
// friction model after Constantinou et al. (1990). In the velocity dependent
// friction model the friction force is given in terms of the friction coefficients
// at low and high velocities and a constant describing the rate of transition.
// If the normal force N is negative the friction force is zero.

#include "FrictionModel.h"

class VelDependent : public FrictionModel
{
public:
    // constructor
    VelDependent();
    VelDependent(int tag, double muSlow, double muFast, double transRate);
    
    // destructor
    ~VelDependent();
    
    const char *getClassType() const {return "VelDependent";};
    
    // public methods to set and obtain response
    int setTrial(double normalForce, double velocity = 0.0);
    double getFrictionForce();
    double getFrictionCoeff();
    double getDFFrcDNFrc();
    double getDFFrcDVel();
    
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    
    FrictionModel *getCopy();
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
protected:

private:
    double muSlow;     // coefficient of friction at low velocity
    double muFast;     // coefficient of friction at high velocity
    double transRate;  // transition rate from low to high velocity
    
    double mu;         // current coefficient of friction (COF)
    double DmuDvel;    // derivative of COF wrt to velocity
};

#endif
