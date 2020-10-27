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

#ifndef VelPressureDep_h
#define VelPressureDep_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class definition for the VelPressureDep
// friction model after Constantinou et al. (1996). In the velocity and pressure
// dependent friction model the friction force is given in terms of the friction
// coefficients at low and high velocities with the latter one being a function of
// pressure. If the normal force N is negative the friction force is zero.

#include "FrictionModel.h"

class VelPressureDep : public FrictionModel
{
public:
    // constructor
    VelPressureDep();
    VelPressureDep(int tag, double muSlow, double muFast0, double A,
        double deltaMu, double alpha, double transRate);
    
    // destructor
    ~VelPressureDep();
    
    const char *getClassType() const {return "VelPressureDep";};
    
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
    double muFast0;    // initial coefficient of friction at high velocity
    double A;          // nominal contact area
    double deltaMu;    // pressure parameter
    double alpha;      // pressure parameter
    double transRate;  // transition rate from low to high velocity
    
    double mu;         // current coefficient of friction (COF)
    double DmuDn;      // derivative of COF wrt to normal force
    double DmuDvel;    // derivative of COF wrt to velocity
};

#endif
