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

#ifndef VelDepMultiLinear_h
#define VelDepMultiLinear_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 04/12
// Revision: A
//
// Description: This file contains the class definition for the VelDepMultiLinear
// friction model. In the velocity dependent multi-linear friction model the
// friction force is given in terms of a multi-linear curve that is defined by
// a number of non-negative velocity/friction pairs. If the normal force N is
// negative the friction force is zero.

#include "FrictionModel.h"

class VelDepMultiLinear : public FrictionModel
{
public:
    // constructor
    VelDepMultiLinear();
    VelDepMultiLinear(int tag,
        const Vector &velocityPoints,
        const Vector &frictionPoints);
    
    // destructor
    ~VelDepMultiLinear();
    
    const char *getClassType() const {return "VelDepMultiLinear";};
    
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
    Vector velocityPoints;  // velocity points on multi-linear curve
    Vector frictionPoints;  // friction points on multi-linear curve
    int trialID;            // trial ID into velocity, friction arrays
    int trialIDmin;         // minimum of trial ID
    int trialIDmax;         // maximum of trial ID
    int numDataPoints;      // number of data points defining curve
    
    double mu;              // current coefficient of friction (COF)
    double DmuDvel;         // derivative of COF wrt to velocity
};

#endif
