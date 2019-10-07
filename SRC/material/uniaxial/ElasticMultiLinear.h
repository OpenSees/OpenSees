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

#ifndef ElasticMultiLinear_h
#define ElasticMultiLinear_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 12/11
// Revision: A
//
// Description: This file contains the class definition for 
// ElasticMultiLinear. ElasticMultiLinear provides the abstraction
// of an elastic multilinear uniaxial material.

#include <UniaxialMaterial.h>

class ElasticMultiLinear : public UniaxialMaterial
{
public:
    // constructor
    ElasticMultiLinear(int tag,
		       const Vector &strainPoints, 
		       const Vector &stressPoints,
               double eta = 0.0);    
    ElasticMultiLinear();    

    // destructor
    ~ElasticMultiLinear();

    const char *getClassType() const {return "ElasticMultiLinear";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain() {return trialStrain;};
    double getStrainRate() {return trialStrainRate;};
    double getStress() {return trialStress;};
    double getTangent() {return trialTangent;};
    double getInitialTangent() {return initTangent;};
    double getDampTangent() {return eta;};

    int commitState();
    int revertToLastCommit();
    int revertToStart();

    UniaxialMaterial *getCopy();

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
        FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

protected:

private:
    Vector strainPoints;     // strain points on multi-linear curve
    Vector stressPoints;     // stress points on multi-linear curve
    double eta;              // damping tangent modulus
    int trialID;             // trial ID into strain, stress arrays
    int trialIDmin;          // minimum of trial ID
    int trialIDmax;          // maximum of trial ID
    int numDataPoints;       // number of data points defining curve
    double initTangent;      // initial tangent (on positive side)
    double trialStrain;      // trial strain
    double trialStrainRate;  // trial strain rate
    double trialStress;      // trial stress
    double trialTangent;     // trial tangent
};

#endif
