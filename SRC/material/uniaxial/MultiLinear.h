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
                                                                        
#ifndef MultiLinear_h
#define MultiLinear_h

// Written: fmk 
// Created: 05/12
// Revision: A
//
// Description: This file contains the class definition for 
// MultiLinear. MultiLinear provides the abstraction
// of an elastic perfectly plastic uniaxial material, 
//
// What: "@(#) MultiLinear.h, revA"

#include <UniaxialMaterial.h>
#include <Matrix.h>

class MultiLinear : public UniaxialMaterial
{
  public:
    MultiLinear(int tag, const Vector &s, const Vector &e);
    MultiLinear();
    virtual ~MultiLinear();

    const char *getClassType(void) const {return "MultiLinear";};

    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void) {return data(0,4);};

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    UniaxialMaterial *getCopy(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
           FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
  
 protected:
  
 private:
    Matrix data;
    int numSlope;
    int tSlope;

    double tStrain;     // current t strain
    double tStress;     // current t stress
    double tTangent;    // current t tangent
    double cStrain;     // last ced strain
    double cStress;     // last ced stress
    double cTangent;    // last cted  tangent
};

#endif
