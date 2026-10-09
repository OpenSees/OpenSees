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

#ifndef FrictionSpringDamper_h
#define FrictionSpringDamper_h

// Written: Mark D. Denavit

#include <UniaxialMaterial.h>

class FrictionSpringDamper : public UniaxialMaterial
{
  public:
    FrictionSpringDamper(int tag, double Ke, double K, double Kp, double preload, double initStrain);
    FrictionSpringDamper();

    ~FrictionSpringDamper();

    const char *getClassType(void) const {return "FrictionSpringDamper";};

    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    UniaxialMaterial *getCopy(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
		 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);

  protected:

  private:
    // Material parameters (from input)
    double Ke;              // Elastic Modulus
    double K;               //
    double Kp;              // 
    double preload;         // 
    double initStrain;      //

    // Trial and last committed state variables
    double trialStrain;	    // current trial strain
    double trialStress;     // current trial stress
    double trialTangent;    // current trial tangent
    double commitStrain;    // last committed strain
    double commitStress;    // last committed stress
    double commitTangent;   // last committed tangent
};


#endif



