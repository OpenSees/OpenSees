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
                                                                        
// $Revision: 1.9 $
// $Date: 2010-04-06 20:18:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ENTMaterial.h,v $
                                                                        
                                                                        
#ifndef ENTMaterial_h
#define ENTMaterial_h

// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for 
// ENTMaterial. ENTMaterial provides the abstraction
// of an elastic  uniaxial material under compression
// i.e. stress = E*strain e
// under tension however it exhbits the following
//     stress = a*(tanh(strain*b))
//     tangent = a*(1-tanh(strain*b)*tanh(strain*b));
//
//
// What: "@(#) ENTMaterial.h, revA"

#include <UniaxialMaterial.h>

class ENTMaterial : public UniaxialMaterial
{
  public:
  ENTMaterial(int tag, double E, double a = 0.0, double b = 1.0);    
    ENTMaterial();    
    ~ENTMaterial();

    const char *getClassType(void) const {return "ENTMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return E;};

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

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int activateParameter(int parameterID);
    double getStressSensitivity(int gradIndex, bool conditional);
    double getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(double strainGradient, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    double E;
    double trialStrain;

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////

    double a;
    double b;
};


#endif

