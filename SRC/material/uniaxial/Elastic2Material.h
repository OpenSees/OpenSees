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
                                                                        
// $Revision: 1.2 $
// $Date: 2007-06-15 20:22:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Elastic2Material.h,v $
                                                                        
                                                                        
#ifndef Elastic2Material_h
#define Elastic2Material_h

// Written: ZHY
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for 
// Elastic2Material. Elastic2Material provides the abstraction
// of an viscoelastic uniaxial material,
// i.e. stress = E*strain + eta*strainrate
//
//
// What: "@(#) Elastic2Material.h, revA"


#include <UniaxialMaterial.h>

class Elastic2Material : public UniaxialMaterial
{
  public:
    Elastic2Material(int tag, double E, double eta = 0.0);    
    Elastic2Material();    
    ~Elastic2Material();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0); 
    double getStrain(void) {return trialStrain;};
    double getStrainRate(void) {return trialStrainRate;};
    double getStress(void);
    double getTangent(void) {return (1-zeroE)*E;};
    double getDampTangent(void) {return eta;};
    double getInitialTangent(void) {return (1-zeroE)*E;};

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
    double trialStrain;
    double trialStrainRate;
    double E;
    double eta;
    double initialStrain;
    static int zeroE;
};


#endif

