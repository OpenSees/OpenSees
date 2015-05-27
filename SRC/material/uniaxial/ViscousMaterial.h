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
                                                                        
// $Revision: 1.6 $
// $Date: 2008-10-17 23:35:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ViscousMaterial.h,v $
                                                                        
// Written: Mehrdad Sasani 
// Created: June 2000
// Revision: A
//
// Description: This file contains the class interface for 
// ViscousMaterial.  ViscousMaterial defines a force(F)-velocity(v)
// relationship of the form F = C*pow(v,a), where C is a prescribed
// constant and a is a real number.

#ifndef ViscousMaterial_h
#define ViscousMaterial_h

#include <UniaxialMaterial.h>

class ViscousMaterial : public UniaxialMaterial
{
  public:
    ViscousMaterial(int tag, double C, double Alpha, double minVel = 1.0e-11);    
    ViscousMaterial(); 
    ~ViscousMaterial();

    const char *getClassType(void) const {return "ViscousMaterial";};

    int setTrialStrain(double velocity, double strainRate = 0.0); 
    double getStrain(void); 
    double getStrainRate(void);
    double getStress(void);

    double getTangent(void);
    double getInitialTangent(void);
    double getDampTangent(void);


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
    double trialRate;
    double C;
    double Alpha;
    double minVel;
    double commitStrain;
    double commitRate;
};


#endif

