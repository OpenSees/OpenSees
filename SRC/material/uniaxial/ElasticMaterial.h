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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-02-25 23:33:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticMaterial.h,v $
                                                                        
                                                                        
#ifndef ElasticMaterial_h
#define ElasticMaterial_h

// File: ~/material/ElasticMaterial.h
//
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for 
// ElasticMaterial. ElasticMaterial provides the abstraction
// of an viscoelastic uniaxial material,
// i.e. stress = E*strain + eta*strainrate
//
//
// What: "@(#) ElasticMaterial.h, revA"


#include <UniaxialMaterial.h>

class ElasticMaterial : public UniaxialMaterial
{
  public:
    ElasticMaterial(int tag, double E, double eta = 0.0);    
    ElasticMaterial();    
    ~ElasticMaterial();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0); 
    double getStrain(void) {return trialStrain;};
    double getStrainRate(void) {return trialStrainRate;};
    double getStress(void);
    double getTangent(void) {return E;};
    double getDampTangent(void) {return eta;};
    double getInitialTangent(void) {return E;};

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
    int setParameter(const char **argv, int argc, Information &info);
    int updateParameter(int parameterID, Information &info);

  protected:
    
  private:
    double trialStrain;
    double trialStrainRate;
    double E;
    double eta;
};


#endif

