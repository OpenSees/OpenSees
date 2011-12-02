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
                                                                        
// $Revision: 1.1 $
// $Date: 2010-08-17 00:23:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticBilin.h,v $
                                                                        
#ifndef ElasticBilin_h
#define ElasticBilin_h

// Written: fmk 
// Created: 07/98
//
// Description: This file contains the class definition for 
// ElasticBilin. ElasticBilin provides the abstraction
// of an elastic bilinear material uniaxial material, 
//
// What: "@(#) ElasticBilin.h, revA"

#include <UniaxialMaterial.h>

class ElasticBilin : public UniaxialMaterial
{
  public:
    ElasticBilin(int tag, double E1, double E2, double eps2);    
    ElasticBilin(int tag, double E1P, double E2P, double epsP, double E1N, double E2N, double eps2N);    
    ElasticBilin();    

    ~ElasticBilin();

    const char *getClassType(void) const {return "ElasticBilin";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void) {return E1P;};

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
    double E1P, E1N, E2P, E2N;   // elastic modulus
    double eps2P;	       // strain at which E2P takes place	
    double eps2N;	       // strain at which E2P takes place	

    double trialStrain, trialStress, trialTangent, commitStrain; 
};


#endif



