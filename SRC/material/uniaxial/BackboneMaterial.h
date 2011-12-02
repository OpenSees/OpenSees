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
// $Date: 2008-11-04 22:20:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BackboneMaterial.h,v $
                                                      
// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the class definition for 
// BackboneMaterial.  BackboneMaterial uses a HystereticBackbone
// object to represent a path-independent uniaxial material.  Since
// it is path-independent, no state information is stored by
// BackboneMaterial.

#ifndef BackboneMaterial_h
#define BackboneMaterial_h

#include <UniaxialMaterial.h>

class HystereticBackbone;

class BackboneMaterial : public UniaxialMaterial
{
  public:
    BackboneMaterial(int tag, HystereticBackbone &backbone); 
    BackboneMaterial();
    ~BackboneMaterial();

    const char *getClassType(void) const {return "BackboneMaterial";}; 

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStrainRate(void);
    double getStress(void);
    double getTangent(void);
    double getDampTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
  private:
    HystereticBackbone *theBackbone;
    double strain;
};


#endif

