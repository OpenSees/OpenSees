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
// $Date: 2006-08-03 23:42:19 $: 2001/07/16 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/CableMaterial.h,v $
                                                                        
                                                                        
#ifndef CableMaterial_h
#define CableMaterial_h

// Written: Charles Chadwell
// Created: 07/01
//
// Description: This file contains the class definition for 
// CableMaterial. CableMaterial provides the abstraction
// of an elastic uniaxial material,
//
// The input parameters are the Prestress, E, Effective Self Weight (gravity component of 
// Weight per volume transverse to the cable), and Length of the cable.
//
// The cable Force Displacement is converted to Stress Strain material for use 
// with the truss element.  The stress strain ranges from slack (large strain at zero 
// stress) to taught (linear with modulus E).  The material has no history and is not
// path dependent.
//
//
// What: "@(#) CableMaterial.h, revA"


#include <UniaxialMaterial.h>

class CableMaterial : public UniaxialMaterial
{
  public:
    CableMaterial(int tag, double Prestress, double E, double unitWeightEff, double L_Element);    
    CableMaterial();    
    ~CableMaterial();

    const char *getClassType(void) const {return "CableMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void) {return trialStrain;};
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return 1.0e-8;}; 

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
    //int setParameter(const char **argv, int argc, Information &info);
    //int updateParameter(int parameterID, Information &info);

  protected:
    
  private:

    double Ps;
    double E;
    double Mue;
    double L;
    double trialStrain;
    
    double evalStress(double stress);
    double abs(double value);

    double trialStress;      // current trial stress
    double trialTangent;     // current trial tangent
};


#endif

