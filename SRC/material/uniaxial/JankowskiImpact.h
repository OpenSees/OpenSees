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
//                                                                     
// Revision: 1.0
// Date: 05/2019
// Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/JankowskiImpact.h
//
// Written: Patrick J. Hughes, University of California - San Diego
// Created: 05/2019
//
// Description: This file contains the class definition for the
// JankowskiImpact uniaxialMaterial.
//
// References:
// Jankowski, R. (2005). "Non-linear Viscoelastic Modelling of Earthquake-Induced Structural Pounding." 
//   Earthquake Engineering and Structural Dynamics, 34, 595-611.
// Jankowski, R. (2006). "Analytical Expression Between the Impact Damping Ratio and the Coefficient of Restitution in the Non-linear Viscoelastic Modelling of Structural Pounding." 
//   Earthquake Engineering and Structural Dynamics, 35, 517-524.
//
// Variables:
// Kh: nonlinear Hertz contact stiffness
// xi: impact damping ratio
// mEff: effective mass of the colliding bodies
// gap: initial gap (must be a negative value)
// n: displacement exponent (default is 1.5)


#ifndef JankowskiImpact_h
#define JankowskiImpact_h

#include <UniaxialMaterial.h>

class JankowskiImpact : public UniaxialMaterial
{
  public: 
    JankowskiImpact(int tag, double Kh, double xi, double mEff, double gap, double n);   
    JankowskiImpact(); 
    ~JankowskiImpact();

    const char *getClassType(void) const {return "JankowskiImpact";};

    int setTrialStrain(double strain, double strainRate); 
	
    double getStrain(void); 
    double getStrainRate(void);
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
    // input variables
    double Kh; // nonlinear Hertz contact stiffness
    double xi; // impact damping ratio
    double mEff; //effective mass of the colliding bodies
    double gap; // initial gap distance (must be a negative value)
    double n; // displacement exponent (default is 1.5)
    
    // state variables
	double commitStrain;
	double commitStrainRate;
	double commitStress;
	double commitTangent;
	double trialStrain;
	double trialStrainRate;
	double trialStress;
	double trialTangent;

	// other variables
	int printFlag; // print flag for impact event
	
};

#endif