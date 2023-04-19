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
// Date: 03/2019
// Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Hertzdamp.h
//
// Written: Patrick J. Hughes, University of California - San Diego
// Created: 03/2019
//
// Description: This file contains the class definition for the
// Hertzdamp uniaxialMaterial.
//
// References:
// Muthukumar, S., and DesRoches, R. (2006). "A Hertz Contact Model with Non-linear Damping for Pounding Simulation." 
//   Earthquake Engineering and Structural Dynamics, 35, 811-828.
// Ye, Kun., Li, L., and Zhu, H. (2008) "A Note on the Hertz Contact Model with Nonlinear Damping for Pounding Simulation."
//	 Earthquake Engineering and Structural Dynamics, 38, 1135-1142.
//
// Variables:
// Kh: nonlinear Hertz contact stiffness
// xiNorm: normalized damping coefficient
// gap: initial gap (must be a negative value)
// n: displacement exponent (default is 1.5)


#ifndef Hertzdamp_h
#define Hertzdamp_h

#include <UniaxialMaterial.h>

class Hertzdamp : public UniaxialMaterial
{
  public: 
    Hertzdamp(int tag, double Kh, double xiNorm, double gap, double n);   
    Hertzdamp(); 
    ~Hertzdamp();

    const char *getClassType(void) const {return "Hertzdamp";};

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
    double xiNorm; // normalized damping coefficient
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
	double xi; // impact damping coefficient (depends on material properties and impact velocity)
	int printFlag; // print flag for impact event
	
};

#endif