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
// $Date: 2007-02-02 22:58:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ECC01.h,v $                                                                                                                                              
                                                                        
#ifndef ECC01_h
#define ECC01_h

// File: ECC01.h
//
// Written: Won Lee of Stanford University 
// Created: 09/04
// Revision: A
//
// Description: This file contains the class definition for 
// ECC01.h
//   - ECC model based on Han et al. model
//      (Han TS, Feenstra PH, Billington SL, ACI Structural Journal,
//			Nov-Dec 2003, "Simulation of Highly Ductile Fiber Reinforced
//			Cement-Based Composite Components Under Cyclic Loading")
//

#include <UniaxialMaterial.h>

class ECC01 : public UniaxialMaterial
{
 public:
  ECC01 (int tag, double SIGT0, double EPST0, double SIGT1, double EPST1, double EPST2, double SIGC0, 
	 double EPSC0, double EPSC1, double ALPHAT1, double ALPHAT2, double ALPHAC, double ALPHACU, double BETAT, double BETAC);
  ECC01 ();
  ~ECC01();

  int setTrialStrain(double strain, double strainRate = 0.0); 
  int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
  double getStrain(void);      
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) {return sigc0/epsc0;}
  
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
  /*** Material Properties ***/
  double sigt0;		// Tensile cracking stress
  double epst0;		// Strain at tensile cracking
  double sigt1;		// Peak tensile stress
  double epst1;		// Peak tensile strain
  double epst2;		// Ultimate tensile strain
  double sigc0;		// Peak compressive stress
  double epsc0;		// Peak compressive strain
  double epsc1;		// Ultimate compressive strain
  double alphaT1;	// Constant parameter for unloading equation in tensile strain-hardening region
  double alphaT2;	// Constant parameter for unloading equation in tensile softening region (=1 for linear unloading)
  double alphaC;	// Constant parameter for unloading equation in compressive softening region
  double alphaCU;	// Constant parameter for envelope compression softening equation
  double betaT;		// Constant parameter for permanent strain in tension
  double betaC;		// Constant parameter for permanent strain in compression
  
  /*** CONVERGED History Variables ***/
  double CminStrain;   // Smallest (most negative) previous concrete strain (compression)
  double CmaxStrain;   // Largest previous conrete strain (tension)
  double Cstmp;			// temporary stress value, used to compute stresses and strains in re/unloading
  double Cetmp;			// temporary strain value, used to compute stresses and strains in re/unloading
  int Cindex;			// Index that tells you where you are on the stress-strain curve
  
  /*** CONVERGED State Variables ***/
  double Cstrain;
  double Cstress;   
  double Ctangent;	// Don't need Ctangent other than for revert and sendSelf/recvSelf
  // Storing it is better than recomputing it!!!
  
  /*** TRIAL History Variables ***/
  double TminStrain;
  double TmaxStrain;   
  double Tstmp; 
  double Tetmp; 
  int Tindex;
  
  /*** TRIAL State Variables ***/
  double Tstrain;
  double Tstress;
  double Ttangent; // Not really a state variable, but declared here
  // for convenience
  
  //void determineTrialState (double dStrain);
  
  void envelope();
  void ECCGetStressAndStiffness(int index, double sigmax, double epstul, double sigmin, double epscul);
};


#endif



