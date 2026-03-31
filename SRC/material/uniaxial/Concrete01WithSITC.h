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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete01WithSITC.h,v $
                                                                        
#ifndef Concrete01WithSITC_h
#define Concrete01WithSITC_h

// Modified by: Won Lee
// Created: 10/3/05
// Modified from Concrete01.h (see details below)
// Description: This file contains the class definition for Concrete01WithSITC
// Description: Concrete01 model modified to include SITC effect (ref. Prof. 
//	John Stanton of Univ. of Washington).  Use modified rules from his paper to include this 
//	effect (J.F. Stanton and H.D. McNiven, "The Development of a Mathematical
//	Model to Predict the Flexural Response of Reinforced Concrete Beams to Cyclic
//	Loads, Using System Identification", EERC Report Number 79/02, January 1979.
//

// BASED ON FILE:
// File: Concrete01.h
//
// Written: MHS 
// Created: 06/99
// Revision: A
//
// Description: This file contains the class definition for 
// Concrete01.h adapted from Concr1.f90 (Filippou)
//   - Modified Kent-Park envelope
//   - No tension
//   - Linear unloading/reloading
//
// What: "@(#) Concrete01.h, revA"


#include <UniaxialMaterial.h>

class Concrete01WithSITC : public UniaxialMaterial
{
 public:
  Concrete01WithSITC (int tag, double fpc, double eco, double fpcu, double ecu, double endStrainSITC = 0.03);
  Concrete01WithSITC ();
  ~Concrete01WithSITC();

  const char *getClassType(void) const {return "Concrete01WithSITC";}

      int setTrialStrain(double strain, double strainRate = 0.0); 
      int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
      double getStrain(void);      
      double getStress(void);
      double getTangent(void);
      double getInitialTangent(void) {return 2.0*fpc/epsc0;}

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
      double fpc;    // Compressive strength
      double epsc0;  // Strain at compressive strength
      double fpcu;   // Crushing strength
      double epscu;  // Strain at crushing strength

      /*** CONVERGED History Variables ***/
      double CminStrain;   // Smallest previous concrete strain (compression)
      double CunloadSlope; // Unloading (reloading) slope from CminStrain
      double CendStrain;   // Strain at the end of unloading from CminStrain
      double CmaxStrain;   // Largest previous concrete strain (tension) 
      double CslopeSITC; 
      double CendStrainSITC; 
      int Cindex; 
      int CsmallStrainIndex;
      
      
      /*** CONVERGED State Variables ***/
      double Cstrain;
      double Cstress;   
      double Ctangent;	// Don't need Ctangent other than for revert and sendSelf/recvSelf
      // Storing it is better than recomputing it!!!

      /*** TRIAL History Variables ***/
      double TminStrain;
      double TunloadSlope;
      double TendStrain;
      double TmaxStrain; 
      double TslopeSITC; 

      int Tindex;
      int TsmallStrainIndex;
      
      /*** TRIAL State Variables ***/
      double Tstrain;
      double Tstress;
      double Ttangent; // Not really a state variable, but declared here
                       // for convenience
      
      void determineTrialState (double dStrain);
      
      void reload();
      void unload();
      void envelope();
      void getSITCslope();
};


#endif


