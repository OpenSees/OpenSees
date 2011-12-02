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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete01.h,v $
                                                                        
                                                                        
#ifndef Concrete01_h
#define Concrete01_h

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

class Concrete01 : public UniaxialMaterial
{
   public:
      Concrete01 (int tag, double fpc, double eco, double fpcu, double ecu,
                  double epsmin=NEG_INF_STRAIN, double epsmax=POS_INF_STRAIN);
		Concrete01 ();
      ~Concrete01();

      int setTrialStrain(double strain, double strainRate = 0.0); 
      double getStrain(void);      
      double getStress(void);
      double getTangent(void);
      double getSecant (void);

      int commitState(void);
      int revertToLastCommit(void);    
      int revertToStart(void);        

      UniaxialMaterial *getCopy(void);
    
      int sendSelf(int commitTag, Channel &theChannel);  
      int recvSelf(int commitTag, Channel &theChannel, 
	  	 FEM_ObjectBroker &theBroker);    
    
      void Print(ostream &s, int flag =0);

   protected:

   private:
      /*** Material Properties ***/
      double fpc;    // Compressive strength
      double epsc0;  // Strain at compressive strength
	  double Ec0;
      double fpcu;   // Crushing strength
      double epscu;  // Strain at crushing strength
      double epsmin; // Strain at compressive failure
      double epsmax; // Strain at tensile failure

      /*** CONVERGED History Variables ***/
      double CminStrain;   // Smallest previous concrete strain (compression)
      double CunloadSlope; // Unloading (reloading) slope from CminStrain
      double CendStrain;   // Strain at the end of unloading from CminStrain
      int Cfailed;         // Flag indicating compressive failure

      /*** CONVERGED State Variables ***/
      double Cstrain;
      double Cstress;   
	  double Ctangent;	// Don't need Ctangent other than for revert and sendSelf/recvSelf
						// Storing it is better than recomputing it!!!

      /*** TRIAL History Variables ***/
      double TminStrain;
      double TunloadSlope;
      double TendStrain;
      int Tfailed;

      /*** TRIAL State Variables ***/
      double Tstrain;
      double Tstress;
      double Ttangent; // Not really a state variable, but declared here
                       // for convenience

      void determineTrialState (double dStrain);

      void reload();
      void unload();
      void envelope();
};


#endif


