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
                                                                        
// $Revision: 1.2 $
// $Date: 2007-06-28 21:46:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete07.h,v $
                                                                        
// Written: Jon Waugh, Iowa State University
// Created: 11/2005
//
// Description: This file contains the class definition for Uniaxial material Concrete07 
//				A simplified form of Chang & Mander Concrrete model from 1994.

#ifndef Concrete07_h
#define Concrete07_h

#include <UniaxialMaterial.h>


class Concrete07 : public UniaxialMaterial
{
	public:
		Concrete07 (int tag, double FPC, double EPSC0, double EC, double FPT, double ESPT0, double XCRP, double XCRN, double R);
		Concrete07 ();
		~Concrete07();

		int setTrialStrain(double strain, double strainRate = 0.0); 
		int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
		double getStrain(void);      
		double getStress(void);
		double getTangent(void);
		double getInitialTangent(void) {return Ec;}

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
      double fpc;				// Compressive strength
      double epsc0;				// Strain at compressive strength
	  double Ec;				// Initial Young's modulus;
	  double fpt;				// Tensile strength
	  double epst0;				// Strain at tensile strength

	  /*** Model Variables ***/
      double xcrn;				// nondimensional critical strain used to determine spalling
	  double xsp;				// nondimensional spalling strain
	  double xcrp;				// nondimensional critical strain used to determine cracking
	  double xcrk;				// nondimensional cracking strain
	  double nn;				// parameter to control stress-strain curve in compression
	  double np;				// parameter to control stress-strain curve in tension
	  double r;					// factor from Tsai's equation
	  double e0;				// origin of the shifted tension side.

      /*** CONVERGED History Variables ***/
      double CminStrain;		// Smallest previous concrete strain (compression)
	  double CminStress;		// Stress that goes with CminStrain
	  double CUnloadStrain;		// Previous partial unloading strain
	  double CUnloadStress;		// Previous partial unloading stress
	  double CUnloadStiffness;	// Previous partial unloading stiffness
	  double CmaxStrain;		// Largest previous concrete strain (tension)
	  double CmaxStress;		// Stress that goes with CmaxStrain (tension)
	  double CReloadStrain;		// Previous partial reloading strain
	  double CReloadStress;		// Previous partial reloading stress
	  double CReloadStiffness;	// Previous partial reloading stiffness
	  double C13Zero;			// Strain where we reversed onto rule 13.
	  double C13Strain;			// Strain where we left reversed on rule 13
	  double C13Stress;			// Stress where we left reversed on rule 13

	  int Cloading;				// Flag for stressing/unstressing
								// 1 = compressing the concrete
								// -1 = uncompressing the concrete
								// 0 = initial state
	  bool Ccracked;			// Flag for if the concrete has cracked
								// true the concrete has cracked and assumed to have no tension capacity
								// false the concrete has not cracked yet
	  int Crule;				// The rule we are following for this step

      /*** CONVERGED State Variables ***/
      double Cstrain;
      double Cstress;   
      double Ctangent;			// Don't need Ctangent other than for revert and sendSelf/recvSelf
								// Storing it is better than recomputing it!!!

      /*** TRIAL History Variables ***/
      double TminStrain;
	  double TminStress;
	  double TUnloadStrain;
	  double TUnloadStress;
	  double TUnloadStiffness;
	  double TmaxStrain;
	  double TmaxStress;
	  double TReloadStrain;
	  double TReloadStress;
	  double TReloadStiffness;
	  double T13Zero;
	  double T13Strain;
	  double T13Stress;
	  int Tloading;
	  bool Tcracked;
	  int Trule;

      /*** TRIAL State Variables ***/
      double Tstrain;
      double Tstress;
      double Ttangent; // Not really a state variable, but declared here
                       // for convenience

	  /*** Private functions to make life easier ***/
	  // Calculates trial state variables based on the trail strain
      void determineTrialState (double dStrain);

	  // Determines is a strain reversal has occurred based on the trial strain
	  void detectStrainReversal (double dStrain);

	  // Calculates the value of y(x) and z(x)
	  void calculateYandZ(double x, double& y, double& z, double n);

	  // Calculate the stress of the transition curve;
	  void calculateStressTransition(double &fc, double &Et, double ec, double eI, double fI, double EI, double eF, double fF, double EF, int rule);
      
	  // Calculate the envelope stress for a given nondimensional strain.
	  void envelope(double x, double& fc, double& Et, int flag);

	  // Calculate the stress-strain curve for Rule 13
	  void calculate13Stress(double &fc, double &Et, double ec, double eI, double eF, double fF, double EF);


};

#endif

