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
                                                                        
// Function contributed by  Quan Gu & Michele Barbato

// $Revision: 1.2 $
// $Date: 2009-03-27 19:19:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SmoothPSConcrete.h,v $
                                                                        
                                                                        
#ifndef SmoothPSConcrete_h
#define SmoothPSConcrete_h

// File: SmoothPSConcrete.h
//
//
// Description: This file contains the class definition for 
// SmoothPSConcrete.h 
//   - Popovics-Saenz envelope
//   - No tension
//   - Smoooth transition to  a linear unloading/reloading
//
// What: "@(#) SmoothPSConcrete.h, "


#include <UniaxialMaterial.h>

class SmoothPSConcrete : public UniaxialMaterial
{
 public:
  SmoothPSConcrete (int tag, double fpc, double fpcu,double Ec,double eco,  double ecu,   double  eta );
  SmoothPSConcrete ();
  ~SmoothPSConcrete();

      int setTrialStrain(double strain, double strainRate = 0.0); 
      int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
      double getStrain(void);      
      double getStress(void);
      double getTangent(void);
      double getInitialTangent(void) {
		//  return TEt;  Quan 2008 Jan 30.  
			return Ec;
	  }
 	double getInitialTangentSensitivity(int gradNumber);
 

      int commitState(void);
      int revertToLastCommit(void);    
      int revertToStart(void);        

      UniaxialMaterial *getCopy(void);
    
      int sendSelf(int commitTag, Channel &theChannel);  
      int recvSelf(int commitTag, Channel &theChannel, 
	  	 FEM_ObjectBroker &theBroker);    
    
      void Print(OPS_Stream &s, int flag =0);
		
//	  Response* setResponse(const char **argv, int argc, OPS_Stream &output);
//	  int getResponse(int responseID, Information &matInfo);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int    setParameter             (const char **argv, int argc, Parameter &param);
    int    updateParameter          (int parameterID, Information &info);
	int    activateParameter        (int parameterID);
	double getStressSensitivity     (int gradNumber, bool conditional);
	double getStrainSensitivity     (int gradNumber);

	int    commitSensitivity        (double strainGradient, int gradNumber, int numGrads);
// AddingSensitivity:END ///////////////////////////////////////////

   protected:

   private:
      /*** Material Properties ***/
      double fc;    // Compressive strength
      double eps0;  // Strain at compressive strength
      double fcu;   // curvature invertion concrete stress
      double epsu;  // curvature invertion concrete strain
      double  Ec ;       // initial stiffness modulus
      double  eta ;      // smoothing coefficient
        

	  int Monotonic_Envelope (double epsc, double * sig, double * Et);
      int Compute_epsp (void);
      int Transition_r(double epsc, double e1,double e2,double s1,double s2,double e1th,double e2th,double Et1,double Et2, double * sig, double * Et);
	  int Transition_p (double Delta);
      int compute_epsmax(double * epsmax,double * sigmax);

// AddingSensitivity:BEGIN //////////////////////////////////////////
	//  int Monotonic_Envelope_sens (double dfcdh, double deps0dh,double depsudh,double dfudh,double dEcdh,double depsdh){return 0;};
    //  int Compute_depspdh (double depsrdh,double dsigrdh,double deps0dh,double dfcdh,double dEcdh){return 0;};
    //  int Transition_r_sens(double e1,double e2,double s1,double s2,double e1th,double e2th,double Et1,double Et2,double depsdh,double de1dh,double de2dh,double ds1dh,double ds2dh,double dEt1dh,double dEt2dh,double dfcdh,double deps0dh,double depsudh,double dfudh,double dEcdh){return 0;};
	 // int Transition_p_sens (double Delta, double depsdh,double depspdh,double dDeltadh,double dEurdh){return 0;};
//	  int Monotonic_Envelope_Et_sens (double depsdh,double dfcdh,double deps0dh,double depsudh,double dfudh,double dEcdh){return 0;};
	  double Monotonic_Envelope_sens(double epsilon,double dfcdh,double deps0dh,double depsudh,double dfudh,double dEcdh,double depsdh);
	  double Compute_depspdh (double epsr,double sigr,double depsrdh,double dsigrdh,double deps0dh,double dfcdh,double dEcdh);
	  double Transition_r_sens(double epsc,double e1,double e2,double s1,double s2,double e1th,double e2th,double Et1,double Et2,
										   double depsdh,double de1dh,double de2dh,double ds1dh,double ds2dh,double dEt1dh,double dEt2dh,
										   double dfcdh,double deps0dh,double depsudh,double dfudh,double dEcdh);
	  double Transition_p_sens(double epsc,double epsp,double Delta,double Eur,double depsdh,double depspdh,double dDeltadh,double dEurdh);
	  double Monotonic_Envelope_Et_sens (double epsilon,double depsdh,double dfcdh,double deps0dh,double depsudh,double dfudh,double dEcdh);
// AddingSensitivity:END ///////////////////////////////////////////	  

	  // ===========================  new members =================   
		/*** CONVERGED History Variables ***/

		double Cepsr;
		double Csigr;
		double Cepsp;
		double CEur;
		int    Cflag;
		double Cepsr1;
		double Cepsr2;
		double Csigr1;
		double Csigr2;
		double CEt2;

       /*** CONVERGED State Variables ***/
		double Csig;
		double CEt;
	    double Cepsc; // total strain
 

		/*** TRIAL History Variables ***/
 
		double Tepsr;
		double Tsigr;
		double Tepsp;
		double TEur;
		int    Tflag;
		double Tepsr1;
		double Tepsr2;
		double Tsigr1;
		double Tsigr2;
		double TEt2;

		/*** TRIAL State Variables ***/
		double Tsig;
		double TEt;
	    double Tepsc;   // total strain
	    double Tdepsc;// total strain increment


		double epsmax;
		double sigmax;

 
       







// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////
};


#endif


