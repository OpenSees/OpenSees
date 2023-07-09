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


#ifndef ConcretewBeta_h
#define ConcretewBeta_h

// Written: 2013, Yuan Lu & Marios Panagiotou (UC Berkeley) 

#include <UniaxialMaterial.h>

class ConcretewBeta : public UniaxialMaterial {
 public:
	ConcretewBeta(int tag, double fpc, double ec0, double fcint, double ecint, double fcres, double ecres, double fct, double ftint, double etint, double ftres, double etres, double lambda, double alpha, double bint, double etbint, double bres, double etbres, double M, double E0, double fcc, double ecc);
	ConcretewBeta();
	~ConcretewBeta();

	const char *getClassType(void) const {return "ConcretewBeta";};

	int setTrialStrain(double strain, double strainRate = 0.0); 

	// extra thing, only call if you know what you're doing.
	int setTrialStrainwBeta(double strain, double et = 0.0, double strainRate = 0.0);
	
	double getStrain(void);      
	double getStress(void);
	double getTangent(void);
	double getBeta(void);
	double getInitialTangent(void);

	int commitState(void);
	int revertToLastCommit(void);    
	int revertToStart(void);  

	UniaxialMaterial *getCopy(void);

	int sendSelf(int commitTag, Channel &theChannel);  
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
  
	void Print(OPS_Stream &s, int flag =0);

	Response *setResponse (const char **argv, int argc, 
				   OPS_Stream &theOutputStream);
	int getResponse (int responseID, Information &matInformation);    


private:
	// material properties
	double fpc;
	double ec0;
	double fcint;
	double ecint;
	double fcres;
	double ecres;
	double fct;
	double ftint;
	double etint;
	double ftres;
	double etres;
	double lambda;
	double alpha;
	double M;
	double fcc;
	double ecc;

	//For computation of beta
	double bint;
	double etbint;
	double bres;
	double etbres;

	// stored values, derived but not to be re-derived. :)
	void updateStoredValues();
	double et0; 
	double alphafct; // alpha*fct, place to unload to.
	double eaftc; // strain associated with alpha*fct
	double ElinearsoftcP1;
	double ElinearsoftcP2;
	double E0;
	double lambdaM;
	//double aLoad;
	//double bLoad;
	//double cLoad;
	//double dLoad;

	/*** CONVERGED History Variables ***/
	double CMaxStrainCompr;
	double CMaxStressCompr; // with beta factored in already.
	double CMaxStressComprPure;  // without beta factored in already.
	double CMaxStrainTens;
	double CMaxStressTens;

	/*** CONVERGED State Variables ***/
	double Cstrain;
	double Cstress;   // with beta factored in already.
	double Ctangent;   // derived value
	double Cbeta;
	  
	/*** TRIAL History Variables ***/
	double TMaxStrainCompr;
	double TMaxStressCompr; // with beta factored in already.
	double TMaxStressComprPure;  // without beta factored in already.
	double TMaxStrainTens;
	double TMaxStressTens;
	  
	/*** TRIAL State Variables ***/
	double Tstrain;
	double Tstress;  // with beta factored in already.
	double Ttangent;  // derived value
	double Tbeta;

	int setValues(double newStrain, double beta, double & newStress, double & newStressPure, double & newTangent);
	double computeBeta(double newStrain, double et);
};

#endif
