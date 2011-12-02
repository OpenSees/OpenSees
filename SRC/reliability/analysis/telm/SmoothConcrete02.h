/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/SmoothConcrete02.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef SmoothConcrete02_h
#define SmoothConcrete02_h

#include <UniaxialMaterial.h>
#include <Parameter.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::scientific;
using std::showpoint;

class SmoothConcrete02 : public UniaxialMaterial
{
	public:
		SmoothConcrete02 
		(int tag, double fpc, double eco, double fpcu, double ecu, 
		 double etension, double gamma, double eta);
		SmoothConcrete02 ();
		~SmoothConcrete02();

		int setTrialStrain(double strain, double strainRate = 0.0); 

		double getStress () { return Tstress; }
		double getStrain () { return Tstrain; }
		double getTangent () { return Ttangent; }
		double getInitialTangent(){ return 2.0*fpc/epsc0; }
		int revertToLastCommit(void);    
		int revertToStart(void);        
		int commitState(void);

		UniaxialMaterial *getCopy(void);
    
		int sendSelf(int commitTag, Channel &theChannel);  
		int recvSelf(int commitTag, Channel &theChannel, 
	  				 FEM_ObjectBroker &theBroker);    
    
		void Print(OPS_Stream &s, int flag =0);

		// Reliability and sensitivity stuff
		int setParameter (const char **argv, int argc, Parameter &param);
		int updateParameter(int parameterID, Information &info);
		int activateParameter(int parameterID);
		double getStressSensitivity (int gradNumber, bool conditional);
		double getStrainSensitivity (int gradNumber);
		double getTangentSensitivity (int gradNumber);
		double getDampTangentSensitivity(int gradNumber);
		double getRhoSensitivity (int gradNumber);
		int commitSensitivity (double strainGradient, int gradNumber, int numGrads);

   protected:

   private:
		/*** Material Properties ***/
		double fpc;    // Compressive strength
		double epsc0;  // Strain at compressive strength
		double fpcu;   // Crushing strength
		double epscu;  // Strain at crushing strength
		double etension; // Strain limit with stress in tension region
		// Additional material parameters
		double gamma;
		double beta;
		double epsilon;

		/*** CONVERGED State ***/
		double Cstrain;	// Strain at previous convergence
		double Cstress;   // Stress at previous convergence
		double CminStrain;   // E_min at previsou cSmallest previous concrete strain (compression)
		double CendStrain;   // Strain at the end of unloading from CminStrain
		double Celastic;   // Elastic stiffness
		double Ctangent;   // Elastic stiffness
		bool   Cloading	;	 // loading indicator;
		int    Cstate;		// state indicator;

		/*** TRIAL History Variables ***/
		double Tstrain;	//  Trial Strain
		double Tstress;	//  Trial Stress
		double TminStrain;   // E_min at previsou cSmallest previous concrete strain (compression)
		double TendStrain;   // Strain at the end of unloading from CminStrain
		double Telastic;   // Strain at the end of unloading from CminStrain
		double Ttangent;   // Elastic stiffness
		bool   Tloading	;	 // loading indicator;
		int    Tstate;		// state indicator;

		// Sensitivity stuff
		int parameterID;
		double DTstressDhTstrainFix;
		double DeminDhUpdated;
		int UnconditionalDerivative;
		Matrix *SHVs;
//		ofstream output;
//		ofstream conc02_output2;
//		int iopen;


		/*** Functions ***/
		void envelope 
		(const double& eps, double& sig, double& tan);
		void envelopewithSens1 
		(const double& eps, double& sig, double& tan, double& dtandeps);
		void envelopeSens 
		(const double& eps, 
		 const double& dfpcdh, const double& dfpcudh, 
		 const double& depsc0dh, const double& depscudh,
		 double& dsigdh, double& dtandh);
//		void envlopewithSens2 
//		(const double& eps, const double& sig, const double& tan, 
//		 const double& dfpcdh, const double& dfpcudh, 
//		 const double& depsc0dh, const double& depscudh,
//		 double& dsigdh double dtandh);

		void unload 
		(const double& emin, const double& sigma,
		 double& eend, double& elastic);
		void unloadwithSens1 
		(const double& emin, const double& sigmin, const double& dsigmindemin,
		 double& eend, double& elastic, double& deenddemin, double& delasticdemin);
		void unloadwithSens2 
		(const double& emin, const double& sigmin, 
		 const double& dsigmindemin, const double& dsigmindh,
		 const double& dfpcdh, const double& dfpcudh, 
		 const double& depsc0dh, const double& depscudh,
		 double& eend, double& elastic, 
		 double& deenddemin, double& delasticdemin,double& deenddh, double& delasticdh);

		void polinomial 
		(const double& eps1, const double& sig1, const double& tan1, 
		 const double& eps2, const double& sig2, const double& tan2,
		 double& a, double& b, double& c, double& d); 
		void polinomialwithSens
		(const double& eps1, const double& sig1, const double& tan1, 
		 const double& eps2, const double& sig2, const double& tan2,
		 const double& a, const double& b, const double& c, const double& d, 
		 const double& deps1dh, const double& dsig1dh, const double& dtan1dh,
		 const double& deps2dh, const double& dsig2dh, const double& dtan2dh,
		 double& dadh, double& dbdh, double& dcdh, double& dddh);

		void solemin1
		(const double& sigma, const double& eps, 
		 double& emin, double& eend, double& elastic, 
		 double& amin, double& bmin, double& cmin, double& dmin);
		void solemin1Sens
		(const double& emin, const double eps, const double& dsigmadh, const double& depsdh,
		 const double& dfpcdh, const double& dfpcudh, 
		 const double& depsc0dh, const double& depscudh, 
		 double& demindh, double& dadh, double& dbdh, double& dcdh, double& dddh); 

		void solemin2
		(const double& sigma, const double& eps, 
		double& emin, double& eend, double& elastic); 
		void solemin2Sens
		(const double& emin, const double& eps, const double& dsigmadh, const double& depsdh,
		 const double& dfpcdh, const double& dfpcudh, 
		 const double& depsc0dh, const double& depscudh, 
		 double& demindh, double& deenddh, double& delasticdh);


		double envelopeSelf(double);
		double envelopeDeriv(double);

		bool path1;
		bool path2;
		bool path3;
		bool path4;
		bool path5;
		bool path6;
		bool path7;
		bool path8;
		bool path9;
		bool path10;
		bool path11;
		bool path12;
		bool path13;
		bool path14;
		bool path15;
		bool path16;
		bool path17;
		bool path18;
		bool path19;
		bool path20;
		bool path21;
		bool path22;
		bool path23;
		bool path24;
		bool path25;
		bool path26;
		bool path27;
		bool path28;
		bool path29;
		bool path30;


};


#endif


