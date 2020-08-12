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
                                                                        
// Written: KJE
// Created: July 2001
// Modified July 2002
// Based on HystereticMaterial by MHS
//
// Description: This file contains the class definition for 
// LimitStateMaterial.  LimitStateMaterial is
// a one-dimensional hysteretic model with pinching of both
// force and deformation, damage due to deformation and energy, and
// degraded unloading stiffness based on maximum ductility.  Based on 
// HystereticMaterial which is a modified implementation of Hyster2.f90 
// by Filippou.  
//
// Modified July 2002 by KJE to include the use of limitCurve to 
// detect failure in hysteretic material (see commitState).
// Option only available for 3 point specification. Behaves 
// same as HystereticMaterial if no limit curve specified.
//
// All code specific to LimitStateMaterial separated by ////////////////

#ifndef LimitStateMaterial_h
#define LimitStateMaterial_h

#include <UniaxialMaterial.h>
class LimitCurve;


class LimitStateMaterial : public UniaxialMaterial
{
	public:
		LimitStateMaterial(int tag,
			double mom1p, double rot1p, double mom2p, double rot2p,
			double mom3p, double rot3p,
			double mom1n, double rot1n, double mom2n, double rot2n,
			double mom3n, double rot3n,
			double pinchX, double pinchY,
			double damfc1 = 0.0, double damfc2 = 0.0,
			double beta = 0.0);
		LimitStateMaterial(int tag,
			double mom1p, double rot1p, double mom2p, double rot2p,
			double mom1n, double rot1n, double mom2n, double rot2n,
			double pinchX, double pinchY,
			double damfc1 = 0.0, double damfc2 = 0.0,
			double beta = 0.0);
		LimitStateMaterial(int tag,
			double mom1p, double rot1p, double mom2p, double rot2p,
			double mom3p, double rot3p,
			double mom1n, double rot1n, double mom2n, double rot2n,
			double mom3n, double rot3n,
			double pinchX, double pinchY,
			double damfc1, double damfc2,
			double beta, LimitCurve &theCurve, 
			int curveType, int degrade);
		LimitStateMaterial();
		~LimitStateMaterial();

		const char *getClassType(void) const {return "LimitStateMaterial";};

	    int setTrialStrain(double strain, double strainRate = 0.0);
	    double getStrain(void);
	    double getStress(void);
	    double getTangent(void);
	    double getInitialTangent(void) {return E1p;};

	    int commitState(void);
	    int revertToLastCommit(void);
	    int revertToStart(void);

	    UniaxialMaterial *getCopy(void);
		
	    int sendSelf(int commitTag, Channel &theChannel);  
	    int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);    
    
	    void Print(OPS_Stream &s, int flag =0);

	    int setParameter (const char **argv, int argc, Parameter &param);



	    Response *setResponse (const char **argv, int argc, 
					   OPS_Stream &theOutputStream);
	    int getResponse (int responseID, Information &matInformation);    

	protected:

	private:
		// Pinching parameters
		double pinchX;		// Deformation pinching
		double pinchY;		// Force pinching

		// Damage parameters
		double damfc1;		// Deformation
		double damfc2;		// Energy

		// Unloading parameter
		double beta;

		// Trial history variables
		double TrotMax;
		double TrotMin;
		double TrotPu;
		double TrotNu;
		double TenergyD;
		int TloadIndicator;

		// Trial state variables
		double Ttangent;
		double Tstress;
		double Tstrain;

		// Converged history variables
		double CrotMax;
		double CrotMin;
		double CrotPu;
		double CrotNu;
		double CenergyD;
		int CloadIndicator;

		// Converged state variables
		double Cstress;
		double Cstrain;

		// Backbone parameters
		double mom1p, rot1p;
		double mom2p, rot2p;
		double mom3p, rot3p;
		double mom1n, rot1n;
		double mom2n, rot2n;
		double mom3n, rot3n;

		double E1p, E1n;
		double E2p, E2n;
		double E3p, E3n;

	double pinchX_orig;
	double pinchY_orig;
	double damfc1_orig;
	double damfc2_orig;
	double beta_orig;
	double mom1p_orig;
	double rot1p_orig;
	double mom2p_orig;
	double rot2p_orig;
	double mom3p_orig;
	double rot3p_orig;
	double mom1n_orig;
	double rot1n_orig;
	double mom2n_orig;
	double rot2n_orig;
	double mom3n_orig;
	double rot3n_orig;

	int constructorType;


		double energyA;

		//////////////////
		// variables for LimitCurve option
		int degrade;			// option to force degradation in opposite direction due to degradation
								// in current direction (requires testing, not documented)
								// (0 = off, 1 = on)
		int curveType;			// type of limit curve 
								// (0 = no curve, 1 = axial curve, 
								// all other curves can be any other integer)
		LimitCurve *theCurve;	// curve defining limit state
		int CstateFlag;			// Flag indicating if the limit state surface has been reached
								// (not reached = 0, reached for first time = 1, 
								// on limit curve = 2, off limit curve = 3, at residual capacity = 4)
		double Ploss;			// loss of axial load capacity.
		double Kdeg;			// degrading slope after failure
		double Eelasticp;		// elastic stiffness used if off failure surface (for axial curve only)
		double Eelasticn;
		//////////////////

		void setEnvelope(void);

		double posEnvlpStress(double strain);
		double negEnvlpStress(double strain);

		double posEnvlpTangent(double strain);
		double negEnvlpTangent(double strain);

		double posEnvlpRotlim(double strain);
		double negEnvlpRotlim(double strain);

		void positiveIncrement(double dStrain);
		void negativeIncrement(double dStrain);

		///////////////
		int getNewBackbone(int flag);
		int mirrorBackbone(void);
		///////////////
};

#endif
