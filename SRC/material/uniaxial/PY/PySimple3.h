/* *********************************************************************
**    Module:	PySimple3.h 
**
**    Purpose:	Provide a p-y spring for OpenSees to better capture
**              small strain nonlinear and viscoelastic behavior.            
**
**    Developed by Benjamin Turner and Scott J. Brandenberg
**    (C) Copyright 2013 and 2017, All Rights Reserved.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2012/12/31
// $Source: /OpenSees/SRC/material/uniaxial/PY/PySimple3.h

#ifndef PYSIMPLE3_H
#define PYSIMPLE3_H

// Written: SJB
// Created: Dec 2013
// tested and checked: BJT, 2015
//
// Description: This file contains the class definition for PySimple3.
// 

#include <UniaxialMaterial.h>

class PySimple3 : public UniaxialMaterial
{
  public:
    PySimple3(int tag, int classtag, double pult, double pyield, double kmax, double Hmod, double dashpot);
    PySimple3();
    ~PySimple3();


    const char *getClassType(void) const {return "PySimple3";};

    int setTrialStrain(double y, double yRate); 
    double getStrain(void);
	double getStrainRate(void);
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);
	double getDampTangent(void);
	double getResidual(double ke, double Cp, double Tp, double dy, double pu, double C, double Tpin, double dashpot, double tstepCurrent, double dyELast, double CyeTotal, double tstepLast, double Pveguess, double bump);
	int sign(double val);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

   
  protected:

	  // Material parameters
    double pult;		// Ultimate capacity
	double pyield;		// Yield force
    double ke;			// Elastic modulus
	double C;			// Plastc modulus parameter
	double dashpot;		// Dashpot coefficient for viscoelastic radiation damping

  private:

	
    // Committed history variables for entire p-y material
	
	double Cy;
	double Cp;
	double CpinF;		// p at start of current plastic loading cycle in forward dir.
	double CpinB;		// p at start of current plastic loading cycle in backward dir.
	double CpinLast;	// value of pin for last converged step.
	double Cyin;		// y at start of current plastic loading cycle
	double Cpalpha;		// p at center of elastic region
	double Ctangent;
	double dyLast;		//dy for converged load step; will be used to keep track of change in velocity during next load step
	double dyELast;		//dyE for converged load step; will be used to keep track of change in elastic velocity during next load step
	double CyeTotal;	//total elastic displacement (not incremental, from beginning of analysis) at start of load step
	double CtstepLast;	//last committed timestep, which will be the EQ record time increment except for the first yielding step when it will get sub-stepped
	int CLastYieldDir;	// Direction of loading increment during previous plastic loading cycle
	
	// Trial history variables for entire p-y material
    double Tp;			// Trial p
    double Ty;			// Trial y
	double TpinUse;		// Trial value of pin that will be used for current step
    double TpinF;		// Trial value of pin in forward dir.
	double TpinB;		// Trial value of pin in backward dir.
	double Tyin;		// Trial value of yin
	double Tpalpha;		// Trial value of p at center of elastic region
    double Ttangent;	// Trial tangent
	double TyRate;		// Trial disp rate
	double dy;
	double dyP;
	double dyE;
	double pn1_a, pn1_b, R1, R2, pn1_guess, R_guess;
	double yRate, tstep, ypRate, TdyP, yeRate, TdyE;		// Trial values of plastic and elastic displacement increments and displacement increment rates
	double TyeTotal;	//Trial value of total elastic displacement
	
	double lam;			//lambda, plastic multiplier
	double lamLB;
	double lamUB;		//lower- and upper-bound guesses for lambda
	double lamIn;		//trial value of lambda passed to getPlasticForce function
	double ypDot;		//plastic displacement rate
	double P1;			//Force (dP1/dt) in viscoelastic spring component
	double P2;			//Force (dP2/dt) in plastic spring component
	double yLast;		//y from the timestep before the currently committed timestep, i.e. y(i)=y, y(i-1)=Cy, y(i-2)=yLast
	double residualLam;	//residual when iterating to determine correct value of plastic multiplier
	double Rlam1;		//residual for lambda using lower- and upper-bound guesses
	double Rlam2;
	double sysTimeStep;	//system time
	double signdyLast;	//direction of last converged load step
	double bumped;		//flag keeps track of whether or not trial state was shifted to yield surface
	double tstepBump;		//timestep over shift
	double pP1;
	double dyELastUse;
	double signPalphaNew;
	double signBump;
	double P1veGuess;	//the purely visco-elastic force prediction
	double tstep1;		//sub-increments of timestep used when first yield occurs
	double tstep2;

	
	
	double f;			// Value of yield function
	int TLastYieldDir;	// Direction of loading increment during current plastic loading cycle
	int signdy;
	double Residual;	// Residual error in stress function
	double dResidual_dTp; // Gradient of residual expression with respect to trial stress

	double initialTangent;
};


#endif
