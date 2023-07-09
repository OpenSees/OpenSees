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
// $Date: 2006-08-03 23:44:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/CloughDamage.h,v $
//
//
// CloughDamage.h: implementation of the CloughDamage class from Fortran version.
// Originally from SNAP PROGRAM by Prof H.K. Krawinkler
//
// Written: A. Altoontash & Prof. G. Deierlein 12/01
// Revised: 03/02
//
// Purpose: This file contains the implementation for the CloughDamage class.
//
//////////////////////////////////////////////////////////////////////

// CloughDamage.h: interface for the CloughDamage class.
//
//////////////////////////////////////////////////////////////////////

#ifndef CloughDamage_H
#define CloughDamage_H

#include <UniaxialMaterial.h>
#include <Vector.h>
#include <stdio.h>

class DamageModel;

class CloughDamage : public UniaxialMaterial
{
public:
	CloughDamage();
	CloughDamage(int tag, Vector inputParam, DamageModel *strength, DamageModel *stiffness,DamageModel *accelerated, DamageModel *capping);
	virtual ~CloughDamage();

	const char *getClassType(void) const {return "CloughDamage";};
	
	int setTrialStrain(double d, double strainRate = 0.0);
	double getStrain(void);
	
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);
	int commitState(void);
	int revertToLastCommit(void);    
	int revertToStart(void);  
	
	//virtual
	UniaxialMaterial *getCopy(void);
	
	int sendSelf(int commitTag, Channel &theChannel);  
	int recvSelf(int commitTag, Channel &theChannel,FEM_ObjectBroker &theBroker);
	
	void Print(OPS_Stream &s, int flag =0);
	
protected:
	void envelPosCap(double fy, double alphaPos, double alphaCap,
		double cpDsp, double d, double *f, double *ek );
	
	void envelNegCap(double fy, double alphaNeg, double alphaCap,
		double cpDsp, double d, double *f, double *ek);
	
	void recordInfo(int cond =0);
	
	
private:
	// Input parameters
	double elstk,fyieldPos,fyieldNeg,alpha,Resfac;		//	Properties
	double capSlope,capDispPos,capDispNeg;	 // Cap

	// Parameters calculated from input data
	double dyieldPos,dyieldNeg;
	DamageModel *StrDamage;
	DamageModel *StfDamage;
	DamageModel *AccDamage;
	DamageModel *CapDamage;

	double hsTrial[24], hsCommit[24], hsLastCommit[24];
	
};

#endif
