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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/Bilinear.h,v $
//
//
// Bilinear.h: implementation of the Bilinear class from Fortran version.
// Originally from SNAP PROGRAM by Luis Ibarra and Prof H.K. Krawinkler
//
// Written: A. Altoontash & Prof. G. Deierlein 12/01
// Revised: 03/02
//
// Purpose: This file contains the implementation for the Bilinear class.
//
//////////////////////////////////////////////////////////////////////

// Bilinear.h: interface for the Bilinear class.
//
//////////////////////////////////////////////////////////////////////

#ifndef BILINEAR_H
#define BILINEAR_H

#include <UniaxialMaterial.h>
#include <Vector.h>
#include <stdio.h>
#include <DamageModel.h>
#include <MaterialResponse.h>

class Bilinear : public UniaxialMaterial  
{
	public:
	Bilinear();
	Bilinear(int tag, Vector inputParam ,DamageModel *strength,DamageModel *stiffness,DamageModel *capping);
	virtual ~Bilinear();
	
	const char *getClassType(void) const {return "Bilinear";};
  
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
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  
	Response* setResponse(const char **argv, int argc, Information &matInfo);
	int getResponse(int responseID, Information &matInfo);

	void Print(OPS_Stream &s, int flag =0);
    int    setParameter             (const char **argv, int argc, Information &info);
    int    updateParameter          (int parameterID, Information &info);
	int    activateParameter        (int parameterID);

/*
  	// Reliability and sensitivity stuff
    double getInitialTangent        (void);

	double getStressSensitivity     (int gradNumber, bool conditional);
	double getStrainSensitivity     (int gradNumber);
	double getTangentSensitivity    (int gradNumber);
	double getDampTangentSensitivity(int gradNumber);
	double getRhoSensitivity        (int gradNumber);
	int    commitSensitivity        (double strainGradient, int gradNumber, int numGrads);
  */

 protected:
	void envelPosCap( double ekelstk, double fy, double ekhard, double dcap,
						   double ekcap, double fRes, double *fuPos, double d, double *f, double *ek );
	void envelNegCap( double ekelstk, double fy, double ekhard, double dcap,
						   double ekcap, double fRes, double *fuNeg, double d, double *f, double *ek );
	void recordInfo(int cond =0);
 
 private:
	 // Input parameters
	double elstk, fyieldPos ,fyieldNeg, alfa;			// Main properties
	double alfaCap, capDispPos, capDispNeg, Resfac;		// Cap properties
	int flagCapenv;
	DamageModel *StrDamage;
	DamageModel *StfDamage;
	DamageModel *CapDamage;
	
	// Hystory data
	double hsTrial[17], hsCommit[17], hsLastCommit[17];
	
	FILE *OutputFile;		// For debugging

	// Sensitivity related variables
    int parameterID;
	Matrix *SHVs;
};

#endif
