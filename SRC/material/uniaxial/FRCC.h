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
// This code is developed and written by FERAS L. KHLEF (feraskhlef@gmail.com)
// $Revision: 1.0 $
// $Date: 2019-03-24$
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FRCC.h,v1 $                                                                                                                                              

#ifndef FRCC_h
#define FRCC_h
#include <UniaxialMaterial.h>
// File: FRCC.h
// Description: This file contains the class definition for FRCC.h.

class FRCC : public UniaxialMaterial
{
public:
	FRCC(int tag, double EC, double ET1, double FT1, double ET2, double FT2, double ETU, double FTU, double STU,
		double PARAT1, double PARAT2, double PARAT3, double PARAT4, double PARAT5, double PARAT6,
		double ECP, double FCP, double ECU, double PARAC1, double PARAC2, double PARAC3, double PARAC4, double PARAC5, double PARAC6, double RC);
	FRCC();
	~FRCC();

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void) { return Ec; }

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	UniaxialMaterial* getCopy(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);

protected:

private:
	/*** Input Parameters ***/
		//Tension inputs--------------------------------------------------------------------------------------------------
	double Ec;      // Initial modulus (tension and compression)
	double et1;	    // Strain at beginning of the hardening stage
	double ft1;		// stress at beginning of the hardening stage
	double et2;  	// Strain at the ending point of the hatdening/softening stage
	double ft2;  	// Stress at the ending point of the hatdening/softening stage
	double etu;	    // Ultimate tensile strain
	double ftu;	    // Stress at ultimate tensile strain
	double stu;	    // Tangent slope at ultimate tensile strain
	double paraT1;  // A parameter to define strain at the returning point on the envelope after reloading (tension)
	double paraT2;  // A parameter to define the amount of degridation in stress due to unloading/reloading (tension)
	double paraT3;  // A parameter to define the regression model for stiffness degridation (tension)
	double paraT4;  // A parameter to define the tangent slope of unloading curve at the intersection point with strain axis(tension)
	double paraT5;  // A parameter to define the shifting of the compression zone into the tension zone
	double paraT6;  // A parameter to define the tangent slope of the reloading curve after the transition zone 4(tension)
	// Compression inputs----------------------------------------------------------------------------------------------
	double ecp;	    // strain at Peak compression stress
	double fcp;  	// Peak compression stress
	double ecu;     // Ultimate compression strain
	double paraC1;  // A parameter to define strain at the returning point on the envelope after reloading (compression)
	double paraC2;  // A parameter to define the amount of degridation in stress due to unloading/reloading (compression)
	double paraC3;  // A parameter to define the regression model for stiffness degridation (compression)
	double paraC4;  // A parameter to define the tangent slope of unloading curve at the intersection point with strain axis(compression)
	double paraC5;  // A parameter to define the shifting of the tension zone into the compression zone
	double paraC6;  // A parameter to define the tangent slope of the reloading curve after the transition zone 2(compression)
	double rc;	    // A parameter to control the descending branch position and shape

/*** Model Variables ***/
	double Ycr;
	double Zcr;
	double Xsp;
	double Eh;
	double nc;
	double D;
	double Esec;
	double fnew;
	double Enew;
	double Eul;
	double Ea1;
	double Ea2;
	double epsca;
	double epsta;
	double epscMax;
	double sigcMax;
	double epstMax;
	double sigtMax;
	double epscRe;
	double epstRe;
	int FlagT0;
	int FlagC0;
	double RRR, AAA, EsecDesce, Xcrk;
	/*** CONVERGED History Variables ***/
	double Cepstul;  
	double CsigtMax; 
	double Cepscul;  
	double CsigcMax; 
	double CepstMax;     
	double CepscMax;     
	double CepstMaxStar; 
	double CepscMaxStar; 
	double CsigtMaxStar; 
	double CsigcMaxStar; 
	double CepstulStar;  
	double CepsculStar;  
	double CsigtulStar;  
	double CsigculStar;  
	double CepsculStar1;
	double CsigculStar1;
	double CepsculStar2;
	double CsigculStar2;
	double CepsculStar3;
	double CsigculStar3;
	double CepstulStar1;
	double CsigtulStar1;
	double CepstulStar2;
	double CsigtulStar2;
	double CepstulStar3;
	double CsigtulStar3;
	double CepstRe;
	double CepscRe;
	double Cepst0;
	double Cepsc0;
	double CTI;
	int CFlagT;
	int CFlagC;
	int Cid1;
	int Cid2;
	int Cid3;
	int Cid4;
	int Cid4_1;
	int Cid2_1;
	/*** CONVERGED State Variables ***/
	double Cstrain;
	double Cstress;
	double Ctangent;


	/*** TRIAL History Variables ***/
	double Tepstul;		 // Trial strain at reloading point (tension) 
	double TsigtMax;	 // Trial Unloading stress (tension)
	double Tepscul;		 // Trial strain at reloading point (compression)
	double TsigcMax;	 // Trial Unloading stress (compression)
	double TepstMax;     // Trial unloading strain (tension)
	double TepscMax;     // Trial unloading strain (compression)
	double TepstMaxStar; // Trial unloading strain for partitial unloading (tension)
	double TepscMaxStar; // Trial unloading strain for partitial unloading (compression)
	double TsigtMaxStar; // Trial unloading stress for partitial unloading (tension)
	double TsigcMaxStar; // Trial unloading stress for partitial unloading (tension) (compression)
	double TepstulStar;  // Trial strain at reloading point for partitial reloading (tension)
	double TepsculStar;  // Trial strain at reloading point for partitial reloading (compression)
	double TsigtulStar;  // Trial stress at reloading point for partitial reloading (tension)
	double TsigculStar;  // Trial stress at reloading point for partitial reloading (compression)
	double TepsculStar1;
	double TsigculStar1;
	double TepsculStar2;
	double TsigculStar2;
	double TepsculStar3;
	double TsigculStar3;
	double TepstulStar1;
	double TsigtulStar1;
	double TepstulStar2;
	double TsigtulStar2;
	double TepstulStar3;
	double TsigtulStar3;
	double TepstRe;
	double TepscRe;
	double Tepst0;
	double Tepsc0;
	double TTI;
	double TTTI;
	int TFlagT;
	int TFlagC;
	int Tid1;
	int Tid2;
	int Tid3;
	int Tid4;
	int Tid4_1;
	int Tid2_1;
	/*** TRIAL State Variables ***/
	double Tstrain;
	double Tstress;
	double Ttangent;

	/*** Defined Private functions ***/
	void RegExpr3(double epsun, double epsu, double EE, double PARA3);
	void RegExpr4(double epsun, double epsu, double EE, double PARA4);
	void RegExpr6(double epsun, double epsu, double EE, double PARA6);
	void TstressAndTtangent(double epsI, double fI, double EI, double epsF, double fF, double EF, double eps);
	void Dfunc(double r, double x, double n);
	void fnewEnew(double epsMax, double epsul, double sigMax, double mlsig);
	void EnvTen(double eps, double epst0);
	void EnvComp(double eps, double epsc0);
};
#endif