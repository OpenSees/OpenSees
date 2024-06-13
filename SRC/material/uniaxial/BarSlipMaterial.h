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
// $Date: 2004-10-06 19:21:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BarSlipMaterial.h,v $
                                                                        
                                                                        
// Written: NM (nmitra@u.washington.edu) 
// Created: January 2002
// Updated: September 2004
//
// Description: This file contains the class definition for 
// bar-slip material. The file generates the 4 point envelope for both positive 
// and negative loading and is basically a wrapper for the Pinching4 material at it's outset.
// Updates: new damage calculations and user interfaces


#ifndef BarSlipMaterial_h
#define BarSlipMaterial_h

#include <UniaxialMaterial.h>
//#include <Pinching4Material.h>
#include <math.h>
#include <Matrix.h>
#include <Vector.h>
#include <FileStream.h>

class BarSlipMaterial : public UniaxialMaterial
{
public :
	BarSlipMaterial(int tag,
		double fc, double fy, double Es, double fu,
		double Eh, double db, double ld, int nbars, double width, double depth,
		int bsflag, int type);

	BarSlipMaterial(int tag,
		double fc, double fy, double Es, double fu,
		double Eh, double db, double ld, int nbars, double width, double depth,
		int bsflag, int type, int damage, int unit);

	BarSlipMaterial();
	~BarSlipMaterial();

    const char *getClassType(void) const {return "BarSlipMaterial";}

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	UniaxialMaterial *getCopy(void);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);

	void Print(OPS_Stream &s, int flag = 0);

protected:

private:
	// bond strength flag
   int bsflag;   // 1 --- weak or 0 --- strong bond strength

	// unit used
	int unit;     // 1 --- mpa or 2 --- psi (the coefficients of shear stress will change accordingly)

	// type of j to be used for beam top barslip it is "1" or it is "0" for beam bottom bar slip or column
	int type;

	int damage; // type of damage in the material;

	// dimensions
	double width;  // width of the member
	double depth;  // depth of the member

	// material used at onset
//	Pinching4Material* material;
	Vector envlpPosStress; Vector envlpPosStrain; 
	Vector envlpNegStress; Vector envlpNegStrain;


	// concrete properties
	double fc;   // compressive strength of concrete
	
	// steel properties
	double fy;    // yield strength of steel
	double Es;    // modulus of elasticity of steel
	double fu;    // ultimate strength of steel
	double Eh;    // modulus of hardening of steel
	double db;    // reinforcing steel bar diameter
	int nbars;  // no of reinforcing bars

	// anchorage length
	double ld;  // anchorage length of the reinforcing steel

		// bond strengths
	double tauET;  // bond strength of steel that is elastic and carries tensile load
	double tauYT;  // bond strength of steel that has yielded in tension
	double tauEC;  // bond strength of steel that is elastic and carries compression load
	double tauYC;  // bond strength of steel that has yielded in compression
	double tauR;   // residual bond strength of steel

	// unloading-reloading parameters 
	double rDispP; double rForceP; double uForceP;
	double rDispN; double rForceN; double uForceN;

	// Damage parameters
	double gammaK1; double gammaK2; double gammaK3; double gammaK4; double gammaKLimit;
	double gammaD1; double gammaD2; double gammaD3; double gammaD4; double gammaDLimit;
	double gammaF1; double gammaF2; double gammaF3; double gammaF4; double gammaFLimit;
	double gammaE;


	// positive and negative envelopes
	Matrix eP;
	Matrix eN;

	//**************************************************************************
	// Trial Set variables
	double Tstrain; double Ttangent; double Tstress;

	// Converged Material History parameters
	int Cstate;
	double Cstrain;
	double Cstress;
	double CstrainRate;
	double lowCstateStrain;
	double lowCstateStress;
	double hghCstateStrain;
	double hghCstateStress;
	double CminStrainDmnd;
	double CmaxStrainDmnd;
	double Cenergy;
	double CgammaK;
	double CgammaD;
	double CgammaF;
	double gammaKUsed;
	double gammaFUsed;

	// Trial Material History Parameters
	int Tstate;
	double dstrain;
	double TstrainRate;
	double lowTstateStrain;
	double lowTstateStress;
	double hghTstateStrain;
	double hghTstateStress;
	double TminStrainDmnd;
	double TmaxStrainDmnd;
	double Tenergy;
	double TgammaK;
	double TgammaD;
	double TgammaF;

	// strength and stiffness parameters;
	double kElasticPos;
	double kElasticNeg;
	double kElasticPosDamgd;
	double kElasticNegDamgd;
	double uMaxDamgd;
	double uMinDamgd;

	Vector envlpPosDamgdStress; Vector envlpNegDamgdStress;

	double kunload;
	Vector state3Stress; Vector state3Strain; Vector state4Stress; Vector state4Strain;
	double elasticStrainEnergy;
	double energyCapacity;
//****************************************************************

	void getBondStrength(void);
	void getBarSlipEnvelope(void);
	void createMaterial(void);


		void SetEnvelope(void);
	void getstate(double, double);
	double posEnvlpStress(double);
	double posEnvlpTangent(double);
	double negEnvlpStress(double);
	double negEnvlpTangent(double);

	void getState3(Vector& , Vector& , double);
	void getState4(Vector& , Vector& , double);
	double Envlp3Tangent(Vector , Vector , double);
	double Envlp3Stress(Vector , Vector , double);
	double Envlp4Tangent(Vector , Vector , double);
	double Envlp4Stress(Vector , Vector , double);

	void updateDmg(double);

//FileStream* fn;
//FileStream* fg;

};
#endif
