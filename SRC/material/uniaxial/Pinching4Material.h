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
                                                                        
// $Revision: 1.4 $
// $Date: 2010-04-06 20:18:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Pinching4Material.h,v $
                                                                        
// Written: NM (nmitra@u.washington.edu) 
// Created: December 2001
// Updated: September 2004
//
// Description: This file contains the class definition for 
// Pinching material which is defined by 4 points on the positive and 
// negative envelopes and a bunch of damage parameters. The material accounts for
// 3 types of damage rules : Strength degradation, Stiffness degradation, 
// unloading stiffness degradation. 
// Updates: damage calculations and several bug fixes


#ifndef Pinching4Material_h
#define Pinching4Material_h

#include <UniaxialMaterial.h>
#include <FileStream.h>
#include <OPS_Stream.h>
#include <Vector.h>

class Pinching4Material : public UniaxialMaterial
{
public :
	Pinching4Material(int tag,
		double stress1p, double strain1p, double stress2p, double strain2p,
		double stress3p, double strain3p, double stress4p, double strain4p,
		double stress1n, double strain1n, double stress2n, double strain2n,
		double stress3n, double strain3n, double stress4n, double strain4n,
		double rDispP, double rForceP, double uForceP,
		double rDispN, double rForceN, double uForceN,
		double gammaK1, double gammaK2, double gammaK3,
		double gammaK4, double gammaKLimit,
		double gammaD1, double gammaD2, double gammaD3,
		double gammaD4, double gammaDLimit,
		double gammaF1, double gammaF2, double gammaF3,
		double gammaF4, double gammaFLimit, double gammaE, int DmgCyc);

	Pinching4Material(int tag,
		double stress1p, double strain1p, double stress2p, double strain2p,
		double stress3p, double strain3p, double stress4p, double strain4p,
		double rDispP, double rForceP, double uForceP,
		double gammaK1, double gammaK2, double gammaK3,
		double gammaK4, double gammaKLimit,
		double gammaD1, double gammaD2, double gammaD3,
		double gammaD4, double gammaDLimit,
		double gammaF1, double gammaF2, double gammaF3,
		double gammaF4, double gammaFLimit, double gammaE, int DmgCyc);

	Pinching4Material();
	~Pinching4Material();

        const char *getClassType(void) const {return "Pinching4Material";}

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

	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);

	double getEnergy(void) { return Tenergy; }

protected:

private:
	// Backbone parameters
	    double stress1p; double strain1p; double stress2p; double strain2p;
		double stress3p; double strain3p; double stress4p; double strain4p;
		double stress1n; double strain1n; double stress2n; double strain2n;
		double stress3n; double strain3n; double stress4n; double strain4n;
		Vector envlpPosStress; Vector envlpPosStrain; 
		Vector envlpNegStress; Vector envlpNegStrain;

	// Damage parameters
	double gammaK1; double gammaK2; double gammaK3; double gammaK4; double gammaKLimit;
	double gammaD1; double gammaD2; double gammaD3; double gammaD4; double gammaDLimit;
	double gammaF1; double gammaF2; double gammaF3; double gammaF4; double gammaFLimit;
	double gammaE;
	double TnCycle, CnCycle; // number of cycles contributing to damage calculation
	int DmgCyc; // flag for indicating whether no. of cycles are to be used for damage calculation

	// unloading-reloading parameters
	double rDispP; double rForceP; double uForceP;
	double rDispN; double rForceN; double uForceN;

	Vector state3Stress; Vector state3Strain; Vector state4Stress; Vector state4Strain;

	Vector envlpPosDamgdStress; Vector envlpNegDamgdStress;

	// Trial State Variables
	double Tstress;
	double Tstrain;
	double Ttangent;

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

	// strength and stiffness parameters
	double kElasticPos;
	double kElasticNeg;
	double kElasticPosDamgd;
	double kElasticNegDamgd;
	double uMaxDamgd;
	double uMinDamgd;


	// energy parameters
	double energyCapacity;
	double kunload;
	double elasticStrainEnergy;

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
	void updateDmg(double, double);


};
#endif
