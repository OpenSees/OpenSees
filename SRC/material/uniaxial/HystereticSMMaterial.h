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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-12-18 23:40:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Hysteretic7ptMaterial.h,v $

// Written: MHS
// Created: July 2000
// Modified to 7 points, etc: Silvia Mazzoni
//	April 2023
//
// Description: This file contains the class definition for 
// Hysteretic7ptMaterial.  Hysteretic7ptMaterial provides the implementation
// of a one-dimensional hysteretic model with pinching of both
// force and deformation, damage due to deformation and energy, and
// degraded unloading stiffness based on maximum ductility.  This
// is a modified implementation of Hyster2.f90 by Filippou.

#ifndef HystereticSMMaterial_h
#define HystereticSMMaterial_h

#include <UniaxialMaterial.h>

class HystereticSMMaterial : public UniaxialMaterial
{
 public:

	HystereticSMMaterial(int tag, const Vector& posEnv, const Vector& negEnv, const Vector& pinchArray, const Vector& damageArray, double beta ,
		const Vector& degEnvArray, const Vector& forceLimitStates, const Vector& defoLimitStates, Vector& internalValues, int YXorder = 1, int printInput = 0);
	//HystereticSMMaterial(int tag,
	//	double m1p, double r1p, double m2p, double r2p, double m3p, double r3p,
	//	double mom4p, double rot4p, double mom5p, double rot5p, double mom6p, double rot6p, double mom7p, double rot7p,
	//	double m1n, double r1n, double m2n, double r2n, double m3n, double r3n,
	//	double mom4n, double rot4n, double mom5n, double rot5n, double mom6n, double rot6n, double mom7n, double rot7n,
	//	double pinchX, double pinchY, const Vector& theLSforce, const Vector& theLSdefo, double degEnvFactor,
	//	double damfc1 = 0.0, double damfc2 = 0.0,
	//	double beta = 0.0);

  HystereticSMMaterial();
  ~HystereticSMMaterial();

  const char *getClassType(void) const {return "HystereticSMMaterial";};
  
  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void) { return E1p; };
  double getMUy(void);
  double getThetaP(void);



  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  
  UniaxialMaterial *getCopy(void);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  

  Response* setResponse(const char** argv, int argc, OPS_Stream& theOutput);
  int getResponse(int responseID, Information& matInfo);

  //by SAJalali
  double getEnergy() { return CenergyD; }

 protected:
  
 private:
  // Pinching parameters
  double pinchX;		// Deformation pinching
  double pinchY;		// Force pinching
  
  // Damage parameters
  double damfc1;		// Deformation
  double damfc2;		// Energy

  // Degrading envelope parameters
  double degEnvp;
  double degEnvn;
  
  
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
  
  double mom4p, rot4p;
  double mom5p, rot5p;
  double mom6p, rot6p;
  double mom7p, rot7p;
  double mom4n, rot4n;
  double mom5n, rot5n;
  double mom6n, rot6n;
  double mom7n, rot7n;


  double E1p, E1n;
  double E2p, E2n;
  double E3p, E3n;
  double E4p, E4n;
  double E5p, E5n;
  double E6p, E6n;
  double E7p, E7n;
  double Eup, Eun;

  double energyA;

  Vector posEnv;
  Vector negEnv;
  int nposEnv;
  int nnegEnv;

  Vector pinchArray;
  Vector damageArray;
  int npinchArray;
  int ndamageArray;

  // Unloading parameter
  double beta;

  int printInput;
  int YXorder;
  int xIndexIncr;
  int yIndexIncr;


  Vector defoLimitStates;
  Vector forceLimitStates;
  int nDefoLimitStates;
  int nForceLimitStates;
  //double degEnvFactor;
  Vector degEnvArray;
  int ndegEnvArray;


  void setEnvelope(void);
  
  double posEnvlpStress(double strain);
  double negEnvlpStress(double strain);
  
  double posEnvlpTangent(double strain);
  double negEnvlpTangent(double strain);
  
  double posEnvlpRotlim(double strain);
  double negEnvlpRotlim(double strain);
  
  void positiveIncrement(double dStrain);
  void negativeIncrement(double dStrain);
};

#endif
