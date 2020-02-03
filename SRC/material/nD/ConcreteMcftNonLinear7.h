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
                                                                        
// $Revision: 1.19 $
// $Date: 2007/11/30 23:34:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete01.cpp,v $
                                                                        
// Written: Murat  
// Created: 03/08
//
// 
// Description: This file contains the class definition for 
// ConcreteMcft.h adapted from Modified Compression Field Theory (Vecchio and Collins) and
// MultiStrain section model of Petrangeli et al. 1991
//

#include <NDMaterial.h>
#include <Information.h>
#include <OPS_Globals.h>
#include <Matrix.h>
#include <Vector.h>
#include <MaterialResponse.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <TaggedObject.h>
#include <OPS_Globals.h>


#ifndef ConcreteMcftNonLinear7_h
#define ConcreteMcftNonLinear7_h


class ConcreteMcftNonLinear7  : public NDMaterial
{


public:
  
 //constructor
  ConcreteMcftNonLinear7 (int tag, double fcu, double ecu, double Ec, double fcr, double Esv, double fyv, double alphaV, double RoV);
  ConcreteMcftNonLinear7 (); 
  ~ConcreteMcftNonLinear7 ();

    NDMaterial* getCopy (void);
    NDMaterial* getCopy (const char * type);

	int setTrialStrain (const Vector &v);
    int setTrialStrain (const Vector &v, const Vector &r);
    int setTrialStrainIncr (const Vector &v);
    int setTrialStrainIncr (const Vector &v, const Vector &r);
    
	const Vector &getStress (void);
    const Vector &getStrain (void);   
	const Matrix &getTangent (void);
    const Matrix &getInitialTangent (void);

    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);

  virtual const char *getType (void) const ;
  virtual int getOrder (void) const ;

   void Print(OPS_Stream &s, int flag = 0);
   int sendSelf(int commitTag, Channel &theChannel);
   int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    Response *setResponse (const char **argv, int argc, 
				   OPS_Stream &s);
    int getResponse (int responseID, Information &matInformation);


  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int    setParameter             (const char **argv, int argc, Parameter &param);
  int    updateParameter          (int parameterID, Information &info);
  int    activateParameter        (int parameterID);
  const Vector  &getStressSensitivity   (int gradNumber, bool conditional);
  const Matrix  &getInitialTangentSensitivity (int gradNumber);
    int commitSensitivity        (const Vector & strainGradient, int gradNumber, int numGrads);


  // AddingSensitivity:END ///////////////////////////////////////////



private:
 
	//defining parameters
	double fcu;		//concrete conmpressive strength
	double ecu;		//concrete compress strain at ultimate strength
	double Ec;      //concrete modulus of elasticity
	double fcr;		//concrete strength at cracking
	double Esv;		//vertical reinforcement Youngs Modulus
	double fyv;     //vertical reinforcement yield stress
	double alphaV;	//vertical steel reinforcement hardening ratio
	double RoV;		//vertical steel reinforcement ratio to concrete area
	double nE;
	//Internal Parameters
	//Principle strains
	 Vector epsf;
	 Vector sigf;

	 Vector sigfsens;
	 Vector sigfsensP;

	double ex;
	double exy;
	double e1;
	double e2;
	double ey;
	double theta;
	double sig1;
	double sig2;
	double exmin;
	double exmax;
	double exymin;
	double exymax;
	double eymax;
	double eymin;
	double exminc;
	double exmaxc;
	double exyminc;
	double exymaxc;
	double eymaxc;
	double eyminc;

	double crackLabelP;
	double e1max;
	double e2min;
	double fe1max;
	double fe2min;
	double e1P;
	double e2P;
	double thetaP;

	double exP;
	double eyP;
	double exyP;
	double fxP;
	double fxyP;
	double loadpath;
	double loadpathP;

	double exminP;
	double exmaxP;
	double eyminP;
	double eymaxP;
	double exyminP;
	double exymaxP;
	double e1maxP;
	double e2minP;

	int parameterID;
	Matrix Dr;
	Matrix Dri;

	//Matrix ComSensMat;
	Matrix DrP;

    double fx;
	double fxy;
	double fy;

	//MCFT Crack check: Tensile stresses
	double InitCrackAngle;
	double FinalAnglex;
	double epsy;
	double crackLabel;
	double Strain1;
	double Strain2;
	double Sigma1;
	double Sigma2;

    double sensParaID;
    Matrix *SHVs;
  
	void ForwardAngleSearch(void);
	void ForwardAngleSearch01(void);
	void ForwardAngleSearch89(void);
	
	void StressEnvelope(double e1, double e2, double e1P, double e2P, double e1max, double e2min);
	void Loadf(void);

	int searchDirection;

	double tangentstifness00(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esvvv, double RoV, double e1P, double e2P, double fe1max, double e1max, double fe2min, double e2min);
	double tangentstifness01(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esvvv, double RoV, double e1P, double e2P, double fe1max, double e1max, double fe2min, double e2min);
	double tangentstifness10(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esvvv, double RoV, double e1P, double e2P, double fe1max, double e1max, double fe2min, double e2min);
	double tangentstifness11(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esvvv, double RoV, double e1P, double e2P, double fe1max, double e1max, double fe2min, double e2min);
	
};

#endif
