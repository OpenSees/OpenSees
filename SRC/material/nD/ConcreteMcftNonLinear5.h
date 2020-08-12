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


#ifndef ConcreteMcftNonLinear5_h
#define ConcreteMcftNonLinear5_h


class ConcreteMcftNonLinear5 : public NDMaterial
{


public:
  
 //constructor
  ConcreteMcftNonLinear5(int tag, double fcu, double ecu, double Ec, double fcr, double Esv, double fyv, double alphaV, double RoV);
  ConcreteMcftNonLinear5(); 
  ~ConcreteMcftNonLinear5();

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
	
	//Internal Parameters
	//Principle strains
	 Vector epsf;
	 Vector sigf;
	 Vector dsigfdfcu;
	 Vector dsigfdRoV;

	double e1;
	double e2;
	double ey;
	double radius;
	double sig1;
	double sig2;

	double exP;
	double exyP;
	double fxP;
	double fxyP;
///sens
	double dfxdfcuP;
	double dfxydfcuP;
	double dfxdRoVP;
	double dfxydRoVP;
///
	int parameterID;

	Matrix Dr;
	Matrix dDri;
	Matrix TT;
	Matrix dDdfcu;      
    Matrix dDdRoV;


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
	double t00;
	double t01;
	double t10;
	double t11;
    double sensParaID;




	int searchDirection;

	double c1tmd00(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c1tmd01(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c1tmd10(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c1tmd11(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);

	double c2tmd00(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c2tmd01(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c2tmd10(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c2tmd11(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);

	//fcu sensitivity of tangent stiffness
	double c1dd00dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c1dd01dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c1dd10dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c1dd11dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);

	double c2dd00dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c2dd01dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c2dd10dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c2dd11dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);

	//RoV sensitivity of tangent stiffness
	double c1dd00dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c1dd01dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c1dd10dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c1dd11dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);

	double c2dd00dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c2dd01dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c2dd10dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);
	double c2dd11dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV);



	double tm3(double e2, double ecu, double fcu, double nE);
};

#endif
