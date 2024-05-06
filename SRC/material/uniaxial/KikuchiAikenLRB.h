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

// $Revision: 1.0 $
// $Date: 2012-08-09 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/KikuchiAikenLRB.h,v $

// Written: Masaru Kikuchi
// Created: August 9, 2012
//
// Kikuchi&Aiken model for lead rubber bearing
//
// Description: This file contains the class definition for KikuchiAikenLRB.
//
// (input)                             (output)
// deformation -> strain -> stress  -> force
//                strain -> modulus -> spring constant
//

#ifndef KikuchiAikenLRB_h
#define KikuchiAikenLRB_h

#include <UniaxialMaterial.h>


class KikuchiAikenLRB : public UniaxialMaterial
{
 public:
  KikuchiAikenLRB(int tag, int type, double ar, double hr, double gr, double ap, double tp,
		  double alph, double beta, double temp, double rk, double rq, double rs, double rf);

  KikuchiAikenLRB();
  ~KikuchiAikenLRB();

  const char *getClassType(void) const {return "KikuchiAikenLRB";};

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

  void Print(OPS_Stream &s, int flag =0);

 protected:

 private:

  void setType(int Tp);
  
  // input values
  int    Type; // type  =1:LRB
  double Ar;   // area (rubber) [m^2]
  double Hr;   // total rubber thickness [m]
  double Gr;   // shear modulus (rubber) [N/m^2]
  double Ap;   // area (lead plug) [m^2]
  double Tp;   // yield stress (lead plug) [N/m^2]
  double Alph; // shear modulus (lead plug) [N/m^2]
  double Beta; // stiffness rate (initial/yield)
  double Temp; // temperature [C]
  double Rk;   // coefficient for Kd
  double Rq;   // coefficient for Qd
  double Rs;   // reduction rate of stiffness (for MSS model)
  double Rf;   // reduction rate of force (for MSS model)

  //
  double qd100; // reference value for Qd
  double kd100; // reference value for Kd
  double ku100; // reference value for Ku
  double qd;    // Qd (force at zero displacement)
  double kd;    // Kd (yielding stiffness)
  double ku;    // Ku (elastic stiffness)

  //
  double trgStrain; // elastic limit strain
  double lmtStrain; // application limit strain
  double initialStiff;   //initial stiffness

  // temporary values for Kikuchi&Aiken model
  double tmpStrain;
  double tmpDeform;
  double keq; //keq
  double heq; //heq
  double u;   //u
  double n;   //n
  double p;   //p
  double a;   //a
  double b;   //b
  double c;   //c
  double xm;  //maximum strain
  double fm;  //maximum stress
  double x;   //normarized strain
  double alpha; //
  double q1Stf; //tangent of Q1
  double q2Stf; //tangent of Q2

  // trial values
  double trialDeform;  //deformation (input)
  double trialForce;   //force (output)
  double trialStiff;   //stiffness (output)
  double trialStrain;       //strain
  bool   trialIfElastic;    //flag
  double trialQ1;           //stress Q1 component (nonlinear elastic)
  double trialQ2;           //stress Q2 component (hysteretic)
  double trialMaxStrain;       //maximum (absolute)
  double trialDDeform;         //incremental strain
  int    trialDDeformLastSign; //+1 (if DDeform > 0), -1(if DDeform < 0)
  int    trialIdxRev;          //index of reversal point


  // commit values
  double commitDeform;
  double commitForce;
  double commitStiff;
  double commitStrain;
  bool   commitIfElastic;
  double commitQ1;
  double commitQ2;
  double commitMaxStrain;
  double commitDDeform;
  int    commitDDeformLastSign;
  int    commitIdxRev;



  //memorize reversal points of hysteretic curve
  // XX[0] skeleton curve
  // XX[1] unload
  // XX[2,3,...] reversal load
  int numIdx; //array size of reversal points
  double *revXBgn;  //x  (reversal point)
  double *revQ2Bgn; //q2 (reversal point)
  double *revXEnd;  //x  (directing point)
  double *revQ2End; //q2 (directing point)
  double *revB;     //b
  double *revAlpha; //alpha


  //formulae of Kikuchi&Aiken model
  double compQ1(double u, double n, double p, double fm, double x);
  double compQ2Unload(double u, double a, double b, double c, double fm, double x);
  double compQ2Masing(double u, double a, double b, double c, double fm, double x1, double x2, double q2i, double alpha);
  double compAlpha(double a, double b1, double b2, double c, double x1, double x2, double alpha0);

  double compQ1Derivertive(double u, double n, double p, double keq, double x);
  double compQ2UnloadDerivertive(double u, double a, double b, double c, double keq, double x);
  double compQ2MasingDerivertive(double u, double a, double b, double c, double keq, double x1, double x2, double alpha);

  //bisection method to calculate parameter "a"
  static double compABisection(double heq, double u, double min, double max, double lim, double tol);

  //Keq, Heq
  static double compKeq(double xm, double qd, double kd);
  static double compHeq(double xm, double qd, double kd, double ku);

  //formulae to calculate parameters
  double (*calcN)(double gm);
  double (*calcP)(double gm);
  double (*calcA)(double gm, double heq, double u);
  double (*calcB)(double gm, double a, double c,double heq, double u);
  double (*calcC)(double gm);
  double (*calcCQd)(double gm);
  double (*calcCKd)(double gm);
  double (*calcCHeq)(double gm);


  //parameters for each rubber
  //LRB
  static double calcNType1(double gm);
  static double calcPType1(double gm);
  static double calcAType1(double gm, double heq, double u);
  static double calcBType1(double gm, double a, double c,double heq, double u);
  static double calcCType1(double gm);
  static double calcCQdType1(double gm);
  static double calcCKdType1(double gm);
  static double calcCHeqType1(double gm);

};

#endif

