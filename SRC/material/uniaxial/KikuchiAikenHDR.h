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
// $Date: 2013-05-29 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/KikuchiAikenHDR.h,v $

// Written: Ken Ishii
// Created: June 2012
//
// Jan 31, 2017: mkiku
//    X0.6 standard and zero vertical stress Geq eqs are modified
//    New compound X0.4 and X0.3
//    compABisection is fixed
//
// Kikuchi&Aiken model for high-damping rubber bearing
//
// Description: This file contains the class definition for KikuchiAikenHDR.
//
// (input)                             (output)
// deformation -> strain -> stress  -> force
//                strain -> modulus -> spring constant
//

#ifndef KikuchiAikenHDR_h
#define KikuchiAikenHDR_h

#include <UniaxialMaterial.h>

class KikuchiAikenHDR : public UniaxialMaterial
{
 public:
  KikuchiAikenHDR(int tag, int tp, double ar, double hr,
		  double cg, double ch, double cu, double rs, double rf);
  KikuchiAikenHDR();
  ~KikuchiAikenHDR();

  const char *getClassType(void) const {return "KikuchiAikenHDR";};

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
  int    Tp; // type = 1;X0.6(standard compressive stress), 2;X0.6(zero compressive stress), 3;X0.4(standard), 4;X0.4(zero), 5;X0.3(standard), 6;X0.3(zero)
  double Ar; // area [m^2]
  double Hr; // total rubber thickness [m]
  double Cg; // correction coefficient of Geq
  double Ch; // correction coefficient of Heq
  double Cu; // correction coefficient of U
  double Rs; // reduction rate of stiffness (used for MSS model)
  double Rf; // reduction rate of force (used for MSS model)

  //
  double trgStrain; // elastic limit strain
  double lmtStrain; // application limit strain
  double initialStiff;   // initial stiffness

  // temporary values for Kikuchi&Aiken model
  //double tmpStrain; //
  double geq; //Geq
  //double heq; //heq
  double u;   //u
  double n;   //n
  double a;   //a
  //double b;   //b
  double c;   //c
  double xm;  //maximum strain
  double fm;  //maximum stress
  //double x;   //normarized strain
  //double alpha; //
  //double q1Tan; //tangent of Q1
  //double q2Tan; //tangent of Q2


  // trial values
  double trialDeform; // deformation (input)
  double trialForce;  // force (output)
  double trialStiff;  // stiffness (output)
  double trialStrain;    //strain
  double trialStress;    //stress
  double trialTangent;   //tangent
  bool   trialIfElastic; //flag
  double trialQ1; //stress Q1 component (nonlinear elastic)
  double trialQ2; //stress Q2 component (hysteretic)
  double trialMaxStrain;       //maximum strain (absolute)
  double trialDStrain;         //incremental strain
  int    trialDStrainLastSign; //+1(if Dstrain > 0), -1(if Dstrain < 0)
  int    trialIdxRev;  //index of reversal point


  // commit values
  double commitDeform;
  double commitForce;
  double commitStiff;
  double commitStrain;
  double commitStress;
  double commitTangent;
  bool   commitIfElastic;
  double commitQ1;
  double commitQ2;
  double commitMaxStrain;
  double commitDStrain;
  int    commitDStrainLastSign;
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
  double compQ1(double u, double n, double fm, double x);
  double compQ2Unload(double u, double a, double b, double c, double fm, double x);
  double compQ2Masing(double u, double a, double b, double c, double fm, double x1, double x2, double q2i, double alpha);
  double compAlpha(double a, double b1, double b2, double c, double x1, double x2, double alpha0);

  double compQ1Derivertive(double u, double n, double geq, double x);
  double compQ2UnloadDerivertive(double u, double a, double b, double c, double geq, double x);
  double compQ2MasingDerivertive(double u, double a, double b, double c, double geq, double x1, double x2, double alpha);

  //bisection method to calculate parameter "a"
  static double compABisection(double heq, double u, double min, double max, double lim, double tol);


  //formulae to calculate parameters
  double (*calcGeq)(double gm);
  double (*calcHeq)(double gm);
  double (*calcU)(double gm);
  double (*calcN)(double gm);
  double (*calcA)(double gm, double heq, double u);
  double (*calcB)(double gm, double a, double c,double heq, double u);
  double (*calcC)(double gm);


  //parameters for each rubber
  //X0.6 standard compressive stress
  static double calcGeqTp1(double gm);
  static double calcHeqTp1(double gm);
  static double calcUTp1(double gm);
  static double calcNTp1(double gm);
  static double calcATp1(double gm, double heq, double u);
  static double calcBTp1(double gm, double a, double c,double heq, double u);
  static double calcCTp1(double gm);
  //X0.6 zero compressive stress
  static double calcGeqTp2(double gm);
  static double calcHeqTp2(double gm);
  static double calcUTp2(double gm);
  static double calcNTp2(double gm);
  static double calcATp2(double gm, double heq, double u);
  static double calcBTp2(double gm, double a, double c,double heq, double u);
  static double calcCTp2(double gm);

  //X0.4 standard compressive stress
  static double calcGeqTp3(double gm);
  static double calcHeqTp3(double gm);
  static double calcUTp3(double gm);
  static double calcNTp3(double gm);
  static double calcATp3(double gm, double heq, double u);
  static double calcBTp3(double gm, double a, double c,double heq, double u);
  static double calcCTp3(double gm);
  //X0.4 zero compressive stress
  static double calcGeqTp4(double gm);
  static double calcHeqTp4(double gm);
  static double calcUTp4(double gm);
  static double calcNTp4(double gm);
  static double calcATp4(double gm, double heq, double u);
  static double calcBTp4(double gm, double a, double c,double heq, double u);
  static double calcCTp4(double gm);

  //X0.3 standard compressive stress
  static double calcGeqTp5(double gm);
  static double calcHeqTp5(double gm);
  static double calcUTp5(double gm);
  static double calcNTp5(double gm);
  static double calcATp5(double gm, double heq, double u);
  static double calcBTp5(double gm, double a, double c,double heq, double u);
  static double calcCTp5(double gm);
  //X0.3 zero compressive stress
  static double calcGeqTp6(double gm);
  static double calcHeqTp6(double gm);
  static double calcUTp6(double gm);
  static double calcNTp6(double gm);
  static double calcATp6(double gm, double heq, double u);
  static double calcBTp6(double gm, double a, double c,double heq, double u);
  static double calcCTp6(double gm);

};

#endif

