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
// $Date: 2013-06-03 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/AxialSpHD.h,v $

// Written: Kazuki Tanimoto
// Created: June 2012
//
// Top and bottom series of axial springs for three-dimensonal multi-spring mechanical model.
// This model is considering hardening at tension side.
//
// Description: This file contains the class definition for AxialSpHD.
//
//   strain  --> disp.
//   stress  --> force
//   modulus --> stiffness
//

#ifndef AxialSpHD_h
#define AxialSpHD_h

#include <UniaxialMaterial.h>

class AxialSpHD : public UniaxialMaterial
{
 public:
  AxialSpHD(int tag, double sce, double fty, double fcy, double bte,
	  double bty, double bth, double bcy, double fcr, double ath);

  AxialSpHD();

  ~AxialSpHD();

  const char *getClassType(void) const {return "AxialSpHD";};

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

  double trialDeformation;   // trial strain
  double trialForce;         // trial stress
  double trialStiffness;     // trial tangent
  double commitDeformation;  // commit strain
  double commitForce;        // commit stress
  double commitStiffness;    // commit tangent

  // input
  double sce;   // compressive modulus
  double fty;   // tensile yield stress (0.0 < fty)
  double fcy;   // compressive yield stress (fcy < 0.0)
  double bte;   // reduction rate for tensile modulus (ratio to sce, 0.0 <= bty <= bth <= bte <= 1.0)
  double bty;   // reduction rate after tensile yield point (ratio to sce)
  double bth;   // reduction rate after tensile hardening point (ratio to sce)
  double bcy;   // reduction rate after compressive yield point (ratio to sce, 0.0 <= bcy <= 1.0)
  double fcr;   // reversal point stress (fcy <= fcr <= 0.0)
  double ath;   // tensile hardening point strain (ratio to tensile yield strain, 1.0 <= ath)

  // data for calculation
  double ste;   // tensile modulus
  double sty;   // tensile modulus after yield point
  double stp;
  double sth;   // tensile modulus after hardening point
  double scy;   // compressive modulus after yield point
  double uty;   // tensile yield strain
  double ucy;   // compressive yield strain
  double ucr;   // reversal point strain
  double utr;   // strain (related to reversal point)
  double ftr;   // stress (related to reversal point)
  double uth;   // hardening point strain
  double fth;   // hardening point stress
  double uch;   // strain (related to hardening point)
  double fch;   // stress (related to hardening point)

  double uc0;   // strain (related to compressive yield point)
  double ur1;   // ref. point no.1
  double fr1;   // ref. point no.1
  double ur2;   // ref. point no.2
  double fr2;   // ref. point no.2
  double ur3;   // ref. point no.3
  double fr3;   // ref. point no.3
  double ur4;   // ref. point no.4
  double fr4;   // ref. point no.4
  double ur5;   // ref. point no.5
  double fr5;   // ref. point no.5
  double ur6;   // ref. point no.6
  double fr6;   // ref. point no.6
  double ur7;   // ref. point no.7
  double fr7;   // ref. point no.7

  int trialStg; // index for Stg
  int commitStg;//

};

#endif

