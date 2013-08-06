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
// $Date: 2012-06-04 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/AxialSp.h,v $

// Written: Kazuki Tanimoto
// Created: June 2012
//
// Top and bottom series of axial springs for three-dimensonal multi-spring mechanical model.
//
// Description: This file contains the class definition for AxialSp.
//
//   strain  --> disp.
//   stress  --> force
//   modulus --> stiffness
//

#ifndef AxialSp_h
#define AxialSp_h

#include <UniaxialMaterial.h>

class AxialSp : public UniaxialMaterial
{
 public:
  AxialSp(int tag, double sce, double fty, double fcy,
	  double bte, double bty, double bcy, double fcr);

  AxialSp();

  ~AxialSp();

  const char *getClassType(void) const {return "AxialSp";};

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
  double trialForce;   // trial stress
  double trialStiffness;  // trial tangent
  double commitDeformation;  // commit strain
  double commitForce;  // commit stress
  double commitStiffness; // commit tangent

  // input
  double sce;
  double fty;
  double fcy;
  double bte;
  double bty;
  double bcy;
  double fcr;

  // data for calculation
  double ste;
  double sty;
  double scy;
  double uty;
  double ucy;
  double ucr;

  double uc0;
  double ur1;
  double fr1;
  double ur2;
  double fr2;
  double ur3;
  double fr3;
  double ur4;
  double fr4;
  double ur5;
  double fr5;

  int trialStg;
  int commitStg;

};

#endif

