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
// $Date: 2019-01-28 17:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/damping/UniformDamping.h,v $

// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 02/2020
// Revision: A
//
// Description: This file contains the definition for the UniformDamping class.
// UniformDamping provides the abstraction of an elemental damping imposition
// providing uniform damping over a frequency range
//
// Reference:
// Yuli Huang, Richard Sturt, Michael Willford,
// A damping model for nonlinear dynamic analysis providing uniform damping over a frequency range,
// Computers & Structures,
// Volume 212,
// 2019,
// Pages 101-109,
// ISSN 0045-7949,
// https://doi.org/10.1016/j.compstruc.2018.10.016
//
// Yuan Tian, Yuli Huang, Zhe Qu, Yifan Fei, Xinzheng Lu,
// High-performance uniform damping model for response history analysis in OpenSees,
// Journal of Earthquake Engineering,
// 2022,
// https://doi.org/10.1080/13632469.2022.2124557
//
// What: "@(#) UniformDamping.h, revA"

#ifndef UniformDamping_h
#define UniformDamping_h

#include <Damping.h>
#include <Vector.h>
#include <TimeSeries.h>

class UniformDamping: public Damping
{
public:
  UniformDamping(int tag, double eta, double freq1, double freq2, double ta, double td, TimeSeries *fac);
  UniformDamping(int tag, double eta, double freq1, double freq2, double ta, double td, TimeSeries *fac, int nFilter, Vector *alpha, Vector *omegac);
  
  UniformDamping();
  ~UniformDamping();
  
  const char *getClassType() const {return "UniformDamping";};

  int Initialize(void);
  
  int setDomain(Domain *domain, int nComp);
  int update(Vector q);
  
  int commitState(void);
  int revertToLastCommit(void);    
  int revertToStart(void);
  
  const Vector &getDampingForce(void);
  double getStiffnessMultiplier(void);
  
  Damping *getCopy(void);
  
  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  
  void Print(OPS_Stream &s, int flag = 0);
  
private:
  
  // internal data
  int nComp, nFilter;
  double eta, freq1, freq2, ta, td;
  TimeSeries *fac;
  Vector *alpha, *omegac;
  Matrix *qL, *qLC;
  Vector *qd, *qdC, *q0, *q0C;
  Domain *theDomain;
};

#endif
