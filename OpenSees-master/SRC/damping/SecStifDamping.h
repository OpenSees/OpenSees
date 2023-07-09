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
// $Source: /usr/local/cvs/OpenSees/SRC/damping/SecStifDamping.h,v $

// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 02/2020
// Revision: A
//
// Description: This file contains the definition for the SecStifDamping class.
// FRDamping provides the abstraction of an elemental damping imposition
// providing quasi tangent stiffness proportional damping
//
// What: "@(#) FRDamping.h, revA"

#ifndef SecStifDamping_h
#define SecStifDamping_h

#include <Damping.h>
#include <Vector.h>
#include <TimeSeries.h>

class SecStifDamping: public Damping
{
public:
  SecStifDamping(int tag, double beta, double ta, double td, TimeSeries *fac);
  
  SecStifDamping();
  ~SecStifDamping();
  
  const char *getClassType() const {return "SecStifDamping";};
  
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
  double beta, ta, td;
  TimeSeries *fac;
  Vector *qd, *q0, *q0C;
  Domain *theDomain;
};

#endif
