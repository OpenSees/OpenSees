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
// $Date: 2021-07-02 14:29:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/damping/URDDamping.h,v $

// Revised by: Y Tian
// Created: 02/2020
// Revision: A
//
// Description: This file contains the definition for the URDDamping class.
// URDDamping provides the abstraction of an elemental damping imposition
// providing user-define damping over a frequency range
//
// Reference:
// 

// What: "@(#) URDDamping.h, revA"

#ifndef URDDamping_h
#define URDDamping_h

#include <Damping.h>
#include <Vector.h>
#include <TimeSeries.h>

class URDDamping: public Damping
{
public:
  URDDamping(int tag, int numfreq, Matrix *etaFreq, double dptol, double ta, double td, TimeSeries *fac, int prttag, int maxiter);
  URDDamping(int tag, int numfreq, Matrix* etaFreq, double dptol, double ta, double td, TimeSeries *fac, int nFilter, Vector *alpha, Vector *omegac, Vector *omegaetaf, int prttag, int maxiter);
  
  URDDamping();
  ~URDDamping();
  
  const char *getClassType() const {return "URDDamping";};

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
  int numfreq,prttag,maxiter;
  int nComp, nFilter;
  double ta, td, dptol;
  TimeSeries *fac;
  Vector *alpha, *omegac, *omegaetaf;
  Vector *Freqlog, *Fredif, *Freqk, *Freqb; 
  Matrix *etaFreq;
  Matrix *qL, *qLC;
  Vector *qd, *qdC, *q0, *q0C;
  Domain *theDomain;
};

#endif
