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
// $Source: /usr/local/cvs/OpenSees/SRC/damping/UniformDamping.cpp,v $

// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 02/2020
// Revision: A
// 
// Description: This file contains the implementation for the UniformDamping class.
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
// https://doi.org/10.1016/j.compstruc.2018.10.016.
// (http://www.sciencedirect.com/science/article/pii/S0045794918300555)

// What: "@(#) UniformDamping.cpp, revA"

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <Domain.h>
#include <Channel.h>
#include <elementAPI.h>
#include <UniformDamping.h>

// constructor:
UniformDamping::UniformDamping(int tag, double cd, double f1, double f2):
Damping(tag, DMP_TAG_UniformDamping),
nComp(0), nFilter(0),
eta(cd), freq1(f1), freq2(f2),
alpha(0), omegac(0), qL(0), qLC(0), qd(0), qdC(0), q0(0), q0C(0)
{
  if (eta <= 0.0) opserr << "UniformDamping::UniformDamping:  Invalid damping ratio\n";
  if (freq1 <= 0.0 || freq2 <= 0.0 || freq1 >= freq2)
    opserr << "UniformDamping::UniformDamping:  Invalid frequency range\n";

  Initialize();
}

// constructor:
UniformDamping::UniformDamping(int tag, double cd, double f1, double f2, int nF, Vector *a, Vector *w):
Damping(tag, DMP_TAG_UniformDamping),
nComp(0), nFilter(0),
eta(cd), freq1(f1), freq2(f2),
alpha(0), omegac(0), qL(0), qLC(0), qd(0), qdC(0), q0(0), q0C(0)
{
  if (eta <= 0.0) opserr << "UniformDamping::UniformDamping:  Invalid damping ratio\n";
  if (freq1 <= 0.0 || freq2 <= 0.0 || freq1 >= freq2)
    opserr << "UniformDamping::UniformDamping:  Invalid frequency range\n";
  if (nF > 0 && a->Size() == nF && w->Size() == nF)
  {
    nFilter = nF;
    alpha = new Vector(*a);
    omegac = new Vector(*w);
  }
  else
  {
    Initialize();
  }
}


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
UniformDamping::UniformDamping():
Damping(0, DMP_TAG_UniformDamping),
nComp(0), nFilter(0),
eta(0.0), freq1(0.0), freq2(0.0),
alpha(0), omegac(0), qL(0), qLC(0), qd(0), qdC(0), q0(0), q0C(0)
{

}


// destructor:
UniformDamping::~UniformDamping() 
{
  if (alpha) delete alpha;
  if (omegac) delete omegac;
  if (qL) delete qL;
  if (qLC) delete qLC;
  if (qd) delete qd;
  if (qdC) delete qdC;
  if (q0) delete q0;
  if (q0C) delete q0C;
}

int
UniformDamping::Initialize(void)
{
  double f1log = log10(freq1);
  double f2log = log10(freq2);
  nFilter = static_cast<int>(floor(f2log - f1log)) + 2;
  double dfreq = (f2log - f1log) / (nFilter - 1);
  alpha = new Vector(nFilter);
  omegac = new Vector(nFilter);

  for (int i = 0; i < nFilter; ++i)
  {
    (*omegac)(i) = 6.28318530718 * pow(10.0, f1log + i * dfreq);
  }

  int nf = 100 * nFilter;
  double df = (f2log - f1log) / (nf - 1);
  Vector *y = new Vector(nFilter);
  Matrix *X = new Matrix(nFilter, nFilter);
  for (int i = 0; i < nf; ++i)
  {
    double omega = 6.28318530718 * pow(10.0, f1log + i * df);
    for (int j = 0; j < nFilter; ++j)
    {
      double wjn = omega / (*omegac)[j];
      double phij = 2.0 * wjn / (1.0 + wjn * wjn);
      (*y)(j) += phij;
      for (int k = 0; k < nFilter; ++k)
      {
        double wkn = omega / (*omegac)[k];
        double phik = 2.0 * wkn / (1.0 + wkn * wkn);
        (*X)(j, k) += phij * phik;
      }
    }
  }

  *alpha = (*y) / (*X);

  delete y;
  delete X;
  
  return 0;
}

int
UniformDamping::commitState(void)
{
  *qdC = *qd;
  *q0C = *q0;
  *qLC = *qL;
  return 0;
}


int
UniformDamping::revertToLastCommit(void)
{
  *qd = *qdC;
  *q0 = *q0C;
  *qL = *qLC;
  return 0;
}


int
UniformDamping::revertToStart(void)
{
  (*qd).Zero();
  (*qdC).Zero();
  (*q0).Zero();
  (*q0C).Zero();
  (*qL).Zero();
  (*qLC).Zero();
  return 0;
}


int 
UniformDamping::setDomain(Domain *domain, int nC)
{
  theDomain = domain;
  nComp = nC;
  
  qd = new Vector(nComp);
  qdC = new Vector(nComp);
  q0 = new Vector(nComp);
  q0C = new Vector(nComp);
  qL = new Matrix(nComp, nFilter);
  qLC = new Matrix(nComp, nFilter);
  
  return 0;
}


int
UniformDamping::update(Vector q)
{
  double dT = theDomain->getDT();
  if (dT > 0.0)
  {
    *q0 = q;
    (*qd).Zero();
    for (int i = 0; i < nFilter; ++i)
    {
      double dTomegac = dT * (*omegac)(i);
      double cqd = 2.0 * (*alpha)(i) * eta / (2.0 + dTomegac);
      double cqL1 = dTomegac / (2.0 + dTomegac);
      double cqL2 = 1.0 - 2.0 * cqL1;
      for (int j = 0; j < nComp; ++j)
      {
        (*qd)(j) += cqd * ((*q0C)(j) + (*q0)(j) - 2.0 * (*qLC)(j,i));
        (*qL)(j,i) = cqL1 * ((*q0C)(j) + (*q0)(j)) + cqL2 * (*qLC)(j,i);
      }
    }
    *qd = *qd * 2.0 - *qdC;
  }
  return 0;
}

const Vector &
UniformDamping::getDampingForce(void)
{
  return (*qd);
}

double UniformDamping::getStiffnessMultiplier(void)
{
  double dT = theDomain->getDT();
  double km = 1.0;
  if (dT > 0.0)
  {
    for (int i = 0; i < nFilter; ++i)
    {
      km += 4.0 * (*alpha)(i) * eta / (2.0 + (*omegac)(i) * dT);
    }
  }
  return km;
}

Damping *UniformDamping::getCopy(void)
{
  // create a new instance of UniformDamping 

  UniformDamping *theCopy;

  theCopy = new UniformDamping(this->getTag(), eta, freq1, freq2, nFilter, alpha, omegac);

  return theCopy;
}


int 
UniformDamping::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(4);
  data(3) = this->getTag();
  data(0) = eta;
  data(1) = freq1;
  data(2) = freq2;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)
  {
    opserr << " UniformDamping::sendSelf() - data could not be sent\n" ;
    return -1;
  }
  return 0;
}


int 
UniformDamping::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static Vector data(4);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)
  {
    opserr << " UniformDamping::recvSelf() - data could not be received\n" ;
    return -1;
  }
    
  this->setTag((int)data(3));
  eta = data(0);
  freq1 = data(1);
  freq2 = data(2);

  if (eta <= 0.0) opserr << "UniformDamping::recvSelf:  Invalid damping ratio\n";
  if (freq1 <= 0.0 || freq2 <= 0.0 || freq1 >= freq2)
    opserr << "UniformDamping::recvSelf:  Invalid frequency range\n";
  
  Initialize();
  return 0;
}

void
UniformDamping::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE)
  {
    s << "\nDamping: " << this->getTag() << " Type: UniformDamping";
    s << "\tdamping ratio: " << 0.5 * eta << endln;
    s << "\tlower bound frequency: " << freq1 << endln;
    s << "\tupper bound frequency: " << freq2 << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON)
  {
    s << "\t\t\t{\"name\": \"" << this->getTag() << "\", \"type\": \"UniformDamping\"";
    s << ", \"damping ratio\": [" << 0.5 * eta << "]";
    s << ", \"lower bound frequency\": [" << freq1 << "]";
    s << ", \"upper bound frequency\": [" << freq2 << "]";
    s << "}";
  }
}
