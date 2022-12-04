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
// $Source: /usr/local/cvs/OpenSees/SRC/damping/URDDampingbeta.h,v $

// Revised by: Y Tian
// Created: 02/2020
// Revision: A
//
// Description: This file contains the definition for the URDDampingbeta class.
// URDDampingbeta provides the abstraction of an elemental damping imposition
// providing user-define damping over a frequency range
//
// Reference:
//

// What: "@(#) URDDampingbeta.cpp, revA"

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <Domain.h>
#include <Channel.h>
#include <elementAPI.h>
#include <URDDampingbeta.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>
//#include <algorithm>
//#include <iostream>
//#include <fstream>
//#include <iomanip>
//using namespace std;

//extern StaticAnalysis *theStaticAnalysis;

// constructor:
URDDampingbeta::URDDampingbeta(int tag, int nfreq, Vector *tmpomegac, Vector *tmpbeta, double t1, double t2, TimeSeries *f):
Damping(tag, DMP_TAG_URDDampingbeta),
nComp(0), nFilter(nfreq),
ta(t1), td(t2), fac(f),
qL(0), qLC(0), qd(0), qdC(0), q0(0), q0C(0)
{
  beta = new Vector(*tmpbeta);
  omegac = new Vector(*tmpomegac);
}



// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
URDDampingbeta::URDDampingbeta():
Damping(0, DMP_TAG_URDDampingbeta),
nComp(0), nFilter(0),
beta(0), omegac(0), ta(0.0), td(0.0), fac(0),
qL(0), qLC(0), qd(0), qdC(0), q0(0), q0C(0)
{

}


// destructor:
URDDampingbeta::~URDDampingbeta() 
{
  if (fac) delete fac;
  if (beta) delete beta;
  if (omegac) delete omegac;
  if (qL) delete qL;
  if (qLC) delete qLC;
  if (qd) delete qd;
  if (qdC) delete qdC;
  if (q0) delete q0;
  if (q0C) delete q0C;
}

int
URDDampingbeta::Initialize(void)
{
    return 0;
}

int
URDDampingbeta::commitState(void)
{
  *qdC = *qd;
  *q0C = *q0;
  *qLC = *qL;
  return 0;
}


int
URDDampingbeta::revertToLastCommit(void)
{
  *qd = *qdC;
  *q0 = *q0C;
  *qL = *qLC;
  return 0;
}


int
URDDampingbeta::revertToStart(void)
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
URDDampingbeta::setDomain(Domain *domain, int nC)
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
URDDampingbeta::update(Vector q)
{
  double t = theDomain->getCurrentTime();
  double dT = theDomain->getDT();
  StaticAnalysis **theStaticAnalysis = OPS_GetStaticAnalysis();  
  if (*theStaticAnalysis)
  {
    *q0 = q;
    (*qd).Zero();
    for (int i = 0; i < nFilter; ++i)
      for (int j = 0; j < nComp; ++j)
        (*qL)(j,i) = q(j);
  }
  else if (dT > 0.0)
  {
    *q0 = q;
    (*qd).Zero();
    if (t < td)
    {
      if (t > ta)
      {
        for (int i = 0; i < nFilter; ++i)
        {
          double dTomegac = dT * (*omegac)(i);
          double cd = 4.0 * (*beta)(i) / (2.0 + dTomegac);
          double c0 = dTomegac / (2.0 + dTomegac);
          double cL = (2.0 - dTomegac) / (2.0 + dTomegac);
          for (int j = 0; j < nComp; ++j)
          {
            (*qd)(j) += cd * ((*q0C)(j) + (*q0)(j) - 2.0 * (*qLC)(j,i));
            (*qL)(j,i) = c0 * ((*q0C)(j) + (*q0)(j)) + cL * (*qLC)(j,i);
          }
        }
        *qd -= *qdC;
      }
      else
      {
        for (int i = 0; i < nFilter; ++i)
          for (int j = 0; j < nComp; ++j)
            (*qL)(j,i) = q(j);
      }
      if (fac) *qd *= fac->getFactor(t);
    }
  }
  return 0;
}

const Vector &
URDDampingbeta::getDampingForce(void)
{
  return (*qd);
}

double URDDampingbeta::getStiffnessMultiplier(void)
{
  double t = theDomain->getCurrentTime();
  double dT = theDomain->getDT();
  double km = 0.0;
  StaticAnalysis **theStaticAnalysis = OPS_GetStaticAnalysis();
  if (!*theStaticAnalysis && dT > 0.0 && t > ta && t < td)
  {
    for (int i = 0; i < nFilter; ++i)
    {
        km += 4.0 * (*beta)(i) / (2.0 + (*omegac)(i) * dT); 
    }
    if (fac) km *= fac->getFactor(t);
  }
  return 1.0 + km;
}

Damping *URDDampingbeta::getCopy(void)
{
  // create a new instance of URDDampingbeta 

  URDDampingbeta *theCopy;

  theCopy = new URDDampingbeta(this->getTag(), nFilter, omegac, beta, ta, td, fac);

  return theCopy;
}


int 
URDDampingbeta::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  static ID idData(2);
  static Vector data(4);
  static Vector dataomegac(nFilter);
  static Vector databeta(nFilter);

  if (fac)
  {
    idData(0) = fac->getClassTag();
    int seriesDbTag = fac->getDbTag();
    if (seriesDbTag == 0)
    {
      seriesDbTag = theChannel.getDbTag();
      fac->setDbTag(seriesDbTag);
    }
    idData(1) = seriesDbTag;
  }
  else
    idData(0) = -1;

  

  data(0) = this->getTag();
  data(1) = nFilter;
  data(2) = ta;
  data(3) = td;
  dataomegac = *omegac;
  databeta = *beta;

  int res = theChannel.sendID(dbTag, cTag, idData);
  res += theChannel.sendVector(dbTag, cTag, data);
  res += theChannel.sendVector(dbTag, cTag, dataomegac);
  res += theChannel.sendVector(dbTag, cTag, databeta);
  if (res < 0)
  {
    opserr << " URDDampingbeta::sendSelf() - data could not be sent\n" ;
    return -1;
  }

  if (fac)
  {
    res = fac->sendSelf(cTag, theChannel);
    if (res < 0)
    {
      opserr << " URDDampingbeta::sendSelf() - failed to send factor series\n";
      return res;
    }
  }

  return 0;
}


int 
URDDampingbeta::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  static ID idData(2);
  static Vector data(4);
  static Vector dataomegac(nFilter);
  static Vector databeta(nFilter);

  int res = theChannel.recvID(dbTag, cTag, idData);
  res += theChannel.recvVector(dbTag, cTag, data);
  res += theChannel.recvVector(dbTag, cTag, dataomegac);
  res += theChannel.recvVector(dbTag, cTag, databeta);
  if (res < 0) {opserr << " URDDampingbeta::recvSelf() - data could not be received\n" ;
    return -1;
  }

  int seriesClassTag = idData(0);
  if (seriesClassTag != -1)
  {
    int seriesDbTag = idData(1);
    if (fac == 0 || fac->getClassTag() != seriesClassTag)
    {
      if (fac != 0)
        delete fac;
      fac = theBroker.getNewTimeSeries(seriesClassTag);
      if (fac == 0)
      {
        opserr << "GroundMotion::recvSelf - could not create a Series object\n";
        return -2;
      }
    }
    fac->setDbTag(seriesDbTag);
    res = fac->recvSelf(cTag, theChannel, theBroker);
    if (res < 0)
    {
      opserr << "URDDampingbeta::recvSelf() - factor series could not be received\n";
      return res;
    }
  }
    
  this->setTag((int)data(0));
  nFilter = data(1);
  ta = data(2);
  td = data(3);
  *omegac = dataomegac;
  *beta = databeta;

  Initialize();
  return 0;
}

void
URDDampingbeta::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE)
  {
    s << "\nDamping: " << this->getTag() << " Type: URDDampingbeta";
    s << "\tnumber of filters: " << nFilter << endln;
    s << "\tfrequency: " << (*omegac) / 6.28318530718 << endln;
    s << "\tbeta: " << (*beta) << endln;
    s << "\tactivation time: " << ta << endln;
    s << "\tdeactivation time: " << td << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON)
  {
    s << "\t\t\t{\"name\": \"" << this->getTag() << "\", \"type\": \"URDDampingbeta\"";
    s << "\tnumber of filters: " << nFilter << endln;
    s << ", \"frequency\": [" << (*omegac) / 6.28318530718 << "]";
    s << ", \"beta\": [" << (*beta) << "]";
    s << ", \"activation time\": [" << ta << "]";
    s << ", \"deactivation time\": [" << td << "]";
    s << "}";
  }
}
