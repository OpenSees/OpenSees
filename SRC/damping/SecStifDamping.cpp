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
// $Source: /usr/local/cvs/OpenSees/SRC/damping/SecStifDamping.cpp,v $

// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 02/2020
// Revision: A
// 
// Description: This file contains the implementation for the SecStifDamping class.
// FRDamping provides the abstraction of an elemental damping imposition
// providing quasi tangent stiffness proportional damping
//
// What: "@(#) SecStifDamping.cpp, revA"

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <Domain.h>
#include <Channel.h>
#include <elementAPI.h>
#include <SecStifDamping.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>

//extern StaticAnalysis *theStaticAnalysis;

// constructor:
SecStifDamping::SecStifDamping(int tag, double b, double t1, double t2, TimeSeries *f):
Damping(tag, DMP_TAG_SecStifDamping),
beta(b), ta(t1), td(t2), fac(f), qd(0), q0(0), q0C(0)
{
  if (beta <= 0.0) opserr << "SecStifDamping::SecStifDamping:  Invalid damping factor\n";
}


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
SecStifDamping::SecStifDamping():
Damping(0, DMP_TAG_SecStifDamping),
beta(0.0), ta(0.0), td(0.0), fac(0), qd(0), q0(0), q0C(0)
{
    
}


// destructor:
SecStifDamping::~SecStifDamping() 
{
  if (fac) delete fac;
  if(qd) delete qd;
  if(q0) delete q0;
  if(q0C) delete q0C;
}


int
SecStifDamping::commitState(void)
{
  *q0C = *q0;
  return 0;
}


int
SecStifDamping::revertToLastCommit(void)
{
  *q0 = *q0C;
  return 0;
}


int
SecStifDamping::revertToStart(void)
{
  (*q0).Zero();
  (*q0C).Zero();
  return 0;
}


int 
SecStifDamping::setDomain(Domain *domain, int nComp)
{       
  theDomain = domain;
  
  qd = new Vector(nComp);
  q0 = new Vector(nComp);
  q0C = new Vector(nComp);
    
  return 0;
}


int
SecStifDamping::update(Vector q)
{       
  double t = theDomain->getCurrentTime();
  double dT = theDomain->getDT();
  StaticAnalysis **theStaticAnalysis = OPS_GetStaticAnalysis();
  if (*theStaticAnalysis)
  {
    (*qd).Zero();
  }
  else if (dT > 0.0)
  {
    *q0 = q;
    if (t > ta && t < td)
    {
      *qd = (beta / dT) * (*q0 - *q0C);
      if (fac) *qd *= fac->getFactor(t);
    }
    else
    {
      (*qd).Zero();
    }
  }
  return 0;
}

const Vector &
SecStifDamping::getDampingForce(void)
{
  return *qd;
}

double SecStifDamping::getStiffnessMultiplier(void)
{
  double t = theDomain->getCurrentTime();
  double dT = theDomain->getDT();
  double km = 0.0;
  StaticAnalysis **theStaticAnalysis = OPS_GetStaticAnalysis();
  if (!*theStaticAnalysis && dT > 0.0 && t > ta && t < td) km = beta / dT;
  return 1.0 + km;
}


Damping *SecStifDamping::getCopy(void)
{
    // create a new instance of SecStifDamping 
    
    SecStifDamping *theCopy;
    
    theCopy = new SecStifDamping(this->getTag(), beta, ta, td, fac);
    
    return theCopy;
}


int 
SecStifDamping::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  static ID idData(2);
  static Vector data(4);

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
  data(1) = beta;
  data(2) = ta;
  data(3) = td;

  int res = theChannel.sendID(dbTag, cTag, idData);
  res += theChannel.sendVector(dbTag, cTag, data);
  if (res < 0)
  {
    opserr << " UniformDamping::sendSelf() - data could not be sent\n" ;
    return -1;
  }

  if (fac)
  {
    res = fac->sendSelf(cTag, theChannel);
    if (res < 0)
    {
      opserr << " UniformDamping::sendSelf() - failed to send factor series\n";
      return res;
    }
  }

  return 0;
}


int 
SecStifDamping::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  static ID idData(2);
  static Vector data(4);

  int res = theChannel.recvID(dbTag, cTag, idData);
  res += theChannel.recvVector(dbTag, cTag, data);
  if (res < 0) {opserr << " UniformDamping::recvSelf() - data could not be received\n" ;
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
      opserr << "UniformDamping::recvSelf() - factor series could not be received\n";
      return res;
    }
  }
    
  this->setTag((int)data(0));
  beta = data(1);
  ta = data(2);
  td = data(3);

  if (beta <= 0.0) opserr << "SecStifDamping::recvSelf:  Invalid damping factor\n";

  return 0;
}

void
SecStifDamping::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_CURRENTSTATE)
	{
    s << "\nDamping: " << this->getTag() << " Type: SecStifDamping";
    s << "\tdamping factor: " << beta << endln;
    s << "\tactivation time: " << ta << endln;
    s << "\tdeactivation time: " << td << endln;
	}
	
	if (flag == OPS_PRINT_PRINTMODEL_JSON)
  {
    s << "\t\t\t{\"name\": \"" << this->getTag() << "\", \"type\": \"SecStifDamping\"";
    s << ", \"damping factor\": [" << beta << "]";
    s << ", \"activation time\": [" << ta << "]";
    s << ", \"deactivation time\": [" << td << "]";
    s << "}";
  }
}
