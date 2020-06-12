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

void* OPS_SecStifDamping()
{
  if(OPS_GetNumRemainingInputArgs() < 2)
  {
    opserr<<"insufficient arguments for SecStifDamping\n";
    return 0;
  }

  // get tag
  int tag;
  int numData = 1;
  if(OPS_GetIntInput(&numData,&tag) < 0) return 0;
  
  double beta;
  if(OPS_GetDoubleInput(&numData,&beta) < 0) return 0;
      
  return new SecStifDamping(tag,beta);
}

// constructor:
SecStifDamping::SecStifDamping(int tag, double b):
Damping(tag, DMP_TAG_SecStifDamping),
beta(b), qd(0), q0(0), q0C(0)
{
  if (beta <= 0.0) opserr << "SecStifDamping::SecStifDamping:  Invalid damping factor\n";
}


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
SecStifDamping::SecStifDamping():
Damping(0, DMP_TAG_SecStifDamping),
beta(0.0), qd(0), q0(0), q0C(0)
{
    
}


// destructor:
SecStifDamping::~SecStifDamping() 
{
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
  double dT = theDomain->getDT();
  if (dT > 0.0)
  {
    *q0 = q;
    *qd = (beta / dT) * (*q0 - *q0C);
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
  double dT = theDomain->getDT();
  double km = 1.0;
  if (dT > 0.0) km += beta / dT;
  return km;
}


Damping *SecStifDamping::getCopy(void)
{
    // create a new instance of SecStifDamping 
    
    SecStifDamping *theCopy;
    
    theCopy = new SecStifDamping(this->getTag(), beta);
    
    return theCopy;
}


int 
SecStifDamping::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
    
  return res;
}


int 
SecStifDamping::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  return res;
}

void
SecStifDamping::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_CURRENTSTATE)
	{
    s << "\nDamping: " << this->getTag() << " Type: SecStifDamping";
    s << "\tdamping factor: " << beta << endln;
	}
	
	if (flag == OPS_PRINT_PRINTMODEL_JSON)
  {
    s << "\t\t\t{\"name\": \"" << this->getTag() << "\", \"type\": \"SecStifDamping\"";
    s << ", \"damping factor\": [" << beta << "]";
    s << "}";
  }
}
