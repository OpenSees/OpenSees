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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the implementation of 
// MultilinearBackbone, which is a backbone defined by
// many points

#include <MultilinearBackbone.h>
#include <Vector.h>
#include <Channel.h>

#include <elementAPI.h>

void *
OPS_MultilinearBackbone(void)
{
  HystereticBackbone *theBackbone = 0;

  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "Invalid number of args, want: hystereticBackbone Multilinear tag? e1? s1? e2? s2? ..." << endln;
    return 0;
  }

  int iData[1];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for hystereticBackbone Multilinear" << endln;
    return 0;
  }

  int numPoints = OPS_GetNumRemainingInputArgs()/2;
  numData = 2*numPoints;
  Vector e(numPoints);
  Vector s(numPoints);
  double *dData = new double[numData];
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for hystereticBackbone Multilinear" << endln;
    return 0;
  }
  for (int i = 0; i < numPoints; i++) {
    e(i) = dData[2*i];
    s(i) = dData[2*i+1];
  }
  
  theBackbone = new MultilinearBackbone(iData[0], numPoints, e, s);
  if (theBackbone == 0) {
    opserr << "WARNING could not create MultilinearBackbone\n";
    return 0;
  }

  delete [] dData;
  
  return theBackbone;
}

MultilinearBackbone::MultilinearBackbone(int tag, int num,
					 const Vector &def, const Vector &force):
  HystereticBackbone(tag,BACKBONE_TAG_Multilinear),
  E(0), e(0), s(0), c(0), numPoints(num)
{
  E = new double [numPoints];
  if (E == 0)
    opserr << "MultilinearBackbone::MultilinearBackbone -- failed to allocate tangent array" << endln;
  
  e = new double [numPoints+1];
  if (e == 0)
    opserr << "MultilinearBackbone::MultilinearBackbone -- failed to allocate strain array" << endln;
  
  s = new double [numPoints+1];
  if (s == 0)
    opserr << "MultilinearBackbone::MultilinearBackbone -- failed to allocate stress array" << endln;
  
  c = new double [numPoints+1];
  if (c == 0)
    opserr << "MultilinearBackbone::MultilinearBackbone -- failed to allocate energy array" << endln;
  
  e[0] = s[0] = c[0] = 0.0;
  
  bool error = false;
  
  int i;
  
  for (i = 1; i <= numPoints; i++) {
    e[i] = def(i-1);
    s[i] = force(i-1);
  }
  
  for (i = 1; i <= numPoints; i++)
    if (e[i] < e[i-1])
      error = true;
  
  if (error) {
    delete [] E;
    delete [] e;
    delete [] s;
    delete [] c;
    
    opserr << "MultilinearBackbone::MultilinearBackbone -- input backbone is not unique (one-to-one)" << endln;
  }
  
  for (i = 1; i <= numPoints; i++) {
    E[i-1] = (s[i]-s[i-1])/(e[i]-e[i-1]);
    c[i] = c[i-1] + 0.5*(s[i]-s[i-1])*(e[i]-e[i-1]);
  }
}

MultilinearBackbone::MultilinearBackbone():
  HystereticBackbone(0,BACKBONE_TAG_Multilinear), 
  E(0), e(0), s(0), c(0), numPoints(0)
{

}

MultilinearBackbone::~MultilinearBackbone()
{
  if (E != 0)
    delete [] E;

  if (e != 0)
    delete [] e;

  if (s != 0)
    delete [] s;

  if (c != 0)
    delete [] c;
}

double
MultilinearBackbone::getTangent (double strain)
{
  for (int i = 1; i <= numPoints; i++) {
    if (strain < e[i])
      return E[i-1];
  }
  
  return E[0]*1.0e-9;
}

double
MultilinearBackbone::getStress (double strain)
{
  for (int i = 1; i <= numPoints; i++) {
    if (strain < e[i])
      return s[i-1] + E[i-1]*(strain-e[i-1]);
  }
  
  return s[numPoints];
}

double
MultilinearBackbone::getEnergy (double strain)
{
  for (int i = 1; i <= numPoints; i++) {
    if (strain < e[i])
      return c[i-1] + 0.5*E[i-1]*(strain-e[i-1])*(strain-e[i-1]);
  }
  
  return c[numPoints] + s[numPoints]*(strain-e[numPoints]);
}

double
MultilinearBackbone::getYieldStrain(void)
{
  return e[1];
}

HystereticBackbone*
MultilinearBackbone::getCopy(void)
{
  Vector def(&e[1], numPoints);
  Vector force(&s[1], numPoints);
  
  MultilinearBackbone *theCopy = 
    new MultilinearBackbone (this->getTag(), numPoints, def, force);
  
  return theCopy;
}

void
MultilinearBackbone::Print (OPS_Stream &o, int flag)
{
  Vector def(&e[1], numPoints);
  Vector force(&s[1], numPoints);
  
  o << "MultilinearBackbone, tag: " << this->getTag() << endln;
  o << "\tStrains: " << def << endln;
  o << "\tStresses: " << force << endln;
}

int
MultilinearBackbone::setVariable (char *argv)
{
  if (strcmp(argv,"yieldStrain") == 0)
    return 1;
  else
    return -1;
}

int
MultilinearBackbone::getVariable (int varID, double &theValue)
{
  switch (varID) {
  case 1:
    theValue = e[1];
    return 1;
  default:
    return -1;
  }
}

int
MultilinearBackbone::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
MultilinearBackbone::recvSelf(int commitTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}
