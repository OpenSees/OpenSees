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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-04-02 22:02:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/FeapMaterial.cpp,v $
                                                                        
// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// FeapMaterial. FeapMaterial wraps a Feap material subroutine.

#include <OPS_Globals.h>
#include <FeapMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

double FeapMaterial::d[200];
double FeapMaterial::sig[6];
double FeapMaterial::dd[36];

Vector FeapMaterial::strain3(3);
Vector FeapMaterial::strain4(4);
Vector FeapMaterial::strain6(6);

Vector FeapMaterial::sigma3(3);
Vector FeapMaterial::sigma4(4);
Vector FeapMaterial::sigma6(sig,6);

Matrix FeapMaterial::tangent3(3,3);
Matrix FeapMaterial::tangent4(4,4);
Matrix FeapMaterial::tangent6(dd,6,6);

FeapMaterial::FeapMaterial(int tag, int classTag, int nhv, int ndata, double r)
  :NDMaterial(tag,classTag), ud(0), hstv(0), rho(r),
   numHV(nhv), numData(ndata), myFormulation(ThreeDimensional)
{
  if (numHV < 0)
    numHV = 0;
  
  if (numHV > 0) {
    // Allocate history array
    hstv = new double[2*numHV];
    if (hstv == 0) {
      opserr <<  "FeapMaterial::FeapMaterial -- failed to allocate history array -- type: " << classTag << endln;
      exit(-1);
    }
			    
    
    // Zero out the history variables
    for (int i = 0; i < 2*numHV; i++)
      hstv[i] = 0.0;
  }
  
  if (numData < 0)
    numData = 0;
  
  if (numData > 0) {
    // Allocate material parameter array
    ud = new double[numData];
    if (ud == 0) {
      opserr << "FeapMaterial::FeapMaterial -- failed to allocate ud array -- type: " << classTag << endln;
      exit(-1);
    }
  }
  
  // Zero the strain vector
  for (int i = 0; i < 6; i++)
    eps[i] = 0.0;
}

FeapMaterial::FeapMaterial(int classTag)
  :NDMaterial(0,classTag), ud(0), hstv(0), rho(0.0),
   numHV(0), numData(0), myFormulation(Unknown)
{
  // Zero the strain vector
  for (int i = 0; i < 6; i++)
    eps[i] = 0.0;
}

FeapMaterial::~FeapMaterial()
{
  if (ud != 0)
    delete [] ud;
  
  if (hstv != 0)
    delete [] hstv;
}

// Here we are assuming the following order on strains
// \epsilon = {11, 22, 33, 12, 23, 31}
int
FeapMaterial::setTrialStrain(const Vector &strain)
{
  switch (myFormulation) {
  case ThreeDimensional:
    eps[0] = strain(0);
    eps[1] = strain(1);
    eps[2] = strain(2);
    eps[3] = strain(3);
    eps[4] = strain(4);
    eps[5] = strain(5);
    break;
  case PlaneStrain:
    eps[0] = strain(0);
    eps[1] = strain(1);
    eps[3] = strain(2);
    break;
  case AxiSymmetric:
    eps[0] = strain(0);
    eps[1] = strain(1);
    eps[2] = strain(2);
    eps[3] = strain(3);
    break;
  default:
    opserr << "FeapMaterial::FeapMaterial -- unknown material formulation\n";
    exit(-1);
    break;
  }

  return 0;
}

// Here we are assuming the following order on strains
// \epsilon = {11, 22, 33, 12, 23, 31}
const Vector&
FeapMaterial::getStrain(void)
{
  switch (myFormulation) {
  case ThreeDimensional:
    strain6(0) = eps[0];
    strain6(1) = eps[1];
    strain6(2) = eps[2];
    strain6(3) = eps[3];
    strain6(4) = eps[4];
    strain6(5) = eps[5];
    return strain6;
  case PlaneStrain:
    strain3(0) = eps[0];
    strain3(1) = eps[1];
    strain3(2) = eps[3];
    return strain3;
  case AxiSymmetric:
    strain4(0) = eps[0];
    strain4(1) = eps[1];
    strain4(2) = eps[2];
    strain4(3) = eps[3];
    return strain4;
  default:
    opserr << "FeapMaterial::getSTrain -- unknown material formulation\n";
    exit(-1);
      
    return strain6;
  }
}

// Here we are assuming the following order on stresses
// \sigma = {11, 22, 33, 12, 23, 31}
const Vector&
FeapMaterial::getStress(void)
{
  int isw = 3;

  // Invoke Feap subroutine
  this->invokeSubroutine(isw);
  
  switch (myFormulation) {
  case ThreeDimensional:
    return sigma6;
  case PlaneStrain:
    sigma3(0) = sig[0];
    sigma3(1) = sig[1];
    sigma3(2) = sig[3];
    return sigma3;
  case AxiSymmetric:
    sigma4(0) = sig[0];
    sigma4(1) = sig[1];
    sigma4(2) = sig[2];
    sigma4(3) = sig[3];
    return sigma4;
  default:
    opserr << "FeapMaterial::getStress -- unknown material formulation\n";
      exit(-1);
    return sigma6;
  }
}

const Matrix&
FeapMaterial::getTangent(void)
{
  int isw = 6;

  // Invoke Feap subroutine
  this->invokeSubroutine(isw);

  int i,j;

  switch (myFormulation) {
  case ThreeDimensional:
    return tangent6;
  case PlaneStrain:
    tangent3(0,0) = tangent6(0,0);
    tangent3(0,1) = tangent6(0,1);
    tangent3(0,2) = tangent6(0,3);
    tangent3(1,0) = tangent6(1,0);
    tangent3(1,1) = tangent6(1,1);
    tangent3(1,2) = tangent6(1,3);
    tangent3(2,0) = tangent6(3,0);
    tangent3(2,1) = tangent6(3,1);
    tangent3(2,2) = tangent6(3,3);
    return tangent3;
  case AxiSymmetric:
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
	tangent4(i,j) = tangent6(i,j);
    return tangent4;
  default:
    opserr << "FeapMaterial::getTangent -- unknown material formulation\n";
    exit(-1);
    return tangent6;
  }
}

double
FeapMaterial::getRho(void)
{
  return rho;
}

int
FeapMaterial::commitState(void)
{
  // Set committed values equal to corresponding trial values
  for (int i = 0; i < numHV; i++)
    hstv[i] = hstv[i+numHV];
  
  return 0;
}

int
FeapMaterial::revertToLastCommit(void)
{
  // Set trial values equal to corresponding committed values
  for (int i = 0; i < numHV; i++)
    hstv[i+numHV] = hstv[i];
  
  return 0;
}

int
FeapMaterial::revertToStart(void)
{
  // Set all trial and committed values to zero
  for (int i = 0; i < 2*numHV; i++)
    hstv[i] = 0.0;
  
  for (int i = 0; i < 6; i++)
    eps[i] = 0.0;

  return 0;
}

NDMaterial*
FeapMaterial::getCopy(void)
{
  FeapMaterial *theCopy = 
    new FeapMaterial(this->getTag(), this->getClassTag(),
		     numHV, numData, rho);
  
  int i;
  for (i = 0; i < 2*numHV; i++)
    theCopy->hstv[i] = hstv[i];
  
  for (i = 0; i < numData; i++)
    theCopy->ud[i] = ud[i];
  
  theCopy->myFormulation = myFormulation;
  
  return theCopy;
}

NDMaterial*
FeapMaterial::getCopy(const char *type)
{
  FeapMaterial *theCopy = (FeapMaterial*)this->getCopy();
  
  if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0)
    theCopy->myFormulation = ThreeDimensional;
  
  else if (strcmp(type, "PlaneStrain") == 0 || strcmp(type, "PlaneStrain2D") == 0)
    theCopy->myFormulation = PlaneStrain;
  
  else if (strcmp(type, "AxiSymmetric") == 0 || strcmp(type, "AxiSymmetric2D") == 0)
    theCopy->myFormulation = AxiSymmetric;
  
  else {
    opserr << "FeapMaterial::getCopy -- Invalid type (" << type << ") for FeapMaterial\n";
    return 0;
  }
  
  return theCopy;
}

const char*
FeapMaterial::getType(void) const
{
  switch (myFormulation) {
  case ThreeDimensional:
    return "ThreeDimensional";
  case PlaneStrain:
    return "PlaneStrain";
  case AxiSymmetric:
    return "AxiSymmetric";
  default:
    opserr << "FeapMaterial::getTYpe -- unknown material formulation\n";
    return "Unknown";
  }
}

int
FeapMaterial::getOrder(void) const
{
  switch (myFormulation) {
  case ThreeDimensional:
    return 6;
  case PlaneStrain:
    return 3;
  case AxiSymmetric:
    return 4;
  default:
    opserr << "FeapMaterial::getOrder -- unknown material formulation\n";
    return 0;
  }
}

int 
FeapMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static ID idData(4);
  
  idData(0) = this->getTag();
  idData(1) = numHV;
  idData(2) = numData;
  idData(3) = myFormulation;
  
  res += theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) 
    opserr << "FeapMaterial::sendSelf() - failed to send ID data\n";
  
  Vector vecData(numHV+numData);
  
  int i, j;
  // Copy history variables into vector
  for (i = 0; i < numHV; i++)
    vecData(i) = hstv[i];
  
  // Copy material properties into vector
  for (i = 0, j = numHV; i < numData; i++, j++)
    vecData(j) = ud[i];
  
  res += theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) 
    opserr << "FeapMaterial::sendSelf() - failed to send Vector data\n";
  
  return res;
}

int
FeapMaterial::recvSelf(int commitTag, Channel &theChannel,
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static ID idData(4);
  
  res += theChannel.recvID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "FeapMaterial::recvSelf() - failed to receive ID data\n";
    return res;
  }
  
  this->setTag(idData(0));
  numHV   = idData(1);
  numData = idData(2);
  myFormulation  = idData(3);
  
  Vector vecData(numHV+numData);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "FeapMaterial::recvSelf() - failed to receive Vector data\n";
    return res;
  }
  
  int i, j;
  // Copy history variables from vector
  for (i = 0; i < numHV; i++)
    hstv[i] = vecData(i);
  
  // Copy material properties from vector
  for (i = 0, j = numHV; i < numData; i++, j++)
    ud[i] = vecData(j);
  
  return res;
}
    
void
FeapMaterial::Print(OPS_Stream &s, int flag)
{
  s << "FeapMaterial, tag: " << this->getTag() << endln;
  s << "Material formulation: " << this->getType() << endln;
  s << "Material subroutine: ";
  
  switch (this->getClassTag()) {
  case ND_TAG_FeapMaterial01:
    s << "matl01" << endln;
    break;

  case ND_TAG_FeapMaterial02:
    s << "matl02" << endln;
    break;

  case ND_TAG_FeapMaterial03:
    s << "matl03" << endln;
    break;
    
    // Add more cases as needed
  default:
    s << this->getClassTag() << endln;
    break;
  }
  
  s << "Material density: " << rho << endln;
}

#ifdef _WIN32

extern "C" int FEAPCOMMON(double *dt, int *niter);

extern "C" int MATL01(double *eps, double *trace, double *td, double *d,
			       double *ud, double *hn, double *h1, int *nh,
			       double *sig, double *dd, int *isw);

extern "C" int MATL02(double *eps, double *trace, double *td, double *d,
			       double *ud, double *hn, double *h1, int *nh,
			       double *sig, double *dd, int *isw);

extern "C" int MATL03(double *eps, double *trace, double *td, double *d,
			       double *ud, double *hn, double *h1, int *nh,
			       double *sig, double *dd, int *isw);

#define feapcommon_ FEAPCOMMON
#define matl01_ MATL01
#define matl02_ MATL02
#define matl03_ MATL03

#else

extern "C" int feapcommon_(double *dt, int *niter);

extern "C" int matl01_(double *eps, double *trace, double *td, double *d,
		       double *ud, double *hn, double *h1, int *nh,
		       double *sig, double *dd, int *isw);

extern "C" int matl02_(double *eps, double *trace, double *td, double *d,
		       double *ud, double *hn, double *h1, int *nh,
		       double *sig, double *dd, int *isw);

extern "C" int matl03_(double *eps, double *trace, double *td, double *d,
		       double *ud, double *hn, double *h1, int *nh,
		       double *sig, double *dd, int *isw);

#endif

int
FeapMaterial::invokeSubroutine(int isw)
{
  // Trace of strain vector
  double trace = eps[0] + eps[1] + eps[2];

  // Temperature change (currently not used)
  double td = 0.0;

  // Zero out stress and tangent arrays as required to do so
  int i;
  for (i = 0; i < 6; i++) {
    sig[i] = 0.0;
    dd[i] = 0.0;
  }
  for ( ; i < 36; i++)
    dd[i] = 0.0;

  // Populate the FEAP material array for this model
  this->fillDArray();

  // Fill in the common blocks
  double dt = ops_Dt;  // From G3Globals.h
  int niter = 1;    // Need to count the number of global iterations!
  feapcommon_(&dt, &niter);

  switch (this->getClassTag()) {
  case ND_TAG_FeapMaterial01:
    matl01_(eps, &trace, &td, d, ud, hstv, &hstv[numHV], &numHV, sig, dd, &isw);
    break;

  case ND_TAG_FeapMaterial02:
    matl02_(eps, &trace, &td, d, ud, hstv, &hstv[numHV], &numHV, sig, dd, &isw);
    break;

  case ND_TAG_FeapMaterial03:
    matl03_(eps, &trace, &td, d, ud, hstv, &hstv[numHV], &numHV, sig, dd, &isw);
    break;
    
    // Add more cases as needed
  default:
    opserr << "FeapMaterial::invokeSubroutine -- unknown material type\n";
    return -1;
  }
  
  return 0;
}

int
FeapMaterial::fillDArray(void)
{
  // Must be done by subclasses!
  return 0;
}
