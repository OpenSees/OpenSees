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

#include <Elliptical2.h>           
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <elementAPI.h>

#include <MaterialResponse.h>
#include <Information.h>
#include <Parameter.h>

Vector Elliptical2::s(2);
Matrix Elliptical2::ks(2,2);
ID Elliptical2::code(2);

void* OPS_Elliptical2()
{
    if (OPS_GetNumRemainingInputArgs() < 8) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section Elliptical tag? E1? E2? sigY1? sigY2? Hiso? Hkin1? Hkin2? <code1? code2?>" << endln;
	return 0;
    }    

    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid Elliptical tag" << endln;
	return 0;
    }

    numdata = 7;
    double data[7];
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid double inputs\n";
	opserr << "section Elliptical: " << tag << endln;
	return 0;
    }
    double E1 = data[0];
    double E2 = data[1];    
    double sigY1 = data[2];
    double sigY2 = data[3];    
    double Hi = data[4];
    double Hk1 = data[5];
    double Hk2 = data[6];    

    if (OPS_GetNumRemainingInputArgs() > 1) {
      int code1, code2;
      const char* type1 = OPS_GetString();
      const char* type2 = OPS_GetString();
      if (strcmp(type1,"Mz") == 0)
	code1 = SECTION_RESPONSE_MZ;
      else if (strcmp(type1,"P") == 0)
	code1 = SECTION_RESPONSE_P;
      else if (strcmp(type1,"Vy") == 0)
	code1 = SECTION_RESPONSE_VY;
      else if (strcmp(type1,"My") == 0)
	code1 = SECTION_RESPONSE_MY;
      else if (strcmp(type1,"Vz") == 0)
	code1 = SECTION_RESPONSE_VZ;
      else if (strcmp(type1,"T") == 0)
	code1 = SECTION_RESPONSE_T;
      else {
	opserr << "WARNING invalid code 1 " << type1 << endln;
	opserr << "section Elliptical: " << tag << endln;
	return 0;
      }
      
      if (strcmp(type2,"Mz") == 0)
	code2 = SECTION_RESPONSE_MZ;
      else if (strcmp(type2,"P") == 0)
	code2 = SECTION_RESPONSE_P;
      else if (strcmp(type2,"Vy") == 0)
	code2 = SECTION_RESPONSE_VY;
      else if (strcmp(type2,"My") == 0)
	code2 = SECTION_RESPONSE_MY;
      else if (strcmp(type2,"Vz") == 0)
	code2 = SECTION_RESPONSE_VZ;
      else if (strcmp(type2,"T") == 0)
	code2 = SECTION_RESPONSE_T;
      else {
	opserr << "WARNING invalid code 2 " << type2 << endln;
	opserr << "section Elliptical: " << tag << endln;
	return 0;
      }
      return new Elliptical2(tag, E1, E2, sigY1, sigY2, Hi, Hk1, Hk2, code1, code2);
    } else {
      return new Elliptical2(tag, E1, E2, sigY1, sigY2, Hi, Hk1, Hk2);
    }

}

Elliptical2::Elliptical2
(int tag, double e1, double e2, double sy1, double sy2, double Hi, 
 double Hk1, double Hk2, int c1, int c2) :
  SectionForceDeformation(tag, SEC_TAG_Elliptical2),
  Hiso(Hi), code1(c1), code2(c2), parameterID(0), SHVs(0)
{
  E[0] = e1;
  E[1] = e2;

  Hkin[0] = Hk1;
  Hkin[1] = Hk2;

  sigY[0] = sy1;
  sigY[1] = sy2;

  for (int i = 0; i < 2; i++) {
    eP_n[i]  = 0.0;
    eP_n1[i] = 0.0;
    e_n1[i] = 0.0;
  }

  dg_n1 = 0.0;
  
  alpha_n  = 0.0;
  alpha_n1 = 0.0;
}

Elliptical2::Elliptical2() :
  SectionForceDeformation(0, SEC_TAG_Elliptical2),
  Hiso(0.0), parameterID(0), SHVs(0)
{
  for (int i = 0; i < 2; i++) {
    E[i] = 0.0;
    Hkin[i] = 0.0;
    sigY[i] = 0.0;

    eP_n[i]  = 0.0;
    eP_n1[i] = 0.0;
    e_n1[i] = 0.0;
  }
  
  dg_n1 = 0.0;

  alpha_n  = 0.0;
  alpha_n1 = 0.0;
  
  code1 = SECTION_RESPONSE_MZ;
  code2 = SECTION_RESPONSE_VY;
}

Elliptical2::~Elliptical2()
{
  if (SHVs != 0)
    delete SHVs;
}

int
Elliptical2::setTrialSectionDeformation(const Vector &e)
{
  e_n1[0] = e(0);
  e_n1[1] = e(1);
  
  return 0;
}

const Matrix&
Elliptical2::getSectionTangent(void)
{
  static Vector str(2);

  // Compute trial stress using elastic tangent
  str(0) = E[0]*(e_n1[0]-eP_n[0]);
  str(1) = E[1]*(e_n1[1]-eP_n[1]);

  double xsi[2];
  xsi[0] = str(0) - Hkin[0]*eP_n[0];
  xsi[1] = str(1) - Hkin[1]*eP_n[1];

  double Q[2];
  Q[0] = 1.0/(sigY[0]*sigY[0]);
  Q[1] = 1.0/(sigY[1]*sigY[1]);

  double q = sqrt(xsi[0]*Q[0]*xsi[0] + xsi[1]*Q[1]*xsi[1]);
  double F = q - (1.0 + Hiso*alpha_n);

  //if (F <= -100*DBL_EPSILON) {
  if (F < -10*DBL_EPSILON) {
    ks(0,0) = E[0];
    ks(1,1) = E[1];
    ks(0,1) = ks(1,0) = 0.0;

    return ks;
  }

  double dg = 0.0;

  static Vector x(3);
  x(0) = xsi[0];
  x(1) = xsi[1];
  x(2) = 0.0;

  static Vector dx(3);

  double n[2];
  n[0] = Q[0]*xsi[0]/q;
  n[1] = Q[1]*xsi[1]/q;

  static Vector R(3);
  R(0) = 0.0;
  R(1) = 0.0;
  R(2) = F;

  static Matrix J(3,3);

  int numIter = 0; int maxIter = 25;
  while (R.Norm() > 1.0e-14 && numIter < maxIter) {

    J(0,0) = 1.0 + dg/q*(E[0]+Hkin[0])*(Q[0]-n[0]*n[0]);
    J(0,1) =       dg/q*(E[0]+Hkin[0])*(    -n[0]*n[1]);
    J(0,2) = (E[0]+Hkin[0])*n[0];

    J(1,0) =       dg/q*(E[1]+Hkin[1])*(    -n[1]*n[0]);
    J(1,1) = 1.0 + dg/q*(E[1]+Hkin[1])*(Q[1]-n[1]*n[1]);
    J(1,2) = (E[1]+Hkin[1])*n[1];

    J(2,0) = n[0];
    J(2,1) = n[1];
    J(2,2) = -Hiso;

    J.Solve(R, dx);

    x = x - dx;
    dg = x(2);
    dg_n1 = dg;

    q = sqrt(x(0)*Q[0]*x(0) + x(1)*Q[1]*x(1));

    n[0] = Q[0]*x(0)/q;
    n[1] = Q[1]*x(1)/q;

    R(0) = x(0) - xsi[0] + dg*(E[0]+Hkin[0])*n[0];
    R(1) = x(1) - xsi[1] + dg*(E[1]+Hkin[1])*n[1];
    R(2) = q - (1.0 + Hiso*(alpha_n+dg));

    numIter++;
  }

  alpha_n1 = alpha_n + dg;

  eP_n1[0] = eP_n[0] + dg*n[0];
  eP_n1[1] = eP_n[1] + dg*n[1];
  

  static Matrix A(2,2);
  static Matrix invA(2,2);
  static Matrix B(2,2);

  A(0,0) = 1.0 + dg/q*Hkin[0]*(Q[0]-n[0]*n[0]);
  A(0,1) =       dg/q*Hkin[0]*(    -n[0]*n[1]);
  A(1,0) =       dg/q*Hkin[1]*(    -n[1]*n[0]);
  A(1,1) = 1.0 + dg/q*Hkin[1]*(Q[1]-n[1]*n[1]);

  A.Invert(invA);

  B(0,0) = dg/q*E[0]*(Q[0]-n[0]*n[0]);
  B(0,1) = dg/q*E[0]*(    -n[0]*n[1]);
  B(1,0) = dg/q*E[1]*(    -n[1]*n[0]);
  B(1,1) = dg/q*E[1]*(Q[1]-n[1]*n[1]);
  
  J(0,0) = 1.0 + B(0,0)*invA(0,0) + B(0,1)*invA(1,0);
  J(0,1) =       B(0,0)*invA(0,1) + B(0,1)*invA(1,1);
  J(1,0) =       B(1,0)*invA(0,0) + B(1,1)*invA(1,0);
  J(1,1) = 1.0 + B(1,0)*invA(0,1) + B(1,1)*invA(1,1);

  J(2,0) = n[0]*invA(0,0) + n[1]*invA(1,0);
  J(2,1) = n[0]*invA(0,1) + n[1]*invA(1,1);

  double b[2];
  b[0] = invA(0,0)*Hkin[0]*n[0] + invA(0,1)*Hkin[1]*n[1];
  b[1] = invA(1,0)*Hkin[0]*n[0] + invA(1,1)*Hkin[1]*n[1];

  J(0,2) = E[0]*n[0] - (B(0,0)*b[0] + B(0,1)*b[1]);
  J(1,2) = E[1]*n[1] - (B(1,0)*b[0] + B(1,1)*b[1]);

  J(2,2) = -(n[0]*b[0]+n[1]*b[1]) - Hiso;


  static Matrix C(3,3);
  J.Invert(C);

  ks(0,0) = E[0]*C(0,0);
  ks(1,0) = E[0]*C(1,0);
  ks(0,1) = E[1]*C(0,1);
  ks(1,1) = E[1]*C(1,1);

  return ks;
}

const Matrix&
Elliptical2::getInitialTangent(void)
{
  ks(0,0) = E[0];
  ks(1,1) = E[1];
  ks(0,1) = ks(1,0) = 0.0;

  return ks;
}

const Vector&
Elliptical2::getStressResultant(void)
{
  // Compute trial stress using elastic tangent
  s(0) = E[0]*(e_n1[0]-eP_n[0]);
  s(1) = E[1]*(e_n1[1]-eP_n[1]);

  double xsi[2];
  xsi[0] = s(0) - Hkin[0]*eP_n[0];
  xsi[1] = s(1) - Hkin[1]*eP_n[1];

  double Q[2];
  Q[0] = 1.0/(sigY[0]*sigY[0]);
  Q[1] = 1.0/(sigY[1]*sigY[1]);

  double q = sqrt(xsi[0]*Q[0]*xsi[0] + xsi[1]*Q[1]*xsi[1]);
  double F = q - (1.0 + Hiso*alpha_n);

  //if (F <= -100*DBL_EPSILON) {
  if (F < -10*DBL_EPSILON) {
    eP_n1[0] = eP_n[0];
    eP_n1[1] = eP_n[1];

    alpha_n1 = alpha_n;

    return s;
  }

  double dg = 0.0;

  static Vector dx(3);

  double n[2];
  n[0] = Q[0]*xsi[0]/q;
  n[1] = Q[1]*xsi[1]/q;

  static Vector x(3);
  x(0) = xsi[0];
  x(1) = xsi[1];
  x(2) = 0.0;

  static Vector R(3);
  R(0) = 0.0;
  R(1) = 0.0;
  R(2) = F;

  static Matrix J(3,3);

  int numIter = 0; int maxIter = 25;
  while (R.Norm() > 1.0e-14 && numIter < maxIter) {

    J(0,0) = 1.0 + dg/q*(E[0]+Hkin[0])*(Q[0]-n[0]*n[0]);
    J(0,1) =       dg/q*(E[0]+Hkin[0])*(    -n[0]*n[1]);
    J(0,2) = (E[0]+Hkin[0])*n[0];

    J(1,0) =       dg/q*(E[1]+Hkin[1])*(    -n[1]*n[0]);
    J(1,1) = 1.0 + dg/q*(E[1]+Hkin[1])*(Q[1]-n[1]*n[1]);
    J(1,2) = (E[1]+Hkin[1])*n[1];

    J(2,0) = n[0];
    J(2,1) = n[1];
    J(2,2) = -Hiso;

    J.Solve(R, dx);

    x = x - dx;
    dg = x(2);
    dg_n1 = dg;

    q = sqrt(x(0)*Q[0]*x(0) + x(1)*Q[1]*x(1));

    n[0] = Q[0]*x(0)/q;
    n[1] = Q[1]*x(1)/q;

    R(0) = x(0) - xsi[0] + dg*(E[0]+Hkin[0])*n[0];
    R(1) = x(1) - xsi[1] + dg*(E[1]+Hkin[1])*n[1];
    R(2) = q - (1.0 + Hiso*(alpha_n+dg));

    numIter++;
  }

  if (numIter == maxIter) {
    //opserr << "Elliptical2::getStressResultant norm(R)=" << R.Norm() << endln;
  }

  alpha_n1 = alpha_n + dg;

  eP_n1[0] = eP_n[0] + dg*n[0];
  eP_n1[1] = eP_n[1] + dg*n[1];

  s(0) = x(0) + Hkin[0]*eP_n1[0];
  s(1) = x(1) + Hkin[1]*eP_n1[1];
  
  return s;
}

const Vector&
Elliptical2::getSectionDeformation(void)
{
  static Vector e(2);

  // Write to static variable for return
  e(0) = e_n1[0];
  e(1) = e_n1[1];
  
  return e;
}

int
Elliptical2::commitState(void)
{
  eP_n[0] = eP_n1[0];
  eP_n[1] = eP_n1[1];
  
  alpha_n = alpha_n1;
  
  return 0;
}

int
Elliptical2::revertToLastCommit(void)
{
  return 0;
}

int
Elliptical2::revertToStart(void)
{
  for (int i = 0; i < 2; i++) {
    eP_n[i]  = 0.0;
    eP_n1[i] = 0.0;
    e_n1[i] = 0.0;
  }

  dg_n1 = 0.0;
  
  alpha_n  = 0.0;
  alpha_n1 = 0.0;
  
  if (SHVs != 0) {
    delete SHVs;
    SHVs = 0;
  }

  return 0;
}

SectionForceDeformation*
Elliptical2::getCopy(void)
{
  Elliptical2 *theCopy =
    new Elliptical2 (this->getTag(), E[0], E[1], sigY[0], sigY[1],
		    Hiso, Hkin[0], Hkin[1], code1, code2);
  
  for (int i = 0; i < 2; i++) {
    theCopy->eP_n[i]  = eP_n[i];
    theCopy->eP_n1[i] = eP_n1[i];
    theCopy->e_n1[i] = e_n1[i];
  }
  
  theCopy->alpha_n  = alpha_n;
  theCopy->alpha_n1 = alpha_n1;
  theCopy->dg_n1 = dg_n1;
  
  theCopy->parameterID = parameterID;

  return theCopy;
}

const ID&
Elliptical2::getType(void)
{
  code(0) = code1;
  code(1) = code2;

  return code;
}

int
Elliptical2::getOrder(void) const
{
  return 2;
}

int 
Elliptical2::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(13);
  
  data(0) = this->getTag();
  data(1) = E[0];
  data(2) = E[1];
  data(3) = sigY[0];
  data(4) = sigY[1];
  data(5) = Hiso;
  data(6) = Hkin[0];
  data(7) = Hkin[1];
  data(8) = code1;
  data(9) = code2;
  data(10) = eP_n[0];
  data(11) = eP_n[1];
  data(12) = alpha_n;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Elliptical2::sendSelf() - failed to send data\n";

  return res;
}

int 
Elliptical2::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(13);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
    opserr << "Elliptical2::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    E[0] = data(1);
    E[1] = data(2);
    sigY[0] = data(3);
    sigY[1] = data(4);
    Hiso = data(5);
    Hkin[0] = data(6);
    Hkin[1] = data(7);
    code1 = (int)data(8);
    code2 = (int)data(9);
    eP_n[0] = data(10);
    eP_n[1] = data(11);
    alpha_n = data(12);
    
    // Set the trial state variables
    revertToLastCommit();
  }
  
  return res;
}

void
Elliptical2::Print(OPS_Stream &s, int flag)
{
  s << "Elliptical2, tag: " << this->getTag() << endln;
  s << "\tE1, E2:    " << E[0] << ", " << E[1] << endln;
  s << "\tsigY1, sigY2:    " << sigY[0] << ", " << sigY[1] << endln;
  s << "\tHiso: " << Hiso << endln;
  s << "\tHkin1, Hkin2: " << Hkin[0] << ", " << Hkin[1] << endln;
  s << "\talpha_n: " << alpha_n << endln;
  s << "\tep_n: " << eP_n[0] << ' ' << eP_n[1] << endln;
}

Response*
Elliptical2::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  if (strcmp(argv[0],"plasticDeformation") == 0)
    return new MaterialResponse(this, 123, Vector(2));
  else
    return SectionForceDeformation::setResponse(argv, argc, output);
}

int
Elliptical2::getResponse(int responseID, Information &eleInformation)
{
  if (responseID == 123) {
    Vector &theVec = *(eleInformation.theVector);
    theVec(0) = eP_n[0];
    theVec(1) = eP_n[1];
    return eleInformation.setVector(theVec);
  }
  else
    return -1;
}

int 
Elliptical2::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"Fy1") == 0) {
    param.setValue(sigY[0]);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"Fy2") == 0) {
    param.setValue(sigY[1]);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"Fy") == 0) {
    param.setValue(sigY[0]);
    return param.addObject(12, this);
  }
  if (strcmp(argv[0],"k1") == 0 || strcmp(argv[0],"E1") == 0) {
    param.setValue(E[0]);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"k2") == 0 || strcmp(argv[0],"E2") == 0) {
    param.setValue(E[1]);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"k") == 0 || strcmp(argv[0],"E") == 0) {
    param.setValue(E[0]);
    return param.addObject(34, this);
  }
  if (strcmp(argv[0],"Hkin1") == 0) {
    param.setValue(Hkin[0]);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"Hkin2") == 0) {
    param.setValue(Hkin[1]);
    return param.addObject(6, this);
  }
  if (strcmp(argv[0],"Hkin") == 0) {
    param.setValue(Hkin[0]);
    return param.addObject(56, this);
  }
  if (strcmp(argv[0],"Hiso") == 0) {
    param.setValue(Hiso);
    return param.addObject(7, this);
  }

  return -1;
}

int 
Elliptical2::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    sigY[0] = info.theDouble;
    return 0;
  case 2:
    sigY[1] = info.theDouble;
    return 0;
  case 12:
    sigY[0] = info.theDouble;
    sigY[1] = info.theDouble;
    return 0;
  case 3:
    E[0] = info.theDouble;
    return 0;
  case 4:
    E[1] = info.theDouble;
    return 0;
  case 34:
    E[0] = info.theDouble;
    E[1] = info.theDouble;
    return 0;
  case 5:
    Hkin[0] = info.theDouble;
    return 0;
  case 6:
    Hkin[1] = info.theDouble;
    return 0;
  case 56:
    Hkin[0] = info.theDouble;
    Hkin[1] = info.theDouble;
    return 0;
  case 7:
    Hiso = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
Elliptical2::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector& 
Elliptical2::getStressResultantSensitivity(int gradIndex,
                                          bool conditional)
{
  s.Zero();

  double dEdh[2];
  dEdh[0] = 0.0;
  dEdh[1] = 0.0;
  double dHkindh[2];
  dHkindh[0] = 0.0;
  dHkindh[1] = 0.0;
  double dHisodh = 0.0;
  double dFydh[2];
  dFydh[0] = 0.0;
  dFydh[1] = 0.0;

  if (parameterID == 1 || parameterID == 12)
    dFydh[0] = 1.0;
  if (parameterID == 2 || parameterID == 12)
    dFydh[1] = 1.0;
  if (parameterID == 3 || parameterID == 34)
    dEdh[0] = 1.0;
  if (parameterID == 4 || parameterID == 34)
    dEdh[1] = 1.0;
  if (parameterID == 5 || parameterID == 56)
    dHkindh[0] = 1.0;
  if (parameterID == 6 || parameterID == 56)
    dHkindh[1] = 1.0;
  if (parameterID == 7)
    dHisodh = 1.0;

  double dePdh[2]; dePdh[0] = 0.0; dePdh[1] = 0.0;
  double dalphadh = 0.0;
  if (SHVs != 0) {
    dePdh[0] = (*SHVs)(0,gradIndex);
    dePdh[1] = (*SHVs)(1,gradIndex);
    dalphadh = (*SHVs)(2,gradIndex);
  }

  double xsi[2];
  xsi[0] = E[0]*e_n1[0] - (E[0]+Hkin[0])*eP_n1[0];
  xsi[1] = E[1]*e_n1[1] - (E[1]+Hkin[1])*eP_n1[1];

  double Q[2];
  Q[0] = 1.0/(sigY[0]*sigY[0]);
  Q[1] = 1.0/(sigY[1]*sigY[1]);

  double q = sqrt(xsi[0]*Q[0]*xsi[0] + xsi[1]*Q[1]*xsi[1]);
  double F = q - (1 + Hiso*alpha_n1);

  s(0) = dEdh[0]*(e_n1[0]-eP_n1[0]) - E[0]*dePdh[0];
  s(1) = dEdh[1]*(e_n1[1]-eP_n1[1]) - E[1]*dePdh[1];

  double dzdh[2];
  dzdh[0] = s(0) - dHkindh[0]*eP_n1[0] - Hkin[0]*dePdh[0];
  dzdh[1] = s(1) - dHkindh[1]*eP_n1[1] - Hkin[1]*dePdh[1];

  double dg = dg_n1;

  //if (F > -100*DBL_EPSILON) {
  if (F >= -10*DBL_EPSILON) {

    double n[2];
    n[0] = Q[0]*xsi[0]/q;
    n[1] = Q[1]*xsi[1]/q;

    static Matrix J(3,3);

    J(0,0) = 1.0 + dg/q*(E[0]+Hkin[0])*(Q[0]-n[0]*n[0]);
    J(0,1) =       dg/q*(E[0]+Hkin[0])*(    -n[0]*n[1]);
    J(0,2) = (E[0]+Hkin[0])*n[0];

    J(1,0) =       dg/q*(E[1]+Hkin[1])*(    -n[1]*n[0]);
    J(1,1) = 1.0 + dg/q*(E[1]+Hkin[1])*(Q[1]-n[1]*n[1]);
    J(1,2) = (E[1]+Hkin[1])*n[1];

    J(2,0) = n[0];
    J(2,1) = n[1];
    J(2,2) = -Hiso;

    double dQ[2];
    dQ[0] = -2*Q[0]/sigY[0]*dFydh[0];
    dQ[1] = -2*Q[1]/sigY[1]*dFydh[1];

    static Matrix B(2,2);
    B(0,0) = 1.0-0.5/q*n[0]*xsi[0];
    B(0,1) =    -0.5/q*n[0]*xsi[1];
    B(1,0) =    -0.5/q*n[1]*xsi[0];
    B(1,1) = 1.0-0.5/q*n[1]*xsi[1];

    static Vector c(3);
    c(0) = dzdh[0] - (E[0]+Hkin[0])*dg/q*(B(0,0)*dQ[0]*xsi[0] + B(0,1)*dQ[1]*xsi[1]);
    c(1) = dzdh[1] - (E[1]+Hkin[1])*dg/q*(B(1,0)*dQ[0]*xsi[0] + B(1,1)*dQ[1]*xsi[1]);
    c(2) = Hiso*dalphadh + dHisodh*alpha_n1 - 0.5/q*(xsi[0]*dQ[0]*xsi[0]+xsi[1]*dQ[1]*xsi[1]);

    static Vector dx(3);
    J.Solve(c, dx);

    dzdh[0] = dx(0);
    dzdh[1] = dx(1);
    double ddgdh = dx(2);

    double dndh[2];
    dndh[0] = (Q[0]-n[0]*n[0])/q*dx(0) - (n[0]*n[1])/q*dx(1) + B(0,0)/q*dQ[0]*xsi[0] + B(0,1)/q*dQ[1]*xsi[1];
    dndh[1] = (Q[1]-n[1]*n[1])/q*dx(1) - (n[1]*n[0])/q*dx(0) + B(1,0)/q*dQ[0]*xsi[0] + B(1,1)/q*dQ[1]*xsi[1];

    dePdh[0] += ddgdh*n[0] + dg*dndh[0];
    dePdh[1] += ddgdh*n[1] + dg*dndh[1];

    s(0) = dzdh[0] + Hkin[0]*dePdh[0] + dHkindh[0]*eP_n1[0];
    s(1) = dzdh[1] + Hkin[1]*dePdh[1] + dHkindh[1]*eP_n1[1];
  }
  return s;
}
    
int 
Elliptical2::commitSensitivity(const Vector &dedh, int gradIndex, int numGrads)
{
  if (SHVs == 0) {
    SHVs = new Matrix(3,numGrads);
  }

  if (gradIndex >= SHVs->noCols()) {
    //opserr << gradIndex << ' ' << SHVs->noCols() << endln;
    return 0;
  }

  //return 0;

  double dEdh[2];
  dEdh[0] = 0.0;
  dEdh[1] = 0.0;
  double dHkindh[2];
  dHkindh[0] = 0.0;
  dHkindh[1] = 0.0;
  double dHisodh = 0.0;
  double dFydh[2];
  dFydh[0] = 0.0;
  dFydh[1] = 0.0;

  if (parameterID == 1 || parameterID == 12)
    dFydh[0] = 1.0;
  if (parameterID == 2 || parameterID == 12)
    dFydh[1] = 1.0;
  if (parameterID == 3 || parameterID == 34)
    dEdh[0] = 1.0;
  if (parameterID == 4 || parameterID == 34)
    dEdh[1] = 1.0;
  if (parameterID == 5 || parameterID == 56)
    dHkindh[0] = 1.0;
  if (parameterID == 6 || parameterID == 56)
    dHkindh[1] = 1.0;
  if (parameterID == 7)
    dHisodh = 1.0;

  double dePdh[2]; dePdh[0] = 0.0; dePdh[1] = 0.0;
  double dalphadh = 0.0;
  if (SHVs != 0) {
    dePdh[0] = (*SHVs)(0,gradIndex);
    dePdh[1] = (*SHVs)(1,gradIndex);
    dalphadh = (*SHVs)(2,gradIndex);
  }

  double xsi[2];
  xsi[0] = E[0]*e_n1[0] - (E[0]+Hkin[0])*eP_n1[0];
  xsi[1] = E[1]*e_n1[1] - (E[1]+Hkin[1])*eP_n1[1];

  double Q[2];
  Q[0] = 1.0/(sigY[0]*sigY[0]);
  Q[1] = 1.0/(sigY[1]*sigY[1]);

  double q = sqrt(xsi[0]*Q[0]*xsi[0] + xsi[1]*Q[1]*xsi[1]);
  double F = q - (1 + Hiso*alpha_n1);

  double dzdh[2];
  dzdh[0] = dEdh[0]*e_n1[0] + E[0]*dedh(0) - (dEdh[0]+dHkindh[0])*eP_n1[0] - (E[0]+Hkin[0])*dePdh[0];
  dzdh[1] = dEdh[1]*e_n1[1] + E[1]*dedh(1) - (dEdh[1]+dHkindh[1])*eP_n1[1] - (E[1]+Hkin[1])*dePdh[1];

  double dg = dg_n1;

  //if (F > -100*DBL_EPSILON) {
  if (F >= -10*DBL_EPSILON) {

    double n[2];
    n[0] = Q[0]*xsi[0]/q;
    n[1] = Q[1]*xsi[1]/q;

    static Matrix J(3,3);

    J(0,0) = 1.0 + dg/q*(E[0]+Hkin[0])*(Q[0]-n[0]*n[0]);
    J(0,1) =       dg/q*(E[0]+Hkin[0])*(    -n[0]*n[1]);
    J(0,2) = (E[0]+Hkin[0])*n[0];

    J(1,0) =       dg/q*(E[1]+Hkin[1])*(    -n[1]*n[0]);
    J(1,1) = 1.0 + dg/q*(E[1]+Hkin[1])*(Q[1]-n[1]*n[1]);
    J(1,2) = (E[1]+Hkin[1])*n[1];

    J(2,0) = n[0];
    J(2,1) = n[1];
    J(2,2) = -Hiso;
    
    double dQ[2];
    dQ[0] = -2*Q[0]/sigY[0]*dFydh[0];
    dQ[1] = -2*Q[1]/sigY[1]*dFydh[1];

    static Matrix B(2,2);
    B(0,0) = 1.0-0.5/q*n[0]*xsi[0];
    B(0,1) =    -0.5/q*n[0]*xsi[1];
    B(1,0) =    -0.5/q*n[1]*xsi[0];
    B(1,1) = 1.0-0.5/q*n[1]*xsi[1];

    static Vector c(3);
    c(0) = dzdh[0] - (E[0]+Hkin[0])*dg/q*(B(0,0)*dQ[0]*xsi[0] + B(0,1)*dQ[1]*xsi[1]);
    c(1) = dzdh[1] - (E[1]+Hkin[1])*dg/q*(B(1,0)*dQ[0]*xsi[0] + B(1,1)*dQ[1]*xsi[1]);
    c(2) = Hiso*dalphadh + dHisodh*alpha_n1 - 0.5/q*(xsi[0]*dQ[0]*xsi[0]+xsi[1]*dQ[1]*xsi[1]);

    static Vector dx(3);
    J.Solve(c, dx);

    double ddgdh = dx(2);

    double dndh[2];
    dndh[0] = (Q[0]-n[0]*n[0])/q*dx(0) - (n[0]*n[1])/q*dx(1) + B(0,0)/q*dQ[0]*xsi[0] + B(0,1)/q*dQ[1]*xsi[1];
    dndh[1] = (Q[1]-n[1]*n[1])/q*dx(1) - (n[1]*n[0])/q*dx(0) + B(1,0)/q*dQ[0]*xsi[0] + B(1,1)/q*dQ[1]*xsi[1];

    dalphadh += ddgdh;
    dePdh[0] += ddgdh*n[0] + dg*dndh[0];
    dePdh[1] += ddgdh*n[1] + dg*dndh[1];

    (*SHVs)(0,gradIndex) = dePdh[0];
    (*SHVs)(1,gradIndex) = dePdh[1];
    (*SHVs)(2,gradIndex) = dalphadh;
  }

  return 0;
}
