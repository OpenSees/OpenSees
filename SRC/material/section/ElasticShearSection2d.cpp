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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-08-26 16:49:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticShearSection2d.cpp,v $

#include <ElasticShearSection2d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <stdlib.h>
#include <Information.h>
#include <Parameter.h>

#include <classTags.h>

Vector ElasticShearSection2d::s(3);
Matrix ElasticShearSection2d::ks(3,3);
ID ElasticShearSection2d::code(3);

ElasticShearSection2d::ElasticShearSection2d(void)
:SectionForceDeformation(0, SEC_TAG_ElasticShear2d),
 E(0.0), A(0.0), I(0.0), G(0.0), alpha(0.0),
 e(3), eCommit(3), parameterID(0)
{
    if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
	code(2) = SECTION_RESPONSE_VY;	// Vy is the third
    }    
}

ElasticShearSection2d::ElasticShearSection2d
(int tag, double E_in, double A_in, double I_in, double G_in, double alpha_in)
:SectionForceDeformation(tag, SEC_TAG_ElasticShear2d),
 E(E_in), A(A_in), I(I_in), G(G_in), alpha(alpha_in),
 e(3), eCommit(3), parameterID(0)
{
    if (E <= 0.0)  {
		opserr << "ElasticShearSection2d::ElasticShearSection2d -- Input E <= 0.0 ... setting E to 1.0\n";
		E = 1.0;
  }
	
    if (A <= 0.0)  {
		opserr << "ElasticShearSection2d::ElasticShearSection2d -- Input A <= 0.0 ... setting A to 1.0\n";
		A = 1.0;
    }
    
    if (I <= 0.0)  {
		opserr << "ElasticShearSection2d::ElasticShearSection2d -- Input I <= 0.0 ... setting I to 1.0\n";
		I = 1.0;
    }    
	
    if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
	code(2) = SECTION_RESPONSE_VY;	// Vy is the third
    }
}

ElasticShearSection2d::~ElasticShearSection2d(void)
{
  return;
}

int 
ElasticShearSection2d::commitState(void)
{
  eCommit = e;

  return 0;
}

int 
ElasticShearSection2d::revertToLastCommit(void)
{
  e = eCommit;

  return 0;
}

int 
ElasticShearSection2d::revertToStart(void)
{
  eCommit.Zero();

  return 0;
}

int
ElasticShearSection2d::setTrialSectionDeformation (const Vector &def)
{
  e = def;

  return 0;
}

const Vector &
ElasticShearSection2d::getSectionDeformation (void)
{
  return e;
}

const Vector &
ElasticShearSection2d::getStressResultant (void)
{
  s(0) = E*A*e(0);
  s(1) = E*I*e(1);    
  s(2) = G*A*alpha*e(2);

  return s;
}

const Matrix &
ElasticShearSection2d::getSectionTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*I;
  ks(2,2) = G*A*alpha;
  
  return ks;
}

const Matrix &
ElasticShearSection2d::getInitialTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*I;
  ks(2,2) = G*A*alpha;
  
  return ks;
}

const Matrix &
ElasticShearSection2d::getSectionFlexibility (void)
{
  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*I);
  ks(2,2) = 1.0/(G*A*alpha);
  
  return ks;
}

const Matrix &
ElasticShearSection2d::getInitialFlexibility(void)
{
  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*I);
  ks(2,2) = 1.0/(G*A*alpha);
  
  return ks;
}

SectionForceDeformation*
ElasticShearSection2d::getCopy(void)
{
  ElasticShearSection2d *theCopy =
    new ElasticShearSection2d (this->getTag(), E, A, I, G, alpha);
  
  theCopy->eCommit = eCommit;
  theCopy->parameterID = parameterID;
  
  return theCopy;
}

const ID&
ElasticShearSection2d::getType(void)
{
  return code;
}

int
ElasticShearSection2d::getOrder(void) const
{
  return 3;
}

int
ElasticShearSection2d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(9);
  
  int dataTag = this->getDbTag();
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = A;
  data(3) = I;    
  data(4) = G;    
  data(5) = alpha;    
  data(6) = eCommit(0);
  data(7) = eCommit(1);
  data(8) = eCommit(2);
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res<0) {
    opserr << "ElasticShearSection2d::sendSelf -- failed to send data\n";
    return res;
  }
  
  return res;
}

int
ElasticShearSection2d::recvSelf(int commitTag, Channel &theChannel,
				FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(9);
  
  int dataTag = this->getDbTag();
  
  res += theChannel.recvVector(dataTag, commitTag, data);
  if(res < 0) {
    opserr << "ElasticShearSection2d::recvSelf -- failed to receive data\n";
    return res;
  }
  
  this->setTag((int)data(0));
  E = data(1);
  A = data(2);
  I = data(3);
  G = data(4);
  alpha = data(5);
  eCommit(0) = data(6);
  eCommit(1) = data(7);
  eCommit(2) = data(8);
  
  return res;
}
 
void
ElasticShearSection2d::Print(OPS_Stream &s, int flag)
{
  s << "ElasticShearSection2d, tag: " << this->getTag() << endln;
  s << "\tE: " << E << endln;
  s << "\tA: " << A << endln;
  s << "\tI: " << I << endln;
  s << "\tG: " << G << endln;
  s << "\talpha: " << alpha << endln;
}

int
ElasticShearSection2d::setParameter(const char **argv, int argc,
				    Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"E") == 0) {
	  param.setValue(E);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"A") == 0) {
	  param.setValue(A);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"I") == 0) {
	  param.setValue(I);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"G") == 0) {
	  param.setValue(G);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"alpha") == 0) {
	  param.setValue(alpha);
    return param.addObject(5, this);
  }
  return -1;
}

int
ElasticShearSection2d::updateParameter(int paramID, Information &info)
{
  if (paramID == 1)
    E = info.theDouble;
  if (paramID == 2)
    A = info.theDouble;
  if (paramID == 3)
    I = info.theDouble;
  if (paramID == 4)
    G = info.theDouble;
  if (paramID == 5)
    alpha = info.theDouble;

  return 0;
}

int
ElasticShearSection2d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
ElasticShearSection2d::getStressResultantSensitivity(int gradIndex,
						     bool conditional)
{
  s.Zero();

  if (parameterID == 1) { // E
    s(0) = A*e(0);
    s(1) = I*e(1);
  }
  if (parameterID == 2) { // A
    s(0) = E*e(0);
    s(2) = G*alpha*e(2);
  }
  if (parameterID == 3) // I
    s(1) = E*e(1);
  if (parameterID == 4) // G
    s(2) = A*alpha*e(2);
  if (parameterID == 5) // alpha
    s(2) = G*A*e(2);

  return s;
}

const Matrix&
ElasticShearSection2d::getInitialTangentSensitivity(int gradIndex)
{
  ks.Zero();

  return ks;
}
