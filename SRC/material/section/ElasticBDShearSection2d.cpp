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
                                                                        
// $Revision: 1.2 $
// $Date: 2008-08-26 16:46:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticBDShearSection2d.cpp,v $

#include <ElasticBDShearSection2d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <stdlib.h>
#include <Information.h>
#include <Parameter.h>

#include <classTags.h>

Vector ElasticBDShearSection2d::s(3);
Matrix ElasticBDShearSection2d::ks(3,3);
ID ElasticBDShearSection2d::code(3);

ElasticBDShearSection2d::ElasticBDShearSection2d(void)
:SectionForceDeformation(0, SEC_TAG_ElasticBDShear2d),
 E(0.0), b(0.0), d(0.0), G(0.0), alpha(0.0),
 e(3), parameterID(0)
{
    if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
	code(2) = SECTION_RESPONSE_VY;	// Vy is the third
    }    
}

ElasticBDShearSection2d::ElasticBDShearSection2d
(int tag, double E_in, double b_in, double d_in, double G_in, double alpha_in)
:SectionForceDeformation(tag, SEC_TAG_ElasticBDShear2d),
 E(E_in), b(b_in), d(d_in), G(G_in), alpha(alpha_in),
 e(3), parameterID(0)
{
    if (E <= 0.0)  {
		opserr << "ElasticBDShearSection2d::ElasticBDShearSection2d -- Input E <= 0.0 ... setting E to 1.0\n";
		E = 1.0;
  }
	
    if (b <= 0.0)  {
		opserr << "ElasticBDShearSection2d::ElasticBDShearSection2d -- Input b <= 0.0 ... setting b to 1.0\n";
		b = 1.0;
    }
    
    if (d <= 0.0)  {
		opserr << "ElasticBDShearSection2d::ElasticBDShearSection2d -- Input d <= 0.0 ... setting d to 1.0\n";
		d = 1.0;
    }    
	
    if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
	code(2) = SECTION_RESPONSE_VY;	// Vy is the third
    }
}

ElasticBDShearSection2d::~ElasticBDShearSection2d(void)
{

}

int 
ElasticBDShearSection2d::commitState(void)
{
  return 0;
}

int 
ElasticBDShearSection2d::revertToLastCommit(void)
{
  return 0;
}

int 
ElasticBDShearSection2d::revertToStart(void)
{
  return 0;
}

int
ElasticBDShearSection2d::setTrialSectionDeformation (const Vector &def)
{
  e = def;

  return 0;
}

const Vector &
ElasticBDShearSection2d::getSectionDeformation (void)
{
  return e;
}

const Vector &
ElasticBDShearSection2d::getStressResultant (void)
{
  double A = b*d;
  double I = b*d*d*d/12.0;

  s(0) = E*A*e(0);
  s(1) = E*I*e(1);    
  s(2) = G*A*alpha*e(2);

  return s;
}

const Matrix &
ElasticBDShearSection2d::getSectionTangent(void)
{
  double A = b*d;
  double I = b*d*d*d/12.0;

  ks(0,0) = E*A;
  ks(1,1) = E*I;
  ks(2,2) = G*A*alpha;
  
  return ks;
}

const Matrix &
ElasticBDShearSection2d::getInitialTangent(void)
{
  double A = b*d;
  double I = b*d*d*d/12.0;

  ks(0,0) = E*A;
  ks(1,1) = E*I;
  ks(2,2) = G*A*alpha;
  
  return ks;
}

const Matrix &
ElasticBDShearSection2d::getSectionFlexibility (void)
{
  double A = b*d;
  double I = b*d*d*d/12.0;

  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*I);
  ks(2,2) = 1.0/(G*A*alpha);
  
  return ks;
}

const Matrix &
ElasticBDShearSection2d::getInitialFlexibility(void)
{
  double A = b*d;
  double I = b*d*d*d/12.0;

  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*I);
  ks(2,2) = 1.0/(G*A*alpha);
  
  return ks;
}

SectionForceDeformation*
ElasticBDShearSection2d::getCopy(void)
{
  ElasticBDShearSection2d *theCopy =
    new ElasticBDShearSection2d (this->getTag(), E, b, d, G, alpha);
  
  theCopy->parameterID = parameterID;
  
  return theCopy;
}

const ID&
ElasticBDShearSection2d::getType(void)
{
  return code;
}

int
ElasticBDShearSection2d::getOrder(void) const
{
  return 3;
}

int
ElasticBDShearSection2d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(6);
  
  int dataTag = this->getDbTag();
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = b;
  data(3) = d;    
  data(4) = G;    
  data(5) = alpha;    
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res<0) {
    opserr << "ElasticBDShearSection2d::sendSelf -- failed to send data\n";
    return res;
  }
  
  return res;
}

int
ElasticBDShearSection2d::recvSelf(int commitTag, Channel &theChannel,
				FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(6);
  
  int dataTag = this->getDbTag();
  
  res += theChannel.recvVector(dataTag, commitTag, data);
  if(res < 0) {
    opserr << "ElasticBDShearSection2d::recvSelf -- failed to receive data\n";
    return res;
  }
  
  this->setTag((int)data(0));
  E = data(1);
  b = data(2);
  d = data(3);
  G = data(4);
  alpha = data(5);
  
  return res;
}
 
void
ElasticBDShearSection2d::Print(OPS_Stream &s, int flag)
{
  s << "ElasticBDShearSection2d, tag: " << this->getTag() << endln;
  s << "\tE: " << E << endln;
  s << "\tb: " << b << endln;
  s << "\td: " << d << endln;
  s << "\tG: " << G << endln;
  s << "\talpha: " << alpha << endln;
}

int
ElasticBDShearSection2d::setParameter(const char **argv, int argc,
				      Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"b") == 0) {
    param.setValue(b);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"d") == 0) {
    param.setValue(d);
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
ElasticBDShearSection2d::updateParameter(int paramID, Information &info)
{
  if (paramID == 1)
    E = info.theDouble;
  if (paramID == 2)
    b = info.theDouble;
  if (paramID == 3)
    d = info.theDouble;
  if (paramID == 4)
    G = info.theDouble;
  if (paramID == 5)
    alpha = info.theDouble;

  return 0;
}

int
ElasticBDShearSection2d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
ElasticBDShearSection2d::getStressResultantSensitivity(int gradIndex,
						       bool conditional)
{
  static Vector dsdh(3);
  dsdh.Zero();

  double A = b*d;
  double I = b*d*d*d/12.0;

  if (parameterID == 1) { // E
    dsdh(0) = A*e(0);
    dsdh(1) = I*e(1);
  }
  if (parameterID == 2) { // b
    dsdh(0) = E*d*e(0);
    dsdh(1) = E*d*d*d/12.0*e(1);
    dsdh(2) = G*alpha*d*e(2);
  }
  if (parameterID == 3) { // d
    dsdh(0) = E*b*e(0);
    dsdh(1) = E*b*d*d/4.0*e(1);
    dsdh(2) = G*alpha*b*e(2);
  }
  if (parameterID == 4) // G
    dsdh(2) = A*alpha*e(2);
  if (parameterID == 5) // alpha
    dsdh(2) = G*A*e(2);

  //opserr << "EBD dsdh - " << parameterID << ' ' << dsdh << ' ' << e;

  return dsdh;
}

const Matrix&
ElasticBDShearSection2d::getInitialTangentSensitivity(int gradIndex)
{
  ks.Zero();

  return ks;
}
