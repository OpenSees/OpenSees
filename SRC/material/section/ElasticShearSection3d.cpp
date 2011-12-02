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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticShearSection3d.cpp,v $

#include <ElasticShearSection3d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <Parameter.h>
#include <stdlib.h>
#include <string.h>

#include <classTags.h>

Vector ElasticShearSection3d::s(6);
Matrix ElasticShearSection3d::ks(6,6);
ID ElasticShearSection3d::code(6);

ElasticShearSection3d::ElasticShearSection3d(void)
:SectionForceDeformation(0, SEC_TAG_Elastic3d),
 E(0.0), A(0.0), Iz(0.0), Iy(0.0), G(0.0), J(0.0), alpha(0.0),
 e(6), eCommit(6)
{
    if (code(0) != SECTION_RESPONSE_P) {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
	code(2) = SECTION_RESPONSE_VY;	// Vy is the third 
	code(3) = SECTION_RESPONSE_MY;	// My is the fourth
	code(4) = SECTION_RESPONSE_VZ;	// Vz is the fifth
	code(5) = SECTION_RESPONSE_T;	// T is the sixth
    }
}

ElasticShearSection3d::ElasticShearSection3d
(int tag, double E_in, double A_in, double Iz_in, double Iy_in,
 double G_in, double J_in, double alpha_in)
:SectionForceDeformation(tag, SEC_TAG_Elastic3d),
 E(E_in), A(A_in), Iz(Iz_in), Iy(Iy_in), G(G_in), J(J_in), alpha(alpha_in),
 e(6), eCommit(6)
{
    if (E <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input E <= 0.0 ... setting E to 1.0\n";
      E = 1.0;
    }
    
    if (A <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input A <= 0.0 ... setting A to 1.0\n";
      A = 1.0;
    }

    if (Iz <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input Iz <= 0.0 ... setting Iz to 1.0\n";
      Iz = 1.0;
    }
    
    if (Iy <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input Iy <= 0.0 ... setting Iy to 1.0\n";
      Iy = 1.0;
    }

    if (G <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input G <= 0.0 ... setting G to 1.0\n";
      G = 1.0;
    }
    
    if (J <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input J <= 0.0 ... setting J to 1.0\n";
      J = 1.0;
    }
    
    if (code(0) != SECTION_RESPONSE_P) {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
	code(2) = SECTION_RESPONSE_VY;	// Vy is the third 
	code(3) = SECTION_RESPONSE_MY;	// My is the fourth
	code(4) = SECTION_RESPONSE_VZ;	// Vz is the fifth
	code(5) = SECTION_RESPONSE_T;	// T is the sixth
    }
}

ElasticShearSection3d::~ElasticShearSection3d(void)
{
    return;
}

int 
ElasticShearSection3d::commitState(void)
{
  eCommit = e;

  return 0;
}

int 
ElasticShearSection3d::revertToLastCommit(void)
{
  e = eCommit;

  return 0;
}

int 
ElasticShearSection3d::revertToStart(void)
{
  eCommit.Zero();

  return 0;
}

int
ElasticShearSection3d::setTrialSectionDeformation (const Vector &def)
{
  e = def;
    
  return 0;
}

const Vector &
ElasticShearSection3d::getSectionDeformation (void)
{
  return e;
}

const Vector &
ElasticShearSection3d::getStressResultant (void)
{
  s(0) = E*A*e(0);
  s(1) = E*Iz*e(1);
  s(3) = E*Iy*e(3);
  s(5) = G*J*e(5);

  double GA = G*A*alpha;
  s(2) = GA*e(2);
  s(4) = GA*e(4);
  
  return s;
}

const Matrix &
ElasticShearSection3d::getSectionTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*Iz;
  ks(3,3) = E*Iy;
  ks(5,5) = G*J;
  
  double GA = G*A*alpha;
  ks(2,2) = GA;
  ks(4,4) = GA;
  
  return ks;
}

const Matrix &
ElasticShearSection3d::getInitialTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*Iz;
  ks(3,3) = E*Iy;
  ks(5,5) = G*J;
  
  double GA = G*A*alpha;
  ks(2,2) = GA;
  ks(4,4) = GA;
  
  return ks;
}

const Matrix &
ElasticShearSection3d::getSectionFlexibility (void)
{
  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*Iz);
  ks(3,3) = 1.0/(E*Iy);
  ks(5,5) = 1.0/(G*J);

  double GA = 1.0/(G*A*alpha);
  ks(2,2) = GA;
  ks(4,4) = GA;

  return ks;
}

const Matrix &
ElasticShearSection3d::getInitialFlexibility (void)
{
  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*Iz);
  ks(3,3) = 1.0/(E*Iy);
  ks(5,5) = 1.0/(G*J);

  double GA = 1.0/(G*A*alpha);
  ks(2,2) = GA;
  ks(4,4) = GA;

  return ks;
}

SectionForceDeformation*
ElasticShearSection3d::getCopy ()
{
    // Make a copy of the hinge
    ElasticShearSection3d *theCopy =
	new ElasticShearSection3d (this->getTag(), E, A, Iz, Iy, G, J, alpha);

    theCopy->eCommit = eCommit;

    return theCopy;
}

const ID&
ElasticShearSection3d::getType ()
{
    return code;
}

int
ElasticShearSection3d::getOrder () const
{
    return 6;
}

int
ElasticShearSection3d::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(14);

    int dataTag = this->getDbTag();
    
    data(0) = this->getTag();
    data(1) = E;
    data(2) = A;    
    data(3) = Iz;
    data(4) = Iy;
    data(5) = G;
    data(6) = J;
    data(7) = alpha;
    data(8) = eCommit(0);
    data(9) = eCommit(1);
    data(10) = eCommit(2);
    data(11) = eCommit(3);
    data(12) = eCommit(4);
    data(13) = eCommit(5);
    
    res += theChannel.sendVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticShearSection3d::sendSelf -- failed to send data\n";
      return res;
    }
    
    return res;
}

int
ElasticShearSection3d::recvSelf(int commitTag, Channel &theChannel,
					 FEM_ObjectBroker &theBroker)
{
    int res = 0;
    
    static Vector data(14);

    int dataTag = this->getDbTag();

    res += theChannel.recvVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticShearSection3d::recvSelf -- failed to receive data\n";
      return res;
    }

    this->setTag((int)data(0));
    E = data(1);
    A = data(2);    
    Iz = data(3);
    Iy = data(4);
    G = data(5);
    J = data(6); 
    alpha = data(7);   
    eCommit(0) = data(8);
    eCommit(1) = data(9);
    eCommit(2) = data(10);
    eCommit(3) = data(11);
    eCommit(4) = data(12);
    eCommit(5) = data(13);

    return res;
}
 
void
ElasticShearSection3d::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {

  } else {
    s << "ElasticShearSection3d, tag: " << this->getTag() << endln;
    s << "\t E: " << E << endln;
    s << "\t A: " << A << endln;
    s << "\tIz: " << Iz << endln;
    s << "\tIy: " << Iy << endln;
    s << "\t G: " << G << endln;
    s << "\t J: " << J << endln;
    s << "\talpha: " << alpha << endln;
  }
}

int
ElasticShearSection3d::setParameter(const char **argv, int argc, Parameter &param)
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
  if (strcmp(argv[0],"Iz") == 0) {
	  param.setValue(Iz);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"Iy") == 0) {
	  param.setValue(Iy);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"G") == 0) {
	  param.setValue(G);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"J") == 0) {
	  param.setValue(J);
    return param.addObject(6, this);
  }
  if (strcmp(argv[0],"alpha") == 0) {
	  param.setValue(alpha);
    return param.addObject(7, this);
  }
  return -1;
}

int
ElasticShearSection3d::updateParameter(int paramID, Information &info)
{
  if (paramID == 1)
    E = info.theDouble;
  if (paramID == 2)
    A = info.theDouble;
  if (paramID == 3)
    Iz = info.theDouble;
  if (paramID == 4)
    Iy = info.theDouble;
  if (paramID == 5)
    G = info.theDouble;
  if (paramID == 6)
    J = info.theDouble;
  if (paramID == 7)
    alpha = info.theDouble;

  return 0;
}

int
ElasticShearSection3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
ElasticShearSection3d::getStressResultantSensitivity(int gradIndex,
						     bool conditional)
{
  s.Zero();

  if (parameterID == 1) { // E
    s(0) = A*e(0);
    s(1) = Iz*e(1);
    s(3) = Iy*e(3);
  }
  if (parameterID == 2) { // A
    s(0) = E*e(0);
    double Galpha = G*alpha;
    s(2) = Galpha*e(2);
    s(4) = Galpha*e(4);
  }
  if (parameterID == 3) // Iz
    s(1) = E*e(1);
  if (parameterID == 4) // Iy
    s(3) = E*e(3);
  if (parameterID == 5) { // G
    double Aalpha = A*alpha;
    s(2) = Aalpha*e(2);
    s(4) = Aalpha*e(4);
    s(5) = J*e(5);
  }
  if (parameterID == 6) // J
    s(5) = G*e(5);
  if (parameterID == 7) { // alpha
    double GA = G*A;
    s(2) = GA*e(2);
    s(4) = GA*e(4);
  }

  return s;
}

const Matrix&
ElasticShearSection3d::getInitialTangentSensitivity(int gradIndex)
{
  ks.Zero();

  return ks;
}
