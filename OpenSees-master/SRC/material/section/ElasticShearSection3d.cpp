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
#include <elementAPI.h>

Vector ElasticShearSection3d::s(6);
Matrix ElasticShearSection3d::ks(6,6);
ID ElasticShearSection3d::code(6);

void* OPS_ElasticShearSection3d()
{
    if(OPS_GetNumRemainingInputArgs() < 9) {
	opserr<<"insufficient arguments for ealstic shear 3d section\n";
	return 0;
    }

    // get tag
    int tag;
    int numData = 1;
    if(OPS_GetIntInput(&numData,&tag) < 0) return 0;

    // get data
    numData = 8;
    double data[8];
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;

    return new ElasticShearSection3d(tag,data[0],data[1],data[2],
				     data[3],data[4],data[5],data[6],
				     data[7]);
}


ElasticShearSection3d::ElasticShearSection3d(void)
:SectionForceDeformation(0, SEC_TAG_ElasticShear3d),
 E(0.0), A(0.0), Iz(0.0), Iy(0.0), G(0.0), J(0.0),
 alphaY(0.0), alphaZ(0.0), e(6)
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
 double G_in, double J_in, double alphaY_in, double alphaZ_in)
:SectionForceDeformation(tag, SEC_TAG_ElasticShear3d),
 E(E_in), A(A_in), Iz(Iz_in), Iy(Iy_in), G(G_in), J(J_in),
 alphaY(alphaY_in), alphaZ(alphaZ_in), e(6)
{
    if (E <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input E <= 0.0\n";
    }
    
    if (A <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input A <= 0.0\n";
    }

    if (Iz <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input Iz <= 0.0\n";
    }
    
    if (Iy <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input Iy <= 0.0\n";
    }

    if (G <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input G <= 0.0\n";
    }
    
    if (J <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input J <= 0.0\n";
    }
    
    if (alphaY <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input alphaY <= 0.0\n";
    }

    if (alphaZ <= 0.0)  {
      opserr << "ElasticShearSection3d::ElasticShearSection3d -- Input alphaZ <= 0.0\n";
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
  return 0;
}

int 
ElasticShearSection3d::revertToLastCommit(void)
{
  return 0;
}

int 
ElasticShearSection3d::revertToStart(void)
{
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

  double GA = G*A;
  s(2) = GA*alphaY*e(2);
  s(4) = GA*alphaZ*e(4);
  
  return s;
}

const Matrix &
ElasticShearSection3d::getSectionTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*Iz;
  ks(3,3) = E*Iy;
  ks(5,5) = G*J;
  
  double GA = G*A;
  ks(2,2) = GA*alphaY;
  ks(4,4) = GA*alphaZ;
  
  return ks;
}

const Matrix &
ElasticShearSection3d::getInitialTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*Iz;
  ks(3,3) = E*Iy;
  ks(5,5) = G*J;
  
  double GA = G*A;
  ks(2,2) = GA*alphaY;
  ks(4,4) = GA*alphaZ;
  
  return ks;
}

const Matrix &
ElasticShearSection3d::getSectionFlexibility (void)
{
  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*Iz);
  ks(3,3) = 1.0/(E*Iy);
  ks(5,5) = 1.0/(G*J);

  double oneOverGA = 1.0/(G*A);
  ks(2,2) = oneOverGA/alphaY;
  ks(4,4) = oneOverGA/alphaZ;

  return ks;
}

const Matrix &
ElasticShearSection3d::getInitialFlexibility (void)
{
  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*Iz);
  ks(3,3) = 1.0/(E*Iy);
  ks(5,5) = 1.0/(G*J);

  double oneOverGA = 1.0/(G*A);
  ks(2,2) = oneOverGA/alphaY;
  ks(4,4) = oneOverGA/alphaZ;

  return ks;
}

SectionForceDeformation*
ElasticShearSection3d::getCopy ()
{
    // Make a copy of the hinge
    ElasticShearSection3d *theCopy =
	new ElasticShearSection3d (this->getTag(), E, A, Iz, Iy,
                               G, J, alphaY, alphaZ);

    theCopy->parameterID = parameterID;

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

    static Vector data(9);

    int dataTag = this->getDbTag();
    
    data(0) = this->getTag();
    data(1) = E;
    data(2) = A;    
    data(3) = Iz;
    data(4) = Iy;
    data(5) = G;
    data(6) = J;
    data(7) = alphaY;
    data(8) = alphaZ;
    
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
    
    static Vector data(9);

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
    alphaY = data(7);
    alphaZ = data(8);

    return res;
}
 
void
ElasticShearSection3d::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
		s << "ElasticShearSection3d, tag: " << this->getTag() << endln;
		s << "\t E: " << E << endln;
		s << "\t A: " << A << endln;
		s << "\tIz: " << Iz << endln;
		s << "\tIy: " << Iy << endln;
		s << "\t G: " << G << endln;
		s << "\t J: " << J << endln;
		s << "\talphaY: " << alphaY << endln;
		s << "\talphaZ: " << alphaZ << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"ElasticShearSection3d\", ";
		s << "\"E\": " << E << ", ";
		s << "\"G\": " << G << ", ";
		s << "\"A\": " << A << ", ";
		s << "\"Avy\": " << alphaY*A << ", ";
		s << "\"Avz\": " << alphaZ*A << ", ";
		s << "\"Jx\": " << J << ", ";
		s << "\"Iy\": " << Iy << ", ";
		s << "\"Iz\": " << Iz << "}";
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
  if (strcmp(argv[0],"alphaY") == 0) {
	  param.setValue(alphaY);
    return param.addObject(7, this);
  }
  if (strcmp(argv[0],"alphaZ") == 0) {
	  param.setValue(alphaZ);
    return param.addObject(8, this);
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
    alphaY = info.theDouble;
  if (paramID == 8)
    alphaZ = info.theDouble;

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
    s(2) = G*alphaY*e(2);
    s(4) = G*alphaZ*e(4);
  }
  if (parameterID == 3) // Iz
    s(1) = E*e(1);
  if (parameterID == 4) // Iy
    s(3) = E*e(3);
  if (parameterID == 5) { // G
    s(2) = A*alphaY*e(2);
    s(4) = A*alphaZ*e(4);
    s(5) = J*e(5);
  }
  if (parameterID == 6) // J
    s(5) = G*e(5);
  if (parameterID == 7) // alphaY
    s(2) = G*A*e(2);
  if (parameterID == 8) // alphaZ
    s(4) = G*A*e(4);

  return s;
}

const Matrix&
ElasticShearSection3d::getInitialTangentSensitivity(int gradIndex)
{
  ks.Zero();

  return ks;
}
