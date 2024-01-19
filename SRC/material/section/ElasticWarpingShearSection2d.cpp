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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticWarpingShearSection2d.cpp,v $

#include <ElasticWarpingShearSection2d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <stdlib.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <classTags.h>
#include <elementAPI.h>

void* OPS_ElasticWarpingShearSection2d()
{
    if (OPS_GetNumRemainingInputArgs() < 9) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section ElasticWarpingShear tag? E? A? Iz? G? alpha? J? B? C?>" << endln;
	return 0;
    }
	
    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid section ElasticWarpingShearSection2d tag" << endln;
	return 0;
    }

    numdata = 8;
    double data[8];
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid double inputs" << endln;
	opserr << "ElasticWarpingShearSection2d section: " << tag << endln;	    
	return 0;
    }
    double E = data[0];
    double A = data[1];
    double Iz = data[2];
    double G = data[3];
    double alpha = data[4];
    double J = data[5];
    double B = data[6];
    double C = data[7];
      
    return new ElasticWarpingShearSection2d(tag, E, A, Iz, G, alpha, J, B, C);

}

Vector ElasticWarpingShearSection2d::s(5);
Matrix ElasticWarpingShearSection2d::ks(5,5);
ID ElasticWarpingShearSection2d::code(5);

ElasticWarpingShearSection2d::ElasticWarpingShearSection2d(void)
:SectionForceDeformation(0, SEC_TAG_ElasticWarpingShear2d),
 E(0.0), A(0.0), I(0.0), G(0.0), alpha(0.0), J(0.0), B(0.0), C(0.0),
 e(5), eCommit(5), parameterID(0)
{
    if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
	code(2) = SECTION_RESPONSE_VY;	// Vy is the third
	code(3) = SECTION_RESPONSE_R;	// R (warping Shear) is the fourth
	code(4) = SECTION_RESPONSE_Q;	// Q (warping Moment) is the fifth

    }    
}

ElasticWarpingShearSection2d::ElasticWarpingShearSection2d
(int tag, double E_in, double A_in, double I_in, double G_in, double alpha_in, double J_in, double B_in, double C_in)
:SectionForceDeformation(tag, SEC_TAG_ElasticWarpingShear2d),
 E(E_in), A(A_in), I(I_in), G(G_in), alpha(alpha_in), J(J_in), B(B_in), C(C_in),
 e(5), eCommit(5), parameterID(0)
{

    if (E <= 0.0)  {
		opserr << "ElasticWarpingShearSection2d::ElasticWarpingShearSection2d -- Input E <= 0.0";
  }
	
    if (A <= 0.0)  {
		opserr << "ElasticWarpingShearSection2d::ElasticWarpingShearSection2d -- Input A <= 0.0";
    }
    
    if (I <= 0.0)  {
		opserr << "ElasticWarpingShearSection2d::ElasticWarpingShearSection2d -- Input I <= 0.0";
    }    
	
    if (G <= 0.0)  {
		opserr << "ElasticWarpingShearSection2d::ElasticWarpingShearSection2d -- Input G <= 0.0";
    }    

	if (alpha <= 0.0)  {
		opserr << "ElasticWarpingShearSection2d::ElasticWarpingShearSection2d -- Input alpha <= 0.0";
	}
	
	if (J <= 0.0)  {
		opserr << "ElasticWarpingShearSection2d::ElasticWarpingShearSection2d -- Input J <= 0.0";
	}

	if (B <= 0.0)  {
		opserr << "ElasticWarpingShearSection2d::ElasticWarpingShearSection2d -- Input B <= 0.0";
	}

	if (C <= 0.0)  {
		opserr << "ElasticWarpingShearSection2d::ElasticWarpingShearSection2d -- Input C <= 0.0";
    }    

		if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
	code(2) = SECTION_RESPONSE_VY;	// Vy is the third
	code(3) = SECTION_RESPONSE_R;	// R (warping Shear) is the fourth
	code(4) = SECTION_RESPONSE_Q;	// Q (warping Moment) is the fifth
    }
}

ElasticWarpingShearSection2d::~ElasticWarpingShearSection2d(void)
{
  return;
}

int 
ElasticWarpingShearSection2d::commitState(void)
{
  eCommit = e;

  return 0;
}

int 
ElasticWarpingShearSection2d::revertToLastCommit(void)
{
  e = eCommit;

  return 0;
}

int 
ElasticWarpingShearSection2d::revertToStart(void)
{
  eCommit.Zero();

  return 0;
}

int
ElasticWarpingShearSection2d::setTrialSectionDeformation (const Vector &def)
{
  e = def;

  return 0;
}

const Vector &
ElasticWarpingShearSection2d::getSectionDeformation (void)
{
  return e;
}

const Vector &
ElasticWarpingShearSection2d::getStressResultant (void)
{
  s(0) = E*A*e(0);
  s(1) = E*I*e(1);    
  s(2) = G*A*alpha*e(2) + G*B*e(3);
  s(3) = G*B*e(2) + G*C*e(3);
  s(4) = E*J*e(4);
  return s;
}

const Matrix &
ElasticWarpingShearSection2d::getSectionTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*I;
  ks(2,2) = G*A*alpha;
  ks(2,3) = G*B;
  ks(3,2) = G*B;
  ks(3,3) = G*C;
  ks(4,4) = E*J;
  
 return ks;
}

const Matrix &
ElasticWarpingShearSection2d::getInitialTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*I;
  ks(2,2) = G*A*alpha;
  ks(2,3) = G*B;
  ks(3,2) = G*B;
  ks(3,3) = G*C;
  ks(4,4) = E*J;
  
 return ks;
}

const Matrix &
ElasticWarpingShearSection2d::getSectionFlexibility (void)
{
	double dtrmnt=0.0;

	dtrmnt =   G * (A * C * alpha - B*B);

    ks(0,0) =  1.0 / (E*A);
	ks(1,1) =  1.0 / (E*I);
    ks(2,2) =  C / dtrmnt;
	ks(2,3) = -B / dtrmnt;
	ks(3,2) = -B / dtrmnt;
	ks(3,3) =  A*alpha / dtrmnt;
	ks(4,4) =  1.0 / (E*J);

  return ks;
}

const Matrix &
ElasticWarpingShearSection2d::getInitialFlexibility(void)
{
	double dtrmnt=0.0;

	dtrmnt =   G * (A * C * alpha - B*B);

    ks(0,0) =  1.0 / (E*A);
	ks(1,1) =  1.0 / (E*I);
    ks(2,2) =  C / dtrmnt;
	ks(2,3) = -B / dtrmnt;
	ks(3,2) = -B / dtrmnt;
	ks(3,3) =  A*alpha / dtrmnt;
	ks(4,4) =  1.0 / (E*J);

  return ks;
}

SectionForceDeformation*
ElasticWarpingShearSection2d::getCopy(void)
{
  ElasticWarpingShearSection2d *theCopy =
    new ElasticWarpingShearSection2d (this->getTag(), E, A, I, G, alpha, J, B, C);
  
  theCopy->eCommit = eCommit;
  theCopy->parameterID = parameterID;
  
  return theCopy;
}

const ID&
ElasticWarpingShearSection2d::getType(void)
{
  return code;
}

int
ElasticWarpingShearSection2d::getOrder(void) const
{
  return 5;
}

int
ElasticWarpingShearSection2d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(14);
  
  int dataTag = this->getDbTag();
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = A;
  data(3) = I;    
  data(4) = G;    
  data(5) = alpha;  
  data(6) = J;
  data(7) = B;
  data(8) = C;
  data(9) =  eCommit(0);
  data(10) = eCommit(1);
  data(11) = eCommit(2);
  data(12) = eCommit(3);
  data(13) = eCommit(4);
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res<0) {
    opserr << "ElasticWarpingShearSection2d::sendSelf -- failed to send data\n";
    return res;
  }
  
  return res;
}

int
ElasticWarpingShearSection2d::recvSelf(int commitTag, Channel &theChannel,
				FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(14);
  
  int dataTag = this->getDbTag();
  
  res += theChannel.recvVector(dataTag, commitTag, data);
  if(res < 0) {
    opserr << "ElasticWarpingShearSection2d::recvSelf -- failed to receive data\n";
    return res;
  }
  
  this->setTag((int)data(0));
  E = data(1);
  A = data(2);
  I = data(3);
  G = data(4);
  alpha = data(5);
  J = data(6);
  B = data(7);
  C= data(8);
  eCommit(0) = data(9);
  eCommit(1) = data(10);
  eCommit(2) = data(11);
  eCommit(3) = data(12);
  eCommit(4) = data(13);
  
  return res;
}
 
void
ElasticWarpingShearSection2d::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
        s << "ElasticWarpingShearSection2d, tag: " << this->getTag() << endln;
        s << "\tE: " << E << endln;
        s << "\tA: " << A << endln;
        s << "\tI: " << I << endln;
        s << "\tG: " << G << endln;
        s << "\talpha: " << alpha << endln;
        s << "\tJ: " << J << endln;
        s << "\tB: " << B << endln;
        s << "\tC: " << C << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ElasticWarpingShearSection2d\", ";
        s << "\"E\": " << E << ", ";
        s << "\"G\": " << G << ", ";
        s << "\"A\": " << A << ", ";
        s << "\"I\": " << I << ", ";
        s << "\"J\": " << J << ", ";
        s << "\"B\": " << B << ", ";
        s << "\"C\": " << C << ", ";
        s << "\"alpha\": " << alpha << "}";
    }
}

int
ElasticWarpingShearSection2d::setParameter(const char **argv, int argc,
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
  if (strcmp(argv[0],"J") == 0) {
	  param.setValue(J);
    return param.addObject(6, this);
  }
  if (strcmp(argv[0],"B") == 0) {
	  param.setValue(B);
    return param.addObject(7, this);
  }
  if (strcmp(argv[0],"C") == 0) {
	  param.setValue(C);
    return param.addObject(8, this);
  }
  return -1;
}

int
ElasticWarpingShearSection2d::updateParameter(int paramID, Information &info)
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
  if (paramID == 6)
    J = info.theDouble;
  if (paramID == 7)
    B = info.theDouble;
  if (paramID == 8)
    C = info.theDouble;

  return 0;
}

int
ElasticWarpingShearSection2d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
ElasticWarpingShearSection2d::getStressResultantSensitivity(int gradIndex,
						bool conditional)
{
  s.Zero();

  if (parameterID == 1) { // E
    s(0) = A*e(0);
    s(1) = I*e(1);
  }
  if (parameterID == 2) // A
    s(0) = E*e(0);
  if (parameterID == 3) // I
    s(1) = E*e(1);

  return s;
}

const Matrix&
ElasticWarpingShearSection2d::getInitialTangentSensitivity(int gradIndex)
{
  ks.Zero();

  return ks;
}
