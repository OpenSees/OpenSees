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
                                                                        
// $Revision: 1.12 $
// $Date: 2008-08-26 16:46:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticSection2d.cpp,v $

#include <ElasticSection2d.h>
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

Vector ElasticSection2d::s(2);
Matrix ElasticSection2d::ks(2,2);
ID ElasticSection2d::code(2);

void* OPS_ElasticSection2d()
{
    if(OPS_GetNumRemainingInputArgs() < 4) {
	opserr<<"insufficient arguments for ealstic section\n";
	return 0;
    }

    // get tag
    int tag;
    int numData = 1;
    if(OPS_GetIntInput(&numData,&tag) < 0) return 0;

    // get data
    numData = 3;
    double data[3];
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;

    return new ElasticSection2d(tag,data[0],data[1],data[2]);
}

ElasticSection2d::ElasticSection2d(void)
:SectionForceDeformation(0, SEC_TAG_Elastic2d),
 E(0.0), A(0.0), I(0.0), e(2)
{
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;	// P is the first quantity
    code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
  }
}

ElasticSection2d::ElasticSection2d
(int tag, double E_in, double A_in, double I_in)
:SectionForceDeformation(tag, SEC_TAG_Elastic2d),
 E(E_in), A(A_in), I(I_in), e(2)
{
  if (E <= 0.0)  {
    //opserr << "ElasticSection2d::ElasticSection2d -- Input E <= 0.0\n";
  }
  
  if (A <= 0.0)  {
    //opserr << "ElasticSection2d::ElasticSection2d -- Input A <= 0.0\n";
  }
  
  if (I <= 0.0)  {
    //opserr << "ElasticSection2d::ElasticSection2d -- Input I <= 0.0\n";
  }    
  
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;	// P is the first quantity
    code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
  }
}

ElasticSection2d::~ElasticSection2d(void)
{
  return;
}

int 
ElasticSection2d::commitState(void)
{
  return 0;
}

int 
ElasticSection2d::revertToLastCommit(void)
{
  return 0;
}

int 
ElasticSection2d::revertToStart(void)
{
  return 0;
}

int
ElasticSection2d::setTrialSectionDeformation(const Vector &def)
{
    e = def;

    return 0;
}

const Vector &
ElasticSection2d::getSectionDeformation(void)
{
    return e;
}

const Vector &
ElasticSection2d::getStressResultant(void)
{
  s(0) = E*A*e(0);
  s(1) = E*I*e(1);    

  return s;
}

const Matrix &
ElasticSection2d::getSectionTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*I;
  
  return ks;
}

const Matrix &
ElasticSection2d::getInitialTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*I;
  
  return ks;
}

const Matrix &
ElasticSection2d::getSectionFlexibility(void)
{
  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*I);
  
  return ks;
}

const Matrix &
ElasticSection2d::getInitialFlexibility(void)
{
  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*I);
  
  return ks;
}

SectionForceDeformation*
ElasticSection2d::getCopy(void)
{
    // Make a copy of the hinge
    ElasticSection2d *theCopy =
	new ElasticSection2d (this->getTag(), E, A, I);

    theCopy->parameterID = parameterID;

    return theCopy;
}

const ID&
ElasticSection2d::getType(void)
{
    return code;
}

int
ElasticSection2d::getOrder(void) const
{
    return 2;
}

int
ElasticSection2d::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(4);
    
    int dataTag = this->getDbTag();
    
    data(0) = this->getTag();
    data(1) = E;
    data(2) = A;
    data(3) = I;    

    res += theChannel.sendVector(dataTag, commitTag, data);
    if (res<0) {
      opserr << "ElasticSection2d::sendSelf -- failed to send data\n";
      return res;
    }
    
    return res;
}

int
ElasticSection2d::recvSelf(int commitTag, Channel &theChannel,
			   FEM_ObjectBroker &theBroker)
{
	int res = 0;

    static Vector data(4);

    int dataTag = this->getDbTag();

    res += theChannel.recvVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticSection2d::recvSelf -- failed to receive data\n";
      return res;
    }
    
	this->setTag((int)data(0));
    E = data(1);
    A = data(2);
    I = data(3);

    return res;
}
 
void
ElasticSection2d::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
    s << "ElasticSection2d, tag: " << this->getTag() << endln;
    s << "\tE: " << E << endln;
    s << "\tA: " << A << endln;
    s << "\tI: " << I << endln;
  }
  
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": \"" << this->getTag() << "\", ";
	s << "\"type\": \"ElasticSection2d\", ";
    s << "\"E\": " << E << ", ";
    s << "\"A\": " << A << ", ";
    s << "\"Iz\": " << I << "}";
  }
}  

int
ElasticSection2d::setParameter(const char **argv, int argc, Parameter &param)
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

  return -1;
}

int
ElasticSection2d::updateParameter(int paramID, Information &info)
{
  if (paramID == 1)
    E = info.theDouble;
  if (paramID == 2)
    A = info.theDouble;
  if (paramID == 3)
    I = info.theDouble;

  return 0;
}

int
ElasticSection2d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
ElasticSection2d::getStressResultantSensitivity(int gradIndex,
						bool conditional)
{
  s.Zero();

  //opserr << "ElasticSection parameterID: " << parameterID << endln;

  if (parameterID == 1) { // E
    s(0) = A*e(0);
    s(1) = I*e(1);
  }
  if (parameterID == 2) // A
    s(0) = E*e(0);
  if (parameterID == 3) // I
    s(1) = E*e(1);

  //opserr << "ElasticSection2d::dsdh: " << s << endln;

  return s;
}

const Matrix&
ElasticSection2d::getSectionTangentSensitivity(int gradIndex)
{
  ks.Zero();

  if (parameterID == 1) { // E
    ks(0,0) = A;
    ks(1,1) = I;
  }
  if (parameterID == 2) // A
    ks(0,0) = E;
  if (parameterID == 3) // I
    ks(1,1) = E;

  return ks;
}

const Matrix&
ElasticSection2d::getInitialTangentSensitivity(int gradIndex)
{
  return this->getSectionTangentSensitivity(gradIndex);
}

const Matrix&
ElasticSection2d::getSectionFlexibilitySensitivity(int gradIndex)
{
  ks.Zero();

  if (parameterID == 1) { // E
    ks(0,0) = -1.0/(E*E*A);
    ks(1,1) = -1.0/(E*E*I);
  }
  if (parameterID == 2) // A
    ks(0,0) = -1.0/(E*A*A);
  if (parameterID == 3) // I
    ks(1,1) = -1.0/(E*I*I);

  return ks;
}

const Matrix& 
ElasticSection2d::getInitialFlexibilitySensitivity(int gradIndex)
{
  return this->getSectionFlexibilitySensitivity(gradIndex);
}
