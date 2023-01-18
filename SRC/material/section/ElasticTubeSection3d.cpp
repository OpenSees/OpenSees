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

#include <ElasticTubeSection3d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <stdlib.h>
#include <Information.h>
#include <Parameter.h>

#include <classTags.h>
#include <elementAPI.h>

void* OPS_ElasticTubeSection3d()
{
    if (OPS_GetNumRemainingInputArgs() < 5) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section ElasticTube tag? E? nu? d? tw? <shape?>" << endln;
	return 0;
    }
	
    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid section ElasticTube tag" << endln;
	return 0;
    }

    numdata = 4;
    double data[4];
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid double inputs" << endln;
	opserr << "ElasticTube section: " << tag << endln;	    
	return 0;
    }

    double E = data[0];
    double nu = data[1];
    double d = data[2];
    double tw = data[3];

    if (OPS_GetNumRemainingInputArgs() > 0) {
      double shape;
      numdata = 1;
      if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid shape factor input" << endln;
	opserr << "ElasticTube section: " << tag << endln;	    
	return 0;
      }
      return new ElasticTubeSection3d(tag, E, nu, d, tw, shape);
    }
    else
      return new ElasticTubeSection3d(tag, E, nu, d, tw);
}

Vector ElasticTubeSection3d::s(6);
Matrix ElasticTubeSection3d::ks(6,6);
ID ElasticTubeSection3d::code(6);

ElasticTubeSection3d::ElasticTubeSection3d(void)
:SectionForceDeformation(0, SEC_TAG_ElasticTube3d),
 E(0.0), nu(0.0), d(0.0), tw(0.0), shape(0.0),
 e(6), parameterID(0)
{
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;	// P is the first quantity
    code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
    code(2) = SECTION_RESPONSE_MY;	// My is the third
    code(3) = SECTION_RESPONSE_T;   // T is the fourth
    code(4) = SECTION_RESPONSE_VY;
    code(5) = SECTION_RESPONSE_VZ;    
  }    
}

ElasticTubeSection3d::ElasticTubeSection3d
(int tag, double E_in, double nu_in, double d_in, double tw_in, double shape_in)
:SectionForceDeformation(tag, SEC_TAG_ElasticTube3d),
 E(E_in), nu(nu_in), d(d_in), tw(tw_in), shape(shape_in),
 e(6), parameterID(0)
{
  if (E <= 0.0)  {
    opserr << "ElasticTubeSection3d::ElasticTubeSection3d -- Input E <= 0.0" << endln;
  }
	
  if (nu <= 0.0)  {
    opserr << "ElasticTubeSection3d::ElasticTubeSection3d -- Input nu <= 0.0" << endln;
  }

  if (nu > 0.5)  {
    opserr << "ElasticTubeSection3d::ElasticTubeSection3d -- Input nu > 0.5 "
	   << "You may have entered G instead of nu" << endln;
  }  
	
  if (d <= 0.0)  {
    opserr << "ElasticTubeSection3d::ElasticTubeSection3d -- Input d <= 0.0" << endln;
  }
    
  if (tw <= 0.0)  {
    opserr << "ElasticTubeSection3d::ElasticTubeSection3d -- Input tw <= 0.0" << endln;
  }    
  
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;	// P is the first quantity
    code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
    code(2) = SECTION_RESPONSE_MY;	// My is the third
    code(3) = SECTION_RESPONSE_T;	// My is the fourth
    code(4) = SECTION_RESPONSE_VY;
    code(5) = SECTION_RESPONSE_VZ;        
  }
}

ElasticTubeSection3d::~ElasticTubeSection3d(void)
{

}

int 
ElasticTubeSection3d::commitState(void)
{
  return 0;
}

int 
ElasticTubeSection3d::revertToLastCommit(void)
{
  return 0;
}

int 
ElasticTubeSection3d::revertToStart(void)
{
  return 0;
}

int
ElasticTubeSection3d::setTrialSectionDeformation (const Vector &def)
{
  e = def;

  return 0;
}

const Vector &
ElasticTubeSection3d::getSectionDeformation (void)
{
  return e;
}

const Vector &
ElasticTubeSection3d::getStressResultant (void)
{
  static const double pi = 3.1415926535897932;
  
  double r1 = 0.5*d;
  double r2 = r1-tw;
  double A = pi*(r1*r1 - r2*r2);
  double I = 0.25*pi*(r1*r1*r1*r1 - r2*r2*r2*r2);
  double J = 2*I;

  double G = 0.5*E/(1+nu);
  
  s(0) = E*A*e(0);
  s(1) = E*I*e(1);    
  s(2) = E*I*e(2);
  s(3) = G*J*e(3);
  s(4) = shape*G*A*e(4);
  s(5) = shape*G*A*e(5);  

  return s;
}

const Matrix &
ElasticTubeSection3d::getSectionTangent(void)
{
  static const double pi = 3.1415926535897932;
  
  double r1 = 0.5*d;
  double r2 = r1-tw;
  double A = pi*(r1*r1 - r2*r2);
  double I = 0.25*pi*(r1*r1*r1*r1 - r2*r2*r2*r2);
  double J = 2*I;

  double G = 0.5*E/(1+nu);
  
  ks(0,0) = E*A;
  ks(1,1) = E*I;
  ks(2,2) = ks(1,1);
  ks(3,3) = G*J;
  ks(4,4) = shape*G*A;
  ks(5,5) = ks(4,4);
    
  return ks;
}

const Matrix &
ElasticTubeSection3d::getInitialTangent(void)
{
  return this->getSectionTangent();
}

const Matrix &
ElasticTubeSection3d::getSectionFlexibility (void)
{
  static const double pi = 3.1415926535897932;
  
  double r1 = 0.5*d;
  double r2 = r1-tw;
  double A = pi*(r1*r1 - r2*r2);
  double I = 0.25*pi*(r1*r1*r1*r1 - r2*r2*r2*r2);
  double J = 2*I;

  double G = 0.5*E/(1+nu);
  
  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*I);
  ks(2,2) = ks(1,1);
  ks(3,3) = 1.0/(G*J);
  ks(4,4) = 1.0/(shape*G*A);
  ks(5,5) = ks(4,4);
  
  return ks;
}

const Matrix &
ElasticTubeSection3d::getInitialFlexibility(void)
{
  return this->getSectionFlexibility();
}

SectionForceDeformation*
ElasticTubeSection3d::getCopy(void)
{
  ElasticTubeSection3d *theCopy =
    new ElasticTubeSection3d (this->getTag(), E, nu, d, tw, shape);
  
  theCopy->parameterID = parameterID;
  theCopy->e = e;
  
  return theCopy;
}

const ID&
ElasticTubeSection3d::getType(void)
{
  return code;
}

int
ElasticTubeSection3d::getOrder(void) const
{
  return 6;
}

int
ElasticTubeSection3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(6);
  
  int dataTag = this->getDbTag();
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = d;
  data(3) = tw;    
  data(4) = nu;
  data(5) = shape;    

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res<0) {
    opserr << "ElasticTubeSection3d::sendSelf -- failed to send data\n";
    return res;
  }
  
  return res;
}

int
ElasticTubeSection3d::recvSelf(int commitTag, Channel &theChannel,
				FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(6);
  
  int dataTag = this->getDbTag();
  
  res += theChannel.recvVector(dataTag, commitTag, data);
  if(res < 0) {
    opserr << "ElasticTubeSection3d::recvSelf -- failed to receive data\n";
    return res;
  }
  
  this->setTag((int)data(0));
  E =  data(1);
  d =  data(2);
  tw = data(3);
  nu =  data(4);
  shape =  data(5);  
  
  return res;
}
 
void
ElasticTubeSection3d::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_CURRENTSTATE) {
        s << "ElasticTubeSection3d, tag: " << this->getTag() << endln;
        s << "\tE: " << E << endln;
        s << "\tnu: " << nu << endln;
        s << "\td: " << d << endln;
        s << "\ttw: " << tw << endln;
        s << "\tshape: " << shape << endln;	
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ElasticTubeSection3d\", ";
        s << "\"E\": " << E << ", ";
        s << "\"nu\": " << nu << ", ";
        s << "\"diameter\": " << d << ", ";
        s << "\"thickness\": " << tw << ", ";
	s << "\"shape\": " << shape << "}";
    }
}

int
ElasticTubeSection3d::setParameter(const char **argv, int argc,
				    Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }

  if (strcmp(argv[0],"tw") == 0 || strcmp(argv[0],"t") == 0) {
    param.setValue(tw);    
    return param.addObject(2, this);
  }

  if (strcmp(argv[0],"d") == 0 || strcmp(argv[0],"D") == 0) {
    param.setValue(d);    
    return param.addObject(3, this);
  }
  
  if (strcmp(argv[0],"nu") == 0) {
    param.setValue(nu);    
    return param.addObject(4, this);
  }

  if (strcmp(argv[0],"shape") == 0) {
    param.setValue(shape);    
    return param.addObject(5, this);
  }  
    
  return -1;
}

int
ElasticTubeSection3d::updateParameter(int paramID, Information &info)
{
  if (paramID == 1)
    E = info.theDouble;
  if (paramID == 2)
    tw = info.theDouble;
  if (paramID == 3)
    d = info.theDouble;
  if (paramID == 4)
    nu = info.theDouble;
  if (paramID == 5)
    shape = info.theDouble;  

  return 0;
}

int
ElasticTubeSection3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
ElasticTubeSection3d::getStressResultantSensitivity(int gradIndex,
						       bool conditional)
{
  static Vector dsdh(6);
  dsdh.Zero();

  static const double pi = 3.1415926535897932;
  
  double r1 = 0.5*d;
  double r2 = r1-tw;
  double A = pi*(r1*r1 - r2*r2);
  double I = 0.25*pi*(r1*r1*r1*r1 - r2*r2*r2*r2);
  double J = 2*I;

  double G = 0.5*E/(1+nu);
  
  if (parameterID == 1) { // E
    dsdh(0) = A*e(0);
    dsdh(1) = I*e(1);
    dsdh(2) = I*e(2);
    double dGdh = 0.5/(1+nu);
    dsdh(3) = dGdh*J*e(3);
    dsdh(4) = dGdh*shape*A*e(4);
    dsdh(5) = dGdh*shape*A*e(5);    
  }

  if (parameterID == 2) { // tw
    double dr1dh =  0.0;
    double dr2dh = -1.0;
    
    double dAdh = pi*(2*r1*dr1dh - 2*r2*dr2dh);
    dsdh(0) = E*dAdh*e(0);

    double dIdh = 0.25*pi*(4*r1*r1*r1*dr1dh - 4*r2*r2*r2*dr2dh);
    dsdh(1) = E*dIdh*e(1);
    dsdh(2) = E*dIdh*e(2);

    double dJdh = 2*dIdh;
    dsdh(3) = G*dJdh*e(3);

    dsdh(4) = G*shape*dAdh*e(4);
    dsdh(5) = G*shape*dAdh*e(5);    
  }
  
  if (parameterID == 3) { // D
    double dr1dh = 0.5;
    double dr2dh = 0.5;
    
    double dAdh = pi*(2*r1*dr1dh - 2*r2*dr2dh);
    dsdh(0) = E*dAdh*e(0);

    double dIdh = 0.25*pi*(4*r1*r1*r1*dr1dh - 4*r2*r2*r2*dr2dh);
    dsdh(1) = E*dIdh*e(1);
    dsdh(2) = E*dIdh*e(2);

    double dJdh = 2*dIdh;
    dsdh(3) = G*dJdh*e(3);

    dsdh(4) = G*shape*dAdh*e(4);
    dsdh(5) = G*shape*dAdh*e(5);        
  }

  if (parameterID == 4) { // nu
    double dGdh = -0.5*E/(1 + 2*nu + nu*nu);    
    dsdh(3) = dGdh*J*e(3);
    dsdh(4) = dGdh*shape*A*e(4);
    dsdh(5) = dGdh*shape*A*e(5);    
  }

  if (parameterID == 5) { // shape
    dsdh(4) = G*A*e(4);
    dsdh(5) = G*A*e(5);    
  }

  return dsdh;
}

const Matrix&
ElasticTubeSection3d::getInitialTangentSensitivity(int gradIndex)
{
  ks.Zero();

  return ks;
}
