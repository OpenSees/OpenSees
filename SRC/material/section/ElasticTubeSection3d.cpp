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
	opserr << "Want: section ElasticTube tag? E? d? tw? G?" << endln;
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
    double d = data[1];
    double tw = data[2];
    double G = data[3];

    return new ElasticTubeSection3d(tag, E, d, tw, G);	
}

Vector ElasticTubeSection3d::s(4);
Matrix ElasticTubeSection3d::ks(4,4);
ID ElasticTubeSection3d::code(4);

ElasticTubeSection3d::ElasticTubeSection3d(void)
:SectionForceDeformation(0, SEC_TAG_ElasticTube3d),
 E(0.0), d(0.0), tw(0.0), G(0.0),
 e(4), parameterID(0)
{
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;	// P is the first quantity
    code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
    code(2) = SECTION_RESPONSE_MY;	// My is the third
    code(3) = SECTION_RESPONSE_T;   // T is the fourth
  }    
}

ElasticTubeSection3d::ElasticTubeSection3d
(int tag, double E_in, double d_in, double tw_in, double G_in)
:SectionForceDeformation(tag, SEC_TAG_ElasticTube3d),
 E(E_in), d(d_in), tw(tw_in), G(G_in),
 e(4), parameterID(0)
{
  if (E <= 0.0)  {
    opserr << "ElasticTubeSection3d::ElasticTubeSection3d -- Input E <= 0.0\n";
  }
	
  if (G <= 0.0)  {
    opserr << "ElasticTubeSection3d::ElasticTubeSection3d -- Input G <= 0.0\n";
  }
	
  if (d <= 0.0)  {
    opserr << "ElasticTubeSection3d::ElasticTubeSection3d -- Input d <= 0.0\n";
  }
    
  if (tw <= 0.0)  {
    opserr << "ElasticTubeSection3d::ElasticTubeSection3d -- Input tw <= 0.0\n";
  }    
  
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;	// P is the first quantity
    code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
    code(2) = SECTION_RESPONSE_MY;	// My is the third
    code(3) = SECTION_RESPONSE_T;	// My is the fourth
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
  double r1 = 0.5*d;
  double r2 = r1-tw;
  double A = 3.14159*(r1*r1 - r2*r2);
  double I = 0.25*3.14159*(r1*r1*r1*r1 - r2*r2*r2*r2);
  double J = 2*I;

  s(0) = E*A*e(0);
  s(1) = E*I*e(1);    
  s(2) = E*I*e(2);
  s(3) = G*J*e(3);

  return s;
}

const Matrix &
ElasticTubeSection3d::getSectionTangent(void)
{
  double r1 = 0.5*d;
  double r2 = r1-tw;
  double A = 3.14159*(r1*r1 - r2*r2);
  double I = 0.25*3.14159*(r1*r1*r1*r1 - r2*r2*r2*r2);
  double J = 2*I;

  ks(0,0) = E*A;
  ks(1,1) = E*I;
  ks(2,2) = E*I;
  ks(3,3) = G*J;
  
  return ks;
}

const Matrix &
ElasticTubeSection3d::getInitialTangent(void)
{
  double r1 = 0.5*d;
  double r2 = r1-tw;
  double A = 3.14159*(r1*r1 - r2*r2);
  double I = 0.25*3.14159*(r1*r1*r1*r1 - r2*r2*r2*r2);
  double J = 2*I;
  
  ks(0,0) = E*A;
  ks(1,1) = E*I;
  ks(2,2) = E*I;
  ks(3,3) = G*J;
  
  return ks;
}

const Matrix &
ElasticTubeSection3d::getSectionFlexibility (void)
{
  double r1 = 0.5*d;
  double r2 = r1-tw;
  double A = 3.14159*(r1*r1 - r2*r2);
  double I = 0.25*3.14159*(r1*r1*r1*r1 - r2*r2*r2*r2);
  double J = 2*I;

  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*I);
  ks(2,2) = 1.0/(E*I);
  ks(3,3) = 1.0/(G*J);
  
  return ks;
}

const Matrix &
ElasticTubeSection3d::getInitialFlexibility(void)
{
  double r1 = 0.5*d;
  double r2 = r1-tw;
  double A = 3.14159*(r1*r1 - r2*r2);
  double I = 0.25*3.14159*(r1*r1*r1*r1 - r2*r2*r2*r2);
  double J = 2*I;

  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*I);
  ks(2,2) = 1.0/(E*I);
  ks(3,3) = 1.0/(G*J);
  
  return ks;
}

SectionForceDeformation*
ElasticTubeSection3d::getCopy(void)
{
  ElasticTubeSection3d *theCopy =
    new ElasticTubeSection3d (this->getTag(), E, d, tw, G);
  
  theCopy->parameterID = parameterID;
  
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
  return 4;
}

int
ElasticTubeSection3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(5);
  
  int dataTag = this->getDbTag();
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = d;
  data(3) = tw;    
  data(4) = G;    

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

  static Vector data(5);
  
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
  G =  data(4);
  
  return res;
}
 
void
ElasticTubeSection3d::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_CURRENTSTATE) {
        s << "ElasticTubeSection3d, tag: " << this->getTag() << endln;
        s << "\tE: " << E << endln;
        s << "\td: " << d << endln;
        s << "\ttw: " << tw << endln;
        s << "\tG: " << G << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ElasticTubeSection3d\", ";
        s << "\"E\": " << E << ", ";
        s << "\"G\": " << G << ", ";
        s << "\"diameter\": " << d << ", ";
        s << "\"thickness\": " << tw << "}";
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
  
  if (strcmp(argv[0],"G") == 0) {
    param.setValue(G);    
    return param.addObject(4, this);
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
    G = info.theDouble;

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
  static Vector dsdh(4);
  dsdh.Zero();

  return dsdh;
}

const Matrix&
ElasticTubeSection3d::getInitialTangentSensitivity(int gradIndex)
{
  ks.Zero();

  return ks;
}
