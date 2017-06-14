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
                                                                        
// $Revision: 1.25 $                                                              
// $Date: 2009-01-29 00:42:03 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicMaterial.cpp,v $                                                                
// Written: MHS 
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for ElasticIsotropicMaterial.
//
// What: "@(#) ElasticIsotropicMaterial.C, revA"

#include <string.h>

#include <ElasticIsotropicMaterial.h>
#include <ElasticIsotropicPlaneStress2D.h>
#include <ElasticIsotropicPlaneStrain2D.h>
#include <ElasticIsotropicAxiSymm.h>
#include <ElasticIsotropicThreeDimensional.h>
#include <ElasticIsotropicPlateFiber.h>
#include <ElasticIsotropicBeamFiber.h>
#include <ElasticIsotropicBeamFiber2d.h>

#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>


void *
OPS_ElasticIsotropicMaterial(void)
{
  NDMaterial *theMaterial = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 3) {
    opserr << "Want: nDMaterial ElasticIsotropic $tag $E $V <$rho>" << endln;
    return 0;	
  }
  
  int iData[1];
  double dData[3];
  dData[2] = 0.0;
  
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial EasticIsotropic \n";
    return 0;
  }
  
  if (numArgs > 3) 
    numData = 3;
  else
    numData = 2;
  
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial EasticIsotropic : " << iData[0] <<"\n";
    return 0;
  }  
  
  theMaterial = new ElasticIsotropicMaterial(iData[0], dData[0], dData[1], dData[2]);
  
  return theMaterial;
}



ElasticIsotropicMaterial::ElasticIsotropicMaterial
(int tag, int classTag, double e, double nu, double r)
  :NDMaterial(tag, classTag), E(e), v(nu), rho(r), parameterID(0)
{

}

ElasticIsotropicMaterial::ElasticIsotropicMaterial
(int tag, double e, double nu, double r)
  :NDMaterial(tag, ND_TAG_ElasticIsotropic), E(e), v(nu), rho(r), parameterID(0)
{

}

ElasticIsotropicMaterial::~ElasticIsotropicMaterial()
{
	
}

double
ElasticIsotropicMaterial::getRho() 
{ 
  return rho ;
}

NDMaterial*
ElasticIsotropicMaterial::getCopy (const char *type)
{
  if (strcmp(type,"PlaneStress2D") == 0 || strcmp(type,"PlaneStress") == 0) {
    ElasticIsotropicPlaneStress2D *theModel;
    theModel = new ElasticIsotropicPlaneStress2D (this->getTag(), E, v, rho);
    return theModel;
  } 

  else if (strcmp(type,"PlaneStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0) {
    ElasticIsotropicPlaneStrain2D *theModel;
    theModel = new ElasticIsotropicPlaneStrain2D (this->getTag(), E, v, rho);
    return theModel;
  }

  else if (strcmp(type,"AxiSymmetric2D") == 0 || strcmp(type,"AxiSymmetric") == 0) {
    ElasticIsotropicAxiSymm *theModel;
    theModel = new ElasticIsotropicAxiSymm(this->getTag(), E, v, rho);
    return theModel;
  }
  
  else if (strcmp(type,"ThreeDimensional") == 0 || strcmp(type,"3D") == 0) {
    ElasticIsotropicThreeDimensional *theModel;
    theModel = new ElasticIsotropicThreeDimensional (this->getTag(), E, v, rho);
    return theModel;
  }

  else if (strcmp(type,"PlateFiber") == 0) {
    ElasticIsotropicPlateFiber *theModel;
    theModel = new ElasticIsotropicPlateFiber(this->getTag(), E, v, rho);
    return theModel;
  }

  else if (strcmp(type,"BeamFiber") == 0) {
    ElasticIsotropicBeamFiber *theModel;
    theModel = new ElasticIsotropicBeamFiber(this->getTag(), E, v, rho);
    return theModel;
  }

  else if (strcmp(type,"BeamFiber2d") == 0) {
    ElasticIsotropicBeamFiber2d *theModel;
    theModel = new ElasticIsotropicBeamFiber2d(this->getTag(), E, v, rho);
    return theModel;
  }

  // Handle other cases
  else
    return NDMaterial::getCopy(type);
}

int
ElasticIsotropicMaterial::setTrialStrain (const Vector &v)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticIsotropicMaterial::setTrialStrain (const Vector &v, const Vector &rate)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticIsotropicMaterial::setTrialStrainIncr (const Vector &v)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticIsotropicMaterial::setTrialStrainIncr (const Vector &v, const Vector &rate)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

const Matrix&
ElasticIsotropicMaterial::getTangent (void)
{
  opserr << "ElasticIsotropicMaterial::getTangent -- subclass responsibility\n";
  exit(-1);

  // Just to make it compile
  Matrix *ret = new Matrix();
  return *ret;
}

const Matrix&
ElasticIsotropicMaterial::getInitialTangent (void)
{
  return this->getTangent();
}

const Vector&
ElasticIsotropicMaterial::getStress (void)
{
  opserr << "ElasticIsotropicMaterial::getStress -- subclass responsibility\n";
  exit(-1);
    
  // Just to make it compile
  Vector *ret = new Vector();
  return *ret;
}

const Vector&
ElasticIsotropicMaterial::getStrain (void)
{
  opserr << "ElasticIsotropicMaterial::getStrain -- subclass responsibility\n";
  exit(-1);

  // Just to make it compile
  Vector *ret = new Vector();
  return *ret;
}

int
ElasticIsotropicMaterial::commitState (void)
{
  opserr << "ElasticIsotropicMaterial::commitState -- subclass responsibility\n";
  exit(-1);
  return -1;
}

int
ElasticIsotropicMaterial::revertToLastCommit (void)
{
  opserr << "ElasticIsotropicMaterial::revertToLastCommit -- subclass responsibility\n";
  exit(-1);
    
  return -1;
}

int
ElasticIsotropicMaterial::revertToStart (void)
{
  opserr << "ElasticIsotropicMaterial::revertToStart -- subclass responsibility\n";
  exit(-1);
  return -1;
}

NDMaterial*
ElasticIsotropicMaterial::getCopy (void)
{
  opserr << "ElasticIsotropicMaterial::getCopy -- subclass responsibility\n";
  exit(-1);
  return 0;
}

const char*
ElasticIsotropicMaterial::getType (void) const
{
  opserr << "ElasticIsotropicMaterial::getType -- subclass responsibility\n";
  exit(-1);	

  return 0;
}

int
ElasticIsotropicMaterial::getOrder (void) const
{
  opserr << "ElasticIsotropicMaterial::getOrder -- subclass responsibility\n";
  exit(-1);
  return -1;
}

int
ElasticIsotropicMaterial::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(4);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = v;
  data(3) = rho;
  
 res += theChannel.sendVector(this->getDbTag(), commitTag, data);
 if (res < 0) {
   opserr << "ElasticIsotropicMaterial::sendSelf -- could not send Vector\n";
   return res;
 }

 return res;
}

int
ElasticIsotropicMaterial::recvSelf (int commitTag, Channel &theChannel, 
				    FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(4);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
   opserr << "ElasticIsotropicMaterial::recvSelf -- could not recv Vector\n";
   return res;
  }
    
  this->setTag((int)data(0));
  E = data(1);
  v = data(2);
  rho = data(3);
  
  return res;
}

void
ElasticIsotropicMaterial::Print (OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "Elastic Isotropic Material Model" << endln;
        s << "\tE:  " << E << endln;
        s << "\tv:  " << v << endln;
        s << "\trho:  " << rho << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ElasticIsotropicMaterial\", ";
        s << "\"E\": " << E << ", ";
        s << "\"nu\": " << v << ", ";
        s << "\"rho\": " << rho << "}";
    }
}

int
ElasticIsotropicMaterial::setParameter(const char **argv, int argc,
				      Parameter &param)
{
  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  else if (strcmp(argv[0],"nu") == 0 || strcmp(argv[0],"v") == 0) {
    param.setValue(v);
    return param.addObject(2, this);
  }
  else if (strcmp(argv[0],"rho") == 0) {
    param.setValue(rho);
    return param.addObject(3, this);
  }

  return -1;
}

int 
ElasticIsotropicMaterial::updateParameter(int parameterID, Information &info)
{ 
  switch(parameterID) {
  case 1:
    E = info.theDouble;
    return 0;
  case 2:
    v = info.theDouble;
    return 0;
  case 3:
    rho = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
ElasticIsotropicMaterial::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}
