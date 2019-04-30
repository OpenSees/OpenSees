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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticOrthotropicMaterial.cpp,v $                                                                
// Written: MHS 
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for ElasticOrthotropicMaterial.
//
// What: "@(#) ElasticOrthotropicMaterial.C, revA"

#include <string.h>

#include <ElasticOrthotropicMaterial.h>
//#include <ElasticOrthotropicPlaneStress2D.h>
//#include <ElasticOrthotropicPlaneStrain2D.h>
//#include <ElasticOrthotropicAxiSymm.h>
#include <ElasticOrthotropicThreeDimensional.h>
//#include <ElasticOrthotropicPlateFiber.h>
//#include <ElasticOrthotropicBeamFiber.h>
//#include <ElasticOrthotropicBeamFiber2d.h>

#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>


void *
OPS_ElasticOrthotropicMaterial(void)
{
  NDMaterial *theMaterial = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 10) {
    opserr << "Want: nDMaterial ElasticOrthotropic $tag $Ex $Ey $Ez $vxy $vyz $vzx $Gxy $Gyz $Gzx <$rho>" << endln;
    return 0;	
  }
  
  int iData[1];
  double dData[10];
  dData[9] = 0.0;
  
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial ElasticOrthotropic \n";
    return 0;
  }
  
  if (numArgs > 10) 
    numData = 10;
  else
    numData = 9;
  
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial EasticIsotropic : " << iData[0] <<"\n";
    return 0;
  }  
  
  theMaterial = new ElasticOrthotropicMaterial(iData[0], 
	dData[0], dData[1], dData[2],
	dData[3], dData[4], dData[5],
	dData[6], dData[7], dData[8], dData[9]);
  
  return theMaterial;
}



ElasticOrthotropicMaterial::ElasticOrthotropicMaterial
(int tag, int classTag, double ex, double ey, double ez,
double nuxy, double nuyz, double nuzx,
double gxy, double gyz, double gzx, double r)
  :NDMaterial(tag, classTag), 
Ex(ex), Ey(ey), Ez(ez),
vxy(nuxy), vyz(nuyz), vzx(nuzx),
Gxy(gxy), Gyz(gyz), Gzx(gzx), 
rho(r),  parameterID(0)
{

}

ElasticOrthotropicMaterial::ElasticOrthotropicMaterial
(int tag, double ex, double ey, double ez,
double nuxy, double nuyz, double nuzx,
double gxy, double gyz, double gzx, double r)
  :NDMaterial(tag, ND_TAG_ElasticOrthotropic), 
Ex(ex), Ey(ey), Ez(ez),
vxy(nuxy), vyz(nuyz), vzx(nuzx),
Gxy(gxy), Gyz(gyz), Gzx(gzx), 
rho(r),  parameterID(0)
{

}

ElasticOrthotropicMaterial::~ElasticOrthotropicMaterial()
{
	
}

double
ElasticOrthotropicMaterial::getRho() 
{ 
  return rho ;
}

NDMaterial*
ElasticOrthotropicMaterial::getCopy (const char *type)
{
  if (strcmp(type,"ThreeDimensional") == 0 || strcmp(type,"3D") == 0) {
    ElasticOrthotropicThreeDimensional *theModel;
      theModel = new ElasticOrthotropicThreeDimensional (this->getTag(), Ex, Ey, Ez,
                                                         vxy, vyz, vzx, Gxy, Gyz, Gzx, rho);

    return theModel;
  }

  // Handle other cases
  else
    return NDMaterial::getCopy(type);
}

int
ElasticOrthotropicMaterial::setTrialStrain (const Vector &v)
{
    opserr << "ElasticOrthotropicMaterial::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticOrthotropicMaterial::setTrialStrain (const Vector &v, const Vector &rate)
{
    opserr << "ElasticOrthotropicMaterial::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticOrthotropicMaterial::setTrialStrainIncr (const Vector &v)
{
    opserr << "ElasticOrthotropicMaterial::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticOrthotropicMaterial::setTrialStrainIncr (const Vector &v, const Vector &rate)
{
    opserr << "ElasticOrthotropicMaterial::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

const Matrix&
ElasticOrthotropicMaterial::getTangent (void)
{
  opserr << "ElasticOrthotropicMaterial::getTangent -- subclass responsibility\n";
  exit(-1);

  // Just to make it compile
  Matrix *ret = new Matrix();
  return *ret;
}

const Matrix&
ElasticOrthotropicMaterial::getInitialTangent (void)
{
  return this->getTangent();
}

const Vector&
ElasticOrthotropicMaterial::getStress (void)
{
  opserr << "ElasticOrthotropicMaterial::getStress -- subclass responsibility\n";
  exit(-1);
    
  // Just to make it compile
  Vector *ret = new Vector();
  return *ret;
}

const Vector&
ElasticOrthotropicMaterial::getStrain (void)
{
  opserr << "ElasticOrthotropicMaterial::getStrain -- subclass responsibility\n";
  exit(-1);

  // Just to make it compile
  Vector *ret = new Vector();
  return *ret;
}

int
ElasticOrthotropicMaterial::commitState (void)
{
  opserr << "ElasticOrthotropicMaterial::commitState -- subclass responsibility\n";
  exit(-1);
  return -1;
}

int
ElasticOrthotropicMaterial::revertToLastCommit (void)
{
  opserr << "ElasticOrthotropicMaterial::revertToLastCommit -- subclass responsibility\n";
  exit(-1);
    
  return -1;
}

int
ElasticOrthotropicMaterial::revertToStart (void)
{
  opserr << "ElasticOrthotropicMaterial::revertToStart -- subclass responsibility\n";
  exit(-1);
  return -1;
}

NDMaterial*
ElasticOrthotropicMaterial::getCopy (void)
{
  opserr << "ElasticOrthotropicMaterial::getCopy -- subclass responsibility\n";
  exit(-1);
  return 0;
}

const char*
ElasticOrthotropicMaterial::getType (void) const
{
  opserr << "ElasticOrthotropicMaterial::getType -- subclass responsibility\n";
  exit(-1);	

  return 0;
}

int
ElasticOrthotropicMaterial::getOrder (void) const
{
  opserr << "ElasticOrthotropicMaterial::getOrder -- subclass responsibility\n";
  exit(-1);
  return -1;
}

int
ElasticOrthotropicMaterial::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(11);
  
  data(0) = this->getTag();
  data(1) = Ex;
  data(2) = Ey;
  data(3) = Ez;
  data(4) = vxy;
  data(5) = vyz;
  data(6) = vzx;
  data(7) = Gxy;
  data(8) = Gyz;
  data(9) = Gzx;
  data(10) = rho;
  
 res += theChannel.sendVector(this->getDbTag(), commitTag, data);
 if (res < 0) {
   opserr << "ElasticOrthotropicMaterial::sendSelf -- could not send Vector\n";
   return res;
 }

 return res;
}

int
ElasticOrthotropicMaterial::recvSelf (int commitTag, Channel &theChannel, 
				    FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(11);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
   opserr << "ElasticOrthotropicMaterial::recvSelf -- could not recv Vector\n";
   return res;
  }
    
  this->setTag((int)data(0));
  Ex = data(1);
  Ey = data(2);
  Ez = data(3);
  vxy = data(4);
  vyz = data(5);
  vzx = data(6);
  Gxy = data(7);
  Gyz = data(8);
  Gzx = data(9);
  rho = data(10);
  
  return res;
}

void
ElasticOrthotropicMaterial::Print (OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "Elastic Isotropic Material Model" << endln;
        s << "\tEx:  " << Ex << endln;
        s << "\tEy:  " << Ey << endln;
        s << "\tEz:  " << Ez << endln;
        s << "\tvxy:  " << vxy << endln;
        s << "\tvyz:  " << vyz << endln;
        s << "\tvzx:  " << vzx << endln;
        s << "\tGxy:  " << Gxy << endln;
        s << "\tGyz:  " << Gyz << endln;
        s << "\tGzx:  " << Gzx << endln;
        s << "\trho:  " << rho << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ElasticOrthotropicMaterial\", ";
        s << "\"Ex\": " << Ex << ", ";
        s << "\"Ey\": " << Ey << ", ";
        s << "\"Ez\": " << Ez << ", ";
        s << "\"nuxy\": " << vxy << ", ";
        s << "\"nuyz\": " << vyz << ", ";
        s << "\"nuzx\": " << vzx << ", ";
        s << "\"Gxy\": " << Gxy << ", ";
        s << "\"Gyz\": " << Gyz << ", ";
        s << "\"Gzx\": " << Gzx << ", ";
        s << "\"rho\": " << rho << "}";
    }
}

int
ElasticOrthotropicMaterial::setParameter(const char **argv, int argc,
				      Parameter &param)
{
  if (strcmp(argv[0],"Ex") == 0) {
    param.setValue(Ex);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"Ey") == 0) {
    param.setValue(Ey);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"Ez") == 0) {
    param.setValue(Ez);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"vxy") == 0 || strcmp(argv[0],"vyx") == 0) {
    param.setValue(vxy);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"vyz") == 0 || strcmp(argv[0],"vzy") == 0) {
    param.setValue(vyz);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"vzx") == 0 || strcmp(argv[0],"vxz") == 0) {
    param.setValue(vzx);
    return param.addObject(6, this);
  }
  if (strcmp(argv[0],"Gxy") == 0 || strcmp(argv[0],"Gyx") == 0) {
    param.setValue(Gxy);
    return param.addObject(7, this);
  }
  if (strcmp(argv[0],"Gyz") == 0 || strcmp(argv[0],"Gzy") == 0) {
    param.setValue(Gyz);
    return param.addObject(8, this);
  }
  if (strcmp(argv[0],"Gzx") == 0 || strcmp(argv[0],"Gxz") == 0) {
    param.setValue(Gzx);
    return param.addObject(9, this);
  }
  if (strcmp(argv[0],"rho") == 0) {
    param.setValue(rho);
    return param.addObject(10, this);
  }

  return -1;
}

int 
ElasticOrthotropicMaterial::updateParameter(int parameterID, Information &info)
{ 
  switch(parameterID) {
  case 1:
    Ex = info.theDouble;
    return 0;
  case 2:
    Ey = info.theDouble;
    return 0;
  case 3:
    Ez = info.theDouble;
    return 0;
  case 4:
    vxy = info.theDouble;
    return 0;
  case 5:
    vyz = info.theDouble;
    return 0;
  case 6:
    vzx = info.theDouble;
    return 0;
  case 7:
    Gxy = info.theDouble;
    return 0;
  case 8:
    Gyz = info.theDouble;
    return 0;
  case 9:
    Gzx = info.theDouble;
    return 0;
  case 10:
    rho = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
ElasticOrthotropicMaterial::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}
