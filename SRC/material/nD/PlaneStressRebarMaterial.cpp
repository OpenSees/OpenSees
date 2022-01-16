/* ****************************************************************** **
**    OpenSees System for Earthquake Engineering Simulation    **
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

// written: fmk derived from code in PlateRebarMaterial from
// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)

#include <PlaneStressRebarMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>   //Antonios Vytiniotis used for the recorder
#include <math.h>
#include <elementAPI.h>

//static vector and matrices
Vector  PlaneStressRebarMaterial::stress(3) ;
Matrix  PlaneStressRebarMaterial::tangent(3,3) ;

//null constructor
PlaneStressRebarMaterial::PlaneStressRebarMaterial( ) : 
NDMaterial(0, ND_TAG_PlaneStressRebarMaterial ), 
strain(5) 
{ }

void * OPS_ADD_RUNTIME_VPV(OPS_PlaneStressRebarMaterial)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 3) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial PlaneStressRebarMaterial tag? matTag? angle?" << endln;
	return 0;
    }

    int tag[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata, tag)<0) {
	opserr << "WARNING invalid nDMaterial PlaneStressRebarMaterial tag or matTag" << endln;
	return 0;
    }

    UniaxialMaterial *theMaterial = OPS_getUniaxialMaterial(tag[1]);
    if (theMaterial == 0) {
      opserr << "WARNING uniaxialmaterial does not exist\n";
      opserr << "UniaxialMaterial: " << tag[1];
      opserr << "\nPlaneStressRebarMaterial nDMaterial: " << tag[0] << endln;
      return 0;
    }
    
    double angle;
    numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &angle)<0) {
      opserr << "WARNING invalid angle" << endln;
      return 0;
    }
    
    NDMaterial* mat = new PlaneStressRebarMaterial( tag[0], *theMaterial, angle);

    if (mat == 0) {
	opserr << "WARNING: failed to create PlaneStressRebarMaterial material\n";
	return 0;
    }
    return mat;
}


//full constructor
PlaneStressRebarMaterial::PlaneStressRebarMaterial(int tag,
						   UniaxialMaterial &uniMat,
						   double ang)
:NDMaterial( tag, ND_TAG_PlaneStressRebarMaterial),
 strain(3),angle(ang)
{
  theMat = uniMat.getCopy() ;
  double rang = ang * 4.0 * asin(1.0)/360.0;
  c = cos(rang);
  s = sin(rang);
}


//destructor
PlaneStressRebarMaterial::~PlaneStressRebarMaterial( ) 
{ 
  if (theMat != 0) delete theMat ;
} 


//make a clone of this material
NDMaterial*
PlaneStressRebarMaterial::getCopy( ) 
{
  PlaneStressRebarMaterial *clone = 0;   //new instance of this class

  clone = new PlaneStressRebarMaterial( this->getTag(), 
					*theMat,
					angle ) ; //make the copy
  return clone ;
}


//make a clone of this material
NDMaterial* 
PlaneStressRebarMaterial::getCopy( const char *type ) 
{
  if (strcmp(type,this->getType()) == 0)
    return this->getCopy( ) ;
  else
    return 0;
}


//send back order of strain in vector form
int 
PlaneStressRebarMaterial::getOrder( ) const
{
  return 3;
}


const char*
PlaneStressRebarMaterial::getType( ) const 
{
  return "PlaneStress2D" ; 
}



//swap history variables
int 
PlaneStressRebarMaterial::commitState( ) 
{
  return theMat->commitState( ) ;
}



//revert to last saved state
int 
PlaneStressRebarMaterial::revertToLastCommit( )
{
  return theMat->revertToLastCommit( ) ;
}


//revert to start
int
PlaneStressRebarMaterial::revertToStart( )
{
  strain.Zero();
  return theMat->revertToStart( ) ;
}


//mass per unit volume
double
PlaneStressRebarMaterial::getRho( )
{
  return theMat->getRho( ) ;
}


//receive the strain
int 
PlaneStressRebarMaterial::setTrialStrain( const Vector &strainFromElement )
{
  strain = strainFromElement;

  if (angle == 0)
    return theMat->setTrialStrain(strain(0));
  else if (angle == 90)
    return theMat->setTrialStrain(strain(1));
  
  return theMat->setTrialStrain(strain(0) * c * c
				+ strain(1) * s * s
				+ strain(2) * c * s,
				0) ;
}


//send back the strain
const Vector& 
PlaneStressRebarMaterial::getStrain( )
{
  return strain ;
}


//send back the stress 
const Vector&  
PlaneStressRebarMaterial::getStress( )
{
  double sig = theMat->getStress();

  stress.Zero();
  if (angle == 0) 
    stress(0) = sig; 
  else if (angle == 90)
    stress(1)= sig;
  else {
    stress(0) = sig * c * c;
    stress(1) = sig * s * s;
    stress(2) = sig * c * s;
  }

  return stress ;
}


//send back the tangent 
const Matrix&  
PlaneStressRebarMaterial::getTangent( )
{
  double tan = theMat->getTangent( ) ;

  tangent.Zero();
  if (angle == 0) 
    tangent(0,0) = tan;
  else if (angle == 90)
    tangent(1,1) = tan;
  else {
    tangent(0,0) = tan * c * c * c * c ;
    tangent(0,2) = tan * c * c * c * s ;
    tangent(0,1) = tan * c * c * s * s ;
    tangent(2,0) = tangent(0,1) ;
    tangent(2,2) = tangent(0,2) ;
    tangent(2,1) = tan * c * s * s * s ;
    tangent(1,0) = tangent(0,2) ;
    tangent(1,2) = tangent(1,2) ;
    tangent(1,1) = tan * s * s * s * s ;
  }

  return tangent ;
}

const Matrix&  
PlaneStressRebarMaterial::getInitialTangent
( )
{
  double tan = theMat->getInitialTangent( ) ;
  tangent.Zero();
  if (angle == 0) 
    tangent(0,0) = tan;
  else if (angle == 90)
    tangent(1,1) = tan;
  else {
    tangent(0,0) = tan * c * c * c * c ;
    tangent(0,1) = tan * c * c * c * s ;
    tangent(0,2) = tan * c * c * s * s ;
    tangent(1,0) = tangent(0,1) ;
    tangent(1,1) = tangent(0,2) ;
    tangent(1,2) = tan * c * s * s * s ;
    tangent(2,0) = tangent(0,2) ;
    tangent(2,1) = tangent(1,2) ;
    tangent(2,2) = tan * s * s * s * s ;
  }

  return tangent ;
}


//print out data
void  
PlaneStressRebarMaterial::Print( OPS_Stream &s, int flag )
{
  s << "PlaneStressPlateRebar Material tag: " << this->getTag() << endln ; 
  s << "using uniaxialmaterials : " << endln ;

  if (theMat != 0)
    theMat->Print( s, flag ) ;

  return ;
}


int 
PlaneStressRebarMaterial::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  int dataTag = this->getDbTag();

  int matDbTag;
  
  static ID idData(3);
  idData(0) = dataTag;
  idData(1) = theMat->getClassTag();
  matDbTag = theMat->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMat->setDbTag(matDbTag);
  }
  idData(2) = matDbTag;

  res = theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlaneStressRebarMaterial::sendSelf() - failed to send data" << endln;
    return res;
  }

  static Vector vecData(1);
  vecData(0) = angle;

  res = theChannel.sendVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlaneStressRebarMaterial::sendSelf() - failed to send data" << endln;
    return res;
  }

  // now send the materials data
  res += theMat->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlaneStressRebarMaterial::sendSelf() - failed to send material1" << endln;

  return res;
}

int 
PlaneStressRebarMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // recv an id containing the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlaneStressRebarMaterial::sendSelf() - failed to receive id data" << endln;
    return res;
  }

  this->setTag(idData(0));
  int matClassTag = idData(1);
  if (theMat->getClassTag() != matClassTag) {
    if (theMat != 0) delete theMat;
    theMat = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMat == 0) {
      opserr << "PlaneStressRebarMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMat->setDbTag(idData(2));

  static Vector vecData(1);
  res = theChannel.recvVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlaneStressRebarMaterial::sendSelf() - failed to receive vector data" << endln;
    return res;
  }
  angle = vecData(0);
  double rang = angle * 0.0174532925;
  c = cos(rang);
  s = sin(rang);

  // now receive the materials data
  res = theMat->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "PlaneStressRebarMaterial::sendSelf() - failed to receive material1" << endln;
  
  return res;
}
 
