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
                                                                        
// $Revision: 1.0 $
// $Date: 2012-05-26 22:03:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateRebarMaterial.cpp,v $

//
// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Generic Plate Rebar Material
//
/* Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
concrete high-rise building induced by extreme earthquakes, 
Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723*/

#include <PlateRebarMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>   //Antonios Vytiniotis used for the recorder
#include <math.h>

//static vector and matrices
Vector  PlateRebarMaterial::stress(5) ;
Matrix  PlateRebarMaterial::tangent(5,5) ;

//null constructor
PlateRebarMaterial::PlateRebarMaterial( ) : 
NDMaterial(0, ND_TAG_PlateRebarMaterial ), 
strain(5) 
{ }


//full constructor
PlateRebarMaterial::PlateRebarMaterial(int tag,
                                       UniaxialMaterial &uniMat,
                                       double ang) :
NDMaterial( tag, ND_TAG_PlateRebarMaterial ),
strain(5),angle(ang)
{
  theMat = uniMat.getCopy() ;
  double rang = ang * 0.0174532925;
  c = cos(rang);
  s = sin(rang);
}


//destructor
PlateRebarMaterial::~PlateRebarMaterial( ) 
{ 
  if (theMat != 0) delete theMat ;
} 


//make a clone of this material
NDMaterial*
PlateRebarMaterial::getCopy( ) 
{
  PlateRebarMaterial *clone ;   //new instance of this class

  clone = new PlateRebarMaterial( this->getTag(), 
                                  *theMat,
                                  angle ) ; //make the copy

  return clone ;
}


//make a clone of this material
NDMaterial* 
PlateRebarMaterial::getCopy( const char *type ) 
{
  return this->getCopy( ) ;
}


//send back order of strain in vector form
int 
PlateRebarMaterial::getOrder( ) const
{
  return 5 ;
}


const char*
PlateRebarMaterial::getType( ) const 
{
  return "PlateRebarMaterial" ; 
}



//swap history variables
int 
PlateRebarMaterial::commitState( ) 
{
  return theMat->commitState( ) ;
}



//revert to last saved state
int 
PlateRebarMaterial::revertToLastCommit( )
{
  return theMat->revertToLastCommit( ) ;
}


//revert to start
int
PlateRebarMaterial::revertToStart( )
{
  strain.Zero();
  return theMat->revertToStart( ) ;
}


//mass per unit volume
double
PlateRebarMaterial::getRho( )
{
  return theMat->getRho( ) ;
}


//receive the strain
int 
PlateRebarMaterial::setTrialStrain( const Vector &strainFromElement )
{
  strain(0) = strainFromElement(0) ;
  strain(1) = strainFromElement(1) ;
  strain(2) = strainFromElement(2) ;
  strain(3) = strainFromElement(3) ;
  strain(4) = strainFromElement(4) ;

  return theMat->setTrialStrain(   strain(0) * c * c
                                 + strain(1) * s * s
                                 + strain(2) * c * s,
                                 0) ;
}


//send back the strain
const Vector& 
PlateRebarMaterial::getStrain( )
{
  return strain ;
}


//send back the stress 
const Vector&  
PlateRebarMaterial::getStress( )
{
  double sig = theMat->getStress();
  stress(0) = sig * c * c;
  stress(1) = sig * s * s;
  stress(2) = sig * c * s;
  stress(3) = 0.0;
  stress(4) = 0.0;

  return stress ;
}


//send back the tangent 
const Matrix&  
PlateRebarMaterial::getTangent( )
{
  double tan = theMat->getTangent( ) ;

  tangent(0,0) = tan * c * c * c * c ;
  tangent(0,1) = tan * c * c * c * s ;
  tangent(0,2) = tan * c * c * s * s ;
  tangent(1,0) = tangent(0,1) ;
  tangent(1,1) = tangent(0,2) ;
  tangent(1,2) = tan * c * s * s * s ;
  tangent(2,0) = tangent(0,2) ;
  tangent(2,1) = tangent(1,2) ;
  tangent(2,2) = tan * s * s * s * s ;

  return tangent ;
}

const Matrix&  
PlateRebarMaterial::getInitialTangent
( )
{
  double tan = theMat->getInitialTangent( ) ;

  tangent(0,0) = tan * c * c * c * c ;
  tangent(0,1) = tan * c * c * c * s ;
  tangent(0,2) = tan * c * c * s * s ;
  tangent(1,0) = tangent(0,1) ;
  tangent(1,1) = tangent(0,2) ;
  tangent(1,2) = tan * c * s * s * s ;
  tangent(2,0) = tangent(0,2) ;
  tangent(2,1) = tangent(1,2) ;
  tangent(2,2) = tan * s * s * s * s ;

  return tangent ;
}


//print out data
void  
PlateRebarMaterial::Print( OPS_Stream &s, int flag )
{
  s << "PlateRebar Material tag: " << this->getTag() << endln ; 
  s << "using uniaxialmaterials : " << endln ;

  theMat->Print( s, flag ) ;

  return ;
}


int 
PlateRebarMaterial::sendSelf(int commitTag, Channel &theChannel) 
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
    opserr << "PlateRebarMaterial::sendSelf() - failed to send data" << endln;
    return res;
  }

  static Vector vecData(1);
  vecData(0) = angle;

  res = theChannel.sendVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateRebarMaterial::sendSelf() - failed to send data" << endln;
    return res;
  }

  // now send the materials data
  res += theMat->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlateRebarMaterial::sendSelf() - failed to send material1" << endln;

  return res;
}

int 
PlateRebarMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // recv an id containg the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlateRebarMaterial::sendSelf() - failed to receive id data" << endln;
    return res;
  }

  this->setTag(idData(0));
  int matClassTag = idData(1);
  if (theMat->getClassTag() != matClassTag) {
    if (theMat != 0) delete theMat;
    theMat = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMat == 0) {
      opserr << "PlateRebarMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMat->setDbTag(idData(2));

  static Vector vecData(1);
  res = theChannel.recvVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateRebarMaterial::sendSelf() - failed to receive vector data" << endln;
    return res;
  }
  angle = vecData(0);
  double rang = angle * 0.0174532925;
  c = cos(rang);
  s = sin(rang);

  // now receive the materials data
  res = theMat->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "PlateRebarMaterial::sendSelf() - failed to receive material1" << endln;
  
  return res;
}
 
