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
// $Date: 2012-05-28 22:03:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateFromPlaneStressMaterialThermal.cpp,v $

//
// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Generic Plate Material from Plane Stress Material
//
/* Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
concrete high-rise building induced by extreme earthquakes, 
Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723*/
//Modified by Liming Jiang for considering elevated temperature



#include <PlateFromPlaneStressMaterialThermal.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>

//static vector and matrices
Vector  PlateFromPlaneStressMaterialThermal::stress(5) ;
Matrix  PlateFromPlaneStressMaterialThermal::tangent(5,5) ;

//null constructor
PlateFromPlaneStressMaterialThermal::PlateFromPlaneStressMaterialThermal( ) : 
NDMaterial(0, ND_TAG_PlateFromPlaneStressMaterialThermal ), 
strain(5) 
{ }


//full constructor
PlateFromPlaneStressMaterialThermal::PlateFromPlaneStressMaterialThermal(    
				   int tag, NDMaterial &ndMat, double g ) :
NDMaterial( tag, ND_TAG_PlateFromPlaneStressMaterialThermal ),
strain(5),gmod(g)
{
  theMat = ndMat.getCopy("PlaneStress") ;
}


//destructor
PlateFromPlaneStressMaterialThermal::~PlateFromPlaneStressMaterialThermal( ) 
{ 
  if (theMat != 0) delete theMat ;
} 



//make a clone of this material
NDMaterial*
PlateFromPlaneStressMaterialThermal::getCopy( ) 
{
  PlateFromPlaneStressMaterialThermal *clone ;   //new instance of this class

  clone = new PlateFromPlaneStressMaterialThermal( this->getTag(), *theMat, gmod ) ;

  return clone ;
}


//make a clone of this material
NDMaterial* 
PlateFromPlaneStressMaterialThermal::getCopy( const char *type ) 
{
  if (strcmp(type, this->getType()) == 0)
    return this->getCopy( ) ;
  else
    return 0;
}


//send back order of strain in vector form
int 
PlateFromPlaneStressMaterialThermal::getOrder( ) const
{
  return 5 ;
}


const char*
PlateFromPlaneStressMaterialThermal::getType( ) const 
{
  return "PlateFiberThermal" ; 
}



//swap history variables
int 
PlateFromPlaneStressMaterialThermal::commitState( ) 
{
  return theMat->commitState( ) ;
}



//revert to last saved state
int 
PlateFromPlaneStressMaterialThermal::revertToLastCommit( )
{
  return theMat->revertToLastCommit( )  ;
}


//revert to start
int
PlateFromPlaneStressMaterialThermal::revertToStart( )
{

  strain.Zero();

  return theMat->revertToStart( ) ;
}


//mass per unit volume
double
PlateFromPlaneStressMaterialThermal::getRho( )
{
  return theMat->getRho( ) ;
}


//receive the strain
int 
PlateFromPlaneStressMaterialThermal::setTrialStrain( const Vector &strainFromElement )
{
  strain(0) = strainFromElement(0) ;
  strain(1) = strainFromElement(1) ;
  strain(2) = strainFromElement(2) ;
  strain(3) = strainFromElement(3) ;
  strain(4) = strainFromElement(4) ;

  static Vector PSStrain(3) ;
  
  PSStrain(0) = strain(0);
  PSStrain(1) = strain(1);
  PSStrain(2) = strain(2);
  
  return theMat->setTrialStrain(PSStrain);
}


//send back the strain
const Vector& 
PlateFromPlaneStressMaterialThermal::getStrain( )
{
  return strain ;
}


//send back the stress 
const Vector&  
PlateFromPlaneStressMaterialThermal::getStress( )
{
  //three dimensional stress
  const Vector &PSStress = theMat->getStress();

  stress(0) = PSStress(0);
  stress(1) = PSStress(1);
  stress(2) = PSStress(2);
  stress(3) = gmod * strain(3);
  stress(4) = gmod * strain(4);

  return stress ;
}


//send back the tangent 
const Matrix&  
PlateFromPlaneStressMaterialThermal::getTangent( )
{
  const Matrix PSTangent = theMat->getTangent( ) ;

  tangent.Zero();
  
  tangent(0,0) = PSTangent(0,0);
  tangent(0,1) = PSTangent(0,1);
  tangent(0,2) = PSTangent(0,2);
  tangent(1,0) = PSTangent(1,0);
  tangent(1,1) = PSTangent(1,1);
  tangent(1,2) = PSTangent(1,2);
  tangent(2,0) = PSTangent(2,0);
  tangent(2,1) = PSTangent(2,1);
  tangent(2,2) = PSTangent(2,2);
  tangent(3,3) = gmod;
  tangent(4,4) = gmod;

  return tangent ;
}

//send back the tangent 
const Matrix&  
PlateFromPlaneStressMaterialThermal::getInitialTangent
( )
{
  const Matrix PSTangent = theMat->getInitialTangent( ) ;

  tangent.Zero();
  
  tangent(0,0) = PSTangent(0,0);
  tangent(0,1) = PSTangent(0,1);
  tangent(0,2) = PSTangent(0,2);
  tangent(1,0) = PSTangent(1,0);
  tangent(1,1) = PSTangent(1,1);
  tangent(1,2) = PSTangent(1,2);
  tangent(2,0) = PSTangent(2,0);
  tangent(2,1) = PSTangent(2,1);
  tangent(2,2) = PSTangent(2,2);
  tangent(3,3) = gmod;
  tangent(4,4) = gmod;

  return tangent ;
}


//print out data
void  
PlateFromPlaneStressMaterialThermal::Print( OPS_Stream &s, int flag )
{
  s << "PlateFromPlaneStress Material tag: " << this->getTag() << "" << endln ; 
  s << "using PlaneStress material : " << endln ;

  theMat->Print( s, flag ) ;

}


int 
PlateFromPlaneStressMaterialThermal::sendSelf(int commitTag, Channel &theChannel) 
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
    opserr << "PlateFromPlaneStressMaterialThermal::sendSelf() - failed to send data" << endln;
    return res;
  }

  static Vector vecData(1);
  vecData(0) = gmod;

  res = theChannel.sendVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateFromPlaneStressMaterialThermal::sendSelf() - failed to send data" << endln;
    return res;
  }

  // now send the materials data
  res += theMat->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlateFromPlaneStressMaterialThermal::sendSelf() - failed to send material1" << endln;

  return res;
}

int 
PlateFromPlaneStressMaterialThermal::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // recv an id containing the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlateFromPlaneStressMaterialThermal::sendSelf() - failed to receive id data" << endln;
    return res;
  }

  this->setTag(idData(0));
  int matClassTag = idData(1);
  if (theMat->getClassTag() != matClassTag) {
    if (theMat != 0) delete theMat;
    theMat = theBroker.getNewNDMaterial(matClassTag);
    if (theMat == 0) {
      opserr << "PlateFromPlaneStressMaterialThermal::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMat->setDbTag(idData(2));

  static Vector vecData(1);
  res = theChannel.recvVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateFromPlaneStressMaterialThermal::sendSelf() - failed to receive vector data" << endln;
    return res;
  }
  gmod = vecData(0);

  // now receive the materials data
  res = theMat->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "PlateFromPlaneStressMaterialThermal::sendSelf() - failed to receive material1" << endln;
  
  return res;
}
 

double 
PlateFromPlaneStressMaterialThermal::getThermalTangentAndElongation(double &TempT, double&ET, double&Elong)
{

theMat->setThermalTangentAndElongation(TempT,ET,Elong );
return 0;
}

const Vector&
PlateFromPlaneStressMaterialThermal::getTempAndElong()
{
	//return theMaterial->getTempAndElong( );
   static Vector returnedVec = Vector(2);
	returnedVec(0)= theMat->getTempAndElong( )(0);
	returnedVec(1) = theMat->getTempAndElong( )(1);
	return returnedVec;
}