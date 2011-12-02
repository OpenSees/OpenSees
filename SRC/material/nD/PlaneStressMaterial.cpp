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
// $Date: 2003-02-14 23:01:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlaneStressMaterial.cpp,v $

//
// Ed "C++" Love
//
// Generic Plane Stress Material
//


#include <PlaneStressMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static vector and matrices
Vector  PlaneStressMaterial::stress(3) ;
Matrix  PlaneStressMaterial::tangent(3,3) ;



//null constructor
PlaneStressMaterial::PlaneStressMaterial( ) : 
NDMaterial(0, ND_TAG_PlaneStressMaterial ), 
strain(3) 
{ }



//full constructor
PlaneStressMaterial::PlaneStressMaterial(    
				   int tag, 
                                   NDMaterial &the3DMaterial ) :
NDMaterial( tag, ND_TAG_PlaneStressMaterial ),
strain(3)
{
  theMaterial = the3DMaterial.getCopy("ThreeDimensional") ;

  Tstrain22 = 0.0 ;
  Tgamma02 = 0.0 ;
  Tgamma12 = 0.0 ;
  Cstrain22 = 0.0 ;
  Cgamma02 = 0.0 ;
  Cgamma12 = 0.0 ;
}



//destructor
PlaneStressMaterial::~PlaneStressMaterial( ) 
{ 
  delete theMaterial ;
} 



//make a clone of this material
NDMaterial*
PlaneStressMaterial::getCopy( ) 
{
  PlaneStressMaterial *clone ;   //new instance of this class

  clone = new PlaneStressMaterial( this->getTag(), 
                                   *theMaterial ) ; //make the copy

  clone->Tstrain22 = this->Tstrain22 ;
  clone->Tgamma02  = this->Tgamma02 ;
  clone->Tgamma12  = this->Tgamma12 ;
  clone->Cstrain22 = this->Cstrain22 ;
  clone->Cgamma02  = this->Cgamma02 ;
  clone->Cgamma12  = this->Cgamma12 ;

  return clone ;
}


//make a clone of this material
NDMaterial* 
PlaneStressMaterial::getCopy( const char *type ) 
{
  return this->getCopy( ) ;
}


//send back order of strain in vector form
int 
PlaneStressMaterial::getOrder( ) const
{
  return 3 ;
}


const char*
PlaneStressMaterial::getType( ) const 
{
  return "PlaneStress" ; 
}



//swap history variables
int 
PlaneStressMaterial::commitState( ) 
{
  Cstrain22 = Tstrain22;
  Cgamma02 = Tgamma02;
  Cgamma12 = Tgamma12;

  return theMaterial->commitState( ) ;
}



//revert to last saved state
int 
PlaneStressMaterial::revertToLastCommit( )
{
  Tstrain22 = Cstrain22;
  Tgamma02 = Cgamma02;
  Tgamma12 = Cgamma12;

  return theMaterial->revertToLastCommit( )  ;
}


//revert to start
int
PlaneStressMaterial::revertToStart( )
{
  this->Tstrain22 = 0.0 ;
  this->Tgamma12  = 0.0 ;
  this->Tgamma02  = 0.0 ;
  this->Cstrain22 = 0.0 ;
  this->Cgamma12  = 0.0 ;
  this->Cgamma02  = 0.0 ;
  
  strain.Zero();

  return theMaterial->revertToStart( ) ;
}


//mass per unit volume
double
PlaneStressMaterial::getRho( )
{
  return theMaterial->getRho( ) ;
}


//receive the strain
int 
PlaneStressMaterial::setTrialStrain( const Vector &strainFromElement )
{
  static const double tolerance = 1.0e-08 ;

  this->strain(0) = strainFromElement(0) ;
  this->strain(1) = strainFromElement(1) ;
  this->strain(2) = strainFromElement(2) ;

  // return theMaterial->setTrialStrain( threeDstrain ) ;
  double norm ;

  static Vector outOfPlaneStress(3) ;
  static Vector strainIncrement(3) ;
  static Vector threeDstress(6) ;
  static Vector threeDstrain(6) ;
  static Matrix threeDtangent(6,6) ;
  static Vector threeDstressCopy(6) ; 
  static Matrix threeDtangentCopy(6,6) ;

  static Matrix dd22(3,3) ;

  int i, j ;
  int ii, jj ;

  //newton loop to solve for out-of-plane strains
  do {

    //set three dimensional strain
    threeDstrain(0) = this->strain(0) ;
    threeDstrain(1) = this->strain(1) ;
    threeDstrain(2) = this->Tstrain22 ;
    threeDstrain(3) = this->strain(2) ; 
    threeDstrain(4) = this->Tgamma12 ;
    threeDstrain(5) = this->Tgamma02 ;

    if (theMaterial->setTrialStrain( threeDstrain ) < 0) {
      opserr << "PlaneStressMaterial::setTrialStrain() - setTrialStrain in material failed with strain " << threeDstrain;
      return -1;
    }

    //three dimensional stress
    threeDstress = theMaterial->getStress( ) ;

    //three dimensional tangent 
    threeDtangent = theMaterial->getTangent( ) ;

    //NDmaterial strain order          = 11, 22, 33, 12, 23, 31 
    //PlaneStressMaterial strain order = 11, 22, 12, 33, 23, 31 

    //swap matrix indices to sort out-of-plane components 
    for ( i=0; i<6; i++ ) {

      ii = this->indexMap(i) ;

      threeDstressCopy(ii) = threeDstress(i) ;

      for ( j=0; j<6; j++ ) {

	jj = this->indexMap(j) ;
	
	threeDtangentCopy(ii,jj) = threeDtangent(i,j) ;

      }//end for j
       
    }//end for i


    //partitioned stresses and tangent
    for ( i=0; i<3; i++ ) {

      outOfPlaneStress(i) = threeDstressCopy(i+3) ;

      for ( j=0; j<3; j++ ) 
	dd22(i,j) = threeDtangentCopy(i+3,j+3) ;

    }//end for i


    //set norm
    norm = outOfPlaneStress.Norm( ) ;

    //int Solve(const Vector &V, Vector &res) const;
    //int Solve(const Matrix &M, Matrix &res) const;
    //condensation 
    dd22.Solve( outOfPlaneStress, strainIncrement ) ;

    //update out of plane strains
    this->Tstrain22 -= strainIncrement(0) ;
    this->Tgamma12  -= strainIncrement(1) ;
    this->Tgamma02  -= strainIncrement(2) ;

  } while ( norm > tolerance ) ;

  return 0;
}


//send back the strain
const Vector& 
PlaneStressMaterial::getStrain( )
{
  return this->strain ;
}


//send back the stress 
const Vector&  
PlaneStressMaterial::getStress( )
{
  //three dimensional stress
  const Vector &threeDstress = theMaterial->getStress();
  static Vector threeDstressCopy(6);

  //partitioned stresses and tangent
  //swap matrix indices to sort out-of-plane components 
  int i, ii;
  for ( i=0; i<6; i++ ) {

    ii = this->indexMap(i) ;

    threeDstressCopy(ii) = threeDstress(i) ;
  }

  for ( i=0; i<3; i++ ) 
    this->stress(i)     = threeDstressCopy(i) ;
  
  return this->stress ;
}


//send back the tangent 
const Matrix&  
PlaneStressMaterial::getTangent( )
{
  static Matrix dd11(3,3) ;
  static Matrix dd12(3,3) ;
  static Matrix dd21(3,3) ;
  static Matrix dd22(3,3) ;

  static Matrix dd22invdd21(3,3) ;
  static Matrix threeDtangentCopy(6,6);

  //three dimensional tangent 
  const Matrix &threeDtangent = theMaterial->getTangent( ) ;

  //NDmaterial strain order          = 11, 22, 33, 12, 23, 31 
  //PlaneStressMaterial strain order = 11, 22, 12, 33, 23, 31 

  //swap matrix indices to sort out-of-plane components 
  int i,j, ii, jj;

  for ( i=0; i<6; i++ ) {

    ii = this->indexMap(i) ;

    for ( j=0; j<6; j++ ) {
      jj = this->indexMap(j) ;
      threeDtangentCopy(ii,jj) = threeDtangent(i,j) ;
    }//end for j

  }//end for i


  //out of plane stress and tangents
  for ( i=0; i<3; i++ ) {
    for ( j=0; j<3; j++ ) {
	
      dd11(i,j) = threeDtangentCopy(i,  j  ) ;
      dd12(i,j) = threeDtangentCopy(i,  j+3) ;
      dd21(i,j) = threeDtangentCopy(i+3,j  ) ;
      dd22(i,j) = threeDtangentCopy(i+3,j+3) ;

    }//end for j
  }//end for i

  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  dd22.Solve( dd21, dd22invdd21 ) ;
  this->tangent   = dd11 ; 
  this->tangent  -= ( dd12*dd22invdd21 ) ;

  return this->tangent ;
}



int 
PlaneStressMaterial::indexMap( int i )
{
  int ii ;

  if ( i == 2 ) 
    ii = 3 ;
  else if ( i == 3 )
    ii = 2 ;
  else 
    ii = i ;

  return ii ;
}



//print out data
void  
PlaneStressMaterial::Print( OPS_Stream &s, int flag )
{
  s << "General Plane Stress Material \n" ;
  s << " Tag: " << this->getTag() << "\n" ; 
  s << "using the 3D material : \n" ;

  theMaterial->Print( s, flag ) ;

  return ;
}


int 
PlaneStressMaterial::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  // put tag and assocaited materials class and database tags into an id and send it
  static ID idData(3);
  idData(0) = this->getTag();
  idData(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  idData(2) = matDbTag;

  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PlaneStressMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  // put the strains in a vector and send it
  static Vector vecData(3);
  vecData(0) = Cstrain22;
  vecData(1) = Cgamma02;
  vecData(2) = Cgamma12;

  res = theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "PlaneStressMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlaneStressMaterial::sendSelf() - failed to send vector material\n";

  return res;
}

int 
PlaneStressMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id containg the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PlaneStressMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  this->setTag(idData(0));
  int matClassTag = idData(1);
  
  // if the associated material has not yet been created or is of the wrong type
  // create a new material for recvSelf later
  if (theMaterial == 0 || theMaterial->getClassTag() != matClassTag) {
    if (theMaterial != 0)
      delete theMaterial;
    theMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "PlaneStressMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(3);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "PlaneStressMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  Cstrain22 = vecData(0);
  Cgamma02 = vecData(1);
  Cgamma12  = vecData(2);

  Tstrain22 = Cstrain22;
  Tgamma02 = Cgamma02;
  Tgamma12  = Cgamma12;

  // now receive the materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "PlaneStressMaterial::sendSelf() - failed to send vector material\n";
  
  return res;
}
 


