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
** ****************************************************************** */


/*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                          |
 |              MultiaxialCyclicPlasticity  NDMaterial                      |
 +                                                                          +
 |--------------------------------------------------------------------------|
 |                                                                          |
 +             Authors: Gang Wang  AND  Professor Nicholas Sitar            +
 |                                                                          |
 |             Department of Civil and Environmental Engineering            |
 +             University of California, Berkeley, CA 94720, USA            +
 |                                                                          |
 |             Email: wang@ce.berkeley.edu (G.W.)                           |
 |                                                                          |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/


#include <MultiaxialCyclicPlasticity3D.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static vectors and matrices
Vector MultiaxialCyclicPlasticity3D :: strain_vec(6) ;
Vector MultiaxialCyclicPlasticity3D :: stress_vec(6) ;
Matrix MultiaxialCyclicPlasticity3D :: tangent_matrix(6,6) ;


//null constructor
MultiaxialCyclicPlasticity3D ::  MultiaxialCyclicPlasticity3D( ) : 
MultiaxialCyclicPlasticity( ) 
{  }


//full constructor
MultiaxialCyclicPlasticity3D :: 
MultiaxialCyclicPlasticity3D(   int    tag,
				double rho,
				double K,
				double G,
				double Su,
				double Ho_kin,
				double Parameter_h,
				double Parameter_m,
				double Parameter_beta,
				double Kcoeff,
				double viscosity ) :
MultiaxialCyclicPlasticity( tag, ND_TAG_MultiaxialCyclicPlasticity3D, rho, K, G, 
		Su, Ho_kin, Parameter_h, Parameter_m, Parameter_beta, Kcoeff, viscosity)
{ 

}


//elastic constructor
MultiaxialCyclicPlasticity3D :: 
MultiaxialCyclicPlasticity3D(   int    tag, 
				 double rho,
                 double K, 
                 double G ) :
MultiaxialCyclicPlasticity( tag, ND_TAG_MultiaxialCyclicPlasticity3D, rho, K, G )
{ 

}



//destructor
MultiaxialCyclicPlasticity3D :: ~MultiaxialCyclicPlasticity3D( ) 
{ } 


//make a clone of this material
NDMaterial* MultiaxialCyclicPlasticity3D :: getCopy( ) 
{ 
  MultiaxialCyclicPlasticity3D  *clone;
  clone = new MultiaxialCyclicPlasticity3D( ) ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}


//send back type of material
const char* MultiaxialCyclicPlasticity3D :: getType( ) const 
{
  return "ThreeDimensional" ;
}


//send back order of strain in vector form
int MultiaxialCyclicPlasticity3D :: getOrder( ) const 
{ 
  return 6 ; 
} 


//get the strain and integrate plasticity equations
int MultiaxialCyclicPlasticity3D :: setTrialStrain( const Vector &strain_from_element) 
{
  strain.Zero( ) ;

  strain(0,0) =        strain_from_element(0) ;
  strain(1,1) =        strain_from_element(1) ;
  strain(2,2) =        strain_from_element(2) ;

  strain(0,1) = 0.50 * strain_from_element(3) ;
  strain(1,0) =        strain(0,1) ;

  strain(1,2) = 0.50 * strain_from_element(4) ;
  strain(2,1) =        strain(1,2) ;
  
  strain(2,0) = 0.50 * strain_from_element(5) ;
  strain(0,2) =        strain(2,0) ;

  
  //opserr<<"DP:setTrialStrain:MaterialStageID= "<<this->MaterialStageID<<endln;

  if (this->MaterialStageID ==1) { 
	  //opserr<<"DP:setTrialStrain:elastic K0"<<endln;
	  this->elastic_integrator( ) ;
  } else if (this->MaterialStageID ==2) {	  
	  //opserr<<"DP:setTrialStrain:plastic"<<endln;
	  // initialize history variables here in first call???
	  
	  this->plastic_integrator( ) ;
  }
  


  return 0 ;
}


//unused trial strain functions
int MultiaxialCyclicPlasticity3D :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

int MultiaxialCyclicPlasticity3D :: setTrialStrainIncr( const Vector &v ) 
{
  static Vector newStrain(6);
  newStrain(0) = strain(0,0) + v(0);
  newStrain(1) = strain(1,1) + v(1);
  newStrain(2) = strain(2,2) + v(2);
  newStrain(3) = 2.0*strain(0,1) + v(3);
  newStrain(4) = 2.0*strain(1,2) + v(4);
  newStrain(5) = 2.0*strain(2,0) + v(5);
  
  return this->setTrialStrain(newStrain);
}

int MultiaxialCyclicPlasticity3D :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
  return this->setTrialStrainIncr(v);
}



//send back the strain
const Vector& MultiaxialCyclicPlasticity3D :: getStrain( ) 
{
  strain_vec(0) =       strain(0,0) ;
  strain_vec(1) =       strain(1,1) ;
  strain_vec(2) =       strain(2,2) ;

  strain_vec(3) = 2.0 * strain(0,1) ;

  strain_vec(4) = 2.0 * strain(1,2) ;

  strain_vec(5) = 2.0 * strain(2,0) ;

  return strain_vec ;
} 


//send back the stress 
const Vector& MultiaxialCyclicPlasticity3D :: getStress( ) 
{
  stress_vec(0) = stress(0,0) ;
  stress_vec(1) = stress(1,1) ;
  stress_vec(2) = stress(2,2) ;

  stress_vec(3) = stress(0,1) ;

  stress_vec(4) = stress(1,2) ;
  
  stress_vec(5) = stress(2,0) ;

  return stress_vec ;
}

//send back the tangent 
const Matrix& MultiaxialCyclicPlasticity3D :: getTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           2 2   
  //   3           0 1  ( or 1 0 )
  //   4           1 2  ( or 2 1 )
  //   5           2 0  ( or 0 2 ) 
    
  int ii, jj ;
  int i, j, k, l ;

  for ( ii = 0; ii < 6; ii++ ) {
    for ( jj = 0; jj < 6; jj++ ) {

      index_map( ii, i, j ) ;
      index_map( jj, k, l ) ;

      tangent_matrix(ii,jj) = tangent[i][j][k][l] ;

    } //end for j
  } //end for i

  return tangent_matrix ;
} 

//send back the tangent 
const Matrix& MultiaxialCyclicPlasticity3D :: getInitialTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           2 2   
  //   3           0 1  ( or 1 0 )
  //   4           1 2  ( or 2 1 )
  //   5           2 0  ( or 0 2 ) 
    
  int ii, jj ;
  int i, j, k, l ;

  this->doInitialTangent();

  for ( ii = 0; ii < 6; ii++ ) {
    for ( jj = 0; jj < 6; jj++ ) {

      index_map( ii, i, j ) ;
      index_map( jj, k, l ) ;

      tangent_matrix(ii,jj) = initialTangent[i][j][k][l] ;

    } //end for j
  } //end for i

  return tangent_matrix ;
} 





