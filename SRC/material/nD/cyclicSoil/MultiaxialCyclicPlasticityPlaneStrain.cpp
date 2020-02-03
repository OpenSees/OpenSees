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

#include <MultiaxialCyclicPlasticityPlaneStrain.h>
#include <Channel.h>
#include  <FEM_ObjectBroker.h>

//static vectors and matrices
Vector MultiaxialCyclicPlasticityPlaneStrain :: strain_vec(3) ;
Vector MultiaxialCyclicPlasticityPlaneStrain :: stress_vec(3) ;
Matrix MultiaxialCyclicPlasticityPlaneStrain :: tangent_matrix(3,3) ;


//null constructor
MultiaxialCyclicPlasticityPlaneStrain ::  MultiaxialCyclicPlasticityPlaneStrain( ) : 
MultiaxialCyclicPlasticity( ) 
{  }


//full constructor
MultiaxialCyclicPlasticityPlaneStrain :: 
MultiaxialCyclicPlasticityPlaneStrain (int    tag, 
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
MultiaxialCyclicPlasticity( tag, ND_TAG_MultiaxialCyclicPlasticityPlaneStrain, rho, K, G, 
      Su, Ho_kin, Parameter_h, Parameter_m, Parameter_beta, Kcoeff, viscosity)
{ 

}



//elastic constructor
MultiaxialCyclicPlasticityPlaneStrain :: 
MultiaxialCyclicPlasticityPlaneStrain(   int    tag, double rho,    // add density by Gang Wang
                 double K, 
                 double G ) :
MultiaxialCyclicPlasticity( tag, ND_TAG_MultiaxialCyclicPlasticityPlaneStrain, rho, K, G )
{ 

}



//destructor
MultiaxialCyclicPlasticityPlaneStrain :: ~MultiaxialCyclicPlasticityPlaneStrain( ) 
{ } 




//make a clone of this material
NDMaterial* MultiaxialCyclicPlasticityPlaneStrain :: getCopy( ) 
{ 
  MultiaxialCyclicPlasticityPlaneStrain  *clone;
  clone = new MultiaxialCyclicPlasticityPlaneStrain() ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}



//send back type of material
const char* MultiaxialCyclicPlasticityPlaneStrain :: getType( ) const 
{
  return "PlaneStrain2D" ;
}


//send back order of strain in vector form
int MultiaxialCyclicPlasticityPlaneStrain :: getOrder( ) const 
{ 
  return 3 ; 
} 


//get the strain and integrate plasticity equations
int MultiaxialCyclicPlasticityPlaneStrain :: setTrialStrain( const Vector &strain_from_element) 
{
  strain.Zero( ) ;

  strain(0,0) =        strain_from_element(0) ;
  strain(1,1) =        strain_from_element(1) ;
  strain(0,1) = 0.50 * strain_from_element(2) ;
  strain(1,0) =        strain(0,1) ;

  if (this->MaterialStageID ==1) { 
	  this->elastic_integrator( ) ;
  } else if (this->MaterialStageID ==2) {	  
	  this->plastic_integrator( ) ;
  }
  return 0 ;
}


//unused trial strain functions
int MultiaxialCyclicPlasticityPlaneStrain :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

int MultiaxialCyclicPlasticityPlaneStrain :: setTrialStrainIncr( const Vector &v ) 
{
  static Vector newStrain(3);
  newStrain(0) = strain(0,0) + v(0);
  newStrain(1) = strain(1,1) + v(1);
  newStrain(2) = 2.0 * strain(0,1) + v(2);

  return this->setTrialStrain(newStrain);  
}

 

int MultiaxialCyclicPlasticityPlaneStrain :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
 
   // implemented April. 4, 2004, Gang Wang
	static Vector newStrain(3);
	newStrain(0) = strain(0,0) + v(0);
	newStrain(1) = strain(1,1) + v(1);
 	newStrain(2) = 2.0*strain(0,1) + v(2);
    	
	return this->setTrialStrain(newStrain);

}


//send back the strain
const Vector& MultiaxialCyclicPlasticityPlaneStrain :: getStrain( ) 
{
  strain_vec(0) =       strain(0,0) ;
  strain_vec(1) =       strain(1,1) ;
  strain_vec(2) = 2.0 * strain(0,1) ;

  return strain_vec ;
} 


//send back the stress 
const Vector& MultiaxialCyclicPlasticityPlaneStrain :: getStress( ) 
{
  stress_vec(0) = stress(0,0) ;
  stress_vec(1) = stress(1,1) ;
  stress_vec(2) = stress(0,1) ;

  return stress_vec ;
}

//send back the tangent 
const Matrix& MultiaxialCyclicPlasticityPlaneStrain :: getTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 ) 
  // 
       
  tangent_matrix(0,0) = tangent [0][0] [0][0] ;
  tangent_matrix(1,1) = tangent [1][1] [1][1] ;
  tangent_matrix(2,2) = tangent [0][1] [0][1] ;

  tangent_matrix(0,1) = tangent [0][0] [1][1] ;
  tangent_matrix(1,0) = tangent [1][1] [0][0] ;

  tangent_matrix(0,2) = tangent [0][0] [0][1] ;
  tangent_matrix(2,0) = tangent [0][1] [0][0] ;

  tangent_matrix(1,2) = tangent [1][1] [0][1] ;
  tangent_matrix(2,1) = tangent [0][1] [1][1] ;

  return tangent_matrix ;
} 


//send back the tangent 
const Matrix& MultiaxialCyclicPlasticityPlaneStrain :: getInitialTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 ) 
  // 

  this->doInitialTangent();

  tangent_matrix(0,0) = initialTangent [0][0] [0][0] ;
  tangent_matrix(1,1) = initialTangent [1][1] [1][1] ;
  tangent_matrix(2,2) = initialTangent [0][1] [0][1] ;

  tangent_matrix(0,1) = initialTangent [0][0] [1][1] ;
  tangent_matrix(1,0) = initialTangent [1][1] [0][0] ;

  tangent_matrix(0,2) = initialTangent [0][0] [0][1] ;
  tangent_matrix(2,0) = initialTangent [0][1] [0][0] ;

  tangent_matrix(1,2) = initialTangent [1][1] [0][1] ;
  tangent_matrix(2,1) = initialTangent [0][1] [1][1] ;

  return tangent_matrix ;
} 

/*
int 
MultiaxialCyclicPlasticityPlaneStrain::commitState( ) 
{
  epsilon_p_n = epsilon_p_nplus1;
  xi_n        = xi_nplus1;

  return 0;
}

int 
MultiaxialCyclicPlasticityPlaneStrain::revertToLastCommit( ) {

  return 0;
}


int 
MultiaxialCyclicPlasticityPlaneStrain::revertToStart( ) 
{
  this->zero( ) ;

  return 0;
}
 

int
MultiaxialCyclicPlasticityPlaneStrain::sendSelf (int commitTag, Channel &theChannel)
{
 
  return -1;
}

int
MultiaxialCyclicPlasticityPlaneStrain::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  
  return -1;
}


*/





