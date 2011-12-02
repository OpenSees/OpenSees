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
                                                                        
// $Revision: 1.1 $
// $Date: 2000-12-13 08:12:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2PlaneStrain.cpp,v $

// Written: Ed "C++" Love
//
// J2PlaneStrain isotropic hardening material class
// 
//  Elastic Model
//  sigma = K*trace(epsilion_elastic) + (2*G)*dev(epsilon_elastic)
//
//  Yield Function
//  phi(sigma,q) = || dev(sigma) ||  - sqrt(2/3)*q(xi) 
//
//  Saturation Isotropic Hardening with linear term
//  q(xi) = simga_infty + (sigma_0 - sigma_infty)*exp(-delta*xi) + H*xi 
//
//  Flow Rules
//  \dot{epsilon_p} =  gamma * d_phi/d_sigma
//  \dot{xi}        = -gamma * d_phi/d_q 
//
//  Linear Viscosity 
//  gamma = phi / eta  ( if phi > 0 ) 
//
//  Backward Euler Integration Routine 
//  Yield condition enforced at time n+1 
//
//  Send strains in following format :
// 
//     strain_vec = {   eps_00
//                      eps_11
//                    2 eps_01   }   <--- note the 2
// 
//  set eta := 0 for rate independent case
//

#include <J2PlaneStrain.h>

//static vectors and matrices
Vector J2PlaneStrain :: strain_vec(3) ;
Vector J2PlaneStrain :: stress_vec(3) ;
Matrix J2PlaneStrain :: tangent_matrix(3,3) ;


//null constructor
J2PlaneStrain ::  J2PlaneStrain( ) : 
J2Plasticity( ) 
{  }


//full constructor
J2PlaneStrain :: 
J2PlaneStrain(   int    tag, 
                 double K,
                 double G,
                 double yield0,
                 double yield_infty,
                 double d,
                 double H,
                 double viscosity ) : 
J2Plasticity( tag, ND_TAG_J2PlaneStrain, 
             K, G, yield0, yield_infty, d, H, viscosity )
{ }


//elastic constructor
J2PlaneStrain :: 
J2PlaneStrain(   int    tag, 
                 double K, 
                 double G ) :
J2Plasticity( tag, ND_TAG_J2PlaneStrain, K, G )
{ }



//destructor
J2PlaneStrain :: ~J2PlaneStrain( ) 
{ } 


//make a clone of this material
NDMaterial* J2PlaneStrain :: getCopy( ) 
{ 
  J2PlaneStrain  *clone;
  clone = new J2PlaneStrain() ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}



//send back type of material
const char* J2PlaneStrain :: getType( ) const 
{
  return "PlaneStrain2D" ;
}


//send back order of strain in vector form
int J2PlaneStrain :: getOrder( ) const 
{ 
  return 3 ; 
} 


//get the strain and integrate plasticity equations
int J2PlaneStrain :: setTrialStrain( const Vector &strain_from_element) 
{
  strain.Zero( ) ;

  strain(0,0) =        strain_from_element(0) ;
  strain(1,1) =        strain_from_element(1) ;
  strain(0,1) = 0.50 * strain_from_element(2) ;
  strain(1,0) =        strain(0,1) ;

  this->plastic_integrator( ) ;

  return 0 ;
}


//unused trial strain functions
int J2PlaneStrain :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

int J2PlaneStrain :: setTrialStrainIncr( const Vector &v ) 
{
    return -1 ;
}

int J2PlaneStrain :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
    return -1 ;
}



//send back the strain
const Vector& J2PlaneStrain :: getStrain( ) 
{
  strain_vec(0) =       strain(0,0) ;
  strain_vec(1) =       strain(1,1) ;
  strain_vec(2) = 2.0 * strain(0,1) ;

  return strain_vec ;
} 


//send back the stress 
const Vector& J2PlaneStrain :: getStress( ) 
{
  stress_vec(0) = stress(0,0) ;
  stress_vec(1) = stress(1,1) ;
  stress_vec(2) = stress(0,1) ;

  return stress_vec ;
}

//send back the tangent 
const Matrix& J2PlaneStrain :: getTangent( ) 
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

//this is mike's problem
int J2PlaneStrain :: setTrialStrain(const Tensor &v) 
{
  return -1 ;
}

int J2PlaneStrain :: setTrialStrain(const Tensor &v, const Tensor &r)     
{
  return -1 ;
}

int J2PlaneStrain :: setTrialStrainIncr(const Tensor &v) 
{
  return -1 ;
}

int J2PlaneStrain :: setTrialStrainIncr(const Tensor &v, const Tensor &r) 
{
  return -1 ;
}

const Tensor& J2PlaneStrain :: getTangentTensor( ) 
{
  return rank4 ;
}

const Tensor& J2PlaneStrain :: getStressTensor( ) 
{
  return rank2 ;
}

const Tensor& J2PlaneStrain :: getStrainTensor( ) 
{
  return rank2 ;
}


//this is frank's problem
int J2PlaneStrain :: sendSelf(int commitTag, Channel &theChannel)
{
  return -1 ;
}

int J2PlaneStrain :: recvSelf(int commitTag, Channel &theChannel, 
	                      FEM_ObjectBroker &theBroker)
{
  return -1 ;
}









