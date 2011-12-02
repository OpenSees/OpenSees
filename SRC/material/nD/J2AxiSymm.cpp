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
// $Date: 2000-12-13 08:10:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2AxiSymm.cpp,v $

// Written: Ed "C++" Love
//
// J2AxiSymmetric isotropic hardening material class
// 
//  Elastic Model
//  sigma = K*trace(epsilion_elastic) + (2*G)*dev(epsilon_elastic)
//
//  Yield Function
//  phi(sigma,q) = || dev(sigma) ||  - sqrt(2/3)*q(xi) 
//
//  Saturation Isotropic Hardening with linear term
//  q(xi) = simga_0 + (sigma_infty - sigma_0)*exp(-delta*xi) + H*xi 
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
//	                eps_22 		      
//                    2 eps_01   }   <--- note the 2
// 
//  set eta := 0 for rate independent case
//

#include <J2AxiSymm.h>

//static vectors and matrices
Vector J2AxiSymm :: strain_vec(4) ;
Vector J2AxiSymm :: stress_vec(4) ;
Matrix J2AxiSymm :: tangent_matrix(4,4) ;


//null constructor
J2AxiSymm ::  J2AxiSymm( ) : 
J2Plasticity( ) 
{  }


//full constructor
J2AxiSymm :: 
J2AxiSymm(   int    tag, 
                 double K,
                 double G,
                 double yield0,
                 double yield_infty,
                 double d,
                 double H,
                 double viscosity ) : 
J2Plasticity( tag, ND_TAG_J2AxiSymm, 
             K, G, yield0, yield_infty, d, H, viscosity )
{ }


//elastic constructor
J2AxiSymm :: 
J2AxiSymm(   int    tag, 
                 double K, 
                 double G ) :
J2Plasticity( tag, ND_TAG_J2AxiSymm, K, G )
{ }



//destructor
J2AxiSymm :: ~J2AxiSymm( ) 
{ } 


//make a clone of this material
NDMaterial* J2AxiSymm :: getCopy( ) 
{ 
  J2AxiSymm  *clone;
  clone = new J2AxiSymm() ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}


//send back type of material
const char* J2AxiSymm :: getType( ) const 
{
  return "AxiSymmetric2D" ;
}


//send back order of strain in vector form
int J2AxiSymm :: getOrder( ) const 
{ 
  return 4 ; 
} 


//get the strain and integrate plasticity equations
int J2AxiSymm :: setTrialStrain( const Vector &strain_from_element) 
{
  strain.Zero( ) ;

  strain(0,0) =        strain_from_element(0) ;
  strain(1,1) =        strain_from_element(1) ;
  strain(2,2) =        strain_from_element(2) ;

  strain(0,1) = 0.50 * strain_from_element(3) ;
  strain(1,0) =        strain(0,1) ;

  this->plastic_integrator( ) ;

  return 0 ;
}


//unused trial strain functions
int J2AxiSymm :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

int J2AxiSymm :: setTrialStrainIncr( const Vector &v ) 
{
    return -1 ;
}

int J2AxiSymm :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
    return -1 ;
}



//send back the strain
const Vector& J2AxiSymm :: getStrain( ) 
{
  strain_vec(0) =       strain(0,0) ;
  strain_vec(1) =       strain(1,1) ;
  strain_vec(2) =       strain(2,2) ;

  strain_vec(3) = 2.0 * strain(0,1) ;

  return strain_vec ;
} 


//send back the stress 
const Vector& J2AxiSymm :: getStress( ) 
{
  stress_vec(0) = stress(0,0) ;
  stress_vec(1) = stress(1,1) ;
  stress_vec(2) = stress(2,2) ;

  stress_vec(3) = stress(0,1) ;

  return stress_vec ;
}

//send back the tangent 
const Matrix& J2AxiSymm :: getTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           2 2   
  //   3           0 1  ( or 1 0 )
  
  int ii, jj ;
  int i, j, k, l ;

  for ( ii = 0; ii < 4; ii++ ) {
    for ( jj = 0; jj < 4; jj++ ) {

      index_map( ii, i, j ) ;
      index_map( jj, k, l ) ;

      tangent_matrix(ii,jj) = tangent[i][j][k][l] ;

    } //end for j
  } //end for i

  return tangent_matrix ;
} 

//this is mike's problem
int J2AxiSymm :: setTrialStrain(const Tensor &v) 
{
  return -1 ;
}

int J2AxiSymm :: setTrialStrain(const Tensor &v, const Tensor &r)     
{
  return -1 ;
}

int J2AxiSymm :: setTrialStrainIncr(const Tensor &v) 
{
  return -1 ;
}

int J2AxiSymm :: setTrialStrainIncr(const Tensor &v, const Tensor &r) 
{
  return -1 ;
}

const Tensor& J2AxiSymm :: getTangentTensor( ) 
{
  return rank4 ;
}

const Tensor& J2AxiSymm :: getStressTensor( ) 
{
  return rank2 ;
}

const Tensor& J2AxiSymm :: getStrainTensor( ) 
{
  return rank2 ;
}


//this is frank's problem
int J2AxiSymm :: sendSelf(int commitTag, Channel &theChannel)
{
  return -1 ;
}

int J2AxiSymm :: recvSelf(int commitTag, Channel &theChannel, 
	                      FEM_ObjectBroker &theBroker)
{
  return -1 ;
}




