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
                                                                        
// $Revision: 1.7 $
// $Date: 2008-10-20 22:23:03 $
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
#include <Channel.h>
#include <FEM_ObjectBroker.h>

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
	     double viscosity,
	     double rho) : 
J2Plasticity( tag, ND_TAG_J2AxiSymm, 
	      K, G, yield0, yield_infty, d, H, viscosity, rho )
{ 

}


//elastic constructor
J2AxiSymm :: 
J2AxiSymm(   int    tag, 
                 double K, 
                 double G ) :
J2Plasticity( tag, ND_TAG_J2AxiSymm, K, G )
{ 

}



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

//send back the tangent 
const Matrix& J2AxiSymm :: getInitialTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           2 2   
  //   3           0 1  ( or 1 0 )

  this->doInitialTangent();

  int ii, jj ;
  int i, j, k, l ;

  for ( ii = 0; ii < 4; ii++ ) {
    for ( jj = 0; jj < 4; jj++ ) {

      index_map( ii, i, j ) ;
      index_map( jj, k, l ) ;

      tangent_matrix(ii,jj) = initialTangent[i][j][k][l] ;

    } //end for j
  } //end for i

  return tangent_matrix ;
} 

//swap history variables
int J2AxiSymm :: commitState( )  
{
  epsilon_p_n = epsilon_p_nplus1 ;
  xi_n        = xi_nplus1 ;

  return 0 ;
}


//revert to last saved state
int J2AxiSymm :: revertToLastCommit( )
{ 
  return 0 ;
} 

//revert to start
int J2AxiSymm :: revertToStart( ) 

{  
  this->zero( ) ;
  return 0 ;
}

int
J2AxiSymm::sendSelf (int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(10+9);
  int cnt = 0;
  data(cnt++) = this->getTag();
  data(cnt++) = bulk;
  data(cnt++) = shear;
  data(cnt++) = sigma_0;
  data(cnt++) = sigma_infty;
  data(cnt++) = delta;
  data(cnt++) = Hard;
  data(cnt++) = eta;
  data(cnt++) = rho;

  data(cnt++) = xi_n;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      data(cnt++) = epsilon_p_n(i,j);

  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "J2AxiSymm::recvSelf - failed to send vector to channel\n";
    return -1;
  }

  return 0;
}

int
J2AxiSymm::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{

  // recv the vector object from the channel which defines material param and state
  static Vector data(10+9);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "J2AxiSymm::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag(data(cnt++));
  bulk = data(cnt++);
  shear = data(cnt++);
  sigma_0 = data(cnt++);
  sigma_infty = data(cnt++);
  delta = data(cnt++);
  Hard = data(cnt++);
  eta = data(cnt++);
  rho = data(cnt++);

  xi_n = data(cnt++);

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++) 
      epsilon_p_n(i,j) = data(cnt++);

  epsilon_p_nplus1 = epsilon_p_n;
  xi_nplus1        = xi_n;

  return 0;
}












