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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2PlaneStress.cpp,v $

// Written: Ed "C++" Love
//
// J2PlaneStress isotropic hardening material class
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

#include <J2PlaneStress.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector J2PlaneStress :: strain_vec(3) ;
Vector J2PlaneStress :: stress_vec(3) ;
Matrix J2PlaneStress :: tangent_matrix(3,3) ;

//null constructor
J2PlaneStress ::  J2PlaneStress( ) : 
J2Plasticity( ) 
{  }


//full constructor
J2PlaneStress :: 
J2PlaneStress(   int    tag, 
                 double K,
                 double G,
                 double yield0,
                 double yield_infty,
                 double d,
                 double H,
                 double viscosity,
		 double rho) : 
J2Plasticity(tag, ND_TAG_J2PlaneStress, 
             K, G, yield0, yield_infty, d, H, viscosity, rho)
{ 

}


//elastic constructor
J2PlaneStress :: 
J2PlaneStress(   int    tag, 
                 double K, 
                 double G ) :
J2Plasticity(tag, ND_TAG_J2PlaneStress, K, G )
{ 

}


//destructor
J2PlaneStress :: ~J2PlaneStress( ) 
{  } 


//make a clone of this material
NDMaterial* J2PlaneStress :: getCopy( ) 
{ 
  J2PlaneStress  *clone;
  clone = new J2PlaneStress( ) ;   //new instance of this class
  *clone = *this ;                 //asignment to make copy
  return clone ;
}


//send back type of material
const char* J2PlaneStress :: getType( ) const 
{
  return "PlaneStress" ;
}


//send back order of strain in vector form
int J2PlaneStress :: getOrder( ) const 
{ 
  return 3 ; 
} 

//get the strain and integrate plasticity equations
int J2PlaneStress :: setTrialStrain( const Vector &strain_from_element ) 
{
  const double tolerance = (1.0e-8)*sigma_0 ;

  const int max_iterations = 25 ;
  int iteration_counter  = 0 ;

  int i, j, k, l ;
  int ii, jj ;

  double eps22  =  strain(2,2) ;
  strain.Zero( ) ;

  strain(0,0) =        strain_from_element(0) ;
  strain(1,1) =        strain_from_element(1) ;
  strain(0,1) = 0.50 * strain_from_element(2) ;
  strain(1,0) =        strain(0,1) ;

  strain(2,2) =        eps22 ; 

  //enforce the plane stress condition sigma_22 = 0 
  //solve for epsilon_22 
  iteration_counter = 0 ;  
  do {

     this->plastic_integrator( ) ;
    
     strain(2,2) -= stress(2,2) / tangent[2][2][2][2] ;

     iteration_counter++ ;
     if ( iteration_counter > max_iterations ) {
       opserr << "More than " << max_iterations ;
       opserr << " iterations in setTrialStrain of J2PlaneStress \n" ;
       break ;
     }// end if 

  } while ( fabs(stress(2,2)) > tolerance ) ;

  //modify tangent for plane stress 
  for ( ii = 0; ii < 3; ii++ ) {
    for ( jj = 0; jj < 3; jj++ )  {

          index_map( ii, i, j ) ;
          index_map( jj, k, l ) ;

          tangent[i][j][k][l] -=   tangent[i][j][2][2] 
                                 * tangent[2][2][k][l] 
                                 / tangent[2][2][2][2] ;

          //minor symmetries 
          tangent [j][i][k][l] = tangent[i][j][k][l] ;
          tangent [i][j][l][k] = tangent[i][j][k][l] ;
          tangent [j][i][l][k] = tangent[i][j][k][l] ;

    } // end for jj
  } // end for ii 

  return 0 ;
}


//unused trial strain functions
int J2PlaneStress :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

int J2PlaneStress :: setTrialStrainIncr( const Vector &v ) 
{
  static Vector newStrain(3);
  newStrain(0) = strain(0,0) + v(0);
  newStrain(1) = strain(1,1) + v(1);
  newStrain(2) = 2.0 * strain(0,1) + v(2);

  return this->setTrialStrain(newStrain);  
}

int J2PlaneStress :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
    return this->setTrialStrainIncr(v);
}


//send back the strain
const Vector& J2PlaneStress :: getStrain( ) 
{
  strain_vec(0) =       strain(0,0) ;
  strain_vec(1) =       strain(1,1) ;
  strain_vec(2) = 2.0 * strain(0,1) ;

  return strain_vec ;
} 


//send back the stress 
const Vector& J2PlaneStress :: getStress( ) 
{
  stress_vec(0) = stress(0,0) ;
  stress_vec(1) = stress(1,1) ;
  stress_vec(2) = stress(0,1) ;

  return stress_vec ;
}

//send back the tangent 
const Matrix& J2PlaneStress :: getTangent( ) 
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
const Matrix& J2PlaneStress :: getInitialTangent( ) 
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

int 
J2PlaneStress::commitState( ) 
{
  epsilon_p_n = epsilon_p_nplus1 ;
  xi_n        = xi_nplus1 ;

  commitEps22 =       strain(2,2);

  return 0;
}

int 
J2PlaneStress::revertToLastCommit( ) {

  strain(2,2)  = commitEps22;

  return 0;
}


int 
J2PlaneStress::revertToStart( ) {
  commitEps22 = 0.0;

  this->zero( ) ;

  return 0;
}

int
J2PlaneStress::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(11+9);
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
  data(cnt++) = commitEps22;

  for (int i=0; i<3; i++) 
    for (int j=0; j<3; j++) 
      data(cnt++) = epsilon_p_n(i,j);

  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "J2PlaneStress::sendSelf - failed to send vector to channel\n";
    return -1;
  }

  return 0;
}

int
J2PlaneStress::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{

  // recv the vector object from the channel which defines material param and state
  static Vector data(11+9);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "J2PlaneStress::recvSelf - failed to recv vector from channel\n";
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
  commitEps22 = data(cnt++);

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++) 
      epsilon_p_n(i,j) = data(cnt++);

  epsilon_p_nplus1 = epsilon_p_n;
  xi_nplus1        = xi_n;
  strain(2,2) = commitEps22;

  return 0;
}


//matrix_index ---> tensor indices i,j
// plane stress different because of condensation on tangent
// case 3 switched to 1-2 and case 4 to 3-3 
void 
J2PlaneStress :: index_map( int matrix_index, int &i, int &j )
{
  switch ( matrix_index+1 ) { //add 1 for standard tensor indices

    case 1 :
      i = 1 ; 
      j = 1 ;
      break ;
 
    case 2 :
      i = 2 ;
      j = 2 ; 
      break ;

    case 3 :
      i = 1 ;
      j = 2 ;
      break ;

    case 4 :
      i = 3 ;
      j = 3 ;
      break ;

    case 5 :
      i = 2 ;
      j = 3 ;
      break ;

    case 6 :
      i = 3 ;
      j = 1 ;
      break ;


    default :
      i = 1 ;
      j = 1 ;
      break ;

  } //end switch

i-- ; //subtract 1 for C-indexing
j-- ;

return ; 
}


