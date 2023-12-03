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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticPlaneStress.cpp,v $

// Written: Ed "C++" Love
//
// ElasticPlaneStress isotropic hardening material class
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

#include <ElasticPlaneStress.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

// Vector ElasticPlaneStress :: strain_vec(3) ;
Vector ElasticPlaneStress :: stress_vec(3) ;
Matrix ElasticPlaneStress :: tangent_matrix(3,3) ;


void* OPS_ElasticPlaneStress(void) 
{
  
  opserr << "ndMaterial ElasticPlaneStress tag E nu rho\n";

  int tag;
  double E, nu, rho;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 4) {
    opserr << "ndMaterial ElasticPlaneStress tag E nu rho\n";
    return 0;
  }

  int iData[1];
  double dData[3];

  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer values: nDMaterial ElasticPlaneStress \n";
    return 0;
  }  
  tag = iData[0]; 

  numData = 3;  
  if (OPS_GetDouble(&numData, dData) != 0) {
      opserr << "WARNING invalid double values: nDMaterial ElasticPlaneStress " << tag << endln;
    return 0;
  }  
  E = dData[0];
  nu = dData[1];
  rho = dData[2];

  opserr << "Creating new ElasticPlaneStress with \n" 
    << "tag  = " << tag << endln
    << "E    = " << E << endln
    << "nu   = " << nu << endln
    << "rho  = " << rho << endln;



  NDMaterial *theMaterial =  new ElasticPlaneStress (tag, E, nu, rho);
  
  return theMaterial;
}






//null constructor
ElasticPlaneStress ::  ElasticPlaneStress( ) : 
NDMaterial(0, ND_TAG_ElasticPlaneStress), 
strain_vec(3),
E(0),
nu(0),
rho(0)
{  }


//full constructor
ElasticPlaneStress :: 
ElasticPlaneStress(int tag, 
                   double E_,
                   double nu_,
                   double rho_) : 
NDMaterial(tag, ND_TAG_ElasticPlaneStress), 
strain_vec(3),
E(E_),
nu(nu_),
rho(rho_)
{ 

}


//destructor
ElasticPlaneStress :: ~ElasticPlaneStress( ) 
{  } 


//make a clone of this material
NDMaterial* ElasticPlaneStress :: getCopy( ) 
{ 
  ElasticPlaneStress  *clone;
  clone = new ElasticPlaneStress( ) ;   //new instance of this class
  *clone = *this ;                 //asignment to make copy
  return clone ;
}


//send back type of material
const char* ElasticPlaneStress :: getType( ) const 
{
  return "PlaneStress" ;
}


//send back order of strain in vector form
int ElasticPlaneStress :: getOrder( ) const 
{ 
  return 3 ; 
} 

//mass per unit volume
double
ElasticPlaneStress::getRho( )
{
  return rho ;
}

//get the strain and integrate plasticity equations
int ElasticPlaneStress :: setTrialStrain( const Vector &strain_from_element ) 
{
  strain_vec = strain_from_element;
  
  return 0 ;
}


//unused trial strain functions
int ElasticPlaneStress :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   opserr << "ElasticPlaneStress :: setTrialStrain( const Vector &v, const Vector &r ) -- should not be used! \n";
   return this->setTrialStrain( v ) ;
} 

int ElasticPlaneStress :: setTrialStrainIncr( const Vector &v ) 
{
   opserr << "ElasticPlaneStress :: setTrialStrainIncr( const Vector &v ) -- should not be used! \n";
   return -1 ;
}

int ElasticPlaneStress :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
   opserr << "ElasticPlaneStress :: setTrialStrainIncr( const Vector &v, const Vector &r ) -- should not be used! \n";
    return this->setTrialStrainIncr(v);
}


//send back the strain
const Vector& ElasticPlaneStress :: getStrain( ) 
{
  return strain_vec ;
} 


//send back the stress 
const Vector& ElasticPlaneStress :: getStress( ) 
{

  //stress_vec = this->getTangent() * strain_vec;
  double den = 1 - nu*nu;
  double G = E / (2*(1+nu));

  stress_vec(0) = E/den*(strain_vec(0) + nu*strain_vec(1));
  stress_vec(1) = E/den*(strain_vec(1) + nu*strain_vec(0));
  stress_vec(2) = G*strain_vec(2);
  
  return stress_vec ;
}

//send back the tangent 
const Matrix& ElasticPlaneStress :: getTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 ) 
  // 
       
  double den = 1 - nu*nu;
  double G = E / (2*(1+nu));

  tangent_matrix(0,0) = E / den ;
  tangent_matrix(1,1) = E / den ;
  tangent_matrix(2,2) = G ;

  tangent_matrix(0,1) = nu * E / den ;
  tangent_matrix(1,0) = nu * E / den ;

  tangent_matrix(0,2) = 0. ;
  tangent_matrix(2,0) = 0. ;

  tangent_matrix(1,2) = 0. ;
  tangent_matrix(2,1) = 0. ;

  return tangent_matrix ;
} 


//send back the tangent 
const Matrix& ElasticPlaneStress :: getInitialTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 ) 
  // 


  return this->getTangent() ;
} 

int 
ElasticPlaneStress::commitState( ) 
{

  return 0;
}

int 
ElasticPlaneStress::revertToLastCommit( ) {


  return 0;
}


int 
ElasticPlaneStress::revertToStart( ) {
 

  return 0;
}

int
ElasticPlaneStress::sendSelf(int commitTag, Channel &theChannel)
{
  

  return -1;
}

int
ElasticPlaneStress::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  

  return -1;
}


//print out material data
void ElasticPlaneStress :: Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "ElasticPlaneStress : " ; 
  s << this->getType( ) << endln ;
  s << "Elastic Modulus =   " << E        << endln ;
  s << "Poisson's ratio =  " << nu       << endln ;
  s << "mass density =        " << rho     << endln ;
  s << endln ;
}


//matrix_index ---> tensor indices i,j
// plane stress different because of condensation on tangent
// case 3 switched to 1-2 and case 4 to 3-3 
void 
ElasticPlaneStress :: index_map( int matrix_index, int &i, int &j )
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


