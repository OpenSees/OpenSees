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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticOrthotropicPlaneStress.cpp,v $

// Written: Ed "C++" Love
//
// ElasticOrthotropicPlaneStress isotropic hardening material class
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

#include <ElasticOrthotropicPlaneStress.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

// Vector ElasticOrthotropicPlaneStress :: strain_vec(3) ;
Vector ElasticOrthotropicPlaneStress :: stress_vec(3) ;
Matrix ElasticOrthotropicPlaneStress :: tangent_matrix(3,3) ;


void* OPS_ElasticOrthotropicPlaneStress(void) 
{
  
  opserr << "ndMaterial ElasticOrthotropicPlaneStress tag E1, E2, nu12, nu21, G12,  rho\n";

  int tag;
  double E1, E2, nu12, nu21, G12, rho;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 7) {
    opserr << "ndMaterial ElasticOrthotropicPlaneStress tag E1, E2, nu12, nu21, G12,  rho\n";
    return 0;
  }

  int iData[1];
  double dData[6];

  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer values: nDMaterial ElasticOrthotropicPlaneStress \n";
    return 0;
  }  
  tag = iData[0]; 

  numData = 6;  
  if (OPS_GetDouble(&numData, dData) != 0) {
      opserr << "WARNING invalid double values: nDMaterial ElasticOrthotropicPlaneStress " << tag << endln;
    return 0;
  }  
  E1 = dData[0];
  E2 = dData[1];
  nu12 = dData[2];
  nu21 = dData[3];
  G12 = dData[4];
  rho = dData[5];

  opserr << "Creating new ElasticOrthotropicPlaneStress with \n" 
    << "tag  = " << tag << endln
    << "E1    = " << E1 << endln
    << "E2    = " << E2 << endln
    << "nu12   = " << nu12 << endln
    << "nu21   = " << nu21 << endln
    << "G12   = " << G12 << endln
    << "rho  = " << rho << endln;



  NDMaterial *theMaterial =  new ElasticOrthotropicPlaneStress (tag, E1, E2, nu12, nu21, G12, rho);
  
  return theMaterial;
}






//null constructor
ElasticOrthotropicPlaneStress ::  ElasticOrthotropicPlaneStress( ) : 
NDMaterial(0, ND_TAG_ElasticOrthotropicPlaneStress), 
strain_vec(3),
E1(0),
E2(0),
nu12(0),
nu21(0),
G12(0),
rho(0)
{  }


//full constructor
ElasticOrthotropicPlaneStress :: 
ElasticOrthotropicPlaneStress(int tag, 
                   double E1_,
                   double E2_,
                   double nu12_,
                   double nu21_,
                   double G12_,
                   double rho_) : 
NDMaterial(tag, ND_TAG_ElasticOrthotropicPlaneStress), 
strain_vec(3),
E1(E1_),
E2(E2_),
nu12(nu12_),
nu21(nu21_),
G12(G12_),
rho(rho_)
{ 

}


//destructor
ElasticOrthotropicPlaneStress :: ~ElasticOrthotropicPlaneStress( ) 
{  } 


//make a clone of this material
NDMaterial* ElasticOrthotropicPlaneStress :: getCopy( ) 
{ 
  ElasticOrthotropicPlaneStress  *clone;
  clone = new ElasticOrthotropicPlaneStress( ) ;   //new instance of this class
  *clone = *this ;                 //asignment to make copy
  return clone ;
}


//send back type of material
const char* ElasticOrthotropicPlaneStress :: getType( ) const 
{
  return "OrthotropicPlaneStress" ;
}


//send back order of strain in vector form
int ElasticOrthotropicPlaneStress :: getOrder( ) const 
{ 
  return 3 ; 
} 

//mass per unit volume
double
ElasticOrthotropicPlaneStress::getRho( )
{
  return rho ;
}

//get the strain and integrate plasticity equations
int ElasticOrthotropicPlaneStress :: setTrialStrain( const Vector &strain_from_element ) 
{
  strain_vec = strain_from_element;
  
  return 0 ;
}


//unused trial strain functions
int ElasticOrthotropicPlaneStress :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   opserr << "ElasticOrthotropicPlaneStress :: setTrialStrain( const Vector &v, const Vector &r ) -- should not be used! \n";
   return this->setTrialStrain( v ) ;
} 

int ElasticOrthotropicPlaneStress :: setTrialStrainIncr( const Vector &v ) 
{
   opserr << "ElasticOrthotropicPlaneStress :: setTrialStrainIncr( const Vector &v ) -- should not be used! \n";
   return -1 ;
}

int ElasticOrthotropicPlaneStress :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
   opserr << "ElasticOrthotropicPlaneStress :: setTrialStrainIncr( const Vector &v, const Vector &r ) -- should not be used! \n";
    return this->setTrialStrainIncr(v);
}


//send back the strain
const Vector& ElasticOrthotropicPlaneStress :: getStrain( ) 
{
  return strain_vec ;
} 


//send back the stress 
const Vector& ElasticOrthotropicPlaneStress :: getStress( ) 
{

  stress_vec = this->getTangent() * strain_vec;

  return stress_vec ;
}

//send back the tangent 
const Matrix& ElasticOrthotropicPlaneStress :: getTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 ) 
  // 
       
  double den = 1 - nu12*nu21;
  // double G12 = E1 / (2*(1+nu21));

  tangent_matrix(0,0) = E1 / den ;
  tangent_matrix(1,1) = E2 / den ;
  tangent_matrix(2,2) = G12 ;

  tangent_matrix(0,1) = nu21 * E1 / den ;
  tangent_matrix(1,0) = nu12 * E2 / den ;

  tangent_matrix(0,2) = 0. ;
  tangent_matrix(2,0) = 0. ;

  tangent_matrix(1,2) = 0. ;
  tangent_matrix(2,1) = 0. ;

  return tangent_matrix ;
} 


//send back the tangent 
const Matrix& ElasticOrthotropicPlaneStress :: getInitialTangent( ) 
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
ElasticOrthotropicPlaneStress::commitState( ) 
{

  return 0;
}

int 
ElasticOrthotropicPlaneStress::revertToLastCommit( ) {


  return 0;
}


int 
ElasticOrthotropicPlaneStress::revertToStart( ) {
 

  return 0;
}

int
ElasticOrthotropicPlaneStress::sendSelf(int commitTag, Channel &theChannel)
{
  

  return -1;
}

int
ElasticOrthotropicPlaneStress::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  

  return -1;
}


//print out material data
void ElasticOrthotropicPlaneStress :: Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "ElasticOrthotropicPlaneStress : " ; 
  s << this->getType( ) << endln ;
  s << "Elastic Modulus 1 =   " << E1        << endln ;
  s << "Elastic Modulus 2 =   " << E2        << endln ;
  s << "Poisson's ratio 12=  " << nu12       << endln ;
  s << "Poisson's ratio 21=  " << nu21       << endln ;
  s << "Shear constant G12=  " << G12       << endln ;
  s << "mass density =        " << rho     << endln ;
  s << endln ;
}


//matrix_index ---> tensor indices i,j
// plane stress different because of condensation on tangent
// case 3 switched to 1-2 and case 4 to 3-3 
void 
ElasticOrthotropicPlaneStress :: index_map( int matrix_index, int &i, int &j )
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


