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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2ThreeDimensionalThermal.cpp,v $

// Written: Ed "C++" Love
// Do not ask Prashant about this code.  He has no clue. 
//
// J2ThreeDimensionalThermal isotropic hardening material class
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
//                    2 eps_01   
//            	      2 eps_12   
//		      2 eps_20    }   <--- note the 2
// 
//  set eta := 0 for rate independent case
//Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io]

#include <J2ThreeDimensionalThermal.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static vectors and matrices
Vector J2ThreeDimensionalThermal :: strain_vec(6) ;
Vector J2ThreeDimensionalThermal :: stress_vec(6) ;
Matrix J2ThreeDimensionalThermal :: tangent_matrix(6,6) ;


//null constructor
J2ThreeDimensionalThermal ::  J2ThreeDimensionalThermal( ) : 
	J2PlasticityThermal( )
{ }


//full constructor
J2ThreeDimensionalThermal :: 
J2ThreeDimensionalThermal(   int    tag, 
		      double K,
		      double G,
		      double yield0,
		      double yield_infty,
		      double d,
		      double H,
		      double viscosity,
		      double rho) : 
  J2PlasticityThermal( tag, ND_TAG_J2ThreeDimensionalThermal, 
		K, G, yield0, yield_infty, d, H, viscosity, rho)
{ 
}


//elastic constructor
J2ThreeDimensionalThermal :: 
J2ThreeDimensionalThermal(   int    tag, 
                 double K, 
                 double G ) :
J2PlasticityThermal( tag, ND_TAG_J2ThreeDimensionalThermal, K, G )
{ 
}



//destructor
J2ThreeDimensionalThermal :: ~J2ThreeDimensionalThermal( ) 
{} 


//make a clone of this material
NDMaterial* J2ThreeDimensionalThermal :: getCopy( ) 
{ 
  J2ThreeDimensionalThermal  *clone;
  clone = new J2ThreeDimensionalThermal( ) ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}


//send back type of material
const char* J2ThreeDimensionalThermal :: getType( ) const 
{
  return "ThreeDimensionalThermal" ;
}


//send back order of strain in vector form
int J2ThreeDimensionalThermal :: getOrder( ) const 
{ 
  return 6 ; 
} 


//get the strain and integrate plasticity equations
int J2ThreeDimensionalThermal :: setTrialStrain( const Vector &strain_from_element) 
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

  this->plastic_integrator( ) ;

  return 0 ;
}


//unused trial strain functions
int J2ThreeDimensionalThermal :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

int J2ThreeDimensionalThermal :: setTrialStrainIncr( const Vector &v ) 
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

int J2ThreeDimensionalThermal :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
  return this->setTrialStrainIncr(v);
}



//send back the strain
const Vector& J2ThreeDimensionalThermal :: getStrain( ) 
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
const Vector& J2ThreeDimensionalThermal :: getStress( ) 
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
const Matrix& J2ThreeDimensionalThermal :: getTangent( ) 
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
const Matrix& J2ThreeDimensionalThermal :: getInitialTangent( ) 
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

//send back TempAndElong(Liming)
const Vector& J2ThreeDimensionalThermal::getTempAndElong() {

	return TempAndElong;
}
