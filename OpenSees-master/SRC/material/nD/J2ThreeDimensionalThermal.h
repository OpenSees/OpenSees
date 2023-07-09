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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2ThreeDimensionalThermal.h,v $

// Written: Ed "C++" Love

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

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>

#include <J2PlasticityThermal.h>

class J2ThreeDimensionalThermal : public J2PlasticityThermal {

//-------------------Declarations-------------------------------

  public : 

  //null constructor
  J2ThreeDimensionalThermal( ) ;

  //full constructor
  J2ThreeDimensionalThermal(   int    tag, 
                   double K,
                   double G,
                   double yield0,
                   double yield_infty,
                   double d,
                   double H,
			double viscosity = 0,
			double rho = 0) ;


  //elastic constructor
  J2ThreeDimensionalThermal( int tag, double K, double G ) ;

  //destructor
  ~J2ThreeDimensionalThermal( ) ;

  const char *getClassType(void) const {return "J2ThreeDimensionalThermal";};

  //make a clone of this material
  NDMaterial* getCopy( ) ;

  //send back type of material
  const char* getType( ) const ;

  //send back order of strain in vector form
  int getOrder( ) const ;

  //get the strain and integrate plasticity equations
  int setTrialStrain( const Vector &strain_from_element) ;

  //unused trial strain functions
  int setTrialStrain( const Vector &v, const Vector &r ) ;
  int setTrialStrainIncr( const Vector &v ) ;
  int setTrialStrainIncr( const Vector &v, const Vector &r ) ;

  //send back the strain
  const Vector& getStrain( ) ;

  //send back the stress 
  const Vector& getStress( ) ;

  //send back the tangent 
  const Matrix& getTangent( ) ;
  const Matrix& getInitialTangent( ) ;

  //added by Limin
  const Vector& getTempAndElong();

  private :

  //static vectors and matrices
  static Vector strain_vec ;     //strain in vector notation
  static Vector stress_vec ;     //stress in vector notation
  static Matrix tangent_matrix ; //material tangent in matrix notation

} ; //end of J2ThreeDimensionalThermal declarations


