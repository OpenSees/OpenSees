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
                                                                        
// $Revision: 1.12 $
// $Date: 2008-10-20 22:23:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2PlasticityThermal.cpp,v $

// Written: Ed "C++" Love
//
// J2 isotropic hardening material class
// 
//  Elastic Model
//  sigma = K*trace(epsilion_elastic) + (2*G)*dev(epsilon_elastic)
//
//  Yield Function
//  phi(sigma,q) = || dev(sigma) ||  - sqrt(2/3)*q(xi) 
//
//  Saturation Isotropic Hardening with linear term
//  q(xi) = simga_infty + (sigma_y - sigma_infty)*exp(-delta*xi) + H*xi 
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
//  set eta := 0 for rate independent case

//Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io]
//

#include <J2PlasticityThermal.h>
#include <J2PlaneStress.h>
#include <J2PlaneStrain.h>
#include <J2AxiSymm.h>
#include <J2PlateFiber.h>

#include <J2ThreeDimensional.h> 
#include <J2ThreeDimensionalThermal.h> 
#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//parameters
const double J2PlasticityThermal :: one3   = 1.0 / 3.0 ;
const double J2PlasticityThermal :: two3   = 2.0 / 3.0 ;
const double J2PlasticityThermal :: four3  = 4.0 / 3.0 ;
const double J2PlasticityThermal :: root23 = sqrt( 2.0 / 3.0 ) ;

double J2PlasticityThermal::initialTangent[3][3][3][3] ;   //material tangent
double J2PlasticityThermal::IIdev[3][3][3][3] ; //rank 4 deviatoric 
double J2PlasticityThermal::IbunI[3][3][3][3] ; //rank 4 I bun I 

//zero internal variables
void J2PlasticityThermal :: zero ( ) 
{
  xi_n = 0.0 ;
  xi_nplus1 = 0.0 ;
  
  epsilon_p_n.Zero( ) ;
  epsilon_p_nplus1.Zero( ) ;

  stress.Zero();
  strain.Zero();
}


//null constructor
J2PlasticityThermal ::  J2PlasticityThermal( ) : 
NDMaterial( ),
epsilon_p_n(3,3),
epsilon_p_nplus1(3,3),
stress(3,3),
strain(3,3),
TempAndElong(2)
{ 
  bulk        = 0.0 ;
  shear       = 0.0 ;
  sigma_y     = 0.0 ;
  bulk_0 = 0.0;
  shear_0 = 0.0;
  sigma_0 = 0.0;

  sigma_infty = 0.0 ;
  delta       = 0.0 ;
  Hard        = 0.0 ;
  eta         = 0.0 ;
  rho = 0.0;

  this->zero( ) ;     // or (*this).zero( ) 

  int i, j, k, l ;

  //zero rank4 IIdev and IbunI 
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ )  {
      for ( k = 0; k < 3; k++ ) {
	for ( l = 0; l < 3; l++)  { 

	  IbunI[i][j][k][l] = 0.0 ;

	  IIdev[i][j][k][l] = 0.0 ;

	} // end for l
      } // end for k
    } // end for j
  } // end for i


  //form rank4 IbunI 

  IbunI [0][0] [0][0] = 1.0 ;
  IbunI [0][0] [1][1] = 1.0 ;
  IbunI [0][0] [2][2] = 1.0 ;
  IbunI [1][1] [0][0] = 1.0 ;
  IbunI [1][1] [1][1] = 1.0 ;
  IbunI [1][1] [2][2] = 1.0 ;
  IbunI [2][2] [0][0] = 1.0 ;
  IbunI [2][2] [1][1] = 1.0 ;
  IbunI [2][2] [2][2] = 1.0 ;

  //form rank4 IIdev

  IIdev [0][0] [0][0] =  two3 ; // 0.666667 
  IIdev [0][0] [1][1] = -one3 ; //-0.333333 
  IIdev [0][0] [2][2] = -one3 ; //-0.333333 
  IIdev [0][1] [0][1] = 0.5 ;
  IIdev [0][1] [1][0] = 0.5 ;
  IIdev [0][2] [0][2] = 0.5 ;
  IIdev [0][2] [2][0] = 0.5 ;
  IIdev [1][0] [0][1] = 0.5 ;
  IIdev [1][0] [1][0] = 0.5 ;
  IIdev [1][1] [0][0] = -one3 ; //-0.333333 
  IIdev [1][1] [1][1] =  two3 ; // 0.666667 
  IIdev [1][1] [2][2] = -one3 ; //-0.333333 
  IIdev [1][2] [1][2] = 0.5 ;
  IIdev [1][2] [2][1] = 0.5 ;
  IIdev [2][0] [0][2] = 0.5 ;
  IIdev [2][0] [2][0] = 0.5 ;
  IIdev [2][1] [1][2] = 0.5 ;
  IIdev [2][1] [2][1] = 0.5 ;
  IIdev [2][2] [0][0] = -one3 ; //-0.333333 
  IIdev [2][2] [1][1] = -one3 ; //-0.333333 
  IIdev [2][2] [2][2] =  two3 ; // 0.666667 

  ThermalElongation = 0.0;
  plastic_integrator();
}


//full constructor
J2PlasticityThermal :: J2PlasticityThermal(int    tag,
			     int classTag,
			     double K,
			     double G,
			     double yield0,
			     double yield_infty,
			     double d,
			     double H,
			     double viscosity,
			     double r) 
: 
  NDMaterial(tag, classTag),
  epsilon_p_n(3,3),
  epsilon_p_nplus1(3,3),
  stress(3,3),
  strain(3,3),
	TempAndElong(2)
{
  bulk        = K ;
  shear       = G ;
  sigma_y     = yield0 ;
  bulk_0 = K;
  shear_0 = G;
  sigma_0 = yield0;

  sigma_infty = yield_infty ;
  delta       = d ;
  Hard        = H ;
  eta         = viscosity ;
  rho = r;

  this->zero( ) ;

  int i, j, k, l ;

  //zero rank4 IIdev and IbunI 
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ )  {
      for ( k = 0; k < 3; k++ ) {
	for ( l = 0; l < 3; l++)  { 

	  IbunI[i][j][k][l] = 0.0 ;

	  IIdev[i][j][k][l] = 0.0 ;

	} // end for l
      } // end for k
    } // end for j
  } // end for i


  //form rank4 IbunI 

  IbunI [0][0] [0][0] = 1.0 ;
  IbunI [0][0] [1][1] = 1.0 ;
  IbunI [0][0] [2][2] = 1.0 ;
  IbunI [1][1] [0][0] = 1.0 ;
  IbunI [1][1] [1][1] = 1.0 ;
  IbunI [1][1] [2][2] = 1.0 ;
  IbunI [2][2] [0][0] = 1.0 ;
  IbunI [2][2] [1][1] = 1.0 ;
  IbunI [2][2] [2][2] = 1.0 ;

  //form rank4 IIdev

  IIdev [0][0] [0][0] =  two3 ; // 0.666667 
  IIdev [0][0] [1][1] = -one3 ; //-0.333333 
  IIdev [0][0] [2][2] = -one3 ; //-0.333333 
  IIdev [0][1] [0][1] = 0.5 ;
  IIdev [0][1] [1][0] = 0.5 ;
  IIdev [0][2] [0][2] = 0.5 ;
  IIdev [0][2] [2][0] = 0.5 ;
  IIdev [1][0] [0][1] = 0.5 ;
  IIdev [1][0] [1][0] = 0.5 ;
  IIdev [1][1] [0][0] = -one3 ; //-0.333333 
  IIdev [1][1] [1][1] =  two3 ; // 0.666667 
  IIdev [1][1] [2][2] = -one3 ; //-0.333333 
  IIdev [1][2] [1][2] = 0.5 ;
  IIdev [1][2] [2][1] = 0.5 ;
  IIdev [2][0] [0][2] = 0.5 ;
  IIdev [2][0] [2][0] = 0.5 ;
  IIdev [2][1] [1][2] = 0.5 ;
  IIdev [2][1] [2][1] = 0.5 ;
  IIdev [2][2] [0][0] = -one3 ; //-0.333333 
  IIdev [2][2] [1][1] = -one3 ; //-0.333333 
  IIdev [2][2] [2][2] =  two3 ; // 0.666667 

  ThermalElongation = 0;
  plastic_integrator();
}


//elastic constructor
J2PlasticityThermal :: 
J2PlasticityThermal(   int    tag, 
                int  classTag,
                double K, 
                double G ) :
NDMaterial(tag, classTag),
epsilon_p_n(3,3),
epsilon_p_nplus1(3,3),
stress(3,3),
strain(3,3),
TempAndElong(2)
{
  bulk        = K ;
  shear       = G ; 
  bulk_0 = K;
  shear_0 = G;
  sigma_y     = 1.0e16*shear ;
  sigma_0 = 1.0e16*shear;
  sigma_infty = sigma_y ;
  delta       = 0.0 ;
  Hard        = 0.0 ;
  eta         = 0.0 ;

  this->zero( ) ;

  int i, j, k, l ;

  //zero rank4 IIdev and IbunI 
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ )  {
      for ( k = 0; k < 3; k++ ) {
	for ( l = 0; l < 3; l++)  { 

	  IbunI[i][j][k][l] = 0.0 ;

	  IIdev[i][j][k][l] = 0.0 ;

	} // end for l
      } // end for k
    } // end for j
  } // end for i


  //form rank4 IbunI 

  IbunI [0][0] [0][0] = 1.0 ;
  IbunI [0][0] [1][1] = 1.0 ;
  IbunI [0][0] [2][2] = 1.0 ;
  IbunI [1][1] [0][0] = 1.0 ;
  IbunI [1][1] [1][1] = 1.0 ;
  IbunI [1][1] [2][2] = 1.0 ;
  IbunI [2][2] [0][0] = 1.0 ;
  IbunI [2][2] [1][1] = 1.0 ;
  IbunI [2][2] [2][2] = 1.0 ;

  //form rank4 IIdev

  IIdev [0][0] [0][0] =  two3 ; // 0.666667 
  IIdev [0][0] [1][1] = -one3 ; //-0.333333 
  IIdev [0][0] [2][2] = -one3 ; //-0.333333 
  IIdev [0][1] [0][1] = 0.5 ;
  IIdev [0][1] [1][0] = 0.5 ;
  IIdev [0][2] [0][2] = 0.5 ;
  IIdev [0][2] [2][0] = 0.5 ;
  IIdev [1][0] [0][1] = 0.5 ;
  IIdev [1][0] [1][0] = 0.5 ;
  IIdev [1][1] [0][0] = -one3 ; //-0.333333 
  IIdev [1][1] [1][1] =  two3 ; // 0.666667 
  IIdev [1][1] [2][2] = -one3 ; //-0.333333 
  IIdev [1][2] [1][2] = 0.5 ;
  IIdev [1][2] [2][1] = 0.5 ;
  IIdev [2][0] [0][2] = 0.5 ;
  IIdev [2][0] [2][0] = 0.5 ;
  IIdev [2][1] [1][2] = 0.5 ;
  IIdev [2][1] [2][1] = 0.5 ;
  IIdev [2][2] [0][0] = -one3 ; //-0.333333 
  IIdev [2][2] [1][1] = -one3 ; //-0.333333 
  IIdev [2][2] [2][2] =  two3 ; // 0.666667 

  ThermalElongation = 0.0;
}


//destructor
J2PlasticityThermal :: ~J2PlasticityThermal( ) 
{  } 



NDMaterial*
J2PlasticityThermal :: getCopy (const char *type)
{
    if (strcmp(type,"PlaneStress2D") == 0 || strcmp(type,"PlaneStress") == 0)
    {
	J2PlaneStress  *clone ;
	clone = new J2PlaneStress(this->getTag(), bulk, shear, sigma_y,
				  sigma_infty, delta, Hard, eta, rho) ;
	return clone ;
    }
    else if ((strcmp(type,"ThreeDimensional") == 0) ||
	     (strcmp(type,"3D") == 0))
    {
	J2ThreeDimensional  *clone ;
	clone = new J2ThreeDimensional(this->getTag(), bulk, shear, sigma_y,
				       sigma_infty, delta, Hard, eta, rho) ;
	return clone ;	
    }
	else if ((strcmp(type, "ThreeDimensionalThermal") == 0) ||
		(strcmp(type, "3DThermal") == 0))
	{
		J2ThreeDimensionalThermal  *clone;
		clone = new J2ThreeDimensionalThermal(this->getTag(), bulk, shear, sigma_y,
			sigma_infty, delta, Hard, eta, rho);
		return clone;
	}
    // Handle other cases
    else
    {
      return NDMaterial::getCopy(type);
    }
}

//print out material data
void J2PlasticityThermal :: Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "J2-Plasticity : " ; 
  s << this->getType( ) << endln ;
  s << "Bulk Modulus =   " << bulk        << endln ;
  s << "Shear Modulus =  " << shear       << endln ;
  s << "sigma_y =        " << sigma_y     << endln ;
  s << "Sigma_infty =    " << sigma_infty << endln ;
  s << "Delta =          " << delta       << endln ;
  s << "H =              " << Hard        << endln ;
  s << "Eta =            " << eta         << endln ;
  s << "Rho =            " << rho         << endln ;
  s << endln ;
}


//--------------------Plasticity-------------------------------------

//plasticity integration routine
void J2PlasticityThermal :: plastic_integrator( )
{
  const double tolerance = (1.0e-8)*sigma_y ;

  const double dt = ops_Dt ; //time step

  static Matrix dev_strain(3,3) ; //deviatoric strain

  static Matrix dev_stress(3,3) ; //deviatoric stress
 
  static Matrix normal(3,3) ;     //normal to yield surface

  double NbunN ; //normal bun normal   

  double norm_tau = 0.0 ;   //norm of deviatoric stress 
  double inv_norm_tau = 0.0 ;

  double phi = 0.0 ; //trial value of yield function

  double trace = 0.0 ; //trace of strain

  double gamma = 0.0 ; //consistency parameter

  double resid = 1.0 ; 
  double tang  = 0.0 ;
  
  double theta = 0.0 ; 
  double theta_inv = 0.0 ;

  double c1 = 0.0 ; 
  double c2 = 0.0 ;
  double c3 = 0.0 ;

  int i,j,k,l;
  int ii, jj ; 

  int iteration_counter ;
  const int max_iterations = 25 ;

  //compute the deviatoric strains

  trace = strain(0,0) + strain(1,1) + strain(2,2) ;
 
  dev_strain = strain ;
  for ( i = 0; i < 3; i++ )
    dev_strain(i,i) -= ( one3*trace ) ;
   
  //compute the trial deviatoric stresses

  //   dev_stress = (2.0*shear) * ( dev_strain - epsilon_p_n ) ;
  dev_stress = dev_strain;
  dev_stress -= epsilon_p_n;
  dev_stress *= 2.0 * shear;

  //compute norm of deviatoric stress

  norm_tau = 0.0 ;
  for ( i = 0; i < 3; i++ ){
    for ( j = 0; j < 3; j++ ) 
      norm_tau += dev_stress(i,j)*dev_stress(i,j) ;
  } //end for i 
   
  norm_tau = sqrt( norm_tau ) ;

  if ( norm_tau > tolerance ) {
    inv_norm_tau = 1.0 / norm_tau ;
    normal =  inv_norm_tau * dev_stress ;
  }
  else {
    normal.Zero( ) ;
    inv_norm_tau = 0.0 ;
  } //end if 

  //compute trial value of yield function

  phi = norm_tau -  root23 * q(xi_n) ;

  // check if phi > 0 
  
  if ( phi > 0.0 ) { //plastic

     //solve for gamma 
     gamma = 0.0 ;
     resid = 1.0 ;
     iteration_counter = 0 ;
     while ( fabs(resid) > tolerance ) {

        resid = norm_tau 
              - (2.0*shear) * gamma 
              - root23 * q( xi_n + root23*gamma ) 
              - (eta/dt) * gamma ;

        tang =  - (2.0*shear)  
                - two3 * qprime( xi_n + root23*gamma )
                - (eta/dt) ;

        gamma -= ( resid / tang ) ;

	iteration_counter++ ;

	if ( iteration_counter > max_iterations ) {
	    opserr << "More than " << max_iterations ;
 	    opserr << " iterations in constituive subroutine J2-plasticity \n" ;
	    break ;
	} //end if 
	
     } //end while resid

     gamma *= (1.0 - 1e-08) ;

     //update plastic internal variables

     epsilon_p_nplus1 = epsilon_p_n + gamma*normal ;

     xi_nplus1 = xi_n + root23*gamma ;

     //recompute deviatoric stresses 

     dev_stress = (2.0*shear) * ( dev_strain - epsilon_p_nplus1 ) ;

     //compute the terms for plastic part of tangent

     theta =  (2.0*shear)  
           +  two3 * qprime( xi_nplus1 )
           +  (eta/dt) ;

     theta_inv = 1.0/theta ;

  }
  else { //elastic 

    //update history variables -- they remain unchanged

    epsilon_p_nplus1 = epsilon_p_n ;

    xi_nplus1 = xi_n ;

    //no extra tangent terms to compute 
    
    gamma = 0.0 ; 
    theta = 0.0 ;
    theta_inv = 0.0 ;

  } //end if phi > 0


  //add on bulk part of stress

  stress = dev_stress ;
  for ( i = 0; i < 3; i++ )
     stress(i,i) += bulk*trace ;

  //compute the tangent

  c1 = -4.0 * shear * shear ;
  c2 = c1 * theta_inv ;
  c3 = c1 * gamma * inv_norm_tau ;

  for ( ii = 0; ii < 6; ii++ ) {
    for ( jj = 0; jj < 6; jj++ )  {

          index_map( ii, i, j ) ;
          index_map( jj, k, l ) ;

          NbunN  = normal(i,j)*normal(k,l) ; 

          //elastic terms
          tangent[i][j][k][l]  = bulk * IbunI[i][j][k][l] ;

          tangent[i][j][k][l] += (2.0*shear) * IIdev[i][j][k][l] ;

          //plastic terms 
          tangent[i][j][k][l] += c2 * NbunN ;

	  tangent[i][j][k][l] += c3 * (  IIdev[i][j][k][l] - NbunN ) ;

          //minor symmetries 
          tangent [j][i][k][l] = tangent[i][j][k][l] ;
          tangent [i][j][l][k] = tangent[i][j][k][l] ;
          tangent [j][i][l][k] = tangent[i][j][k][l] ;

    } // end for jj
  } // end for ii

  return ;
} 




double
J2PlasticityThermal::setThermalTangentAndElongation(double &tempT, double&ET, double&Elong)
{
	double TempT = tempT + 20;
	double E00; //Initial tangent 
	ET = 2E11;
	E00 = 2E11;

	// EN 1992 pt 1-2-1. Class N hot rolled  reinforcing steel at elevated temperatures
	
	if (TempT <= 100) {

	}
	else if (TempT <= 200) {

		bulk = bulk_0*(1 - (TempT - 100)*0.1 / 100);
		shear = shear_0*(1 - (TempT - 100)*0.1 / 100);
		sigma_y = sigma_0;
		ET = E00*(1 - (TempT - 100)*0.1 / 100);
		Hard = 0.01*ET / 2.8;

	}
	else if (TempT <= 300) {

		bulk = bulk_0*(0.9 - (TempT - 200)*0.1 / 100);
		shear = shear_0*(0.9 - (TempT - 200)*0.1 / 100);

		sigma_y = sigma_0;
		ET = E00*(0.9 - (TempT - 200)*0.1 / 100);
		Hard = 0.01*ET / 2.8;

	}
	else if (TempT <= 400) {
		bulk = bulk_0*(0.8 - (TempT - 300)*0.1 / 100);
		shear = shear_0*(0.8 - (TempT - 300)*0.1 / 100);

		sigma_y = sigma_0;
		ET = E00*(0.8 - (TempT - 300)*0.1 / 100);
		Hard = 0.01*ET / 2.8;
	}
	else if (TempT <= 500) {
		bulk = bulk_0*(0.7 - (TempT - 400)*0.1 / 100);
		shear = shear_0*(0.7 - (TempT - 400)*0.1 / 100);

		sigma_y = sigma_0*(1 - (TempT - 400)*0.22 / 100);
		ET = E00*(0.7 - (TempT - 400)*0.1 / 100);
		Hard = 0.01*ET / 2.8;
	}
	else if (TempT <= 600) {
		bulk = bulk_0*(0.6 - (TempT - 500)*0.29 / 100);
		shear = shear_0*(0.6 - (TempT - 500)*0.29 / 100);
		sigma_y = sigma_0*(0.78 - (TempT - 500)*0.31 / 100);
		ET = E00*(0.6 - (TempT - 500)*0.29 / 100);
		Hard = 0.01*ET / 2.8;
	}
	else if (TempT <= 700) {
		bulk = bulk_0*(0.31 - (TempT - 600)*0.18 / 100);
		shear = shear_0*(0.31 - (TempT - 600)*0.18 / 100);

		sigma_y = sigma_0*(0.47 - (TempT - 600)*0.24 / 100);
		ET = E00*(0.31 - (TempT - 600)*0.18 / 100);
		Hard = 0.01*ET / 2.8;
	}
	else if (TempT <= 800) {
		bulk = bulk_0*(0.13 - (TempT - 700)*0.04 / 100);
		shear = shear_0*(0.13 - (TempT - 700)*0.04 / 100);

		sigma_y = sigma_0*(0.23 - (TempT - 700)*0.12 / 100);
		ET = E00*(0.13 - (TempT - 700)*0.04 / 100);
		Hard = 0.01*ET / 2.8;
	}
	else if (TempT <= 900) {
		bulk = bulk_0*(0.09 - (TempT - 800)*0.02 / 100);
		shear = shear_0*(0.09 - (TempT - 800)*0.02 / 100);

		sigma_y = sigma_0*(0.11 - (TempT - 800)*0.05 / 100);
		ET = E00*(0.09 - (TempT - 800)*0.02 / 100);
		Hard = 0.01*ET / 2.8;
	}
	else if (TempT <= 1000) {
		bulk = bulk_0*(0.0675 - (TempT - 900)*(0.00675 - 0.0045) / 100);
		shear = shear_0*(0.0675 - (TempT - 900)*(0.00675 - 0.0045) / 100);

		sigma_y = sigma_0*(0.06 - (TempT - 900)*0.02 / 100);
		ET = E00*(0.0675 - (TempT - 900)*(0.00675 - 0.0045) / 100);
		Hard = 0.01*ET / 2.8;
	}

	else {
		opserr << "the temperature is invalid\n";
	}
	
	// Calculate thermal elongation 
	if (TempT <= 20) {
		ThermalElongation = 0.0;
	}
	else if (TempT <= 750) {
		ThermalElongation = -2.416e-4 + 1.2e-5 *TempT + 0.4e-8 *TempT*TempT;
		
	}
	else if (TempT <= 860) {
		ThermalElongation = 11e-3;
		
	}
	else if (TempT <= 1200) {
		ThermalElongation = -6.2e-3 + 2e-5*TempT;

	}
	else {
		opserr << "the temperature is invalid\n";
	}
	
	//ThermalElongation = 12e-6*(tempT);
	TempAndElong(0) = TempT - 20;
	TempAndElong(1) = ThermalElongation;
	//bulk = bulk_0;
	//shear = shear_0;
	//sigma_y = sigma_0;

	//ET = 2E11;
	//ET = E;  
	//ET = 3.84e10;
	Elong = ThermalElongation;
	this->plastic_integrator();
	return 0;
}




// set up for initial elastic
void J2PlasticityThermal :: doInitialTangent( )
{
  int ii,jj,i,j,k,l;

  //compute the deviatoric strains
  for ( ii = 0; ii < 6; ii++ ) {
    for ( jj = 0; jj < 6; jj++ )  {

          index_map( ii, i, j ) ;
          index_map( jj, k, l ) ;

          //elastic terms
          initialTangent[i][j][k][l]  = bulk * IbunI[i][j][k][l] ;
          initialTangent[i][j][k][l] += (2.0*shear) * IIdev[i][j][k][l] ;

          //minor symmetries 
          //minor symmetries 
          initialTangent [j][i][k][l] = initialTangent[i][j][k][l] ;
          initialTangent [i][j][l][k] = initialTangent[i][j][k][l] ;
          initialTangent [j][i][l][k] = initialTangent[i][j][k][l] ;

    } // end for jj
  } // end for ii

  return ;
} 



//hardening function
double J2PlasticityThermal :: q( double xi ) 
{
//  q(xi) = simga_infty + (sigma_y - sigma_infty)*exp(-delta*xi) + H*xi 

 return    sigma_infty
         + (sigma_y - sigma_infty)*exp(-delta*xi)
         + Hard*xi ;
}


//hardening function derivative
double J2PlasticityThermal :: qprime( double xi )
{
  return  (sigma_y - sigma_infty) * (-delta) * exp(-delta*xi)
         + Hard ;
}


//matrix_index ---> tensor indices i,j
void J2PlasticityThermal :: index_map( int matrix_index, int &i, int &j )
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
      i = 3 ;
      j = 3 ;
      break ;

    case 4 :
      i = 1 ;
      j = 2 ;
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


NDMaterial*
J2PlasticityThermal::getCopy (void)
{
  opserr << "J2PlasticityThermal::getCopy -- subclass responsibility\n"; 
  exit(-1);
  return 0;
}

const char*
J2PlasticityThermal::getType (void) const
{
    opserr << "J2PlasticityThermal::getType -- subclass responsibility\n";
    exit(-1);
    return 0;
}

int
J2PlasticityThermal::getOrder (void) const
{
    opserr << "J2PlasticityThermal::getOrder -- subclass responsibility\n";
    exit(-1);
    return 0;
}


int 
J2PlasticityThermal::commitState( ) 
{
  epsilon_p_n = epsilon_p_nplus1 ;
  xi_n        = xi_nplus1 ;

  return 0;
}

int 
J2PlasticityThermal::revertToLastCommit( ) 
{
  return 0;
}


int 
J2PlasticityThermal::revertToStart( ) {

  // added: C.McGann, U.Washington for InitialStateAnalysis
  if (ops_InitialStateAnalysis) {
	// do nothing, keep state variables from last step
  } else {
	// normal call for revertToStart (not initialStateAnalysis)
    this->zero( ) ;
  }

  return 0;
}

int
J2PlasticityThermal::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(10+9);
  int cnt = 0;
  data(cnt++) = this->getTag();
  data(cnt++) = bulk;
  data(cnt++) = shear;
  data(cnt++) = sigma_y;
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
    opserr << "J2PlasticityThermal::sendSelf - failed to send vector to channel\n";
    return -1;
  }

  return 0;
}

int
J2PlasticityThermal::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  // recv the vector object from the channel which defines material param and state
  static Vector data(10+9);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "J2PlasticityThermal::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag(data(cnt++));
  bulk = data(cnt++);
  shear = data(cnt++);
  sigma_y = data(cnt++);
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
