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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2Plasticity.cpp,v $

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
//  set eta := 0 for rate independent case
//

#include <J2Plasticity.h>
#include <J2PlaneStress.h>
#include <J2PlaneStrain.h>
#include <J2AxiSymm.h>
#include <J2ThreeDimensional.h> 

//this is mike's problem
Tensor J2Plasticity :: rank2(2, def_dim_2, 0.0 ) ;
Tensor J2Plasticity :: rank4(2, def_dim_2, 0.0 ) ;

//parameters
const double J2Plasticity :: one3   = 1.0 / 3.0 ;
const double J2Plasticity :: two3   = 2.0 / 3.0 ;
const double J2Plasticity :: four3  = 4.0 / 3.0 ;
const double J2Plasticity :: root23 = sqrt( 2.0 / 3.0 ) ;


//zero internal variables
void J2Plasticity :: zero ( ) 
{
  xi_n = 0.0 ;
  xi_nplus1 = 0.0 ;
  
  epsilon_p_n.Zero( ) ;
  epsilon_p_nplus1.Zero( ) ;
}


//null constructor
J2Plasticity ::  J2Plasticity( ) : 
NDMaterial( ),
epsilon_p_n(3,3),
epsilon_p_nplus1(3,3),
stress(3,3),
strain(3,3)
{ 
  bulk        = 0.0 ;
  shear       = 0.0 ;
  sigma_0     = 0.0 ;
  sigma_infty = 0.0 ;
  delta       = 0.0 ;
  Hard        = 0.0 ;
  eta         = 0.0 ;
  this->zero( ) ;     // or (*this).zero( ) 
}


//full constructor
J2Plasticity :: J2Plasticity(int    tag,
			     int classTag,
			     double K,
			     double G,
			     double yield0,
			     double yield_infty,
			     double d,
			     double H,
			     double viscosity) 
: 
  NDMaterial(tag, classTag),
  epsilon_p_n(3,3),
  epsilon_p_nplus1(3,3),
  stress(3,3),
  strain(3,3)
{
  bulk        = K ;
  shear       = G ;
  sigma_0     = yield0 ;
  sigma_infty = yield_infty ;
  delta       = d ;
  Hard        = H ;
  eta         = viscosity ;
  this->zero( ) ;
}


//elastic constructor
J2Plasticity :: 
J2Plasticity(   int    tag, 
                int  classTag,
                double K, 
                double G ) :
NDMaterial(tag, classTag),
epsilon_p_n(3,3),
epsilon_p_nplus1(3,3),
stress(3,3),
strain(3,3)
{
  bulk        = K ;
  shear       = G ; 
  sigma_0     = 1.0e16*shear ;
  sigma_infty = sigma_0 ;
  delta       = 0.0 ;
  Hard        = 0.0 ;
  eta         = 0.0 ;
  this->zero( ) ;
}


//destructor
J2Plasticity :: ~J2Plasticity( ) 
{  } 



NDMaterial*
J2Plasticity :: getCopy (const char *type)
{
    if (strcmp(type,"PlaneStress2D") == 0)
    {
	J2PlaneStress  *clone ;
	clone = new J2PlaneStress(this->getTag(), bulk, shear, sigma_0,
				  sigma_infty, delta, Hard, eta) ;
	return clone ;
    }
    else if (strcmp(type,"PlaneStrain2D") == 0)
    {
	J2PlaneStrain  *clone ;
	clone = new J2PlaneStrain(this->getTag(), bulk, shear, sigma_0,
				  sigma_infty, delta, Hard, eta) ;
	return clone ;
    }
    else if (strcmp(type,"AxiSymmetric2D") == 0)
    {
	J2AxiSymm  *clone ;
	clone = new J2AxiSymm(this->getTag(), bulk, shear, sigma_0,
			      sigma_infty, delta, Hard, eta) ;
	return clone ;	
    }
    else if ((strcmp(type,"ThreeDimensional") == 0) ||
	     (strcmp(type,"3D") == 0))
    {
	J2ThreeDimensional  *clone ;
	clone = new J2ThreeDimensional(this->getTag(), bulk, shear, sigma_0,
			      sigma_infty, delta, Hard, eta) ;
	return clone ;	
    }
 
    // Handle other cases
    else
    {
	g3ErrorHandler->fatal("J2Plasticity::getModel failed to get model %s",
			      type);

	return 0;
    }
}

//swap history variables
int J2Plasticity :: commitState( )  
{
  epsilon_p_n = epsilon_p_nplus1 ;
  xi_n        = xi_nplus1 ;

  return 0 ;
}


//revert to last saved state
int J2Plasticity :: revertToLastCommit( )
{ 
  return 0 ;
} 

//revert to start
int J2Plasticity :: revertToStart( ) 
{
  this->zero( ) ;
  return 0 ;
}


//print out material data
void J2Plasticity :: Print( ostream &s, int flag )
{
  s << '\n' ;
  s << "J2-Plasticity : " ; 
  s << this->getType( ) << '\n' ;
  s << "Bulk Modulus =   " << bulk        << '\n' ;
  s << "Shear Modulus =  " << shear       << '\n' ;
  s << "Sigma_0 =        " << sigma_0     << '\n' ;
  s << "Sigma_infty =    " << sigma_infty << '\n' ;
  s << "Delta =          " << delta       << '\n' ;
  s << "H =              " << Hard        << '\n' ;
  s << "Eta =            " << eta         << '\n' ;
  s << endl ;
}


//--------------------Plasticity-------------------------------------

//plasticity integration routine
void J2Plasticity :: plastic_integrator( )
{
  const double tolerance = 1.0e-12 ;

  const double dt = 1.0 ; //time step = 1 for now

  static Matrix dev_strain(3,3) ; //deviatoric strain

  static Matrix dev_stress(3,3) ; //deviatoric stress
 
  static Matrix normal(3,3) ;     //normal to yield surface

  double IIdev[3][3][3][3] ; //rank 4 deviatoric 

  double IbunI[3][3][3][3] ; //rank 4 I bun I 

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

  int ii, jj ; 
  int i, j, k, l ;

  int iteration_counter ;
  const int max_iterations = 25 ;

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



  //compute the deviatoric strains

  trace = strain(0,0) + strain(1,1) + strain(2,2) ;
 
  dev_strain = strain ;
  for ( i = 0; i < 3; i++ )
    dev_strain(i,i) -= ( one3*trace ) ;
   
  //compute the trial deviatoric stresses

  dev_stress = (2.0*shear) * ( dev_strain - epsilon_p_n ) ;

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
	    cerr << "More than " << max_iterations ;
	    cerr << " iterations in constituive subroutine J2-plasticity \n" ;
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



//hardening function
double J2Plasticity :: q( double xi ) 
{
//  q(xi) = simga_infty + (sigma_0 - sigma_infty)*exp(-delta*xi) + H*xi 

 return    sigma_infty
         + (sigma_0 - sigma_infty)*exp(-delta*xi)
         + Hard*xi ;
}


//hardening function derivative
double J2Plasticity :: qprime( double xi )
{
  return  (sigma_0 - sigma_infty) * (-delta) * exp(-delta*xi)
         + Hard ;
}


//matrix_index ---> tensor indices i,j
void J2Plasticity :: index_map( int matrix_index, int &i, int &j )
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
J2Plasticity::getCopy (void)
{
    g3ErrorHandler->fatal("J2Plasticity::getCopy -- subclass responsibility");
    return 0;
}

const char*
J2Plasticity::getType (void) const
{
    g3ErrorHandler->fatal("J2Plasticity::getType -- subclass responsibility");
    return 0;
}

int
J2Plasticity::getOrder (void) const
{
    g3ErrorHandler->fatal("J2Plasticity::getOrder -- subclass responsibility");
    return 0;
}

int
J2Plasticity::sendSelf (int commitTag, Channel &theChannel)
{
    g3ErrorHandler->fatal("J2Plasticity::sendSelf -- subclass responsibility");
    return 0;
}

int
J2Plasticity::recvSelf (int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
    g3ErrorHandler->fatal("J2Plasticity::recvSelf -- subclass responsibility");
    return 0;
}
