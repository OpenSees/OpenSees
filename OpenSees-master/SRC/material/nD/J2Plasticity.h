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
// $Date: 2008-10-21 18:58:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2Plasticity.h,v $

#ifndef J2Plasticity_h
#define J2Plasticity_h

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
//  set eta := 0 for rate independent case
//


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>


class J2Plasticity : public NDMaterial {

//-------------------Declarations-------------------------------

  public : 

  //null constructor
  J2Plasticity() ;

  //full constructor
  J2Plasticity(    int    tag,
		   int    classTag,
                   double K,
                   double G,
                   double yield0,
                   double yield_infty,
                   double d,
                   double H,
                   double viscosity = 0,
		   double rho = 0.0) ;

  //elastic constructor
  J2Plasticity( int tag, int classTag, double K, double G ) ;

  //destructor
  virtual ~J2Plasticity( ) ;

  virtual const char *getClassType(void) const {return "J2Plasticity";};

  virtual NDMaterial* getCopy (const char *type);

  //swap history variables
  virtual int commitState( ) ; 

  //revert to last saved state
  virtual int revertToLastCommit( ) ;

  //revert to start
  virtual int revertToStart( ) ;

  //sending and receiving
  virtual int sendSelf(int commitTag, Channel &theChannel) ;  
  virtual int recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker ) ;

  //print out material data
  void Print(OPS_Stream &s, int flag = 0) ;

  virtual NDMaterial *getCopy (void) ;
  virtual const char *getType (void) const ;
  virtual int getOrder (void) const ;

  double getRho(void) {return rho;}

  virtual int setParameter(const char **argv, int argc, Parameter &param);
  virtual int updateParameter(int parameterID, Information &info);
  virtual int activateParameter(int paramID);

  protected :

  //material parameters
  double bulk ;        //bulk modulus
  double shear ;       //shear modulus
  double sigma_0 ;     //initial yield stress
  double sigma_infty ; //final saturation yield stress
  double delta ;       //exponential hardening parameter
  double Hard ;        //linear hardening parameter
  double eta ;         //viscosity

  //internal variables
  Matrix epsilon_p_n ;       // plastic strain time n
  Matrix epsilon_p_nplus1 ;  // plastic strain time n+1
  double xi_n ;              // xi time n
  double xi_nplus1 ;         // xi time n+1

  //material response 
  Matrix stress ;                //stress tensor
  double tangent[3][3][3][3] ;   //material tangent
  static double initialTangent[3][3][3][3] ;   //material tangent
  static double IIdev[3][3][3][3] ; //rank 4 deviatoric 
  static double IbunI[3][3][3][3] ; //rank 4 I bun I 

  //material input
  Matrix strain ;               //strain tensor

  //parameters
  static const double one3 ;
  static const double two3 ;
  static const double four3 ;
  static const double root23 ;

  //zero internal variables
  void zero( ) ;

  //plasticity integration routine
  void plastic_integrator( ) ;

  void doInitialTangent( ) ;

  //hardening function
  double q( double xi ) ;

  //hardening function derivative
  double qprime( double xi ) ;
  
  //matrix index to tensor index mapping
  virtual void index_map( int matrix_index, int &i, int &j ) ;

  double rho;
  
  int parameterID;

} ; //end of J2Plasticity declarations

#endif
