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

#ifndef ElasticPlaneStress_h
#define ElasticPlaneStress_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>

#include <NDMaterial.h>

class ElasticPlaneStress : public NDMaterial {

//-------------------Declarations-------------------------------

  public : 

  //null constructor
  ElasticPlaneStress( ) ;

  //full constructor
  ElasticPlaneStress(int tag, 
                   double E,
                   double nu,
                   double rho) ;


  //destructor
  ~ElasticPlaneStress( ) ;

  const char *getClassType(void) const {return "ElasticPlaneStress";};

  //make a clone of this material
  NDMaterial* getCopy( ) ;

  //send back type of material
  const char* getType( ) const ;

  //send back order of strain in vector form
  int getOrder( ) const ;

  //mass per unit volume
  double getRho();

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

  //swap history variables
  int commitState( ) ; 
  int revertToLastCommit( ) ;
  int revertToStart( ) ;

  //sending and receiving
  int sendSelf(int commitTag, Channel &theChannel) ;  
  int recvSelf(int commitTag, Channel &theChannel, 
               FEM_ObjectBroker &theBroker ) ;
  //print out material data
  void Print(OPS_Stream &s, int flag = 0) ;
  
  private : 
  
  //static vectors and matrices
  Vector strain_vec ;     //strain in vector notation
  static Vector stress_vec ;     //stress in vector notation
  static Matrix tangent_matrix ; //material tangent in matrix notation

  double E, nu, rho;

  //index mapping special for plane stress because of 
  // condensation on tangent
  void index_map( int matrix_index, int &i, int &j ) ;

} ; //end of ElasticPlaneStress declarations


#endif
