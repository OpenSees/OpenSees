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


/*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                          |
 |              MultiaxialCyclicPlasticity  NDMaterial                      |
 +                                                                          +
 |--------------------------------------------------------------------------|
 |                                                                          |
 +             Authors: Gang Wang  AND  Professor Nicholas Sitar            +
 |                                                                          |
 |             Department of Civil and Environmental Engineering            |
 +             University of California, Berkeley, CA 94720, USA            +
 |                                                                          |
 |             Email: wang@ce.berkeley.edu (G.W.)                           |
 |                                                                          |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/


#ifndef MultiaxialCyclicPlasticity3D_h
#define MultiaxialCyclicPlasticity3D_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>

#include <MultiaxialCyclicPlasticity.h>

class MultiaxialCyclicPlasticity3D : public MultiaxialCyclicPlasticity {

//-------------------Declarations-------------------------------

  public : 

  //null constructor
  MultiaxialCyclicPlasticity3D( ) ;

  //full constructor
  MultiaxialCyclicPlasticity3D(   int    tag,
				  double rho,
				  double K,
				  double G,
				  double Su,
				  double Ho_kin,
				  double Parameter_h,
				  double Parameter_m,
				  double Parameter_beta,
				  double Kcoeff,
				  double viscosity = 0
				  ) ;
  

  //elastic constructor
  MultiaxialCyclicPlasticity3D( int tag, double rho, double K, double G ) ;

  //destructor
  ~MultiaxialCyclicPlasticity3D( ) ;

  const char *getClassType(void) const {return "MultiaxialCyclicPlasticity3D";};

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

  private :

  //static vectors and matrices
  static Vector strain_vec ;     //strain in vector notation
  static Vector stress_vec ;     //stress in vector notation
  static Matrix tangent_matrix ; //material tangent in matrix notation

} ; //end of MultiaxialCyclicPlasticity3D declarations

#endif
