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
                                                                        
// $Revision: 1 $
// $Date: 2012-09-01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/CycLiqCPSP3D.h,v $

// Written: Rui Wang, Tsinghua University, August, 2013
//
// Cyclic constitutive model for post-liquefaction shear deformationof sand 
// 
//
//  Cutting Plane Integration Scheme 
//
//

#ifndef CycLiqCPSP3D_h
#define CycLiqCPSP3D_h



#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
//#include <T2Vector.h>
#include <NDMaterial.h>

#include <CycLiqCPSP.h>

class CycLiqCPSP3D : public CycLiqCPSP {

//-------------------Declarations-------------------------------

  public : 

  //null constructor
  CycLiqCPSP3D() ;

  //full constructor
  CycLiqCPSP3D(    int    tag,
	           double G01,
	           double kappa1,
	           double h1,
	           double Mfc1,       //critical state
	           double dre11,
	           double Mdc1,
	           double dre21,
	           double rdr1,
	           double eta1,
	           double dir1,
			   double lamdac1,
			   double ksi1,
			   double e01,
			   double nb1,
			   double nd1,
	           double ein1,      //initial void ratio
		       double rho1=0.0) ;

  //destructor
  ~CycLiqCPSP3D( ) ;

  const char *getClassType(void) const {return "CycLiqCPSP3D";};
  
  //make a clone of this material
  NDMaterial* getCopy ();

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
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);  




  private:

  //static vectors and matrices sent back in get functions
  static Vector strain_vec ;     //strain in vector notation
  static Vector stress_vec ;     //stress in vector notation
  static Matrix tangent_matrix ; //material tangent in matrix notation

} ; //end of CycLiqCPSP declarations

#endif
