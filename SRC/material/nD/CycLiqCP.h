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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/CycLiqCP.h,v $

// Written: Rui Wang, Tsinghua University, September, 2012
//
// Cyclic constitutive model for liquefaction of sand 
// 
//
// Please refer to Zhang and Wang, 2012 "Large post-liquefaction deformation of sand, part I: physical mechanism, constitutive description and numerical algorithm"
// Acta Geotechnica
//
//  Cutting Plane Integration Scheme 
//
//

#ifndef CycLiqCP_h
#define CycLiqCP_h



#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
//#include <T2Vector.h>
#include <NDMaterial.h>


class CycLiqCP : public NDMaterial {

//-------------------Declarations-------------------------------

  public : 

  //null constructor
  CycLiqCP() ;

  //full constructor
  CycLiqCP(int    tag,
	   int classTag,
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
	   double ein1,      //initial void ratio
	   double rho1=0.0) ;
  
  //destructor
  ~CycLiqCP( ) ;
  
  virtual const char *getClassType(void) const {return "CycLiqCP";};
  
  //make a clone of this material
  virtual NDMaterial *getCopy(const char *type);
  virtual NDMaterial* getCopy (void);
  
  //send back type of material
  virtual const char* getType( void ) const ;
  
  //send back order of strain in vector form
  virtual int getOrder( void ) const ;
  
  //swap history variables
  int commitState( ) ; 
  
  //revert to last saved state
  int revertToLastCommit( ) ;
  
  //revert to start
  int revertToStart( ) ;

  //sending and receiving
  int sendSelf(int commitTag, Channel &theChannel) ;  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker ) ;
  
  //print out material data
  void Print(OPS_Stream &s, int flag = 0) ;
  
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int responseID, Information &eleInformation);
  
  double getRho(void) {return rho;}
  
  protected :

  //cutoff pressure
  //double p_cut=
  
  //material parameters
  double G0;
  double kappa;
  double h;
  double Mfc;       //critical state
  double dre1;
  double Mdc;
  double dre2;
  double rdr;
  double eta;
  double dir;
  double ein;      //initial void ratio
  double rho;      //saturated mass density
  
  
  
  static const double pat;
  
  double G;
  double K;
  double en; //void ratio at step n
  double enplus1;
  double D;
  double lambda;
  double dlambda;
  double N;
  
  //commit variables
  
  double epsvir_n;
  double epsvir_nplus1;
  double epsvre_n;
  double epsvre_nplus1;
  double epsvc_n;
  double epsvc_nplus1;
  double gammarem;
  double etam;
  double etamplus1;
  double epsvc0;
  Matrix alpha_n;            // back stress tensor time n
  Matrix alpha_nplus1;            // back stress tensor time n+1
  Matrix strain_n;           // strain time n (previous step)
  Matrix strain_nplus1;           // strain time n (previous step)
  Matrix stress_n ;                //stress tensor at n
  Matrix stress_nplus1;            //stress tensor at n+1
  double epsvc_ns0,epsvc_ns01;
  
  //internal variables
  Matrix R ;       // n+1/3*D*I
  Matrix L ;  // n-1/3*n:r*I
  Matrix r;                 // dev_stress/p_n
  Matrix rbar;
  Matrix r_nplus1;
  
  double p_n ;              // confining stress time n
  double p_nplus1 ;         // confining stress time n+1
  double p0;
  double gammamono;
  double gammamonos;
  double loadindex;
  double roubar;
  double rou;
  double H;
  double eta_n;
  double eta_nplus1;
  double Dre_n;
  double Dir_n;
  double prem;

  static const double pcut;    // cut off confining stress
  static const double pmin;      // cut off confining stress
  static double mElastFlag;    // elastic/plastic stage
  bool initializeState;
                // (dev_stress-alpha_n)/p_n/sqrt(2/3)/m


  //material response 
  double tangent[3][3][3][3] ;   //material tangent
  static double initialTangent[3][3][3][3] ;   //material tangent
  static double IIdev[3][3][3][3] ; //rank 4 deviatoric 
  static double IbunI[3][3][3][3] ; //rank 4 I bun I 
  static Matrix I; //rank 2 I


  //parameters
  static const double one3 ;
  static const double two3 ;
  static const double four3 ;
  static const double root23 ;
  static const double root32 ;

  //zero internal variables
  void zero( ) ;

  //plasticity integration routine
  void plastic_integrator( ) ;

  void doInitialTangent( ) ;

  //double contraction function for two second tensors
  double doublecontraction(Matrix a, Matrix b);
  Matrix doublecontraction42(double tangent[][3][3][3], Matrix a);
  Matrix doublecontraction24(Matrix a, double tangent[][3][3][3]);
  
  //matrix index to tensor index mapping
  void index_map( int matrix_index, int &i, int &j ) ;


  private:

  //static vectors and matrices sent back in get functions
  static Vector strain_vec ;     //strain in vector notation
  static Vector stress_vec ;     //stress in vector notation
  static Matrix tangent_matrix ; //material tangent in matrix notation

} ; //end of CycLiqCP declarations

#endif
