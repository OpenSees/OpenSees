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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/CycLiqCPSPPlaneStrain.cpp,v $

// Written: Rui Wang, Tsinghua University, August, 2013
//
// Cyclic constitutive model for post-liquefaction shear deformationof sand 
// 
//
//
//  Cutting Plane Integration Scheme 
//
//

#include <CycLiqCPSPPlaneStrain.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>


//static vectors and matrices
Vector CycLiqCPSPPlaneStrain :: strain_vec(3) ;
Vector CycLiqCPSPPlaneStrain :: stress_vec(3) ;
Matrix CycLiqCPSPPlaneStrain :: tangent_matrix(3,3) ;

//null constructor
CycLiqCPSPPlaneStrain :: CycLiqCPSPPlaneStrain() :
CycLiqCPSP()
{}

CycLiqCPSPPlaneStrain :: CycLiqCPSPPlaneStrain(    int    tag,
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
		       double rho1):
CycLiqCPSP(tag, ND_TAG_CycLiqCPSPPlaneStrain, G01, kappa1, h1, Mfc1, dre11,  dre21, rdr1, eta1, dir1, lamdac1, ksi1, e01, nb1, nd1, ein1, rho1 )
{

}

CycLiqCPSPPlaneStrain :: ~CycLiqCPSPPlaneStrain( ) 
{
}

//make a clone of this material
NDMaterial* CycLiqCPSPPlaneStrain :: getCopy( ) 
{ 
  CycLiqCPSPPlaneStrain  *clone;
  clone = new CycLiqCPSPPlaneStrain( ) ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}


//send back type of material
const char* CycLiqCPSPPlaneStrain :: getType( ) const 
{
  return "PlaneStrain" ;
}


//send back order of strain in vector form
int CycLiqCPSPPlaneStrain :: getOrder( ) const 
{ 
  return 3 ; 
} 


//get the strain and integrate plasticity equations
int CycLiqCPSPPlaneStrain :: setTrialStrain( const Vector &strain_from_element) 
{
  strain_nplus1.Zero( ) ;

  strain_nplus1(0,0) =        strain_from_element(0) ;
  strain_nplus1(1,1) =        strain_from_element(1) ;

  strain_nplus1(0,1) = 0.50 * strain_from_element(2) ;
  strain_nplus1(1,0) =        strain_nplus1(0,1) ;

  this->plastic_integrator( ) ;

  return 0 ;
}


//unused trial strain functions
int CycLiqCPSPPlaneStrain :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

int CycLiqCPSPPlaneStrain :: setTrialStrainIncr( const Vector &v ) 
{
  static Vector newStrain(3);
  newStrain(0) = strain_nplus1(0,0) + v(0);
  newStrain(1) = strain_nplus1(1,1) + v(1);
  newStrain(2) = 2.0*strain_nplus1(0,1) + v(2);
  return this->setTrialStrain(newStrain);
}

int CycLiqCPSPPlaneStrain :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
  return this->setTrialStrainIncr(v);
}



//send back the strain
const Vector& CycLiqCPSPPlaneStrain :: getStrain( ) 
{
  strain_vec(0) =       strain_nplus1(0,0) ;
  strain_vec(1) =       strain_nplus1(1,1) ;

  strain_vec(2) = 2.0 * strain_nplus1(0,1) ;


  return strain_vec ;
} 


//send back the stress 
const Vector& CycLiqCPSPPlaneStrain :: getStress( ) 
{
  stress_vec(0) = stress_nplus1(0,0) ;
  stress_vec(1) = stress_nplus1(1,1) ;

  stress_vec(2) = stress_nplus1(0,1) ;

  return stress_vec ;
}

//send back the tangent 
const Matrix& CycLiqCPSPPlaneStrain :: getTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 ) 
  // 
       
  tangent_matrix(0,0) = tangent [0][0] [0][0] ;
  tangent_matrix(1,1) = tangent [1][1] [1][1] ;
  tangent_matrix(2,2) = tangent [0][1] [0][1] ;

  tangent_matrix(0,1) = tangent [0][0] [1][1] ;
  tangent_matrix(1,0) = tangent [1][1] [0][0] ;

  tangent_matrix(0,2) = tangent [0][0] [0][1] ;
  tangent_matrix(2,0) = tangent [0][1] [0][0] ;

  tangent_matrix(1,2) = tangent [1][1] [0][1] ;
  tangent_matrix(2,1) = tangent [0][1] [1][1] ;

  return tangent_matrix ;
} 

//send back the tangent 
const Matrix& CycLiqCPSPPlaneStrain :: getInitialTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 ) 
  // 

  this->doInitialTangent();

  tangent_matrix(0,0) = initialTangent [0][0] [0][0] ;
  tangent_matrix(1,1) = initialTangent [1][1] [1][1] ;
  tangent_matrix(2,2) = initialTangent [0][1] [0][1] ;

  tangent_matrix(0,1) = initialTangent [0][0] [1][1] ;
  tangent_matrix(1,0) = initialTangent [1][1] [0][0] ;

  tangent_matrix(0,2) = initialTangent [0][0] [0][1] ;
  tangent_matrix(2,0) = initialTangent [0][1] [0][0] ;

  tangent_matrix(1,2) = initialTangent [1][1] [0][1] ;
  tangent_matrix(2,1) = initialTangent [0][1] [1][1] ;

  return tangent_matrix ;
} 

int
CycLiqCPSPPlaneStrain::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
	int res = 0;
  static Vector data(22+9*3);
  int cnt = 0;
  data(cnt++) = this->getTag();

  data(cnt++) =	   G0;
  data(cnt++) =	   kappa;
  data(cnt++) =	   h;
  data(cnt++) =	   Mfc;   
  data(cnt++) =	   dre1;
  data(cnt++) =	   Mdc;
  data(cnt++) =	   dre2;
  data(cnt++) =	   rdr;
  data(cnt++) =	   eta;
  data(cnt++) =	   dir;
  data(cnt++) =	   lamdac;
  data(cnt++) =	   ksi;
  data(cnt++) =	   e0;
  data(cnt++) =	   nb;
  data(cnt++) =	   nd;
  data(cnt++) =	   ein;   
  data(cnt++) =	   rho;
  data(cnt++) =	   epsvir_nplus1;
  data(cnt++) =	   epsvre_nplus1;
  data(cnt++) =	   gammamonos;   
  data(cnt++) =	   epsvc_nplus1;
  data(cnt++) =	   etam;


  for (int i=0; i<3; i++) 
  {
    for (int j=0; j<3; j++) 
	{
	  data(cnt+9)   = strain_nplus1(i,j);
	  data(cnt+9*2) = alpha_nplus1(i,j);
	  data(cnt+9*3) = stress_nplus1(i,j);
	}
  }

	res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  // send the vector object to the channel
  if (res < 0) {
    opserr << "CycLiqCPSP::sendSelf - failed to send vector to channel\n";
    return res;
  }

  return res;
}

int
CycLiqCPSPPlaneStrain::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
	int res = 0;
  // recv the vector object from the channel which defines material param and state
  static Vector data(22+9*3);
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "CycLiqCPSP::recvSelf - failed to recv vector from channel\n";
    return res;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag(data(cnt++));
  G0 =	data(cnt++);    
  kappa =	data(cnt++);    
  h =	data(cnt++);    
  Mfc    =	data(cnt++);    
  dre1 =	data(cnt++);    
  Mdc =	data(cnt++);    
  dre2 =	data(cnt++);    
  rdr =	data(cnt++);    
  eta =	data(cnt++);    
  dir =	data(cnt++); 
  lamdac =	data(cnt++);    
  ksi =	data(cnt++);    
  e0 =	data(cnt++);    
  nb =	data(cnt++);    
  nd =	data(cnt++);
  ein    =	data(cnt++);    
  rho =	data(cnt++);    
  epsvir_n =	data(cnt++);    
  epsvre_n =	data(cnt++);    
  gammamono  =	data(cnt++);      
  epsvc_n =	data(cnt++);    
  etam =	data(cnt++);    


  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++) 
	{
      strain_n(i,j) = data(cnt+9);
	  alpha_n(i,j) = data(cnt+9*2);
	  stress_n(i,j) = data(cnt+9*3);
	}
  }

  return res;
}
