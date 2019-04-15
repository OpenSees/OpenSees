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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/CycLiqCPSP3D.cpp,v $

// Written: Rui Wang, Tsinghua University, August, 2013
//
// Cyclic constitutive model for post-liquefaction shear deformationof sand 
// 
//
//
//  Cutting Plane Integration Scheme 
//
//

#include <CycLiqCPSP3D.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>


//static vectors and matrices
Vector CycLiqCPSP3D :: strain_vec(6) ;
Vector CycLiqCPSP3D :: stress_vec(6) ;
Matrix CycLiqCPSP3D :: tangent_matrix(6,6) ;

//null constructor
CycLiqCPSP3D :: CycLiqCPSP3D() :
CycLiqCPSP()
{}

CycLiqCPSP3D :: CycLiqCPSP3D(    int    tag,
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
CycLiqCPSP(tag, ND_TAG_CycLiqCPSP3D, G01, kappa1, h1, Mfc1, dre11, dre21, rdr1, eta1, dir1, lamdac1, ksi1, e01, nb1, nd1, ein1, rho1 )
{

}

CycLiqCPSP3D :: ~CycLiqCPSP3D( ) 
{
}

//make a clone of this material
NDMaterial* CycLiqCPSP3D :: getCopy( ) 
{ 
  CycLiqCPSP3D  *clone;
  clone = new CycLiqCPSP3D( ) ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}


//send back type of material
const char* CycLiqCPSP3D :: getType( ) const 
{
  return "ThreeDimensional" ;
}


//send back order of strain in vector form
int CycLiqCPSP3D :: getOrder( ) const 
{ 
  return 6 ; 
} 


//get the strain and integrate plasticity equations
int CycLiqCPSP3D :: setTrialStrain( const Vector &strain_from_element) 
{
  strain_nplus1.Zero( ) ;

  strain_nplus1(0,0) =        strain_from_element(0) ;
  strain_nplus1(1,1) =        strain_from_element(1) ;
  strain_nplus1(2,2) =        strain_from_element(2) ;

  strain_nplus1(0,1) = 0.50 * strain_from_element(3) ;
  strain_nplus1(1,0) =        strain_nplus1(0,1) ;

  strain_nplus1(1,2) = 0.50 * strain_from_element(4) ;
  strain_nplus1(2,1) =        strain_nplus1(1,2) ;
  
  strain_nplus1(2,0) = 0.50 * strain_from_element(5) ;
  strain_nplus1(0,2) =        strain_nplus1(2,0) ;

  this->plastic_integrator( ) ;

  return 0 ;
}


//unused trial strain functions
int CycLiqCPSP3D :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

int CycLiqCPSP3D :: setTrialStrainIncr( const Vector &v ) 
{
  static Vector newStrain(6);
  newStrain(0) = strain_nplus1(0,0) + v(0);
  newStrain(1) = strain_nplus1(1,1) + v(1);
  newStrain(2) = strain_nplus1(2,2) + v(2);
  newStrain(3) = 2.0*strain_nplus1(0,1) + v(3);
  newStrain(4) = 2.0*strain_nplus1(1,2) + v(4);
  newStrain(5) = 2.0*strain_nplus1(2,0) + v(5);
  return this->setTrialStrain(newStrain);
}

int CycLiqCPSP3D :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
  return this->setTrialStrainIncr(v);
}



//send back the strain
const Vector& CycLiqCPSP3D :: getStrain( ) 
{
  strain_vec(0) =       strain_nplus1(0,0) ;
  strain_vec(1) =       strain_nplus1(1,1) ;
  strain_vec(2) =       strain_nplus1(2,2) ;

  strain_vec(3) = 2.0 * strain_nplus1(0,1) ;

  strain_vec(4) = 2.0 * strain_nplus1(1,2) ;

  strain_vec(5) = 2.0 * strain_nplus1(2,0) ;

  return strain_vec ;
} 


//send back the stress 
const Vector& CycLiqCPSP3D :: getStress( ) 
{
  stress_vec(0) = stress_nplus1(0,0) ;
  stress_vec(1) = stress_nplus1(1,1) ;
  stress_vec(2) = stress_nplus1(2,2) ;

  stress_vec(3) = stress_nplus1(0,1) ;

  stress_vec(4) = stress_nplus1(1,2) ;
  
  stress_vec(5) = stress_nplus1(2,0) ;

  return stress_vec ;
}

//send back the tangent 
const Matrix& CycLiqCPSP3D :: getTangent( ) 
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
	  if (jj>2)
	  {
		  tangent_matrix(ii,jj) = tangent_matrix(ii,jj);
	  }

    } //end for j
  } //end for i


  return tangent_matrix ;
} 

//send back the tangent 
const Matrix& CycLiqCPSP3D :: getInitialTangent( ) 
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

	  if (jj>2)
	  {
		  tangent_matrix(ii,jj) = tangent_matrix(ii,jj);
	  }

    } //end for j
  } //end for i

  return tangent_matrix ;
} 

int
CycLiqCPSP3D::sendSelf(int commitTag, Channel &theChannel)
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
CycLiqCPSP3D::recvSelf (int commitTag, Channel &theChannel, 
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
