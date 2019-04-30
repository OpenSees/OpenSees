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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/CycLiqCP3D.cpp,v $

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

#include <CycLiqCP3D.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>


//static vectors and matrices
Vector CycLiqCP3D :: strain_vec(6) ;
Vector CycLiqCP3D :: stress_vec(6) ;
Matrix CycLiqCP3D :: tangent_matrix(6,6) ;

//null constructor
CycLiqCP3D :: CycLiqCP3D() :
CycLiqCP()
{}

CycLiqCP3D :: CycLiqCP3D(    int    tag,
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
		       double rho1):
CycLiqCP(tag, ND_TAG_CycLiqCP3D, G01, kappa1, h1, Mfc1, dre11, Mdc1, dre21, rdr1, eta1, dir1, ein1, rho1 )
{

}

CycLiqCP3D :: ~CycLiqCP3D( ) 
{
}

//make a clone of this material
NDMaterial* CycLiqCP3D :: getCopy( ) 
{ 
  CycLiqCP3D  *clone;
  clone = new CycLiqCP3D( ) ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}


//send back type of material
const char* CycLiqCP3D :: getType( ) const 
{
  return "ThreeDimensional" ;
}


//send back order of strain in vector form
int CycLiqCP3D :: getOrder( ) const 
{ 
  return 6 ; 
} 


//get the strain and integrate plasticity equations
int CycLiqCP3D :: setTrialStrain( const Vector &strain_from_element) 
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
int CycLiqCP3D :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   return this->setTrialStrain( v ) ;
} 

int CycLiqCP3D :: setTrialStrainIncr( const Vector &v ) 
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

int CycLiqCP3D :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
  return this->setTrialStrainIncr(v);
}



//send back the strain
const Vector& CycLiqCP3D :: getStrain( ) 
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
const Vector& CycLiqCP3D :: getStress( ) 
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
const Matrix& CycLiqCP3D :: getTangent( ) 
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

    } //end for j
  } //end for i


  return tangent_matrix ;
} 

//send back the tangent 
const Matrix& CycLiqCP3D :: getInitialTangent( ) 
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

    } //end for j
  } //end for i

  return tangent_matrix ;
} 

int
CycLiqCP3D::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(20+9*3);
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
  data(cnt++) =	   ein;   
  data(cnt++) =	   rho;
  data(cnt++) =	   epsvir_nplus1;
  data(cnt++) =	   epsvre_nplus1;
  data(cnt++) =	   gammamonos;   
  data(cnt++) =	   epsvc_nplus1;
  data(cnt++) =	   etam;
  data(cnt++) =	   epsvc0;
  data(cnt++) =	   p0;


  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++) 
	{
	  data(cnt+9)   = strain_nplus1(i,j);
	  data(cnt+9*2) = alpha_nplus1(i,j);
	  data(cnt+9*3) = stress_nplus1(i,j);
	}
  }

  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "CycLiqCP::sendSelf - failed to send vector to channel\n";
    return -1;
  }

  return 0;
}

int
CycLiqCP3D::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{

  // recv the vector object from the channel which defines material param and state
  static Vector data(20+9*3);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "CycLiqCP::recvSelf - failed to recv vector from channel\n";
    return -1;
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
  ein    =	data(cnt++);    
  rho =	data(cnt++);    
  epsvir_n =	data(cnt++);    
  epsvre_n =	data(cnt++);    
  gammamono  =	data(cnt++);      
  epsvc_n =	data(cnt++);    
  etam =	data(cnt++);    
  epsvc0 =	data(cnt++);    
  p0 =	data(cnt++);  

  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++) 
	{
      strain_n(i,j) = data(cnt+9);
	  alpha_n(i,j) = data(cnt+9*2);
	  stress_n(i,j) = data(cnt+9*3);
	}
  }


  return 0;
}