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

#include <MultiaxialCyclicPlasticityAxiSymm.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static vectors and matrices
Vector MultiaxialCyclicPlasticityAxiSymm :: strain_vec(4) ;
Vector MultiaxialCyclicPlasticityAxiSymm :: stress_vec(4) ;
Matrix MultiaxialCyclicPlasticityAxiSymm :: tangent_matrix(4,4) ;

//null constructor
MultiaxialCyclicPlasticityAxiSymm ::  MultiaxialCyclicPlasticityAxiSymm( ) : 
MultiaxialCyclicPlasticity( ) 
{  }


//full constructor
MultiaxialCyclicPlasticityAxiSymm :: 
MultiaxialCyclicPlasticityAxiSymm(       int    tag, 
				 double rho,
                 double K,
                 double G,
				 double Su,
 			     double Ho_kin,
                 double Parameter_h,
                 double Parameter_m,
                 double Parameter_beta,
                 double Kcoeff,
                 double viscosity ) : 
MultiaxialCyclicPlasticity( tag, ND_TAG_MultiaxialCyclicPlasticityAxiSymm, rho, K, G, 
      Su, Ho_kin, Parameter_h, Parameter_m, Parameter_beta, Kcoeff, viscosity)
{ 

}



//elastic constructor
MultiaxialCyclicPlasticityAxiSymm :: 
MultiaxialCyclicPlasticityAxiSymm(   int    tag, double rho,    // add density by Gang Wang
                 double K, 
                 double G ) :
MultiaxialCyclicPlasticity( tag, ND_TAG_MultiaxialCyclicPlasticityAxiSymm, rho, K, G )
{ 

}



//destructor
MultiaxialCyclicPlasticityAxiSymm :: ~MultiaxialCyclicPlasticityAxiSymm( ) 
{ } 


//make a clone of this material
NDMaterial* MultiaxialCyclicPlasticityAxiSymm :: getCopy( ) 
{ 
  MultiaxialCyclicPlasticityAxiSymm  *clone;
  clone = new MultiaxialCyclicPlasticityAxiSymm() ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}


//send back type of material
const char* MultiaxialCyclicPlasticityAxiSymm :: getType( ) const 
{
  return "AxiSymmetric" ;
}


//send back order of strain in vector form
int MultiaxialCyclicPlasticityAxiSymm :: getOrder( ) const 
{ 
  return 4 ; 
} 


//get the strain and integrate plasticity equations
int MultiaxialCyclicPlasticityAxiSymm :: setTrialStrain( const Vector &strain_from_element) 
{
  strain.Zero( ) ;

  strain(0,0) =        strain_from_element(0) ;
  strain(1,1) =        strain_from_element(1) ;
  strain(2,2) =        strain_from_element(2) ;

  strain(0,1) = 0.50 * strain_from_element(3) ;
  strain(1,0) =        strain(0,1) ;

//  this->plastic_integrator( ) ;

   //opserr<<"MCP:integrator"<<endln;

   if (this->MaterialStageID ==1) { 
 	  this->elastic_integrator( ) ;
  } else if (this->MaterialStageID ==2) {	  
 	  this->plastic_integrator( ) ;
  }

  return 0 ;
}


//unused trial strain functions
int MultiaxialCyclicPlasticityAxiSymm :: setTrialStrain( const Vector &v, const Vector &r )
{ 
  //opserr<<"MCP::setTrialStrain v"<<v<<endln;
   return this->setTrialStrain( v ) ;
} 

int MultiaxialCyclicPlasticityAxiSymm :: setTrialStrainIncr( const Vector &v ) 
{
	// implemented Dec. 5, 2003, Gang Wang
	static Vector newStrain(4);
	newStrain(0) = strain(0,0) + v(0);
	newStrain(1) = strain(1,1) + v(1);
	newStrain(2) = strain(2,2) + v(2);
	newStrain(3) = 2.0*strain(0,1) + v(3);
	
	//opserr<<"MCP::setTrialStrainIncr "<<newStrain<<endln;
	return this->setTrialStrain(newStrain);
	// return -1 ;  
}

int MultiaxialCyclicPlasticityAxiSymm :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
    // implemented Dec. 5, 2003, Gang Wang
	//return this->setTrialStrainIncr(v);

   // implemented April. 4, 2004, Gang Wang
	static Vector newStrain(4);
	newStrain(0) = strain(0,0) + v(0);
	newStrain(1) = strain(1,1) + v(1);
	newStrain(2) = strain(2,2) + v(2);
	newStrain(3) = 2.0*strain(0,1) + v(3);
    	
	opserr<<"MCP::setTrialStrainIncr"<<strain; 

    //this->ResidStress(0,0)= r(0);
    //this->ResidStress(1,1)= r(1);
    //this->ResidStress(2,2)= r(2);
    //this->ResidStress(0,1)= r(3);
    //this->ResidStress(1,0)= r(3);

	return this->setTrialStrain(newStrain);

    // return -1 ;
}



//send back the strain
const Vector& MultiaxialCyclicPlasticityAxiSymm :: getStrain( ) 
{
  strain_vec(0) =       strain(0,0) ;
  strain_vec(1) =       strain(1,1) ;
  strain_vec(2) =       strain(2,2) ;

  strain_vec(3) = 2.0 * strain(0,1) ;

  //opserr<<"MCP::getStrain() "<<strain_vec;
  return strain_vec ;
  
} 


//send back the stress 
const Vector& MultiaxialCyclicPlasticityAxiSymm :: getStress( ) 
{
  stress_vec(0) = stress(0,0) ;
  stress_vec(1) = stress(1,1) ;
  stress_vec(2) = stress(2,2) ;

  stress_vec(3) = stress(0,1) ;

  //opserr<<"MCP::getStress()= "<<stress_vec;
  return stress_vec ;
}

//send back the tangent 
const Matrix& MultiaxialCyclicPlasticityAxiSymm :: getTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           2 2   
  //   3           0 1  ( or 1 0 )
  
  int ii, jj ;
  int i, j, k, l ;

  for ( ii = 0; ii < 4; ii++ ) {
    for ( jj = 0; jj < 4; jj++ ) {

      index_map( ii, i, j ) ;
      index_map( jj, k, l ) ;

      tangent_matrix(ii,jj) = tangent[i][j][k][l] ;

    } //end for j
  } //end for i

  //opserr<<"MCP::return tangent"<<endln;
  //opserr<<"MCP:"<<tangent_matrix<<endln;

  return tangent_matrix ;
} 

//send back the tangent 
const Matrix& MultiaxialCyclicPlasticityAxiSymm :: getInitialTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           2 2   
  //   3           0 1  ( or 1 0 )

  this->doInitialTangent();

  int ii, jj ;
  int i, j, k, l ;

  for ( ii = 0; ii < 4; ii++ ) {
    for ( jj = 0; jj < 4; jj++ ) {

      index_map( ii, i, j ) ;
      index_map( jj, k, l ) ;

      tangent_matrix(ii,jj) = initialTangent[i][j][k][l] ;

    } //end for j
  } //end for i

  return tangent_matrix ;
} 

//swap history variables
/* 
int MultiaxialCyclicPlasticityAxiSymm :: commitState( )  
{

  e_p_n        = e_p_nplus1 ;       // dev plastic part
  theta_p_n   = theta_p_nplus1 ;  // add vol plastic part
  xi_n        = xi_nplus1 ;

  return 0;

}


//revert to last saved state
int MultiaxialCyclicPlasticityAxiSymm :: revertToLastCommit( )
{ 
  return 0 ;
} 

//revert to start
int MultiaxialCyclicPlasticityAxiSymm :: revertToStart( ) 

{  
  this->zero( ) ;
  return 0 ;
}

 


int
MultiaxialCyclicPlasticityAxiSymm::sendSelf (int commitTag, Channel &theChannel)
{
  return -1;
}

int
MultiaxialCyclicPlasticityAxiSymm::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  return -1;
}


*/
