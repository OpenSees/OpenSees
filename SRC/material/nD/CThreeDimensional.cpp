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
                                                                        
// $Revision: 1.6 $
// $Date: 2005/03/25 00:34:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/CThreeDimensional.cpp,v $

//
//  Send strains in following format :
// 
//     strain_vec = {   eps_00
//                      eps_11
//	                	eps_22
 		      
//                    2 eps_01   
//            	      2 eps_12   
//		              2 eps_20    }   <--- note the 2
// 

#include <CThreeDimensional.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static vectors and matrices
Vector CThreeDimensional :: strain_vec(6) ;
Vector CThreeDimensional :: stress_vec(6) ;
Matrix CThreeDimensional :: tangent_matrix(6,6) ;

//null constructor
CThreeDimensional ::  CThreeDimensional( ) : 
Concrete( ) 
{  }

//full constructor
CThreeDimensional :: 
CThreeDimensional(	int    tag,
					double E_i,
					double ni_i,
					double f01d_i,
					double f02d_i,
					double f0p_i,
					double beta_i,
					double An_i,
					double Bn_i,
					double Ap_i,
					double gammaC_i, // 14/12/2011 Diego Talledo: Added Shear Retention Factor
					double dchem_i,  // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
					double dchemp_i,  // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
					int    def_i,    // 10/01/2013 Diego Talledo: Added different equivalent tension definitions
					double dpMax_i,  // 11/01/2013 Diego Talledo: Added Positive Damage Limit
					double dnMax_i,  // 11/01/2013 Diego Talledo: Added Negative Damage Limit
					bool srfCompr_i, // 11/03/2013 Diego Talledo: Apply SRF also to compression.
					bool isDchemVar_i,  // 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
					double eps_u_i  // 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
				):
Concrete( tag, ND_TAG_CThreeDimensional, 
             E_i,ni_i,f01d_i,f02d_i,f0p_i,beta_i,An_i,Bn_i,Ap_i,gammaC_i, dchem_i, dchemp_i, def_i, dpMax_i, dnMax_i, srfCompr_i, isDchemVar_i, eps_u_i )
{ 

}

//elastic constructor
CThreeDimensional :: 
CThreeDimensional(   int    tag, 
                 double E_i, 
                 double ni_i ) :
Concrete( tag, ND_TAG_CThreeDimensional, E_i,ni_i )
{ 

}

//destructor
CThreeDimensional :: ~CThreeDimensional( ) 
{ 

} 

//make a clone of this material
NDMaterial* CThreeDimensional :: getCopy( ) 
{ 
	CThreeDimensional  *clone;
	clone = new CThreeDimensional( ) ;   //new instance of this class
	*clone = *this ;          //asignment to make copy
	return clone ;
}


//send back type of material
const char* CThreeDimensional :: getType( ) const 
{
	return "ThreeDimensional" ;
}


//send back order of strain in vector form
int CThreeDimensional :: getOrder( ) const 
{ 
	return 6 ; 
} 


//get the strain and integrate plasticity equations
int CThreeDimensional :: setTrialStrain( const Vector &strain_from_element) 
{
	strain.Zero( ) ;

	strain(0,0) =        strain_from_element(0) ;
	strain(1,1) =        strain_from_element(1) ;
	strain(2,2) =        strain_from_element(2) ;

	strain(0,1) = 0.50 * strain_from_element(3) ;
	strain(1,0) =        strain(0,1) ;

	strain(1,2) = 0.50 * strain_from_element(4) ;
	strain(2,1) =        strain(1,2) ;
  
	strain(2,0) = 0.50 * strain_from_element(5) ;
	strain(0,2) =        strain(2,0) ;

	this->plastic_integrator( ) ;

	return 0 ;
}


//unused trial strain functions
int CThreeDimensional :: setTrialStrain( const Vector &v, const Vector &r )
{ 
	return this->setTrialStrain( v ) ;
} 

int CThreeDimensional :: setTrialStrainIncr( const Vector &v ) 
{
	static Vector newStrain(6);
	newStrain(0) = strain(0,0) + v(0);
	newStrain(1) = strain(1,1) + v(1);
	newStrain(2) = strain(2,2) + v(2);
	newStrain(3) = 2.0*strain(0,1) + v(3);
	newStrain(4) = 2.0*strain(1,2) + v(4);
	newStrain(5) = 2.0*strain(2,0) + v(5);
  
	return this->setTrialStrain(newStrain);
}

int CThreeDimensional :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
	return this->setTrialStrainIncr(v);
}

//send back the strain
const Vector& CThreeDimensional :: getStrain( ) 
{
	strain_vec(0) =       strain(0,0) ;
	strain_vec(1) =       strain(1,1) ;
	strain_vec(2) =       strain(2,2) ;

	strain_vec(3) = 2.0 * strain(0,1) ;

	strain_vec(4) = 2.0 * strain(1,2) ;

	strain_vec(5) = 2.0 * strain(2,0) ;

	return strain_vec ;
} 

//send back the stress 
const Vector& CThreeDimensional :: getStress( ) 
{
	stress_vec(0) = stress(0,0) ;
	stress_vec(1) = stress(1,1) ;
	stress_vec(2) = stress(2,2) ;

	stress_vec(3) = stress(0,1) ;

	stress_vec(4) = stress(1,2) ;
  
	stress_vec(5) = stress(2,0) ;

	return stress_vec ;
}

//send back the tangent 
const Matrix& CThreeDimensional :: getTangent( ) 
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
const Matrix& CThreeDimensional :: getInitialTangent( ) 
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