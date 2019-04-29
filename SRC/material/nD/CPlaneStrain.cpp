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
// $Date: 2008/10/20 22:23:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2PlaneStrain.cpp,v $

//
//  Send strains in following format :
// 
//     strain_vec = {   eps_00
//                      eps_11
//                    2 eps_01   }   <--- note the 2
// 
//  set eta := 0 for rate independent case
//

#include <CPlaneStrain.h>
#include <Channel.h>
#include  <FEM_ObjectBroker.h>

//static vectors and matrices
Vector CPlaneStrain :: strain_vec(3) ;
Vector CPlaneStrain :: stress_vec(3) ;
Matrix CPlaneStrain :: tangent_matrix(3,3) ;

//null constructor
CPlaneStrain ::  CPlaneStrain( ) : 
Concrete( ) 
{  }

//full constructor
CPlaneStrain :: CPlaneStrain(	int    tag, 
								double E_i,
								double ni_i,
								double f01d_i,
								double f02d_i,
								double f0p_i,
								double beta_i,
								double An_i,
								double Bn_i,
								double Ap_i,
								double GammaC_i,
								double dchemn_i,   // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
								double dchemp_i,   // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
								int    def_i,    
								double dpMax_i,  
								double dnMax_i,
								bool srfCompr_i,  // 11/03/2013 Diego Talledo: Apply SRF also to compression.
								bool isDchemVar_i,  // 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
								double eps_u_i  // 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
				):
Concrete( tag, ND_TAG_CPlaneStrain, 
	E_i,ni_i,f01d_i,f02d_i,f0p_i,beta_i,An_i,Bn_i,Ap_i, GammaC_i, dchemn_i, dchemp_i, def_i, dpMax_i, dnMax_i,srfCompr_i,isDchemVar_i, eps_u_i)  // Talledo: da sistemare lo 0
{ 

}

//elastic constructor
CPlaneStrain :: 
CPlaneStrain(	int    tag,
				double E_i,
				double ni_i ) :
Concrete( tag, ND_TAG_CPlaneStrain, E_i, ni_i )
{ 

}

//destructor
CPlaneStrain :: ~CPlaneStrain( ) 
{ 

} 

//make a clone of this material
NDMaterial* CPlaneStrain :: getCopy( ) 
{ 
	CPlaneStrain  *clone;
	clone = new CPlaneStrain() ;   //new instance of this class
	*clone = *this ;          //asignment to make copy
	return clone ;
}

//send back type of material
const char* CPlaneStrain :: getType( ) const 
{
	return "PlaneStrain" ;
}

//send back order of strain in vector form
int CPlaneStrain :: getOrder( ) const 
{ 
	return 3 ; 
}

//get the strain and integrate plasticity equations
int CPlaneStrain :: setTrialStrain( const Vector &strain_from_element) 
{
	strain.Zero( ) ;

	strain(0,0) =        strain_from_element(0) ;
	strain(1,1) =        strain_from_element(1) ;
	strain(0,1) = 0.50 * strain_from_element(2) ;
	strain(1,0) =        strain(0,1) ;

	this->plastic_integrator( ) ;

	return 0 ;
}

//unused trial strain functions
int CPlaneStrain :: setTrialStrain( const Vector &v, const Vector &r )
{ 
	return this->setTrialStrain( v ) ;
} 

int CPlaneStrain :: setTrialStrainIncr( const Vector &v ) 
{
	static Vector newStrain(3);
	newStrain(0) = strain(0,0) + v(0);
	newStrain(1) = strain(1,1) + v(1);
	newStrain(2) = 2.0 * strain(0,1) + v(2);

	return this->setTrialStrain(newStrain);  
}

int CPlaneStrain :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
	return this->setTrialStrainIncr(v);
}

//send back the strain
const Vector& CPlaneStrain :: getStrain( ) 
{
	strain_vec(0) =       strain(0,0) ;
	strain_vec(1) =       strain(1,1) ;
	strain_vec(2) = 2.0 * strain(0,1) ;

	return strain_vec ;
} 

//send back the stress 
const Vector& CPlaneStrain :: getStress( ) 
{
	stress_vec(0) = stress(0,0) ;
	stress_vec(1) = stress(1,1) ;
	stress_vec(2) = stress(0,1) ;

	return stress_vec ;
}

//send back the tangent 
const Matrix& CPlaneStrain :: getTangent( ) 
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
const Matrix& CPlaneStrain :: getInitialTangent( ) 
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
CPlaneStrain::commitState( ) 
{
	epsilon_p_n = epsilon_p_nplus1 ;
	rnp=rnp1p;
	rnn=rnp1n;

	strain_n = strain;
	stress_n = stress;

	return 0;
}

int 
CPlaneStrain::revertToLastCommit( ) {

  return 0;
}


int 
CPlaneStrain::revertToStart( ) 
{
  this->zero( ) ;

  return 0;
}

int
CPlaneStrain::sendSelf (int commitTag, Channel &theChannel)
{
	// we place all the data needed to define material and it's state
	// int a vector object
	static Vector data(14);
	int cnt = 0;
	data(cnt++) = this->getTag();
	data(cnt++) = E    ;
	data(cnt++) = ni   ;
	data(cnt++) = f01d ;
	data(cnt++) = f02d ;
	data(cnt++) = f0p  ;
	data(cnt++) = beta ;
	data(cnt++) = An   ;
	data(cnt++) = Bn   ;
	data(cnt++) = Ap   ;
	data(cnt++) = rnn  ;
	data(cnt++) = rnp  ;
	data(cnt++) = rnp1n; 
	data(cnt++) = rnp1p;


	// send the vector object to the channel
	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "CPlaneStrain::sendSelf - failed to send vector to channel\n";
		return -1;
	}
	return 0;
}

int
CPlaneStrain::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{

	// recv the vector object from the channel which defines material param and state
	static Vector data(14);
	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "CPlaneStrain::recvSelf - failed to recv vector from channel\n";
		return -1;
	}

	// set the material parameters and state variables
	int cnt = 0;
	this->setTag(data(cnt++));
	E    = data(cnt++);
	ni   = data(cnt++);
	f01d = data(cnt++);
	f02d = data(cnt++);
	f0p  = data(cnt++);
	beta = data(cnt++);
	An   = data(cnt++);
	Bn   = data(cnt++);
	Ap   = data(cnt++);
	rnn  = data(cnt++);
	rnp  = data(cnt++);
	rnp1n= data(cnt++);
	rnp1p= data(cnt++);
	//rnp1n=rnn;
	//rnp1p=rnp;
	epsilon_p_nplus1 = epsilon_p_n;

	return 0;
}








