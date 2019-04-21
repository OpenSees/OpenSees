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

#include <CPlaneStress.h>
#include <Channel.h>
#include  <FEM_ObjectBroker.h>

//static vectors and matrices
Vector CPlaneStress :: strain_vec(3) ;
Vector CPlaneStress :: stress_vec(3) ;
Matrix CPlaneStress :: tangent_matrix(3,3) ;


//null constructor
CPlaneStress ::  CPlaneStress( ) : 
Concrete( ) 
{  }


//full constructor
CPlaneStress :: CPlaneStress(	int    tag, 
								double E_i,
								double ni_i,
								double f01d_i,
								double f02d_i,
								double f0p_i,
								double beta_i,
								double An_i,
								double Bn_i,
								double Ap_i,
								double gammaC_i,
								double dchemn_i,
								double dchemp_i,   // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
								int    def_i,    
								double dpMax_i,  
								double dnMax_i,
								bool srfCompr_i,  // 11/03/2013 Diego Talledo: Apply SRF also to compression.
								bool isDchemVar_i,  // 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
								double eps_u_i  // 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
				):
Concrete( tag, ND_TAG_CPlaneStress, 
	E_i,ni_i,f01d_i,f02d_i,f0p_i,beta_i,An_i,Bn_i,Ap_i,gammaC_i,dchemn_i,dchemp_i,def_i,dpMax_i,dnMax_i,srfCompr_i, isDchemVar_i, eps_u_i)  //Talledo: da sistemare lo 0
{ 

}

//elastic constructor
CPlaneStress :: 
CPlaneStress(   int    tag, 
                 double E_i, 
                 double ni_i ) :
Concrete( tag, ND_TAG_CPlaneStress, E_i, ni_i )
{ 

}

//destructor
CPlaneStress :: ~CPlaneStress( ) 
{ 

} 

//make a clone of this material
NDMaterial* CPlaneStress :: getCopy( ) 
{ 
	CPlaneStress  *clone;
	clone = new CPlaneStress() ;   //new instance of this class
	*clone = *this ;          //asignment to make copy
	return clone ;
}

//send back type of material
const char* CPlaneStress :: getType( ) const 
{
	return "PlaneStress" ;
}


//send back order of strain in vector form
int CPlaneStress :: getOrder( ) const 
{ 
	return 3 ; 
} 

//get the strain and integrate plasticity equations
int CPlaneStress :: setTrialStrain( const Vector &strain_from_element) 
{
	const double tolerance = 1e-12 ;

	const int max_iterations = 25 ;
	int iteration_counter  = 0 ;

	int i, j, k, l ;
	int ii, jj ;

	double eps22  =  strain(2,2) ;
	strain.Zero( ) ;

	strain(0,0) =        strain_from_element(0) ;
	strain(1,1) =        strain_from_element(1) ;
	strain(0,1) = 0.50 * strain_from_element(2) ; // eps_12 = 0.5 * gamma_12
	strain(1,0) =        strain(0,1) ;

	strain(2,2) =        eps22 ; 

	//enforce the plane stress condition sigma_22 = 0 
	//solve for epsilon_22 
	iteration_counter = 0 ;  
	do {

		this->plastic_integrator( ) ;
    
		strain(2,2) -= stress(2,2) / tangent[2][2][2][2] ;

		iteration_counter++ ;
		if ( iteration_counter > max_iterations ) {
			/*opserr << "More than " << max_iterations ;
			opserr << " iterations in setTrialStrain of CPlaneStress \n"; */
			break ;
		}// end if 

	} while ( fabs(stress(2,2)) > tolerance ) ;

	//modify tangent for plane stress 
	for ( ii = 0; ii < 3; ii++ ) {
		for ( jj = 0; jj < 3; jj++ )  {

			index_map( ii, i, j ) ;
			index_map( jj, k, l ) ;

			tangent[i][j][k][l] -=   tangent[i][j][2][2] 
									* tangent[2][2][k][l] 
									/ tangent[2][2][2][2] ;

			//minor symmetries 
			tangent [j][i][k][l] = tangent[i][j][k][l] ;
			tangent [i][j][l][k] = tangent[i][j][k][l] ;
			tangent [j][i][l][k] = tangent[i][j][k][l] ;

		} // end for jj
	} // end for ii 
	return 0 ;
}

//unused trial strain functions
int CPlaneStress :: setTrialStrain( const Vector &v, const Vector &r )
{ 
	return this->setTrialStrain( v ) ;
} 

int CPlaneStress :: setTrialStrainIncr( const Vector &v ) 
{
	static Vector newStrain(3);
	newStrain(0) = strain(0,0) + v(0);
	newStrain(1) = strain(1,1) + v(1);
	newStrain(2) = 2.0 * strain(0,1) + v(2);

	return this->setTrialStrain(newStrain);  
}

int CPlaneStress :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
	return this->setTrialStrainIncr(v);
}

//send back the strain
const Vector& CPlaneStress :: getStrain( ) 
{
	strain_vec(0) =       strain(0,0) ;
	strain_vec(1) =       strain(1,1) ;
	strain_vec(2) = 2.0 * strain(0,1) ;

	return strain_vec ;
} 

//send back the stress 
const Vector& CPlaneStress :: getStress( ) 
{
	stress_vec(0) = stress(0,0) ;
	stress_vec(1) = stress(1,1) ;
	stress_vec(2) = stress(0,1) ;

	return stress_vec ;
}

//send back the tangent 
const Matrix& CPlaneStress :: getTangent( ) 
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
const Matrix& CPlaneStress :: getInitialTangent( ) 
{
	// matrix to tensor mapping
	//  Matrix      Tensor
	// -------     -------
	//   0           0 0
	//   1           1 1
	//   2           0 1  ( or 1 0 ) 
	// 

	/* // modifica del 26/11/2012
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
	*/
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



int 
CPlaneStress::commitState( ) 
{
	commitEps22 = strain(2,2);
	epsilon_p_n = epsilon_p_nplus1 ;
	rnp=rnp1p;
	rnn=rnp1n;

	strain_n = strain;
	stress_n = stress;

	return 0;
}

int 
CPlaneStress::revertToLastCommit( ) {
	strain(2,2)  = commitEps22;

	return 0;
}


int 
CPlaneStress::revertToStart( ) 
{

	commitEps22 = 0.0;
  
	this->zero( ) ;

	return 0;
}

int
CPlaneStress::sendSelf (int commitTag, Channel &theChannel)
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
		opserr << "CPlaneStress::sendSelf - failed to send vector to channel\n";
		return -1;
	}
	return 0;
}

int
CPlaneStress::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{

  // recv the vector object from the channel which defines material param and state
	static Vector data(14);
	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "CPlaneStress::recvSelf - failed to recv vector from channel\n";
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

//matrix_index ---> tensor indices i,j
// plane stress different because of condensation on tangent
// case 3 switched to 1-2 and case 4 to 3-3 
void 
CPlaneStress :: index_map( int matrix_index, int &i, int &j )
{
	switch ( matrix_index+1 ) { //add 1 for standard tensor indices

		case 1 :
			i = 1 ; 
			j = 1 ;
			break ;
 
		case 2 :
			i = 2 ;
			j = 2 ; 
			break ;

		case 3 :
			i = 1 ;
			j = 2 ;
			break ;

		case 4 :
			i = 3 ;
			j = 3 ;
			break ;

		case 5 :
			i = 2 ;
			j = 3 ;
			break ;

		case 6 :
			i = 3 ;
			j = 1 ;
			break ;

		default :
			i = 1 ;
			j = 1 ;
			break ;
	
	} //end switch

	i-- ; //subtract 1 for C-indexing
	j-- ;

	return ; 
}







