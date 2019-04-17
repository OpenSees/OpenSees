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
                                                                        
// $Revision: 1.20 $
// $Date: 2012/06/06 12:38:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/Concrete.cpp,v $
// Written: Tesser L., Talledo D.A.
//  
/* ****************************************************************** **
**                                                                    **
** Concrete: Implementation of Concrete class as Tesser et al. (2012) **
**                                                                    **
** Two parameters damage model for isotropic material.                **
** Plasticity simplified as Faria et al. (1997)                       **
** Environmental damage as Saetta et al. (1998)                       **
**                                                                    **
** ****************************************************************** */

#include <Concrete.h>
#include <CThreeDimensional.h> 
#include <CPlaneStrain.h>
#include <CPlaneStress.h>
#include <CPlaneStress2d.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

//#include <stresst.h>  // Talledo

#include <fstream>  // Tesser
#include <iostream>  // Talledo
using namespace std;
//using std::cerr;
//using std::ios;

//parameters
const double Concrete :: one3   = 1.0 / 3.0 ;
const double Concrete :: two3   = 2.0 / 3.0 ;
const double Concrete :: four3  = 4.0 / 3.0 ;
const double Concrete :: root23 = sqrt( 2.0 / 3.0 ) ;

double Concrete::initialTangent[3][3][3][3] ;   //material tangent
double Concrete::IIdev[3][3][3][3] ; //rank 4 deviatoric
double Concrete::IbunII[3][3][3][3]; //rank 4 general
double Concrete::IbunI[3][3][3][3] ; //rank 4 I bun I 

#include <elementAPI.h>
#define OPS_Export 
#define PI 3.14159

static int numConcreteMaterials = 0;

void *
OPS_NewConcreteMaterial(void)
{
	if (numConcreteMaterials == 0) {
	  numConcreteMaterials++;
	  opserr << "nDmaterial Concrete - Written by Leopoldo Tesser, Univ. Padua - Italy";
	  opserr << "contact: leopoldo.tesser@dicea.unipd.it\n";
	//to be completed with reference
	}
	
	// Pointer to a NDMaterial material that will be returned
	NDMaterial *theMaterial = 0;
	
	int numArgs = OPS_GetNumRemainingInputArgs();

	// Mandatory parameters
	if (numArgs < 10)
	{
		opserr << "WARNING invalid number of input arguments for nDMaterial Concrete" << endln;
		opserr << "Want: nDMaterial Concrete tag? E? ni? f01d? f02d? f0p? beta? An? Bn? Ap? <-srf gammaC?> <-dchem dchem?> <-eqTensDef def?> <-dpMax max?> <-dnMax max?>" << endln;
		return 0;
	}
	
	int tag = 0;
	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid material tag for nDMaterial Concrete" << endln;
		return 0;
	}

	double dData[9];
	numData = 9;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid material data for nDMaterial Concrete with tag: " << tag << endln;
		return 0;
	}

	char *Param = 0;
	numData = 1;

	// Optional Parameteres
	//default parameters
	// 14/12/2011 Diego Talledo: Added Shear Retention Factor
	double gammaC = -0.001; // (negative by default -> no influence)
	// 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
	double dchemn = 0.0;  // (zero by default -> no influence)
	double dchemp = 0.0;  // (zero by default -> no influence)
	// 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
	bool dchemVariable = false;  // by default dchem is constant
	double eps_u = 0.0; // By default eps_u is not used
	// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
	int def = 2; // // (2 by default -> COMPDYN definition)
	// 11/01/2013 Diego Talledo: Added Negative and Positive Damage Limit
	double dpMax = 0.99999; // (0.99999 by default )
	double dnMax = 0.99999; // (0.99999 by default -> as usal)
	// 11/03/2013 Diego Talledo: Apply SRF also to compression.
	bool srfCompr = false; // Apply only in tension by Default

	for (int ii = 1; ii < numArgs-9; ii++) {
	    const char *Param = OPS_GetString();
	    if (strcmp(Param,"-srf") == 0 && ++ii < numArgs) {
			if (OPS_GetDouble(&numData, &gammaC) != 0) {
				opserr << "WARNING invalid eps_ref for nDMaterial Concrete material with tag: " << tag << endln;
				return 0;		
			}
		} else if ((strcmp(Param,"-dchem") == 0 || strcmp(Param,"-dchemn") == 0) && ++ii < numArgs) {
			if (OPS_GetDouble(&numData, &dchemn) != 0) {
				opserr << "WARNING invalid dchem for nDMaterial Concrete material with tag: " << tag << endln;
				return 0;		
			}
			else dchemp = dchemn;
		} else if (strcmp(Param,"-dchemp") == 0 && ++ii < numArgs) {
			if (OPS_GetDouble(&numData, &dchemp) != 0) {
				opserr << "WARNING invalid dchemp for nDMaterial Concrete material with tag: " << tag << endln;
				return 0;		
			}
		} else if (strcmp(Param,"-dchemvar") == 0 && ++ii < numArgs) { // 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
			if (OPS_GetDouble(&numData, &eps_u) != 0) {
				opserr << "WARNING invalid eps_u for nDMaterial Concrete material with tag: " << tag << endln;
				return 0;		
			} else dchemVariable = true;
		} else if (strcmp(Param,"-eqTensDef") == 0 && ++ii < numArgs) {
			if (OPS_GetInt(&numData, &def) != 0) {
				opserr << "WARNING invalid number for Eq Tens Definition for nDMaterial Concrete material with tag: " << tag << endln;
				return 0;		
			}
		} else if (strcmp(Param,"-dpMax") == 0 && ++ii < numArgs) {
			if (OPS_GetDouble(&numData, &dpMax) != 0) {
				opserr << "WARNING invalid dpMax for nDMaterial Concrete material with tag: " << tag << endln;
				return 0;		
			}
		} else if (strcmp(Param,"-dnMax") == 0 && ++ii < numArgs) {
			if (OPS_GetDouble(&numData, &dnMax) != 0) {
				opserr << "WARNING invalid dnMax for nDMaterial Concrete material with tag: " << tag << endln;
				return 0;		
			}
		} else if (strcmp(Param,"-srfCompr") == 0) {
			srfCompr = true;
		}
	}
	
	theMaterial = new Concrete(tag,0,dData[0],dData[1],dData[2],dData[3],dData[4],dData[5],dData[6],dData[7],dData[8],gammaC,dchemn,dchemp,def,dpMax,dnMax,srfCompr,dchemVariable,eps_u);	
	
	if (theMaterial == 0) {
		opserr << "WARNING ran out of memory for nDMaterial Concrete with tag: " << tag << endln;
	}

	return theMaterial;
}

//zero internal variables
void Concrete :: zero ( ) 
{
	epsilon_p_n.Zero( ) ;
	epsilon_p_nplus1.Zero( ) ;

	stress.Zero();
	strain.Zero();

	stress_n.Zero();  // Diego Talledo: a cosa servono stress_n e strain_n ?
	strain_n.Zero();
}

//null constructor
Concrete ::  Concrete( ) : 
NDMaterial( ),
epsilon_p_n(3,3),
epsilon_p_nplus1(3,3),
stress(3,3),
strain(3,3),
stress_n(3,3),
strain_n(3,3)
{ 
	E    = 0.0;
	ni   = 0.0;
	f01d = 0.0;
	f02d = 0.0;
	f0p  = 0.0;
	beta = 0.0;
	An   = 0.0;
	Bn   = 0.0;
	Ap   = 0.0;
	gammaC = 0.0; // 14/12/2011 Diego Talledo: Added Shear Retention Factor
  
	lambda  = 0.0;
	shear_p = 0.0;
	bulk_p  = 0.0;
	k_p     = 0.0;
	r0n     = 0.0;
	r0p     = 0.0;

	// 14/12/2011 Diego Talledo: Added Shear Retention Factor
	SRF12 = 0.0;
	SRF12 = 0.0;
	SRF12 = 0.0;
  
	dn  = 0.0;
	dp  = 0.0;
	rnn  = r0n;
	rnp1n= rnn;
	rnp  = r0p;
	rnp1p= rnp;

	// 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
	dchemn = 0.0;
	dchemp = 0.0;
	dpstar = 0.0;
	dnstar = 0.0;

	// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
	usedEquivalentTensionDefinition = COMPDYN; // By default we use the COMPDYN implemented definition

	// 11/01/2013 Diego Talledo: Added Positive and Negative Damage Limit
	dpMax = 0.99999;   // Default value
	dnMax = 0.99999; // Default value

	// 11/03/2013 Diego Talledo: Apply SRF also to compression.
	srfCompr = false; // Default value
	
	this->zero( ) ;     // or (*this).zero( ) 

	int i, j, k, l ;

	//zero rank4 IIdev, IbunI and IbunII 
	for ( i = 0; i < 3; i++ ) {
		for ( j = 0; j < 3; j++ )  {
			for ( k = 0; k < 3; k++ ) {
				for ( l = 0; l < 3; l++)  { 
					IbunI[i][j][k][l] = 0.0 ;
					IbunII[i][j][k][l] = 0.0;
					IIdev[i][j][k][l] = 0.0 ;
				} // end for l
			} // end for k
		} // end for j
	} // end for i

	//form rank4 IbunI 

	IbunI [0][0] [0][0] = 1.0 ;
	IbunI [0][0] [1][1] = 1.0 ;
	IbunI [0][0] [2][2] = 1.0 ;
	IbunI [1][1] [0][0] = 1.0 ;
	IbunI [1][1] [1][1] = 1.0 ;
	IbunI [1][1] [2][2] = 1.0 ;
	IbunI [2][2] [0][0] = 1.0 ;
	IbunI [2][2] [1][1] = 1.0 ;
	IbunI [2][2] [2][2] = 1.0 ;

	//form rank4 IIdev

	IIdev [0][0] [0][0] =  two3 ; // 0.666667 
	IIdev [0][0] [1][1] = -one3 ; //-0.333333 
	IIdev [0][0] [2][2] = -one3 ; //-0.333333 
	IIdev [0][1] [0][1] = 0.5 ;
	IIdev [0][1] [1][0] = 0.5 ;
	IIdev [0][2] [0][2] = 0.5 ;
	IIdev [0][2] [2][0] = 0.5 ;
	IIdev [1][0] [0][1] = 0.5 ;
	IIdev [1][0] [1][0] = 0.5 ;
	IIdev [1][1] [0][0] = -one3 ; //-0.333333 
	IIdev [1][1] [1][1] =  two3 ; // 0.666667 
	IIdev [1][1] [2][2] = -one3 ; //-0.333333 
	IIdev [1][2] [1][2] = 0.5 ;
	IIdev [1][2] [2][1] = 0.5 ;
	IIdev [2][0] [0][2] = 0.5 ;
	IIdev [2][0] [2][0] = 0.5 ;
	IIdev [2][1] [1][2] = 0.5 ;
	IIdev [2][1] [2][1] = 0.5 ;
	IIdev [2][2] [0][0] = -one3 ; //-0.333333 
	IIdev [2][2] [1][1] = -one3 ; //-0.333333 
	IIdev [2][2] [2][2] =  two3 ; // 0.666667 

	//form rank4 IbunII

	IbunII [0][0] [0][0] = 1 ;
	IbunII [0][0] [1][1] = 1 ;  
	IbunII [0][0] [2][2] = 1 ; 
	IbunII [0][1] [0][1] = 1 ;
	IbunII [0][1] [1][0] = 1 ;
	IbunII [0][2] [0][2] = 1 ;
	IbunII [0][2] [2][0] = 1 ;
	IbunII [1][0] [0][1] = 1 ;
	IbunII [1][0] [1][0] = 1 ;
	IbunII [1][1] [0][0] = 1 ; 
	IbunII [1][1] [1][1] = 1 ; 
	IbunII [1][1] [2][2] = 1 ; 
	IbunII [1][2] [1][2] = 1 ;
	IbunII [1][2] [2][1] = 1 ;
	IbunII [2][0] [0][2] = 1 ;
	IbunII [2][0] [2][0] = 1 ;
	IbunII [2][1] [1][2] = 1 ;
	IbunII [2][1] [2][1] = 1 ;
	IbunII [2][2] [0][0] = 1 ; 
	IbunII [2][2] [1][1] = 1 ; 
	IbunII [2][2] [2][2] = 1 ; 
	
	plastic_integrator();
}

//full constructor
Concrete :: Concrete(	int      tag,
						int classTag,
						double E_i,
						double ni_i,
						double f01d_i,
						double f02d_i,
						double f0p_i,
						double beta_i,
						double An_i,
						double Bn_i,
						double Ap_i,
						double gammaC_i,    // 14/12/2011 Diego Talledo: Added Shear Retention Factor
						double dchemn_i,    // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
						double dchemp_i,    // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
						int    def_i,       // 10/01/2013 Diego Talledo: Added different equivalent tension definitions
						double dpMax_i,     // 11/01/2013 Diego Talledo: Added Positive Damage Limit
						double dnMax_i,     // 11/01/2013 Diego Talledo: Added Negative Damage Limit
						bool srfCompr_i,    // 11/03/2013 Diego Talledo: Apply SRF also to compression.
						bool isDchemVar_i,  // 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
						double eps_u_i      // 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
						) 
:
NDMaterial(tag, classTag),
epsilon_p_n(3,3),
epsilon_p_nplus1(3,3),
stress(3,3),
strain(3,3),
stress_n(3,3),
strain_n(3,3)
{
	// Set INPUT Parameters

	E    = E_i;
	ni   = ni_i;
	f01d = f01d_i;
	f02d = f02d_i;
	f0p  = f0p_i;
	beta = beta_i;
	An   = An_i;
	Bn   = Bn_i;
	Ap   = Ap_i;
	gammaC = gammaC_i; // 14/12/2011 Diego Talledo: Added Shear Retention Factor
	// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
	if (def_i == 1)
		usedEquivalentTensionDefinition = HOMOGENEOUS;
	else if (def_i == 2)
		usedEquivalentTensionDefinition = COMPDYN;
	else
		usedEquivalentTensionDefinition = ORIGINAL;

	// Calculate initial Parameters
	lambda  = E*ni/((1+ni)*(1-2*ni));
	shear_p = E/(2*(1+ni));
	bulk_p  = lambda+2/3.0*shear_p;
	k_p     = sqrt(2.0)*(f02d-f01d)/(2*f02d-f01d);
	if (usedEquivalentTensionDefinition == ORIGINAL)
	{
		r0n     = sqrt(sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3);
		r0p     = sqrt(f0p/sqrt(E));
	}
	else if (usedEquivalentTensionDefinition == HOMOGENEOUS)
	{
		r0n     = sqrt(sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3);
		r0p     = sqrt(f0p);
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		r0n     = sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3;
		r0p     = f0p;
	}

	// Environmental damage
	dchemn = dchemn_i; // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
	dchemp = dchemp_i; // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
	eps_u = eps_u_i;  // 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
	isDchemVariableSelected = isDchemVar_i;  // 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
	if (isDchemVariableSelected) {
		if ((usedEquivalentTensionDefinition == ORIGINAL) || (usedEquivalentTensionDefinition == HOMOGENEOUS))
		{
			tau_n_eps_u = sqrt(sqrt(3.0)*(k_p-sqrt(2.0))*(E*(eps_u-((eps_u - f01d/E) * beta)))/3);
		}
		else if (usedEquivalentTensionDefinition == COMPDYN)
		{
			tau_n_eps_u = sqrt(3.0)*(k_p-sqrt(2.0))*(E*(eps_u-((eps_u - f01d/E) * beta)))/3;
		}
	}

	// Initialization of damage parameters
	dn  = 0.0;
	dp  = 0.0;
	rnn = r0n;
	rnp1n=rnn;
	rnp = r0p;
	rnp1p=rnp;
	// 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
	dnstar = 0.0;
	dpstar = 0.0;
	// 14/12/2011 Diego Talledo: Added Shear Retention Factor
	SRF12 = 0.0;
	SRF23 = 0.0;
	SRF13 = 0.0;

	// 11/01/2013 Diego Talledo: Added Positive and Negative Damage Limit
	dpMax = dpMax_i;
	dnMax = dnMax_i;
	
	// 11/03/2013 Diego Talledo: Apply SRF also to compression.
	srfCompr = srfCompr_i;

	this->zero( ) ;     // or (*this).zero( ) 

	int i, j, k, l ;

	//zero rank4 IIdev, IbunI and IbunII 
	for ( i = 0; i < 3; i++ ) {
		for ( j = 0; j < 3; j++ )  {
			for ( k = 0; k < 3; k++ ) {
				for ( l = 0; l < 3; l++)  { 
					IbunI[i][j][k][l] = 0.0 ;
					IbunII[i][j][k][l] = 0.0;
					IIdev[i][j][k][l] = 0.0 ;
				} // end for l
			} // end for k
		} // end for j
	} // end for i

	//form rank4 IbunI 

	IbunI [0][0] [0][0] = 1.0 ;
	IbunI [0][0] [1][1] = 1.0 ;
	IbunI [0][0] [2][2] = 1.0 ;
	IbunI [1][1] [0][0] = 1.0 ;
	IbunI [1][1] [1][1] = 1.0 ;
	IbunI [1][1] [2][2] = 1.0 ;
	IbunI [2][2] [0][0] = 1.0 ;
	IbunI [2][2] [1][1] = 1.0 ;
	IbunI [2][2] [2][2] = 1.0 ;

	//form rank4 IIdev

	IIdev [0][0] [0][0] =  two3 ; // 0.666667 
	IIdev [0][0] [1][1] = -one3 ; //-0.333333 
	IIdev [0][0] [2][2] = -one3 ; //-0.333333 
	IIdev [0][1] [0][1] = 0.5 ;
	IIdev [0][1] [1][0] = 0.5 ;
	IIdev [0][2] [0][2] = 0.5 ;
	IIdev [0][2] [2][0] = 0.5 ;
	IIdev [1][0] [0][1] = 0.5 ;
	IIdev [1][0] [1][0] = 0.5 ;
	IIdev [1][1] [0][0] = -one3 ; //-0.333333 
	IIdev [1][1] [1][1] =  two3 ; // 0.666667 
	IIdev [1][1] [2][2] = -one3 ; //-0.333333 
	IIdev [1][2] [1][2] = 0.5 ;
	IIdev [1][2] [2][1] = 0.5 ;
	IIdev [2][0] [0][2] = 0.5 ;
	IIdev [2][0] [2][0] = 0.5 ;
	IIdev [2][1] [1][2] = 0.5 ;
	IIdev [2][1] [2][1] = 0.5 ;
	IIdev [2][2] [0][0] = -one3 ; //-0.333333 
	IIdev [2][2] [1][1] = -one3 ; //-0.333333 
	IIdev [2][2] [2][2] =  two3 ; // 0.666667 

	//form rank4 IbunII

	IbunII [0][0] [0][0] = 1 ;
	IbunII [0][0] [1][1] = 1 ;  
	IbunII [0][0] [2][2] = 1 ; 
	IbunII [0][1] [0][1] = 1 ;
	IbunII [0][1] [1][0] = 1 ;
	IbunII [0][2] [0][2] = 1 ;
	IbunII [0][2] [2][0] = 1 ;
	IbunII [1][0] [0][1] = 1 ;
	IbunII [1][0] [1][0] = 1 ;
	IbunII [1][1] [0][0] = 1 ; 
	IbunII [1][1] [1][1] = 1 ; 
	IbunII [1][1] [2][2] = 1 ; 
	IbunII [1][2] [1][2] = 1 ;
	IbunII [1][2] [2][1] = 1 ;
	IbunII [2][0] [0][2] = 1 ;
	IbunII [2][0] [2][0] = 1 ;
	IbunII [2][1] [1][2] = 1 ;
	IbunII [2][1] [2][1] = 1 ;
	IbunII [2][2] [0][0] = 1 ; 
	IbunII [2][2] [1][1] = 1 ; 
	IbunII [2][2] [2][2] = 1 ; 
	
	plastic_integrator();
}

//elastic constructor
Concrete :: Concrete(	int    tag, 
						int  classTag,
						double E_i, 
						double ni_i ) :
NDMaterial(tag, classTag),
epsilon_p_n(3,3),
epsilon_p_nplus1(3,3),
stress(3,3),
strain(3,3),
stress_n(3,3),
strain_n(3,3)
{

	// Set INPUT Parameters

	E    = E_i;
	ni   = ni_i;
	f01d = -1.0e16*E;
	f02d = -1.0e16*E;
	f0p  = 1.0e16*E;
	beta = 0.0;
	An   = 1.0;
	Bn   = 1.0;
	Ap   = 1.0;
	gammaC = 0.0; // 14/12/2011 Diego Talledo: Added Shear Retention Factor
	dchemn = 0.0; // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
	dchemp = dchemn;
	// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
	usedEquivalentTensionDefinition = COMPDYN; // By default we use the COMPDYN implemented definition

	// Calculate initial Parameters
	lambda  = E*ni/((1+ni)*(1-2*ni));
	shear_p = E/(2*(1+ni));
	bulk_p  = lambda+2/3.0*shear_p;
	k_p     = sqrt(2.0)*(f02d-f01d)/(2*f02d-f01d);
	if (usedEquivalentTensionDefinition == ORIGINAL)
	{
		r0n     = sqrt(sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3);
		r0p     = sqrt(f0p/sqrt(E));
	}
	else if (usedEquivalentTensionDefinition == HOMOGENEOUS)
	{
		r0n     = sqrt(sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3);
		r0p     = sqrt(f0p);
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		r0n     = sqrt(3.0)*(k_p-sqrt(2.0))*f01d/3;
		r0p     = f0p;
	}

	// Initialization of damage parameters
	dn  = 0.0;
	dp  = 0.0;
	rnn = r0n;
	rnp1n=rnn;
	rnp = r0p;
	rnp1p=rnp;
	// 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
	dnstar = 0.0;
	dpstar = 0.0;
	// 14/12/2011 Diego Talledo: Added Shear Retention Factor
	SRF12 = 0.0;
	SRF23 = 0.0;
	SRF13 = 0.0;

	// 11/01/2013 Diego Talledo: Added Positive and Negative Damage Limit
	dpMax = 0.99999;   // Default value
	dnMax = 0.99999; // Default value

	// 11/03/2013 Diego Talledo: Apply SRF also to compression.
	srfCompr = false;  // Default value
	
	this->zero( ) ;     // or (*this).zero( ) 

	int i, j, k, l ;

	//zero rank4 IIdev, IbunI and IbunII 
	for ( i = 0; i < 3; i++ ) {
		for ( j = 0; j < 3; j++ )  {
			for ( k = 0; k < 3; k++ ) {
				for ( l = 0; l < 3; l++)  { 
					IbunI[i][j][k][l] = 0.0 ;
					IbunII[i][j][k][l] = 0.0;
					IIdev[i][j][k][l] = 0.0 ;
				} // end for l
			} // end for k
		} // end for j
	} // end for i

	//form rank4 IbunI 

	IbunI [0][0] [0][0] = 1.0 ;
	IbunI [0][0] [1][1] = 1.0 ;
	IbunI [0][0] [2][2] = 1.0 ;
	IbunI [1][1] [0][0] = 1.0 ;
	IbunI [1][1] [1][1] = 1.0 ;
	IbunI [1][1] [2][2] = 1.0 ;
	IbunI [2][2] [0][0] = 1.0 ;
	IbunI [2][2] [1][1] = 1.0 ;
	IbunI [2][2] [2][2] = 1.0 ;

	//form rank4 IIdev

	IIdev [0][0] [0][0] =  two3 ; // 0.666667 
	IIdev [0][0] [1][1] = -one3 ; //-0.333333 
	IIdev [0][0] [2][2] = -one3 ; //-0.333333 
	IIdev [0][1] [0][1] = 0.5 ;
	IIdev [0][1] [1][0] = 0.5 ;
	IIdev [0][2] [0][2] = 0.5 ;
	IIdev [0][2] [2][0] = 0.5 ;
	IIdev [1][0] [0][1] = 0.5 ;
	IIdev [1][0] [1][0] = 0.5 ;
	IIdev [1][1] [0][0] = -one3 ; //-0.333333 
	IIdev [1][1] [1][1] =  two3 ; // 0.666667 
	IIdev [1][1] [2][2] = -one3 ; //-0.333333 
	IIdev [1][2] [1][2] = 0.5 ;
	IIdev [1][2] [2][1] = 0.5 ;
	IIdev [2][0] [0][2] = 0.5 ;
	IIdev [2][0] [2][0] = 0.5 ;
	IIdev [2][1] [1][2] = 0.5 ;
	IIdev [2][1] [2][1] = 0.5 ;
	IIdev [2][2] [0][0] = -one3 ; //-0.333333 
	IIdev [2][2] [1][1] = -one3 ; //-0.333333 
	IIdev [2][2] [2][2] =  two3 ; // 0.666667 

	//form rank4 IbunII

	IbunII [0][0] [0][0] = 1 ;
	IbunII [0][0] [1][1] = 1 ;  
	IbunII [0][0] [2][2] = 1 ; 
	IbunII [0][1] [0][1] = 1 ;
	IbunII [0][1] [1][0] = 1 ;
	IbunII [0][2] [0][2] = 1 ;
	IbunII [0][2] [2][0] = 1 ;
	IbunII [1][0] [0][1] = 1 ;
	IbunII [1][0] [1][0] = 1 ;
	IbunII [1][1] [0][0] = 1 ; 
	IbunII [1][1] [1][1] = 1 ; 
	IbunII [1][1] [2][2] = 1 ; 
	IbunII [1][2] [1][2] = 1 ;
	IbunII [1][2] [2][1] = 1 ;
	IbunII [2][0] [0][2] = 1 ;
	IbunII [2][0] [2][0] = 1 ;
	IbunII [2][1] [1][2] = 1 ;
	IbunII [2][1] [2][1] = 1 ;
	IbunII [2][2] [0][0] = 1 ; 
	IbunII [2][2] [1][1] = 1 ; 
	IbunII [2][2] [2][2] = 1 ; 
}

//destructor
Concrete :: ~Concrete( ) 
{  } 

NDMaterial*
Concrete :: getCopy (const char *type)
{
	if (strcmp(type,"PlaneStress") == 0)
	{
		CPlaneStress  *clone ;
		clone = new CPlaneStress(this->getTag(), E, ni, f01d, f02d,
			f0p, beta, An, Bn, Ap, gammaC, dchemn, dchemp, usedEquivalentTensionDefinition, dpMax, dnMax, srfCompr, isDchemVariableSelected, eps_u) ;
		return clone ;
	}
	else if (strcmp(type,"PlaneStress2D") == 0)
	{
		CPlaneStress2d  *clone ;
		clone = new CPlaneStress2d(this->getTag(), E, ni, f01d, f02d,
			f0p, beta, An, Bn, Ap, gammaC, dchemn, dchemp, usedEquivalentTensionDefinition, dpMax, dnMax, srfCompr, isDchemVariableSelected, eps_u) ;
		return clone ;
	}
	else if (strcmp(type,"PlaneStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0)
	{
		CPlaneStrain  *clone ;
		clone = new CPlaneStrain(this->getTag(), E, ni, f01d, f02d,
			f0p, beta, An, Bn, Ap, gammaC, dchemn, dchemp, usedEquivalentTensionDefinition, dpMax, dnMax, srfCompr, isDchemVariableSelected, eps_u) ;
		return clone ;
	}
	else if (strcmp(type,"AxiSymmetric2D") == 0 || strcmp(type,"AxiSymmetric") == 0)
	{
		opserr << "Concrete::getModel failed to get model: " << type << endln;
		return 0;
	}
    else if ((strcmp(type,"ThreeDimensional") == 0) || (strcmp(type,"3D") == 0))
	{
		CThreeDimensional  *clone ;
		clone = new CThreeDimensional(this->getTag(), E, ni, f01d, f02d,
			f0p, beta, An, Bn, Ap, gammaC,  dchemn, dchemp, usedEquivalentTensionDefinition, dpMax, dnMax, srfCompr, isDchemVariableSelected, eps_u) ;
		return clone ;	
    }
	else if ((strcmp(type,"PlateFiber") == 0))
	{
		opserr << "Concrete::getModel failed to get model: " << type << endln;
		return 0;
	}
    // Handle other cases
	else
	{
		opserr << "Concrete::getModel failed to get model: " << type << endln;
		return 0;
	}
}

//print out material data
void Concrete :: Print( OPS_Stream &s, int flag )
{
	s << endln ;
	s << "Concrete : " ; 
	s << this->getType( ) << endln ;
	s << "Young Modulus   =  " << E    << endln ;
	s << "Poisson Modulus =  " << ni   << endln ;
	s << "f01d            =  " << f01d << endln ;
	s << "f02d            =  " << f02d << endln ;
	s << "f0p             =  " << f0p  << endln ;
	s << "beta            =  " << beta << endln ;
	s << "An              =  " << An   << endln ;
	s << "Bn              =  " << Bn   << endln ;
	s << "Ap              =  " << Ap   << endln ;
	// 14/12/2011 Diego Talledo: Added Shear Retention Factor
	s << "gammaC			=  " << gammaC << endln;
	// 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
	s << "dchemn          = " << dchemn << endln;
	s << "dchemp          = " << dchemp << endln;
	if (isDchemVariableSelected) 
	{
		s << "dchem linearly variable " << endln;
		s << "eps_u       = " << eps_u << endln;
	}
	// 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
	// 10/01/2013 Diego Talledo: Added different equivalent tension definitions
	if (usedEquivalentTensionDefinition == ORIGINAL)
	{
		s << "Using original equivalent tension definition ("; 
	}
	else if (usedEquivalentTensionDefinition == HOMOGENEOUS)
	{
		s << "Using homogeneous equivalent tension definition ("; 
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		s << "Using compdyn equivalent tension definition ("; 
	}
	s << usedEquivalentTensionDefinition << ")" << endln;
}

//--------------------Plasticity-------------------------------------

//plasticity integration routine
void Concrete :: plastic_integrator( )
{
	const double tolerance = (1.0e-14)*f0p ;

	static Matrix dev_strain(3,3) ; //deviatoric strain
	static Matrix dev_stress(3,3) ; //deviatoric stress

	static Matrix tmp_senp1t(3,3) ; //Sigma_eff n+1 trial

	double trace = 0.0 ; //trace of strain

	int i,j,k,l;
	int ii, jj ; 

	double senp1t[3][3];
	double V[3][3],EE[3][3];
	double dg[3],dgn[3],dgp[3];
	double sigoct,tauoct,taun,taup;
	double alfa,diag,norms,pint,lambdap;
	double st[3][3],dea[3][3],ls[3][3],ses[3][3];
	double trans[3][3];
	double supp1[3][3],supp2[3][3];

	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
			{
				if (i==j)
					EE[i][j]=1.0;
				else
					EE[i][j]=0.0;
			};
	};

	double g, rhoQ, thetaQ, rhoP, alfasp, alfasn, rhoL, thetaL;

	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
			dep[i][j]=0.0;
	};
	
	// evaluate trace of elastic strain
	trace=0.0;
	for (i=0;i<3;i++)
    {
		trace+=strain(i,i)-epsilon_p_n(i,i);
    };

	// Deviatoric strain tensor
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
			dev_strain(i,j)=strain(i,j)-epsilon_p_n(i,j);
    };
	for (i=0;i<3;i++)
    {
		dev_strain(i,i)-= 1*trace/3;
    };

	// Deviatoric Stress
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
			dev_stress(i,j)=dev_strain(i,j)*2.0*shear_p;  // Strain(1,2) è epsilon12, per 2 = gamma12
    };

	// Sigma Eff (n+1,Trial) = DevStress + I * K * trace(StrainElastic)
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
			senp1t[i][j]=dev_stress(i,j);
    };
	for (i=0;i<3;i++)
		senp1t[i][i]=senp1t[i][i]+bulk_p*trace;

    if (beta>0.0)
	{
		// Splitting of Trial Effective Stress in Principal System

		// Project in principal system
		eigen_decomposition(senp1t,EE,V,dg);

		// Negative part of principal Trial Effective Stress
		for (i=0;i<3;i++)
		    dgn[i]=(dg[i]-fabs(dg[i]))/2;

		// Compute tau,n
		sigoct=(dgn[0]+dgn[1]+dgn[2])/3;
		tauoct=sqrt((dgn[0]-dgn[1])*(dgn[0]-dgn[1])+(dgn[0]-dgn[2])*(dgn[0]-dgn[2])+(dgn[1]-dgn[2])*(dgn[1]-dgn[2]))/3;
		if ((usedEquivalentTensionDefinition == ORIGINAL)  || (usedEquivalentTensionDefinition == HOMOGENEOUS))
		{
			taun=sqrt(3.0)*(k_p*sigoct+tauoct);
			if (taun >= 0)
				taun = sqrt(taun);
			else
				taun = 0;
		}
		else if (usedEquivalentTensionDefinition == COMPDYN)
		{
			taun=sqrt(3.0)*(k_p*sigoct+tauoct);
		}
		
		// Positive part of principal Trial Effective Stress
		for (i=0;i<3;i++)
			dgp[i]=(dg[i]+fabs(dg[i]))/2; 

		// Compute tau,p
        for (i=0;i<3;i++)
			dgp[i]=(dg[i]+fabs(dg[i]))/2; 
		diag=(dgp[0]+dgp[1]+dgp[2])*ni/(-E);
		for (i=0;i<3;i++)
		{
			for (j=0;j<3;j++)
			{
				trans[i][j]=0.0;
				if (i==j)
					trans[i][i]=(dgp[i]*(1+ni)/E)+diag;
			};
		};
		taup=0.0;
		for (i=0;i<3;i++)
		{
			for (j=0;j<3;j++)
				if (i==j)
					taup+=trans[i][j]*dgp[i];
		};
		if (usedEquivalentTensionDefinition == ORIGINAL)
		{
			taup=sqrt(sqrt(taup));
		}
		else if (usedEquivalentTensionDefinition == HOMOGENEOUS)
		{
			taup=sqrt(sqrt(taup*E));
		}
		else if (usedEquivalentTensionDefinition == COMPDYN)
		{
			taup=sqrt(taup*E);
		}

		// Compute if inside damage surface
        g=((taup/rnp)*(taup/rnp))+((taun/rnn)*(taun/rnn))-1;
		if (g>tolerance)
		{
			//Q: trial
			//P: intersection ellipse with line OQ
			rhoQ = sqrt(taup*taup+taun*taun);
			if (taun > 1e-15) {
				thetaQ=atan(taup/taun);
			} else {
				thetaQ=PI/2.0;
			}
			// Modifica 26 Marzo 2012: Diego - Cambio formula per il calcolo di rhoP
			//rhoP=sqrt((rnp*rnp*rnn*rnn)/(rnn*rnn*sin(thetaQ)*sin(thetaQ)+rnp*rnp*cos(thetaQ)*cos(thetaQ)));
			rhoP = rnp*rnn*sqrt((taun*taun+taup*taup)/((taun*rnp)*(taun*rnp)+(taup*rnn)*(taup*rnn)));
			// Diego 23 Marzo 2012: piccola modifica per problemi numerici 
			// (può essere che rhoP venga un pelino più piccolo di rnp o un pelino più grande di rnn)
			if (rnn >= rnp)
			{
				if (rhoP<rnp)
					rhoP =rnp;
				if (rhoP>rnn)
					rhoP = rnn;
			}
			else if (rnn < rnp)
			{
				if (rhoP>rnp)
					rhoP =rnp;
				if (rhoP<rnn)
					rhoP = rnn;
			}
			// Fine modifica 23 Marzo 2012
			alfa=rhoQ/rhoP;
			// Valida solo per definizioni 2 e 4:
			if ((usedEquivalentTensionDefinition == ORIGINAL) || (usedEquivalentTensionDefinition == HOMOGENEOUS))
				alfa*=alfa;

			// Compute dea as usual
			for (i=0;i<3;i++)
			{
				for (j=0;j<3;j++)
				{
					st[i][j]=senp1t[i][j]*(1-1/alfa);
					dea[i][j]=st[i][j]*(1+ni)/E;
				};
			};
			diag=(st[0][0]+st[1][1]+st[2][2])*ni/(-E);
			for (i=0;i<3;i++)
				dea[i][i]=dea[i][i]+diag;
			
            norms=0;
	        for (i=0;i<3;i++)
			{
				for (j=0;j<3;j++)
				    norms+=senp1t[i][j]*senp1t[i][j];
			};
			norms=sqrt(norms);
	        for (i=0;i<3;i++)
			{
		        for (j=0;j<3;j++)
					ls[i][j]=senp1t[i][j]/norms;
			};
			
			pint=0;
			for (i=0;i<3;i++)
			{
				for (j=0;j<3;j++)
				{
					pint+=ls[i][j]*dea[i][j];
				};
			};

			if (pint>0)
			{
				lambdap=1-beta*E*pint/norms;

				for (i=0;i<3;i++)
				{
					for (j=0;j<3;j++)
						ses[i][j]=senp1t[i][j]*lambdap;

				};
				if ((usedEquivalentTensionDefinition == ORIGINAL)  || (usedEquivalentTensionDefinition == HOMOGENEOUS))
				{
					taun*=sqrt(lambdap);
					taup*=sqrt(lambdap);
				}
				else if (usedEquivalentTensionDefinition == COMPDYN)
				{
					taun*=lambdap;
					taup*=lambdap;
				}
				
				// check if inside of ellipse in taup,taun space
				g=(taup/rnp)*(taup/rnp)+(taun/rnn)*(taun/rnn)-1;
				if (g>tolerance)
				{
					for (i=0;i<3;i++)
					{
						for (j=0;j<3;j++)
							senp1t[i][j]=ses[i][j];
					};
					for (i=0;i<3;i++)
					{
						for (j=0;j<3;j++)
							trans[i][j]=ls[i][j]*(1+ni)/E;
					};
					diag=(ls[0][0]+ls[1][1]+ls[2][2])*ni/(-E);
					
					for (i=0;i<3;i++)
						trans[i][i]=trans[i][i]+diag;
					for (i=0;i<3;i++)
					{
						for (j=0;j<3;j++)
						{
							dep[i][j]=beta*E*pint*trans[i][j];
							epsilon_p_nplus1(i,j)=epsilon_p_n(i,j)+dep[i][j];
						};
					};
				}
				else
				{
					for (i=0;i<3;i++)
					{
						for (j=0;j<3;j++)
							epsilon_p_nplus1(i,j)=epsilon_p_n(i,j);
					};
				};
			}
			else {
				for (i=0;i<3;i++)
				{
					for (j=0;j<3;j++)
						epsilon_p_nplus1(i,j)=epsilon_p_n(i,j);
				};
			};
		}
		else
		{
			for (i=0;i<3;i++)
			{
				for (j=0;j<3;j++)
					epsilon_p_nplus1(i,j)=epsilon_p_n(i,j);
			};
		};
	};
	// Diego: aggiunto l'aggiornamento della matrice diagonale dg
	// e il calcolo di taun
	eigen_decomposition(senp1t,EE,V,dg);
	for (i=0;i<3;i++)
		dgn[i]=(dg[i]-fabs(dg[i]))/2.0;
	sigoct=(dgn[0]+dgn[1]+dgn[2])/3.0;
	tauoct=sqrt((dgn[0]-dgn[1])*(dgn[0]-dgn[1])+(dgn[0]-dgn[2])*(dgn[0]-dgn[2])+(dgn[1]-dgn[2])*(dgn[1]-dgn[2]))/3.0;
	if ((usedEquivalentTensionDefinition == ORIGINAL)  || (usedEquivalentTensionDefinition == HOMOGENEOUS))
	{
		taun=sqrt(3.0)*(k_p*sigoct+tauoct);
		if (taun >= 0)
			taun = sqrt(taun);
		else
			taun = 0;
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		taun=sqrt(3.0)*(k_p*sigoct+tauoct);
	}

	// Calcolo di taup
	for (i=0;i<3;i++)
		dgp[i]=(dg[i]+fabs(dg[i]))/2.0;
	diag=(dgp[0]+dgp[1]+dgp[2])*ni/(-E);
	for (i=0;i<3;i++)
	{
	    for (j=0;j<3;j++)
		{
			trans[i][j]=0.0;
			if (i==j)  // Diego!!!!
				trans[i][i]=dgp[i]*(1+ni)/E+diag;
		};
	};
	taup=0.0;
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
			if (i==j)   // Diego!!!!
				taup+=trans[i][j]*dgp[i];
	};
	if (usedEquivalentTensionDefinition == ORIGINAL)
	{
		taup=sqrt(sqrt(taup));
	}
	else if (usedEquivalentTensionDefinition == HOMOGENEOUS)
	{
		taup=sqrt(sqrt(taup*E));
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		taup=sqrt(taup*E);
	}

	// Aggiornamento di G
	g=(taup/rnp)*(taup/rnp)+(taun/rnn)*(taun/rnn)-1;
	if (g>tolerance)
	{
		// Evolution of damage
		//compute alfa
		//Q: trial
		//P: intersection ellipse with line OQ
		rhoQ = sqrt(taup*taup+taun*taun); 
		// Diego 23 Marzo 2012: prova modifica calcolo tangente.
		// Modifica 26 Marzo 2012: Diego - Cambio formula per il calcolo di rhoP
		if (taun > 1e-15) {
			thetaQ=atan(taup/taun);
		} else {
			thetaQ=PI/2.0;
		}
		// Fine modifica 23 Marzo 2012
		// Modifica 26 Marzo 2012: Diego - Cambio formula per il calcolo di rhoP
		//rhoP=sqrt((rnp*rnp*rnn*rnn)/(rnn*rnn*sin(thetaQ)*sin(thetaQ)+rnp*rnp*cos(thetaQ)*cos(thetaQ)));
		rhoP = rnp*rnn*sqrt((taun*taun+taup*taup)/((taun*rnp)*(taun*rnp)+(taup*rnn)*(taup*rnn)));
		// diego 23 Marzo 2012: piccola modifica per problemi numerici 
		// (può essere che rhoP venga un pelino più piccolo di rnp o un pelino più grande di rnn)
		if (rnn >= rnp)
		{
			if (rhoP<rnp)
				rhoP =rnp;
			if (rhoP>rnn)
				rhoP = rnn;
		}
		else if (rnn < rnp)
		{
			if (rhoP>rnp)
				rhoP =rnp;
			if (rhoP<rnn)
				rhoP = rnn;
		}
		// Fine modifica 23 Marzo 2012
		alfa=rhoQ/rhoP;

		//compute coordinates of point L
		// point where normal has 1 scope
		thetaL = atan((rnp*rnp)/(rnn*rnn));
		rhoL=sqrt((rnp*rnp*rnn*rnn)/(rnn*rnn*sin(thetaL)*sin(thetaL)+rnp*rnp*cos(thetaL)*cos(thetaL)));
       
		// Diego: 13 marzo 2012 - modificato IF per garantire soluzione in caso si inverta ellisse!
		//if (rhoP>rhoL) { // Diego: Versione vecchia ed eliminata il 13 marzo 2012 dell'IF
		if (((rhoP>rhoL) && (rhoP<=rnn)) || ((rhoP>=rnn) && (rhoP<rhoL))) {
			alfasp=1+(alfa-1)*(rnn-rhoP)/(rnn-rhoL);
			rnp1p=rnp*alfasp;
			rnp1n=sqrt((rnp1p*rnp1p*taun*taun)/(rnp1p*rnp1p-taup*taup));
		//} else {  // // Diego: Versione vecchia ed eliminata il 13 marzo 2012 dell'IF
		} else if (((rhoP>rhoL) && (rhoP<=rnp)) || ((rhoP>=rnp) && (rhoP<rhoL))) {
			alfasn=1+(alfa-1)*(rhoP-rnp)/(rhoL-rnp);
			rnp1n=rnn*alfasn;
			rnp1p=sqrt((rnp1n*rnp1n*taup*taup)/(rnp1n*rnp1n-taun*taun));
		} else {
			// Caso di Ellisse degenerata in Circonferenza
			alfasn=1+(alfa-1)*(rhoP-rnp)/(rhoL-rnp);
			rnp1n=rnn*alfasn;
			rnp1p=sqrt((rnp1n*rnp1n*taup*taup)/(rnp1n*rnp1n-taun*taun));
		}
	} 
	else
	{
		rnp1n=rnn;
		rnp1p=rnp;
    };
	
	// COMPUTE DAMAGE VARIABLES
	// compute negative damage
	if ((usedEquivalentTensionDefinition == ORIGINAL)  || (usedEquivalentTensionDefinition == HOMOGENEOUS))
	{
		dn=1-r0n/rnp1n*(1-An)-An*exp(Bn*(1-rnp1n/r0n));
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		dn=1-(sqrt(r0n))/(sqrt(rnp1n))*(1-An)-An*exp(Bn*(1-(sqrt(rnp1n)/(sqrt(r0n)))));
	}
		
	// Limiti su DN 
	// 11/01/2013 Diego Talledo: Added Negative Damage Limit
    if (dn>=dnMax)  
	{
		dn=dnMax;
		//dn = 0.99999; // togliere
	};
	// Diego: controllo che il danno negativo sia maggiore di zero. 13 marzo 2012 blocco dn
	if (dn<0.0)
	{
		dn=0.0;
	};

	// compute positive damage
	if ((usedEquivalentTensionDefinition == ORIGINAL)  || (usedEquivalentTensionDefinition == HOMOGENEOUS))
	{
		dp=1-((r0p*r0p)/(rnp1p*rnp1p))*exp(Ap*(1-(rnp1p*rnp1p)/(r0p*r0p)));
	}
	else if (usedEquivalentTensionDefinition == COMPDYN)
	{
		dp=1-((r0p)/(rnp1p))*exp(Ap*(1-(rnp1p)/(r0p)));
	}
		
	// Limiti su DP -> necessaria ulteriore ricerca
	// 11/01/2013 Diego Talledo: Added Positive Damage Limit
	if (dp>=dpMax)
	{
		dp=dpMax;
	};
    
	// 14/12/2011 Diego Talledo: Added Shear Retention Factor
	if (gammaC <= 0.0) {
		SRF12 = 0.0;
		SRF23 = 0.0;
		SRF13 = 0.0;
	}
	else {
		SRF12 = 1-abs(strain(0,1))/gammaC;
		SRF23 = 1-abs(strain(1,2))/gammaC;
		SRF13 = 1-abs(strain(0,2))/gammaC;
		if (SRF12 <= 0.0) {
			SRF12 = 0.0;
		}
		if (SRF23 <= 0.0) {
			SRF23 = 0.0;
		}
		if (SRF13 <= 0.0) {
			SRF13 = 0.0;
		}
	}
	
	// Reproject in local system
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
		{
			supp1[i][j]=0.0;
			supp2[i][j]=0.0;
			supp1[i][j]=supp1[i][j]+V[i][j]*dgn[j];
			supp2[i][j]=supp2[i][j]+V[i][j]*dgp[j];
		};
	};
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
		{
			Dn[i][j]=0.0;
			Dp[i][j]=0.0;
			for (k=0;k<3;k++)
			{
				Dn[i][j]=Dn[i][j]+supp1[i][k]*V[j][k];
				Dp[i][j]=Dp[i][j]+supp2[i][k]*V[j][k];
			};
		};
	};
	
	// 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
	// Compute dnstar and dpstar
	dnstar = dn + dchemn - (dchemn * dn);
	dpstar = dp + dchemp - (dchemp * dp);

	// 14/07/2014 Diego Talledo: Aggiunto dchemn variabile
	// dchem variabile: retta tra dchemn e 1. Con 1 raggiunto per epsilon pari a eps_u
	if (isDchemVariableSelected) {
		double dchem_var = (1-dchemn)/tau_n_eps_u*rnp1n+dchemn;
		dnstar = dn + dchem_var - (dchem_var * dn);
	}
	// Fine prova dchem variabile
	
	// Compute Cauchy Damaged Stress
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
		{
			//stress(i,j)=(1-dn)*Dn[i][j]+(1-dp)*Dp[i][j];
			stress(i,j)=(1-dnstar)*Dn[i][j]+(1-dpstar)*Dp[i][j]; // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage

			// 14/12/2011 Diego Talledo: Added Shear Retention Factor
			if ((i==1 && j==0) || (i==0 && j==1)) {
				//stress(i,j)=(1-(1-SRF12)*dnstar)*Dn[i][j]+(1-(1-SRF12)*dpstar)*Dp[i][j];
				if (!srfCompr)
					stress(i,j)=(1-dnstar)*Dn[i][j]+(1-(1-SRF12)*dpstar)*Dp[i][j];
				else
					stress(i,j)=(1-(1-SRF12)*dnstar)*Dn[i][j]+(1-(1-SRF12)*dpstar)*Dp[i][j];
			}
			if ((i==1 && j==2) || (i==2 && j==1)) {
				//stress(i,j)=(1-(1-SRF23)*dnstar)*Dn[i][j]+(1-(1-SRF23)*dpstar)*Dp[i][j];
				if (!srfCompr)
					stress(i,j)=(1-dnstar)*Dn[i][j]+(1-(1-SRF23)*dpstar)*Dp[i][j];
				else
					stress(i,j)=(1-(1-SRF23)*dnstar)*Dn[i][j]+(1-(1-SRF23)*dpstar)*Dp[i][j];
			}
			if ((i==2 && j==0) || (i==0 && j==2)) {
				//stress(i,j)=(1-(1-SRF13)*dnstar)*Dn[i][j]+(1-(1-SRF13)*dpstar)*Dp[i][j];
				if (!srfCompr)
					stress(i,j)=(1-dnstar)*Dn[i][j]+(1-(1-SRF13)*dpstar)*Dp[i][j];
				else
					stress(i,j)=(1-(1-SRF13)*dnstar)*Dn[i][j]+(1-(1-SRF13)*dpstar)*Dp[i][j];
			}
		};
	};
	
	//compute the tangent

	for ( ii = 0; ii < 6; ii++ ) {
		for ( jj = 0; jj < 6; jj++ )  {
			
			index_map( ii, i, j ) ;
			index_map( jj, k, l ) ;
			
			//elastic terms
			//tangent[i][j][k][l]  = (1-dn)*bulk_p * IbunI[i][j][k][l] ;
			//tangent[i][j][k][l] += (1-dn)*(2.0*shear_p) * IIdev[i][j][k][l] ;
			tangent[i][j][k][l]  = (1-dnstar)*bulk_p * IbunI[i][j][k][l] ;  // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage
			tangent[i][j][k][l] += (1-dnstar)*(2.0*shear_p) * IIdev[i][j][k][l] ;  // 20/06/2012 Diego Talledo: Added Environmental Chemical Damage

			//minor symmetries 
			tangent [j][i][k][l] = tangent[i][j][k][l] ;
			tangent [i][j][l][k] = tangent[i][j][k][l] ;
			tangent [j][i][l][k] = tangent[i][j][k][l] ;
		} // end for jj
	} // end for ii

	return ;
} 

// set up for initial elastic
void Concrete :: doInitialTangent( )
{
	int ii,jj,i,j,k,l;
	
	//compute the deviatoric strains
	for ( ii = 0; ii < 6; ii++ ) {
		for ( jj = 0; jj < 6; jj++ )  {

			index_map( ii, i, j ) ;
			index_map( jj, k, l ) ;

			//elastic terms
			initialTangent[i][j][k][l]  = bulk_p * IbunI[i][j][k][l] ;
			initialTangent[i][j][k][l] += (2.0*shear_p) * IIdev[i][j][k][l] ;

			//minor symmetries 
			initialTangent [j][i][k][l] = initialTangent[i][j][k][l] ;
			initialTangent [i][j][l][k] = initialTangent[i][j][k][l] ;
			initialTangent [j][i][l][k] = initialTangent[i][j][k][l] ;
		} // end for jj
	} // end for ii
	
	return ;
} 

//matrix_index ---> tensor indices i,j
void Concrete :: index_map( int matrix_index, int &i, int &j )
{
	switch ( matrix_index+1 ) { //add 1 for standard tensor indices
		
		case 1 :
			i = 1 ; 
			j = 1 ;
			break ;

		case 2:
			i = 2 ;
			j = 2 ; 
			break ;

		case 3 :
			i = 3 ;
			j = 3 ;
			break ;

		case 4 :
			i = 1 ;
			j = 2 ;
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

NDMaterial*
Concrete::getCopy (void)
{
	opserr << "Concrete::getCopy -- subclass responsibility\n"; 
	exit(-1);
	return 0;
}

const char*
Concrete::getType (void) const
{
	opserr << "Concrete::getType -- subclass responsibility\n";
	exit(-1);
	return 0;
}

int
Concrete::getOrder (void) const
{
	opserr << "Concrete::getOrder -- subclass responsibility\n";
	exit(-1);
	return 0;
}

int 
Concrete::commitState( ) 
{
	epsilon_p_n = epsilon_p_nplus1 ;
	rnp=rnp1p;
	rnn=rnp1n;

	strain_n = strain;
	stress_n = stress;

	return 0;
}


int 
Concrete::revertToLastCommit( ) 
{
	return 0;
}


int 
Concrete::revertToStart( )
{
	this->zero( ) ;
	
	return 0;
}

int
Concrete::sendSelf(int commitTag, Channel &theChannel)
{
	//opserr<<"sendSelf"<<endln;
	// we place all the data needed to define material and it's state
	// int a vector object
	static Vector data(14);
	int cnt = 0;
	data(cnt++) = this->getTag();
	data(cnt++) = E;
	data(cnt++) = ni;
	data(cnt++) = f01d;
	data(cnt++) = f02d;
	data(cnt++) = f0p;
	data(cnt++) = beta;
	data(cnt++) = An;
	data(cnt++) = Bn;
	data(cnt++) = Ap;
	data(cnt++) = rnn;
	data(cnt++) = rnp;
	data(cnt++) = rnp1n;
	data(cnt++) = rnp1p;
  
	// send the vector object to the channel
	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "Concrete::sendSelf - failed to send vector to channel\n";
		return -1;
	}

	return 0;
}

int
Concrete::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
	//opserr<<"recvSelf"<<endln;
	// recv the vector object from the channel which defines material param and state
	static Vector data(14);
	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "Concrete::recvSelf - failed to recv vector from channel\n";
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
	//xi_nplus1        = xi_n;

	return 0;
}
	
const Vector & Concrete::getFiberDamage()
{
	static Vector tmp(2);
	tmp(0) = this->dn;
	tmp(1) = this->dp;

	return tmp;
}

const Vector & Concrete::getDamage()
{
	static Vector tmp(2);
	tmp(0) = this->dn;
	tmp(1) = this->dp;

	return tmp;
}

// Diego Gennaio 2015 : Per implementazione Danno Globale
const Vector & Concrete::getElasticFreeEnergy()
{
	static Vector tmp(2); // Psi0- Psi0+

	// Valutazione sigma effetiva elastica e decomposizione spettrale
	static Matrix dev_strain(3,3) ; //deviatoric strain
	static Matrix dev_stress(3,3) ; //deviatoric stress

	double trace = 0.0 ; //trace of strain
	double senp1t[3][3];
	double V[3][3],EE[3][3];
	double dg[3],dgn[3],dgp[3];
	double diag;
	double trans[3][3];
	double supp1[3][3],supp2[3][3];

	// evaluate trace of elastic strain
	trace=0.0;
	for (int i=0;i<3;i++)
    {
		trace+=strain(i,i)-epsilon_p_n(i,i);  // è corretta, non trial perché stato già COMMITED
    };
	// Deviatoric strain tensor
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
			dev_strain(i,j)=strain(i,j)-epsilon_p_n(i,j);
    };
	for (int i=0;i<3;i++)
    {
		dev_strain(i,i)-= 1*trace/3;
    };

	// Deviatoric Stress
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
			dev_stress(i,j)=dev_strain(i,j)*2.0*shear_p;  // Strain(1,2) è epsilon12, per 2 = gamma12
    };

	// Sigma Eff (n+1,Trial) = DevStress + I * K * trace(StrainElastic)
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
			senp1t[i][j]=dev_stress(i,j);
    };
	for (int i=0;i<3;i++)
		senp1t[i][i]=senp1t[i][i]+bulk_p*trace;

	// Project in principal system
	eigen_decomposition(senp1t,EE,V,dg);

	// Negative and Positive part of principal Effective Stress
	for (int i=0;i<3;i++) {
		dgn[i]=(dg[i]-fabs(dg[i]))/2;
		dgp[i]=(dg[i]+fabs(dg[i]))/2;
	}
	diag=(dgp[0]+dgp[1]+dgp[2])*ni/(-E);
	for (int i=0;i<3;i++)
	{
	    for (int j=0;j<3;j++)
		{
			trans[i][j]=0.0;
			if (i==j)  // Diego!!!!
				trans[i][i]=dgp[i]*(1+ni)/E+diag;
		};
	};

	// Reproject in local system
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			supp1[i][j]=0.0;
			supp2[i][j]=0.0;
			supp1[i][j]=supp1[i][j]+V[i][j]*dgn[j];
			supp2[i][j]=supp2[i][j]+V[i][j]*dgp[j];
		};
	};
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			Dn[i][j]=0.0;
			Dp[i][j]=0.0;
			for (int k=0;k<3;k++)
			{
				Dn[i][j]=Dn[i][j]+supp1[i][k]*V[j][k];
				Dp[i][j]=Dp[i][j]+supp2[i][k]*V[j][k];
			};
		};
	};

	tmp(0) = tmp(1) = 0.0;
	for (int i = 0; i <3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			// Psi0-
			tmp(0)+=Dn[i][j]*(strain(i,j)-epsilon_p_n(i,j));
			// Psi0+
			tmp(1)+=Dp[i][j]*(strain(i,j)-epsilon_p_n(i,j));
		}
	}
	tmp(0)*=0.5;
	tmp(1)*=0.5;

	return tmp;
}

const Vector & Concrete::getDamagedFreeEnergy()
{
	static Vector tmp(2);
	const Vector &ElFreeEnergy = this->getElasticFreeEnergy();

	tmp(0) = ElFreeEnergy(0) * (1-this->dn);
	tmp(1) = ElFreeEnergy(1) * (1-this->dp);

	return tmp;
}
