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
 +                                                                          +
 |  Disclaimers:                                                            |
 |  (1) This is implementation of MultiaxialCyclicPlasticity for clays      |
 +      Model References:                                                   +
 |      Borja R.I, Amies, A.P. Multiaxial Cyclic Plasticity Model for       |
 |            Clays, ASCE J. Geotech. Eng. Vol 120, No 6, 1051-1070         |
 +      Montans F.J, Borja R.I. Implicit J2-bounding Surface Plasticity     +
 |            using Prager's translation rule. Int. J. Numer. Meth. Engng.  |
 |            55:1129-1166, 2002                                            |
 +      Code References:                                                    +
 |      Ignacio Romero and Adrian Rodriguez Marek, Brick element model with |
 |            a Multiaxial Cyclic Plasticity Model, in GEOFEAP, UC Berkeley |
 +  (2) Questions regarding this code should be directed to Gang Wang       +
 |  (3) Documentation could be found at                                     |
 |        www.ce.berkeley.edu/~wang/papers/MultiaxialCyclicPlasticity.pdf   |
 +                                                                          +
 |  Development History:                                                    |
 |  First Draft   -- April 2004                                             |
 +  Rewrite       -- Nov   2004                                             +
 |  Final Release --                                                        |
 |                                                                          |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                              User Command

   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    nDMaterial MultiaxialCyclicPlasticity $tag, $rho, $K, $G,
	                                      $Su , $Ho , $h, $m, $beta, $KCoeff
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   where:
      tag   : tag for this material
	  rho   : density
	  K     : buck modulus
	  G     : maximum (small strain) shear modulus
	  Su    : undrained shear strength, size of bounding surface R=sqrt(8/3)*Su
	  Ho    : linear kinematic hardening modulus of bounding surface
	  h     : hardening parameter
	  m     : hardening parameter
	  beta  : integration parameter, usually beta=0.5
	  KCoeff: coefficient of earth pressure, K0

  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#include <math.h>
#include <MultiaxialCyclicPlasticity.h>
#include <MultiaxialCyclicPlasticity3D.h> 
#include <MultiaxialCyclicPlasticityAxiSymm.h>
#include <MultiaxialCyclicPlasticityPlaneStrain.h>
#include <Information.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <elementAPI.h>

//parameters
const double MultiaxialCyclicPlasticity :: one3   = 1.0 / 3.0 ;
const double MultiaxialCyclicPlasticity :: two3   = 2.0 / 3.0 ;
const double MultiaxialCyclicPlasticity :: four3  = 4.0 / 3.0 ;
const double MultiaxialCyclicPlasticity :: root23 = sqrt( 2.0 / 3.0 ) ;
const double MultiaxialCyclicPlasticity :: infinity  = 1.0e12;

double MultiaxialCyclicPlasticity::initialTangent[3][3][3][3] ;   //material tangent
double MultiaxialCyclicPlasticity::IIdev[3][3][3][3] ; //rank 4 deviatoric 
double MultiaxialCyclicPlasticity::IbunI[3][3][3][3] ; //rank 4 I bun I 

int MultiaxialCyclicPlasticity::MaterialStageID = 2;   // classwide load stage tag, 
                                                       // elasto-plastic by default
int MultiaxialCyclicPlasticity::IncrFormulationFlag=1;

void * OPS_ADD_RUNTIME_VPV(OPS_MultiaxialCyclicPlasticity)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 10) {
	opserr << "WARNING: Insufficient arguments\n";
	opserr << "Want: nDMaterial MultiaxialCyclicPlasticity tag? rho? K? G? Su? Ho? h? m? beta? KCoeff? <eta?>" << endln;
	return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING invalid MultiaxialCyclicPlasticity tag\n";
	return 0;
    }

    double data[10];
    data[9] = 0.0;
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 10) {
	numdata = 10;
    }
    if (OPS_GetDoubleInput(&numdata,data)) {
	opserr << "WARNING invalid MultiaxialCyclicPlasticity double inputs\n";
	return 0;
    }

    NDMaterial* mat = new MultiaxialCyclicPlasticity(tag,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9]);
    if (mat == 0) {
	opserr << "WARNING: failed to create Multiaxialcyclicplasticity material\n";
	return 0;
    }

    return mat;
}

//zero internal variables
void MultiaxialCyclicPlasticity :: initialize ( ) 
{
  stress.Zero();
  strain.Zero();
  stress_n.Zero();
  strain_n.Zero();
  backs.Zero();                  // back stress fopr BS
  backs_n.Zero(); 
  so.Zero();                     // reversal deviatoric back stress
  so_n.Zero();
  

  // some flags

  flagjustunload=0;
  flagfirstload=0;
  icounter=0;
  iternum=0;
  plasticflag=0;
  plasticflag_n=0;
 
  // initialize s0 and kappa
  kappa = infinity;
  Psi   = 2 * shear;
  load  = 0.0; 
  //opserr<<"MCP::debug128: kappa = "<< kappa <<endln;
  X[1] = 2 * shear;
  X[2] = kappa;
  alp=0;
}


//null constructor
MultiaxialCyclicPlasticity ::  MultiaxialCyclicPlasticity( ) : 
NDMaterial( ), stress(3,3), strain(3,3),
stress_n(3,3), strain_n(3,3), backs(3,3), backs_n(3,3), so(3,3), so_n(3,3)
{ 
  bulk        = 0.0 ;
  shear       = 0.0 ;
  bulk_K0     = 0.0 ;
  shear_K0    = 0.0 ;
  eta         = 0.0 ;

  density     = 0.0 ; 

  this->initialize( ) ;     // or (*this).zero( ) 

  int i, j, k, l ;

  //zero rank4 IIdev and IbunI 
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ )  {
      for ( k = 0; k < 3; k++ ) {
		for ( l = 0; l < 3; l++)  { 

		  IbunI[i][j][k][l] = 0.0 ;

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
}


//full constructor
MultiaxialCyclicPlasticity::MultiaxialCyclicPlasticity(int    tag,
						       int    classTag,
						       double rho,
						       double K,
						       double G,
						       double Su,
						       double Ho_kin,
						       double Parameter_h,
						       double Parameter_m,
						       double Parameter_beta,
						       double Kcoeff,
						       double viscosity)
  : NDMaterial(tag, ND_TAG_MultiaxialCyclicPlasticity),
    stress(3,3), strain(3,3), stress_n(3,3), strain_n(3,3), 
    backs(3,3), backs_n(3,3), so(3,3), so_n(3,3)
{
  density     = rho; 
  bulk        = K ;
  shear       = G ;
  R           = sqrt(8.0/3.0)*Su;
  //opserr<<"MCP::235 R= "<<R <<endln; 
  Ho          = Ho_kin;        // kinematic hardening modulus
  h			  = Parameter_h; 
  m           = Parameter_m ;
  beta        = Parameter_beta ;
  eta         = viscosity ;
 
  ///////////// compute modulus for stage1 K0 consolidation ////////////
  K0 = Kcoeff;
  // Poisson's ratio v=K0/(1+K0)
  double v = K0 /(1+K0); 
  double E = 9 * bulk * shear / (3*bulk+shear);
  // compute shear and bulk modulus for K0 consolidation, cf. Fung, p.141
  shear_K0 =  E /(2*(1+v));
  bulk_K0  =  E/(3*(1-2*v));
  ///////////////////////////////////////////////////////////////
  
  if (tag==200) { // to be a fluid
	  v = 0.499;
      shear_K0 =  1;
      bulk_K0  *=  1000;
  }

  this->initialize( ) ;

  int i, j, k, l ;

  //zero rank4 IIdev and IbunI 
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      for ( k = 0; k < 3; k++ ) {
		for ( l = 0; l < 3; l++)  { 

		  IbunI[i][j][k][l] = 0.0 ;

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
}


//elastic constructor
MultiaxialCyclicPlasticity :: 
MultiaxialCyclicPlasticity(  int tag, int classTag, 
			     double rho, double K, double G ) 
:NDMaterial(tag, classTag),
stress(3,3), strain(3,3), stress_n(3,3), strain_n(3,3), 
backs(3,3), backs_n(3,3), so(3,3), so_n(3,3)
{ 
  density     = rho; 
  bulk        = K ;
  shear       = G ; 
  bulk_K0     = K ;
  shear_K0    = G ; 
  eta         = 0.0 ;

  this->initialize( ) ;

  int i, j, k, l ;

  //zero rank4 IIdev and IbunI 
  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ )  {
      for ( k = 0; k < 3; k++ ) {
		for ( l = 0; l < 3; l++)  { 

		  IbunI[i][j][k][l] = 0.0 ;

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
}


//destructor
MultiaxialCyclicPlasticity :: ~MultiaxialCyclicPlasticity( ) 
{  } 


NDMaterial*
MultiaxialCyclicPlasticity :: getCopy (const char *type)
{
    if (strcmp(type,"PlaneStress2D") == 0 || strcmp(type,"PlaneStress") == 0)
    {
    opserr << "MultiaxialCyclicPlasticity type plane stress material is NOT available now....";
	return 0;
    }
    else if (strcmp(type,"PlaneStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0)
    {
	MultiaxialCyclicPlasticityPlaneStrain  *clone ;
	clone = new MultiaxialCyclicPlasticityPlaneStrain(this->getTag(), density, bulk, shear, sqrt(3.0/8.0)*R,
		          Ho, h, m, beta, K0, eta) ; 
	return clone ;	
	}
    else if (strcmp(type,"AxiSymmetric2D") == 0 || strcmp(type,"AxiSymmetric") == 0)
    {
	MultiaxialCyclicPlasticityAxiSymm  *clone ;
	clone = new MultiaxialCyclicPlasticityAxiSymm(this->getTag(), density, bulk, shear, sqrt(3.0/8.0)*R,
		          Ho, h, m, beta, K0, eta) ; 
	return clone ;	
	}
    else if ((strcmp(type,"ThreeDimensional") == 0) || (strcmp(type,"3D") == 0))
    {
	MultiaxialCyclicPlasticity3D  *clone ;
	clone = new MultiaxialCyclicPlasticity3D(this->getTag(), density, bulk, shear, sqrt(3.0/8.0)*R,
						 Ho, h, m, beta, K0, eta) ; 
	return clone ;	
    }
    else if ( (strcmp(type,"PlateFiber") == 0) )
    {
    opserr << "MultiaxialCyclicPlasticity type plate fiber material is NOT available now....";
	return 0;
    }
    // Handle other cases
    else
    {
      opserr << "MultiaxialCyclicPlasticity::getModel failed to get model: " << type << endln;
      return 0;
    }
}

//print out material data
void MultiaxialCyclicPlasticity :: Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "MultiaxialCyclicPlasticity : " ; 
  s << this->getType( )  << endln ;
  s << "K    =   " << bulk        << endln ;
  s << "Gmax =   " << shear       << endln ;
  s << "Rho  =   " << density     << endln ;
  s << "R    =   " << R           << endln ;
  s << "Ho   =   " << Ho          << endln ;
  s << "h    =   " << h           << endln ;
  s << "m    =   " << m           << endln ;
  s << "beta =   " << beta        << endln ;
  s << "eta  =   " << eta         << endln ;
  s << endln ;
}


void MultiaxialCyclicPlasticity :: elastic_integrator( )
{
  static Matrix dev_strain(3,3) ; //deviatoric strain  

  static Matrix dev_stress(3,3) ; //deviatoric stress

  // add
  double pressure;                // 1/3 trace(stress) 

  double trace = 0.0 ; //trace of strain

  int i,j,k,l;
  int ii, jj ; 


  if (IncrFormulationFlag==0){

  trace = strain(0,0) + strain(1,1) + strain(2,2) ;
  
  //compute the deviatoric strains
  dev_strain = strain ;

  for ( i = 0; i < 3; i++ )
   dev_strain(i,i) -= ( one3*trace ) ;

  //compute the trial deviatoric stresses
  dev_stress = dev_strain;
  dev_stress *= 2.0 * shear_K0;	
  
  // compute the trial pressure    
  pressure   = trace ;
  pressure  *= bulk_K0 ;						 
  }

  static Matrix IncrStrain(3,3);

  static Matrix DevStress_n(3,3);

  static double pressure_n;

  if (IncrFormulationFlag==1){
   
	  
  // opserr<<"DP::elasticintegrator:stress_n"<<stress_n;
  // opserr<<"DP::elasticintegrator:strain_n"<<strain_n;

  IncrStrain = strain;
  IncrStrain -= strain_n ;

  trace = IncrStrain(0,0) + IncrStrain(1,1) + IncrStrain(2,2) ;

  dev_strain = IncrStrain ;
  for ( i = 0; i < 3; i++ )  dev_strain(i,i) -= ( one3*trace ) ;
 
  pressure_n  = one3 * (stress_n(0,0)+stress_n(1,1)+stress_n(2,2));
  
  DevStress_n = stress_n ;
  
  for ( i = 0; i < 3; i++ )  DevStress_n(i,i) -= pressure_n ;

  // Delta_S = 2*shear*Delta_e
  // Delta_p = bulk * Delta_theta
  
  // incremental formulation, NOTE: now dev_strain and trace are INCREMENTAL strains

  dev_stress = dev_strain;
  dev_stress *= 2.0 * shear_K0;
  dev_stress += DevStress_n;
  
  pressure = trace;
  pressure *= bulk_K0;
  pressure += pressure_n;

  }  
  
  // add on bulk part of stress, compute TOTAL stress at t=n+1

  stress = dev_stress ;
  for ( i = 0; i < 3; i++ )
     stress(i,i) += pressure ;


  //compute the tangent 
  
  for ( ii = 0; ii < 6; ii++ ) {
    for ( jj = 0; jj < 6; jj++ )  {

          index_map( ii, i, j ) ;
          index_map( jj, k, l ) ;

          //elastic terms
          tangent[i][j][k][l]  = bulk_K0 * IbunI[i][j][k][l] ;

          tangent[i][j][k][l] += (2.0*shear_K0) * IIdev[i][j][k][l] ;

          //minor symmetries 
          tangent [j][i][k][l] = tangent[i][j][k][l] ;
          tangent [i][j][l][k] = tangent[i][j][k][l] ;
          tangent [j][i][l][k] = tangent[i][j][k][l] ;

    } // end for jj
  } // end for ii

  flagfirstload=0; // prepare for MCP integrator
  
  return ;
} 




// set up for initial elastic
void MultiaxialCyclicPlasticity :: doInitialTangent( )
{
  int ii,jj,i,j,k,l;

  //compute the deviatoric strains
  for ( ii = 0; ii < 6; ii++ ) {
    for ( jj = 0; jj < 6; jj++ )  {

          index_map( ii, i, j ) ;
          index_map( jj, k, l ) ;

          //elastic terms
          initialTangent[i][j][k][l]  = bulk * IbunI[i][j][k][l] ;
          initialTangent[i][j][k][l] += (2.0*shear) * IIdev[i][j][k][l] ;

          //minor symmetries 
          initialTangent [j][i][k][l] = initialTangent[i][j][k][l] ;
          initialTangent [i][j][l][k] = initialTangent[i][j][k][l] ;
          initialTangent [j][i][l][k] = initialTangent[i][j][k][l] ;

    } // end for jj
  } // end for ii

  return ;
} 



//matrix_index ---> tensor indices i,j
void MultiaxialCyclicPlasticity :: index_map( int matrix_index, int &i, int &j )
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
MultiaxialCyclicPlasticity::getCopy (void)
{
  opserr << "MultiaxialCyclicPlasticity::getCopy -- subclass responsibility\n"; 
  exit(-1);
  return 0;
}

const char*
MultiaxialCyclicPlasticity::getType (void) const
{
    opserr << "MultiaxialCyclicPlasticity::getType -- subclass responsibility\n";
    exit(-1);
    return 0;
}

int
MultiaxialCyclicPlasticity::getOrder (void) const
{
    opserr << "MultiaxialCyclicPlasticity::getOrder -- subclass responsibility\n";
    exit(-1);
    return 0;
}


int 
MultiaxialCyclicPlasticity::commitState( ) 
{
 
  stress_n   = stress;     // april 2, 2004
  strain_n   = strain;     // april 2, 2004

  backs_n    = backs;
  so_n       = so; 
 
  Psi   = X[1];
  kappa = X[2];

  plasticflag_n=plasticflag;
  if (plasticflag==2 ){
	  plasticflag_n=1;
  }

  iternum=0;   // reset global newton iteration counter
  return 0;

}


int 
MultiaxialCyclicPlasticity::revertToLastCommit( ) 
{
  return 0;
}


int 
MultiaxialCyclicPlasticity::revertToStart( ) {

  // added: C.McGann, U.Washington for InitialStateAnalysis
  if (ops_InitialStateAnalysis) {
	// do nothing, keep state variables from last step
  } else {
	// normal call for revertToStart (not initialStateAnalysis)
    this->initialize( ) ;
  }

  return 0;
}

int
MultiaxialCyclicPlasticity::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(10); 
  int cnt = 0;
  data(cnt++) = this->getTag();
  data(cnt++) = density;   //add
  data(cnt++) = bulk;
  data(cnt++) = shear;
  data(cnt++) = R;   //add
  data(cnt++) = Ho;   //add
  data(cnt++) = h;
  data(cnt++) = m;
  data(cnt++) = beta;
  data(cnt++) = eta;
 
  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "MultiaxialCyclicPlasticity::sendSelf - failed to send vector to channel\n";
    return -1;
  }

  return 0;
}

int
MultiaxialCyclicPlasticity::recvSelf (int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{

  // recv the vector object from the channel which defines material param and state
  static Vector data(10);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "MultiaxialCyclicPlasticity::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag(data(cnt++));
  density = data(cnt++);
  bulk    = data(cnt++);
  shear   = data(cnt++);
  R       = data(cnt++);
  Ho      = data(cnt++);
  h       = data(cnt++);
  m       = data(cnt++);
  beta    = data(cnt++);
  eta     = data(cnt++);
  
 
  return 0;
}





double 
MultiaxialCyclicPlasticity::getRho()
{  
   return density; 
}



int MultiaxialCyclicPlasticity::updateParameter(int responseID, Information &info)
{

  // not finished
  this-> MaterialStageID = responseID;
   return 0;
}


// GEOFEAP subroutine model34(d,eps,sig,cc,cee,hn,h1,nhv,isw,n)
void MultiaxialCyclicPlasticity::plastic_integrator()
{ 

  int unloadflag;

  // toleraces

  static double ZERO      = 1.0e-16; 

  static double zerotol   = 1.0e-10;    // tol for zeroloading...||de||<zerotol
  static double tolforload= 1.0e-10;     // tol to determine loading/unloading
  static double tolforX   = 1.0e-6;;       // tol for compute norm g1g2 in iteration X[1], X[2]
  
  double twomu; 
  int ii,jj;              // for loop iterators
  int i,j,k,l;

  static Matrix de(3,3) ;   //incremental deviatoric strain  
  static Matrix  s(3,3) ;   //deviatoric stress
 // add
  static double p;          // pressure, 1/3 trace(stress) 
  static double e;          // incr. vol.strain, trace(incremental strains)

  static Matrix IncrStrain(3,3); // Frank let all Matrix be static
  static Matrix s_n(3,3);        // dev. stress at t_n 
  static double p_n;             // pressure at t_n
  //static Matrix soinit(3,3);   // save s0_n here
  double normchi=0;
  static Matrix strial(3,3);
  static Matrix chitri(3,3);
  static Matrix alpha_n(3,3); // backstress of loading surface at t_n
  static double  Psi_split; // Psi for strain split step
  static Matrix temp6(3,3);
  static double dottemp6;
  double norm = 0;
	
  double Hn =0;
  double Hn1=0;
  double ftrial;
  double temp;
  double normde;
  double fn;
  int zeroloadflag;
  int showdebugInfo;

  showdebugInfo=0;

  Vector debugInfo(10);

  debugInfo.Zero();
  for ( i = 0; i < 10; i++ )  debugInfo(i) = -1.0 ;

  unloadflag=0;
  zeroloadflag=0;
  load = 0;
  icounter=0;
  flagjustunload = 0;
  X[1]=0; X[2]=0;

  iternum++;

  // (0) calculate input incr. strain: (de, e)
  //     retrieve stress and stress of t_n step (s_n, p_n)


  IncrStrain = strain;
  IncrStrain -= strain_n ;       // incremental
  e  = IncrStrain(0,0) + IncrStrain(1,1) + IncrStrain(2,2) ;
  de = IncrStrain ;
  for ( i = 0; i < 3; i++ )  de(i,i) -= ( one3*e ) ;
 
  // compute p and s 
  p_n = one3 * (stress_n(0,0)+stress_n(1,1)+stress_n(2,2));
  s_n = stress_n ;

  for ( i = 0; i < 3; i++ )  s_n(i,i) -= p_n ;
 
  flagfirstload++;

  // initialize so_n if it is a very first call of this routine
  if (flagfirstload==1){
	  //opserr<<"MCP::firstload"<<endln;
	  so_n=s_n;
	  so=so_n;
	  backs_n=s_n;
	  kappa=infinity;
	  Psi=2*shear;
	  plasticflag=0;
	  X[1]=2*shear;
	  X[2]=infinity;
	  debugInfo(1)=1;
	  goto LABEL10;
  }

  so=so_n;


  // normde=||de||
  normde=0;
  for ( i = 0; i < 3; i++ ){
	  for ( j = 0; j < 3; j++ ) {
		  normde += de(i,j)*de(i,j) ;
	  }
  } //end for i 

  if (normde >= 0.0) {
	  normde=sqrt(normde); 
  } else {
	  opserr<<"MCP:1061:problem!!"<<endln;
  }

  if (normde<=zerotol)
  {
	  X[1] = 2*shear;
      X[2] = infinity;
	  plasticflag=0;
	  load = 1000; 
	  zeroloadflag=1;
	  debugInfo(1)=2;
	  //opserr<<"MCP::small strain incr" <<endln;
	  goto LABEL10;
  }
  MCPparameter(9) = normde ;



  //(1) Solve Kappa_n from s_n 
  //    assign plasticflag_n

   double kappasave;
   int plasticflagsave;

   kappasave=kappa;  // this kappa is obtained from previous solution of X[2]
   plasticflagsave=plasticflag_n;

   fn=0;

   chitri=s_n-backs_n;
   normchi=0;
   for ( i = 0; i < 3; i++ ){
	   for ( j = 0; j < 3; j++ ) {
		   normchi += chitri(i,j)*chitri(i,j) ;
	   }
	} //end for i 

    if (normchi >= 0) {
	  normchi =sqrt(normchi); 
	} else {
	  normchi = 0;
      opserr<<"Problem!!"<<endln;
	}

	fn = normchi - R;

	if ((fn >=0)&&(plasticflag_n!=1))
	{
		opserr<<"MCP::995:strange, fn="<<fn<<" plasticflag_n="<<plasticflag_n<<endln;
		showdebugInfo=1;
	}

	if (fn >=0 ){
        kappa=0;
		plasticflag_n=1;
		//opserr<<"MCP::tn on yield surface"<<endln;
	}
	else {

		plasticflag_n=0;
		// compute kappa from s_n, backs_n, so
		double n1 =0;
		double n2 =0;
		double t1dott2=0;

		n1=0.0;

		for ( i = 0; i < 3; i++ ){
			for ( j = 0; j < 3; j++ ) {
				temp=s_n(i,j)-backs_n(i,j);
			    n1 += temp*temp;
			}
		} //end for i 
		
		if (n1 >= 0.0) {
			n1=sqrt(n1); 
		}

		n2=0.0;
		for ( i = 0; i < 3; i++ ){
			for ( j = 0; j < 3; j++ ) {
                temp=s_n(i,j)-so_n(i,j);
				n2 += temp*temp  ;
			}
		} //end for i 

		if (n2 >= 0.0) {
			n2=sqrt(n2); 
		}

		if (n2==0){ // just unload
			kappa=infinity;
		} else{
			t1dott2=0;
			for ( i = 0; i < 3; i++ ){
				for ( j = 0; j < 3; j++ ) {
					t1dott2 += (s_n(i,j)-backs_n(i,j))*(s_n(i,j)-so_n(i,j)) ;
				}
			} //end for i 
		    t1dott2=t1dott2/(n1*n2);


		    // find kappa directly for t_n    cf Montans Eq. (15)
			//kappa =-n1/n2*t1dott2+sqrt(R*R/n2/n2-n1*n1/n2/n2*(1-t1dott2*t1dott2));
	              		
 
			// just another formula, cf Borja   Eq. (29), better than Montans!!
			// because Montans requires n1,n2!=0, borja only requires n2!=0
			// which is just unload case for n2=0

			temp=t1dott2*n1*n2,
			kappa =sqrt(temp*temp+n2*n2*(R*R-normchi*normchi))-t1dott2*n1*n2;
			kappa *=1.0/(n2*n2);

			if ((fabs(kappa-kappasave)>1e-6)&&((fabs(kappasave)<infinity))){
			   // opserr<<"MCP:992 big error in last computed X[2]" <<endln;
			//	opserr<<"MCP:992 kappa-kappaS   ="<<kappa-kappasave   <<endln;	
			//	opserr<<"MCP:992 kappa="<<kappa<<" kappaS="<<kappasave   <<endln;
			}
			
		}
	           
	}


   


	/// note: two options, (1) calculate tn (2) not, use previous value
    /// looks (2) is better!!!
	
	kappa=kappasave;
    //plasticflag_n=plasticflagsave;


	//(2) compute alpha_n

	for ( i = 0; i < 3; i++ ){
		for ( j = 0; j < 3; j++ ) {
			alpha_n(i,j)= so_n(i,j)+ (backs_n(i,j)-so_n(i,j)) /(1+kappa);
		}
	} //end for i 





	// (3) compute for Psi for t_n

	Hn=h*pow(kappa,m)+Ho; 

	Psi=2.0*shear*(1.0-shear/(shear+Hn/3.0));

	// in general Psi is not equal X[1]. since the latter is secant of
	// step n-> n+1, and Psi is secant of step n only

	
    Psi_split=2.0*shear*(1.0-shear/(shear+((1-beta)*Hn+beta*Ho)/3.0)); // Gang
    //Psi_split=2.0*shear/(1.0+3.0*shear*((1-beta)/Hn+beta/Ho)); // Borja



	// (4) check for loading/unloading

  load=0; temp6.Zero(); temp=0;
  for ( i = 0; i < 3; i++ ){
	  for ( j = 0; j < 3; j++ ) {
	    temp6(i,j) = s_n(i,j)-alpha_n(i,j);
        load += temp6(i,j)*de(i,j) ;
	    temp += temp6(i,j)*temp6(i,j);
	  }
  } //end for i 

  load *= 1.0/(sqrt(temp));


   if (load < tolforload*normde)  // unloading
	{   
	  if (plasticflag_n==1){
		  debugInfo(2)=1;
		  //opserr<<"MCP1091::unload checked!! from plastic, load="<<load<<" normde="<<normde<<endln;
		  //opserr<<"MCP1091::plasticflag_n="<<plasticflag_n<<endln;
		  //opserr<<"MCP1091::fn="<<fn<<" R="<<R<<" load="<<load<<endln;
	  } else {
		  debugInfo(2)=2;
		  //opserr<<"MCP1095::unload checked!! within B.S.,  load="<<load<<" normde="<<normde<<endln;
		  //opserr<<"MCP1095::plasticflag_n="<<plasticflag_n<<endln;
		  //opserr<<"MCP1095::fn="<<fn<<" R="<<R<<" load="<<load<<endln;
	  }
		
		kappa = infinity ;
		Psi   = 2*shear;
		so = s_n; 
		flagjustunload = 1;   // not used
		unloadflag = 1;  
		X[1] = 2*shear;
		X[2] = infinity;
		//unloadcounter +=1;
		plasticflag=0;     // force it go to plasticflag=0 tangent
		goto LABEL10;
		//goto LABEL6;
	}


  // (5) set plasticflag


 // begin: check for plastic loading if LAST loading is plastic
  if (plasticflag_n==1){    
  	  strial = s_n + 2.0 * shear * de; 
	  chitri = strial ;
	  chitri-= backs_n;

	  // check yield
	  normchi = 0 ;
	  for ( i = 0; i < 3; i++ ){
		  for ( j = 0; j < 3; j++ ) {
			normchi += chitri(i,j)*chitri(i,j) ;
		  }
		} //end for i 
      
      // compute normchi
	  if (normchi >= 0){
		  normchi = sqrt(normchi);
	  } else {
		  opserr<<"MCP 1060::Problem, normchi<0"<<endln;
	  }

	  ftrial=normchi-R;
	  if (ftrial > 0){
		  plasticflag=1;
		  debugInfo(3)=1;
	  } else {
		  debugInfo(3)=2;
		  opserr<<"MCP1419::unload from plastic" <<endln;
		  // update plasticflag
		  plasticflag = 0;
		  kappa = infinity ;
		  Psi   = 2*shear;
		  so    = s_n;   
		  X[1] = 2*shear;
		  X[2] = infinity;
		  flagjustunload = 1;   // not used
		  unloadflag = 1;
		  goto LABEL10;
		  //goto LABEL6;
	  }
  } else {  // last loading is NOT plastic, assign plasticflag=0 or 2
      debugInfo(3)=3;
	  chitri  = Psi*de;
	  chitri += s_n;
	  chitri -= backs_n;
 
	  normchi=0;
	  for ( i = 0; i < 3; i++ ){
		  for ( j = 0; j < 3; j++ ) {
		   normchi += chitri(i,j)*chitri(i,j) ;
		  }
		} //end for i 
      
      // compute normchi
	  if (normchi >= 0){
		  normchi =sqrt(normchi);
	  } else {
		  opserr<<"MCP 1089::Problem, normchi<0 normchi="<< normchi<<endln;
	  }   

	  plasticflag=0;

	  ftrial=normchi-R;
	  if (ftrial >= 0) {   
   			//Psi_split is calculated in the very beginning
			chitri  = Psi_split*de;
			chitri += s_n;
			chitri -= backs_n;
			
			normchi=0;
			for ( i = 0; i < 3; i++ ){
				for ( j = 0; j < 3; j++ ) {  
					normchi += chitri(i,j)*chitri(i,j) ;
				}
			} //end for i 

			// compute normchi
			if (normchi >= 0){
				normchi =sqrt(normchi);
			} else {
				opserr<<"MCP 1111::Problem, normchi<0"<<endln;
			}   

			ftrial=normchi-R;
			if (ftrial>=0){            
				debugInfo(3)=4;
				plasticflag=2;
			} else {
				debugInfo(3)=5;  // it is really an awkward situation when it comes here....
				plasticflag=0;   // diverge if it set to 1
			}
	  } //end if ftrial>0
  }
   
  // end: check for plastic loading




  // (6) if loading within B.S., solve  Psi_n+1, kappa_n+1 (i.e. X[1] and X[2])
 
LABEL6:
  if (plasticflag==0){
	// initialize

    // converge better
	X[1] = Psi  ;         // stiff at time n
	X[2] = kappa;         // kappa at time n


	//X[1] = 2*shear  ;        // stiff at time n
	//X[2] = infinity;         // kappa at time n
    
	double g1, g2;
	double aa, bb, cc, dd;	

	Hn  = h* pow(kappa, m) + Ho;
	Hn1 = h* pow(X[2],m) + Ho;
	// define g1, g2
	//g1= 2.0*shear-X[1]-3.0*shear*X[1]*((1-beta)/Hn + beta/Hn1);  // borja
	g1=X[1]-2.0*shear*(1.0-shear/(shear+((1-beta)*Hn+beta*Hn1)/3.0));    // Gang

	temp6.Zero();

	//temp6 = s-backs + X[1]*de + X[2]*(s+X[1]*de-so);  
    for ( i = 0; i < 3; i++ ){
		for ( j = 0; j < 3; j++ ) {
			temp6(i,j) = s_n(i,j)-backs_n(i,j) + X[1]*de(i,j) + X[2]*(s_n(i,j)+X[1]*de(i,j)-so(i,j));  
		}
	 } //end for i 


	dottemp6=0;
    for ( i = 0; i < 3; i++ ){
		for ( j = 0; j < 3; j++ ) {
	      dottemp6 += temp6(i,j)*temp6(i,j);
		}
	 } //end for i 

	 if (dottemp6 >=0){
		 g2=R-sqrt(dottemp6);
	 } else {
		 opserr<<"MCP1244::problem, dottemp6<0 "<<dottemp6<<endln;
		 g2=R;
	 }

	 norm = g1*g1+g2*g2;
 	 if (norm >= 0){
		 norm = sqrt(norm);
	 } else { 
		 opserr<<"MCP::problem, norm<0 "<<norm<<endln;
		 norm = 0.0; 
	 }

	
	 icounter = 0;
   // begin iterations to solve Psi_n+1, Kappa_n+1 !!!!	 
	 while ((norm>tolforX)&&(icounter<60)&&(X[2]>0))
	 //while ((icounter<30))
	 { 
	  //   compute Jacobian 
	  //   aa = dg1/dPsi    bb=dg1/dKappa
	  //   cc = dg2/dPsi    dd=dg2/dKappa
	  //  
	  //       | aa  bb |
	  //   J = |        |
	  //       | cc  dd |



		//Hn  = h* pow(kappa, m) + Ho;
		Hn1 = h* pow(X[2],m) + Ho;


		//aa = -1.0 - 3.0*shear*((1.0-beta)/Hn + beta/Hn1); // borja
		//bb = 3.0 * shear * X[1] * beta * m * h * pow(X[2],m-1.0) / pow(Hn1,2); // borja

		// reformulate 
		aa = 1.0;                                                 // Gang
		temp=shear + ((1.0-beta)*Hn+beta*Hn1)/3.0;

		bb= -2.0/3.0*beta*shear*shear/(temp*temp)*h*m*pow(X[2],m-1.0);  // Gang
        // we will have problem if X[2]=0, then bb=NaN. Problem: pow(X[2],m-1)


 		if (sqrt(dottemp6) > ZERO) {
			cc = 0 ;
			for ( i = 0; i < 3; i++ ){
				for ( j = 0; j < 3; j++ ) {
				   cc += temp6(i,j)*de(i,j);
				}
			} 
			cc *= -(1+X[2])/sqrt(dottemp6);

			dd = 0 ;
			for ( i = 0; i < 3; i++ ){
				for ( j = 0; j < 3; j++ ) {
		         dd += temp6(i,j)*(s_n(i,j)+X[1]*de(i,j)-so(i,j));
				}
			 }  
		    dd *= -1.0/sqrt(dottemp6);
		} else {
			opserr<<"MCP:: singularity in Jacobian, dottemp6="<<dottemp6<<endln;
		    //opserr<<"MCP:: icounter="<<icounter<<endln;
			opserr<<"MCP:: X[1]="<<X[1]<<" X[2]="<<X[2]<<endln;
			opserr<<"MCP:: plasticflag_n="<<plasticflag_n<<" kappa="<<kappa<<" Psi="<<Psi <<endln;
			cc = 0;
			dd = 0;
		}

		if (fabs(aa*dd-cc*bb)>=ZERO){
			X[1] += -1.0 /(aa*dd - cc*bb)*( dd*g1 - bb*g2);
			X[2] += -1.0 /(aa*dd - cc*bb)*(-cc*g1 + aa*g2);
		} else {
			opserr<<"MCP:: Fatal error: infinite Jacobian"  <<endln;
		    // remarks: arrive here maybe simply X[2] < 0, s.t. pow() gives NaN -1.#IND
			// so we must exit the loop once X[2]<0 is found
            opserr<<"MCP::pow()="<<pow(X[2],m)<<endln;
			opserr<<"MCP::X[2]="<<X[2]<<endln;  // Gang

			//exit(1);
		}
 
		if (X[2]<=0) {
 		   icounter = 100; // go out of the loop
		   //  exit the loop
		}

        // re-evaluate
		//Hn  = h* pow(kappa, m) + Ho;
		Hn1 = h* pow(X[2],m) + Ho;


 		// update g1, g2
        // g1= 2.0*shear-X[1]-3.0*shear*X[1]*((1.0-beta)/Hn + beta/Hn1);   // borja

		g1=X[1]-2.0*shear*(1.0-shear/(shear+((1.0-beta)*Hn+beta*Hn1)/3.0));    // Gang


 		temp6.Zero();
 
		for ( i = 0; i < 3; i++ ){
			for ( j = 0; j < 3; j++ ) {
				temp6(i,j) = s_n(i,j)-backs_n(i,j) + X[1]*de(i,j) + X[2]*(s_n(i,j)+X[1]*de(i,j)-so(i,j));  
			}
		}


		dottemp6 = 0;
		for ( i = 0; i < 3; i++ ){
			for ( j = 0; j < 3; j++ ) {
				dottemp6 += temp6(i,j)*temp6(i,j);
			}
		 } //end for i 

		if (dottemp6 >=0){
			 g2=R-sqrt(dottemp6);
		} else {
			 // remark: come here wrong because X[2]=0
             opserr<<"MCP1353::aa="<<aa<<" bb="<<bb<<" cc="<<cc<<" dd="<<dd<<endln;
			 opserr<<"MCP1353::X[2]="<<X[2]<<endln;
			 opserr<<"MCP1353::pow(X[2],m-1)"<<pow(X[2],m-1)<<endln;
			 icounter=1353;
		}
 	 	norm = g1*g1+g2*g2;
		if (norm > 0){
			norm = sqrt(norm);
		} else { 
			norm = 0.0; 
		}
 
		icounter +=1;
    
	 } // end of while loop
 
	 // begin: check kappa converged to a positive value
    debugInfo(4)=1;

    if ((norm > tolforX)&&(X[2]!=0)){
	 opserr<<endln<<endln<<"MCP::X[1] X[2] is not converged!! norm = "<<norm <<" icounter="<<icounter<<endln;
     opserr<<"MCP::plasticflag_n= "<<plasticflag_n <<" plasticflag="<< plasticflag<<endln;
	 opserr<<"MCP::X[2] ="<< X[2]<<endln<<endln;
	 opserr<<"MCP::debugInfo= "<<debugInfo<<endln;
	 showdebugInfo=1;
	}

	if (X[2]<=0.0)
	 {   debugInfo(4)=2;
		 
		 //opserr<< endln;
		 //opserr<<"MCP1299::WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endln;
		 //opserr<<"MCP1299::norm = "<<norm <<" icounter="<<icounter<<endln;
		 //opserr<<"MCP1299::kappa= "<<kappa<<" Psi ="<<Psi<<endln;
		 //opserr<<"MCP1299::Negative value for X[2] (kappa_n+1) "<< X[2]<<" set it to 0 "<<endln<<endln;
			 
 
	     X[2]=0.0;
		 
		 if (unloadflag == 1){
			 debugInfo(4)=3;
			 plasticflag=0;
			 X[1] = 2*shear;
		 } else {
			chitri = Psi_split*de;
			chitri += s_n;
			chitri -= backs_n;
			normchi=0;
			for ( i = 0; i < 3; i++ ){
				for ( j = 0; j < 3; j++ ) {
					normchi += chitri(i,j)*chitri(i,j) ;
				}
			} //end for i 

			// compute normchi
			if (normchi >= 0){
				normchi =sqrt(normchi);
			} else {
				opserr<<"MCP 1111::Problem, normchi<0"<<endln;
			}   

			ftrial=normchi-R;
			if (ftrial>=0){
				debugInfo(4)=4;
				plasticflag=2;
				X[1]=Psi_split; //Psi_split=2.0*shear*(1-shear/(shear+((1-beta)*Hn+beta*Ho)/3.0)); // Gang
			} else {
				debugInfo(4)=X[2];
				plasticflag=0;  // diverge if it set to 1
				// question: should plasticflag=0?? or 1?? for this case?
				X[1]=2.0*shear*(1.0-shear/(shear+((1.0-beta)*Hn+beta*Ho)/3.0));
				//X[1]=Psi;  // cannot converge if use this
			}
		 }
			 
	} //end if X[2]<0
	 
 
 } //end if (plasticflag==0)



 // (7) if strain split, calculate ratio alp and secant modulus Psi_split

 if (plasticflag==2){
		   // we are here when the loading is partially inside BS, partially on BS
		   // alp is defined as epstotal=eps(inside)*alp+eps(outside)*(1-alp)
		   // value of Psi between kappa (t=n) inside BS and kappa=0 on BS 	
	        debugInfo(4)=6;
	
 			double chichi;
			double chide ;
			chichi=0; chide=0;
			for ( i = 0; i < 3; i++ ){
				for ( j = 0; j < 3; j++ ) {
					chichi += (s_n(i,j)-backs_n(i,j))*(s_n(i,j)-backs_n(i,j));
					chide  += (s_n(i,j)-backs_n(i,j))*de(i,j);
				} // end for j
			} //end for i 

			// compute strain splitter
			alp = (-chide+sqrt(chide*chide+normde*normde*(R*R-chichi)));
			alp = alp /(Psi_split *normde*normde); 

 
		    double insidesqrt;
			insidesqrt=chide*chide+normde*normde*(R*R-chichi);

			if ((alp > 1.0)||(insidesqrt<0)||(alp < 0.0) ){
				debugInfo(4)=7;
				opserr<<"MCP:1394::WRONG!!! alpha="<<alp<<" chichi-R*R="<<chichi-R*R<<endln;
				opserr<<"MCP::debugInfo= "<<debugInfo<<endln;
				// remark: will have problem if chichi>R*R
				alp=0;
				plasticflag=1;
				showdebugInfo=1;
			}

 }

	
 //(8) Compute consistent tangent, consistency parameter, and update stresses

////////////////////////////// TANGENT /////////////////////

LABEL10:

   debugInfo(0)=plasticflag_n;
   debugInfo(6)=normde;
   debugInfo(7)=plasticflag;
   debugInfo(8)=X[1];
   debugInfo(9)=X[2];

 if (plasticflag==0)  // bounding surface mapping rule
 {  

	      debugInfo(5)=1;
	      // update cauchy stresses using incremental strain
		  twomu=X[1];

		  s  = s_n ;
		  s += twomu * de;
		  p = p_n + bulk  * e ;
          // add on bulk part of stress, compute TOTAL stress at t=n+1
		  stress = s; 
		  for ( i = 0; i < 3; i++ )  stress(i,i)= s(i,i) + p ;
		  
		  backs=backs_n;
		  // check stress_n+1 is on the yield surface
		  
		  temp6=s-backs;
		  temp=0;
		  for ( i = 0; i < 3; i++ ){
			  for ( j = 0; j < 3; j++ ) {
				  temp += temp6(i,j)*temp6(i,j);
			  } // end for j
		  } //end for i 

		 
		  if (sqrt(temp)-R>=0){
			  if (unloadflag!=1){
				  plasticflag=1;
			  }

			// Remark:
			// a big problem for this update is, upon unload from plastic, the 
			// predicted de is still very high, may cause this >=0. And if you
			// just simply let plasticflag = 1, it is not unload at all. Continue
			// to yield at plastic
		    // opserr<<"MCP1535::plasticflag=0: |S|-R = " << sqrt(temp)-R<<endln;
		  }


		  // a consistent tangent

          // compute tangent again
	      for ( ii = 0; ii < 6; ii++ ) {
			  for ( jj = 0; jj < 6; jj++ )  {
				  index_map( ii, i, j ) ;
				  index_map( jj, k, l ) ;
				  tangent[i][j][k][l]  = bulk * IbunI[i][j][k][l] ;
				  tangent[i][j][k][l] += twomu * IIdev[i][j][k][l] ;
				  //minor symmetries 
				  tangent [j][i][k][l] = tangent[i][j][k][l] ;
				  tangent [i][j][l][k] = tangent[i][j][k][l] ;
				  tangent [j][i][l][k] = tangent[i][j][k][l] ;
			  }
		  }
   } 

   if (plasticflag > 0) // plasticflag = 1,2   plastic loading   pp.15
   { 
	    
	    if ((X[1]==2.0*shear)&&(unloadflag!=1)&&(zeroloadflag!=1)){
			  opserr<<"MCP::warning...WHY X[1]=2G at plasticflag>0 and not unload?"<<endln;
			  opserr<<"MCP::debugInfo= "<<debugInfo<<endln;
			  showdebugInfo=1;
		}


		twomu   = 2.0 * shear;
		if (plasticflag==1)
		{	debugInfo(5)=2;
			alp = 0.0;
			strial  = s_n ;
			strial += twomu*de;
			chitri  = strial;
			chitri -= backs_n;
			Psi_split = 0; // just give a number
		}//end plasticfl=1
 
		if (plasticflag==2)
		{
			debugInfo(5)=3;
			strial  = s_n ;
			strial += alp * Psi_split * de;
		    strial += (1-alp)*twomu*de;
			chitri  = strial;
			chitri -= backs_n;
		} // end if plasticflag=2

		normchi = 0.0;
		// double check yield
		for ( i = 0; i < 3; i++ ){
			for ( j = 0; j < 3; j++ ) {
				normchi += chitri(i,j)*chitri(i,j);
			} // end for j
		} //end for i 
		if (normchi >= 0) {
			normchi = sqrt (normchi);
		} 
		else {
			opserr<<"MCP1873:: normchi = " << normchi <<endln;
		    opserr<<"MCP:: plastic loading with zero stresses, problem" <<endln; 
		}

		ftrial =  normchi - R; 

		if (ftrial < ZERO)
		{   // Remark: come here, means, unload, but goes to wrong place
			opserr<<"MCP1472:: ftrial = " <<ftrial<<endln;
			opserr<<"MCP1472:: strange! set ftrial=0, plasticflag="<<plasticflag <<endln;
			opserr<<"MCP1372:: plasticflag_n="<<plasticflag_n<<" fn="<<fn<<" normde="<<normde<<" Psi_split="<<Psi_split<<endln;
			opserr<<"MCP1372:: load="<<load<<endln;
			showdebugInfo=1;
		}

		// compute consistency parameter
		double gamma;
		gamma =  ftrial /(twomu + 2.0*Ho/3.0);
		// leave stress slightly outside yield surface
        gamma *= (1.0 - 1e-08) ; // add July 2, 2004 !!!!!! very important add !!!!!
				
		// update plastic variable for backstresses
		
		backs  = backs_n;
		backs +=  two3 * Ho * gamma * chitri/ normchi ;

        s  = s_n + alp * Psi_split * de;
		s += twomu * ((1-alp) * de - gamma * chitri/normchi );  


		// check stress_n+1 is on the yield surface
  		temp6=s-backs;
		temp=0;
		for ( i = 0; i < 3; i++ ){
			for ( j = 0; j < 3; j++ ) {
				temp += temp6(i,j)*temp6(i,j);
				} // end for j
		} //end for i 

		 
		if ((sqrt(temp)-R>(1.0e-3)*R)||(sqrt(temp)-R<0)){
           opserr<<"MCP1690::plastic: alp   = " << alp <<endln;
		   opserr<<"MCP1690::plastic: |S|-R = " << sqrt(temp)-R<<endln;
		   opserr<<"MCP::debugInfo= "<<debugInfo<<endln;
		   showdebugInfo=1;
		   // Remark: dangerous if <0... 
		}

		//s  = s_n;
		//s += twomu * ( de - gamma * chitri/normchi) ;  
		p = p_n + bulk  * e;
		
		stress = s ;
		for ( i = 0; i < 3; i++ )
			stress(i,i) += p ;   // add on pressure part

 
		// remark: this update is very important, since when commit
		// we assign Psi=X[1], kappa=X[2]
		if (Ho==0){
			X[1] = 0; 
		} else 
		{	
			X[1] = 2.0 * shear /(1.0 + 3.0*shear/Ho);
		}
        X[2] = 0 ;   // kappa = 0


        //compute the tangent
		double theta1 = 1.0 - 2.0 *shear *gamma/normchi;
		double theta2 = 1.0/(1.0+Ho/3.0/shear) - (1.0-theta1);

		double NbunN;

		for ( ii = 0; ii < 6; ii++ ) {
			for ( jj = 0; jj < 6; jj++ )  {

				index_map( ii, i, j ) ;
				index_map( jj, k, l ) ;

				NbunN  = chitri(i,j)*chitri(k,l)/(normchi*normchi) ; 

				// elasto-plastic term  (1-alp) * C_ep
    			// ****** Elasto-plastic consistent modulus, from Simo & Hughes Box 3.2
				// ---------------------------------------------------------------------
                //        C_ep = K IbunI + 2 shear ( theta1* IIdev - theta2 * NbunN)
				// _____________________________________________________________________

                ///opserr<<"MCP tangent! changed July 18"<<endln;
				tangent[i][j][k][l]  = (1-alp)* bulk * IbunI[i][j][k][l] ;
				tangent[i][j][k][l] += (1-alp)* 2.0*shear*theta1 * IIdev[i][j][k][l] ;
				tangent[i][j][k][l] -= (1-alp)* 2.0*shear*theta2 * NbunN ;
				
				// bounding surface terms alp * C_bound
				//****** Consistent modulus from Bounding Surface Plasticity
				//----------------------------------------------------------------------- 
				//       C_bound = K IbunI + Psi_split IIdev
				//_______________________________________________________________________

				tangent[i][j][k][l] += alp * bulk * IbunI[i][j][k][l] ;
                tangent[i][j][k][l] += alp * Psi_split * IIdev[i][j][k][l] ;
              
				// tangent = alp * C_bound + (1-alp) * C_ep;

				//minor symmetries 
				tangent[j][i][k][l] = tangent[i][j][k][l] ;
				tangent[i][j][l][k] = tangent[i][j][k][l] ;
				tangent[j][i][l][k] = tangent[i][j][k][l] ;

			} // end for jj
		} // end for ii
   } // end if(plasticflag)



   if (showdebugInfo==1){
	  opserr<<"END OF INTEGRATOR::debugInfo= "<<debugInfo<<endln;   
   }
	  //if ((iternum>20)&&(this->getEleTag()==151)){
	   //opserr<<"MCP::debugInfo= "<<debugInfo<<endln;
       //opserr<<"MCP::plasticflag= "<<plasticflag<<endln;
	  //}
}  // end of this subroutine


Vector MultiaxialCyclicPlasticity::MCPparameter(10) ;




Vector& 
MultiaxialCyclicPlasticity::getMCPparameter()
{
 //opserr<<"getMCPparameter"<<endln;
 MCPparameter(0) = plasticflag;
 MCPparameter(1) = X[1];
 MCPparameter(2) = X[2]; // kappa
 MCPparameter(3) = alp ;
 MCPparameter(4) = stress(0,1) ;
 MCPparameter(5) = backs(0,1)  ;


 double norm = 0 ;
 double p = one3 * (stress(0,0)+stress(1,1)+stress(2,2));
 Matrix s = stress ;
 int i; 
 for (i = 0; i < 3; i++ )  s(i,i) -= p ;
 
 for (i = 0; i < 3; i++ ){
   for (int j = 0; j < 3; j++ ) {
     norm   += (s(i,j)-backs(i,j)) * (s(i,j)-backs(i,j)) ;
   } // end for j
 } //end for i 

 MCPparameter(6) = sqrt(norm) ;
 MCPparameter(7) = load ;

 norm=0;
 for (i = 0; i < 3; i++ ){
	for (int j = 0; j < 3; j++ ) {
	   norm   += (strain(i,j)) * (strain(i,j)) ;
	} // end for j
  } //end for i 

 MCPparameter(8) = norm ;
 // MCPparameter(9) = normde ;   get directly from subroutine plastic_integrator
 return MCPparameter;

}
