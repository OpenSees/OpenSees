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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/CycLiqCP.cpp,v $

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

#include <CycLiqCP.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include <CycLiqCP3D.h>
#include <CycLiqCPPlaneStrain.h>

//parameters
const double CycLiqCP :: one3   = 1.0 / 3.0 ;
const double CycLiqCP :: two3   = 2.0 / 3.0 ;
const double CycLiqCP :: four3  = 4.0 / 3.0 ;
const double CycLiqCP :: root23 = sqrt( 2.0 / 3.0 ) ;
const double CycLiqCP :: root32 = sqrt( 3.0 / 2.0 ) ;
const double CycLiqCP :: pat=101.0;
const double CycLiqCP :: pcut=0.5;    // cut off confining stress
const double CycLiqCP :: pmin=0.5;    // cut off confining stress

double CycLiqCP::initialTangent[3][3][3][3] ;   //material tangent
double CycLiqCP::IIdev[3][3][3][3] ; //rank 4 deviatoric 
double CycLiqCP::IbunI[3][3][3][3] ; //rank 4 I bun I 
double CycLiqCP::mElastFlag = 0;

#include <elementAPI.h>
#define OPS_Export 

static int numCycLiqCPMaterials = 0;

void *
OPS_CycLiqCPMaterial(void)
{
  if (numCycLiqCPMaterials == 0) {
    numCycLiqCPMaterials=1;
    //OPS_Error("\nCycLiqCP - Written: Rui Wang, Jian-Min Zhang, Gang Wang, Tsinghua University\nPlease refer to: Zhang and Wang, 2012, 'Large post-liquefaction deformation of sand, part I: physical mechanism, constitutive description and numerical algorithm', Acta Geotechnica\n", 1);
    opserr<<"\nCycLiqCP - Written: Rui Wang, Jian-Min Zhang, Gang Wang, Tsinghua University\nPlease refer to: Zhang and Wang, 2012, 'Large post-liquefaction deformation of sand, part I: physical mechanism, constitutive description and numerical algorithm', Acta Geotechnica\n";
  }

  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 12) {
    opserr << "Want: nDmaterial CycLiqCP tag? G0? kappa? h? Mfc? dre1? Mdc? dre2? rdr? eta? dir? ein? <rho?>" << endln;
    return 0;	
  }
  
  int tag;
  double dData[12];

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial CycLiqCP material  tag" << endln;
    return 0;
  }
  if (numArgs == 12) {
  numData = 11;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial CycLiqCP  with tag: " << tag << endln;
    return 0;
  }
    theMaterial = new CycLiqCP(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], 
                                            dData[6], dData[7], dData[8], dData[9], dData[10],  0.0);
  }
  else if (numArgs > 12) {
  numData = 12;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial CycLiqCP  with tag: " << tag << endln;
    return 0;
  }
    theMaterial = new CycLiqCP(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], 
                                            dData[6], dData[7], dData[8], dData[9], dData[10],  dData[11]);
  }
  
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial CycLiqCP  with tag: " << tag << endln;
  }

  return theMaterial;
}

//static vectors and matrices
Vector CycLiqCP :: strain_vec(6) ;
Vector CycLiqCP :: stress_vec(6) ;
Matrix CycLiqCP :: tangent_matrix(6,6) ;
Matrix CycLiqCP :: I(6,6);


//zero internal variables
void CycLiqCP :: zero ( ) 
{
  strain_n.Zero();
  alpha_n.Zero();
  alpha_nplus1.Zero();
  r.Zero();

  epsvir_n=0.;
  epsvir_nplus1=0.;
  epsvre_n=0.;
  epsvre_nplus1=0.;
  epsvc_n=0.;
  epsvc_nplus1=0.;
  gammamono=0.;
  gammamonos=0.;
  etam=0.;
  lambda=0.0;

  epsvc_ns0=0.0;
  epsvc_ns01=0.0;

  stress_n.Zero();

  stress_nplus1.Zero();
  strain_nplus1.Zero();
  initializeState = true;
  p0=pmin;
}


//null constructor
CycLiqCP ::  CycLiqCP( ) : 
NDMaterial( ),
strain_n(3,3),
alpha_n(3,3),
r(3,3),
stress_n(3,3),
strain_nplus1(3,3),
alpha_nplus1(3,3),
stress_nplus1(3,3)
{ 
  G0=0.0;
  kappa=0.0;
  h=0.0;
  Mfc=0.0;    
  dre1=0.0;
  Mdc=0.0;
  dre2=0.0;
  rdr=0.0;
  eta=0.0;
  dir=0.0;    
  ein=0.0;
  rho=0.0;    

  this->zero( ) ;     // or (*this).zero( ) 

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

  I.Zero();
  I(0,0)=I(1,1)=I(2,2)=1.0;

  plastic_integrator();
}


//full constructor
CycLiqCP :: CycLiqCP(
                       int    tag,
					   int classTag,
                       double G01,
                       double kappa1,
                       double h1,
                       double Mfc1,  
                       double dre11,
                       double Mdc1,
                       double dre21,
                       double rdr1,
                       double eta1,
                       double dir1,
                       double ein1,  
                       double rho1

	) 
: NDMaterial(tag,classTag),
strain_n(3,3),
alpha_n(3,3),
r(3,3),
stress_n(3,3),
strain_nplus1(3,3),
alpha_nplus1(3,3),
stress_nplus1(3,3)
{
    G0=G01;
    kappa=kappa1;
    h=h1;
    Mfc=Mfc1;    
    dre1=dre11;
    Mdc=Mdc1;
    dre2=dre21;
    rdr=rdr1;
    eta=eta1;
    dir=dir1;   
	rho=rho1;    
	ein=ein1; 

  this->zero( ) ;

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

  I.Zero();
  I(0,0)=I(1,1)=I(2,2)=1.0;

    plastic_integrator();
}

//destructor
CycLiqCP :: ~CycLiqCP( ) 
{  } 

//send back type of material
const char* CycLiqCP :: getType( ) const 
{
    opserr << "CycLiqCP::getCopy -- subclass responsibility\n"; 
  exit(-1);
  return 0;
}

int
CycLiqCP::getOrder (void) const
{
     opserr << "CycLiqCP::getCopy -- subclass responsibility\n"; 
  exit(-1);
  return 0;
}

NDMaterial *
CycLiqCP::getCopy (const char *type)
{
if ((strcmp(type,"ThreeDimensional") == 0) ||
	     (strcmp(type,"3D") == 0))
    {
	CycLiqCP3D  *clone ;
	clone = new CycLiqCP3D(this->getTag(), G0, kappa, h, Mfc, dre1, Mdc, dre2, rdr, eta, dir, ein, rho) ;
	return clone ;	
  }
    else if (strcmp(type,"PlaneStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0)
    {
	CycLiqCPPlaneStrain  *clone ;
	clone = new CycLiqCPPlaneStrain(this->getTag(), G0, kappa, h, Mfc, dre1, Mdc, dre2, rdr, eta, dir, ein, rho) ;
	return clone ;
    }

    // Handle other cases
    else
    {
      return NDMaterial::getCopy(type);
    }
}

NDMaterial*
CycLiqCP :: getCopy ()
{
	opserr << "CycLiqCP::getCopy -- subclass responsibility\n"; 
  exit(-1);
  return 0;
}



//--------------------Plasticity-------------------------------------

//plasticity integration routine
void CycLiqCP :: plastic_integrator( )
{
  const double tolerance = 1.0e-8*pcut ; //tolerance
  const double tole = 1.0e-6*pcut ; //tolerance
  const double dt = ops_Dt ; //time step

  static Matrix dev_strain(3,3) ; //deviatoric strain

  static Matrix dev_strain_n(3,3) ; //deviatoric strain of last step

  static Matrix ddev_strain_p(3,3); //plastic deviatoric strain increment

  static Matrix dev_stress(3,3) ; //deviatoric stress

  static Matrix dev_stress_n(3,3) ; //deviatoric stress of last step
 
  static Matrix normal(3,3) ;     //(dev_stress-alpha_n)/p_n/sqrt(2/3)/m

  Matrix pass(3,3); //matrix to pass on values
  double G0ep;

  double phi = 0.0 ; //trial value of yield function
  double phi_n=0.0;

  double trace = 0.0 ; //trace of strain
  double trace_n=0.0; //trace of strain from previous step
  double dtracep=0.0; //plastic volumetric strain increment
  double LDR=0.0;
  double iconv=0.0;
  double lambdamax=0.0;
  double lambdamin=0.0;
  double sdila=1.0;
  double fy_pmax=0.0;
  double pmin2;
  double epsvc02;
  double epsvre_ns,epsvir_ns,epsvc_ns,epsvc_ns00=0.0;
  Matrix DR(3,3),LD(3,3),alpha_ns;
  DR.Zero();
  LD.Zero();

  double intm; //intermediate variable
  double ax,bx,cx,beta;
  int isub,sub1,sub2,sub;
  int i,j,k,l,wr,wr1;
  int ii, jj ;

  //compression as positive
	  stress_n=-1.0*stress_n;
	  strain_nplus1=-1.0*strain_nplus1;
	  strain_n=-1.0*strain_n;

  //compute the deviatoric stress of last step
  p_n=one3*(stress_n(0,0) + stress_n(1,1) + stress_n(2,2) );
  
  dev_stress_n=stress_n;
  for ( i = 0; i < 3; i++ )
      dev_stress_n(i,i) -= ( p_n ) ;
  if ((p_n)<pmin)
  {   
	  p_n=pmin;
      //dev_stress_n.Zero();
  }
      

  //compute the deviatoric strains

  trace = strain_nplus1(0,0) + strain_nplus1(1,1) + strain_nplus1(2,2) ;
  trace_n=strain_n(0,0) + strain_n(1,1) + strain_n(2,2) ;

  dev_strain_n = strain_n ;
  for ( i = 0; i < 3; i++ )
    dev_strain_n(i,i) -= ( one3*trace_n ) ;

  dev_strain = strain_nplus1 ;
  for ( i = 0; i < 3; i++ )
    dev_strain(i,i) -= ( one3*trace ) ;


  en=(1+ein)*exp(-trace_n)-1.;//void ratio

  // force elastic response if initialization analysis is designated===================================
	if (mElastFlag == 0) {
		p_nplus1=pat*pow(sqrt(p_n/pat)+(1+en)/2.0/kappa*(trace-trace_n),2);
	    K=(1+en)/kappa*pat*sqrt((p_n+p_nplus1)/2.0/pat);
		G=G0*pat*(pow((2.97-en),2)/(1+en))*sqrt((p_n+p_nplus1)/2.0/pat);
		dev_stress=dev_stress_n+2.0 * G*(dev_strain-dev_strain_n);
		K=(1+en)/kappa*pat*sqrt(p_nplus1/pat);
        G=G0*pat*(pow((2.97-en),2)/(1+en))*sqrt(p_nplus1/pat);
  		for ( ii = 0; ii < 6; ii++ ) {
        for ( jj = 0; jj < 6; jj++ )  {
      
              index_map( ii, i, j ) ;
              index_map( jj, k, l ) ;
      
              //elastic terms
              tangent[i][j][k][l]  = K * IbunI[i][j][k][l] ;
      
              tangent[i][j][k][l] += (2.0*G) * IIdev[i][j][k][l] ;
      
              //minor symmetries 
              tangent [j][i][k][l] = tangent[i][j][k][l] ;
              tangent [i][j][l][k] = tangent[i][j][k][l] ;
              tangent [j][i][l][k] = tangent[i][j][k][l] ;
      
          } // end for jj
        } // end for ii
		initializeState = false;
 // proceed with full algorithm if initialization analysis is not designated==========================
    }else if (mElastFlag == 1) {
        // set initial confining pressure p0
    	if (!initializeState) {
    		p0=p_n;
			epsvc0=-2*kappa/(1+en)*(sqrt(p0/pat)-sqrt(pmin/pat));
			//opserr<<"p0="<<p0<<"\n";
			r=dev_stress_n/(p_n);
			epsvc_ns0=0.0;
			epsvc_ns01=0.0;
	        etam=root32*sqrt(doublecontraction(r,r));
    		initializeState = true;
    	}
  // --------------(I)Initialize-------------------------------------
		alpha_nplus1=alpha_n;
		epsvir_nplus1=epsvir_n;
		epsvre_nplus1=epsvre_n;
		epsvc_nplus1=epsvc_n+trace-trace_n;
		epsvc_ns01=epsvc_ns0+trace-trace_n;
		r=dev_stress_n/(p_n);
		ddev_strain_p.Zero();
        dtracep=0.0;
        lambda=0.0;
		p0=p_n;
		epsvc0=-2*kappa/(1+en)*(sqrt(p0/pat)-sqrt(pmin/pat));

  // --------------Trial Before Substep-------------------------------------
        if (epsvc_nplus1>epsvc0)
		    p_nplus1=pat*pow(sqrt(p0/pat)+(1+en)/2.0/kappa*epsvc_nplus1,2);
		else
		{
			p_nplus1=pmin;
		}
	    K=(1+en)/kappa*pat*sqrt((p_n+p_nplus1)/2.0/pat);
		G=G0*pat*(pow((2.97-en),2)/(1+en))*sqrt((p_n+p_nplus1)/2.0/pat);
		dev_stress=dev_stress_n+2.0 * G*(dev_strain-dev_strain_n);
		r_nplus1=dev_stress/p_nplus1;
		sub1=(int)(root32*sqrt(doublecontraction(r_nplus1-r,r_nplus1-r))/0.05)+1;
		sub2=(int)(sqrt(two3*doublecontraction(dev_strain-dev_strain_n,dev_strain-dev_strain_n))/0.001)+1;
		sub=sub1;
		if (sub2>sub1)
			sub=sub2;
		if (sub>100)
			sub=100;
		alpha_ns=alpha_n;
        epsvir_ns=epsvir_n;
        epsvre_ns=epsvre_n;
        epsvc_ns=epsvc_n;
		epsvc_ns00=epsvc_ns0;
		gammamonos=gammamono;
 // --------------Trial Substep End-------------------------------------

  // --------------(II)Elastic Predict-------------------------------------
	for (isub=0;isub<sub;isub++)
	{
		// --------------(I)Initialize-------------------------------------
		alpha_nplus1=alpha_ns;
		epsvir_nplus1=epsvir_ns;
		epsvre_nplus1=epsvre_ns;
		epsvc_nplus1=epsvc_ns+(trace-trace_n)/sub;
		epsvc_ns01=epsvc_ns0+(trace-trace_n)/sub;
		r=dev_stress_n/(p_n);
		ddev_strain_p.Zero();
        dtracep=0.0;
        lambda=0.0;
		p0=p_n;
		epsvc0=-2*kappa/(1+en)*(sqrt(p0/pat)-sqrt(pmin/pat));
		if (epsvc_nplus1<epsvc0)
		{
			p_nplus1=pmin;
		}
		else
		{
			p_nplus1=pat*pow(sqrt(p0/pat)+(1+en)/2.0/kappa*epsvc_nplus1,2);
			
		}
	    K=(1+en)/kappa*pat*sqrt((p_n+p_nplus1)/2.0/pat);
		G=G0*pat*(pow((2.97-en),2)/(1+en))*sqrt((p_n+p_nplus1)/2.0/pat);
		dev_stress=dev_stress_n+2.0 * G*(dev_strain-dev_strain_n)/sub;
		r_nplus1=dev_stress/p_nplus1;
		eta_n=root32*sqrt(doublecontraction(r,r));
		if (etam>Mfc-tolerance)
		{
			etam=Mfc-tolerance;
		}
		if (eta_n>Mfc-tolerance)
		{
			eta_n=Mfc-tolerance;
		}
		if (eta_n>etam)
		{
			etam=eta_n;
		}
		if (etam-eta_n<=tolerance)
		{
			beta=1.0;
		}
		else
		{
	    	ax=doublecontraction(r-alpha_ns,r-alpha_ns);
	        bx=2*doublecontraction(alpha_ns,r-alpha_ns);
	        cx=doublecontraction(alpha_ns,alpha_ns)-two3*etam*etam;
	        beta=(-bx+sqrt(bx*bx-4*ax*cx))/2.0/ax;
	        if (beta<1.0-tolerance)
	      	  {
	    		  opserr << "beta="<<beta<<"\n";
	    		  beta=(-bx-sqrt(bx*bx-4*ax*cx))/2.0/ax;
	    	}
	        if (beta<1.0-tolerance)
	        {opserr << "beta="<<beta<<"\n";
	        exit(-1);
	        }
		}
		rbar=alpha_ns+beta*(r-alpha_ns);
		normal=3.0/2.0*rbar/root32/sqrt(doublecontraction(rbar,rbar));
		N=doublecontraction(r,normal);
  // --------------(III)Loading/Unloading-------------------------------------
	  
	  phi=doublecontraction(dev_stress-dev_stress_n,normal)-(p_nplus1-p_n)*N;
	  phi_n=doublecontraction(r_nplus1-r,normal);
	    // --------------(IV)Unloading-------------------------------------
	  if (phi<tolerance||phi_n<tolerance)
	  {
		  gammamonos=0.0;
		  alpha_nplus1=r;
	  }
  // --------------(V)Loading-------------------------------------
	  else if (phi>tolerance&&phi_n>tolerance)
	  {
		  epsvc_nplus1=epsvc_ns+(trace-trace_n)/sub-dtracep;
		  loadindex=0.0;
		  lambda=0.0;
		  lambdamin=0.0;
		  lambdamax=0.0;
		  rou=root32*sqrt(doublecontraction(r-alpha_ns,r-alpha_ns));
		  roubar=root32*sqrt(doublecontraction(rbar-alpha_ns,rbar-alpha_ns));
		  if (roubar>tolerance)
		  {
		  H=h*G*(Mfc/etam*roubar/rou-1.0);
		  if (H<tolerance)
		  {
			  H=tolerance;
		  }
		  eta_nplus1=root32*sqrt(doublecontraction(r_nplus1,r_nplus1));
		  sdila=1.0;
		  Dre_n=dre1*(Mdc-sdila*eta_n);
		  if (doublecontraction(r,r_nplus1-r)<0.0)
		  {
			  sdila=-1.0;
		      Dre_n=dre1*(Mdc-sdila*eta_n);
		  }
		  if (Dre_n>0.0)
		  {
			  Dre_n=pow(-dre2*epsvre_ns,2);
			  if (-epsvre_n<tolerance)
				  Dre_n=0.0;
		  }
		  Dir_n=dir*exp(-eta*epsvir_ns)/pow((1.0+gammamonos/rdr),2);
		  D=Dir_n+Dre_n;
		  prem=p_n;
		  wr=1;
		  iconv=0.0;
		  do
		  {
			  if (iconv==0.0)
			  {
				  dlambda=phi/(H+3*G-K*D*N);
				  lambda+=dlambda;
			  }
			  else
			  {
				  lambda=0.5*(lambdamax+lambdamin);
			  }
			  prem=p_nplus1;
			  loadindex=H*lambda;
			  dtracep=lambda*D;
	          ddev_strain_p=lambda*normal;
			  epsvc_nplus1=epsvc_ns+(trace-trace_n)/sub-dtracep;
			  epsvc_ns01=epsvc_ns00+(trace-trace_n)/sub-dtracep;
			  if (epsvc_nplus1<epsvc0)
			  {
				  p_nplus1=pmin;
				  epsvc_nplus1=epsvc0;
			  }
		      else
			  {
				  p_nplus1=pat*pow(sqrt(p0/pat)+(1+en)/2.0/kappa*epsvc_nplus1,2);

			  }
			  G=G0*pat*(pow((2.97-en),2)/(1+en))*sqrt((p_n+p_nplus1)/2.0/pat);
			  dev_stress=dev_stress_n+2*G*((dev_strain-dev_strain_n)/sub-ddev_strain_p);

			  phi=doublecontraction(dev_stress-dev_stress_n,normal)-(p_nplus1-p_n)*N-loadindex;
			  if (phi<-tolerance)
			  {
				  iconv=1.0;
				  lambdamax=lambda;
			  }
			  if (phi>tolerance && iconv==1.0)
				  lambdamin=lambda;
			  if (wr>20000)
			  {
				  opserr<<"lambdamax="<<lambdamax<<"\t";
				  opserr<<"p_nplus1="<<p_nplus1<<"\t";
				  opserr<<"iconv="<<iconv<<"\t";
				  opserr<<"phimax="<<phi<<"\n";
				  break;
			  }
			  wr=wr+1;
				epsvir_nplus1=lambda*Dir_n+epsvir_ns;
		        epsvre_nplus1=lambda*Dre_n+epsvre_ns;
				
		  }while (abs(phi)>tolerance);
		  gammamonos=gammamonos+lambda;
		  }

	  }

	  r_nplus1=dev_stress/p_nplus1;
	  eta_nplus1=root32*sqrt(doublecontraction(r_nplus1,r_nplus1));
	  if ( eta_nplus1>Mfc)
	  {
		  intm=eta_nplus1/Mfc+tolerance;
		  dev_stress=dev_stress/intm;
	  }
	  eta_nplus1=root32*sqrt(doublecontraction(dev_stress,dev_stress))/p_nplus1;
	  r_nplus1=dev_stress/p_nplus1;
	  if ( eta_nplus1>Mfc)
	  {
		  intm=eta_nplus1/Mfc+tolerance;
		  dev_stress=dev_stress/intm;
		  opserr<<"you shouldn't see this="<<intm<<"\n";
	  }
	    alpha_ns=alpha_nplus1;
        epsvir_ns=epsvir_nplus1;
        epsvre_ns=epsvre_nplus1;
		
        epsvc_ns=epsvc_ns-epsvc_nplus1+(trace-trace_n)/sub-dtracep;
		epsvc_ns00=epsvc_ns01;
		epsvc_nplus1=epsvc_ns;
        p_n=p_nplus1;
        dev_stress_n=dev_stress;
	  }
	  K=(1+en)/kappa*pat*sqrt((p_nplus1+p_nplus1)/2.0/pat);
      G=G0*pat*(pow((2.97-en),2)/(1+en))*sqrt((p_nplus1+p_nplus1)/2.0/pat);

  		for ( ii = 0; ii < 6; ii++ ) {
        for ( jj = 0; jj < 6; jj++ )  {
      
              index_map( ii, i, j ) ;
              index_map( jj, k, l ) ;
      
              //elastic terms
              tangent[i][j][k][l]  = K * IbunI[i][j][k][l] ;
      
              tangent[i][j][k][l] += (2.0*G) * IIdev[i][j][k][l] ;
      
              //minor symmetries 
              tangent [j][i][k][l] = tangent[i][j][k][l] ;
              tangent [i][j][l][k] = tangent[i][j][k][l] ;
              tangent [j][i][l][k] = tangent[i][j][k][l] ;
      
          } // end for jj
        } // end for ii

  if (lambda!=0.0)
  {
	  
	  r=dev_stress/(p_nplus1);
      eta_n=root32*sqrt(doublecontraction(r,r));
      if (eta_n>etam)
      {
      	etam=eta_n;
      }
	  rou=root32*sqrt(doublecontraction(r-alpha_nplus1,r-alpha_nplus1));
		  if (etam-eta_n<=tolerance)
		{
			beta=1.0;
		}
		else
		{
	        ax=doublecontraction(r-alpha_nplus1,r-alpha_nplus1);
	        bx=2*doublecontraction(alpha_nplus1,r-alpha_nplus1);
	        cx=doublecontraction(alpha_nplus1,alpha_nplus1)-two3*etam*etam;
	        beta=(-bx+sqrt(bx*bx-4*ax*cx))/2.0/ax;
	        if (beta<1.0-tolerance)
	      	{
	  	  	    opserr << "beta="<<beta<<"\n";
	  	  	    beta=(-bx-sqrt(bx*bx-4*ax*cx))/2.0/ax;
	  	    }
	        if (beta<1.0-tolerance)
	        {
				opserr << "beta="<<beta<<"\n";
	            exit(-1);
	        }
		 }
	  rbar=alpha_nplus1+beta*(r-alpha_nplus1);
	  normal=3.0/2.0*rbar/root32/sqrt(doublecontraction(rbar,rbar));
	  N=doublecontraction(r,normal);

	  L=normal-one3*N*I;

      roubar=root32*sqrt(doublecontraction(rbar-alpha_nplus1,rbar-alpha_nplus1));
      H=h*G*(Mfc/etam*roubar/rou-1.0);

	  if (H<0.01*G)
	  {
		  H=0.01*G;
	  }
	  Dre_n=dre1*(Mdc-sdila*eta_n);
	  if (Dre_n>0.0)
	  {
		  Dre_n=pow(-dre2*epsvre_n,2);
		  if (-epsvre_n<tolerance)
			  Dre_n=0.0;
	  }
      Dir_n=dir*exp(-eta*epsvir_nplus1)/pow((1.0+gammamonos/rdr),2);
      D=Dir_n+Dre_n;
	  if (p_nplus1<(pmin+tolerance))
	  {
		  D=0;
	  }

	  R=normal+one3*D*I;

	  DR=doublecontraction42(tangent,R);
	  LD=doublecontraction24(L,tangent);
	  LDR=4.0*(doublecontraction(LD,R)+doublecontraction(LD,R));

	  for (i =0; i<3; i++) {
      for (j =0; j<3; j++) {
      for ( k = 0; k < 3; k++ ) {
      for ( l = 0; l < 3; l++ )  {
		  tangent[i][j][k][l]+=-(DR(i,j)*LD(k,l))/(H+LDR);
	  }
	  }
	  }
	  }
  }

}


    stress_nplus1=dev_stress;
	for ( i = 0; i < 3; i++ )
    stress_nplus1(i,i) += p_nplus1;

      //tension as positive
	  stress_n=-1.0*stress_n;
	  strain_nplus1=-1.0*strain_nplus1;
	  strain_n=-1.0*strain_n;
	  stress_nplus1=-1.0*stress_nplus1;

  return ;
} 




// set up for initial elastic
void CycLiqCP :: doInitialTangent( )
{
  int ii,jj,i,j,k,l;

  //compute the deviatoric strains
  for ( ii = 0; ii < 6; ii++ ) {
    for ( jj = 0; jj < 6; jj++ )  {

          index_map( ii, i, j ) ;
          index_map( jj, k, l ) ;

          //elastic terms
          initialTangent[i][j][k][l]  = K * IbunI[i][j][k][l] ;
          initialTangent[i][j][k][l] += (2.0*G) * IIdev[i][j][k][l] ;

          //minor symmetries 
          //minor symmetries 
          initialTangent [j][i][k][l] = initialTangent[i][j][k][l] ;
          initialTangent [i][j][l][k] = initialTangent[i][j][k][l] ;
          initialTangent [j][i][l][k] = initialTangent[i][j][k][l] ;

    } // end for jj
  } // end for ii

  return ;
} 

//double contraction of two 2 order tensors
double CycLiqCP :: doublecontraction(Matrix a, Matrix b)
{
	int i,j;
	double contra_result=0.0;
	for (i=0;i<a.noRows();i++) {
		for (j=0;j<a.noCols();j++) {
			contra_result+=a(i,j)*b(i,j);
		}
	}
	return contra_result;
}

//double contraction of 4 and 2 order tensors
Matrix CycLiqCP :: doublecontraction42(double b[][3][3][3], Matrix a)
{
	int i,j,k,l;
	Matrix Da(3,3) ;
	Da.Zero();
      for (i =0; i<3; i++) {
      for (j =0; j<3; j++) {
      for ( k = 0; k < 3; k++ ) {
      for ( l = 0; l < 3; l++ )  {
		  Da(i,j)=Da(i,j)+b[i][j][k][l]*a(k,l);
	  }
	  }
	  }
	  }
	  return Da;
}

//double contraction of 4 and 2 order tensors
Matrix CycLiqCP :: doublecontraction24(Matrix a, double b[][3][3][3])
{ 
	int i,j,k,l;
	Matrix aD(3,3) ;
	aD.Zero();
      for (i =0; i<3; i++) {
      for (j =0; j<3; j++) {
      for ( k = 0; k < 3; k++ ) {
      for ( l = 0; l < 3; l++ )  {
		  aD(i,j)=aD(i,j)+b[k][l][i][j]*a(k,l);
	  }
	  }
	  }
	  }
	  return aD;
}

//matrix_index ---> tensor indices i,j
void CycLiqCP :: index_map( int matrix_index, int &i, int &j )
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

//print out material data
void CycLiqCP :: Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "CycLiqCP : " ; 
  s << endln ;
}


int 
CycLiqCP::commitState( ) 
{
  strain_n        = strain_nplus1 ;
  alpha_n=alpha_nplus1;
  epsvir_n=epsvir_nplus1;
  epsvre_n=epsvre_nplus1;
  gammamono=gammamonos;
  epsvc_n=epsvc_nplus1;
  etam=etam;
  stress_n=stress_nplus1;
  epsvc_ns0=epsvc_ns01;
  //epsvc0=epsvc0;

  return 0;
}

int 
CycLiqCP::revertToLastCommit( ) 
{
  return 0;
}


int 
CycLiqCP::revertToStart( ) {

  // added: C.McGann, U.Washington for InitialStateAnalysis
  if (ops_InitialStateAnalysis) {
	// do nothing, keep state variables from last step
  } else {
	// normal call for revertToStart (not initialStateAnalysis)
    this->zero( ) ;
  }

  return 0;
}

int
CycLiqCP::sendSelf(int commitTag, Channel &theChannel)
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
	  cnt=cnt+1;
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
CycLiqCP::recvSelf (int commitTag, Channel &theChannel, 
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
	  cnt=cnt+1;
	}
  }


  return 0;
}


int
CycLiqCP::setParameter(const char **argv, int argc, Parameter &param)
{
  	if (argc < 2)
    	return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0],"updateMaterialStage") == 0) {
			return param.addObject(1, this);
		}
	}

    return -1;
}

int
CycLiqCP::updateParameter(int responseID, Information &info)
{
	// called updateMaterialStage in tcl file
	if (responseID == 1) {
		mElastFlag = info.theInt;
	}
	return 0;
}
