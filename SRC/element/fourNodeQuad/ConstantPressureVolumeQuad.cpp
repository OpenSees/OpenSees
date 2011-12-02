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
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.17 $
// $Date: 2006-08-04 19:07:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/ConstantPressureVolumeQuad.cpp,v $

// Ed "C++" Love
//
// Constant Presssure/Volume Four Node Quadrilateral
// Plane Strain (NOT PLANE STRESS)


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <ConstantPressureVolumeQuad.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static data
double ConstantPressureVolumeQuad::matrixData[64];
Matrix ConstantPressureVolumeQuad :: stiff(matrixData,8,8)   ;
Vector ConstantPressureVolumeQuad :: resid(8)     ;
Matrix ConstantPressureVolumeQuad :: mass(8,8)    ;
Matrix ConstantPressureVolumeQuad :: damping(8,8) ;
 
//volume-pressure constants
double ConstantPressureVolumeQuad :: one3  = 1.0 / 3.0 ;
double ConstantPressureVolumeQuad :: two3  = 2.0 / 3.0 ;
double ConstantPressureVolumeQuad :: four3 = 4.0 / 3.0 ;
double ConstantPressureVolumeQuad :: one9  = 1.0 / 9.0 ;
    
//quadrature data
double ConstantPressureVolumeQuad :: root3 = sqrt(3.0) ;
double ConstantPressureVolumeQuad :: one_over_root3 = 1.0 / root3 ;

double ConstantPressureVolumeQuad :: sg[] = { -one_over_root3,  
					       one_over_root3, 
					       one_over_root3, 
					      -one_over_root3 } ;

double ConstantPressureVolumeQuad :: tg[] = { -one_over_root3, 
                                              -one_over_root3, 
                                               one_over_root3,  
                                               one_over_root3 } ;

double ConstantPressureVolumeQuad :: wg[] = { 1.0, 1.0, 1.0, 1.0 } ;
  

//null constructor
ConstantPressureVolumeQuad :: ConstantPressureVolumeQuad( ) :
Element( 0, ELE_TAG_ConstantPressureVolumeQuad ),
connectedExternalNodes(4), load(0)
{ 
  for (int i=0; i<9; i++)
    materialPointers[i] = 0;
}


//full constructor
ConstantPressureVolumeQuad :: ConstantPressureVolumeQuad( 
                            int tag, 
                  	    int node1,
			    int node2,
			    int node3,
			    int node4,
			    NDMaterial &theMaterial ) :
Element( tag, ELE_TAG_ConstantPressureVolumeQuad ),
connectedExternalNodes(4), load(0)
{
  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;

  int i ;
  for ( i = 0 ;  i < 4; i++ ) {

      materialPointers[i] = theMaterial.getCopy("AxiSymmetric2D") ;

      if (materialPointers[i] == 0) {
	opserr << "ConstantPressureVolumeQuad::constructor - failed to get a material of type: AxiSymmetric2D\n";
	exit(-1);
      } //end if
      
  } //end for i 

}


//destructor 
ConstantPressureVolumeQuad :: ~ConstantPressureVolumeQuad( )
{
  int i ;
  for ( i = 0 ;  i < 4; i++ ) {

    delete materialPointers[i] ;
    materialPointers[i] = 0 ; 

    nodePointers[i] = 0 ;
  } //end for i

  if (load != 0)
    delete load;
}


//set domain
void ConstantPressureVolumeQuad :: setDomain( Domain *theDomain ) 
{  
  int i ;
  for ( i = 0;  i < 4; i++ ) {

     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i)  ) ;

     if ( nodePointers[i] != 0 ) {
	 const Vector &coor = nodePointers[i]->getCrds( ) ;

	 xl[0][i] = coor(0) ;
	 xl[1][i] = coor(1) ; 
     } // end if 

  } //end for i 
  
  this->DomainComponent::setDomain(theDomain);
}


//get the number of external nodes
int ConstantPressureVolumeQuad :: getNumExternalNodes( ) const
{
  return 4 ;
} 


//return connected external nodes
const ID& ConstantPressureVolumeQuad :: getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 


Node **
ConstantPressureVolumeQuad::getNodePtrs(void) 
{
  return nodePointers;
} 


//return number of dofs
int ConstantPressureVolumeQuad :: getNumDOF( ) 
{
  return 8 ;
}


//commit state
int ConstantPressureVolumeQuad :: commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "ConstantPressureVolumeQuad::commitState () - failed in base class";
  }    

  for (int i = 0; i < 4; i++ ) 
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}
 


//revert to last commit 
int ConstantPressureVolumeQuad :: revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int ConstantPressureVolumeQuad :: revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

int 
ConstantPressureVolumeQuad :: update( ) 
{
  // strains ordered  00, 11, 22, 01  
  //            i.e.  11, 22, 33, 12 
  //
  //            strain(0) =   eps_00
  //            strain(1) =   eps_11
  //            strain(2) =   eps_22
  //            strain(3) = 2*eps_01
  //
  //  same ordering for stresses but no 2 

  int i,  k, l;
  int node ;
  int success = 0;
  
  static double tmp_shp[3][4] ; //shape functions

  static double shp[3][4][4] ; //shape functions at each gauss point

  static double vol_avg_shp[3][4] ; // volume averaged shape functions

  double xsj ;  // determinant jacaobian matrix 

  static Matrix sx(2,2) ; // inverse jacobian matrix 

  double dvol[4] ; //volume elements

  double volume = 0.0 ; //volume of element

  double theta = 0.0 ; //average volume change (trace of strain) 

  static Vector strain(4) ; //strain in vector form 

  double trace = 0.0 ; //trace of the strain 

  static Vector one(4) ; //rank 2 identity as a vector
  
  //one vector
  one(0) = 1.0 ;
  one(1) = 1.0 ;
  one(2) = 1.0 ;
  one(3) = 0.0 ;


  //zero stuff
  volume = 0.0 ;

  for ( k = 0; k < 3; k++ ){
    for ( l = 0; l < 4; l++ ) 
        vol_avg_shp[k][l] = 0.0 ; 
  } //end for k


  //gauss loop to compute volume averaged shape functions

  for ( i = 0; i < 4; i++ ){
    
    shape2d( sg[i], tg[i], xl, tmp_shp, xsj, sx ) ;

    dvol[i] = wg[i] * xsj ;  // multiply by radius for axisymmetry 

    volume += dvol[i] ;

    for ( k = 0; k < 3; k++ ){
      for ( l = 0; l < 4; l++ ) {

	shp[k][l][i] = tmp_shp[k][l] ;

        vol_avg_shp[k][l] += tmp_shp[k][l] * dvol[i] ;

      } // end for l
    } //end for k

  } //end for i 


  //compute volume averaged shape functions
  for ( k = 0; k < 3; k++ ){
    for ( l = 0; l < 4; l++ ) 
        vol_avg_shp[k][l] /= volume ; 
  } //end for k


  //compute theta
  theta = 0.0 ;
  for ( i = 0; i < 4; i++ ) {

    strain.Zero( ) ;

    //node loop to compute strain
    for ( node = 0; node < 4; node++ ) {

      const Vector &ul = nodePointers[node]->getTrialDisp( ) ;

      strain(0) += shp[0][node][i] * ul(0) ;

      strain(1) += shp[1][node][i] * ul(1) ;

      strain(2) = 0.0 ;  // not zero for axisymmetry

    } // end for node

    trace  =  strain(0) + strain(1) + strain(2) ;

    theta +=  trace * dvol[i] ;

  } // end for i
  theta /= volume ;


  //compute strain in materials
  for ( i = 0; i < 4; i++ ) {

    strain.Zero( ) ;

    //node loop to compute strain
    for ( node = 0; node < 4; node++ ) {

      const Vector &ul = nodePointers[node]->getTrialDisp( ) ;

      strain(0) += shp[0][node][i] * ul(0) ;

      strain(1) += shp[1][node][i] * ul(1) ;

      strain(2) = 0.0 ; // not zero for axisymmetry

      strain(3) +=  shp[1][node][i] * ul(0) 	
	          + shp[0][node][i] * ul(1) ; 

    } // end for node

    trace = strain(0) + strain(1) + strain(2) ;

    //strain -= (one3*trace)*one ;
    strain.addVector(1.0,  one, -one3*trace ) ;

    //strain += (one3*theta)*one ;
    strain.addVector(1.0,  one, one3*theta ) ;

    success += materialPointers[i]->setTrialStrain( strain ) ;

  } // end for i

  return success;
}

//print out element data
void ConstantPressureVolumeQuad :: Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "Four Node Quad -- Mixed Pressure/Volume -- Plane Strain \n" ;
  s << "Element Number " << this->getTag() << endln ;
  s << "Node 1 : " << connectedExternalNodes(0) << endln ;
  s << "Node 2 : " << connectedExternalNodes(1) << endln ;
  s << "Node 3 : " << connectedExternalNodes(2) << endln ;
  s << "Node 4 : " << connectedExternalNodes(3) << endln ;
  s << "Material Information : \n " ;

  materialPointers[0]->Print( s, flag ) ;

  s << endln ;
}

//return stiffness matrix 
const Matrix& ConstantPressureVolumeQuad :: getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    

const Matrix& ConstantPressureVolumeQuad :: getInitialStiff( ) 
{
  int i,  j,  k, l, p, q ;
  int jj, kk ;
  
  static double tmp_shp[3][4] ; //shape functions

  static double shp[3][4][4] ; //shape functions at each gauss point

  static double vol_avg_shp[3][4] ; // volume averaged shape functions

  double xsj ;  // determinant jacaobian matrix 

  static Matrix sx(2,2) ; // inverse jacobian matrix 

  double dvol[4] ; //volume elements

  double volume = 0.0 ; //volume of element

  double pressure = 0.0 ; //constitutive pressure  

  static Vector strain(4) ; //strain in vector form 

  // static Vector sigBar(4) ; //stress in vector form
  static Vector sig(4) ; //mixed stress in vector form

  double trace = 0.0 ; //trace of the strain 

  static Matrix BJtran(2,4) ; 
  static Matrix BK(4,2) ;

  static Matrix littleBJtran(2,1) ;
  static Matrix littleBK(1,2) ; 

  static Matrix stiffJK(2,2) ; //nodeJ-nodeK 2x2 stiffness
  static Vector residJ(2) ; //nodeJ residual 
  
  static Vector one(4) ; //rank 2 identity as a vector
  
  static Matrix Pdev(4,4) ; //deviator projector

  //  static Matrix dd(4,4) ;  //material tangent

  static Matrix ddPdev(4,4) ;
  static Matrix PdevDD(4,4) ;

  static double Pdev_dd_Pdev_data[16];
  static double Pdev_dd_one_data[4]; 
  static double one_dd_Pdev_data[4];
  static Matrix Pdev_dd_Pdev(Pdev_dd_Pdev_data, 4, 4);
  static Matrix Pdev_dd_one(Pdev_dd_one_data, 4, 1); 
  static Matrix one_dd_Pdev(one_dd_Pdev_data, 1,4) ;

  double bulk ;
  static Matrix BJtranD(2,4) ;
  static Matrix BJtranDone(2,1) ;

  static Matrix littleBJoneD(2,4) ;
  static Matrix littleBJtranBulk(2,1) ;
  
  //zero stiffness and residual 
  stiff.Zero();

  //one vector
  one(0) = 1.0 ;
  one(1) = 1.0 ;
  one(2) = 1.0 ;
  one(3) = 0.0 ;

  //Pdev matrix
  Pdev.Zero( ) ;

  Pdev(0,0) =  two3 ;
  Pdev(0,1) = -one3 ;
  Pdev(0,2) = -one3 ;

  Pdev(1,0) = -one3 ;
  Pdev(1,1) =  two3 ;
  Pdev(1,2) = -one3 ;

  Pdev(2,0) = -one3 ;
  Pdev(2,1) = -one3 ;
  Pdev(2,2) =  two3 ;

  Pdev(3,3) = 1.0 ;

  //zero stuff
  volume = 0.0 ;

  for ( k = 0; k < 3; k++ ){
    for ( l = 0; l < 4; l++ ) 
        vol_avg_shp[k][l] = 0.0 ; 
  } //end for k


  //gauss loop to compute volume averaged shape functions

  for ( i = 0; i < 4; i++ ){
    
    shape2d( sg[i], tg[i], xl, tmp_shp, xsj, sx ) ;

    dvol[i] = wg[i] * xsj ;  // multiply by radius for axisymmetry 

    volume += dvol[i] ;

    for ( k = 0; k < 3; k++ ){
      for ( l = 0; l < 4; l++ ) {

	shp[k][l][i] = tmp_shp[k][l] ;

        vol_avg_shp[k][l] += tmp_shp[k][l] * dvol[i] ;

      } // end for l
    } //end for k

  } //end for i 


  //compute volume averaged shape functions
  for ( k = 0; k < 3; k++ ){
    for ( l = 0; l < 4; l++ ) 
        vol_avg_shp[k][l] /= volume ; 
  } //end for k

  //residual and tangent calculations gauss loop
  for ( i = 0; i < 4; i++ ) {

    static Matrix dd(4,4);

    dd = materialPointers[i]->getInitialTangent( ) ;

    dd *= dvol[i] ;
    
    //Pdev_dd_Pdev = Pdev * dd * Pdev ;
    Pdev_dd_Pdev.addMatrixTripleProduct(0.0,
					Pdev,
					dd,
					1.0) ;
      
    //Pdev_dd_one  = one3 * ( Pdev * dd * oneMatrix ) ;
    PdevDD.addMatrixProduct(0.0, Pdev, dd, 1.0) ;
    Pdev_dd_one(0,0) = one3 * (PdevDD(0,0) + PdevDD(0,1) + PdevDD(0,2));
    Pdev_dd_one(1,0) = one3 * (PdevDD(1,0) + PdevDD(1,1) + PdevDD(1,2));
    Pdev_dd_one(2,0) = one3 * (PdevDD(2,0) + PdevDD(2,1) + PdevDD(2,2));
    Pdev_dd_one(3,0) = one3 * (PdevDD(3,0) + PdevDD(3,1) + PdevDD(3,2));
    
    //one_dd_Pdev  = one3 * ( oneTran * dd * Pdev ) ;
    ddPdev.addMatrixProduct(0.0, dd, Pdev, 1.0) ;
    one_dd_Pdev(0,0) = one3 * (ddPdev(0,0) + ddPdev(1,0) + ddPdev(2,0));
    one_dd_Pdev(0,1) = one3 * (ddPdev(0,1) + ddPdev(1,1) + ddPdev(2,1));
    one_dd_Pdev(0,2) = one3 * (ddPdev(0,2) + ddPdev(1,2) + ddPdev(2,2));
    one_dd_Pdev(0,3) = one3 * (ddPdev(0,3) + ddPdev(1,3) + ddPdev(2,3));
    
    bulk = one9 * ( dd(0,0) + dd(0,1) + dd(0,2)  
		    + dd(1,0) + dd(1,1) + dd(1,2) 
		    + dd(2,0) + dd(2,1) + dd(2,2) ) ;
    
    jj = 0 ;
    for ( j = 0; j < 4; j++ ) {
      
      double BJ00 = shp[0][j][i];
      double BJ11 = shp[1][j][i]; 
      double BJ30 = shp[1][j][i];
      double BJ31 = shp[0][j][i];
      
      BJtran.Zero( );
      BJtran(0,0) = shp[0][j][i] ;
      BJtran(1,1) = shp[1][j][i] ; 
      
      // BJ(2,0) for axi-symmetry 
      
      BJtran(0,3) = shp[1][j][i]  ;
      BJtran(1,3) = shp[0][j][i]  ;
      
      //compute residual 
      
      double ltBJ00 = vol_avg_shp[0][j] ;
      double ltBJ01 = vol_avg_shp[1][j] ;
      
      //BJtranD          =  BJtran * Pdev_dd_Pdev ;
      // BJtranD.addMatrixProduct(0.0,  BJtran, Pdev_dd_Pdev, 1.0);
      
      //littleBJoneD     =  littleBJtran * one_dd_Pdev ;
      // littleBJoneD.addMatrixProduct(0.0,  littleBJtran, one_dd_Pdev, 1.0);
      
      static double Adata[8];
      static Matrix A(Adata, 2, 4);
      
      // A = BJtranD;
      // A += littleBJoneD;
      
      for (int colA = 0, loc = 0, colPdev = 0; colA<4; colA++, colPdev += 4) {
	double data3colA = Pdev_dd_Pdev_data[3+colPdev];
	Adata[loc++] = BJ00*Pdev_dd_Pdev_data[colPdev] + BJ30*data3colA + ltBJ00*one_dd_Pdev_data[colA];
	Adata[loc++] = BJ11*Pdev_dd_Pdev_data[1+colPdev] + BJ31*data3colA + ltBJ01*one_dd_Pdev_data[colA];
      }
      
      //BJtranDone       =  BJtran * Pdev_dd_one ;
      // BJtranDone.addMatrixProduct(0.0,  BJtran, Pdev_dd_one, 1.0);
      
      //littleBJtranBulk =  bulk * littleBJtran ;
      // littleBJtranBulk = littleBJtran ;
      // littleBJtranBulk *= bulk ;
      
      double B1, B2;
      // B1 = BJtranDone(0,0) + littleBJtranBulk(0,0);
      // B2 = BJtranDone(1,0) + littleBJtranBulk(1,0);
      
      B1 = BJ00*Pdev_dd_one_data[0] + BJ30*Pdev_dd_one_data[3] + ltBJ00 *bulk;
      B2 = BJ11*Pdev_dd_one_data[1] + BJ31*Pdev_dd_one_data[3] + ltBJ01 *bulk;
      
      int colkk, colkkP1;
      for ( k = 0, kk=0, colkk =0, colkkP1 =8; 
	    k < 4; 
	    k++, kk += 2, colkk += 16, colkkP1 += 16 ) {
	
	double BK00 = shp[0][k][i];
	double BK11 = shp[1][k][i];
	double BK30 = shp[1][k][i];
	double BK31 = shp[0][k][i];
	
	double littleBK00 = vol_avg_shp[0][k];
	double littleBK01 = vol_avg_shp[1][k];

	//compute stiffness matrix
        
	stiff( jj,   kk   ) += Adata[0]*BK00 + Adata[6]*BK30 + B1 * littleBK00;
	stiff( jj+1, kk   ) += Adata[1]*BK00 + Adata[7]*BK30 + B2 * littleBK00;
	stiff( jj,   kk+1 ) += Adata[2]*BK11 + Adata[6]*BK31 + B1 * littleBK01;
	stiff( jj+1, kk+1 ) += Adata[3]*BK11 + Adata[7]*BK31 + B2 * littleBK01;
	
      } // end for k
      
      jj += 2 ;
    } // end for j 
  } //end for i

  return stiff ;
}    

//return mass matrix
const Matrix& ConstantPressureVolumeQuad :: getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;
  return mass ;
} 



void ConstantPressureVolumeQuad :: zeroLoad( )
{
  if (load != 0)
    load->Zero();

  return ;
}

int 
ConstantPressureVolumeQuad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "ConstantPressureVolumeQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}

int
ConstantPressureVolumeQuad::addInertiaLoadToUnbalance(const Vector &accel)
{
  static const int numberGauss = 9 ;
  static const int numberNodes = 9 ;
  static const int ndf = 2 ; 

  int i;

  // check to see if have mass
  int haveRho = 0;
  for (i = 0; i < numberGauss; i++) {
    if (materialPointers[i]->getRho() != 0.0)
      haveRho = 1;
  }

  if (haveRho == 0)
    return 0;

  // Compute mass matrix
  int tangFlag = 1 ;
  formInertiaTerms( tangFlag ) ;

  // store computed RV fro nodes in resid vector
  int count = 0;

  for (i=0; i<numberNodes; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int j=0; j<ndf; j++)
      resid(count++) = Raccel(i);
  }

  // create the load vector if one does not exist
  if (load == 0) 
    load = new Vector(numberNodes*ndf);

  // add -M * RV(accel) to the load vector
  load->addMatrixVector(1.0, mass, resid, -1.0);
  
  return 0;
}


//get residual
const Vector& ConstantPressureVolumeQuad :: getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  // subtract external loads 
  if (load != 0)
    resid -= *load;

  return resid ;   
}


//get residual with inertia terms
const Vector& ConstantPressureVolumeQuad :: getResistingForceIncInertia( )
{
  int tang_flag = 0 ; //don't get the tangent

  static Vector res(8);

  //do tangent and residual here 
  formResidAndTangent( tang_flag ) ;

  //inertia terms
  formInertiaTerms( tang_flag ) ;
  res = resid;

  // subtract external loads 
  if (load != 0)
    res -= *load;

  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();

  return res ;
}

//*****************************************************************************
//form inertia terms

void   ConstantPressureVolumeQuad::formInertiaTerms( int tangFlag ) 
{

  static const int ndm = 2 ;

  static const int ndf = 2 ; 

  static const int numberNodes = 4 ;

  static const int numberGauss = 4 ;

  static const int nShape = 3 ;

  static const int massIndex = nShape - 1 ;

  double xsj ;  // determinant jacaobian matrix 

  double dvol ; //volume element

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static Vector momentum(ndf) ;

  static Matrix sx(ndm,ndm) ;

  int i, j, k, p ;
  int jj, kk ;

  double temp, rho, massJK ;


  //zero mass 
  mass.Zero( ) ;

  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) {

    //get shape functions    
    shape2d( sg[i], tg[i], xl, shp, xsj, sx ) ;
    
    //volume element
    dvol = wg[i] * xsj ;

    //node loop to compute acceleration
    momentum.Zero( ) ;
    for ( j = 0; j < numberNodes; j++ ) 
      //momentum += shp[massIndex][j] * ( nodePointers[j]->getTrialAccel()  ) ; 
      momentum.addVector( 1.0,
			  nodePointers[j]->getTrialAccel(),
			  shp[massIndex][j] ) ;

    //density
    rho = materialPointers[i]->getRho() ;

    //multiply acceleration by density to form momentum
    momentum *= rho ;


    //residual and tangent calculations node loops
    jj = 0 ;
    for ( j = 0; j < numberNodes; j++ ) {

      temp = shp[massIndex][j] * dvol ;

      if ( tangFlag == 1 ) {

	 //multiply by density
	 temp *= rho ;

	 //node-node mass
         kk = 0 ;
         for ( k = 0; k < numberNodes; k++ ) {

	    massJK = temp * shp[massIndex][k] ;

            for ( p = 0; p < ndf; p++ )  
	      mass( jj+p, kk+p ) += massJK ;
            
            kk += ndf ;
          } // end for k loop

      } // end if tang_flag 

      else
	for ( p = 0; p < ndf; p++ )
	  resid( jj+p ) += ( temp * momentum(p) )  ;

      jj += ndf ;

    } // end for j loop


  } //end for i gauss loop 



}

//*********************************************************************

//form residual and tangent
void ConstantPressureVolumeQuad ::  formResidAndTangent( int tang_flag ) 
{
  // strains ordered  00, 11, 22, 01  
  //            i.e.  11, 22, 33, 12 
  //
  //            strain(0) =   eps_00
  //            strain(1) =   eps_11
  //            strain(2) =   eps_22
  //            strain(3) = 2*eps_01
  //
  //  same ordering for stresses but no 2 

  int i,  j,  k, l, p, q ;
  int jj, kk ;
  
  static double tmp_shp[3][4] ; //shape functions

  static double shp[3][4][4] ; //shape functions at each gauss point

  static double vol_avg_shp[3][4] ; // volume averaged shape functions

  double xsj ;  // determinant jacaobian matrix 

  static Matrix sx(2,2) ; // inverse jacobian matrix 

  double dvol[4] ; //volume elements

  double volume = 0.0 ; //volume of element

  double pressure = 0.0 ; //constitutive pressure  

  static Vector strain(4) ; //strain in vector form 

  // static Vector sigBar(4) ; //stress in vector form
  static Vector sig(4) ; //mixed stress in vector form

  double trace = 0.0 ; //trace of the strain 

  static Matrix BJtran(2,4) ; 
  static Matrix BK(4,2) ;

  static Matrix littleBJtran(2,1) ;
  static Matrix littleBK(1,2) ; 

  static Matrix stiffJK(2,2) ; //nodeJ-nodeK 2x2 stiffness
  static Vector residJ(2) ; //nodeJ residual 
  
  static Vector one(4) ; //rank 2 identity as a vector
  
  static Matrix Pdev(4,4) ; //deviator projector

  //  static Matrix dd(4,4) ;  //material tangent

  static Matrix ddPdev(4,4) ;
  static Matrix PdevDD(4,4) ;

  static double Pdev_dd_Pdev_data[16];
  static double Pdev_dd_one_data[4]; 
  static double one_dd_Pdev_data[4];
  static Matrix Pdev_dd_Pdev(Pdev_dd_Pdev_data, 4, 4);
  static Matrix Pdev_dd_one(Pdev_dd_one_data, 4, 1); 
  static Matrix one_dd_Pdev(one_dd_Pdev_data, 1,4) ;

  double bulk ;
  static Matrix BJtranD(2,4) ;
  static Matrix BJtranDone(2,1) ;

  static Matrix littleBJoneD(2,4) ;
  static Matrix littleBJtranBulk(2,1) ;
  
  //zero stiffness and residual 
  if ( tang_flag == 1 ) 
    stiff.Zero();
  else
    resid.Zero();

  //one vector
  one(0) = 1.0 ;
  one(1) = 1.0 ;
  one(2) = 1.0 ;
  one(3) = 0.0 ;

  //Pdev matrix
  Pdev.Zero( ) ;

  Pdev(0,0) =  two3 ;
  Pdev(0,1) = -one3 ;
  Pdev(0,2) = -one3 ;

  Pdev(1,0) = -one3 ;
  Pdev(1,1) =  two3 ;
  Pdev(1,2) = -one3 ;

  Pdev(2,0) = -one3 ;
  Pdev(2,1) = -one3 ;
  Pdev(2,2) =  two3 ;

  Pdev(3,3) = 1.0 ;

  //zero stuff
  volume = 0.0 ;

  for ( k = 0; k < 3; k++ ){
    for ( l = 0; l < 4; l++ ) 
        vol_avg_shp[k][l] = 0.0 ; 
  } //end for k


  //gauss loop to compute volume averaged shape functions

  for ( i = 0; i < 4; i++ ){
    
    shape2d( sg[i], tg[i], xl, tmp_shp, xsj, sx ) ;

    dvol[i] = wg[i] * xsj ;  // multiply by radius for axisymmetry 

    volume += dvol[i] ;

    for ( k = 0; k < 3; k++ ){
      for ( l = 0; l < 4; l++ ) {

	shp[k][l][i] = tmp_shp[k][l] ;

        vol_avg_shp[k][l] += tmp_shp[k][l] * dvol[i] ;

      } // end for l
    } //end for k

  } //end for i 


  //compute volume averaged shape functions
  for ( k = 0; k < 3; k++ ){
    for ( l = 0; l < 4; l++ ) 
        vol_avg_shp[k][l] /= volume ; 
  } //end for k

  //compute pressure if residual calculation
  if (tang_flag != 1) {
    pressure = 0.0 ;
    for ( i = 0; i < 4; i++ ) {

      const Vector &sigBar = materialPointers[i]->getStress( ) ;

      pressure +=  one3 * ( sigBar(0) + sigBar(1) + sigBar(2) ) * dvol[i] ;
      
    } // end for i

    pressure /= volume ;
  } // end if != tang_flag

  //residual and tangent calculations gauss loop
  
  for ( i = 0; i < 4; i++ ) {

    if ( tang_flag == 1 ) {    // compute matrices for stiffness calculation
      
      static Matrix dd(4,4);
      dd = materialPointers[i]->getTangent( ) ;
      
      dd *= dvol[i] ;
      
      //Pdev_dd_Pdev = Pdev * dd * Pdev ;
      Pdev_dd_Pdev.addMatrixTripleProduct(0.0,
					  Pdev,
					  dd,
					  1.0) ;
      
      //Pdev_dd_one  = one3 * ( Pdev * dd * oneMatrix ) ;
      PdevDD.addMatrixProduct(0.0, Pdev, dd, 1.0) ;
      Pdev_dd_one(0,0) = one3 * (PdevDD(0,0) + PdevDD(0,1) + PdevDD(0,2));
      Pdev_dd_one(1,0) = one3 * (PdevDD(1,0) + PdevDD(1,1) + PdevDD(1,2));
      Pdev_dd_one(2,0) = one3 * (PdevDD(2,0) + PdevDD(2,1) + PdevDD(2,2));
      Pdev_dd_one(3,0) = one3 * (PdevDD(3,0) + PdevDD(3,1) + PdevDD(3,2));
      
      //one_dd_Pdev  = one3 * ( oneTran * dd * Pdev ) ;
      ddPdev.addMatrixProduct(0.0, dd, Pdev, 1.0) ;
      one_dd_Pdev(0,0) = one3 * (ddPdev(0,0) + ddPdev(1,0) + ddPdev(2,0));
      one_dd_Pdev(0,1) = one3 * (ddPdev(0,1) + ddPdev(1,1) + ddPdev(2,1));
      one_dd_Pdev(0,2) = one3 * (ddPdev(0,2) + ddPdev(1,2) + ddPdev(2,2));
      one_dd_Pdev(0,3) = one3 * (ddPdev(0,3) + ddPdev(1,3) + ddPdev(2,3));
      
      bulk = one9 * ( dd(0,0) + dd(0,1) + dd(0,2)  
		      + dd(1,0) + dd(1,1) + dd(1,2) 
		      + dd(2,0) + dd(2,1) + dd(2,2) ) ;
      
    } else { // compute stress for residual calculation
      //stress for equilibrium
      const Vector &sigBar = materialPointers[i]->getStress( ) ; 
      trace = sigBar(0) + sigBar(1) + sigBar(2) ;
      sig  = sigBar ;

      //sig -= (one3*trace)*one ;
      sig.addVector(1.0,  one, -one3*trace ) ;
      sig.addVector(1.0,  one, pressure ) ;
      
      //multilply by volume elements and compute 
      sig *= dvol[i] ;
    }
        
    //residual and tangent loop over nodes
    jj = 0 ;
    for ( j = 0; j < 4; j++ ) {

      /********** expanding for efficiency the matrix operations that use these
      BJ.Zero( );
      BJ(0,0) = shp[0][j][i] ;
      BJ(1,1) = shp[1][j][i] ; 

      // BJ(2,0) for axi-symmetry 

      BJ(3,0) = shp[1][j][i]  ;
      BJ(3,1) = shp[0][j][i]  ;

      littleBJ(0,0) = vol_avg_shp[0][j] ;
      littleBJ(0,1) = vol_avg_shp[1][j] ;

      // BJtran = this->transpose( 4, 2, BJ ) ;
      for (p=0; p<2; p++) {
	for (q=0; q<4; q++) 
	  BJtran(p,q) = BJ(q,p) ;
      }//end for p

      for (p=0; p<2; p++) {
        for (q=0; q<1; q++) 
          littleBJtran(p,q) = littleBJ(q,p) ;
      }//end for p

      **********************************************************************/


      double BJ00 = shp[0][j][i];
      double BJ11 = shp[1][j][i]; 
      double BJ30 = shp[1][j][i];
      double BJ31 = shp[0][j][i];

      BJtran.Zero( );
      BJtran(0,0) = shp[0][j][i] ;
      BJtran(1,1) = shp[1][j][i] ; 

      // BJ(2,0) for axi-symmetry 

      BJtran(0,3) = shp[1][j][i]  ;
      BJtran(1,3) = shp[0][j][i]  ;

      //compute residual 

      if ( tang_flag == 1 ) { //stiffness matrix

	double ltBJ00 = vol_avg_shp[0][j] ;
	double ltBJ01 = vol_avg_shp[1][j] ;

	//BJtranD          =  BJtran * Pdev_dd_Pdev ;
	// BJtranD.addMatrixProduct(0.0,  BJtran, Pdev_dd_Pdev, 1.0);

	//littleBJoneD     =  littleBJtran * one_dd_Pdev ;
	// littleBJoneD.addMatrixProduct(0.0,  littleBJtran, one_dd_Pdev, 1.0);

	static double Adata[8];
	static Matrix A(Adata, 2, 4);

	// A = BJtranD;
	// A += littleBJoneD;

	for (int colA = 0, loc = 0, colPdev = 0; colA<4; colA++, colPdev += 4) {
	  double data3colA = Pdev_dd_Pdev_data[3+colPdev];
	  Adata[loc++] = BJ00*Pdev_dd_Pdev_data[colPdev] + BJ30*data3colA + ltBJ00*one_dd_Pdev_data[colA];
	  Adata[loc++] = BJ11*Pdev_dd_Pdev_data[1+colPdev] + BJ31*data3colA + ltBJ01*one_dd_Pdev_data[colA];
	}

        //BJtranDone       =  BJtran * Pdev_dd_one ;
	// BJtranDone.addMatrixProduct(0.0,  BJtran, Pdev_dd_one, 1.0);

	//littleBJtranBulk =  bulk * littleBJtran ;
	// littleBJtranBulk = littleBJtran ;
	// littleBJtranBulk *= bulk ;

	double B1, B2;
	// B1 = BJtranDone(0,0) + littleBJtranBulk(0,0);
	// B2 = BJtranDone(1,0) + littleBJtranBulk(1,0);

	B1 = BJ00*Pdev_dd_one_data[0] + BJ30*Pdev_dd_one_data[3] + ltBJ00 *bulk;
	B2 = BJ11*Pdev_dd_one_data[1] + BJ31*Pdev_dd_one_data[3] + ltBJ01 *bulk;

	int colkk, colkkP1;
	for ( k = 0, kk=0, colkk =0, colkkP1 =8; 
	      k < 4; 
	      k++, kk += 2, colkk += 16, colkkP1 += 16 ) {

	  /**************************************************************
	   REPLACING THESE LINES WITH THE 4 BELOW COMMENT FOR EFFICIENCY
	   BK.Zero( );
	   BK(0,0) = shp[0][k][i];
	   BK(1,1) = shp[1][k][i];
	   
	   // BK(2,0) for axi-symmetry 
	   
	   BK(3,0) = shp[1][k][i];
	   BK(3,1) = shp[0][k][i];
	  **************************************************************/

	  double BK00 = shp[0][k][i];
	  double BK11 = shp[1][k][i];
	  double BK30 = shp[1][k][i];
	  double BK31 = shp[0][k][i];

	  double littleBK00 = vol_avg_shp[0][k];
	  double littleBK01 = vol_avg_shp[1][k];

	  //compute stiffness matrix
        
	  // stiffJK =  ( BJtranD + littleBJoneD ) * BK
	  //        +  ( BJtranDone + littleBJtranBulk ) * littleBK ; 

	  /**************************************************************
	    REPLACING THESE LINES WITH THE 4 BELOW COMMENT FOR EFFICIENCY
	  //stiffJK.addMatrixProduct(0.0, A, BK, 1.0);

	  //stiff( jj,   kk   ) += stiffJK(0,0) + B1 * littleBK00;
	  //stiff( jj+1, kk   ) += stiffJK(1,0) + B2 * littleBK00;
	  //stiff( jj,   kk+1 ) += stiffJK(0,1) + B1 * littleBK01;
	  //stiff( jj+1, kk+1 ) += stiffJK(1,1) + B2 * littleBK01;	  
	  ***************************************************************/

	  // matrixData[  colkk +   jj] += Adata[0]*BK00 + Adata[6]*BK30 + B1 * littleBK00;
	  // matrixData[  colkk + jj+1] += Adata[1]*BK00 + Adata[7]*BK30 + B2 * littleBK00;
	  // matrixData[colkkP1 +   jj] += Adata[2]*BK11 + Adata[6]*BK31 + B1 * littleBK01;
	  // matrixData[colkkP1 + jj+1] += Adata[3]*BK11 + Adata[7]*BK31 + B2 * littleBK01;
	  stiff( jj,   kk   ) += Adata[0]*BK00 + Adata[6]*BK30 + B1 * littleBK00;
	  stiff( jj+1, kk   ) += Adata[1]*BK00 + Adata[7]*BK30 + B2 * littleBK00;
	  stiff( jj,   kk+1 ) += Adata[2]*BK11 + Adata[6]*BK31 + B1 * littleBK01;
	  stiff( jj+1, kk+1 ) += Adata[3]*BK11 + Adata[7]*BK31 + B2 * littleBK01;

	} // end for k

      } else { // residual calculation
	
	//residJ = BJtran * sig; 
	residJ.addMatrixVector(0.0,  BJtran, sig, 1.0);
	resid( jj   ) += residJ(0);
	resid( jj+1 ) += residJ(1);
      }
      
      jj += 2 ;
    } // end for j 

  } //end for i

  
  return ;
}


//shape function routine for four node quads
void ConstantPressureVolumeQuad :: shape2d( double ss, double tt, 
		                            const double x[2][4], 
		                            double shp[3][4], 
		                            double &xsj, 
		                            Matrix &sx ) 
{ 

  int i, j, k ;

  double temp ;
     
  static const double s[] = { -0.5,  0.5, 0.5, -0.5 } ;
  static const double t[] = { -0.5, -0.5, 0.5,  0.5 } ;

  static double xs[2][2];

  //  static Matrix xs(2,2) ;

  for ( i = 0; i < 4; i++ ) {
      shp[2][i] = ( 0.5 + s[i]*ss )*( 0.5 + t[i]*tt ) ;
      shp[0][i] = s[i] * ( 0.5 + t[i]*tt ) ;
      shp[1][i] = t[i] * ( 0.5 + s[i]*ss ) ;
  } // end for i

  
  // Construct jacobian and its inverse
  
  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++ ) {

      double value = 0;
      for ( k = 0; k < 4; k++ )
	value +=  x[i][k] * shp[j][k] ;
      // xs(i,j) +=  x[i][k] * shp[j][k] ;
      xs[i][j] = value;
      
    } //end for j
  }  // end for i 

  xsj = xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0] ;

  sx(0,0) =  xs[1][1] / xsj ;
  sx(1,1) =  xs[0][0] / xsj ;
  sx(0,1) = -xs[0][1] / xsj ; 
  sx(1,0) = -xs[1][0] / xsj ; 

  //form global derivatives 
  
  for ( i = 0; i < 4; i++ ) {
    temp      = shp[0][i]*sx(0,0) + shp[1][i]*sx(1,0) ;
    shp[1][i] = shp[0][i]*sx(0,1) + shp[1][i]*sx(1,1) ;
    shp[0][i] = temp ;
  } // end for i

  return ;
}
	
//***********************************************************************

int
ConstantPressureVolumeQuad::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    const Vector &end1Crd = nodePointers[0]->getCrds();
    const Vector &end2Crd = nodePointers[1]->getCrds();	
    const Vector &end3Crd = nodePointers[2]->getCrds();	
    const Vector &end4Crd = nodePointers[3]->getCrds();	

    static Matrix coords(4,3) ;
    static Vector values(4) ;

    coords.Zero( ) ;

    values(0) = 1 ;
    values(1) = 1 ;
    values(2) = 1 ;
    values(3) = 1 ;


    if (displayMode >= 0) {    
      
      const Vector &end1Disp = nodePointers[0]->getDisp();
      const Vector &end2Disp = nodePointers[1]->getDisp();
      const Vector &end3Disp = nodePointers[2]->getDisp();
      const Vector &end4Disp = nodePointers[3]->getDisp();

      for (int i = 0; i < 2; i++) {
	coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
	coords(1,i) = end2Crd(i) + end2Disp(i)*fact;    
	coords(2,i) = end3Crd(i) + end3Disp(i)*fact;    
	coords(3,i) = end4Crd(i) + end4Disp(i)*fact;    
      }
    } else {
      int mode = displayMode  *  -1;
      const Matrix &eigen1 = nodePointers[0]->getEigenvectors();
      const Matrix &eigen2 = nodePointers[1]->getEigenvectors();
      const Matrix &eigen3 = nodePointers[2]->getEigenvectors();
      const Matrix &eigen4 = nodePointers[3]->getEigenvectors();
      if (eigen1.noCols() >= mode) {
	for (int i = 0; i < 2; i++) {
	  coords(0,i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	  coords(1,i) = end2Crd(i) + eigen2(i,mode-1)*fact;
	  coords(2,i) = end3Crd(i) + eigen3(i,mode-1)*fact;
	  coords(3,i) = end4Crd(i) + eigen4(i,mode-1)*fact;
	}    
      } else {
	for (int i = 0; i < 2; i++) {
	  coords(0,i) = end1Crd(i);
	  coords(1,i) = end2Crd(i);
	  coords(2,i) = end3Crd(i);
	  coords(3,i) = end4Crd(i);
	}    
      }
    }

    //opserr << coords;
    int error = 0;

    error += theViewer.drawPolygon (coords, values);

    return error;
}


Response*
ConstantPressureVolumeQuad::setResponse(const char **argv, int argc, 
					Information &eleInfo, OPS_Stream &output)
{
  Response *theResponse =0;

  output.tag("ElementOutput");
  output.attr("eleType","ConstantPressureVolumeQuad");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  output.attr("node3",connectedExternalNodes[2]);
  output.attr("node4",connectedExternalNodes[3]);

  char dataOut[10];
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
    
    for (int i=1; i<=4; i++) {
      sprintf(dataOut,"P1_%d",i);
      output.tag("ResponseType",dataOut);
      sprintf(dataOut,"P2_%d",i);
      output.tag("ResponseType",dataOut);
    }
    
    theResponse =  new ElementResponse(this, 1, resid);
  }   else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",sg[pointNum-1]);
      output.attr("neta",tg[pointNum-1]);

      theResponse =  materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, eleInfo, output);
      
      output.endTag();

  } else if (strcmp(argv[0],"stresses") ==0) {

      for (int i=0; i<4; i++) {
	output.tag("GaussPoint");
	output.attr("number",i+1);
	output.attr("eta",sg[i]);
	output.attr("neta",tg[i]);

	output.tag("NdMaterialOutput");
	output.attr("classType", materialPointers[i]->getClassTag());
	output.attr("tag", materialPointers[i]->getTag());

	output.tag("ResponseType","UnknownStress");
	output.tag("ResponseType","UnknownStress");
	output.tag("ResponseType","UnknownStress");
	output.tag("ResponseType","UnknownStress");

	output.endTag(); // GaussPoint
	output.endTag(); // NdMaterialOutput
      }

      theResponse =  new ElementResponse(this, 3, Vector(16));
    }
  }
	
  output.endTag(); // ElementOutput

  return theResponse;
}

int 
ConstantPressureVolumeQuad::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {

    // Loop over the integration points
    static Vector stresses(16);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStress();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      stresses(cnt+3) = sigma(2);
      cnt += 4;
    }
    return eleInfo.setVector(stresses);
	
  } else

    return -1;
}
   
//-----------------------------------------------------------------------

int ConstantPressureVolumeQuad :: sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(13);
  
  int i;
  for (i = 0; i < 4; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i+4) = matDbTag;
  }
  
  idData(8) = this->getTag();
  idData(9) = connectedExternalNodes(0);
  idData(10) = connectedExternalNodes(1);
  idData(11) = connectedExternalNodes(2);
  idData(12) = connectedExternalNodes(3);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ConstantPressureVolumeQuad::sendSelf() - " << this->getTag() << " failed to send ID\n"; 
			    
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING ConstantPressureVolumeQuad::sendSelf() - " << this->getTag() << "failed to send its Material\n";
      return res;
    }
  }
  return res;
}
    
int 
ConstantPressureVolumeQuad :: recvSelf (int commitTag, 
					Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(13);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ConstantPressureVolumeQuad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(8));
  connectedExternalNodes(0) = idData(9);
  connectedExternalNodes(1) = idData(10);
  connectedExternalNodes(2) = idData(11);
  connectedExternalNodes(3) = idData(12);

  int i;

  if (materialPointers[0] == 0) {
    for (i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
	opserr << "ConstantPressureVolumeQuad::recvSelf() - " <<
	  "Broker could not create NDMaterial of class type" << matClassTag << endln;
	return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "NLBeamColumn3d::recvSelf() - material " << 
	  i << "failed to recv itself\n";
	return res;
      }
    }
  }
  // Number of materials is the same, receive materials into current space
  else {
    for (i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
	delete materialPointers[i];
	materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
	if (materialPointers[i] == 0) {
	  opserr << "ConstantPressureVolumeQuad::recvSelf() - " << 
	    "Broker could not create NDMaterial of class type" << matClassTag << endln;
	exit(-1);
	}
      }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "ConstantPressureVolumeQuad::recvSelf() - material  " << i <<
	  "failed to recv itself\n";
	return res;
      }
    }
  }
  
  return res;
}

