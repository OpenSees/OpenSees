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
                                                                        
// $Revision: 1.14 $
// $Date: 2006-03-21 22:19:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/ShellMITC4.cpp,v $

// Ed "C++" Love
//
// B-bar four node shell element with membrane and drill
//

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <ShellMITC4.h>
#include <R3vectors.h>
#include <Renderer.h>


#include <Channel.h>
#include <FEM_ObjectBroker.h>

#define min(a,b) ( (a)<(b) ? (a):(b) )

//static data
Matrix  ShellMITC4::stiff(24,24) ;
Vector  ShellMITC4::resid(24) ;
Matrix  ShellMITC4::mass(24,24) ;

//intialize pointers to zero using intialization list 
Matrix**  ShellMITC4::GammaB1pointer = 0 ;
Matrix**  ShellMITC4::GammaD1pointer = 0 ;
Matrix**  ShellMITC4::GammaA2pointer = 0 ;
Matrix**  ShellMITC4::GammaC2pointer = 0 ; 
Matrix**  ShellMITC4::Bhat           = 0 ;
    
//quadrature data
const double  ShellMITC4::root3 = sqrt(3.0) ;
const double  ShellMITC4::one_over_root3 = 1.0 / root3 ;

double ShellMITC4::sg[4] ;
double ShellMITC4::tg[4] ;
double ShellMITC4::wg[4] ;

/*
const double  ShellMITC4::sg[] = { -one_over_root3,  
				one_over_root3, 
				one_over_root3, 
	                       -one_over_root3 } ;

const double  ShellMITC4::tg[] = { -one_over_root3, 
                               -one_over_root3, 
                                one_over_root3,  
                                one_over_root3 } ;

const double  ShellMITC4::wg[] = { 1.0, 1.0, 1.0, 1.0 } ;
*/
  

//null constructor
ShellMITC4::ShellMITC4( ) :
Element( 0, ELE_TAG_ShellMITC4 ),
connectedExternalNodes(4), load(0), Ki(0)
{ 
  for (int i = 0 ;  i < 4; i++ ) 
    materialPointers[i] = 0;

  //shear matrix pointers
  if ( GammaB1pointer == 0 ) {
	GammaB1pointer = new Matrix*[4] ;      //four matrix pointers
	GammaB1pointer[0] = new Matrix(1,3) ;  //
	GammaB1pointer[1] = new Matrix(1,3) ;  //    four
	GammaB1pointer[2] = new Matrix(1,3) ;  //  1x3 matrices
	GammaB1pointer[3] = new Matrix(1,3) ;  //
  } //end if B1

  if ( GammaD1pointer == 0 ) {
	GammaD1pointer = new Matrix*[4] ;
	GammaD1pointer[0] = new Matrix(1,3) ;
	GammaD1pointer[1] = new Matrix(1,3) ;
	GammaD1pointer[2] = new Matrix(1,3) ;
	GammaD1pointer[3] = new Matrix(1,3) ;
  } //end if D1

  if ( GammaA2pointer == 0 ) {
	GammaA2pointer = new Matrix*[4] ;
	GammaA2pointer[0] = new Matrix(1,3) ;
	GammaA2pointer[1] = new Matrix(1,3) ;
	GammaA2pointer[2] = new Matrix(1,3) ;
	GammaA2pointer[3] = new Matrix(1,3) ;
  } //end if A2

  if ( GammaC2pointer == 0 ) {
	GammaC2pointer = new Matrix*[4] ;
	GammaC2pointer[0] = new Matrix(1,3) ;
	GammaC2pointer[1] = new Matrix(1,3) ;
	GammaC2pointer[2] = new Matrix(1,3) ;
	GammaC2pointer[3] = new Matrix(1,3) ;
  } //end if C2

  if ( Bhat == 0 ) {
	Bhat = new Matrix*[4] ;
	Bhat[0] = new Matrix(2,3) ;
	Bhat[1] = new Matrix(2,3) ;
	Bhat[2] = new Matrix(2,3) ;
	Bhat[3] = new Matrix(2,3) ;
  } //end if Bhat

  sg[0] = -one_over_root3;
  sg[1] = one_over_root3;
  sg[2] = one_over_root3;
  sg[3] = -one_over_root3;  

  tg[0] = -one_over_root3;
  tg[1] = -one_over_root3;
  tg[2] = one_over_root3;
  tg[3] = one_over_root3;  

  wg[0] = 1.0;
  wg[1] = 1.0;
  wg[2] = 1.0;
  wg[3] = 1.0;
}


//*********************************************************************
//full constructor
ShellMITC4::ShellMITC4(  int tag, 
                         int node1,
                         int node2,
   	                 int node3,
                         int node4,
	                 SectionForceDeformation &theMaterial ) :
Element( tag, ELE_TAG_ShellMITC4 ),
connectedExternalNodes(4), load(0), Ki(0)
{
  int i;

  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;

  for ( i = 0 ;  i < 4; i++ ) {

      materialPointers[i] = theMaterial.getCopy( ) ;

      if (materialPointers[i] == 0) {

	opserr << "ShellMITC4::constructor - failed to get a material of type: ShellSection\n";
      } //end if
      
  } //end for i 

  //shear matrix pointers
  if ( GammaB1pointer == 0 ) {
	GammaB1pointer = new Matrix*[4] ;      //four matrix pointers
	GammaB1pointer[0] = new Matrix(1,3) ;  //
	GammaB1pointer[1] = new Matrix(1,3) ;  //    four
	GammaB1pointer[2] = new Matrix(1,3) ;  //  1x3 matrices
	GammaB1pointer[3] = new Matrix(1,3) ;  //
  } //end if B1

  if ( GammaD1pointer == 0 ) {
	GammaD1pointer = new Matrix*[4] ;
	GammaD1pointer[0] = new Matrix(1,3) ;
	GammaD1pointer[1] = new Matrix(1,3) ;
	GammaD1pointer[2] = new Matrix(1,3) ;
	GammaD1pointer[3] = new Matrix(1,3) ;
  } //end if D1

  if ( GammaA2pointer == 0 ) {
	GammaA2pointer = new Matrix*[4] ;
	GammaA2pointer[0] = new Matrix(1,3) ;
	GammaA2pointer[1] = new Matrix(1,3) ;
	GammaA2pointer[2] = new Matrix(1,3) ;
	GammaA2pointer[3] = new Matrix(1,3) ;
  } //end if A2

  if ( GammaC2pointer == 0 ) {
	GammaC2pointer = new Matrix*[4] ;
	GammaC2pointer[0] = new Matrix(1,3) ;
	GammaC2pointer[1] = new Matrix(1,3) ;
	GammaC2pointer[2] = new Matrix(1,3) ;
	GammaC2pointer[3] = new Matrix(1,3) ;
  } //end if C2

  if ( Bhat == 0 ) {
	Bhat = new Matrix*[4] ;
	Bhat[0] = new Matrix(2,3) ;
	Bhat[1] = new Matrix(2,3) ;
	Bhat[2] = new Matrix(2,3) ;
	Bhat[3] = new Matrix(2,3) ;
  } //end if Bhat

  sg[0] = -one_over_root3;
  sg[1] = one_over_root3;
  sg[2] = one_over_root3;
  sg[3] = -one_over_root3;  

  tg[0] = -one_over_root3;
  tg[1] = -one_over_root3;
  tg[2] = one_over_root3;
  tg[3] = one_over_root3;  

  wg[0] = 1.0;
  wg[1] = 1.0;
  wg[2] = 1.0;
  wg[3] = 1.0;

}
//******************************************************************

//destructor 
ShellMITC4::~ShellMITC4( )
{
  int i ;
  for ( i = 0 ;  i < 4; i++ ) {

    delete materialPointers[i] ;
    materialPointers[i] = 0 ; 

    nodePointers[i] = 0 ;

  } //end for i

  //shear matrix pointers
  if ( GammaB1pointer != 0 ) {
	delete GammaB1pointer[0] ;
	delete GammaB1pointer[1] ;
	delete GammaB1pointer[2] ;
	delete GammaB1pointer[3] ;
	delete [] GammaB1pointer ;
	GammaB1pointer = 0 ;
  } //end if B1

  if ( GammaD1pointer != 0 ) {
	delete GammaD1pointer[0] ;
	delete GammaD1pointer[1] ;
	delete GammaD1pointer[2] ;
	delete GammaD1pointer[3] ;
	delete [] GammaD1pointer ;
	GammaD1pointer = 0 ;
  } //end if D1

  if ( GammaA2pointer != 0 ) {
	delete GammaA2pointer[0] ;
	delete GammaA2pointer[1] ;
	delete GammaA2pointer[2] ;
	delete GammaA2pointer[3] ;
	delete [] GammaA2pointer ;
	GammaA2pointer = 0 ;
  } //end if A2

  if ( GammaC2pointer != 0 ) {
	delete GammaC2pointer[0] ;
	delete GammaC2pointer[1] ;
	delete GammaC2pointer[2] ;
	delete GammaC2pointer[3] ;
	delete [] GammaC2pointer ;
	GammaC2pointer = 0 ;
  } //end if C2

  if ( Bhat != 0 ) {
	delete Bhat[0] ;
	delete Bhat[1] ;
	delete Bhat[2] ;
	delete Bhat[3] ;
	delete [] Bhat ;
	Bhat = 0 ;
  } //end if Bhat

  if (load != 0)
    delete load;

  if (Ki != 0)
    delete Ki;
}
//**************************************************************************


//set domain
void  ShellMITC4::setDomain( Domain *theDomain ) 
{  
  int i, j ;
  static Vector eig(3) ;
  static Matrix ddMembrane(3,3) ;

  //node pointers
  for ( i = 0; i < 4; i++ ) {
     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
     
     if (nodePointers[i] == 0) {
       opserr << "ShellMITC4::setDomain - no node " << connectedExternalNodes(i);
       opserr << " exists in the model\n";
     }
  }

  //compute drilling stiffness penalty parameter
  const Matrix &dd = materialPointers[0]->getInitialTangent( ) ;

  //assemble ddMembrane ;
  for ( i = 0; i < 3; i++ ) {
      for ( j = 0; j < 3; j++ )
         ddMembrane(i,j) = dd(i,j) ;
  } //end for i 

  //eigenvalues of ddMembrane
  eig = LovelyEig( ddMembrane ) ;
  
  //set ktt 
  //Ktt = dd(2,2) ;  //shear modulus 
  Ktt = min( eig(2), min( eig(0), eig(1) ) ) ;
  //Ktt = dd(2,2);

  //basis vectors and local coordinates
  computeBasis( ) ;

  this->DomainComponent::setDomain(theDomain);
}


//get the number of external nodes
int  ShellMITC4::getNumExternalNodes( ) const
{
  return 4 ;
} 
 

//return connected external nodes
const ID&  ShellMITC4::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 


Node **
ShellMITC4::getNodePtrs(void) 
{
  return nodePointers;
} 

//return number of dofs
int  ShellMITC4::getNumDOF( ) 
{
  return 24 ;
}


//commit state
int  ShellMITC4::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "ShellMITC4::commitState () - failed in base class";
  }    

  for (int i = 0; i < 4; i++ ) 
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}
 


//revert to last commit 
int  ShellMITC4::revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int  ShellMITC4::revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

//print out element data
void  ShellMITC4::Print( OPS_Stream &s, int flag )
{
  if (flag == -1) {
    int eleTag = this->getTag();
    s << "EL_QUAD4\t" << eleTag << "\t";
    s << eleTag << "\t" << 1; 
    s  << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
    s  << "\t" << connectedExternalNodes(2) << "\t" << connectedExternalNodes(3) << "\t0.00";
    s << endln;
    s << "PROP_2D\t" << eleTag << "\t";
    s << eleTag << "\t" << 1; 
    s  << "\t" << -1 << "\tSHELL\t1.0\0.0";
    s << endln;
  }  else if (flag < -1) {

     int counter = (flag + 1) * -1;
     int eleTag = this->getTag();
     int i,j;
     for ( i = 0; i < 4; i++ ) {
       const Vector &stress = materialPointers[i]->getStressResultant();
       
       s << "STRESS\t" << eleTag << "\t" << counter << "\t" << i << "\tTOP";
       for (j=0; j<6; j++)
	 s << "\t" << stress(j);
       s << endln;
     }

   } else {
    s << endln ;
    s << "MITC4 Bbar Non-Locking Four Node Shell \n" ;
    s << "Element Number: " << this->getTag() << endln ;
    s << "Node 1 : " << connectedExternalNodes(0) << endln ;
    s << "Node 2 : " << connectedExternalNodes(1) << endln ;
    s << "Node 3 : " << connectedExternalNodes(2) << endln ;
    s << "Node 4 : " << connectedExternalNodes(3) << endln ;
    
    s << "Material Information : \n " ;
    materialPointers[0]->Print( s, flag ) ;
    
    s << endln ;
  }
}

//return stiffness matrix 
const Matrix&  ShellMITC4::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    

//return secant matrix 
const Matrix&  ShellMITC4::getInitialStiff( ) 
{
  if (Ki != 0)
    return *Ki;

  static const int ndf = 6 ; //two membrane plus three bending plus one drill

  static const int nstress = 8 ; //three membrane, three moment, two shear

  static const int ngauss = 4 ;

  static const int numnodes = 4 ;

  int i,  j,  k, p, q ;
  int jj, kk ;
  int node ;

  double volume = 0.0 ;

  static double xsj ;  // determinant jacaobian matrix 

  static double dvol[ngauss] ; //volume element

  static double shp[3][numnodes] ;  //shape functions at a gauss point

  static double Shape[3][numnodes][ngauss] ; //all the shape functions

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 

  static Matrix dd(nstress,nstress) ;  //material tangent

  static Matrix J0(2,2) ;  //Jacobian at center
 
  static Matrix J0inv(2,2) ; //inverse of Jacobian at center

  //---------B-matrices------------------------------------

    static Matrix BJ(nstress,ndf) ;      // B matrix node J

    static Matrix BJtran(ndf,nstress) ;

    static Matrix BK(nstress,ndf) ;      // B matrix node k

    static Matrix BJtranD(ndf,nstress) ;


    static Matrix Bbend(3,3) ;  // bending B matrix

    static Matrix Bshear(2,3) ; // shear B matrix

    static Matrix Bmembrane(3,2) ; // membrane B matrix


    static double BdrillJ[ndf] ; //drill B matrix

    static double BdrillK[ndf] ;  

    double *drillPointer ;

    static double saveB[nstress][ndf][numnodes] ;

  //-------------------------------------------------------

  stiff.Zero( ) ;


  //compute basis vectors and local nodal coordinates
  //computeBasis( ) ;

  //compute Jacobian and inverse at center
  double L1 = 0.0 ;
  double L2 = 0.0 ;
  computeJacobian( L1, L2, xl, J0, J0inv ) ; 

  //compute the gamma's
  computeGamma( xl, J0 ) ;

  //zero Bhat = \frac{1}{volume} \int{  B - \bar{B} } \diff A
  for ( node = 0;  node < numnodes;  node++ ) {
    Bhat[node]->Zero( ) ;
  }

  //gauss loop to compute Bhat's 
  for ( i = 0; i < ngauss; i++ ) {

    //get shape functions    
    shape2d( sg[i], tg[i], xl, shp, xsj ) ;

    //save shape functions
    for ( p = 0; p < 3; p++ ) {
      for ( q = 0; q < numnodes; q++ )
	  Shape[p][q][i] = shp[p][q] ;
    } // end for p

    //volume element to also be saved
    dvol[i] = wg[i] * xsj ;  

    volume += dvol[i] ;

    for ( node = 0; node < numnodes; node++ ) {

      //compute B shear matrix for this node
      //Bhat[node] += (  dvol[i] * computeBshear(node,shp)  ) ;
      Bhat[node]->addMatrix(1.0, 
			   computeBshear(node,shp),
			   dvol[i] ) ;

      //compute B-bar shear matrix for this node
      //Bhat[node] -= ( dvol[i] *
      //            computeBbarShear( node, sg[i], tg[i], J0inv ) 
      //                ) ;
      Bhat[node]->addMatrix(1.0, 
			   computeBbarShear(node,sg[i],tg[i],J0inv),
			   -dvol[i] ) ;

    } //end for node   

  } // end for i gauss loop

  //compute Bhat 
    for ( node = 0;  node < numnodes;  node++ )
      //(*Bhat[node]) /= volume ;
      Bhat[node]->operator/=(volume) ;
    

  //gauss loop 
  for ( i = 0; i < ngauss; i++ ) {

    //extract shape functions from saved array
    for ( p = 0; p < 3; p++ ) {
       for ( q = 0; q < numnodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
    } // end for p


    // j-node loop to compute strain 
    for ( j = 0; j < numnodes; j++ )  {

      //compute B matrix 

      Bmembrane = computeBmembrane( j, shp ) ;

      Bbend = computeBbend( j, shp ) ;

      Bshear = computeBbarShear( j, sg[i], tg[i], J0inv ) ;

      //add Bhat to shear terms 
      Bshear += (*Bhat[j]) ;

      BJ = assembleB( Bmembrane, Bbend, Bshear ) ;

      //save the B-matrix
      for (p=0; p<nstress; p++) {
	for (q=0; q<ndf; q++ )
	  saveB[p][q][j] = BJ(p,q) ;
      }//end for p

      //drilling B matrix
      drillPointer = computeBdrill( j, shp ) ;
      for (p=0; p<ndf; p++ ) {
	//BdrillJ[p] = *drillPointer++ ;
	BdrillJ[p] = *drillPointer ; //set p-th component
	drillPointer++ ;             //pointer arithmetic
      }//end for p
    } // end for j
  

    dd = materialPointers[i]->getInitialTangent( ) ;
    dd *= dvol[i] ;

    //residual and tangent calculations node loops

    jj = 0 ;
    for ( j = 0; j < numnodes; j++ ) {

      //Bmembrane = computeBmembrane( j, shp ) ;

      //Bbend = computeBbend( j, shp ) ;

      //multiply bending terms by (-1.0) for correct statement
      // of equilibrium  
      //Bbend *= (-1.0) ;

      //Bshear = computeBbarShear( j, sg[i], tg[i], J0inv ) ;

      //add Bhat to shear terms
      //Bshear += (*Bhat[j]) ;

      //BJ = assembleB( Bmembrane, Bbend, Bshear ) ;

      //extract BJ
      for (p=0; p<nstress; p++) {
	for (q=0; q<ndf; q++ )
	  BJ(p,q) = saveB[p][q][j]   ;
      }//end for p

      //multiply bending terms by (-1.0) for correct statement
      // of equilibrium  
      for ( p = 3; p < 6; p++ ) {
	for ( q = 3; q < 6; q++ ) 
	  BJ(p,q) *= (-1.0) ;
      } //end for p


      //transpose 
      //BJtran = transpose( 8, ndf, BJ ) ;
      for (p=0; p<ndf; p++) {
	for (q=0; q<nstress; q++) 
	  BJtran(p,q) = BJ(q,p) ;
      }//end for p

      //drilling B matrix
      drillPointer = computeBdrill( j, shp ) ;
      for (p=0; p<ndf; p++ ) {
	BdrillJ[p] = *drillPointer ;
	drillPointer++ ;
      }//end for p

      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixProduct(0.0, BJtran,dd,1.0 ) ;
      
      for (p=0; p<ndf; p++) 
	BdrillJ[p] *= ( Ktt*dvol[i] ) ;
      
      kk = 0 ;
      for ( k = 0; k < numnodes; k++ ) {

	//Bmembrane = computeBmembrane( k, shp ) ;
	
	//Bbend = computeBbend( k, shp ) ;
	
	//Bshear = computeBbarShear( k, sg[i], tg[i], J0inv ) ;

	//Bshear += (*Bhat[k]) ;
	
	//BK = assembleB( Bmembrane, Bbend, Bshear ) ;
	
	//extract BK
	for (p=0; p<nstress; p++) {
	  for (q=0; q<ndf; q++ )
	    BK(p,q) = saveB[p][q][k]   ;
	}//end for p
	
	
	//drilling B matrix
	drillPointer = computeBdrill( k, shp ) ;
	for (p=0; p<ndf; p++ ) {
	  BdrillK[p] = *drillPointer ;
	  drillPointer++ ;
	}//end for p
	
	//stiffJK = BJtranD * BK  ;
	// +  transpose( 1,ndf,BdrillJ ) * BdrillK ; 
	stiffJK.addMatrixProduct(0.0, BJtranD,BK,1.0 ) ;
	
	for ( p = 0; p < ndf; p++ )  {
	  for ( q = 0; q < ndf; q++ ) {
	    stiff( jj+p, kk+q ) += stiffJK(p,q) 
	      + ( BdrillJ[p]*BdrillK[q] ) ;
	  }//end for q
	}//end for p
	
	kk += ndf ;
      } // end for k loop
      
      jj += ndf ;
    } // end for j loop
    
  } //end for i gauss loop 

  Ki = new Matrix(stiff);
  
  return stiff ;
}
    

//return mass matrix
const Matrix&  ShellMITC4::getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 


void  ShellMITC4::zeroLoad( )
{
  if (load != 0)
    load->Zero();

  return ;
}


int 
ShellMITC4::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "ShellMITC4::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}



int 
ShellMITC4::addInertiaLoadToUnbalance(const Vector &accel)
{
  int tangFlag = 1 ;

  int i;

  int allRhoZero = 0;
  for (i=0; i<4; i++) {
    if (materialPointers[i]->getRho() != 0.0)
      allRhoZero = 1;
  }

  if (allRhoZero == 0) 
    return 0;

  int count = 0;
  for (i=0; i<4; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int j=0; j<6; j++)
      resid(count++) = Raccel(i);
  }

  formInertiaTerms( tangFlag ) ;
  if (load == 0) 
    load = new Vector(24);
  load->addMatrixVector(1.0, mass, resid, -1.0);

  return 0;
}



//get residual
const Vector&  ShellMITC4::getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  // subtract external loads 
  if (load != 0)
    resid -= *load;

  return resid ;   
}


//get residual with inertia terms
const Vector&  ShellMITC4::getResistingForceIncInertia( )
{
  static Vector res(24);
  int tang_flag = 0 ; //don't get the tangent

  //do tangent and residual here 
  formResidAndTangent( tang_flag ) ;

  formInertiaTerms( tang_flag ) ;

  res = resid;
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();

  // subtract external loads 
  if (load != 0)
    res -= *load;

  return res;
}

//*********************************************************************
//form inertia terms

void   
ShellMITC4::formInertiaTerms( int tangFlag ) 
{

  //translational mass only
  //rotational inertia terms are neglected


  static const int ndf = 6 ; 

  static const int numberNodes = 4 ;

  static const int numberGauss = 4 ;

  static const int nShape = 3 ;

  static const int massIndex = nShape - 1 ;

  double xsj ;  // determinant jacaobian matrix 

  double dvol ; //volume element

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static Vector momentum(ndf) ;


  int i, j, k, p;
  int jj, kk ;

  double temp, rhoH, massJK ;


  //zero mass 
  mass.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  //computeBasis( ) ;


  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) {

    //get shape functions    
    shape2d( sg[i], tg[i], xl, shp, xsj ) ;

    //volume element to also be saved
    dvol = wg[i] * xsj ;  


    //node loop to compute accelerations
    momentum.Zero( ) ;
    for ( j = 0; j < numberNodes; j++ ) 
      //momentum += ( shp[massIndex][j] * nodePointers[j]->getTrialAccel() ) ;
      momentum.addVector(1.0,  
			 nodePointers[j]->getTrialAccel(),
			 shp[massIndex][j] ) ;

      
    //density
    rhoH = materialPointers[i]->getRho() ;

    //multiply acceleration by density to form momentum
    momentum *= rhoH ;


    //residual and tangent calculations node loops
    //jj = 0 ;
    for ( j=0, jj=0; j<numberNodes; j++, jj+=ndf ) {

      temp = shp[massIndex][j] * dvol ;

      for ( p = 0; p < 3; p++ )
        resid( jj+p ) += ( temp * momentum(p) ) ;

      
      if ( tangFlag == 1 && rhoH != 0.0) {

	 //multiply by density
	 temp *= rhoH ;

	 //node-node translational mass
         //kk = 0 ;
         for ( k=0, kk=0; k<numberNodes; k++, kk+=ndf ) {

	   massJK = temp * shp[massIndex][k] ;

	   for ( p = 0; p < 3; p++ ) 
	      mass( jj+p, kk+p ) +=  massJK ;
            
	   //kk += ndf ;
          } // end for k loop

      } // end if tang_flag 

      //jj += ndf ;
    } // end for j loop


  } //end for i gauss loop 

}

//*********************************************************************

//form residual and tangent
void  
ShellMITC4::formResidAndTangent( int tang_flag ) 
{
  //
  //  six(6) nodal dof's ordered :
  //
  //    -        - 
  //   |    u1    |   <---plate membrane
  //   |    u2    |
  //   |----------|
  //   |  w = u3  |   <---plate bending
  //   |  theta1  | 
  //   |  theta2  | 
  //   |----------|
  //   |  theta3  |   <---drill 
  //    -        -  
  //
  // membrane strains ordered :
  //
  //            strain(0) =   eps00     i.e.   (11)-strain
  //            strain(1) =   eps11     i.e.   (22)-strain
  //            strain(2) =   gamma01   i.e.   (12)-shear
  //
  // curvatures and shear strains ordered  :
  //
  //            strain(3) =     kappa00  i.e.   (11)-curvature
  //            strain(4) =     kappa11  i.e.   (22)-curvature
  //            strain(5) =   2*kappa01  i.e. 2*(12)-curvature 
  //
  //            strain(6) =     gamma02  i.e.   (13)-shear
  //            strain(7) =     gamma12  i.e.   (23)-shear
  //
  //  same ordering for moments/shears but no 2 
  //  
  //  Then, 
  //              epsilon00 = -z * kappa00      +    eps00_membrane
  //              epsilon11 = -z * kappa11      +    eps11_membrane
  //  gamma01 = 2*epsilon01 = -z * (2*kappa01)  +  gamma01_membrane 
  //
  //  Shear strains gamma02, gamma12 constant through cross section
  //

  static const int ndf = 6 ; //two membrane plus three bending plus one drill

  static const int nstress = 8 ; //three membrane, three moment, two shear

  static const int ngauss = 4 ;

  static const int numnodes = 4 ;

  int i,  j,  k, p, q ;
  int jj, kk ;
  int node ;

  int success ;
  
  double volume = 0.0 ;

  static double xsj ;  // determinant jacaobian matrix 

  static double dvol[ngauss] ; //volume element

  static Vector strain(nstress) ;  //strain

  static double shp[3][numnodes] ;  //shape functions at a gauss point

  static double Shape[3][numnodes][ngauss] ; //all the shape functions

  static Vector residJ(ndf) ; //nodeJ residual 

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 

  static Vector stress(nstress) ;  //stress resultants

  static Matrix dd(nstress,nstress) ;  //material tangent

  static Matrix J0(2,2) ;  //Jacobian at center
 
  static Matrix J0inv(2,2) ; //inverse of Jacobian at center

  double epsDrill = 0.0 ;  //drilling "strain"

  double tauDrill = 0.0 ; //drilling "stress"

  //---------B-matrices------------------------------------

    static Matrix BJ(nstress,ndf) ;      // B matrix node J

    static Matrix BJtran(ndf,nstress) ;

    static Matrix BK(nstress,ndf) ;      // B matrix node k

    static Matrix BJtranD(ndf,nstress) ;


    static Matrix Bbend(3,3) ;  // bending B matrix

    static Matrix Bshear(2,3) ; // shear B matrix

    static Matrix Bmembrane(3,2) ; // membrane B matrix


    static double BdrillJ[ndf] ; //drill B matrix

    static double BdrillK[ndf] ;  

    double *drillPointer ;

    static double saveB[nstress][ndf][numnodes] ;

  //-------------------------------------------------------

    

  //zero stiffness and residual 
  stiff.Zero( ) ;
  resid.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  //computeBasis( ) ;

  //compute Jacobian and inverse at center
  double L1 = 0.0 ;
  double L2 = 0.0 ;
  computeJacobian( L1, L2, xl, J0, J0inv ) ; 

  //compute the gamma's
  computeGamma( xl, J0 ) ;

  //zero Bhat = \frac{1}{volume} \int{  B - \bar{B} } \diff A
  for ( node = 0;  node < numnodes;  node++ ) {
    Bhat[node]->Zero( ) ;
  }

  //gauss loop to compute Bhat's 
  for ( i = 0; i < ngauss; i++ ) {

    //get shape functions    
    shape2d( sg[i], tg[i], xl, shp, xsj ) ;

    //save shape functions
    for ( p = 0; p < 3; p++ ) {
      for ( q = 0; q < numnodes; q++ )
	  Shape[p][q][i] = shp[p][q] ;
    } // end for p

    //volume element to also be saved
    dvol[i] = wg[i] * xsj ;  

    volume += dvol[i] ;

    for ( node = 0; node < numnodes; node++ ) {

      //compute B shear matrix for this node
      //Bhat[node] += (  dvol[i] * computeBshear(node,shp)  ) ;
      Bhat[node]->addMatrix(1.0, 
			   computeBshear(node,shp),
			   dvol[i] ) ;

      //compute B-bar shear matrix for this node
      //Bhat[node] -= ( dvol[i] *
      //            computeBbarShear( node, sg[i], tg[i], J0inv ) 
      //                ) ;
      Bhat[node]->addMatrix(1.0, 
			   computeBbarShear(node,sg[i],tg[i],J0inv),
			   -dvol[i] ) ;

    } //end for node   

  } // end for i gauss loop

  //compute Bhat 
    for ( node = 0;  node < numnodes;  node++ )
      //(*Bhat[node]) /= volume ;
      Bhat[node]->operator/=(volume) ;
    

  //gauss loop 
  for ( i = 0; i < ngauss; i++ ) {

    //extract shape functions from saved array
    for ( p = 0; p < 3; p++ ) {
       for ( q = 0; q < numnodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
    } // end for p


    //zero the strains
    strain.Zero( ) ;
    epsDrill = 0.0 ;


    // j-node loop to compute strain 
    for ( j = 0; j < numnodes; j++ )  {

      //compute B matrix 

      Bmembrane = computeBmembrane( j, shp ) ;

      Bbend = computeBbend( j, shp ) ;

      Bshear = computeBbarShear( j, sg[i], tg[i], J0inv ) ;

      //add Bhat to shear terms 
      Bshear += (*Bhat[j]) ;

      BJ = assembleB( Bmembrane, Bbend, Bshear ) ;

      //save the B-matrix
      for (p=0; p<nstress; p++) {
	for (q=0; q<ndf; q++ )
	  saveB[p][q][j] = BJ(p,q) ;
      }//end for p


      //nodal "displacements" 
      const Vector &ul = nodePointers[j]->getTrialDisp( ) ;

      //compute the strain
      //strain += (BJ*ul) ; 
      strain.addMatrixVector(1.0, BJ,ul,1.0 ) ;

      //drilling B matrix
      drillPointer = computeBdrill( j, shp ) ;
      for (p=0; p<ndf; p++ ) {
	//BdrillJ[p] = *drillPointer++ ;
	BdrillJ[p] = *drillPointer ; //set p-th component
	drillPointer++ ;             //pointer arithmetic
      }//end for p

      //drilling "strain" 
      for ( p = 0; p < ndf; p++ )
	epsDrill +=  BdrillJ[p]*ul(p) ;

    } // end for j
  

    //send the strain to the material 
    success = materialPointers[i]->setTrialSectionDeformation( strain ) ;

    //compute the stress
    stress = materialPointers[i]->getStressResultant( ) ;

    //drilling "stress" 
    tauDrill = Ktt * epsDrill ;

    //multiply by volume element
    stress   *= dvol[i] ;
    tauDrill *= dvol[i] ;

    if ( tang_flag == 1 ) {
      dd = materialPointers[i]->getSectionTangent( ) ;
      dd *= dvol[i] ;
    } //end if tang_flag


    //residual and tangent calculations node loops

    jj = 0 ;
    for ( j = 0; j < numnodes; j++ ) {

      //Bmembrane = computeBmembrane( j, shp ) ;

      //Bbend = computeBbend( j, shp ) ;

      //multiply bending terms by (-1.0) for correct statement
      // of equilibrium  
      //Bbend *= (-1.0) ;

      //Bshear = computeBbarShear( j, sg[i], tg[i], J0inv ) ;

      //add Bhat to shear terms
      //Bshear += (*Bhat[j]) ;

      //BJ = assembleB( Bmembrane, Bbend, Bshear ) ;

      //extract BJ
      for (p=0; p<nstress; p++) {
	for (q=0; q<ndf; q++ )
	  BJ(p,q) = saveB[p][q][j]   ;
      }//end for p

      //multiply bending terms by (-1.0) for correct statement
      // of equilibrium  
      for ( p = 3; p < 6; p++ ) {
	for ( q = 3; q < 6; q++ ) 
	  BJ(p,q) *= (-1.0) ;
      } //end for p


      //transpose 
      //BJtran = transpose( 8, ndf, BJ ) ;
      for (p=0; p<ndf; p++) {
	for (q=0; q<nstress; q++) 
	  BJtran(p,q) = BJ(q,p) ;
      }//end for p


      //residJ = BJtran * stress ;
      residJ.addMatrixVector(0.0, BJtran,stress,1.0 ) ;

      //drilling B matrix
      drillPointer = computeBdrill( j, shp ) ;
      for (p=0; p<ndf; p++ ) {
	BdrillJ[p] = *drillPointer ;
	drillPointer++ ;
      }//end for p

      //residual including drill
      for ( p = 0; p < ndf; p++ )
        resid( jj + p ) += ( residJ(p) + BdrillJ[p]*tauDrill ) ;


      if ( tang_flag == 1 ) {

        //BJtranD = BJtran * dd ;
	BJtranD.addMatrixProduct(0.0, BJtran,dd,1.0 ) ;

	for (p=0; p<ndf; p++) 
	  BdrillJ[p] *= ( Ktt*dvol[i] ) ;

        kk = 0 ;
        for ( k = 0; k < numnodes; k++ ) {

          //Bmembrane = computeBmembrane( k, shp ) ;

          //Bbend = computeBbend( k, shp ) ;

          //Bshear = computeBbarShear( k, sg[i], tg[i], J0inv ) ;

          //Bshear += (*Bhat[k]) ;

          //BK = assembleB( Bmembrane, Bbend, Bshear ) ;

	  //extract BK
	  for (p=0; p<nstress; p++) {
	    for (q=0; q<ndf; q++ )
	      BK(p,q) = saveB[p][q][k]   ;
	  }//end for p

	  
	  //drilling B matrix
	  drillPointer = computeBdrill( k, shp ) ;
	  for (p=0; p<ndf; p++ ) {
	    BdrillK[p] = *drillPointer ;
	    drillPointer++ ;
	  }//end for p
  
 
          //stiffJK = BJtranD * BK  ;
	  // +  transpose( 1,ndf,BdrillJ ) * BdrillK ; 
	  stiffJK.addMatrixProduct(0.0, BJtranD,BK,1.0 ) ;

          for ( p = 0; p < ndf; p++ )  {
	    for ( q = 0; q < ndf; q++ ) {
	       stiff( jj+p, kk+q ) += stiffJK(p,q) 
		                   + ( BdrillJ[p]*BdrillK[q] ) ;
	    }//end for q
          }//end for p

          kk += ndf ;
        } // end for k loop

      } // end if tang_flag 

      jj += ndf ;
    } // end for j loop


  } //end for i gauss loop 

  
  return ;
}


//************************************************************************
//compute local coordinates and basis

void   
ShellMITC4::computeBasis( ) 
{
  //could compute derivatives \frac{ \partial {\bf x} }{ \partial L_1 } 
  //                     and  \frac{ \partial {\bf x} }{ \partial L_2 }
  //and use those as basis vectors but this is easier 
  //and the shell is flat anyway.

  static Vector temp(3) ;

  static Vector v1(3) ;
  static Vector v2(3) ;
  static Vector v3(3) ;

  //get two vectors (v1, v2) in plane of shell by 
  // nodal coordinate differences

  const Vector &coor0 = nodePointers[0]->getCrds( ) ;

  const Vector &coor1 = nodePointers[1]->getCrds( ) ;

  const Vector &coor2 = nodePointers[2]->getCrds( ) ;
  
  const Vector &coor3 = nodePointers[3]->getCrds( ) ;

  v1.Zero( ) ;
  //v1 = 0.5 * ( coor2 + coor1 - coor3 - coor0 ) ;
  v1  = coor2 ;
  v1 += coor1 ;
  v1 -= coor3 ;
  v1 -= coor0 ;
  v1 *= 0.50 ;
  
  v2.Zero( ) ;
  //v2 = 0.5 * ( coor3 + coor2 - coor1 - coor0 ) ;
  v2  = coor3 ;
  v2 += coor2 ;
  v2 -= coor1 ;
  v2 -= coor0 ;
  v2 *= 0.50 ;
    

 
  //normalize v1 
  //double length = LovelyNorm( v1 ) ;
  double length = v1.Norm( ) ;
  v1 /= length ;


  //Gram-Schmidt process for v2 

  //double alpha = LovelyInnerProduct( v2, v1 ) ;
  double alpha = v2^v1 ;

  //v2 -= alpha*v1 ;
  temp = v1 ;
  temp *= alpha ;
  v2 -= temp ;


  //normalize v2 
  //length = LovelyNorm( v2 ) ;
  length = v2.Norm( ) ;
  v2 /= length ;


  //cross product for v3  
  v3 = LovelyCrossProduct( v1, v2 ) ;

  
  //local nodal coordinates in plane of shell

  int i ;
  for ( i = 0; i < 4; i++ ) {

       const Vector &coorI = nodePointers[i]->getCrds( ) ;

       xl[0][i] = coorI^v1 ;  
       xl[1][i] = coorI^v2 ;

  }  //end for i 

  //basis vectors stored as array of doubles
  for ( i = 0; i < 3; i++ ) {
      g1[i] = v1(i) ;
      g2[i] = v2(i) ;
      g3[i] = v3(i) ;
  }  //end for i 

}

//*************************************************************************
//compute Bdrill

double*
ShellMITC4::computeBdrill( int node, const double shp[3][4] )
{

  //static Matrix Bdrill(1,6) ;
  static double Bdrill[6] ;

  static double B1 ;
  static double B2 ;
  static double B6 ;


//---Bdrill Matrix in standard {1,2,3} mechanics notation---------
//
//             -                                       -
//   Bdrill = | -0.5*N,2   +0.5*N,1    0    0    0   -N |   (1x6) 
//             -                                       -  
//
//----------------------------------------------------------------

  //  Bdrill.Zero( ) ;

  //Bdrill(0,0) = -0.5*shp[1][node] ;

  //Bdrill(0,1) = +0.5*shp[0][node] ;

  //Bdrill(0,5) =     -shp[2][node] ;


  B1 =  -0.5*shp[1][node] ; 

  B2 =  +0.5*shp[0][node] ;

  B6 =  -shp[2][node] ;

  Bdrill[0] = B1*g1[0] + B2*g2[0] ;
  Bdrill[1] = B1*g1[1] + B2*g2[1] ; 
  Bdrill[2] = B1*g1[2] + B2*g2[2] ;

  Bdrill[3] = B6*g3[0] ;
  Bdrill[4] = B6*g3[1] ; 
  Bdrill[5] = B6*g3[2] ;
 
  return Bdrill ;

}

//*************************************************************************
//compute the gamma's

void   
ShellMITC4::computeGamma( const double xl[2][4], const Matrix &J )
{
  //static Matrix e1(1,2) ;
  //static Matrix e2(1,2) ;

  static double Jtran[2][2] ;  // J transpose

  static Matrix e1Jtran(1,2) ; // e1*Jtran 

  static Matrix e2Jtran(1,2) ; // e2*Jtran

  static Matrix Bshear(2,3) ; // shear strain-displacement matrix

  static double shape[3][4] ; // shape functions 

  double xsj ; 

  double L1, L2 ;

  int node ;

 /*  
  e1(0,0) = 1.0 ;
  e1(0,1) = 0.0 ;

  e2(0,0) = 0.0 ;
  e2(0,1) = 1.0 ;
 */

  //Jtran = transpose( 2, 2, J ) ;
  Jtran[0][0] = J(0,0) ;
  Jtran[1][1] = J(1,1) ;
  Jtran[0][1] = J(1,0) ;
  Jtran[1][0] = J(0,1) ;

  //first row of Jtran
  //e1Jtran = e1*Jtran ; 
  e1Jtran(0,0) = Jtran[0][0] ;
  e1Jtran(0,1) = Jtran[0][1] ;

  //second row of Jtran 
  //e2Jtran = e2*Jtran ;
  e2Jtran(0,0) = Jtran[1][0] ;
  e2Jtran(0,1) = Jtran[1][1] ;


  //-------------GammaB1--------------------------------

    //set natural coordinate point location
    L1 =  0.0 ;
    L2 = -1.0 ;

    //get shape functions    
    shape2d( L1, L2, xl, shape, xsj ) ;

    for ( node = 0; node < 4; node++ ) {

        //compute shear B        
        Bshear = computeBshear( node, shape ) ;

        GammaB1pointer[node]->Zero( ) ;

	//( *GammaB1pointer[node] ) = e1Jtran*Bshear ;
	GammaB1pointer[node]->addMatrixProduct(0.0, e1Jtran,Bshear,1.0 );

    } //end for node
 
  //----------------------------------------------------


  //-------------GammaD1--------------------------------

    //set natural coordinate point location
    L1 =  0.0 ;
    L2 = +1.0 ;

    //get shape functions    
    shape2d( L1, L2, xl, shape, xsj ) ;

    for ( node = 0; node < 4; node++ ) {

        //compute shear B        
        Bshear = computeBshear( node, shape ) ;

        GammaD1pointer[node]->Zero( ) ;

        //( *GammaD1pointer[node] ) = e1Jtran*Bshear ;
	GammaD1pointer[node]->addMatrixProduct(0.0, e1Jtran,Bshear,1.0 );

    } //end for node
 
  //----------------------------------------------------    


  //-------------GammaA2--------------------------------

    //set natural coordinate point location
    L1 = -1.0 ;
    L2 =  0.0 ;

    //get shape functions    
    shape2d( L1, L2, xl, shape, xsj ) ;

    for ( node = 0; node < 4; node++ ) {

        //compute shear B        
        Bshear = computeBshear( node, shape ) ;

        GammaA2pointer[node]->Zero( ) ;

        //( *GammaA2pointer[node] ) = e2Jtran*Bshear ;
	GammaA2pointer[node]->addMatrixProduct(0.0, e2Jtran,Bshear,1.0 );
		   

    } //end for node
 
  //----------------------------------------------------     


  //-------------GammaC2--------------------------------

    //set natural coordinate point location
    L1 = +1.0 ;
    L2 =  0.0 ;

    //get shape functions    
    shape2d( L1, L2, xl, shape, xsj ) ;

    for ( node = 0; node < 4; node++ ) {

        //compute shear B        
        Bshear = computeBshear( node, shape ) ;

        GammaC2pointer[node]->Zero( ) ;

        //( *GammaC2pointer[node] ) = e2Jtran*Bshear ;
	GammaC2pointer[node]->addMatrixProduct(0.0, e2Jtran,Bshear,1.0 );

    } //end for node
 
  //----------------------------------------------------     

 return ;

}

//********************************************************************
//assemble a B matrix

const Matrix&  
ShellMITC4::assembleB( const Matrix &Bmembrane,
                               const Matrix &Bbend, 
                               const Matrix &Bshear ) 
{

  //Matrix Bbend(3,3) ;  // plate bending B matrix

  //Matrix Bshear(2,3) ; // plate shear B matrix

  //Matrix Bmembrane(3,2) ; // plate membrane B matrix


    static Matrix B(8,6) ;

    static Matrix BmembraneShell(3,3) ; 
    
    static Matrix BbendShell(3,3) ; 

    static Matrix BshearShell(2,6) ;
 
    static Matrix Gmem(2,3) ;

    static Matrix Gshear(3,6) ;

    int p, q ;
    int pp ;

//    
// For Shell : 
//
//---B Matrices in standard {1,2,3} mechanics notation---------
//
//            -                     _          
//           | Bmembrane  |     0    |
//           | --------------------- |     
//    B =    |     0      |  Bbend   |   (8x6) 
//           | --------------------- |
//           |         Bshear        |
//            -           -         -
//
//-------------------------------------------------------------
//
//


    //shell modified membrane terms
    
    Gmem(0,0) = g1[0] ;
    Gmem(0,1) = g1[1] ;
    Gmem(0,2) = g1[2] ;

    Gmem(1,0) = g2[0] ;
    Gmem(1,1) = g2[1] ;
    Gmem(1,2) = g2[2] ;

    //BmembraneShell = Bmembrane * Gmem ;
    BmembraneShell.addMatrixProduct(0.0, Bmembrane,Gmem,1.0 ) ;


    //shell modified bending terms 

    Matrix &Gbend = Gmem ;

    //BbendShell = Bbend * Gbend ;
    BbendShell.addMatrixProduct(0.0, Bbend,Gbend,1.0 ) ; 


    //shell modified shear terms 

    Gshear.Zero( ) ;

    Gshear(0,0) = g3[0] ;
    Gshear(0,1) = g3[1] ;
    Gshear(0,2) = g3[2] ;

    Gshear(1,3) = g1[0] ;
    Gshear(1,4) = g1[1] ;
    Gshear(1,5) = g1[2] ;

    Gshear(2,3) = g2[0] ;
    Gshear(2,4) = g2[1] ;
    Gshear(2,5) = g2[2] ;

    //BshearShell = Bshear * Gshear ;
    BshearShell.addMatrixProduct(0.0, Bshear,Gshear,1.0 ) ;
   

  B.Zero( ) ;

  //assemble B from sub-matrices 

  //membrane terms 
  for ( p = 0; p < 3; p++ ) {

      for ( q = 0; q < 3; q++ ) 
          B(p,q) = BmembraneShell(p,q) ;

  } //end for p


  //bending terms
  for ( p = 3; p < 6; p++ ) {
    pp = p - 3 ;
    for ( q = 3; q < 6; q++ ) 
        B(p,q) = BbendShell(pp,q-3) ; 
  } // end for p
    

  //shear terms 
  for ( p = 0; p < 2; p++ ) {
      pp = p + 6 ;
      
      for ( q = 0; q < 6; q++ ) {
          B(pp,q) = BshearShell(p,q) ; 
      } // end for q
 
  } //end for p
  
  return B ;

}

//***********************************************************************
//compute Bmembrane matrix

const Matrix&   
ShellMITC4::computeBmembrane( int node, const double shp[3][4] ) 
{

  static Matrix Bmembrane(3,2) ;

//---Bmembrane Matrix in standard {1,2,3} mechanics notation---------
//
//                -             -
//               | +N,1      0   | 
// Bmembrane =   |   0     +N,2  |    (3x2)
//               | +N,2    +N,1  |
//                -             -  
//
//  three(3) strains and two(2) displacements (for plate)
//-------------------------------------------------------------------

  Bmembrane.Zero( ) ;

  Bmembrane(0,0) = shp[0][node] ;
  Bmembrane(1,1) = shp[1][node] ;
  Bmembrane(2,0) = shp[1][node] ;
  Bmembrane(2,1) = shp[0][node] ;

  return Bmembrane ;

}

//***********************************************************************
//compute Bbend matrix

const Matrix&   
ShellMITC4::computeBbend( int node, const double shp[3][4] )
{

    static Matrix Bbend(3,2) ;

//---Bbend Matrix in standard {1,2,3} mechanics notation---------
//
//            -             -
//   Bbend = |    0    -N,1  | 
//           |  +N,2     0   |    (3x2)
//           |  +N,1   -N,2  |
//            -             -  
//
//  three(3) curvatures and two(2) rotations (for plate)
//----------------------------------------------------------------

    Bbend.Zero( ) ;

    Bbend(0,1) = -shp[0][node] ;
    Bbend(1,0) =  shp[1][node] ;
    Bbend(2,0) =  shp[0][node] ;
    Bbend(2,1) = -shp[1][node] ; 

    return Bbend ;
}

//***********************************************************************
//compute Bbar shear matrix

const Matrix&  
ShellMITC4::computeBbarShear( int node, double L1, double L2,
			      const Matrix &Jinv ) 
{
    static Matrix Bshear(2,3) ;
    static Matrix BshearNat(2,3) ;

    static Matrix JinvTran(2,2) ;  // J-inverse-transpose

    static Matrix Gamma1(1,3) ;
    static Matrix Gamma2(1,3) ; 

    static Matrix temp1(1,3) ;
    static Matrix temp2(1,3) ;


    //JinvTran = transpose( 2, 2, Jinv ) ;
    JinvTran(0,0) = Jinv(0,0) ;
    JinvTran(1,1) = Jinv(1,1) ;
    JinvTran(0,1) = Jinv(1,0) ;
    JinvTran(1,0) = Jinv(0,1) ;


    // compute BShear from Bbar interpolation
    
    //Gamma1 = ( 1.0 - L2 )*( *GammaB1pointer[node] )  +
    //         ( 1.0 + L2 )*( *GammaD1pointer[node] )   ;
    temp1  = *GammaB1pointer[node] ;
    temp1 *= ( 1.0 - L2 ) ;
    
    temp2  = *GammaD1pointer[node] ;
    temp2 *= ( 1.0 + L2 ) ;

    Gamma1  = temp1 ;
    Gamma1 += temp2 ;


    //Gamma2 = ( 1.0 - L1 )*( *GammaA2pointer[node] )  +
    //         ( 1.0 + L1 )*( *GammaC2pointer[node] )   ;
    temp1  = *GammaA2pointer[node] ;
    temp1 *= ( 1.0 - L1 ) ;
    
    temp2  = *GammaC2pointer[node] ;
    temp2 *= ( 1.0 + L1 ) ;

    Gamma2  = temp1 ;
    Gamma2 += temp2 ;

    //don't forget the 1/2
    //Gamma1 *= 0.50 ;
    //Gamma2 *= 0.50 ;


    BshearNat(0,0) = Gamma1(0,0) ;
    BshearNat(0,1) = Gamma1(0,1) ;
    BshearNat(0,2) = Gamma1(0,2) ;

    BshearNat(1,0) = Gamma2(0,0) ;
    BshearNat(1,1) = Gamma2(0,1) ;
    BshearNat(1,2) = Gamma2(0,2) ;


    //strain tensor push on BshearNat
    //Bshear = ( JinvTran * BshearNat ) ;
    Bshear.addMatrixProduct(0.0,  JinvTran,BshearNat,0.5 ) ;
    //note the 1/2 -----------------------------------^

    return Bshear ;
}

//***********************************************************************
//compute standard Bshear matrix

const Matrix&  
ShellMITC4::computeBshear( int node, const double shp[3][4] )
{

  static Matrix Bshear(2,3) ;

//---Bshear Matrix in standard {1,2,3} mechanics notation------
//
//             -                -
//   Bshear = | +N,1      0    +N |  (2x3)
//            | +N,2     -N     0 |
//             -                -  
//
//  two(2) shear strains, one(1) displacement and two(2) rotations for plate 
//-----------------------------------------------------------------------

  Bshear.Zero( ) ;

  Bshear(0,0) =  shp[0][node] ;
  Bshear(0,2) =  shp[2][node] ;
  Bshear(1,0) =  shp[1][node] ;
  Bshear(1,1) = -shp[2][node] ;

  return Bshear ;

}  

//***********************************************************************
//compute Jacobian matrix and inverse at point {L1,L2} 

void   
ShellMITC4::computeJacobian( double L1, double L2,
				    const double x[2][4], 
                                    Matrix &JJ, 
                                    Matrix &JJinv )
{
  int i, j, k ;
     
  static const double s[] = { -0.5,  0.5, 0.5, -0.5 } ;
  static const double t[] = { -0.5, -0.5, 0.5,  0.5 } ;

  static double shp[2][4] ;

  double ss = L1 ;
  double tt = L2 ;

  for ( i = 0; i < 4; i++ ) {
      shp[0][i] = s[i] * ( 0.5 + t[i]*tt ) ;
      shp[1][i] = t[i] * ( 0.5 + s[i]*ss ) ;
  } // end for i

  
  // Construct jacobian and its inverse
  
  JJ.Zero( ) ;
  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++ ) {

      for ( k = 0; k < 4; k++ )
	  JJ(i,j) +=  x[i][k] * shp[j][k] ;

    } //end for j
  }  // end for i 

  double xsj = JJ(0,0)*JJ(1,1) - JJ(0,1)*JJ(1,0) ;

  //inverse jacobian
  double jinv = 1.0 / xsj ;
  JJinv(0,0) =  JJ(1,1) * jinv ;
  JJinv(1,1) =  JJ(0,0) * jinv ;
  JJinv(0,1) = -JJ(0,1) * jinv ; 
  JJinv(1,0) = -JJ(1,0) * jinv ; 

  return ;

}

//************************************************************************
//shape function routine for four node quads

void  
ShellMITC4::shape2d( double ss, double tt, 
		           const double x[2][4], 
		           double shp[3][4], 
		           double &xsj            )

{ 

  int i, j, k ;

  double temp ;
     
  static const double s[] = { -0.5,  0.5, 0.5, -0.5 } ;
  static const double t[] = { -0.5, -0.5, 0.5,  0.5 } ;

  static double xs[2][2] ;
  static double sx[2][2] ;

  for ( i = 0; i < 4; i++ ) {
      shp[2][i] = ( 0.5 + s[i]*ss )*( 0.5 + t[i]*tt ) ;
      shp[0][i] = s[i] * ( 0.5 + t[i]*tt ) ;
      shp[1][i] = t[i] * ( 0.5 + s[i]*ss ) ;
  } // end for i

  
  // Construct jacobian and its inverse
  
  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++ ) {

      xs[i][j] = 0.0 ;

      for ( k = 0; k < 4; k++ )
	  xs[i][j] +=  x[i][k] * shp[j][k] ;

    } //end for j
  }  // end for i 

  xsj = xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0] ;

  //inverse jacobian
  double jinv = 1.0 / xsj ;
  sx[0][0] =  xs[1][1] * jinv ;
  sx[1][1] =  xs[0][0] * jinv ;
  sx[0][1] = -xs[0][1] * jinv ; 
  sx[1][0] = -xs[1][0] * jinv ; 


  //form global derivatives 
  
  for ( i = 0; i < 4; i++ ) {
    temp      = shp[0][i]*sx[0][0] + shp[1][i]*sx[1][0] ;
    shp[1][i] = shp[0][i]*sx[0][1] + shp[1][i]*sx[1][1] ;
    shp[0][i] = temp ;
  } // end for i

  return ;
}
	   
//**********************************************************************

Matrix  
ShellMITC4::transpose( int dim1, 
                                       int dim2, 
		                       const Matrix &M ) 
{
  int i ;
  int j ;

  Matrix Mtran( dim2, dim1 ) ;

  for ( i = 0; i < dim1; i++ ) {
     for ( j = 0; j < dim2; j++ ) 
         Mtran(j,i) = M(i,j) ;
  } // end for i

  return Mtran ;
}

//**********************************************************************

int  ShellMITC4::sendSelf (int commitTag, Channel &theChannel)
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
    opserr << "WARNING ShellMITC4::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  static Vector vectData(1);
  vectData(0) = Ktt;

  res += theChannel.sendVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellMITC4::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING ShellMITC4::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;
}
    
int  ShellMITC4::recvSelf (int commitTag, 
		       Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(13);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellMITC4::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(8));
  connectedExternalNodes(0) = idData(9);
  connectedExternalNodes(1) = idData(10);
  connectedExternalNodes(2) = idData(11);
  connectedExternalNodes(3) = idData(12);

  static Vector vectData(1);
  res += theChannel.recvVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellMITC4::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  Ktt = vectData(0);

  int i;

  if (materialPointers[0] == 0) {
    for (i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewSection(matClassTag);
      if (materialPointers[i] == 0) {
	opserr << "ShellMITC4::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
	return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "NLBeamColumn3d::recvSelf() - material " << i <<
	  "failed to recv itself\n";
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
	materialPointers[i] = theBroker.getNewSection(matClassTag);
	if (materialPointers[i] == 0) {
	  opserr << "ShellMITC4::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;
	  exit(-1);
	}
      }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "ShellMITC4::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }
  
  return res;
}
//**************************************************************************

int
ShellMITC4::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    const Vector &end1Crd = nodePointers[0]->getCrds();
    const Vector &end2Crd = nodePointers[1]->getCrds();	
    const Vector &end3Crd = nodePointers[2]->getCrds();	
    const Vector &end4Crd = nodePointers[3]->getCrds();	

    static Matrix coords(4,3);
    static Vector values(4);
    static Vector P(24) ;

    for (int j=0; j<4; j++)
      values(j) = 0.0;

    if (displayMode >= 0) {
      const Vector &end1Disp = nodePointers[0]->getDisp();
      const Vector &end2Disp = nodePointers[1]->getDisp();
      const Vector &end3Disp = nodePointers[2]->getDisp();
      const Vector &end4Disp = nodePointers[3]->getDisp();
      
      if (displayMode < 8 && displayMode > 0) {
	for (int i=0; i<4; i++) {
	  const Vector &stress = materialPointers[i]->getStressResultant();
	  values(i) = stress(displayMode-1);
	}
      }

      for (int i = 0; i < 3; i++) {
	coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
	coords(1,i) = end2Crd(i) + end2Disp(i)*fact;    
	coords(2,i) = end3Crd(i) + end3Disp(i)*fact;    
	coords(3,i) = end4Crd(i) + end4Disp(i)*fact;    

      }
    } else {
      int mode = displayMode * -1;
      const Matrix &eigen1 = nodePointers[0]->getEigenvectors();
      const Matrix &eigen2 = nodePointers[1]->getEigenvectors();
      const Matrix &eigen3 = nodePointers[2]->getEigenvectors();
      const Matrix &eigen4 = nodePointers[3]->getEigenvectors();
      if (eigen1.noCols() >= mode) {
	for (int i = 0; i < 3; i++) {
	  coords(0,i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	  coords(1,i) = end2Crd(i) + eigen2(i,mode-1)*fact;
	  coords(2,i) = end3Crd(i) + eigen3(i,mode-1)*fact;
	  coords(3,i) = end4Crd(i) + eigen4(i,mode-1)*fact;
	}    
      } else {
	for (int i = 0; i < 3; i++) {
	  coords(0,i) = end1Crd(i);
	  coords(1,i) = end2Crd(i);
	  coords(2,i) = end3Crd(i);
	  coords(3,i) = end4Crd(i);
	}    
      }
    }

    //opserr << coords;
    int error = 0;

    //error += theViewer.drawPolygon (coords, values);
    error += theViewer.drawPolygon (coords, values);

    return error;
}
