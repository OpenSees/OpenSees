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
                                                                        
// $Revision: 1.11 $
// $Date: 2006-08-04 19:07:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/NineNodeMixedQuad.cpp,v $

// Ed "C++" Love
//
// Mixed Presssure/Volume Nine Node Quadrilateral
// Plane Strain (NOT PLANE STRESS)
//

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
#include <NineNodeMixedQuad.h>
#include <ElementResponse.h>

#include <Renderer.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static data
Matrix  NineNodeMixedQuad :: stiff(18,18)   ;
Vector  NineNodeMixedQuad :: resid(18)     ;
Matrix  NineNodeMixedQuad :: mass(18,18)    ;
double  NineNodeMixedQuad::xl[2][9];
 
//quadrature data
double   NineNodeMixedQuad::root06 = sqrt(0.6) ;
double   NineNodeMixedQuad::sg[] = { -root06,   0.0,      root06  } ;
double   NineNodeMixedQuad::wg[] = {  5.0/9.0,  8.0/9.0,  5.0/9.0 } ;
  

//null constructor
NineNodeMixedQuad :: NineNodeMixedQuad( ) :
Element( 0, ELE_TAG_NineNodeMixedQuad ),
connectedExternalNodes(9) , load(0), Ki(0)
{ 
  for (int i=0; i<9; i++)
    materialPointers[i] = 0;
}


//full constructor
NineNodeMixedQuad :: NineNodeMixedQuad( int tag, 
					int node1,
					int node2,
					int node3,
					int node4,
					int node5,
					int node6,
					int node7,
					int node8,
					int node9,
					NDMaterial &theMaterial ) :
Element( tag, ELE_TAG_NineNodeMixedQuad ),
connectedExternalNodes(9) , load(0), Ki(0)
{
  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;
  connectedExternalNodes(4) = node5 ;
  connectedExternalNodes(5) = node6 ;
  connectedExternalNodes(6) = node7 ;
  connectedExternalNodes(7) = node8 ;
  connectedExternalNodes(8) = node9 ;

  int i ;
  for ( i=0 ; i<9; i++ ) {

    materialPointers[i] = theMaterial.getCopy("AxiSymmetric2D") ;
    
    if (materialPointers[i] == 0) {

      opserr << "NineNodeMixedQuad::constructor() - failed to get a material of type: AxiSymmetric2D\n";
    } //end if
      
  } //end for i 
}


//destructor 
NineNodeMixedQuad :: ~NineNodeMixedQuad( )
{
  int i ;
  for ( i=0 ; i<9; i++ ) {

    delete materialPointers[i] ;
    materialPointers[i] = 0 ; 

    nodePointers[i] = 0 ;
  } //end for i

  if (load != 0)
    delete load;

  if (Ki != 0)
    delete Ki;
}


//set domain
void 
NineNodeMixedQuad::setDomain( Domain *theDomain ) 
{  
  for ( int i = 0; i<9; i++ ) 
    nodePointers[i] = theDomain->getNode( connectedExternalNodes(i)  ) ;
  
  this->DomainComponent::setDomain(theDomain);
}


//get the number of external nodes
int 
NineNodeMixedQuad::getNumExternalNodes( ) const
{
  return 9 ;
} 
 

//return connected external nodes
const ID& 
NineNodeMixedQuad::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 

Node **
NineNodeMixedQuad::getNodePtrs(void) 
{
  return nodePointers;
} 


//return number of dofs
int 
NineNodeMixedQuad::getNumDOF( ) 
{
  return 18 ;
}


//commit state
int 
NineNodeMixedQuad::commitState( )
{
  int i ;
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "NineNodeMixedQuad::commitState () - failed in base class\n";
  }    

  for ( i=0; i<9; i++ ) 
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}
 


//revert to last commit 
int 
NineNodeMixedQuad::revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i=0; i<9; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int 
NineNodeMixedQuad::revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i=0; i<9; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

//print out element data
void 
NineNodeMixedQuad::Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "Nine Node Quad -- Mixed Pressure/Volume -- Plane Strain \n" ;
  s << "Element Number " << this->getTag() << endln ;
  s << "Node 1 : " << connectedExternalNodes(0) << endln ;
  s << "Node 2 : " << connectedExternalNodes(1) << endln ;
  s << "Node 3 : " << connectedExternalNodes(2) << endln ;
  s << "Node 4 : " << connectedExternalNodes(3) << endln ;
  s << "Node 5 : " << connectedExternalNodes(4) << endln ;
  s << "Node 6 : " << connectedExternalNodes(5) << endln ;
  s << "Node 7 : " << connectedExternalNodes(6) << endln ;
  s << "Node 8 : " << connectedExternalNodes(7) << endln ;
  s << "Node 9 : " << connectedExternalNodes(8) << endln ;

  s << "Material Information : \n " ;

  materialPointers[0]->Print( s, flag ) ;

  s << endln ;
}


//return stiffness matrix 
const Matrix& 
NineNodeMixedQuad::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    


//return secant matrix 
const Matrix& 
NineNodeMixedQuad::getInitialStiff( ) 
{
  if (Ki != 0)
    return *Ki;

  static const int ndm = 2 ;

  static const int ndf = 2 ; 

  static const int nstress = 4 ;
 
  static const int numberNodes = 9 ;

  static const int numberGauss = 9 ;

  static const int nShape = 3 ;

  static const int nMixed = 3 ;

  int i, j, k, p, q, r, s ;
  int jj, kk ;

  static double volume ;

  static double xsj ;  // determinant jacaobian matrix 

  static double dvol[numberGauss] ; //volume element

  static double gaussPoint[ndm] ;

  static double natCoorArray[ndm][numberGauss] ;

  static Vector strain(nstress) ;  //strain

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

  static double shpBar[nShape][numberNodes][nMixed] ; //mean value of shape functions

  static double rightHandSide[nShape][numberNodes][nMixed] ;

  static Vector residJ(ndf) ; //nodeJ residual 

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 

  static Vector stress(nstress) ;  //stress

  static Matrix dd(nstress,nstress) ;  //material tangent

  static double interp[nMixed] ;

  static Matrix Proj(3,3) ;   //projection matrix 
  static Matrix ProjInv(3,3) ;

  static Matrix Iden(3,3) ;
  Iden(0,0) = 1.0 ;
  Iden(1,1) = 1.0 ;
  Iden(2,2) = 1.0 ;

  //---------B-matrices------------------------------------

    static Matrix BJ(nstress,ndf) ;      // B matrix node J

    static Matrix BJtran(ndf,nstress) ;

    static Matrix BK(nstress,ndf) ;      // B matrix node k

    static Matrix BJtranD(ndf,nstress) ;

  //-------------------------------------------------------

  
  //zero stiffness and residual 
  stiff.Zero( ) ;

  //node coordinates
  computeBasis() ;

  //zero mean shape functions
  for ( p=0; p<nShape; p++ ) {
    for ( q=0; q<numberNodes; q++ ) {

      for (r=0; r<nMixed; r++ ) {
	shpBar[p][q][r] = 0.0 ;
	rightHandSide[p][q][r] = 0.0 ;
      }

    }//end for q
  } // end for p


  //zero volume
  volume = 0.0 ;

  //zero projection matrix  
  Proj.Zero( ) ;
  ProjInv.Zero( ) ;

  //gauss loop to compute and save shape functions 
  int count = 0 ;

  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ ) {

        gaussPoint[0] = sg[i] ;        
	gaussPoint[1] = sg[j] ;        


	//save gauss point locations
	natCoorArray[0][count] = gaussPoint[0] ;
	natCoorArray[1][count] = gaussPoint[1] ;


	//get shape functions    
	shape2dNine( gaussPoint, xl, shp, xsj ) ;


	//save shape functions
	for ( p=0; p<nShape; p++ ) {
	  for ( q=0; q<numberNodes; q++ )
	    Shape[p][q][count] = shp[p][q] ;
	} // end for p

	
	//volume element to also be saved
	dvol[count] = ( wg[i]*wg[j] ) * xsj ;  


        //add to projection matrix
	interp[0] = 1.0 ;
	interp[1] = gaussPoint[0] ;
	interp[2] = gaussPoint[1] ;
	
	for ( r=0; r<nMixed; r++ ) {
	  for ( s=0; s<nMixed; s++ ) 
	    Proj(r,s) += ( interp[r]*interp[s] * dvol[count] ) ;
	}//end for r

	volume += dvol[count] ;
	
	
	//add to mean shape functions
	for ( p=0; p<nShape; p++ ) {
	  for ( q=0; q<numberNodes; q++ ) {

	    for ( s=0; s<nMixed; s++ ) 
	      rightHandSide[p][q][s] += ( shp[p][q] * interp[s] * dvol[count] ) ;

	  }//end for q 
	} // end for p


	//increment gauss point counter
	count++ ;

    } //end for j
  } // end for i 
  


  //invert projection matrix
  //int Solve(const Matrix &M, Matrix &res) const;
  Proj.Solve( Iden, ProjInv ) ;
  
  //mean value of shape functions
  for ( p=0; p<nShape; p++ ) {
    for ( q=0; q<numberNodes; q++ ) {

      for (r=0; r<nMixed; r++ ) {
	for (s=0; s<nMixed; s++ ) 
	  shpBar[p][q][r] += ( ProjInv(r,s) * rightHandSide[p][q][s] ) ;
      }//end for r

    }//end for q
  }//end for p


  //gauss loop 
  for ( i=0; i<numberGauss; i++ ) {
    
    //extract gauss point location
    gaussPoint[0] = natCoorArray[0][i] ;
    gaussPoint[1] = natCoorArray[1][i] ;

    //extract shape functions from saved array
    for ( p=0; p<nShape; p++ ) {
       for ( q=0; q<numberNodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
    } // end for p

    dd = materialPointers[i]->getInitialTangent( ) ;
    dd *= dvol[i] ;


    //residual and tangent calculations node loops

    jj = 0 ;
    for ( j=0; j<numberNodes; j++ ) {

      BJ = computeBbar( j, gaussPoint, shp, shpBar ) ;
   
      //transpose 
      //BJtran = transpose( nstress, ndf, BJ ) ;
      for (p=0; p<ndf; p++) {
	for (q=0; q<nstress; q++) 
	  BJtran(p,q) = BJ(q,p) ;
      }//end for p


      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixProduct(0.0,  BJtran,dd,1.0);

      kk = 0 ;
      for ( k=0; k<numberNodes; k++ ) {
	
	BK = computeBbar( k, gaussPoint, shp, shpBar ) ;
  
	
	//stiffJK =  BJtranD * BK  ;
	stiffJK.addMatrixProduct(0.0,  BJtranD,BK,1.0) ;

	for ( p=0; p<ndf; p++ )  {
	  for ( q=0; q<ndf; q++ )
	    stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
	} //end for p
	
	kk += ndf ;
      }//end for k loop

      jj += ndf ;
    }//end for j loop
  }//end for i gauss loop 

  Ki = new Matrix(stiff);

  return stiff;
}
    

//return mass matrix
const Matrix& 
NineNodeMixedQuad::getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 



void 
NineNodeMixedQuad::zeroLoad( )
{
  if (load != 0)
    load->Zero();
}



int 
NineNodeMixedQuad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "NineNodeMixedQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}

int
NineNodeMixedQuad::addInertiaLoadToUnbalance(const Vector &accel)
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
const Vector& 
NineNodeMixedQuad::getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  // subtract external loads 
  if (load != 0)
    resid -= *load;

  return resid ;   
}


//get residual with inertia terms
const Vector& 
NineNodeMixedQuad::getResistingForceIncInertia( )
{
  int tang_flag = 0 ; //don't get the tangent

  static Vector res(18);

  //do tangent and residual here 
  formResidAndTangent( tang_flag ) ;

  //inertia terms
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

//*****************************************************************************
//form inertia terms

void   
NineNodeMixedQuad::formInertiaTerms( int tangFlag ) 
{

  static const int ndm = 2 ;

  static const int ndf = 2 ; 

  static const int numberNodes = 9 ;

  //  static const int numberGauss = 9 ;

  static const int nShape = 3 ;

  static const int massIndex = nShape - 1 ;

  double xsj ;  // determinant jacaobian matrix 

  double dvol ; //volume element

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static Vector momentum(ndf) ;

  static Matrix sx(ndm,ndm) ;

  static double GaussPoint[2] ;

  int j, k, p, q, r ;
  int jj, kk ;

  double temp, rho, massJK ;


  //zero mass 
  mass.Zero( ) ;

  //node coordinates
  computeBasis() ;

  //gauss loop 
  int count = 0 ;
  for ( p=0; p<3; p++ ) {
    for ( q=0; q<3; q++ ) {
 
    GaussPoint[0] = sg[p] ;
    GaussPoint[1] = sg[q] ;

    //get shape functions    
    shape2dNine( GaussPoint, xl, shp, xsj ) ;

    //volume element
    dvol = ( wg[p] * wg[q] ) * xsj ;

    //node loop to compute acceleration
    momentum.Zero( ) ;
    for ( j = 0; j < numberNodes; j++ ) 
      //momentum += shp[massIndex][j] * ( nodePointers[j]->getTrialAccel()  ) ; 
      momentum.addVector( 1.0,
			  nodePointers[j]->getTrialAccel(),
			  shp[massIndex][j] ) ;


    //density
    rho = materialPointers[count]->getRho() ;

    //multiply acceleration by density to form momentum
    momentum *= rho ;


    //residual and tangent calculations node loops
    for ( jj=0, j=0; j<numberNodes; j++, jj+=ndf ) {

      temp = shp[massIndex][j] * dvol ;

      for ( r=0; r<ndf; r++ )
        resid( jj+r ) += ( temp * momentum(r) )  ;

      
      if ( tangFlag == 1 ) {

	 //multiply by density
	 temp *= rho ;

	 //node-node mass
         for ( kk=0, k=0; k<numberNodes; k++, kk+=ndf ) {

	    massJK = temp * shp[massIndex][k] ;

            for ( r=0; r<ndf; r++ )  
	      mass( jj+r, kk+r ) += massJK ;
            
          } // end for k loop

      } // end if tang_flag 

    } // end for j loop

    
    count++ ;
    }//end for q gauss loop
  } //end for p gauss loop 

}

//*********************************************************************
//form residual and tangent

void  
NineNodeMixedQuad::formResidAndTangent( int tang_flag ) 
{

  //strains ordered : eps11, eps22, eps33, 2*eps12 
  //volumtric strains projected onto {1, \xi, \eta} natural coordinates

  static const int ndm = 2 ;

  static const int ndf = 2 ; 

  static const int nstress = 4 ;
 
  static const int numberNodes = 9 ;

  static const int numberGauss = 9 ;

  static const int nShape = 3 ;

  static const int nMixed = 3 ;

  int i, j, k, p, q, r, s ;
  int jj, kk ;

  int success ;
  
  static double volume ;

  static double xsj ;  // determinant jacaobian matrix 

  static double dvol[numberGauss] ; //volume element

  static double gaussPoint[ndm] ;

  static double natCoorArray[ndm][numberGauss] ;

  static Vector strain(nstress) ;  //strain

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

  static double shpBar[nShape][numberNodes][nMixed] ; //mean value of shape functions

  static double rightHandSide[nShape][numberNodes][nMixed] ;

  static Vector residJ(ndf) ; //nodeJ residual 

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 

  static Vector stress(nstress) ;  //stress

  static Matrix dd(nstress,nstress) ;  //material tangent

  static double interp[nMixed] ;

  static Matrix Proj(3,3) ;   //projection matrix 
  static Matrix ProjInv(3,3) ;

  static Matrix Iden(3,3) ;
  Iden(0,0) = 1.0 ;
  Iden(1,1) = 1.0 ;
  Iden(2,2) = 1.0 ;

  //---------B-matrices------------------------------------

    static Matrix BJ(nstress,ndf) ;      // B matrix node J

    static Matrix BJtran(ndf,nstress) ;

    static Matrix BK(nstress,ndf) ;      // B matrix node k

    static Matrix BJtranD(ndf,nstress) ;

  //-------------------------------------------------------

  
  //zero stiffness and residual 
  stiff.Zero( ) ;
  resid.Zero( ) ;

  //node coordinates
  computeBasis() ;

  //zero mean shape functions
  for ( p=0; p<nShape; p++ ) {
    for ( q=0; q<numberNodes; q++ ) {

      for (r=0; r<nMixed; r++ ) {
	shpBar[p][q][r] = 0.0 ;
	rightHandSide[p][q][r] = 0.0 ;
      }

    }//end for q
  } // end for p


  //zero volume
  volume = 0.0 ;

  //zero projection matrix  
  Proj.Zero( ) ;
  ProjInv.Zero( ) ;


  //gauss loop to compute and save shape functions 
  int count = 0 ;

  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ ) {

        gaussPoint[0] = sg[i] ;        
	gaussPoint[1] = sg[j] ;        


	//save gauss point locations
	natCoorArray[0][count] = gaussPoint[0] ;
	natCoorArray[1][count] = gaussPoint[1] ;


	//get shape functions    
	shape2dNine( gaussPoint, xl, shp, xsj ) ;


	//save shape functions
	for ( p=0; p<nShape; p++ ) {
	  for ( q=0; q<numberNodes; q++ )
	    Shape[p][q][count] = shp[p][q] ;
	} // end for p

	
	//volume element to also be saved
	dvol[count] = ( wg[i]*wg[j] ) * xsj ;  


        //add to projection matrix
	interp[0] = 1.0 ;
	interp[1] = gaussPoint[0] ;
	interp[2] = gaussPoint[1] ;
	
	for ( r=0; r<nMixed; r++ ) {
	  for ( s=0; s<nMixed; s++ ) 
	    Proj(r,s) += ( interp[r]*interp[s] * dvol[count] ) ;
	}//end for r

	volume += dvol[count] ;
	
	
	//add to mean shape functions
	for ( p=0; p<nShape; p++ ) {
	  for ( q=0; q<numberNodes; q++ ) {

	    for ( s=0; s<nMixed; s++ ) 
	      rightHandSide[p][q][s] += ( shp[p][q] * interp[s] * dvol[count] ) ;

	  }//end for q 
	} // end for p


	//increment gauss point counter
	count++ ;

    } //end for j
  } // end for i 
  


  //invert projection matrix
  //int Solve(const Matrix &M, Matrix &res) const;
  Proj.Solve( Iden, ProjInv ) ;
  
  //mean value of shape functions
  for ( p=0; p<nShape; p++ ) {
    for ( q=0; q<numberNodes; q++ ) {

      for (r=0; r<nMixed; r++ ) {
	for (s=0; s<nMixed; s++ ) 
	  shpBar[p][q][r] += ( ProjInv(r,s) * rightHandSide[p][q][s] ) ;
      }//end for r

    }//end for q
  }//end for p


  //gauss loop 
  for ( i=0; i<numberGauss; i++ ) {
    
    //extract gauss point location
    gaussPoint[0] = natCoorArray[0][i] ;
    gaussPoint[1] = natCoorArray[1][i] ;

    //extract shape functions from saved array
    for ( p=0; p<nShape; p++ ) {
       for ( q=0; q<numberNodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
    } // end for p


    //zero the strains
    strain.Zero( ) ;

    // j-node loop to compute strain 
    for ( j=0; j<numberNodes; j++ )  {

      //compute B matrix 

      BJ = computeBbar( j, gaussPoint, shp, shpBar ) ;
      
      //nodal displacements 
      const Vector &ul = nodePointers[j]->getTrialDisp( ) ;

      //compute the strain
      //strain += (BJ*ul) ; 
      strain.addMatrixVector(1.0,  BJ,ul,1.0 ) ;

    } // end for j
  


    //send the strain to the material 
    success = materialPointers[i]->setTrialStrain( strain ) ;

    //compute the stress
    stress = materialPointers[i]->getStress( ) ;


    //multiply by volume element
    stress  *= dvol[i] ;

    if ( tang_flag == 1 ) {
      dd = materialPointers[i]->getTangent( ) ;
      dd *= dvol[i] ;
    } //end if tang_flag


    //residual and tangent calculations node loops

    jj = 0 ;
    for ( j=0; j<numberNodes; j++ ) {

      BJ = computeBbar( j, gaussPoint, shp, shpBar ) ;
   
      //transpose 
      //BJtran = transpose( nstress, ndf, BJ ) ;
      for (p=0; p<ndf; p++) {
	for (q=0; q<nstress; q++) 
	  BJtran(p,q) = BJ(q,p) ;
      }//end for p


      //residual
      //residJ = BJtran * stress ;
      residJ.addMatrixVector(0.0,  BJtran,stress,1.0);

      //residual 
      for ( p=0; p<ndf; p++ )
        resid( jj + p ) += residJ(p)  ;


      if ( tang_flag == 1 ) {

	//BJtranD = BJtran * dd ;
	BJtranD.addMatrixProduct(0.0,  BJtran,dd,1.0);

         kk = 0 ;
         for ( k=0; k<numberNodes; k++ ) {

            BK = computeBbar( k, gaussPoint, shp, shpBar ) ;
  
 
            //stiffJK =  BJtranD * BK  ;
	    stiffJK.addMatrixProduct(0.0,  BJtranD,BK,1.0) ;

            for ( p=0; p<ndf; p++ )  {
               for ( q=0; q<ndf; q++ )
                  stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
            } //end for p

            kk += ndf ;
	 }//end for k loop

      }//end if tang_flag 

      jj += ndf ;
    }//end for j loop


  }//end for i gauss loop 

  
  return ;
}


//*************************************************************************
//compute local coordinates and basis

void   
NineNodeMixedQuad::computeBasis( ) 
{

  //nodal coordinates 

  for ( int i = 0; i<9; i++ ) {

       const Vector &coorI = nodePointers[i]->getCrds( ) ;

       xl[0][i] = coorI(0) ;
       xl[1][i] = coorI(1) ;

  }  //end for i 

}

//*************************************************************************
//compute Bbar

const Matrix&   
NineNodeMixedQuad::computeBbar( int node, 
			    const double natCoor[2],
			    const double shp[3][9], 
			    double shpBar[3][9][3] )
{

  static Matrix Bbar(4,2) ;

  static double Bdev[3][2] ;

  static double BbarVol[3][2] ;

  static const double one3 = 1.0/3.0 ;

  static double interp[3] ;

  static double c0, c1 ;

  int i, j ;

//---B Matrices in standard {1,2,3} mechanics notation---------
//
//                -                  -
//               |  2N,1    -N,2      | 
// Bdev =  (1/3) |  -N,1    2N,2      |  (3x2)
//               |  -N,1    -N,2      |   
//                -                  -
//
//                -                 -
//               |  N,1      N,2     | 
// Bvol =  (1/3) |  N,1      N,2     |  (3x2)
//               |  N,1      N,2     |   
//                -                 -
//
//                -                   -
//               |                     |
//               |    Bdev + Bvol      |
//   B       =   |                     | 
//               |---------------------|   (4x2)
//               |   N,2     N,1       |
//                -                   -       
//
//---------------------------------------------------------------


  Bbar.Zero( ) ;

  //deviatoric
  Bdev[0][0] = 2.0*shp[0][node] ;
  Bdev[0][1] =    -shp[1][node] ;

  Bdev[1][0] =    -shp[0][node] ;
  Bdev[1][1] = 2.0*shp[1][node] ;

  Bdev[2][0] =    -shp[0][node] ;
  Bdev[2][1] =    -shp[1][node] ;

  
  //volume interpolation functions
  interp[0] = 1.0 ;
  interp[1] = natCoor[0] ;
  interp[2] = natCoor[1] ; 

  //volumetric 
  c0 = 0.0 ;
  c1 = 0.0 ;

  for (i=0; i<3; i++) {

    c0 += ( shpBar[0][node][i] * interp[i] ) ; 
    c1 += ( shpBar[1][node][i] * interp[i] ) ; 

  }//end for i

  //standard displacement formulation
  //c0 = shp[0][node] ;
  //c1 = shp[1][node] ;

  BbarVol[0][0] = c0 ;
  BbarVol[0][1] = c1 ;

  BbarVol[1][0] = c0 ;
  BbarVol[1][1] = c1 ;

  BbarVol[2][0] = c0 ;
  BbarVol[2][1] = c1 ;


  //extensional terms
  for ( i=0; i<3; i++ ){
    for ( j=0; j<2; j++ ) 
      Bbar(i,j) = one3*( Bdev[i][j] + BbarVol[i][j] ) ;
  }//end for i


  //shear terms
  Bbar(3,0) = shp[1][node] ;
  Bbar(3,1) = shp[0][node] ;


  return Bbar ;
}

//**************************************************************************
//shape function routine for nine node quads

void 
NineNodeMixedQuad::shape2dNine( double coor[2], 
		            const double x[2][9], 
		            double shp[3][9], 
			    double &xsj )  
{ 
  static const int nNode = 9 ;
  static const int ndm   = 2 ;

  static const int node1[] = { 0,1,1,0, 2,1,2,0, 2 } ;
  static const int node2[] = { 0,0,1,1, 0,2,1,2, 2 } ;
  int n1, n2 ;

  int i, j, k, q ;

  static double xs[ndm][ndm] ;
  static double sx[ndm][ndm] ;

  double ss = coor[0] ;
  double tt = coor[1] ;
  double tempSS, tempTT ;

  //shape functions and derivatives in natural coordinates
  for ( q=0; q<nNode; q++ ) {
    
    n1 = node1[q] ;
    n2 = node2[q] ;

    tempSS = this->shape1d(1,n1,ss) ;
    tempTT = this->shape1d(1,n2,tt) ; 

    //function
    shp[2][q] =  tempSS * tempTT ;

    //ss-derivative
    shp[0][q] =  (this->shape1d(0,n1,ss)) * tempTT ;
    
    //tt-derivative
    shp[1][q] =  tempSS * (this->shape1d(0,n2,tt)) ;  

  }//end for q


  // Construct jacobian and its inverse
  for ( i=0; i<ndm; i++ ) {
    for ( j=0; j<ndm; j++ ) {
      
      xs[i][j] = 0.0 ;

      for ( k=0; k<nNode; k++ )
	  xs[i][j] += ( x[i][k] * shp[j][k] ) ;

    } //end for j
  }  // end for i 


  xsj = xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0] ;

  double xsjInv = 1.0/xsj ;

  sx[0][0] =  xs[1][1] * xsjInv ;
  sx[1][1] =  xs[0][0] * xsjInv ;
  sx[0][1] = -xs[0][1] * xsjInv ; 
  sx[1][0] = -xs[1][0] * xsjInv ; 


  //form global derivatives 
  double temp ;
  for ( i=0; i<nNode; i++ ) {
    temp      = shp[0][i]*sx[0][0] + shp[1][i]*sx[1][0] ;
    shp[1][i] = shp[0][i]*sx[0][1] + shp[1][i]*sx[1][1] ;
    shp[0][i] = temp ;
  } // end for i


  return ;
}
	
//***********************************************************************
//1d quadratic shape functions

double 
NineNodeMixedQuad::shape1d( int code, int node, double xi ) 
{

  //one-dimensional quadratic shape functions
  //
  // o------o------o
  // 0      2      1
  //

  double result ;

  static double oneHalf = 0.50 ;

  //shape functions
  if ( code == 1 ) {

    if ( node == 0 )
      result = oneHalf*xi*(xi-1.0) ;

    if ( node == 1 )  
      result = oneHalf*xi*(xi+1.0) ;

    if ( node == 2 ) 
      result = 1.0 - xi*xi ;

  }
  //shape function derivatives
  else if ( code == 0 ) {

    if ( node == 0 )
      result = oneHalf * ( 2.0*xi - 1.0 ) ;

    if ( node == 1 )  
      result = oneHalf * ( 2.0*xi + 1.0 ) ;

    if ( node == 2 ) 
      result = -2.0*xi ;

  }//end if

  return result ;
}

//***********************************************************************

int
NineNodeMixedQuad::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    const Vector &end1Crd = nodePointers[0]->getCrds();
    const Vector &end2Crd = nodePointers[4]->getCrds();	
    const Vector &end3Crd = nodePointers[1]->getCrds();	
    const Vector &end4Crd = nodePointers[5]->getCrds();	
    const Vector &end5Crd = nodePointers[2]->getCrds();
    const Vector &end6Crd = nodePointers[6]->getCrds();	
    const Vector &end7Crd = nodePointers[3]->getCrds();	
    const Vector &end8Crd = nodePointers[7]->getCrds();	

    const Vector &end1Disp = nodePointers[0]->getDisp();
    const Vector &end2Disp = nodePointers[4]->getDisp();
    const Vector &end3Disp = nodePointers[1]->getDisp();
    const Vector &end4Disp = nodePointers[5]->getDisp();
    const Vector &end5Disp = nodePointers[2]->getDisp();
    const Vector &end6Disp = nodePointers[6]->getDisp();
    const Vector &end7Disp = nodePointers[3]->getDisp();
    const Vector &end8Disp = nodePointers[7]->getDisp();


    static Matrix coords(8,3) ;
    static Vector values(8) ;
    static Vector P(8) ;

    coords.Zero( ) ;

    values(0) = 1 ;
    values(1) = 1 ;
    values(2) = 1 ;
    values(3) = 1 ;
    values(4) = 1 ;
    values(5) = 1 ;
    values(6) = 1 ;
    values(7) = 1 ;


    if (displayMode < 3 && displayMode > 0)
      P = this->getResistingForce();

    for (int i = 0; i < 2; i++) {
      coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
      coords(1,i) = end2Crd(i) + end2Disp(i)*fact;    
      coords(2,i) = end3Crd(i) + end3Disp(i)*fact;    
      coords(3,i) = end4Crd(i) + end4Disp(i)*fact;    
      coords(4,i) = end5Crd(i) + end5Disp(i)*fact;
      coords(5,i) = end6Crd(i) + end6Disp(i)*fact;    
      coords(6,i) = end7Crd(i) + end7Disp(i)*fact;    
      coords(7,i) = end8Crd(i) + end8Disp(i)*fact;    
      /*      if (displayMode < 3 && displayMode > 0)
	values(i) = P(displayMode*2+i);
      else
      values(i) = 1;  */
    }

    //opserr << coords;
    int error = 0;

    error += theViewer.drawPolygon (coords, values);

    return error;
}
   
//**************************************************************************


Response*
NineNodeMixedQuad::setResponse(const char **argv, int argc, 
			  Information &eleInfo, OPS_Stream &output)
{
  Response *theResponse =0;

  output.tag("ElementOutput");
  output.attr("eleType","NineNodeMixedQuad");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  output.attr("node3",connectedExternalNodes[2]);
  output.attr("node4",connectedExternalNodes[3]);
  output.attr("node5",connectedExternalNodes[4]);
  output.attr("node6",connectedExternalNodes[5]);
  output.attr("node7",connectedExternalNodes[6]);
  output.attr("node8",connectedExternalNodes[7]);
  output.attr("node9",connectedExternalNodes[8]);

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
    if (pointNum > 0 && pointNum <= 9) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",sg[pointNum-1]);
      output.attr("neta",sg[pointNum-1]);

      theResponse =  materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, eleInfo, output);
      
      output.endTag();

  } else if (strcmp(argv[0],"stresses") ==0) {

      for (int i=0; i<9; i++) {
	output.tag("GaussPoint");
	output.attr("number",i+1);
	output.attr("eta",sg[i]);
	output.attr("neta",sg[i]);

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

      theResponse =  new ElementResponse(this, 3, Vector(4*9));
    }
  }
	
  output.endTag(); // ElementOutput

  return theResponse;
}

int 
NineNodeMixedQuad::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {

    // Loop over the integration points
    static Vector stresses(4*9);
    int cnt = 0;
    for (int i = 0; i < 9; i++) {

      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStress();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      stresses(cnt+3) = sigma(3);
      cnt += 4;
    }
    return eleInfo.setVector(stresses);
	
  } else

    return -1;
}


int 
NineNodeMixedQuad::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(28);
  
  int i;
  for (i = 0; i < 9; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i+9) = matDbTag;
  }
  
  idData(18) = this->getTag();
  idData(19) = connectedExternalNodes(0);
  idData(20) = connectedExternalNodes(1);
  idData(21) = connectedExternalNodes(2);
  idData(22) = connectedExternalNodes(3);
  idData(23) = connectedExternalNodes(4);
  idData(24) = connectedExternalNodes(5);
  idData(25) = connectedExternalNodes(6);
  idData(26) = connectedExternalNodes(7);
  idData(27) = connectedExternalNodes(8);


  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING NineNodeMixedQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
      
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 9; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING NineNodeMixedQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;
}
    
int 
NineNodeMixedQuad::recvSelf (int commitTag, 
			     Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(28);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING NineNodeMixedQuad::recvSelf() - " << this->getTag() << "  failed to receive ID\n";
    return res;
  }

  this->setTag(idData(18));
  connectedExternalNodes(0) = idData(19);
  connectedExternalNodes(1) = idData(20);
  connectedExternalNodes(2) = idData(21);
  connectedExternalNodes(3) = idData(22);
  connectedExternalNodes(4) = idData(23);
  connectedExternalNodes(5) = idData(24);
  connectedExternalNodes(6) = idData(25);
  connectedExternalNodes(7) = idData(26);
  connectedExternalNodes(8) = idData(27);
  

  int i;

  if (materialPointers[0] == 0) {
    for (i = 0; i < 9; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+9);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
	opserr << "NineNodeMixedQuad::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;
	return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "NineNodeMixedQuad::recvSelf() - material " <<
	  i << "failed to recv itself\n";
	return res;
      }
    }
  }
  // Number of materials is the same, receive materials into current space
  else {
    for (i = 0; i < 9; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+9);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
	delete materialPointers[i];
	materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
	if (materialPointers[i] == 0) {
	  opserr << "NineNodeMixedQuad::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;
	  exit(-1);
	}
      }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "NineNodeMixedQuad::recvSelf() - material " <<
	  i << "failed to recv itself\n";
	return res;
      }
    }
  }
  
  return res;
}

