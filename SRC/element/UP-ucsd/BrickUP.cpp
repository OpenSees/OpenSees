/* *********************************************************************
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
  
///////////////////////////////////////////////////////////////////////////////
// Description: This file contains the class declaration for                 //
// BrickUP, an 8-node cubic element for solid-fluid fully coupled analysis.  //
// This implementation is a simplified u-p formulation of Biot theory        //
// (u - solid displacement, p - fluid pressure). Each element node has three //
// DOFs for u and 1 DOF for p.                                               //
//                                                                           //
// Written by Zhaohui Yang	(March 2004)                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
                                                                           
// $Revision: 1.1 $
// $Date: 2005-09-22 21:28:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/UP-ucsd/BrickUP.cpp,v $

// by Zhaohui Yang (Modified based on Ed "C++" Love's Brick element)
//
// Eight node BrickUP element
//

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <BrickUP.h>
#include <shp3d.h>
#include <Renderer.h>
#include <ElementResponse.h>


#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static data
double  BrickUP::xl[4][8] ;

Matrix  BrickUP::stiff(32,32) ;
Vector  BrickUP::resid(32) ;
Matrix  BrickUP::mass(32,32) ;
Matrix  BrickUP::damp(32,32) ;
    
//quadrature data
const double  BrickUP::root3 = sqrt(3.0) ;
const double  BrickUP::one_over_root3 = 1.0 / root3 ;

const double  BrickUP::sg[] = { -one_over_root3,  
			       one_over_root3  } ;

const double  BrickUP::wg[] = { 1.0, 1.0, 1.0, 1.0, 
                              1.0, 1.0, 1.0, 1.0  } ;

  

//null constructor
BrickUP::BrickUP( ) :
Element( 0, ELE_TAG_BrickUP ),
connectedExternalNodes(8), load(0), Ki(0), kc(0), rho(0)
{ 
  for (int i=0; i<8; i++ ) {
    materialPointers[i] = 0;
    nodePointers[i] = 0;
  }
  b[0] = b[1] = b[2] = 0.;
  perm[0] = perm[1] = perm[2] = 0.;
}


//*********************************************************************
//full constructor
BrickUP::BrickUP(  int tag, 
                         int node1,
                         int node2,
   	                     int node3,
                         int node4,
                         int node5,
                         int node6,
                         int node7,
			             int node8,
			 NDMaterial &theMaterial, double bulk, double rhof, 
			double p1, double p2, double p3,
		   double b1, double b2, double b3) :
Element( tag, ELE_TAG_BrickUP ),
connectedExternalNodes(8), load(0), Ki(0), kc(bulk), rho(rhof)
{
  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;

  connectedExternalNodes(4) = node5 ;
  connectedExternalNodes(5) = node6 ;
  connectedExternalNodes(6) = node7 ;
  connectedExternalNodes(7) = node8 ;

  int i ;
  for ( i=0; i<8; i++ ) {

      materialPointers[i] = theMaterial.getCopy("ThreeDimensional") ;

      if (materialPointers[i] == 0) {
	  opserr <<"BrickUP::constructor - failed to get a material of type: ThreeDimensional\n";
	  exit(-1);
      } //end if
      
  } //end for i 

	// Body forces
	b[0] = b1;
	b[1] = b2;
	b[2] = b3;
	// Permeabilities
  perm[0] = p1;
  perm[1] = p2;
  perm[2] = p3;
}
//******************************************************************


//destructor 
BrickUP::~BrickUP( )
{
  int i ;
  for ( i=0 ; i<8; i++ ) {

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
void  BrickUP::setDomain( Domain *theDomain ) 
{  
  int i,dof ;

  // Check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    for ( i=0; i<8; i++ ) 
    nodePointers[i] = 0;
	return;
  }

  //node pointers
  for ( i=0; i<8; i++ ) {
     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
     if (nodePointers[i] == 0) {
	   opserr << "FATAL ERROR BrickUP ("<<this->getTag()<<"): node not found in domain"<<endln;
	   return;
     }
     
     dof = nodePointers[i]->getNumberDOF();
     if (dof != 4) {
	   opserr << "FATAL ERROR BrickUP ("<<this->getTag()<<"): has differing number of DOFs at its nodes"<<endln;
	   return;
     }
  }

  this->DomainComponent::setDomain(theDomain);
}


//get the number of external nodes
int  BrickUP::getNumExternalNodes( ) const
{
  return 8 ;
} 
 

//return connected external nodes
const ID&  BrickUP::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 

//return connected external node
Node **  
BrickUP::getNodePtrs(void) 
{
  return nodePointers ;
} 


//return number of dofs
int  BrickUP::getNumDOF( ) 
{
  return 32 ;
}


//commit state
int  BrickUP::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "BrickUP::commitState () - failed in base class";
  }    

  for (int i=0; i<8; i++ ) 
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}
 


//revert to last commit 
int  BrickUP::revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i=0; i<8; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int  BrickUP::revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i=0; i<8; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

//print out element data
void  BrickUP::Print( OPS_Stream &s, int flag )
{

  if (flag == 2) {

    s << "#Brick\n";
    
    int i;
    const int numNodes = 8;
    const int nstress = 6 ;
    
    for (i=0; i<numNodes; i++) {
      const Vector &nodeCrd = nodePointers[i]->getCrds();
      const Vector &nodeDisp = nodePointers[i]->getDisp();
      s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << nodeCrd(2)
	<< " " << nodeDisp(0) << " " << nodeDisp(1) << " " << nodeDisp(2) << endln;
     }
    
    // spit out the section location & invoke print on the scetion
    const int numMaterials = 8;

    static Vector avgStress(7);
    static Vector avgStrain(nstress);
    avgStress.Zero();
    avgStrain.Zero();
    for (i=0; i<numMaterials; i++) {
      avgStress += materialPointers[i]->getCommittedStress();
      avgStrain += materialPointers[i]->getCommittedStrain();
    }
    avgStress /= numMaterials;
    avgStrain /= numMaterials;

    s << "#AVERAGE_STRESS ";
    for (i=0; i<7; i++)
      s << avgStress(i) << " " ;
    s << endln;

    s << "#AVERAGE_STRAIN ";
    for (i=0; i<nstress; i++)
      s << avgStrain(i) << " " ;
    s << endln;

    /*
    for (i=0; i<numMaterials; i++) {
      s << "#MATERIAL\n";
      //      materialPointers[i]->Print(s, flag); 
      s << materialPointers[i]->getStress(); 
    }
    */

  } else {

    s << endln ;
    s << "Eight Node BrickUP \n" ;
    s << "Element Number: " << this->getTag() << endln ;
    s << "Node 1 : " << connectedExternalNodes(0) << endln ;
    s << "Node 2 : " << connectedExternalNodes(1) << endln ;
    s << "Node 3 : " << connectedExternalNodes(2) << endln ;
    s << "Node 4 : " << connectedExternalNodes(3) << endln ;
    s << "Node 5 : " << connectedExternalNodes(4) << endln ;
    s << "Node 6 : " << connectedExternalNodes(5) << endln ;
    s << "Node 7 : " << connectedExternalNodes(6) << endln ;
    s << "Node 8 : " << connectedExternalNodes(7) << endln ;
    
    s << "Material Information : \n " ;
    materialPointers[0]->Print( s, flag ) ;
    
    s << endln ;
  }
}
 
 
//return stiffness matrix 
const Matrix&  BrickUP::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    


//return secant matrix 
//const Matrix&  BrickUP::getSecantStiff( ) 

const Matrix&  BrickUP::getInitialStiff( ) 
{
  if (Ki != 0)
    return *Ki;

  //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 
  static const int ndm = 3 ;
  static const int ndf = 3 ; 
  static const int ndff = 4 ; 
  static const int nstress = 6 ;
  static const int numberNodes = 8 ;
  static const int numberGauss = 8 ;
  static const int nShape = 4 ;

  int i, j, k, p, q ;
  int jj, kk ;

  static double volume ;
  static double xsj ;  // determinant jacaobian matrix 
  static double dvol[numberGauss] ; //volume element
  static double gaussPoint[ndm] ;
  static Vector strain(nstress) ;  //strain
  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point
  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions
  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 
  static Matrix dd(nstress,nstress) ;  //material tangent


  //---------B-matrices------------------------------------

    static Matrix BJ(nstress,ndf) ;      // B matrix node J

    static Matrix BJtran(ndf,nstress) ;

    static Matrix BK(nstress,ndf) ;      // B matrix node k

    static Matrix BJtranD(ndf,nstress) ;

  //-------------------------------------------------------

  
  //zero stiffness and residual 
  stiff.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

  //gauss loop to compute and save shape functions 

  int count = 0 ;
  volume = 0.0 ;

  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++ ) {
      for ( k = 0; k < 2; k++ ) {

        gaussPoint[0] = sg[i] ;        
	gaussPoint[1] = sg[j] ;        
	gaussPoint[2] = sg[k] ;

	//get shape functions    
	shp3d( gaussPoint, xsj, shp, xl ) ;

	//save shape functions
	for ( p = 0; p < nShape; p++ ) {
	  for ( q = 0; q < numberNodes; q++ )
	    Shape[p][q][count] = shp[p][q] ;
	} // end for p


	//volume element to also be saved
	dvol[count] = wg[count] * xsj ;  

	//volume += dvol[count] ;

	count++ ;

      } //end for k
    } //end for j
  } // end for i 
  

  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) {

    //extract shape functions from saved array
    for ( p = 0; p < nShape; p++ ) {
       for ( q = 0; q < numberNodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
    } // end for p


    dd = materialPointers[i]->getInitialTangent( ) ;
    dd *= dvol[i] ;
    
    jj = 0;
    for ( j = 0; j < numberNodes; j++ ) {

      BJ = computeB( j, shp ) ;
   
      //transpose 
      //BJtran = transpose( nstress, ndf, BJ ) ;
      for (p=0; p<ndf; p++) {
	for (q=0; q<nstress; q++) 
	  BJtran(p,q) = BJ(q,p) ;
      }//end for p

      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixProduct(0.0,  BJtran, dd, 1.0) ;
      
      kk = 0 ;
      for ( k = 0; k < numberNodes; k++ ) {
	
	BK = computeB( k, shp ) ;
	
	
	//stiffJK =  BJtranD * BK  ;
	stiffJK.addMatrixProduct(0.0,  BJtranD, BK, 1.0) ;
	
	for ( p = 0; p < ndf; p++ )  {
	  for ( q = 0; q < ndf; q++ )
	    stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
	} //end for p

	kk += ndff ;

      } // end for k loop

      jj += ndff ;

    } // end for j loop
  } //end for i gauss loop 

  Ki = new Matrix(stiff);

  return stiff ;
}    


//return mass matrix
const Matrix&  BrickUP::getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 


//return mass matrix
const Matrix&  BrickUP::getDamp( ) 
{
  int tangFlag = 1 ;

  formDampingTerms( tangFlag ) ;

  return damp ;
} 

void BrickUP::formDampingTerms( int tangFlag )
{
  static const int ndm = 3 ;
  static const int ndf = 3 ; 
  static const int ndff = 4 ;
  static const int numberNodes = 8 ;
  static const int numberGauss = 8 ;
  static const int numberDOFs = 32 ;
  static const int nShape = 4 ;
  static double volume ;
  static double xsj ;  // determinant jacaobian matrix 
  static double dvol[numberGauss] ; //volume element
  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point
  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions
  static double gaussPoint[ndm] ;
  static Vector a(ndff*numberNodes) ;

  int i, j, k, p, q, m, i1, j1;
  //int jj, kk ;
  //double temp, rhot, massJK ;


  //zero damp 
  damp.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

  //zero volume
  volume = 0.0 ;

  //gauss loop to compute and save shape functions 

  int count = 0 ;

  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++ ) {
      for ( k = 0; k < 2; k++ ) {

        gaussPoint[0] = sg[i] ;        
	gaussPoint[1] = sg[j] ;        
	gaussPoint[2] = sg[k] ;

	//get shape functions    
	shp3d( gaussPoint, xsj, shp, xl ) ;

	//save shape functions
	for ( p = 0; p < nShape; p++ ) {
	  for ( q = 0; q < numberNodes; q++ )
	    Shape[p][q][count] = shp[p][q] ;
	} // end for p


	//volume element to also be saved
	dvol[count] = wg[count] * xsj ;  

        //add to volume
	volume += dvol[count] ;

	count++ ;

      } //end for k
    } //end for j
  } // end for i 

  if (betaK != 0.0)
    damp.addMatrix(1.0, this->getTangentStiff(), betaK);      
  if (betaK0 != 0.0)
    damp.addMatrix(1.0, this->getInitialStiff(), betaK0);      
  if (betaKc != 0.0)
    damp.addMatrix(1.0, *Kc, betaKc);      


  if (alphaM != 0.0) {
	this->getMass();
    for (i = 0; i < numberDOFs; i += ndff) {
      for (j = 0; j < numberDOFs; j += ndff) {
        damp(i,j) += mass(i,j)*alphaM;
        damp(i+1,j+1) += mass(i+1,j+1)*alphaM;
        damp(i+2,j+2) += mass(i+2,j+2)*alphaM;
	  }
    }  
  }

  // Compute coupling matrix
  for (i = 0; i < numberDOFs; i += ndff) {
    i1 = i / ndff;
    for (j = 3; j < numberDOFs; j += ndff) {
      j1 = (j-3) / ndff;
      for (m = 0; m < numberGauss; m++) {
	    damp(i,j) += -dvol[m]*Shape[0][i1][m]*Shape[3][j1][m];
	    damp(i+1,j) += -dvol[m]*Shape[1][i1][m]*Shape[3][j1][m];
	    damp(i+2,j) += -dvol[m]*Shape[2][i1][m]*Shape[3][j1][m];
	  }
      damp(j,i) = damp(i,j);
      damp(j,i+1) = damp(i+1,j);
      damp(j,i+2) = damp(i+2,j);
    }
  }
  
  // Compute permeability matrix
  for (i = 3; i < numberDOFs; i += ndff) {
    int i1 = (i-3) / ndff;
    for (j = 3; j < numberDOFs; j += ndff) {
      int j1 = (j-3) / ndff;
      for (m = 0; m < numberGauss; m++) {
	    damp(i,j) -= dvol[m]*(perm[0]*Shape[0][i1][m]*Shape[0][j1][m] +
	                          perm[1]*Shape[1][i1][m]*Shape[1][j1][m]+
					          perm[2]*Shape[2][i1][m]*Shape[2][j1][m]);
	  }
    }
  }

  if (tangFlag == 0) {
    for ( k = 0; k < numberNodes; k++ ) {
      const Vector &vel = nodePointers[k]->getTrialVel();
	  for ( p = 0; p < ndff; p++ )
		  a( k*ndff+p ) = vel(p); 
	} // end for k loop

    resid.addMatrixVector(1.0,damp,a,1.0);
  } // end if tang_flag 

  return ;
} 


void  BrickUP::zeroLoad( )
{
  if (load != 0)
    load->Zero();

  return ;
}


int 
BrickUP::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "BrickUP::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  return -1;
}

int
BrickUP::addInertiaLoadToUnbalance(const Vector &accel)
{
  static const int numberNodes = 8 ;
  static const int numberGauss = 8 ;
  static const int ndf = 3 ; 
  static const int ndff = 4 ; 

  int i;

  // Compute mass matrix
  int tangFlag = 1 ;
  formInertiaTerms( tangFlag ) ;

  // store computed RV from nodes in resid vector
  int count = 0;
  for (i=0; i<numberNodes; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int j=0; j<ndf; j++)
      resid(count++) = Raccel(j);

	resid(count++) = 0.0;  
  }

  // create the load vector if one does not exist
  if (load == 0) 
    load = new Vector(numberNodes*ndff);

  // add -M * RV(accel) to the load vector
  load->addMatrixVector(1.0, mass, resid, -1.0);
  
  return 0;
}


//get residual
const Vector&  BrickUP::getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  if (load != 0)
    resid -= *load;

  return resid ;   
}


//get residual with inertia terms
const Vector&  BrickUP::getResistingForceIncInertia( )
{
  static Vector res(32);

  int tang_flag = 0 ; //don't get the tangent

  //do tangent and residual here 
  formResidAndTangent( tang_flag ) ;

  formInertiaTerms( tang_flag ) ;

  formDampingTerms( tang_flag ) ;

  res = resid;

  if (load != 0)
    res -= *load;

  return res;
}


//*********************************************************************
//form inertia terms

void   BrickUP::formInertiaTerms( int tangFlag ) 
{
  static const int ndm = 3 ;
  static const int ndf = 3 ; 
  static const int ndff = 4 ;
  static const int numberNodes = 8 ;
  static const int numberGauss = 8 ;
  static const int nShape = 4 ;
  static const int massIndex = nShape - 1 ;
  static double volume ;
  static double xsj ;  // determinant jacaobian matrix 
  static double dvol[numberGauss] ; //volume element
  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point
  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions
  static double gaussPoint[ndm] ;
  static Vector a(ndff*numberNodes) ;

  int i, j, k, p, q ;
  int jj, kk ;

  double temp, rhot, massJK ;


  //zero mass 
  mass.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

  //zero volume
  volume = 0.0 ;

  //gauss loop to compute and save shape functions 

  int count = 0 ;

  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++ ) {
      for ( k = 0; k < 2; k++ ) {

        gaussPoint[0] = sg[i] ;        
	gaussPoint[1] = sg[j] ;        
	gaussPoint[2] = sg[k] ;

	//get shape functions    
	shp3d( gaussPoint, xsj, shp, xl ) ;

	//save shape functions
	for ( p = 0; p < nShape; p++ ) {
	  for ( q = 0; q < numberNodes; q++ )
	    Shape[p][q][count] = shp[p][q] ;
	} // end for p


	//volume element to also be saved
	dvol[count] = wg[count] * xsj ;  

        //add to volume
	volume += dvol[count] ;

	count++ ;

      } //end for k
    } //end for j
  } // end for i 

  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) {

    //extract shape functions from saved array
    for ( p = 0; p < nShape; p++ ) {
       for ( q = 0; q < numberNodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
    } // end for p

    // average material density 
      rhot = mixtureRho(i);

    //mass and compressibility calculations node loops
    jj = 0 ;
    for ( j = 0; j < numberNodes; j++ ) {

      temp = shp[massIndex][j] * dvol[i] ;

	 //multiply by density
	 temp *= rhot ;

	 //node-node mass
         kk = 0 ;
         for ( k = 0; k < numberNodes; k++ ) {

	    massJK = temp * shp[massIndex][k] ;

            for ( p = 0; p < ndf; p++ )  
	          mass( jj+p, kk+p ) += massJK ;
            
            // Compute compressibility terms
            mass( jj+3, kk+3 ) += -dvol[i]*Shape[3][j][i]*Shape[3][k][i]/kc;

            kk += ndff ;
          } // end for k loop

      jj += ndff ;
    } // end for j loop

  } //end for i gauss loop 

  if ( tangFlag == 0 ) {
    for ( k = 0; k < numberNodes; k++ ) {
      const Vector &acc = nodePointers[k]->getTrialAccel();
	  for ( p = 0; p < ndff; p++ )
		  a( k*ndff+p ) = acc(p); 
	} // end for k loop

    resid.addMatrixVector(1.0,mass,a,1.0);
  } // end if tang_flag 
}

//*********************************************************************
//form residual and tangent
void  BrickUP::formResidAndTangent( int tang_flag ) 
{

  //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 

  static const int ndm = 3 ;
  static const int ndf = 3 ; 
  static const int ndff = 4 ; 
  static const int nstress = 6 ;
  static const int numberNodes = 8 ;
  static const int numberGauss = 8 ;
  static const int nShape = 4 ;

  int i, j, k, p, q ;
  int jj, kk ;

  int success ;
  
  static double volume ;
  static double xsj ;  // determinant jacaobian matrix 
  static double dvol[numberGauss] ; //volume element
  static double gaussPoint[ndm] ;
  static Vector strain(nstress) ;  //strain
  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point
  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions
  static Vector residJ(ndf) ; //nodeJ residual 
  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 
  static Vector stress(nstress) ;  //stress
  static Matrix dd(nstress,nstress) ;  //material tangent


  //---------B-matrices------------------------------------
    static Matrix BJ(nstress,ndf) ;      // B matrix node J
    static Matrix BJtran(ndf,nstress) ;
    static Matrix BK(nstress,ndf) ;      // B matrix node k
    static Matrix BJtranD(ndf,nstress) ;
  //-------------------------------------------------------

  
  //zero stiffness and residual 
  stiff.Zero( ) ;
  resid.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

  //gauss loop to compute and save shape functions 

  int count = 0 ;
  volume = 0.0 ;

  for ( i = 0; i < 2; i++ ) {
    for ( j = 0; j < 2; j++ ) {
      for ( k = 0; k < 2; k++ ) {

        gaussPoint[0] = sg[i] ;        
	gaussPoint[1] = sg[j] ;        
	gaussPoint[2] = sg[k] ;

	//get shape functions    
	shp3d( gaussPoint, xsj, shp, xl ) ;

	//save shape functions
	for ( p = 0; p < nShape; p++ ) {
	  for ( q = 0; q < numberNodes; q++ )
	    Shape[p][q][count] = shp[p][q] ;
	} // end for p


	//volume element to also be saved
	dvol[count] = wg[count] * xsj ;  

	//volume += dvol[count] ;

	count++ ;

      } //end for k
    } //end for j
  } // end for i 
  

  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) {

    //extract shape functions from saved array
    for ( p = 0; p < nShape; p++ ) {
       for ( q = 0; q < numberNodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
    } // end for p


	double rhot;

    //zero the strains
    strain.Zero( ) ;

    // j-node loop to compute strain 
    for ( j = 0; j < numberNodes; j++ )  {

      //compute B matrix 
      BJ = computeB( j, shp ) ;
      
      //nodal displacements 
      const Vector &ul = nodePointers[j]->getTrialDisp( ) ;
      Vector ul3(3);
	  ul3(0) = ul(0);
      ul3(1) = ul(1);
	  ul3(2) = ul(2);
      //compute the strain
      //strain += (BJ*ul) ; 
      strain.addMatrixVector(1.0,BJ,ul3,1.0 ) ;

    } // end for j
  
    //send the strain to the material 
    success = materialPointers[i]->setTrialStrain( strain ) ;  


    //residual and tangent calculations node loops

    if ( tang_flag == 1 ) {
      dd = materialPointers[i]->getTangent( ) ;
      dd *= dvol[i] ;
    } //end if tang_flag


    if ( tang_flag == 0 ) {
      //compute the stress
      stress = materialPointers[i]->getStress( ) ;

      //multiply by volume element
      stress  *= dvol[i] ;

      rhot = mixtureRho(i);
    }

    jj = 0 ;
    for ( j = 0; j < numberNodes; j++ ) {

      BJ = computeB( j, shp ) ;
   
      //transpose 
      //BJtran = transpose( nstress, ndf, BJ ) ;
      for (p=0; p<ndf; p++) {
	    for (q=0; q<nstress; q++) 
	      BJtran(p,q) = BJ(q,p) ;
      }//end for p

      if ( tang_flag == 0 ) {
        //residual
        //residJ = BJtran * stress ;
        residJ.addMatrixVector(0.0,  BJtran,stress,1.0);

        for ( p = 0; p < ndf; p++ ) {
          resid( jj + p ) += residJ(p)  ;

          // Subtract equiv. body forces from the nodes
	      resid( jj + p ) -= dvol[i]*rhot*b[p]*Shape[3][j][i];
		}

        // Subtract fluid body force
		resid( jj + 3 ) += dvol[i]*rho*(perm[0]*b[0]*Shape[0][j][i] +
                                        perm[1]*b[1]*Shape[1][j][i] +
										perm[2]*b[2]*Shape[2][j][i]);
      } // end if tang_flag

      if ( tang_flag == 1 ) {
	     //BJtranD = BJtran * dd ;
	     BJtranD.addMatrixProduct(0.0,  BJtran,dd,1.0) ;

         kk = 0 ;
         for ( k = 0; k < numberNodes; k++ ) {

            BK = computeB( k, shp ) ;
  
            //stiffJK =  BJtranD * BK  ;
	        stiffJK.addMatrixProduct(0.0,  BJtranD,BK,1.0) ;

            for ( p = 0; p < ndf; p++ )  {
               for ( q = 0; q < ndf; q++ )
                  stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
            } //end for p

            kk += ndff ;
          } // end for k loop

      } // end if tang_flag 

      jj += ndff ;
    } // end for j loop

  } //end for i gauss loop 

  
  return ;
}

double BrickUP::mixtureRho(int i)
{
  double rhoi, e, n;

  rhoi= materialPointers[i]->getRho();
  //e = 0.7;  //theMaterial[i]->getVoidRatio();
  //n = e / (1.0 + e);
  //return n * rho + (1.0-n) * rhoi;
  return rhoi;
}

//************************************************************************
//compute local coordinates and basis

void   BrickUP::computeBasis( ) 
{

  //nodal coordinates 

  int i ;
  for ( i = 0; i < 8; i++ ) {

       const Vector &coorI = nodePointers[i]->getCrds( ) ;

       xl[0][i] = coorI(0) ;
       xl[1][i] = coorI(1) ;
       xl[2][i] = coorI(2) ;

  }  //end for i 

}

//*************************************************************************
//compute B

const Matrix&   
BrickUP::computeB( int node, const double shp[4][8] )
{

  static Matrix B(6,3) ;

//---B Matrix in standard {1,2,3} mechanics notation---------
//
//                -                   -
//               | N,1      0     0    | 
//   B       =   |   0     N,2    0    |
//               |   0      0     N,3  |   (6x3)
//               | N,2     N,1     0   |
//               |   0     N,3    N,2  |
//               | N,3      0     N,1  |
//                -                   -       
//
//-------------------------------------------------------------------


  B.Zero( ) ;

  B(0,0) = shp[0][node] ;
  B(1,1) = shp[1][node] ;
  B(2,2) = shp[2][node] ;

  B(3,0) = shp[1][node] ;
  B(3,1) = shp[0][node] ;

  B(4,1) = shp[2][node] ;
  B(4,2) = shp[1][node] ;

  B(5,0) = shp[2][node] ;
  B(5,2) = shp[0][node] ;

  return B ;

}

//***********************************************************************

Matrix  BrickUP::transpose( int dim1, 
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

int  BrickUP::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(25);

  idData(24) = this->getTag();
  
  int i;
  for (i = 0; i < 8; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i+8) = matDbTag;
  }
  
  idData(16) = connectedExternalNodes(0);
  idData(17) = connectedExternalNodes(1);
  idData(18) = connectedExternalNodes(2);
  idData(19) = connectedExternalNodes(3);
  idData(20) = connectedExternalNodes(4);
  idData(21) = connectedExternalNodes(5);
  idData(22) = connectedExternalNodes(6);
  idData(23) = connectedExternalNodes(7);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING BrickUP::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }


  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 8; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING BrickUP::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;

}
    
int  BrickUP::recvSelf (int commitTag, 
		       Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(25);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING BrickUP::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(24));

  connectedExternalNodes(0) = idData(16);
  connectedExternalNodes(1) = idData(17);
  connectedExternalNodes(2) = idData(18);
  connectedExternalNodes(3) = idData(19);
  connectedExternalNodes(4) = idData(20);
  connectedExternalNodes(5) = idData(21);
  connectedExternalNodes(6) = idData(22);
  connectedExternalNodes(7) = idData(23);
  

  int i;
  if (materialPointers[0] == 0) {
    for (i = 0; i < 8; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
	opserr << "BrickUP::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }
  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (i = 0; i < 8; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
	delete materialPointers[i];
	materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
	if (materialPointers[i] == 0) {
	  opserr << "BrickUP::recvSelf() - Broker could not create NDMaterial of class type " <<
	    matClassTag << endln;
	  exit(-1);
	}
      materialPointers[i]->setDbTag(matDbTag);
      }
      // Receive the material

      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "BrickUP::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }

  return res;
}
//**************************************************************************

int
BrickUP::displaySelf(Renderer &theViewer, int displayMode, float fact)
{

    const Vector &end1Crd = nodePointers[0]->getCrds();
    const Vector &end2Crd = nodePointers[1]->getCrds();	
    const Vector &end3Crd = nodePointers[2]->getCrds();	
    const Vector &end4Crd = nodePointers[3]->getCrds();	

    const Vector &end5Crd = nodePointers[4]->getCrds();
    const Vector &end6Crd = nodePointers[5]->getCrds();	
    const Vector &end7Crd = nodePointers[6]->getCrds();	
    const Vector &end8Crd = nodePointers[7]->getCrds();	

    const Vector &end1Disp = nodePointers[0]->getDisp();
    const Vector &end2Disp = nodePointers[1]->getDisp();
    const Vector &end3Disp = nodePointers[2]->getDisp();
    const Vector &end4Disp = nodePointers[3]->getDisp();

    const Vector &end5Disp = nodePointers[4]->getDisp();
    const Vector &end6Disp = nodePointers[5]->getDisp();
    const Vector &end7Disp = nodePointers[6]->getDisp();
    const Vector &end8Disp = nodePointers[7]->getDisp();

    const Vector &stress1 = materialPointers[0]->getStress();
    const Vector &stress2 = materialPointers[1]->getStress();
    const Vector &stress3 = materialPointers[2]->getStress();
    const Vector &stress4 = materialPointers[3]->getStress();
    const Vector &stress5 = materialPointers[4]->getStress();
    const Vector &stress6 = materialPointers[5]->getStress();
    const Vector &stress7 = materialPointers[6]->getStress();
    const Vector &stress8 = materialPointers[7]->getStress();

    static Matrix coords(4,3);
    static Vector values(4);
    static Vector P(24) ;

    values(0) = 1 ;
    values(1) = 1 ;
    values(2) = 1 ;
    values(3) = 1 ;


    // for each face of the BrickUP we:
    //   1) determine the coordinates of the displaced point
    //   2) determine the value to be drawn, the stress at nearest gauss point in displayMode dirn
    //   3) get the renderer to draw the face

    int error = 0;
    int i;
    for (i = 0; i < 3; i++) {
      coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
      coords(1,i) = end2Crd(i) + end2Disp(i)*fact;    
      coords(2,i) = end3Crd(i) + end3Disp(i)*fact;    
      coords(3,i) = end4Crd(i) + end4Disp(i)*fact;
    }

    if (displayMode < 3 && displayMode > 0) {
      int index = displayMode - 1;
      values(0) = stress1(index);
      values(1) = stress2(index);
      values(2) = stress3(index);
      values(3) = stress4(index);
    }

    error += theViewer.drawPolygon (coords, values);

    for (i = 0; i < 3; i++) {
      coords(0,i) = end5Crd(i) + end5Disp(i)*fact;
      coords(1,i) = end6Crd(i) + end6Disp(i)*fact;
      coords(2,i) = end7Crd(i) + end7Disp(i)*fact;
      coords(3,i) = end8Crd(i) + end8Disp(i)*fact;
    }

    if (displayMode < 3 && displayMode > 0) {
      int index = displayMode - 1;
      values(0) = stress5(index);
      values(1) = stress6(index);
      values(2) = stress7(index);
      values(3) = stress8(index);
    }

    error += theViewer.drawPolygon (coords, values);

    for (i = 0; i < 3; i++) {
      coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
      coords(1,i) = end4Crd(i) + end4Disp(i)*fact;
      coords(2,i) = end8Crd(i) + end8Disp(i)*fact;
      coords(3,i) = end5Crd(i) + end5Disp(i)*fact;
    }

    if (displayMode < 3 && displayMode > 0) {
      int index = displayMode - 1;
      values(0) = stress1(index);
      values(1) = stress4(index);
      values(2) = stress8(index);
      values(3) = stress5(index);
    }

    error += theViewer.drawPolygon (coords, values);

    for (i = 0; i < 3; i++) {
      coords(0,i) = end2Crd(i) + end2Disp(i)*fact;
      coords(1,i) = end3Crd(i) + end3Disp(i)*fact;
      coords(2,i) = end7Crd(i) + end7Disp(i)*fact;
      coords(3,i) = end6Crd(i) + end6Disp(i)*fact;
    }
    if (displayMode < 3 && displayMode > 0) {
      int index = displayMode - 1;
      values(0) = stress2(index);
      values(1) = stress3(index);
      values(2) = stress7(index);
      values(3) = stress6(index);
    }

    error += theViewer.drawPolygon (coords, values);


    for (i = 0; i < 3; i++) {
      coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
      coords(1,i) = end2Crd(i) + end2Disp(i)*fact;
      coords(2,i) = end6Crd(i) + end6Disp(i)*fact;
      coords(3,i) = end5Crd(i) + end5Disp(i)*fact;
    }

    if (displayMode < 3 && displayMode > 0) {
      int index = displayMode - 1;
      values(0) = stress1(index);
      values(1) = stress2(index);
      values(2) = stress6(index);
      values(3) = stress5(index);
    }

    error += theViewer.drawPolygon (coords, values);

    for (i = 0; i < 3; i++) {
      coords(0,i) = end4Crd(i) + end4Disp(i)*fact;
      coords(1,i) = end3Crd(i) + end3Disp(i)*fact;
      coords(2,i) = end7Crd(i) + end7Disp(i)*fact;
      coords(3,i) = end8Crd(i) + end8Disp(i)*fact;
    }

    if (displayMode < 3 && displayMode > 0) {
      int index = displayMode - 1;
      values(0) = stress3(index);
      values(1) = stress4(index);
      values(2) = stress7(index);
      values(3) = stress8(index);
    }

    error += theViewer.drawPolygon (coords, values);

    return error;
}

Response*
BrickUP::setResponse(const char **argv, int argc, Information &eleInfo)
{
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
    return new ElementResponse(this, 1, resid);
  
  else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)
    return new ElementResponse(this, 2, stiff);
  
  else if (strcmp(argv[0],"mass") == 0)
	return new ElementResponse(this, 3, mass);

  else if (strcmp(argv[0],"damp") == 0)
	return new ElementResponse(this, 4, damp);

  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 8)
      return materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, eleInfo);
    else 
      return 0;
  } else if (strcmp(argv[0],"stresses") ==0) {
    return new ElementResponse(this, 5, Vector(48));
  }
  
  // otherwise response quantity is unknown for the BrickUP class
  else
    return 0;
}

int 
BrickUP::getResponse(int responseID, Information &eleInfo)
{
  static Vector stresses(48);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2)
    return eleInfo.setMatrix(this->getTangentStiff());

  else if (responseID == 3)
    return eleInfo.setMatrix(this->getMass());

  else if (responseID == 4)
    return eleInfo.setMatrix(this->getDamp());
    
  else if (responseID == 5) {
    
    // Loop over the integration points
    int cnt = 0;
    for (int i = 0; i < 8; i++) {
      
      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStress();
      stresses(cnt++) = sigma(0);
      stresses(cnt++) = sigma(1);
      stresses(cnt++) = sigma(2);
      stresses(cnt++) = sigma(3);
      stresses(cnt++) = sigma(4);
      stresses(cnt++) = sigma(5);
    }
    return eleInfo.setVector(stresses);
    
  }
  else
    
    return -1;
}
