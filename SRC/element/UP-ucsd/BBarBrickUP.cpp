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
// Description: This file contains the class declaration for BBarBrickUP,    //
// an 8-node cubic element for solid-fluid fully coupled analysis.           //
// This implementation is a simplified u-p formulation of Biot theory        //
// (u - solid displacement, p - fluid pressure). Each element node has three //
// DOFs for u and 1 DOF for p.                                               //
// Constant volume/pressure integration (BBar method) is used for integration//
// of the volumetric component of solid phase and the fulid phase.           //
//                                                                           //
// Written by Zhaohui Yang	(September 2009)                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// $Revision: 1.2 $
// $Date: 2009-11-02 21:22:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/UP-ucsd/BBarBrickUP.cpp,v $

// by Zhaohui Yang (Modified based on Ed "C++" Love's Brick element)
//
// Eight node BBarBrickUP element
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
#include <BBarBrickUP.h>
#include <shp3d.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

void* OPS_BBarBrickUP()
{
    if (OPS_GetNDM() != 3 || OPS_GetNDF() != 4) {
	opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with QuadUP element\n";
	return 0;
    }
    if (OPS_GetNumRemainingInputArgs() < 15) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: element brickUP eleTag? N1? N2? N3? N4? N5? N6? N7? N8? matTag? bulk? rhof? perm_x? perm_y? perm_z? <b1? b2? b3?>\n";
	return 0;
    }

    // BBarBrickUPId, Nod[8], matID
    int tags[10];
    int num = 10;
    if (OPS_GetIntInput(&num,tags) < 0) {
	opserr<<"WARNING: invalid integer input\n";
	return 0;
    }

    NDMaterial* mat = OPS_getNDMaterial(tags[9]);
    if (mat == 0) {
	opserr << "WARNING material not found\n";
	opserr << "Material: " << tags[9];
	opserr << "\nBBarBrickUP element: " << tags[0] << endln;
	return 0;
    }

    // bk, r, perm1, perm2, perm3
    double data[5];
    num = 5;
    if (OPS_GetDoubleInput(&num,data) < 0) {
	opserr<<"WARNING: invalid double input\n";
	return 0;
    }

    // b1, b2, b3
    double opt[3] = {0,0,0};
    num = OPS_GetNumRemainingInputArgs();
    if (num > 3) {
	num = 3;
    }
    if (num > 0) {
	if (OPS_GetDoubleInput(&num,opt) < 0) {
	    opserr<<"WARNING: invalid double input\n";
	    return 0;
	}
    }

    return new BBarBrickUP(tags[0],tags[1],tags[2],tags[3],tags[4],
			   tags[5],tags[6],tags[7],tags[8],
			   *mat,data[0],data[1],data[2],data[3],data[4],
			   opt[0],opt[1],opt[2]);
}

//static data
double  BBarBrickUP::xl[4][8] ;
double  BBarBrickUP::Shape[4][8][8];
double  BBarBrickUP::shpBar[3][8];
double  BBarBrickUP::BBar[6][3][8][8];
double  BBarBrickUP::BBarp[3][8][8];
double  BBarBrickUP::dvol[8];

Matrix  BBarBrickUP::stiff(32,32) ;
Vector  BBarBrickUP::resid(32) ;
Matrix  BBarBrickUP::mass(32,32) ;
Matrix  BBarBrickUP::damp(32,32) ;

//quadrature data
const double  BBarBrickUP::root3 = sqrt(3.0) ;
const double  BBarBrickUP::one_over_root3 = 1.0 / root3 ;

const double  BBarBrickUP::sg[] = { -one_over_root3,
			       one_over_root3  } ;

const double  BBarBrickUP::wg[] = { 1.0, 1.0, 1.0, 1.0,
                              1.0, 1.0, 1.0, 1.0  } ;



//null constructor
BBarBrickUP::BBarBrickUP( ) :
Element( 0, ELE_TAG_BBarBrickUP ),
connectedExternalNodes(8), applyLoad(0), load(0), Ki(0), kc(0), rho(0)
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
BBarBrickUP::BBarBrickUP(int tag,
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
Element( tag, ELE_TAG_BBarBrickUP ),
connectedExternalNodes(8), applyLoad(0), load(0), Ki(0), kc(bulk), rho(rhof)
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
	  opserr <<"BBarBrickUP::constructor - failed to get a material of type: ThreeDimensional\n";
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
BBarBrickUP::~BBarBrickUP( )
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
void  BBarBrickUP::setDomain( Domain *theDomain )
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
	   opserr << "FATAL ERROR BBarBrickUP ("<<this->getTag()<<"): node not found in domain"<<endln;
	   return;
     }

     dof = nodePointers[i]->getNumberDOF();
     if (dof != 4) {
	   opserr << "FATAL ERROR BBarBrickUP ("<<this->getTag()<<"): has differing number of DOFs at its nodes"<<endln;
	   return;
     }
  }

  this->DomainComponent::setDomain(theDomain);
}


//get the number of external nodes
int  BBarBrickUP::getNumExternalNodes( ) const
{
  return 8 ;
}


//return connected external nodes
const ID&  BBarBrickUP::getExternalNodes( )
{
  return connectedExternalNodes ;
}

//return connected external node
Node **
BBarBrickUP::getNodePtrs(void)
{
  return nodePointers ;
}


//return number of dofs
int  BBarBrickUP::getNumDOF( )
{
  return 32 ;
}


//commit state
int  BBarBrickUP::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "BBarBrickUP::commitState () - failed in base class";
  }

  for (int i=0; i<8; i++ )
    success += materialPointers[i]->commitState( ) ;

  return success ;
}



//revert to last commit
int  BBarBrickUP::revertToLastCommit( )
{
  int i ;
  int success = 0 ;

  for ( i=0; i<8; i++ )
    success += materialPointers[i]->revertToLastCommit( ) ;

  return success ;
}


//revert to start
int  BBarBrickUP::revertToStart( )
{
  int i ;
  int success = 0 ;

  for ( i=0; i<8; i++ )
    success += materialPointers[i]->revertToStart( ) ;

  return success ;
}

//print out element data
void  BBarBrickUP::Print( OPS_Stream &s, int flag )
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
      avgStress += materialPointers[i]->getStress();
      avgStrain += materialPointers[i]->getStrain();
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
    s << "Eight Node BBarBrickUP \n" ;
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
const Matrix&  BBarBrickUP::getTangentStiff( )
{
  int tang_flag = 1 ; //get the tangent

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;

  return stiff ;
}


//return secant matrix
//const Matrix&  BBarBrickUP::getSecantStiff( )

const Matrix&  BBarBrickUP::getInitialStiff( )
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

  static double xsj ;  // determinant jacaobian matrix
  static double gaussPoint[ndm] ;
  static Vector strain(nstress) ;  //strain
  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness
  static Matrix dd(nstress,nstress) ;  //material tangent
  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  //---------B-matrices------------------------------------
  static Matrix BJ(nstress,ndf) ;      // B matrix node J
  static Matrix BJtran(ndf,nstress) ;
  static Matrix BK(nstress,ndf) ;      // B matrix node k
  static Matrix BJtranD(ndf,nstress) ;


  //zero stiffness and residual
  stiff.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

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

        count++ ;

      } //end for k
    } //end for j
  } // end for i

  computeBBar();

  //gauss loop
  for ( i = 0; i < numberGauss; i++ ) {

    dd = materialPointers[i]->getInitialTangent( ) ;
    dd *= dvol[i] ;

    jj = 0;
    for ( j = 0; j < numberNodes; j++ ) {

      BJ = computeB(j, i) ;

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

	    BK = computeB(k, i) ;

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
const Matrix&  BBarBrickUP::getMass( )
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
}


//return mass matrix
const Matrix&  BBarBrickUP::getDamp( )
{
  int tangFlag = 1 ;

  formDampingTerms( tangFlag ) ;

  return damp ;
}

void BBarBrickUP::formDampingTerms( int tangFlag )
{
  static const int ndm = 3 ;
  static const int ndf = 3 ;
  static const int ndff = 4 ;
  static const int numberNodes = 8 ;
  static const int numberGauss = 8 ;
  static const int numberDOFs = 32 ;
  static const int nShape = 4 ;
  static double xsj ;  // determinant jacaobian matrix
  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point
  static double gaussPoint[ndm] ;
  static Vector a(ndff*numberNodes) ;

  int i, j, k, p, q, m, i1, j1;

  //zero damp
  damp.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

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

	count++ ;

      } //end for k
    } //end for j
  } // end for i

  computeBBar();

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
	    damp(i,j) += -dvol[m]*Shape[3][j1][m]
	                 *(BBar[0][0][i1][m]+BBar[1][0][i1][m]+BBar[2][0][i1][m]);
	    damp(i+1,j) += -dvol[m]*Shape[3][j1][m]
	                 *(BBar[0][1][i1][m]+BBar[1][1][i1][m]+BBar[2][1][i1][m]);
	    damp(i+2,j) += -dvol[m]*Shape[3][j1][m]
	                 *(BBar[0][2][i1][m]+BBar[1][2][i1][m]+BBar[2][2][i1][m]);
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
	    damp(i,j) -= dvol[m]*(perm[0]*BBarp[0][i1][m]*BBarp[0][j1][m] +
	                          perm[1]*BBarp[1][i1][m]*BBarp[1][j1][m]+
					          perm[2]*BBarp[2][i1][m]*BBarp[2][j1][m]);
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


void  BBarBrickUP::zeroLoad( )
{
  if (load != 0)
    load->Zero();

  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;

  return ;
}


int
BBarBrickUP::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // Added option for applying body forces in load pattern: C.McGann, U.Washington
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_BrickSelfWeight) {
      applyLoad = 1;
      appliedB[0] += loadFactor * b[0];
      appliedB[1] += loadFactor * b[1];
      appliedB[2] += loadFactor * b[2];
    return 0;
  } else if (type == LOAD_TAG_SelfWeight) {
      // added compatibility with selfWeight class implemented for all continuum elements, C.McGann, U.W.
      applyLoad = 1;
	  appliedB[0] += loadFactor*data(0)*b[0];
	  appliedB[1] += loadFactor*data(1)*b[1];
	  appliedB[2] += loadFactor*data(2)*b[2];
	  return 0;
  } else {
    opserr << "BBarBrickUP::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
  }

  return -1;
}

int
BBarBrickUP::addInertiaLoadToUnbalance(const Vector &accel)
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
    for (int j=0; j<ndf; j++) resid(count++) = Raccel(j);

    resid(count++) = 0.0;
  }

  // create the load vector if one does not exist
  if (load == 0)  load = new Vector(numberNodes*ndff);

  // add -M * RV(accel) to the load vector
  load->addMatrixVector(1.0, mass, resid, -1.0);

  return 0;
}


//get residual
const Vector&  BBarBrickUP::getResistingForce( )
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  if (load != 0)
    resid -= *load;

  return resid ;
}


//get residual with inertia terms
const Vector&  BBarBrickUP::getResistingForceIncInertia( )
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

void   BBarBrickUP::formInertiaTerms( int tangFlag )
{
  static const int ndm = 3 ;
  static const int ndf = 3 ;
  static const int ndff = 4 ;
  static const int numberNodes = 8 ;
  static const int numberGauss = 8 ;
  static const int nShape = 4 ;
  static const int massIndex = nShape - 1 ;
  static double xsj ;  // determinant jacaobian matrix
  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point
  static double gaussPoint[ndm] ;
  static Vector a(ndff*numberNodes) ;

  int i, j, k, p, q ;
  int jj, kk ;

  double temp, rhot, massJK ;


  //zero mass
  mass.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

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

	    count++ ;

      } //end for k
    } //end for j
  } // end for i

  computeBBar();

  //gauss loop
  for ( i = 0; i < numberGauss; i++ ) {

    // average material density
      rhot = mixtureRho(i);

    //mass and compressibility calculations node loops
    jj = 0 ;
    for ( j = 0; j < numberNodes; j++ ) {

      temp = Shape[3][j][i] * dvol[i] ;

	  //multiply by density
	  temp *= rhot ;

	  //node-node mass
      kk = 0 ;
      for ( k = 0; k < numberNodes; k++ ) {

	    massJK = temp * Shape[3][k][i] ;

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
void  BBarBrickUP::formResidAndTangent( int tang_flag )
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

  static double xsj ;  // determinant jacaobian matrix
  static double gaussPoint[ndm] ;
  static Vector strain(nstress) ;  //strain
  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point
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

    	count++ ;

      } //end for k
    } //end for j
  } // end for i

  computeBBar();

  //gauss loop
  for ( i = 0; i < numberGauss; i++ ) {

	double rhot;

    //zero the strains
    strain.Zero( ) ;

    // j-node loop to compute strain
    for ( j = 0; j < numberNodes; j++ )  {

      //compute B matrix
      BJ = computeB( j, i ) ;

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

      BJ = computeB( j, i ) ;

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
	  if (applyLoad == 0) {
	    resid( jj + p ) -= dvol[i]*rhot*b[p]*Shape[3][j][i];
	  } else {
	    resid( jj + p ) -= dvol[i]*rhot*appliedB[p]*Shape[3][j][i];
	  }
	}

        // Subtract fluid body force
		if (applyLoad == 0) {
			resid( jj + 3 ) += dvol[i]*rho*(perm[0]*b[0]*BBarp[0][j][i] +
                                        perm[1]*b[1]*BBarp[1][j][i] +
										perm[2]*b[2]*BBarp[2][j][i]);
		} else {
			resid( jj + 3 ) += dvol[i]*rho*(perm[0]*appliedB[0]*BBarp[0][j][i] +
                                        perm[1]*appliedB[1]*BBarp[1][j][i] +
										perm[2]*appliedB[2]*BBarp[2][j][i]);
		}
      } // end if tang_flag

      if ( tang_flag == 1 ) {
	     //BJtranD = BJtran * dd ;
	     BJtranD.addMatrixProduct(0.0,  BJtran,dd,1.0) ;

         kk = 0 ;
         for ( k = 0; k < numberNodes; k++ ) {

            BK = computeB( k, i ) ;

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


double BBarBrickUP::mixtureRho(int i)
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

void   BBarBrickUP::computeBasis( )
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


void  BBarBrickUP::computeBBar()
{

  static double volume ;
  int i, j, k;

  volume = 0;

  for (i=0; i<3; i++) {     // Directions iteration
	 for (j=0; j<8; j++) {  // Nodes iteration
		 shpBar[i][j] = 0.;
	 }
  }

  for (k=0; k<8; k++) {        // Gauss points iteration
     for (i=0; i<3; i++) {     // Directions iteration
	    for (j=0; j<8; j++) {  // Nodes iteration
	       shpBar[i][j] += Shape[i][j][k] * dvol[k];
	    }
     }

     volume += dvol[k];
  }

  for (i=0; i<3; i++) {
	 for (j=0; j<8; j++) {
	    shpBar[i][j] /= volume;
	 }
  }


//----- See Tom Hughes, Chapter 4: Mixed and Penalty Methods---------
//
//            -                   -
//   -       | B,5     B,6    B,8  |
//   B   =   | B,4     B,7    B,8  |
//           | B,4     B,6    B,9  |   (6x3)
//           | B,2     B,1     0   |
//           |   0     B,3    B,2  |
//           | B,3      0     B,1  |
//            -                   -
//where
//         -
//  B,4 = (B,1 - B,1) / 3
//
//  B,5 = B,1 + B,4
//         -
//  B,6 = (B,2 - B,2) / 3
//
//  B,7 = B,2 + B,6
//         -
//  B,8 = (B,3 - B,3) / 3
//
//  B,9 = B,3 + B,8
//-------------------------------------------------------------------

  for (k=0; k<8; k++) {
	 for (j=0; j<8; j++) {

       BBar[0][0][j][k] = (2*Shape[0][j][k] + shpBar[0][j]) / 3. ;
       BBar[0][1][j][k] = (shpBar[1][j] - Shape[1][j][k]) / 3. ;
       BBar[0][2][j][k] = (shpBar[2][j] - Shape[2][j][k]) / 3. ;
       BBar[1][0][j][k] = (shpBar[0][j] - Shape[0][j][k]) / 3. ;
       BBar[1][1][j][k] = (2*Shape[1][j][k] + shpBar[1][j]) / 3. ;
       BBar[1][2][j][k] = BBar[0][2][j][k];
       BBar[2][0][j][k] = BBar[1][0][j][k];
       BBar[2][1][j][k] = BBar[0][1][j][k];
       BBar[2][2][j][k] = (2*Shape[2][j][k] + shpBar[2][j]) / 3. ;
       BBar[3][0][j][k] = Shape[1][j][k];
       BBar[3][1][j][k] = Shape[0][j][k];
       BBar[3][2][j][k] = 0;
       BBar[4][0][j][k] = 0;
       BBar[4][1][j][k] = Shape[2][j][k];
       BBar[4][2][j][k] = Shape[1][j][k];
       BBar[5][0][j][k] = Shape[2][j][k];
       BBar[5][1][j][k] = 0;
       BBar[5][2][j][k] = Shape[0][j][k];

       BBarp[0][j][k] = BBar[0][0][j][k];
	   BBarp[1][j][k] = BBar[1][1][j][k];
	   BBarp[2][j][k] = BBar[2][2][j][k];
     }
   }

}


const Matrix&
BBarBrickUP::computeB( int node, int Gauss )
{

  static Matrix B(6,3) ;

  int i, j;

  for (i=0; i<6; i++) {
	 for (j=0; j<3; j++) {

        B(i,j) = BBar[i][j][node][Gauss] ;
	 }
  }

  return B ;

}

//***********************************************************************

Matrix  BBarBrickUP::transpose( int dim1,
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

int  BBarBrickUP::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // BBarBrickUP packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(13);
  data(0) = this->getTag();
  data(1) = rho;
  data(2) = b[0];
  data(3) = b[1];
  data(4) = b[2];

  data(5) = alphaM;
  data(6) = betaK;
  data(7) = betaK0;
  data(8) = betaKc;

  data( 9) = kc;
  data(10) = perm[0];
  data(11) = perm[1];
  data(12) = perm[2];

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING BBarBrickUP::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }

  // Now BBarBrickUP sends the ids of its materials
  int matDbTag;

  static ID idData(24);

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
    opserr << "WARNING BBarBrickUP::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }


  // Finally, BBarBrickUP asks its material objects to send themselves
  for (i = 0; i < 8; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING BBarBrickUP::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  return res;

}

int  BBarBrickUP::recvSelf (int commitTag,
		       Channel &theChannel,
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // BBarBrickUP creates a Vector, receives the Vector and then sets the
  // internal data with the data in the Vector
  static Vector data(13);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING FourNodeQuadUP::recvSelf() - failed to receive Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  rho = data(1);
  b[0] = data(2);
  b[1] = data(3);
  b[2] = data(4);

  alphaM = data(5);
  betaK = data(6);
  betaK0 = data(7);
  betaKc = data(8);

  kc = data(9);
  perm[0] = data(10);
  perm[1] = data(11);
  perm[2] = data(12);

  static ID idData(24);
  // BBarBrickUP now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING BBarBrickUP::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

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
    // Allocate new materials
    //materialPointers = new NDMaterial[8];
    //if (materialPointers == 0) {
    //  opserr << "BBarBrickUP::recvSelf() - Could not allocate NDMaterial array\n";
    //  return -1;
    //}
    for (i = 0; i < 8; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
	opserr << "BBarBrickUP::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "BBarBrickUP::recvSelf() - material " << i << "failed to recv itself\n";
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
	  opserr << "BBarBrickUP::recvSelf() - Broker could not create NDMaterial of class type " <<
	    matClassTag << endln;
	  exit(-1);
	}
      }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "BBarBrickUP::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }

  return res;
}
//**************************************************************************

int
BBarBrickUP::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // get vertex display coordinate vectors
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    static Vector v5(3);
    static Vector v6(3);
    static Vector v7(3);
    static Vector v8(3);
    nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
    nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
    nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
    nodePointers[3]->getDisplayCrds(v4, fact, displayMode);
    nodePointers[4]->getDisplayCrds(v5, fact, displayMode);
    nodePointers[5]->getDisplayCrds(v6, fact, displayMode);
    nodePointers[6]->getDisplayCrds(v7, fact, displayMode);
    nodePointers[7]->getDisplayCrds(v8, fact, displayMode);

    // add to coord matrix
    static Matrix coords(8, 3);
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v2(i);
        coords(2, i) = v3(i);
        coords(3, i) = v4(i);
        coords(4, i) = v5(i);
        coords(5, i) = v6(i);
        coords(6, i) = v7(i);
        coords(7, i) = v8(i);
    }

    // create color vector
    static Vector values(8);
    if (displayMode < 3 && displayMode > 0) {
        // get stress vectors
        const Vector& stress1 = materialPointers[0]->getStress();
        const Vector& stress2 = materialPointers[1]->getStress();
        const Vector& stress3 = materialPointers[2]->getStress();
        const Vector& stress4 = materialPointers[3]->getStress();
        const Vector& stress5 = materialPointers[4]->getStress();
        const Vector& stress6 = materialPointers[5]->getStress();
        const Vector& stress7 = materialPointers[6]->getStress();
        const Vector& stress8 = materialPointers[7]->getStress();
        // apply to color value vector
        int index = displayMode - 1;
        values(0) = stress1(index);
        values(1) = stress2(index);
        values(2) = stress3(index);
        values(3) = stress4(index);
        values(4) = stress5(index);
        values(5) = stress6(index);
        values(6) = stress7(index);
        values(7) = stress8(index);
    }
    else {
        // default color
        for (int i = 0; i < 8; i++)
            values(i) = 1.0;
    }

    // draw cube
    return theViewer.drawCube(coords, values, this->getTag());
}


Response*
BBarBrickUP::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","BBarBrickUP");
  output.attr("eleTag",this->getTag());
  for (int i=1; i<=8; i++) {
    sprintf(outputData,"node%d",i);
    output.attr(outputData,nodePointers[i-1]->getTag());
  }

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=8; i++) {
      sprintf(outputData,"P1_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P2_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P3_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"Pp_%d",i);
      output.tag("ResponseType",outputData);
    }

    theResponse = new ElementResponse(this, 1, resid);

  }   else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0) {
    theResponse = new ElementResponse(this, 2, stiff);

  }   else if (strcmp(argv[0],"mass") == 0) {
    theResponse = new ElementResponse(this, 3, mass);

  }   else if (strcmp(argv[0],"damp") == 0) {
    theResponse = new ElementResponse(this, 4, damp);

  } else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 8) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);

      theResponse =  materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, output);

      output.endTag(); // GaussPoint
    }

  } else if (strcmp(argv[0],"stresses") ==0) {

    for (int i=0; i<8; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());

      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma33");
      output.tag("ResponseType","sigma12");
      output.tag("ResponseType","sigma13");
      output.tag("ResponseType","sigma23");

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }
    theResponse = new ElementResponse(this, 5, Vector(48));
  }

  output.endTag(); // ElementOutput
  return theResponse;

}

int
BBarBrickUP::getResponse(int responseID, Information &eleInfo)
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

int
BBarBrickUP::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  // permeability in horizontal direction
  if (strcmp(argv[0],"hPerm") == 0) {
    return param.addObject(3, this);

  // permeability in vertical direction
  } else if (strcmp(argv[0],"vPerm") == 0) {
    return param.addObject(4, this);

  } else {

    int matRes = res;
    for (int i=0; i<8; i++) {
        matRes =  materialPointers[i]->setParameter(argv, argc, param);
        if (matRes != -1)
            res = matRes;
    }
  }

    return res;
}

int
BBarBrickUP::updateParameter(int parameterID, Information &info)
{
  int res = -1;
  int matRes = res;
  switch( parameterID ) {
	case 3:
		perm[0] = info.theDouble;
		this->getDamp();	// update mass matrix
		return 0;
	case 4:
		perm[1] = info.theDouble;
		perm[2] = info.theDouble;
		this->getDamp();	// update mass matrix
		return 0;
	default:
		return -1;
  }
  return -1;
}

