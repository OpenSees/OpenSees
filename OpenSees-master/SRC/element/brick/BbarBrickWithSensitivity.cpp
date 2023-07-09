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

// $Revision: 1.1 $
// $Date: 2009-08-07 00:17:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/brick/BbarBrickWithSensitivity.cpp,v $

// Sensitivity Part: Quan Gu UCSD
//
// Eight node BbarBrickWithSensitivity element
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
#include <BbarBrickWithSensitivity.h>
#include <shp3d.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

//static data
double  BbarBrickWithSensitivity::xl[3][8] ;

Matrix  BbarBrickWithSensitivity::stiff(24,24) ;
Vector  BbarBrickWithSensitivity::resid(24) ;
Matrix  BbarBrickWithSensitivity::mass(24,24) ;


//quadrature data
const double  BbarBrickWithSensitivity::root3 = sqrt(3.0) ;
const double  BbarBrickWithSensitivity::one_over_root3 = 1.0 / root3 ;

const double  BbarBrickWithSensitivity::sg[] = { -one_over_root3,
			       one_over_root3  } ;

const double  BbarBrickWithSensitivity::wg[] = { 1.0, 1.0, 1.0, 1.0,
                              1.0, 1.0, 1.0, 1.0  } ;


# define ELE_TAG_BbarBrickWithSensitivity 1984587234

void* OPS_BbarBrickWithSensitivity()
{
    if (OPS_GetNumRemainingInputArgs() < 10) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: element Brick eleTag? Node1? Node2? Node3? Node4? Node5? Node6? Node7? Node 8? matTag?\n";
	return 0;
    }

    int idata[10];
    int num = 10;
    if (OPS_GetIntInput(&num,idata)<0) {
	opserr<<"WARNING: invalid integer data\n";
	return 0;
    }

    NDMaterial* mat = OPS_getNDMaterial(idata[9]);
    if (mat == 0) {
	opserr << "WARNING material not found\n";
	opserr << "material tag: " << idata[9];
	opserr << "\nBrick element: " << idata[0] << endln;
    }

    double data[3] = {0,0,0};
    num = OPS_GetNumRemainingInputArgs();
    if (num > 3) {
	num = 3;
    }
    if (num > 0) {
	if (OPS_GetDoubleInput(&num,data) < 0) {
	    opserr<<"WARNING: invalid double data\n";
	    return 0;
	}	
    }

    return new BbarBrickWithSensitivity(idata[0],idata[1],idata[2],idata[3],idata[4],idata[5],idata[6],idata[7],idata[8],*mat,data[0],data[1],data[2]);
}

//null constructor
BbarBrickWithSensitivity::BbarBrickWithSensitivity( ) :
Element( 0, ELE_TAG_BbarBrickWithSensitivity ),
connectedExternalNodes(8), applyLoad(0), load(0), Ki(0)
{
  for (int i=0; i<8; i++ ) {
    materialPointers[i] = 0;
    nodePointers[i] = 0;
  }
  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;
// AddingSensitivity:BEGIN ////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END /////////////////////////////////
}


//*********************************************************************
//full constructor
BbarBrickWithSensitivity::BbarBrickWithSensitivity(  int tag,
                         int node1,
                         int node2,
   	                 int node3,
                         int node4,
                         int node5,
                         int node6,
                         int node7,
			 int node8,
			 NDMaterial &theMaterial,
			 double b1, double b2, double b3) :
Element( tag, ELE_TAG_BbarBrickWithSensitivity ),
connectedExternalNodes(8), applyLoad(0), load(0), Ki(0)
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
	  opserr <<"BbarBrickWithSensitivity::constructor - failed to get a material of type: ThreeDimensional\n";
	  exit(-1);
      } //end if

  } //end for i

  // Body forces
  b[0] = b1;
  b[1] = b2;
  b[2] = b3;

// AddingSensitivity:BEGIN ////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END /////////////////////////////////
}
//******************************************************************


//destructor
BbarBrickWithSensitivity::~BbarBrickWithSensitivity( )
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
void  BbarBrickWithSensitivity::setDomain( Domain *theDomain )
{

  int i ;

  //node pointers
  for ( i=0; i<8; i++ )
     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;

  this->DomainComponent::setDomain(theDomain);

}


//get the number of external nodes
int  BbarBrickWithSensitivity::getNumExternalNodes( ) const
{
  return 8 ;
}


//return connected external nodes
const ID&  BbarBrickWithSensitivity::getExternalNodes( )
{
  return connectedExternalNodes ;
}

//return connected external node
Node **
BbarBrickWithSensitivity::getNodePtrs(void)
{
  return nodePointers ;
}


//return number of dofs
int  BbarBrickWithSensitivity::getNumDOF( )
{
  return 24 ;
}


//commit state
int  BbarBrickWithSensitivity::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "Brick::commitState () - failed in base class";
  }

  for (int i=0; i<8; i++ )
    success += materialPointers[i]->commitState( ) ;

  return success ;
}



//revert to last commit
int  BbarBrickWithSensitivity::revertToLastCommit( )
{
  int i ;
  int success = 0 ;

  for ( i=0; i<8; i++ )
    success += materialPointers[i]->revertToLastCommit( ) ;

  return success ;
}


//revert to start
int  BbarBrickWithSensitivity::revertToStart( )
{
  int i ;
  int success = 0 ;

  for ( i=0; i<8; i++ )
    success += materialPointers[i]->revertToStart( ) ;

  return success ;
}

//print out element data
void  BbarBrickWithSensitivity::Print( OPS_Stream &s, int flag )
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "Element Number: " << this->getTag();
        s << "     Node 1 : " << connectedExternalNodes(0);
        s << "     Node 2 : " << connectedExternalNodes(1);
        s << "     Node 3 : " << connectedExternalNodes(2);
        s << "     Node 4 : " << connectedExternalNodes(3);
        s << "     Node 5 : " << connectedExternalNodes(4);
        s << "     Node 6 : " << connectedExternalNodes(5);
        s << "     Node 7 : " << connectedExternalNodes(6);
        s << "     Node 8 : " << connectedExternalNodes(7) << endln;
        
        //s << "Material Information : \n " ;
        //materialPointers[0]->Print( s, flag ) ;
        
        //s << endln ;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"BbarBrickWithSensitivity\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
        for (int i = 1; i < 7; i++)
            s << connectedExternalNodes(i) << ", ";
        s << connectedExternalNodes(7) << "], ";
        s << "\"bodyForces\": [" << b[0] << ", " << b[1] << ", " << b[2] << "], ";
        s << "\"material\": \"" << materialPointers[0]->getTag() << "\"}";
    }
}

//return stiffness matrix
const Matrix&  BbarBrickWithSensitivity::getTangentStiff( )
{
  int tang_flag = 1 ; //get the tangent

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;

  return stiff ;
}


//return stiffness matrix
const Matrix&  BbarBrickWithSensitivity::getInitialStiff( )
{
  if (Ki != 0)
    return *Ki;

  //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31

  static const int ndm = 3 ;

  static const int ndf = 3 ;

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

  static double shpBar[nShape][numberNodes] ;  //mean value of shape functions

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

  //zero mean shape functions
  for ( p = 0; p < nShape; p++ ) {
    for ( q = 0; q < numberNodes; q++ )
      shpBar[p][q] = 0.0 ;
  } // end for p

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

	//add to mean shape functions
	for ( p = 0; p < nShape; p++ ) {
	  for ( q = 0; q < numberNodes; q++ )
	    shpBar[p][q] += ( dvol[count] * shp[p][q] ) ;
	} // end for p

	count++ ;

      } //end for k
    } //end for j
  } // end for i


  //mean value of shape functions
  for ( p = 0; p < nShape; p++ ) {
    for ( q = 0; q < numberNodes; q++ )
      shpBar[p][q] /= volume ;
  } // end for p


  //gauss loop
  for ( i = 0; i < numberGauss; i++ ) {

    //extract shape functions from saved array
    for ( p = 0; p < nShape; p++ ) {
       for ( q = 0; q < numberNodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
    } // end for p

    dd = materialPointers[i]->getInitialTangent( ) ;
    dd *= dvol[i] ;

    //residual and tangent calculations node loops

    jj = 0 ;
    for ( j = 0; j < numberNodes; j++ ) {

      BJ = computeBbar( j, shp, shpBar ) ;

      //transpose
      //BJtran = transpose( nstress, ndf, BJ ) ;
      for (p=0; p<ndf; p++) {
	for (q=0; q<nstress; q++)
	  BJtran(p,q) = BJ(q,p) ;
      }//end for p

      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixProduct(0.0,  BJtran,dd,1.0);

      kk = 0 ;
      for ( k = 0; k < numberNodes; k++ ) {

	BK = computeBbar( k, shp, shpBar ) ;

	//stiffJK =  BJtranD * BK  ;
	stiffJK.addMatrixProduct(0.0,  BJtranD,BK,1.0) ;

	for ( p = 0; p < ndf; p++ )  {
	  for ( q = 0; q < ndf; q++ )
	    stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
	} //end for p

	kk += ndf ;
      } // end for k loop

      jj += ndf ;

    } // end for j loop
  } //end for i gauss loop

  Ki = new Matrix(stiff);

  return stiff ;
}


//return mass matrix
const Matrix&  BbarBrickWithSensitivity::getMass( )
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
}


void  BbarBrickWithSensitivity::zeroLoad( )
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
BbarBrickWithSensitivity::addLoad(ElementalLoad *theLoad, double loadFactor)
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
    opserr << "BbarBrickWithSensitivity::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
  }

  return -1;
}

int
BbarBrickWithSensitivity::addInertiaLoadToUnbalance(const Vector &accel)
{
  static const int numberNodes = 8 ;
  static const int numberGauss = 8 ;
  static const int ndf = 3 ;

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

  // store computed RV for nodes in resid vector
  int count = 0;
  for (i=0; i<numberNodes; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int j=0; j<ndf; j++)
      resid(count++) = Raccel(j);
  }

  // create the load vector if one does not exist
  if (load == 0)
    load = new Vector(numberNodes*ndf);

  // add -M * RV(accel) to the load vector
  load->addMatrixVector(1.0, mass, resid, -1.0);

  return 0;
}


//get residual
const Vector&  BbarBrickWithSensitivity::getResistingForce( )
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  if (load != 0)
    resid -= *load;

  return resid ;
}


//get residual with inertia terms
const Vector&  BbarBrickWithSensitivity::getResistingForceIncInertia( )
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

  if (load != 0)
    res -= *load;

  return res ;
}


//*********************************************************************
//form inertia terms

void   BbarBrickWithSensitivity::formInertiaTerms( int tangFlag )
{

  static const int ndm = 3 ;

  static const int ndf = 3 ;

  static const int numberNodes = 8 ;

  static const int numberGauss = 8 ;

  static const int nShape = 4 ;

  static const int massIndex = nShape - 1 ;

  double xsj ;  // determinant jacaobian matrix

  double dvol[numberGauss] ; //volume element

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

  static double gaussPoint[ndm] ;

  static Vector momentum(ndf) ;

  int i, j, k, p, q ;
  int jj, kk ;

  double temp, rho, massJK ;


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



  //gauss loop
  for ( i = 0; i < numberGauss; i++ ) {

    //extract shape functions from saved array
    for ( p = 0; p < nShape; p++ ) {
       for ( q = 0; q < numberNodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
    } // end for p


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

      temp = shp[massIndex][j] * dvol[i] ;

      for ( p = 0; p < ndf; p++ )
        resid( jj+p ) += ( temp * momentum(p) )  ;


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

      jj += ndf ;
    } // end for j loop


  } //end for i gauss loop

}

//*********************************************************************
//form residual and tangent
void  BbarBrickWithSensitivity::formResidAndTangent( int tang_flag )
{


// Quan	change it. This code becomes BAD however at least will not run element state determination twice!

if (tang_flag ==0) {// stress only




		//strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31

  static const int ndm = 3 ;

  static const int ndf = 3 ;

  static const int nstress = 6 ;

  static const int numberNodes = 8 ;

  static const int numberGauss = 8 ;

  static const int nShape = 4 ;

  int i, j, k, p, q ;
  int jj ;

  int success ;

  static double volume ;

  static double xsj ;  // determinant jacaobian matrix

  static double dvol[numberGauss] ; //volume element

  static double gaussPoint[ndm] ;

  static Vector strain(nstress) ;  //strain

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

  static double shpBar[nShape][numberNodes] ;  //mean value of shape functions

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


  //zero mean shape functions
  for ( p = 0; p < nShape; p++ ) {
	for ( q = 0; q < numberNodes; q++ )
	  shpBar[p][q] = 0.0 ;
  } // end for p

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

	//add to mean shape functions
	for ( p = 0; p < nShape; p++ ) {
	  for ( q = 0; q < numberNodes; q++ )
		shpBar[p][q] += ( dvol[count] * shp[p][q] ) ;
	} // end for p

	count++ ;

	  } //end for k
	} //end for j
  } // end for i


  //mean value of shape functions
  for ( p = 0; p < nShape; p++ ) {
	for ( q = 0; q < numberNodes; q++ )
	  shpBar[p][q] /= volume ;
  } // end for p


  //gauss loop
  for ( i = 0; i < numberGauss; i++ ) {

	//extract shape functions from saved array
	for ( p = 0; p < nShape; p++ ) {
	   for ( q = 0; q < numberNodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
	} // end for p



	//zero the strains
	strain.Zero( ) ;


    // j-node loop to compute strain
    for ( j = 0; j < numberNodes; j++ )  {

	  //compute B matrix

	  BJ = computeBbar( j, shp, shpBar ) ;

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


	//residual and tangent calculations node loops

	jj = 0 ;
	for ( j = 0; j < numberNodes; j++ ) {

	  BJ = computeBbar( j, shp, shpBar ) ;

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
	  for ( p = 0; p < ndf; p++ ) {
		resid( jj + p ) += residJ(p)  ;
		if (applyLoad == 0) {
		  resid( jj + p ) -= dvol[i]*b[p]*shp[3][j];
		} else {
		  resid( jj + p ) -= dvol[i]*appliedB[p]*shp[3][j];
		}
	  }


	  jj += ndf ;
	} // end for j loop


  } //end for i gauss loop


}   else if (tang_flag ==1){   //Quan tangent only





	//strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31

  static const int ndm = 3 ;

  static const int ndf = 3 ;

  static const int nstress = 6 ;

  static const int numberNodes = 8 ;

  static const int numberGauss = 8 ;

  static const int nShape = 4 ;

  int i, j, k, p, q ;
  int jj, kk ;

//Quan	  int success ;

  static double volume ;

  static double xsj ;  // determinant jacaobian matrix

  static double dvol[numberGauss] ; //volume element

  static double gaussPoint[ndm] ;

//Quan	  static Vector strain(nstress) ;  //strain

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

  static double shpBar[nShape][numberNodes] ;  //mean value of shape functions

//Quan	  static Vector residJ(ndf) ; //nodeJ residual

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness

//Quan	  static Vector stress(nstress) ;  //stress

  static Matrix dd(nstress,nstress) ;  //material tangent


  //---------B-matrices------------------------------------

	static Matrix BJ(nstress,ndf) ;      // B matrix node J

	static Matrix BJtran(ndf,nstress) ;

	static Matrix BK(nstress,ndf) ;      // B matrix node k

	static Matrix BJtranD(ndf,nstress) ;

  //-------------------------------------------------------


  //zero stiffness and residual
  stiff.Zero( ) ;
//Quan	  resid.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;


  //zero mean shape functions
  for ( p = 0; p < nShape; p++ ) {
	for ( q = 0; q < numberNodes; q++ )
	  shpBar[p][q] = 0.0 ;
  } // end for p

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

	//add to mean shape functions
	for ( p = 0; p < nShape; p++ ) {
	  for ( q = 0; q < numberNodes; q++ )
		shpBar[p][q] += ( dvol[count] * shp[p][q] ) ;
	} // end for p

	count++ ;

	  } //end for k
	} //end for j
  } // end for i


  //mean value of shape functions
  for ( p = 0; p < nShape; p++ ) {
	for ( q = 0; q < numberNodes; q++ )
	  shpBar[p][q] /= volume ;
  } // end for p


  //gauss loop
  for ( i = 0; i < numberGauss; i++ ) {

	//extract shape functions from saved array
	for ( p = 0; p < nShape; p++ ) {
	   for ( q = 0; q < numberNodes; q++ )
	  shp[p][q]  = Shape[p][q][i] ;
	} // end for p


	//zero the strains
//Quan		strain.Zero( ) ;


// j-node loop to compute strain
/*Quan		for ( j = 0; j < numberNodes; j++ )  {

	  //compute B matrix

	  BJ = computeBbar( j, shp, shpBar ) ;

	  //nodal displacements
	  const Vector &ul = nodePointers[j]->getTrialDisp( ) ;

	  //compute the strain
	  //strain += (BJ*ul) ;
	  strain.addMatrixVector(1.0,  BJ,ul,1.0 ) ;

	} // end for j

 Quan */

		//send the strain to the material
//Quan		success = materialPointers[i]->setTrialStrain( strain ) ;

		//compute the stress
//Quan		stress = materialPointers[i]->getStress( ) ;


		//multiply by volume element
//Quan		stress  *= dvol[i] ;

//Quan		if ( tang_flag == 1 ) {
		  dd = materialPointers[i]->getTangent( ) ;
		  dd *= dvol[i] ;
//Quan		} //end if tang_flag


		//residual and tangent calculations node loops

		jj = 0 ;
		for ( j = 0; j < numberNodes; j++ ) {

		  BJ = computeBbar( j, shp, shpBar ) ;

		  //transpose
		  //BJtran = transpose( nstress, ndf, BJ ) ;
		  for (p=0; p<ndf; p++) {
		for (q=0; q<nstress; q++)
		  BJtran(p,q) = BJ(q,p) ;
		  }//end for p


		  //residual
		  //residJ = BJtran * stress ;
//Quan		  residJ.addMatrixVector(0.0,  BJtran,stress,1.0);

		  //residual
//Quan		  for ( p = 0; p < ndf; p++ ) {
//Quan			resid( jj + p ) += residJ(p)  ;
//Quan			resid( jj + p ) -= dvol[i]*b[p]*shp[3][j];
//Quan		  }

//Quan		  if ( tang_flag == 1 ) {

		//BJtranD = BJtran * dd ;
		BJtranD.addMatrixProduct(0.0,  BJtran,dd,1.0);

			 kk = 0 ;
			 for ( k = 0; k < numberNodes; k++ ) {

				BK = computeBbar( k, shp, shpBar ) ;


				//stiffJK =  BJtranD * BK  ;
			stiffJK.addMatrixProduct(0.0,  BJtranD,BK,1.0) ;

				for ( p = 0; p < ndf; p++ )  {
				   for ( q = 0; q < ndf; q++ )
					  stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
				} //end for p

				kk += ndf ;
			  } // end for k loop

//Quan		  } // end if tang_flag

		  jj += ndf ;
		} // end for j loop


	  } //end for i gauss loop

} // else Quan

  return ;
}


//************************************************************************
//compute local coordinates and basis

void   BbarBrickWithSensitivity::computeBasis( )
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
BbarBrickWithSensitivity::computeBbar( int node,
				 const double shp[4][8],
				 const double shpBar[4][8] )
{

  static Matrix Bbar(6,3) ;

  //static Matrix Bdev(3,3) ;
  static double Bdev[3][3] ;

  //static Matrix BbarVol(3,3) ;
  static double BbarVol[3][3] ;

  static const double one3 = 1.0/3.0 ;


//---B Matrices in standard {1,2,3} mechanics notation---------
//
//                -                        -
//               |  2N,1    -N,2     -N,3   |
// Bdev =  (1/3) |  -N,1    2N,2     -N,3   |  (3x3)
//               |  -N,1    -N,2     2N,3   |
//                -                        -
//
//                -                       -
//               |  N,1      N,2     N,3   |
// Bvol =  (1/3) |  N,1      N,2     N.3   |  (3x3)
//               |  N,1      N,2     N,3   |
//                -                       -
//
//                -                   -
//               |                     |
//               |    Bdev + Bvol      |
//   B       =   |                     |
//               |---------------------|   (6x3)
//               | N,2     N,1     0   |
//               |   0     N,3    N,2  |
//               | N,3      0     N,1  |
//                -                   -
//
//---------------------------------------------------------------


  Bbar.Zero( ) ;

  //deviatoric
  Bdev[0][0] = 2.0*shp[0][node] ;
  Bdev[0][1] =    -shp[1][node] ;
  Bdev[0][2] =    -shp[2][node] ;

  Bdev[1][0] =    -shp[0][node] ;
  Bdev[1][1] = 2.0*shp[1][node] ;
  Bdev[1][2] =    -shp[2][node] ;

  Bdev[2][0] =    -shp[0][node] ;
  Bdev[2][1] =    -shp[1][node] ;
  Bdev[2][2] = 2.0*shp[2][node] ;

  //volumetric
  BbarVol[0][0] = shpBar[0][node] ;
  BbarVol[0][1] = shpBar[1][node] ;
  BbarVol[0][2] = shpBar[2][node] ;

  BbarVol[1][0] = shpBar[0][node] ;
  BbarVol[1][1] = shpBar[1][node] ;
  BbarVol[1][2] = shpBar[2][node] ;

  BbarVol[2][0] = shpBar[0][node] ;
  BbarVol[2][1] = shpBar[1][node] ;
  BbarVol[2][2] = shpBar[2][node] ;



  //extensional terms
  for ( int i=0; i<3; i++ ){
    for ( int j=0; j<3; j++ )
      Bbar(i,j) = one3*( Bdev[i][j] + BbarVol[i][j] ) ;
  }//end for i


  //shear terms
  Bbar(3,0) = shp[1][node] ;
  Bbar(3,1) = shp[0][node] ;

  Bbar(4,1) = shp[2][node] ;
  Bbar(4,2) = shp[1][node] ;

  Bbar(5,0) = shp[2][node] ;
  Bbar(5,2) = shp[0][node] ;

  return Bbar ;

}

//***********************************************************************

Matrix  BbarBrickWithSensitivity::transpose( int dim1,
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





int  BbarBrickWithSensitivity::sendSelf (int commitTag, Channel &theChannel)
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
    opserr << "WARNING BbarBrickWithSensitivity::sendSelf() - " << this->getTag() << "failed to send ID\n";
    return res;
  }


  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 8; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING BbarBrickWithSensitivity::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  return res;
}

int  BbarBrickWithSensitivity::recvSelf (int commitTag,
		       Channel &theChannel,
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  static ID idData(25);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING BbarBrickWithSensitivity::recvSelf() - " << this->getTag() << " failed to receive ID\n";
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


  if (materialPointers[0] == 0) {
    for (int i = 0; i < 8; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
	  opserr << "BbarBrickWithSensitivity::recvSelf() - Broker could not create NDMaterial of class type" <<
	    matClassTag << endln;
	  exit(-1);
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
  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < 8; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
	delete materialPointers[i];
	materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
	if (materialPointers[i] == 0) {
	  opserr << "BbarBrickWithSensitivity::recvSelf() - Broker could not create NDMaterial of class type" <<
	    matClassTag << endln;
	  exit(-1);
	}
      materialPointers[i]->setDbTag(matDbTag);
      }
      // Receive the material

      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "NLBeamColumn3d::recvSelf() - material " <<
	  i << "failed to recv itself\n";
	return res;
      }
    }
  }

  return res;
}
//**************************************************************************

int
BbarBrickWithSensitivity::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// vertex display coordinate vectors
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	static Vector v5(3);
	static Vector v6(3);
	static Vector v7(3);
	static Vector v8(3);
	static Matrix coords(8, 3); // polygon coordinate matrix
	static Vector values(8); // color vector
	static Vector P(24);
	int i;

	// get display coords
	nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
	nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
	nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
	nodePointers[3]->getDisplayCrds(v4, fact, displayMode);
	nodePointers[4]->getDisplayCrds(v5, fact, displayMode);
	nodePointers[5]->getDisplayCrds(v6, fact, displayMode);
	nodePointers[6]->getDisplayCrds(v7, fact, displayMode);
	nodePointers[7]->getDisplayCrds(v8, fact, displayMode);

	// add to coord matrix
	for (i = 0; i < 3; i++) {
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
	if (displayMode > 0)
		for (i = 0; i < 8; i++)
			values(i) = 1.0;
	else
		for (i = 0; i < 8; i++)
			values(i) = 0.0;

	if (displayMode < 3 && displayMode > 0)
		P = this->getResistingForce();

	return theViewer.drawCube(coords, values, this->getTag());
}

Response*
BbarBrickWithSensitivity::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","BbarBrick");
  output.attr("eleTag",this->getTag());
  for (int i=1; i<=8; i++) {
    sprintf(outputData,"node%d",i);
    output.attr(outputData,nodePointers[i-1]->getTag());
  }

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    char outputData[10];
    for (int i=1; i<=8; i++) {
      sprintf(outputData,"P1_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P2_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P3_%d",i);
      output.tag("ResponseType",outputData);
    }

    theResponse = new ElementResponse(this, 1, resid);

  }   else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {



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
    theResponse =  new ElementResponse(this, 3, Vector(48));
  }

  output.endTag(); // ElementOutput
  return theResponse;
}

int
BbarBrickWithSensitivity::getResponse(int responseID, Information &eleInfo)
{
  static Vector stresses(48);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2)
    return eleInfo.setMatrix(this->getTangentStiff());

  else if (responseID == 3) {

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



//////////////////////// Add sensitivity /////////////////////////////////////////////




int
BbarBrickWithSensitivity::setParameter(const char **argv, int argc, Parameter &param)
{


	int numberGauss=8;

    if (strstr(argv[0],"material") != 0) {
		int ok;
		for (int i=0; i<numberGauss; i++) {
			ok = materialPointers[i]->setParameter(&argv[1], argc-1, param);
			if (ok < 0 ) {
				opserr<<"BbarBrickWithSensitivity::setParameter() can not setParameter for "<<i<<"th Gauss Point\n";
				return -1;
			}
		}
		return ok;
	}
	else{
		opserr<<"BbarBrickWithSensitivity can not setParameter!"<<endln;
    // otherwise parameter is unknown for the Truss class
	return -1;
	}


}


int
BbarBrickWithSensitivity::updateParameter(int parameterID, Information &info)
{

  opserr<<"warnning: BbarBrickWithSensitivity can not updateParameter!"<<endln;
  return -1;
}



int
BbarBrickWithSensitivity::activateParameter(int passedParameterID)
{

	int numberGauss=8;
	parameterID = passedParameterID;

	if (parameterID == 1) {
		// Don't treat this case for now
	}

	else if (parameterID==0) {
		// Go down to the materials and zero out the identifier
		int ok;
		for (int i=0; i<numberGauss; i++) {
			ok = materialPointers[i]->activateParameter(parameterID);
			if (ok < 0 ) {
				return -1;
			}
		}
	}
	else if (parameterID > 100) {
		// In this case the parameter belongs to the material
		int ok;
		for (int i=0; i<numberGauss; i++) {
			ok = materialPointers[i]->activateParameter(parameterID-100);
			if (ok < 0 ) {
				return -1;
			}
		}
	}
	else {
		opserr << "BbarBrickWithSensitivity::activateParameter() -- unknown parameter " << endln;
	}

	return 0;
}



const Matrix &
BbarBrickWithSensitivity::getKiSensitivity(int gradNumber)
{
	stiff.Zero();
	return stiff;
}

const Matrix &
BbarBrickWithSensitivity::getMassSensitivity(int gradNumber)
{
	mass.Zero();
	return mass;
}




// Need to do/////////////////////////////////////////////
/////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


const Vector &
BbarBrickWithSensitivity::getResistingForceSensitivity(int gradNumber)
{


		//strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31

	  static const int ndm = 3 ;

	  static const int ndf = 3 ;

	  static const int nstress = 6 ;

	  static const int numberNodes = 8 ;

	  static const int numberGauss = 8 ;

	  static const int nShape = 4 ;

	  int i, j, k, p, q ;
	  int jj ;


	  static int markGU=0;
//	  int success ;

	  static double volume ;

	  static double xsj ;  // determinant jacaobian matrix

	  static double dvol[numberGauss] ; //volume element

	  static double gaussPoint[ndm] ;

//	  static Vector strain(nstress) ;  //strain

	  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

	  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

	  static double shpBar[nShape][numberNodes] ;  //mean value of shape functions

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

	  resid.Zero( ) ;

	  //compute basis vectors and local nodal coordinates
	  computeBasis( ) ;


	  //zero mean shape functions
	  for ( p = 0; p < nShape; p++ ) {
		for ( q = 0; q < numberNodes; q++ )
		  shpBar[p][q] = 0.0 ;
	  } // end for p

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

		//add to mean shape functions
		for ( p = 0; p < nShape; p++ ) {
		  for ( q = 0; q < numberNodes; q++ )
			shpBar[p][q] += ( dvol[count] * shp[p][q] ) ;
		} // end for p

		count++ ;

		  } //end for k
		} //end for j
	  } // end for i


	  //mean value of shape functions
	  for ( p = 0; p < nShape; p++ ) {
		for ( q = 0; q < numberNodes; q++ )
		  shpBar[p][q] /= volume ;
	  } // end for p


	  //gauss loop
	  for ( i = 0; i < numberGauss; i++ ) {

		//extract shape functions from saved array
		for ( p = 0; p < nShape; p++ ) {
		   for ( q = 0; q < numberNodes; q++ )
		  shp[p][q]  = Shape[p][q][i] ;
		} // end for p



		//compute the stress
		stress = materialPointers[i]->getStressSensitivity(gradNumber,true) ;

        //opserr<<"soil material. gradNumber="<<gradNumber<<". and stressSens("<<i<<") is:" <<stress<<endln;
		//multiply by volume element
		stress  *= dvol[i] ;


		//residual and tangent calculations node loops

		jj = 0 ;
		for ( j = 0; j < numberNodes; j++ ) {

		  BJ = computeBbar( j, shp, shpBar ) ;

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
		  for ( p = 0; p < ndf; p++ ) {
			resid( jj + p ) += residJ(p)  ;
//			resid( jj + p ) -= dvol[i]*b[p]*shp[3][j];
		  }


		  jj += ndf ;
		} // end for j loop


	  } //end for i gauss loop





//opserr<<"BbarBrickWithSensitivity::getResistingForceSensitivity(int gradNumber=" <<gradNumber<<") resid is: "<<resid<<endln;

//if (fabs(resid(0)+6.918032043683422e-007)<1.e-22){
//	opserr<<"mark passed!"<<endln; markGU=1;}

return resid;

}


int
BbarBrickWithSensitivity::commitSensitivity(int gradNumber, int numGrads)
{
/*	const Vector &disp1 = theNodes[0]->getDispSensitivity(gradNumber);
	const Vector &disp2 = theNodes[1]->getDispSensitivity(gradNumber);
	const Vector &disp3 = theNodes[2]->getDispSensitivity(gradNumber);
	const Vector &disp4 = theNodes[3]->getDispSensitivity(gradNumber);

  static double u[2][4];

	u[0][0] = disp1(0);
	u[1][0] = disp1(1);
	u[0][1] = disp2(0);
	u[1][1] = disp2(1);
	u[0][2] = disp3(0);
	u[1][2] = disp3(1);
	u[0][3] = disp4(0);
	u[1][3] = disp4(1);


	refer
		sens(0) = theNodes[0]->getDispSensitivity(1,gradNumber);
		sens(1) = theNodes[0]->getDispSensitivity(2,gradNumber);
		sens(2) = theNodes[1]->getDispSensitivity(1,gradNumber);
		sens(3) = theNodes[1]->getDispSensitivity(2,gradNumber);
		sens(4) = theNodes[2]->getDispSensitivity(1,gradNumber);
		sens(5) = theNodes[2]->getDispSensitivity(2,gradNumber);
		sens(6) = theNodes[3]->getDispSensitivity(1,gradNumber);
		sens(7) = theNodes[3]->getDispSensitivity(2,gradNumber);

*/

		//strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31


	  static int markGUU=0;

	   const int ndm = 3 ;

	   const int ndf = 3 ;

	   const int nstress = 6 ;

	   const int numberNodes = 8 ;

	   const int numberGauss = 8 ;

	   const int nShape = 4 ;

	  int i, j, k, p, q ;
//	  int jj, kk ;

	  int success ;

	   double volume ;

	   double xsj ;  // determinant jacaobian matrix

	   double dvol[numberGauss] ; //volume element

	   double gaussPoint[ndm] ;

	  static Vector strain(nstress) ;  //strain

	   double shp[nShape][numberNodes] ;  //shape functions at a gauss point

	   double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

	   double shpBar[nShape][numberNodes] ;  //mean value of shape functions

//Quan	  static Vector residJ(ndf) ; //nodeJ residual

	  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness

//Quan	  static Vector stress(nstress) ;  //stress

	  static Matrix dd(nstress,nstress) ;  //material tangent


	  //---------B-matrices------------------------------------

		static Matrix BJ(nstress,ndf) ;      // B matrix node J

		static Matrix BJtran(ndf,nstress) ;

		static Matrix BK(nstress,ndf) ;      // B matrix node k

		static Matrix BJtranD(ndf,nstress) ;

	  //-------------------------------------------------------



	  //compute basis vectors and local nodal coordinates
	  computeBasis( ) ;


	  //zero mean shape functions
	  for ( p = 0; p < nShape; p++ ) {
		for ( q = 0; q < numberNodes; q++ )
		  shpBar[p][q] = 0.0 ;
	  } // end for p

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

			//add to mean shape functions
			for ( p = 0; p < nShape; p++ ) {
			  for ( q = 0; q < numberNodes; q++ )
				shpBar[p][q] += ( dvol[count] * shp[p][q] ) ;
			} // end for p

			count++ ;

			  } //end for k
			} //end for j
		  } // end for i


		  //mean value of shape functions
		  for ( p = 0; p < nShape; p++ ) {
			for ( q = 0; q < numberNodes; q++ )
			  shpBar[p][q] /= volume ;
		  } // end for p


	  //gauss loop
	  for ( i = 0; i < numberGauss; i++ ) {

		//extract shape functions from saved array
		for ( p = 0; p < nShape; p++ ) {
		   for ( q = 0; q < numberNodes; q++ )
		  shp[p][q]  = Shape[p][q][i] ;
		} // end for p


		//zero the strains
		strain.Zero( ) ;

		static Vector ul(3);
	// j-node loop to compute strain
		for ( j = 0; j < numberNodes; j++ )  {

			  //compute B matrix

			  BJ = computeBbar( j, shp, shpBar ) ;

			  //nodal displacements
//			  const Vector &ul = nodePointers[j]->getDispSensitivity(gradNumber) ;


     			ul(0) = nodePointers[j]->getDispSensitivity(1,gradNumber);
				ul(1) = nodePointers[j]->getDispSensitivity(2,gradNumber);
				ul(2) = nodePointers[j]->getDispSensitivity(3,gradNumber);


			  //compute the strain
			  //strain += (BJ*ul) ;
			  strain.addMatrixVector(1.0,  BJ,ul,1.0 ) ;



/*			  if (markGUU==1){
			  opserr<<"j is:"<<j<<endln;
			  opserr<<"BJ:"<<BJ<<endln;
			  opserr<<"ul"<<ul<<endln;
			  opserr<<"strain"<<strain<<endln;
			  }


*/

		} // end for j



		//send the strain to the material

//		opserr<< "bbarBrick commitSens. gradNumber="<<gradNumber<<". Gauss (" <<i<<") is:" <<strain<<endln;


//		if (fabs(strain(0)-3.174725487662866e-011)<1.e-25){
//			opserr<<"bbarBrick commitSens mark passed!"<<endln; markGUU=1;}


		success = materialPointers[i]->commitSensitivity(strain,gradNumber,numGrads ) ;


		//residual and tangent calculations node loops



	  } //end for i gauss loop


	return success;
}





