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
                                                                        
// $Revision: 1.15 $
// $Date: 2003-02-25 23:32:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/brick/BbarBrick.cpp,v $

// Ed "C++" Love
//
// Eight node BbarBrick element
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
#include <BbarBrick.h>
#include <shp3d.h>
#include <Renderer.h>
#include <ElementResponse.h>


#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static data
double  BbarBrick::xl[3][8] ;

Matrix  BbarBrick::stiff(24,24) ;
Vector  BbarBrick::resid(24) ;
Matrix  BbarBrick::mass(24,24) ;

    
//quadrature data
const double  BbarBrick::root3 = sqrt(3.0) ;
const double  BbarBrick::one_over_root3 = 1.0 / root3 ;

const double  BbarBrick::sg[] = { -one_over_root3,  
			       one_over_root3  } ;

const double  BbarBrick::wg[] = { 1.0, 1.0, 1.0, 1.0, 
                              1.0, 1.0, 1.0, 1.0  } ;

  

//null constructor
BbarBrick::BbarBrick( ) :
Element( 0, ELE_TAG_BbarBrick ),
connectedExternalNodes(8), load(0), Ki(0)
{ 
  for (int i=0; i<8; i++ ) {
    materialPointers[i] = 0;
    nodePointers[i] = 0;
  }
}


//*********************************************************************
//full constructor
BbarBrick::BbarBrick(  int tag, 
                         int node1,
                         int node2,
   	                 int node3,
                         int node4,
                         int node5,
                         int node6,
                         int node7,
			 int node8,
			 NDMaterial &theMaterial ) :
Element( tag, ELE_TAG_BbarBrick ),
connectedExternalNodes(8), load(0), Ki(0)
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
	  opserr <<"BbarBrick::constructor - failed to get a material of type: ThreeDimensional\n";
	  exit(-1);
      } //end if
      
  } //end for i 


}
//******************************************************************


//destructor 
BbarBrick::~BbarBrick( )
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
void  BbarBrick::setDomain( Domain *theDomain ) 
{  

  int i ;

  //node pointers
  for ( i=0; i<8; i++ ) 
     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;

  this->DomainComponent::setDomain(theDomain);

}


//get the number of external nodes
int  BbarBrick::getNumExternalNodes( ) const
{
  return 8 ;
} 
 

//return connected external nodes
const ID&  BbarBrick::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 

//return connected external node
Node **  
BbarBrick::getNodePtrs(void) 
{
  return nodePointers ;
} 


//return number of dofs
int  BbarBrick::getNumDOF( ) 
{
  return 24 ;
}


//commit state
int  BbarBrick::commitState( )
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
int  BbarBrick::revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i=0; i<8; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int  BbarBrick::revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i=0; i<8; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

//print out element data
void  BbarBrick::Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "Volume/Pressure Eight Node BbarBrick \n" ;
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

//return stiffness matrix 
const Matrix&  BbarBrick::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    


//return stiffness matrix 
const Matrix&  BbarBrick::getInitialStiff( ) 
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
const Matrix&  BbarBrick::getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 


void  BbarBrick::zeroLoad( )
{
  if (load != 0)
    load->Zero();

  return ;
}

int 
BbarBrick::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "BbarBrick::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}

int
BbarBrick::addInertiaLoadToUnbalance(const Vector &accel)
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

  // store computed RV fro nodes in resid vector
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
const Vector&  BbarBrick::getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  if (load != 0)
    resid -= *load;

  return resid ;   
}


//get residual with inertia terms
const Vector&  BbarBrick::getResistingForceIncInertia( )
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

void   BbarBrick::formInertiaTerms( int tangFlag ) 
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
void  BbarBrick::formResidAndTangent( int tang_flag ) 
{

  //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 

  static const int ndm = 3 ;

  static const int ndf = 3 ; 

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

    if ( tang_flag == 1 ) {
      dd = materialPointers[i]->getTangent( ) ;
      dd *= dvol[i] ;
    } //end if tang_flag


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
      for ( p = 0; p < ndf; p++ )
        resid( jj + p ) += residJ(p)  ;


      if ( tang_flag == 1 ) {

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

      } // end if tang_flag 

      jj += ndf ;
    } // end for j loop


  } //end for i gauss loop 

  
  return ;
}


//************************************************************************
//compute local coordinates and basis

void   BbarBrick::computeBasis( ) 
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
BbarBrick::computeBbar( int node, 
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

Matrix  BbarBrick::transpose( int dim1, 
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




    
int  BbarBrick::sendSelf (int commitTag, Channel &theChannel)
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
    opserr << "WARNING BbarBrick::sendSelf() - " << this->getTag() << "failed to send ID\n";
    return res;
  }


  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 8; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING BbarBrick::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;
}

int  BbarBrick::recvSelf (int commitTag, 
		       Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(25);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING BbarBrick::recvSelf() - " << this->getTag() << " failed to receive ID\n";
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
	  opserr << "BbarBrick::recvSelf() - Broker could not create NDMaterial of class type" <<
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
	  opserr << "BbarBrick::recvSelf() - Broker could not create NDMaterial of class type" <<
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
BbarBrick::displaySelf(Renderer &theViewer, int displayMode, float fact)
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

    static Matrix coords(4,3);
    static Vector values(4);
    static Vector P(24) ;

    values(0) = 1 ;
    values(1) = 1 ;
    values(2) = 1 ;
    values(3) = 1 ;

    if (displayMode < 3 && displayMode > 0)
      P = this->getResistingForce();

    int error = 0;

    int i;
    for (i = 0; i < 3; i++) {
      coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
      coords(1,i) = end2Crd(i) + end2Disp(i)*fact;
      coords(2,i) = end3Crd(i) + end3Disp(i)*fact;
      coords(3,i) = end4Crd(i) + end4Disp(i)*fact;
    }
    error += theViewer.drawPolygon (coords, values);

    for (i = 0; i < 3; i++) {
      coords(0,i) = end5Crd(i) + end5Disp(i)*fact;
      coords(1,i) = end6Crd(i) + end6Disp(i)*fact;
      coords(2,i) = end7Crd(i) + end7Disp(i)*fact;
      coords(3,i) = end8Crd(i) + end8Disp(i)*fact;
    }
    error += theViewer.drawPolygon (coords, values);

    for (i = 0; i < 3; i++) {
      coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
      coords(1,i) = end4Crd(i) + end4Disp(i)*fact;
      coords(2,i) = end8Crd(i) + end8Disp(i)*fact;
      coords(3,i) = end5Crd(i) + end5Disp(i)*fact;
    }
    error += theViewer.drawPolygon (coords, values);

    for (i = 0; i < 3; i++) {
      coords(0,i) = end2Crd(i) + end2Disp(i)*fact;
      coords(1,i) = end3Crd(i) + end3Disp(i)*fact;
      coords(2,i) = end7Crd(i) + end7Disp(i)*fact;
      coords(3,i) = end6Crd(i) + end6Disp(i)*fact;
    }
    error += theViewer.drawPolygon (coords, values);


    for (i = 0; i < 3; i++) {
      coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
      coords(1,i) = end2Crd(i) + end2Disp(i)*fact;
      coords(2,i) = end6Crd(i) + end6Disp(i)*fact;
      coords(3,i) = end5Crd(i) + end5Disp(i)*fact;
    }
    error += theViewer.drawPolygon (coords, values);

    for (i = 0; i < 3; i++) {
      coords(0,i) = end4Crd(i) + end4Disp(i)*fact;
      coords(1,i) = end3Crd(i) + end3Disp(i)*fact;
      coords(2,i) = end7Crd(i) + end7Disp(i)*fact;
      coords(3,i) = end8Crd(i) + end8Disp(i)*fact;
    }
    error += theViewer.drawPolygon (coords, values);


    return error;
}

Response*
BbarBrick::setResponse(const char **argv, int argc, Information &eleInfo)
{
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
    return new ElementResponse(this, 1, resid);
  
  else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)
    return new ElementResponse(this, 2, stiff);
  
  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 8)
      return materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, eleInfo);
    else 
      return 0;
  } else if (strcmp(argv[0],"stresses") ==0) {
    return new ElementResponse(this, 3, Vector(48));
  }
  
  // otherwise response quantity is unknown for the brick class
  else
    return 0;
}

int 
BbarBrick::getResponse(int responseID, Information &eleInfo)
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
