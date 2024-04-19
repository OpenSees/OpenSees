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
                                                                        
// $Revision: 1.31 $
// $Date: 2010-04-23 22:56:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/brick/Brick.cpp,v $

// Ed "C++" Love
//
// Eight node Brick element
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
#include <Brick.h>
#include <shp3d.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

void* OPS_Brick()
{
    int dampingTag = 0;
    Damping* theDamping = 0;

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
    //option,Tang.S
    int numData = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        std::string theType = OPS_GetString();
        if (theType == "-damp") {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
                theDamping = OPS_getDamping(dampingTag);
                if (theDamping == 0) {
                    opserr << "damping not found\n";
                    return 0;
                }
            }
        }
    }

    return new Brick(idata[0],idata[1],idata[2],idata[3],idata[4],idata[5],idata[6],idata[7],
		     idata[8],*mat,data[0],data[1],data[2], theDamping);
}

//static data
double  Brick::xl[3][8] ;

Matrix  Brick::stiff(24,24) ;
Vector  Brick::resid(24) ;
Matrix  Brick::mass(24,24) ;

    
//quadrature data
const double  Brick::root3 = sqrt(3.0) ;
const double  Brick::one_over_root3 = 1.0 / root3 ;

const double  Brick::sg[] = { -one_over_root3,  
			       one_over_root3  } ;

const double  Brick::wg[] = { 1.0, 1.0, 1.0, 1.0, 
                              1.0, 1.0, 1.0, 1.0  } ;

  
static Matrix B(6,3) ;

//null constructor
Brick::Brick( ) 
:Element( 0, ELE_TAG_Brick ),
 connectedExternalNodes(8), applyLoad(0), load(0), Ki(0)
{
  B.Zero();

  for (int i=0; i<8; i++ ) {
    materialPointers[i] = 0;
    nodePointers[i] = 0;
  }

  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;

  for (int i = 0 ;  i < 8; i++ ) 
    theDamping[i] = 0;
}


//*********************************************************************
//full constructor
Brick::Brick(int tag, 
	     int node1,
	     int node2,
	     int node3,
	     int node4,
	     int node5,
	     int node6,
	     int node7,
	     int node8,
	     NDMaterial &theMaterial,
	     double b1, double b2, double b3,
       Damping *damping)
  :Element(tag, ELE_TAG_Brick),
   connectedExternalNodes(8), applyLoad(0), load(0), Ki(0)
{
  B.Zero();

  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;

  connectedExternalNodes(4) = node5 ;
  connectedExternalNodes(5) = node6 ;
  connectedExternalNodes(6) = node7 ;
  connectedExternalNodes(7) = node8 ;

  for (int i=0; i<8; i++ ) {
      materialPointers[i] = theMaterial.getCopy("ThreeDimensional") ;
      if (materialPointers[i] == 0) {
	opserr << "Brick::constructor - failed to get a material of type: ThreeDimensional\n";
	exit(-1);
      } //end if
      nodePointers[i] = 0;
  } //end for i 

  // Body forces
  b[0] = b1;
  b[1] = b2;
  b[2] = b3;

  if (damping)
  {
    for (int i = 0; i < 8; i++)
    {
      theDamping[i] =(*damping).getCopy();
    
      if (!theDamping[i]) {
        opserr << "FourNodeQuad::FourNodeQuad -- failed to get copy of damping\n";
        exit(-1);
      }
    }
  }
  else
  {
    for (int i = 0; i < 8; i++) theDamping[i] = 0;
  }
}
//******************************************************************


//destructor 
Brick::~Brick( )
{

  for (int i=0 ; i<8; i++ ) {
    delete materialPointers[i] ;
  } //end for i

  if (load != 0)
    delete load;

  if (Ki != 0)
    delete Ki;

  for (int i = 0; i < 8; i++)
  {
    if (theDamping[i])
    {
      delete theDamping[i];
      theDamping[i] = 0;
    }
  }
  
}


//set domain
void  Brick::setDomain( Domain *theDomain ) 
{  

  int i ;

  //node pointers
  for ( i=0; i<8; i++ ) 
     nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;

    
  for (int i = 0; i < 8; i++)
  {
    if (theDamping[i] && theDamping[i]->setDomain(theDomain, 6)) {
      opserr << "Brick::setDomain -- Error initializing damping\n";
      return;
    }
  }

  this->DomainComponent::setDomain(theDomain);

}


int
Brick::setDamping(Domain *theDomain, Damping *damping)
{
  if (theDomain && damping)
  {
    for (int i = 0; i < 8; i++)
    {
      if (theDamping[i]) delete theDamping[i];

      theDamping[i] =(*damping).getCopy();
    
      if (!theDamping[i]) {
        opserr << "Brick::setDamping -- failed to get copy of damping\n";
        return -1;
      }
      if (theDamping[i] && theDamping[i]->setDomain(theDomain, 6)) {
        opserr << "Brick::setDamping -- Error initializing damping\n";
        return -2;
      }
    }
  }
  
  return 0;
}

//get the number of external nodes
int  Brick::getNumExternalNodes( ) const
{
  return 8 ;
} 
 

//return connected external nodes
const ID&  Brick::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 

Node **  
Brick::getNodePtrs(void) 
{
  return nodePointers ;
} 

//return number of dofs
int  Brick::getNumDOF( ) 
{
  return 24 ;
}


//commit state
int  Brick::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "Brick::commitState () - failed in base class";
  }    


  for (int i=0; i<8; i++ ) 
    success += materialPointers[i]->commitState( ) ;

  for (int i = 0; i < 8; i++ )
    if (theDamping[i]) success += theDamping[i]->commitState();

  
  return success ;
}
 


//revert to last commit 
int  Brick::revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i=0; i<8; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;

  for (int i = 0; i < 8; i++ )
    if (theDamping[i]) success += theDamping[i]->revertToLastCommit();
  
  return success ;
}
    

//revert to start 
int  Brick::revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i=0; i<8; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;

  for (int i = 0; i < 8; i++ )
    if (theDamping[i]) success += theDamping[i]->revertToStart();
  
  return success ;
}

//print out element data
void  Brick::Print(OPS_Stream &s, int flag)
{
    if (flag == 2) {
        
        s << "#Brick\n";
        
        int i;
        const int numNodes = 8;
        const int nstress = 6;

        for (i = 0; i < numNodes; i++) {
            const Vector &nodeCrd = nodePointers[i]->getCrds();
            const Vector &nodeDisp = nodePointers[i]->getDisp();
            s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << nodeCrd(2)
                << " " << nodeDisp(0) << " " << nodeDisp(1) << " " << nodeDisp(2) << endln;
        }
        
        // spit out the section location & invoke print on the scetion
        const int numMaterials = 8;
        
        static Vector avgStress(nstress);
        static Vector avgStrain(nstress);
        avgStress.Zero();
        avgStrain.Zero();
        for (i = 0; i < numMaterials; i++) {
            avgStress += materialPointers[i]->getStress();
            avgStrain += materialPointers[i]->getStrain();
        }
        avgStress /= numMaterials;
        avgStrain /= numMaterials;
        
        s << "#AVERAGE_STRESS ";
        for (i = 0; i < nstress; i++)
            s << avgStress(i) << " ";
        s << endln;
        
        s << "#AVERAGE_STRAIN ";
        for (i = 0; i < nstress; i++)
            s << avgStrain(i) << " ";
        s << endln;
        
        /*
        for (i=0; i<numMaterials; i++) {
          s << "#MATERIAL\n";
          //      materialPointers[i]->Print(s, flag);
          s << materialPointers[i]->getStress();
        }
        */
    }
    
    if (flag == OPS_PRINT_CURRENTSTATE) {
        
        s << "Standard Eight Node Brick \n";
        s << "Element Number: " << this->getTag() << endln;
        s << "Nodes: " << connectedExternalNodes;
        
        s << "Material Information : \n ";
        materialPointers[0]->Print(s, flag);
        
        s << endln;
        s << this->getTag() << " " << connectedExternalNodes(0)
            << " " << connectedExternalNodes(1)
            << " " << connectedExternalNodes(2)
            << " " << connectedExternalNodes(3)
            << " " << connectedExternalNodes(4)
            << " " << connectedExternalNodes(5)
            << " " << connectedExternalNodes(6)
            << " " << connectedExternalNodes(7)
            << endln;
        
        s << "Body Forces: " << b[0] << " " << b[1] << " " << b[2] << endln;
        s << "Resisting Force (no inertia): " << this->getResistingForce();
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"Brick\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
        for (int i = 1; i < 7; i++)
            s << connectedExternalNodes(i) << ", ";
        s << connectedExternalNodes(7) << "], ";
        s << "\"bodyForces\": [" << b[0] << ", " << b[1] << ", " << b[2] << "], ";
        s << "\"material\": \"" << materialPointers[0]->getTag() << "\"}";
    }
}
 
 
//return stiffness matrix 
const Matrix&  Brick::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    


//return secant matrix 
//const Matrix&  Brick::getSecantStiff( ) 

const Matrix&  Brick::getInitialStiff( ) 

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
    if(theDamping[i]) dd *= theDamping[i]->getStiffnessMultiplier();
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

	kk += ndf ;

      } // end for k loop

      jj += ndf ;

    } // end for j loop
  } //end for i gauss loop 

  Ki = new Matrix(stiff);

  return stiff ;
}    


//return mass matrix
const Matrix&  Brick::getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 



void  Brick::zeroLoad( )
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
Brick::addLoad(ElementalLoad *theLoad, double loadFactor)
{
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
    opserr << "Brick::addLoad() - ele with tag: " << this->getTag() << " does not deal with load type: " << type << "\n";
    return -1;
  }

  return -1;
}

int
Brick::addInertiaLoadToUnbalance(const Vector &accel)
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
const Vector&  Brick::getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  if (load != 0)
    resid -= *load;

  return resid ;   
}


//get residual with inertia terms
const Vector&  Brick::getResistingForceIncInertia( )
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

  return res;
}


//*********************************************************************
//form inertia terms

void   Brick::formInertiaTerms( int tangFlag ) 
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
int  
Brick::update(void) 
{

  //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 

  static const int ndm = 3 ;

  static const int ndf = 3 ; 

  static const int nstress = 6 ;
 
  static const int numberNodes = 8 ;

  static const int numberGauss = 8 ;

  static const int nShape = 4 ;

  int i, j, k, p, q ;
  int success ;
  
  static double volume ;

  static double xsj ;  // determinant jacaobian matrix 

  static double dvol[numberGauss] ; //volume element

  static double gaussPoint[ndm] ;

  static Vector strain(nstress) ;  //strain

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

  //---------B-matrices------------------------------------

    static Matrix BJ(nstress,ndf) ;      // B matrix node J

    static Matrix BJtran(ndf,nstress) ;

    static Matrix BK(nstress,ndf) ;      // B matrix node k

    static Matrix BJtranD(ndf,nstress) ;

  //-------------------------------------------------------

  
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


    //zero the strains
    strain.Zero( ) ;

    // j-node loop to compute strain 
    for ( j = 0; j < numberNodes; j++ )  {

      /**************** fmk - unwinding for performance
      //compute B matrix 
      BJ = computeB( j, shp ) ;

      //nodal displacements 
      const Vector &ul = nodePointers[j]->getTrialDisp( ) ;

      //compute the strain
      //strain += (BJ*ul) ; 
      strain.addMatrixVector(1.0,  BJ,ul,1.0 ) ;
      ***************************************************/


      //               | N,1      0     0    | 
      //   B       =   |   0     N,2    0    |
      //               |   0      0     N,3  |   (6x3)
      //               | N,2     N,1     0   |
      //               |   0     N,3    N,2  |
      //               | N,3      0     N,1  |

      //      B(0,0) = shp[0][node] ;
      //      B(1,1) = shp[1][node] ;
      //      B(2,2) = shp[2][node] ;
      //      B(3,0) = shp[1][node] ;
      //      B(3,1) = shp[0][node] ;
      //      B(4,1) = shp[2][node] ;
      //      B(4,2) = shp[1][node] ;
      //      B(5,0) = shp[2][node] ;
      //      B(5,2) = shp[0][node] ;

      double b00 = shp[0][j];
      double b11 = shp[1][j];
      double b22 = shp[2][j];
      double b30 = shp[1][j];
      double b31 = shp[0][j];
      double b41 = shp[2][j];
      double b42 = shp[1][j];
      double b50 = shp[2][j];
      double b52 = shp[0][j];

      const Vector &ul = nodePointers[j]->getTrialDisp();

      double ul0 = ul(0);
      double ul1 = ul(1);
      double ul2 = ul(2);

      strain(0) += b00 * ul0;
      strain(1) += b11 * ul1;
      strain(2) += b22 * ul2;
      strain(3) += b30 * ul0 + b31 * ul1;
      strain(4) += b41 * ul1 + b42 * ul2;
      strain(5) += b50 * ul0 + b52 * ul2;

    } // end for j
    
    //send the strain to the material 
    success = materialPointers[i]->setTrialStrain( strain ) ;

  } //end for i gauss loop 

  return 0;
}


//*********************************************************************
//form residual and tangent
void  Brick::formResidAndTangent( int tang_flag ) 
{

  //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 

  static const int ndm = 3 ;

  static const int ndf = 3 ; 

  static const int nstress = 6 ;
 
  static const int numberNodes = 8 ;

  static const int numberGauss = 8 ;

  static const int nShape = 4 ;

  int i, j, k, p, q ;


  static double volume ;

  static double xsj ;  // determinant jacaobian matrix 

  static double dvol[numberGauss] ; //volume element

  static double gaussPoint[ndm] ;

  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point

  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

  static Vector residJ(ndf) ; //nodeJ residual 

  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 

  static Vector stress(nstress) ;  //stress

  static Vector dampingStress(nstress) ;  //damping stress

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


    //compute the stress
    stress = materialPointers[i]->getStress( ) ;

    if (theDamping[i])
    {
      theDamping[i]->update(stress);
      dampingStress = theDamping[i]->getDampingForce();
      dampingStress *= dvol[i];
    }

    //multiply by volume element
    stress  *= dvol[i] ;

    if ( tang_flag == 1 ) {
      dd = materialPointers[i]->getTangent( ) ;
      if(theDamping[i]) dd *= theDamping[i]->getStiffnessMultiplier();
      dd *= dvol[i] ;
    } //end if tang_flag


    double stress0 = stress(0);
    double stress1 = stress(1);
    double stress2 = stress(2);
    double stress3 = stress(3);
    double stress4 = stress(4);
    double stress5 = stress(5);

    //residual and tangent calculations node loops

    int jj = 0 ;
    for ( j = 0; j < numberNodes; j++ ) {

      /* ************** fmk - unwinding for performance 
      ************************************************* */

      //               | N,1      0     0    | 
      //   B       =   |   0     N,2    0    |
      //               |   0      0     N,3  |   (6x3)
      //               | N,2     N,1     0   |
      //               |   0     N,3    N,2  |
      //               | N,3      0     N,1  |

      //      B(0,0) = shp[0][node] ;
      //      B(1,1) = shp[1][node] ;
      //      B(2,2) = shp[2][node] ;
      //      B(3,0) = shp[1][node] ;
      //      B(3,1) = shp[0][node] ;
      //      B(4,1) = shp[2][node] ;
      //      B(4,2) = shp[1][node] ;
      //      B(5,0) = shp[2][node] ;
      //      B(5,2) = shp[0][node] ;

      double b00 = shp[0][j];
      double b11 = shp[1][j];
      double b22 = shp[2][j];
      double b30 = shp[1][j];
      double b31 = shp[0][j];
      double b41 = shp[2][j];
      double b42 = shp[1][j];
      double b50 = shp[2][j];
      double b52 = shp[0][j];

      residJ(0) = b00 * stress0 + b30 * stress3 + b50 * stress5;
      residJ(1) = b11 * stress1 + b31 * stress3 + b41 * stress4;
      residJ(2) = b22 * stress2 + b42 * stress4 + b52 * stress5;

      residJ(0) += b00 * dampingStress[0] + b30 * dampingStress[3] + b50 * dampingStress[5];
      residJ(1) += b11 * dampingStress[1] + b31 * dampingStress[3] + b41 * dampingStress[4];
      residJ(2) += b22 * dampingStress[2] + b42 * dampingStress[4] + b52 * dampingStress[5];
      
      BJ = computeB( j, shp ) ;
   
      //transpose 
      //BJtran = transpose( nstress, ndf, BJ ) ;
      for (p=0; p<ndf; p++) {
	for (q=0; q<nstress; q++) 
	  BJtran(p,q) = BJ(q,p) ;
      }//end for p


      //residual 
      for ( p = 0; p < ndf; p++ ) {
        resid( jj + p ) += residJ(p)  ;
	if (applyLoad == 0)
	  resid( jj + p ) -= dvol[i]*b[p]*shp[3][j];
	else
	  resid( jj + p ) -= dvol[i]*appliedB[p]*shp[3][j];
      }

      if ( tang_flag == 1 ) {

	//BJtranD = BJtran * dd ;
	BJtranD.addMatrixProduct(0.0,  BJtran,dd,1.0) ;

	int kk = 0 ;
         for ( k = 0; k < numberNodes; k++ ) {

            BK = computeB( k, shp ) ;
  
 
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

void   Brick::computeBasis( ) 
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
Brick::computeB( int node, const double shp[4][8] )
{

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

Matrix  Brick::transpose( int dim1, 
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

int  Brick::sendSelf (int commitTag, Channel &theChannel)
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
  
  static ID idData(28);

  idData(24) = this->getTag();
  if (alphaM != 0 || betaK != 0 || betaK0 != 0 || betaKc != 0) 
    idData(25) = 1;
  else
    idData(25) = 0;
  
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

  idData(26) = 0;
  idData(27) = 0;
  if (theDamping[0]) {
    idData(26) = theDamping[0]->getClassTag();
    int dbTag = theDamping[0]->getDbTag();
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
        for (i = 0 ;  i < 8; i++)
	        theDamping[i]->setDbTag(dbTag);
	  }
    idData(27) = dbTag;
  }

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING Brick::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  static Vector dData(7);
  dData(0) = alphaM;
  dData(1) = betaK;
  dData(2) = betaK0;
  dData(3) = betaKc;
  dData(4) = b[0];
  dData(5) = b[1];
  dData(6) = b[2];

  if (theChannel.sendVector(dataTag, commitTag, dData) < 0) {
    opserr << "Brick::sendSelf() - failed to send double data\n";
    return -1;
  }    

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 8; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING Brick::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  // Ask the Damping to send itself
  if (theDamping[0]) {
    for (i = 0 ;  i < 8; i++) {
      res += theDamping[i]->sendSelf(commitTag, theChannel);
      if (res < 0) {
        opserr << "Brick::sendSelf -- could not send Damping\n";
        return res;
      }
    }
  }

  return res;

}
    
int  Brick::recvSelf (int commitTag, 
		      Channel &theChannel, 
		      FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(28);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING Brick::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(24));

  static Vector dData(7);
  if (theChannel.recvVector(dataTag, commitTag, dData) < 0) {
    opserr << "DispBeamColumn2d::sendSelf() - failed to recv double data\n";
    return -1;
  }    
  alphaM = dData(0);
  betaK = dData(1);
  betaK0 = dData(2);
  betaKc = dData(3);
  b[0] = dData(4);
  b[1] = dData(5);
  b[2] = dData(6);


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
	opserr << "Brick::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
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
    for (int i = 0; i < 8; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
	delete materialPointers[i];
	materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
	if (materialPointers[i] == 0) {
	  opserr << "Brick::recvSelf() - Broker could not create NDMaterial of class type " <<
	    matClassTag << endln;
	  exit(-1);
	}
	materialPointers[i]->setDbTag(matDbTag);
      }
      // Receive the material

      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "Brick::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }

  int dmpTag = (int)idData(26);
  if (dmpTag) {
    for (int i = 0 ;  i < 8; i++) {
      // Check if the Damping is null; if so, get a new one
      if (theDamping[i] == 0) {
        theDamping[i] = theBroker.getNewDamping(dmpTag);
        if (theDamping[i] == 0) {
          opserr << "Brick::recvSelf -- could not get a Damping\n";
          exit(-1);
        }
      }
  
      // Check that the Damping is of the right type; if not, delete
      // the current one and get a new one of the right type
      if (theDamping[i]->getClassTag() != dmpTag) {
        delete theDamping[i];
        theDamping[i] = theBroker.getNewDamping(dmpTag);
        if (theDamping[i] == 0) {
          opserr << "Brick::recvSelf -- could not get a Damping\n";
          exit(-1);
        }
      }
  
      // Now, receive the Damping
      theDamping[i]->setDbTag((int)idData(27));
      res += theDamping[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "Brick::recvSelf -- could not receive Damping\n";
        return res;
      }
    }
  }
  else {
    for (int i = 0; i < 8; i++)
    {
      if (theDamping[i])
      {
        delete theDamping[i];
        theDamping[i] = 0;
      }
    }
  }
    
  return res;
}
//**************************************************************************

int
Brick::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
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
    
    if (displayMode < 3 && displayMode > 0) {
        // determine the value to be drawn, the stress at nearest gauss point in displayMode dirn
        const Vector& stress1 = materialPointers[0]->getStress();
        const Vector& stress2 = materialPointers[1]->getStress();
        const Vector& stress3 = materialPointers[2]->getStress();
        const Vector& stress4 = materialPointers[3]->getStress();
        const Vector& stress5 = materialPointers[4]->getStress();
        const Vector& stress6 = materialPointers[5]->getStress();
        const Vector& stress7 = materialPointers[6]->getStress();
        const Vector& stress8 = materialPointers[7]->getStress();
        for (i = 0; i < 8; i++) {
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
    } else if (displayMode < 0)
        for (i = 0; i < 8; i++)
            values(i) = 0.0;

    return theViewer.drawCube(coords, values, this->getTag());
}

Response*
Brick::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","Brick");
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
      output.tag("ResponseType","sigma23");
      output.tag("ResponseType","sigma13");      

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }
    theResponse =  new ElementResponse(this, 3, Vector(48));

  } else if (strcmp(argv[0],"strains") ==0) {

    for (int i=0; i<8; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());

      output.tag("ResponseType","eps11");
      output.tag("ResponseType","eps22");
      output.tag("ResponseType","eps33");
      output.tag("ResponseType","eps12");
      output.tag("ResponseType","eps23");
      output.tag("ResponseType","eps13");      

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }
    theResponse =  new ElementResponse(this, 4, Vector(48));
    
  } else if (strcmp(argv[0],"dampingStresses") ==0) {

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
      output.tag("ResponseType","sigma23");
      output.tag("ResponseType","sigma13");      

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }
    theResponse =  new ElementResponse(this, 5, Vector(48));

  }

  
  output.endTag(); // ElementOutput
  return theResponse;
}

int 
Brick::getResponse(int responseID, Information &eleInfo)
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

  } else if (responseID == 4) {
    
    // Loop over the integration points
    int cnt = 0;
    for (int i = 0; i < 8; i++) {
      
      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStrain();
      stresses(cnt++) = sigma(0);
      stresses(cnt++) = sigma(1);
      stresses(cnt++) = sigma(2);
      stresses(cnt++) = sigma(3);
      stresses(cnt++) = sigma(4);
      stresses(cnt++) = sigma(5);
    }
    return eleInfo.setVector(stresses);
  }
  else if (responseID == 5) {
    
    // Loop over the integration points
    int cnt = 0;
    for (int i = 0; i < 8; i++) {
      
      // Get material stress response
      const Vector &sigma = theDamping[i]->getDampingForce();
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
Brick::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  if ((strstr(argv[0],"material") != 0) && (strcmp(argv[0],"materialState") != 0)) {

    if (argc < 3)
      return -1;

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 8)
      return materialPointers[pointNum-1]->setParameter(&argv[2], argc-2, param);
    else 
      return -1;
  }
  
  // otherwise it could be just a forall material parameter
  else {
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
Brick::updateParameter(int parameterID, Information &info)
{
    int res = -1;
	int matRes = res;

    if (parameterID == res) {
        return -1;
    } else {
        for (int i = 0; i<8; i++) {
            matRes = materialPointers[i]->updateParameter(parameterID, info);
        }
		if (matRes != -1) {
			res = matRes;
		}
		return res;
    }
}

