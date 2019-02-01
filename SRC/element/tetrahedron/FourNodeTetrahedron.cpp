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
// $Source: /usr/local/cvs/OpenSees/SRC/element/brick/FourNodeTetrahedron.cpp,v $

// Ed "C++" Love
//
// Eight node FourNodeTetrahedron element
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
#include <FourNodeTetrahedron.h>
// #include <shp3d.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <map>

void* OPS_FourNodeTetrahedron()
{
    if (OPS_GetNumRemainingInputArgs() < 6) 
    {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: element FourNodeTetrahedron eleTag? Node1? Node2? Node3? Node4? matTag?\n";
      return 0;
    }

    int idata[6];
    int num = 6;
    if (OPS_GetIntInput(&num,idata)<0) 
    {
      opserr<<"WARNING: invalid integer data\n";
      return 0;
    }

    NDMaterial* mat = OPS_getNDMaterial(idata[5]);
    if (mat == 0) 
    {
      opserr << "WARNING material not found\n";
      opserr << "material tag: " << idata[9];
      opserr << "\nFourNodeTetrahedron element: " << idata[0] << endln;
    }

    double data[3] = {0,0,0};
    num = OPS_GetNumRemainingInputArgs();

    if (num > 3) 
    {
      num = 3;
    }
    if (num > 0) 
    {
      if (OPS_GetDoubleInput(&num,data) < 0) 
      {
        opserr<<"WARNING: invalid double data\n";
        return 0;
      }     
    }
    return new FourNodeTetrahedron(idata[0],idata[1],idata[2],idata[3],idata[4],*mat,data[0],data[1],data[2]);
}

void* OPS_FourNodeTetrahedron(const ID& info)
{
    if (info.Size() == 0) {
	opserr << "WARNING: info is empty -- FourNodeTetrahedron\n";
	return 0;
    }

    // save data
    static std::map<int, Vector> meshdata;
    int idata[6];
    double data[3] = {0,0,0};
    if (info(0) == 1) {

	// check input
	if (info.Size() < 2) {
	    opserr << "WARNING: need info -- inmesh, meshtag\n";
	    return 0;
	}
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING insufficient arguments:\n";
	    opserr << "matTag <b1, b2, b3>\n";
	    return 0;
	}

	// get tag
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &idata[5]) < 0) {
	    opserr << "WARNING: failed to get material tag -- FourNodeTetrahedron\n";
	    return 0;
	}

	// get body forces
	numdata = OPS_GetNumRemainingInputArgs();
	if (numdata > 3) {
	    numdata = 3;
	}
	if (numdata > 0) {
	    if (OPS_GetDoubleInput(&numdata,data) < 0) {
		opserr << "WARNING: failed to get body force -- FourNodeTetrahedron\n";
		return 0;
	    }
	}

	// save data for a mesh
	Vector& mdata = meshdata[info(1)];
	mdata.resize(4);
	mdata(0) = (double)idata[5];
	for (int i=0; i<3; ++i) {
	    mdata(i+1) = data[i];
	}
	return &meshdata;
    }

    // load data
    if (info(0) == 2) {
	if (info.Size() < 7) {
	    opserr << "WARNING: need info -- inmesh, meshtag, eleTag, nd1, nd2, nd3, nd4\n";
	    return 0;
	}

	// get data for a mesh
	Vector& mdata = meshdata[info(1)];
	if (mdata.Size() < 4) {
	    return 0;
	}

	idata[5] = (int)mdata(0);
	for (int i=0; i<3; ++i) {
	    data[i] = mdata(i+1);
	}

	for (int i=2; i<7; ++i) {
	    idata[i-2] = info(i);
	}

	// check material
	NDMaterial* mat = OPS_getNDMaterial(idata[5]);
	if (mat == 0) {
	    opserr << "WARNING material not found\n";
	    opserr << "material tag: " << idata[5];
	    opserr << "\nFourNodeTetrahedron element: " << idata[0] << endln;
	}

	return new FourNodeTetrahedron(idata[0],idata[1],idata[2],idata[3],idata[4],*mat,data[0],data[1],data[2]);
    }

    return 0;
    
}

//static data
double  FourNodeTetrahedron::xl[3][NumNodes] ;
Matrix  FourNodeTetrahedron::stiff(NumDOFsTotal,NumDOFsTotal) ;
Vector  FourNodeTetrahedron::resid(NumDOFsTotal) ;
Matrix  FourNodeTetrahedron::mass(NumDOFsTotal,NumDOFsTotal) ;
Matrix  FourNodeTetrahedron::damping(NumDOFsTotal,NumDOFsTotal) ;

    
//quadrature data
const double  FourNodeTetrahedron::root3 = sqrt(3.0) ;
const double  FourNodeTetrahedron::one_over_root3 = 1.0 / root3 ;

const double  FourNodeTetrahedron::sg[] = { 0.25 } ;

const double  FourNodeTetrahedron::wg[] = { 0.16667 } ;

  
static Matrix B(NumStressComponents,NumDOFsPerNode) ;

//null constructor
FourNodeTetrahedron::FourNodeTetrahedron( ) 
:Element( 0, ELE_TAG_FourNodeTetrahedron ),
 connectedExternalNodes(NumNodes), applyLoad(0), load(0), Ki(0)
{
  B.Zero();

  for (int i=0; i < NumNodes; i++ ) {
    nodePointers[i] = 0;
  }

  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;
  
  materialPointers[0] = 0;

  for (int i = 0; i < NumNodes; ++i)
  {
    initDisp[i] = Vector(3);
    initDisp[i].Zero();
  }
  do_update = 1;
}


//*********************************************************************
//full constructor
FourNodeTetrahedron::FourNodeTetrahedron(int tag, 
       int node1,
       int node2,
       int node3,
       int node4,
       NDMaterial &theMaterial,
       double b1, double b2, double b3)
  :Element(tag, ELE_TAG_FourNodeTetrahedron),
   connectedExternalNodes(4), applyLoad(0), load(0), Ki(0)
{
  B.Zero();
  do_update = 1;
  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;


  for (int i=0; i<NumGaussPoints; i++ ) 
  {
      materialPointers[i] = theMaterial.getCopy("ThreeDimensional") ;
      if (materialPointers[i] == 0) 
      {
         opserr << "FourNodeTetrahedron::constructor - failed to get a material of type: ThreeDimensional\n";
         exit(-1);
      } //end if
      nodePointers[i] = 0;
  } //end for i 

  // Body forces
  b[0] = b1;
  b[1] = b2;
  b[2] = b3;

  for (int i = 0; i < NumNodes; ++i)
  {
    initDisp[i] = Vector(3);
    initDisp[i].Zero();
  }
}
//******************************************************************


//destructor 
FourNodeTetrahedron::~FourNodeTetrahedron( )
{

  for (int i=0 ; i<NumGaussPoints; i++ ) {
    delete materialPointers[i] ;
  } //end for i

  if (load != 0)
    delete load;

  if (Ki != 0)
    delete Ki;
  
}


//set domain
void  FourNodeTetrahedron::setDomain( Domain *theDomain ) 
{  

  int i ;

  //node pointers
  for ( i=0; i<NumNodes; i++ ) 
  {
      nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
      initDisp[i] = nodePointers[i]->getDisp();
  }

  this->DomainComponent::setDomain(theDomain);

}


//get the number of external nodes
int  FourNodeTetrahedron::getNumExternalNodes( ) const
{
  return NumNodes ;
} 
 

//return connected external nodes
const ID&  FourNodeTetrahedron::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 

Node **  
FourNodeTetrahedron::getNodePtrs(void) 
{
  return nodePointers ;
} 

//return number of dofs
int  FourNodeTetrahedron::getNumDOF( ) 
{
  return NumDOFsTotal ;
}


//commit state
int  FourNodeTetrahedron::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "FourNodeTetrahedron::commitState () - failed in base class";
  }    


  for (int i=0; i<NumGaussPoints; i++ ) 
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}
 


//revert to last commit 
int  FourNodeTetrahedron::revertToLastCommit( ) 
{
  int i ;
  int success = 0 ;

  for ( i=0; i<NumGaussPoints; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}
    

//revert to start 
int  FourNodeTetrahedron::revertToStart( ) 
{
  int i ;
  int success = 0 ;

  for ( i=0; i<NumGaussPoints; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

//print out element data
void  FourNodeTetrahedron::Print(OPS_Stream &s, int flag)
{
    if (flag == 2) {
        
        s << "#FourNodeTetrahedron\n";
        
        int i;
        const int numNodes = NumNodes;
        const int nstress = NumStressComponents;

        for (i = 0; i < numNodes; i++) {
            const Vector &nodeCrd = nodePointers[i]->getCrds();
            const Vector &nodeDisp = nodePointers[i]->getDisp();
            s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << nodeCrd(2)
                << " " << nodeDisp(0) << " " << nodeDisp(1) << " " << nodeDisp(2) << endln;
        }
        
        // spit out the section location & invoke print on the scetion
        const int numMaterials = 1;
        
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
        
        s << "Standard FourNodeTetrahedron \n";
        s << "Element Number: " << this->getTag() << endln;
        s << "Nodes: " << connectedExternalNodes;
        
        s << "Material Information : \n ";
        materialPointers[0]->Print(s, flag);
        
        s << endln;
        
        // s << "DEBUGME!" << endln;
        
        s << "Body Forces: " << b[0] << " " << b[1] << " " << b[2] << endln;
        
        // s << "DEBUGME!" << endln;

        s << "Resisting Force (no inertia): " << this->getResistingForce();
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"FourNodeTetrahedron\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
        for (int i = 1; i < 2; i++)
            s << connectedExternalNodes(i) << ", ";
        s << connectedExternalNodes(3) << "], ";
        s << "\"bodyForces\": [" << b[0] << ", " << b[1] << ", " << b[2] << "], ";
        s << "\"material\": \"" << materialPointers[0]->getTag() << "\"}";
    }
}
 
 
//return stiffness matrix 
const Matrix&  FourNodeTetrahedron::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    


//return secant matrix 
//const Matrix&  FourNodeTetrahedron::getSecantStiff( ) 

const Matrix&  FourNodeTetrahedron::getInitialStiff( ) 

{
  if (Ki != 0)
    return *Ki;

  //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 
  static const int ndm = 3 ;
  static const int ndf = NumDOFsPerNode ; 
  static const int nstress = NumStressComponents ;
  static const int numberNodes = NumNodes ;
  static const int numberGauss = NumGaussPoints ;
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

  // for ( i = 0; i < NumGaussPointsm; i++ ) 
  {
    // for ( j = 0; j < NumGaussPointsm; j++ ) 
    {
      // for ( k = 0; k < NumGaussPointsm; k++ ) 
      {

        i = j = k = 0; // Just one Gauss point in a tet

        gaussPoint[0] = sg[i] ;        
        gaussPoint[1] = sg[j] ;        
        gaussPoint[2] = sg[k] ;

        //get shape functions    
        shp3d( gaussPoint, xsj, shp, xl ) ;

        //save shape functions
        for ( p = 0; p < nShape; p++ ) 
        {
          for ( q = 0; q < numberNodes; q++ )
          {
            Shape[p][q][count] = shp[p][q] ;
          }
         } // end for p

        //volume element to also be saved
        dvol[count] = wg[count] * xsj ;  

        //volume += dvol[count] ;

        count++ ;
      } //end for k
    } //end for j
  } // end for i 
  

  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) 
  {

    //extract shape functions from saved array
    for ( p = 0; p < nShape; p++ ) 
    {
       for ( q = 0; q < numberNodes; q++ )
       {
         shp[p][q]  = Shape[p][q][i] ;
       }
    } // end for p


    dd = materialPointers[i]->getInitialTangent( ) ;
    dd *= dvol[i] ;
    
    jj = 0;
    for ( j = 0; j < numberNodes; j++ ) 
    {

      BJ = computeB( j, shp ) ;
   
      //transpose 
      //BJtran = transpose( nstress, ndf, BJ ) ;
      for (p=0; p<ndf; p++) 
      {
        for (q=0; q<nstress; q++) 
        {
          BJtran(p,q) = BJ(q,p) ;
        }
      }//end for p

      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixProduct(0.0,  BJtran, dd, 1.0) ;
      
      kk = 0 ;
      for ( k = 0; k < numberNodes; k++ ) 
      {
  
         BK = computeB( k, shp ) ;
  
  
        //stiffJK =  BJtranD * BK  ;
        stiffJK.addMatrixProduct(0.0,  BJtranD, BK, 1.0) ;
        
        for ( p = 0; p < ndf; p++ )  
        {
          for ( q = 0; q < ndf; q++ )
          {
            stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
          }
        } //end for p
       kk += ndf ;
      } // end for k loop

      jj += ndf ;
    } // end for j loop
  } //end for i gauss loop 


  // opserr << "STIFF = " << stiff << endln;

  Ki = new Matrix(stiff);

  return stiff ;
}    


//return mass matrix
const Matrix&  FourNodeTetrahedron::getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 



void  FourNodeTetrahedron::zeroLoad( )
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
FourNodeTetrahedron::addLoad(ElementalLoad *theLoad, double loadFactor)
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
      // if( loadFactor > 0)
      // {
      //     opserr << "loadfactor = " << loadFactor << endln;
      //       opserr << "      data = " << data;
      //       opserr << "      b    = " << b[0] << " " << b[1] << " " << b[2] << "\n"   ;
      //       opserr << "      appliedB    = " << appliedB[0] << " " << appliedB[1] << " " << appliedB[2] << "\n"   ;
      //     }
      return 0;
  } else {
    opserr << "FourNodeTetrahedron::addLoad() - ele with tag: " << this->getTag() << " does not deal with load type: " << type << "\n";
    return -1;
  }

  return -1;
}

int
FourNodeTetrahedron::addInertiaLoadToUnbalance(const Vector &accel)
{
  static const int numberNodes = 4 ;
  static const int numberGauss = 1 ;
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
  for (i=0; i<numberNodes; i++) 
  {
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
const Vector&  FourNodeTetrahedron::getResistingForce( ) 
{
  int tang_flag = 0 ; //don't get the tangent

  formResidAndTangent( tang_flag ) ;

  if (load != 0)
    resid -= *load;

  return resid ;   
}


//get residual with inertia terms
const Vector&  FourNodeTetrahedron::getResistingForceIncInertia( )
{
  static Vector res(12); res.Zero();

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

void   FourNodeTetrahedron::formInertiaTerms( int tangFlag ) 
{

  static const int ndm = 3 ;

  static const int ndf = NumDOFsPerNode ; 

  static const int numberNodes = NumNodes ;

  static const int numberGauss = NumGaussPoints ;

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

  if(do_update == 0)
  {
    return ;
  }


  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

  //gauss loop to compute and save shape functions 

  int count = 0 ;

  // for ( i = 0; i < 2; i++ ) 
  {
    // for ( j = 0; j < 2; j++ ) 
    {
      // for ( k = 0; k < 2; k++ ) 
      {
        i = j = k = 0;

        gaussPoint[0] = sg[i] ;        
        gaussPoint[1] = sg[j] ;        
        gaussPoint[2] = sg[k] ;

        //get shape functions    
        shp3d( gaussPoint, xsj, shp, xl ) ;

        //save shape functions
        for ( p = 0; p < nShape; p++ ) 
        {
          for ( q = 0; q < numberNodes; q++ )
          {
            Shape[p][q][count] = shp[p][q] ;
          }
        } // end for p

      //volume element to also be saved
      dvol[count] = wg[count] * xsj ;  

      count++ ;

      } //end for k
    } //end for j
  } // end for i 
  


  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) 
  {

    //extract shape functions from saved array
    for ( p = 0; p < nShape; p++ ) 
    {
        for ( q = 0; q < numberNodes; q++ )
        {
          shp[p][q]  = Shape[p][q][i] ;
        }
    } // end for p


    //node loop to compute acceleration
    momentum.Zero( ) ;
    for ( j = 0; j < numberNodes; j++ ) 
    {
      momentum.addVector( 1.0, nodePointers[j]->getTrialAccel(), shp[massIndex][j] ) ;
    }


    //density
    rho = materialPointers[i]->getRho() ;


    //multiply acceleration by density to form momentum
    momentum *= rho ;


    //residual and tangent calculations node loops
    jj = 0 ;
    for ( j = 0; j < numberNodes; j++ ) 
    {

      temp = shp[massIndex][j] * dvol[i] ;

      for ( p = 0; p < ndf; p++ )
      {
        resid( jj+p ) += ( temp * momentum(p) )  ;
      }

      if ( tangFlag == 1 ) 
      {

         //multiply by density
         temp *= rho ;

         //node-node mass
         kk = 0 ;
         for ( k = 0; k < numberNodes; k++ ) 
         {
            massJK = temp * shp[massIndex][k] ;
            for ( p = 0; p < ndf; p++ )  
            {
              mass( jj+p, kk+p ) += massJK ;
            }
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
FourNodeTetrahedron::update(void) 
{

  if(do_update == 0)
  {
    stiff.Zero();
    resid.Zero();
    mass.Zero();
    // damping.Zero();
    return 0;
  }


  //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 

  static const int ndm = 3 ;

  static const int ndf = NumDOFsPerNode ; 

  static const int nstress = NumStressComponents ;
 
  static const int numberNodes = NumNodes ;

  static const int numberGauss = NumGaussPoints ;

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

  // for ( i = 0; i < 2; i++ ) 
  {
    // for ( j = 0; j < 2; j++ ) 
    {
      // for ( k = 0; k < 2; k++ ) 
      {
        i = j = k = 0;

        gaussPoint[0] = sg[i] ;        
        gaussPoint[1] = sg[j] ;        
        gaussPoint[2] = sg[k] ;

        //get shape functions    
        shp3d( gaussPoint, xsj, shp, xl ) ;

        //save shape functions
        for ( p = 0; p < nShape; p++ ) 
        {
          for ( q = 0; q < numberNodes; q++ )
          {
            Shape[p][q][count] = shp[p][q] ;
          }
        } // end for p

        //volume element to also be saved
        dvol[count] = wg[count] * xsj ;  
        count++ ;
      } //end for k
    } //end for j
  } // end for i 
  

  //gauss loop 
  for ( i = 0; i < numberGauss; i++ ) 
  {

    //extract shape functions from saved array
    for ( p = 0; p < nShape; p++ ) 
    {
      for ( q = 0; q < numberNodes; q++ )
      {
        shp[p][q]  = Shape[p][q][i] ;
      }
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

      const Vector &ul = nodePointers[j]->getTrialDisp() - initDisp[j];

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

    // opserr << "8NB> strain = " << strain << endln;

  } //end for i gauss loop 

  return 0;
}


//*********************************************************************
//form residual and tangent
void  FourNodeTetrahedron::formResidAndTangent( int tang_flag ) 
{

  //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 

  static const int ndm = 3 ;

  static const int ndf = NumDOFsPerNode ; 

  static const int nstress = NumStressComponents ;
 
  static const int numberNodes = NumNodes ;

  static const int numberGauss = NumGaussPoints ;

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

  static Matrix dd(nstress,nstress) ;  //material tangent


  //---------B-matrices------------------------------------

    static Matrix BJ(nstress,ndf) ;      // B matrix node J

    static Matrix BJtran(ndf,nstress) ;

    static Matrix BK(nstress,ndf) ;      // B matrix node k

    static Matrix BJtranD(ndf,nstress) ;

  //-------------------------------------------------------

  // opserr << "DEBUGME!" << endln;

  //zero stiffness and residual 
  stiff.Zero( ) ;
  resid.Zero( ) ;

  if (do_update == 0)
  {
    return ;
  }

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

  //gauss loop to compute and save shape functions 

  // opserr << "DEBUGME!" << endln;

  int count = 0 ;
  volume = 0.0 ;

  // for ( i = 0; i < 2; i++ ) 
  {
    // for ( j = 0; j < 2; j++ ) 
    {
      // for ( k = 0; k < 2; k++ ) 
      {
        i = j = k = 0;

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
  for ( i = 0; i < numberGauss; i++ ) 
  {

    //extract shape functions from saved array
    for ( p = 0; p < nShape; p++ ) 
    {
      for ( q = 0; q < numberNodes; q++ )
      {
        shp[p][q]  = Shape[p][q][i] ;
      }
    } // end for p


    //compute the stress
    stress = materialPointers[i]->getStress( ) ;


    //multiply by volume element
    stress  *= dvol[i] ;

    if ( tang_flag == 1 ) {
      dd = materialPointers[i]->getTangent( ) ;
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
    for ( j = 0; j < numberNodes; j++ ) 
    {

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
      
      BJ = computeB( j, shp ) ;
   
      //transpose 
      //BJtran = transpose( nstress, ndf, BJ ) ;
      for (p=0; p<ndf; p++) 
      {
        for (q=0; q<nstress; q++) 
        {
          BJtran(p,q) = BJ(q,p) ;
        }
      }//end for p


      //residual 
      for ( p = 0; p < ndf; p++ ) 
      {
        resid( jj + p ) += residJ(p)  ;
        if (applyLoad == 0)
        {
          // resid( jj + p ) -= dvol[i]*b[p]*shp[3][j];
        }
        else
        {
          resid( jj + p ) -= dvol[i]*appliedB[p]*shp[3][j];
        }
      }

      if ( tang_flag == 1 ) 
      {

        //BJtranD = BJtran * dd ;
        BJtranD.addMatrixProduct(0.0,  BJtran, dd,1.0) ;
        // opserr << "DEBUUG!! BJtranD = " << BJtranD << endln;
        // opserr << "DEBUUG!! dd = " << dd << endln;

        int kk = 0 ;
        for ( k = 0; k < numberNodes; k++ ) 
        {
          // opserr << "DEBUUG!! k = " << k << endln;

          BK = computeB( k, shp ) ;
          stiffJK.addMatrixProduct(0.0,  BJtranD, BK,1.0) ;

          // opserr << "DEBUUG!! BK = " << BK << endln;
          // opserr << "DEBUUG!! stiffJK = " << stiffJK << endln;

          for ( p = 0; p < ndf; p++ )  
          {
            for ( q = 0; q < ndf; q++ )
            {
              stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
            }
          } //end for p
          kk += ndf ;
        } // end for k loop
        // opserr << "STIFF = " << stiff << endln;
      } // end if tang_flag 
      jj += ndf ;
    } // end for j loop
  } //end for i gauss loop 

  
  


  return ;
}


//************************************************************************
//compute local coordinates and basis

void   FourNodeTetrahedron::computeBasis( ) 
{
  //nodal coordinates 
  int i ;
  for ( i = 0; i < 4; i++ ) {
       const Vector &coorI = nodePointers[i]->getCrds( ) ;
       xl[0][i] = coorI(0) ;
       xl[1][i] = coorI(1) ;
       xl[2][i] = coorI(2) ;
  }  //end for i 

}

//*************************************************************************
//compute B

const Matrix&   
FourNodeTetrahedron::computeB( int node, const double shp[4][4] )
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

Matrix  FourNodeTetrahedron::transpose( int dim1, int dim2,  const Matrix &M ) 
{
  int i, j ;

  Matrix Mtran( dim2, dim1 ) ;

  for ( i = 0; i < dim1; i++ ) {
     for ( j = 0; j < dim2; j++ ) 
         Mtran(j,i) = M(i,j) ;
  } // end for i

  return Mtran ;
}

//**********************************************************************

int  FourNodeTetrahedron::sendSelf (int commitTag, Channel &theChannel)
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
  
  static ID idData(27);

  idData(24) = this->getTag();
  if (alphaM != 0 || betaK != 0 || betaK0 != 0 || betaKc != 0) 
    idData(25) = 1;
  else
    idData(25) = 0;
  
  int i;
  for (i = 0; i < NumGaussPoints; i++) 
  {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) 
    {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
      {
        materialPointers[i]->setDbTag(matDbTag);
      }
    }
    idData(i+8) = matDbTag;
  }
  
  idData(16) = connectedExternalNodes(0);
  idData(17) = connectedExternalNodes(1);
  idData(18) = connectedExternalNodes(2);
  idData(19) = connectedExternalNodes(3);
  idData(26) = do_update;
  // idData(20) = connectedExternalNodes(4);
  // idData(21) = connectedExternalNodes(5);
  // idData(22) = connectedExternalNodes(6);
  // idData(23) = connectedExternalNodes(7);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING FourNodeTetrahedron::sendSelf() - " << this->getTag() << " failed to send ID\n";
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
    opserr << "FourNodeTetrahedron::sendSelf() - failed to send double data\n";
    return -1;
  }    

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < NumGaussPoints; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) 
    {
      opserr << "WARNING FourNodeTetrahedron::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;

}
    
int  FourNodeTetrahedron::recvSelf (int commitTag, 
          Channel &theChannel, 
          FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(27);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING FourNodeTetrahedron::recvSelf() - " << this->getTag() << " failed to receive ID\n";
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
  do_update = idData(26);
  // connectedExternalNodes(4) = idData(20);
  // connectedExternalNodes(5) = idData(21);
  // connectedExternalNodes(6) = idData(22);
  // connectedExternalNodes(7) = idData(23);


  if (materialPointers[0] == 0) 
  {
    for (int i = 0; i < NumGaussPoints; i++) 
    {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);

      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      
      if (materialPointers[i] == 0) 
      {
        opserr << "FourNodeTetrahedron::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
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
  else 
  {
    for (int i = 0; i < NumGaussPoints; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) 
      {
          delete materialPointers[i];
          materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
          if (materialPointers[i] == 0) 
          {
            opserr << "FourNodeTetrahedron::recvSelf() - Broker could not create NDMaterial of class type " <<
              matClassTag << endln;
            exit(-1);
          }
          materialPointers[i]->setDbTag(matDbTag);
      }
      // Receive the material

      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) 
      {
        opserr << "FourNodeTetrahedron::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }

  return res;
}
//**************************************************************************

int
FourNodeTetrahedron::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    return -1;
}

Response*
FourNodeTetrahedron::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","FourNodeTetrahedron");
  output.attr("eleTag",this->getTag());
  for (int i=1; i<=4; i++) 
  {
    sprintf(outputData,"node%d",i);
    output.attr(outputData,nodePointers[i-1]->getTag());
  }


  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
  {
    for (int i=1; i<=3; i++) 
    {
      sprintf(outputData,"P1_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P2_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P3_%d",i);
      output.tag("ResponseType",outputData);
    }

    theResponse = new ElementResponse(this, 1, resid);
  }   
  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) 
  {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 1) 
    {
      output.tag("GaussPoint");
      output.attr("number",pointNum);
      theResponse =  materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, output);
      output.endTag(); // GaussPoint
    }


  } 
  else if (strcmp(argv[0],"stresses") ==0) 
  {
    for (int i=0; i<1; i++) 
    {
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
    theResponse =  new ElementResponse(this, 3, Vector(6));

  } 
  else if (strcmp(argv[0],"strains") ==0) 
  {
    for (int i=0; i<1; i++) 
    {
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
    theResponse =  new ElementResponse(this, 4, Vector(6));
  }
  output.endTag(); // ElementOutput
  return theResponse;
}

int 
FourNodeTetrahedron::getResponse(int responseID, Information &eleInfo)
{
  static Vector stresses(48);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2)
    return eleInfo.setMatrix(this->getTangentStiff());
    
  else if (responseID == 3) {
    
    // Loop over the integration points
    int cnt = 0;
    for (int i = 0; i < 1; i++) {
      
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
    for (int i = 0; i < 1; i++) {
      
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

  else
    return -1;
}

int
FourNodeTetrahedron::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  if ((strstr(argv[0],"material") != 0) && (strcmp(argv[0],"materialState") != 0)) 
  {
    if (argc < 3)
    {
      return -1;
    }

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 1)
    {
      return materialPointers[pointNum-1]->setParameter(&argv[2], argc-2, param);
    }
    else 
    {
      return -1;
    }
  }
  else if ((strstr(argv[0],"setDispInit") != 0) && (strcmp(argv[0],"setdispinit") == 0)) 
  {
    return param.addObject(1313, this);
  }
  else if ( strcmp(argv[0],"update") == 0 )
  {
    return param.addObject(1414, this);
  }
  else // otherwise it could be just a forall material parameter
  {
    int matRes = res;
    for (int i=0; i<1; i++) 
    {
      matRes =  materialPointers[i]->setParameter(argv, argc, param);
      if (matRes != -1)
      {
        res = matRes;
      }
    }
  }
  return res;
}
    
int
FourNodeTetrahedron::updateParameter(int parameterID, Information &info)
{
    int res = -1;
    int matRes = res;

    if (parameterID == res) 
    {
        return -1;
    } 
    else if (parameterID == 1313)
    {
      int doit = info.theDouble;
      if (doit == 1)
      {
        Domain * mydomain = this->getDomain();
        opserr << "FourNodeTetrahedron::updateParameter - ele tag = " << this->getTag()  << " - sets init disp ";
        for ( int i = 0; i < NumNodes; i++ ) 
        {
            nodePointers[i] = mydomain->getNode( connectedExternalNodes(i) ) ;
            initDisp[i] = nodePointers[i]->getDisp();
            opserr << " (" << initDisp[i](0) << " " << initDisp[i](1) << " " << initDisp[i](1) << ") ";
        }
        opserr << endln;
      }
      return 0;
    }
    else if (parameterID == 1414)
    {
      int new_do_update = info.theDouble;
      if (do_update == 0 & new_do_update == 1)
      {
        do_update = 1;
        Domain * mydomain = this->getDomain();
        opserr << "4Ntet::updateParameter - ele tag = " << this->getTag()  << " - sets to update and init disp ";
        for ( int i = 0; i < NumNodes; i++ ) 
        {
            nodePointers[i] = mydomain->getNode( connectedExternalNodes(i) ) ;
            initDisp[i] = nodePointers[i]->getDisp();
            opserr << " (" << initDisp[i](0) << " " << initDisp[i](1) << " " << initDisp[i](1) << ") ";
        }
        opserr << endln;
      }
      if(new_do_update == 0)
      {
        opserr << "4Ntet::updateParameter - ele tag = " << this->getTag()  << " - will not update\n";
      }
      do_update = new_do_update;
      return 0;
    }
    else 
    {
      for (int i = 0; i<1; i++) 
      {
        matRes = materialPointers[i]->updateParameter(parameterID, info);
      }
      if (matRes != -1) {
        res = matRes;
      }
      return res;
    }
}

/*
      Inputs:
         ss[4]     - Natural coordinates of point
         xl[3][4]  - Nodal coordinates for element

      Outputs:
         xsj        - Jacobian determinant at point
         shp[4][4]  - Shape functions and derivatives at point
                     shp[0][i] = dN_i/dx
                     shp[1][i] = dN_i/dy
                     shp[2][i] = dN_i/dzc
                     shp[3][i] =  N_i
*/
void  
FourNodeTetrahedron::shp3d( const double ss[4], double &xsj, double shp[4][4], const double xl[3][4]   )
{
    // Mathematica formulation by Carlos Felippa.
    //
    // IsoTet4ShapeFunDer[xyztet_,numer_]:= Module[{
    // x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,
    // x12,x13,x14,x23,x24,x34,x21,x31,x32,x42,x43,
    // y12,y13,y14,y23,y24,y34,y21,y31,y32,y42,y43,
    // z12,z13,z14,z23,z24,z34,z21,z31,z32,z42,z43,
    // a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,Nfx,Nyy,Nfz,Jdet},\
    // x12=x1-x2; x13=x1-x3; x14=x1-x4; x23=x2-x3; x24=x2-x4; x34=x3-x4;
    // x21=-x12; x31=-x13; x41=-x14; x32=-x23; x42=-x24; x43=-x34;
    // y12=y1-y2; y13=y1-y3; y14=y1-y4; y23=y2-y3; y24=y2-y4; y34=y3-y4;
    // y21=-y12; y31=-y13; y41=-y14; y32=-y23; y42=-y24; y43=-y34;
    // z12=z1-z2; z13=z1-z3; z14=z1-z4; z23=z2-z3; z24=z2-z4; z34=z3-z4;
    // z21=-z12; z31=-z13; z41=-z14; z32=-z23; z42=-z24; z43=-z34;
    double x12 = xl[0][0] - xl[0][1];
    double x13 = xl[0][0] - xl[0][2];
    double x14 = xl[0][0] - xl[0][3];
    double x23 = xl[0][1] - xl[0][2];
    double x24 = xl[0][1] - xl[0][3];
    double x34 = xl[0][2] - xl[0][3];
  
    double x21=-x12;
    double x31=-x13;
    double x41=-x14;
    double x32=-x23;
    double x42=-x24;
    double x43=-x34;
    
    double y12 = xl[1][0]-xl[1][1];
    double y13 = xl[1][0]-xl[1][2];
    double y14 = xl[1][0]-xl[1][3];
    double y23 = xl[1][1]-xl[1][2];
    double y24 = xl[1][1]-xl[1][3];
    double y34 = xl[1][2]-xl[1][3];
  
    double y21=-y12;
    double y31=-y13;
    double y41=-y14;
    double y32=-y23;
    double y42=-y24;
    double y43=-y34;
  
    double z12 = xl[2][0]-xl[2][1];
    double z13 = xl[2][0]-xl[2][2];
    double z14 = xl[2][0]-xl[2][3];
    double z23 = xl[2][1]-xl[2][2];
    double z24 = xl[2][1]-xl[2][3];
    double z34 = xl[2][2]-xl[2][3];
  
    double z21=-z12;
    double z31=-z13;
    double z41=-z14;
    double z32=-z23;
    double z42=-z24;
    double z43=-z34;
    
    // Jdet=x21*(y23*z34-y34*z23)+x32*(y34*z12-y12*z34)+x43*(y12*z23-y23*z12);
    double Jdet=x21*(y23*z34-y34*z23)+x32*(y34*z12-y12*z34)+x43*(y12*z23-y23*z12);
 
    // a1=y42*z32-y32*z42; b1=x32*z42-x42*z32; c1=x42*y32-x32*y42;
    // a2=y31*z43-y34*z13; b2=x43*z31-x13*z34; c2=x31*y43-x34*y13;
    // a3=y24*z14-y14*z24; b3=x14*z24-x24*z14; c3=x24*y14-x14*y24;
    // a4=y13*z21-y12*z31; b4=x21*z13-x31*z12; c4=x13*y21-x12*y31;
    double a1=y42*z32-y32*z42; double b1=x32*z42-x42*z32; double c1=x42*y32-x32*y42;
    double a2=y31*z43-y34*z13; double b2=x43*z31-x13*z34; double c2=x31*y43-x34*y13;
    double a3=y24*z14-y14*z24; double b3=x14*z24-x24*z14; double c3=x24*y14-x14*y24;
    double a4=y13*z21-y12*z31; double b4=x21*z13-x31*z12; double c4=x13*y21-x12*y31;

    xsj=Jdet;

    // Nfx={a1,a2,a3,a4}; Nfy={b1,b2,b3,b4}; Nfz={c1,c2,c3,c4};
    shp[0][0] = a1/Jdet;
    shp[0][1] = a2/Jdet;
    shp[0][2] = a3/Jdet;
    shp[0][3] = a4/Jdet;

    shp[1][0] = b1/Jdet;
    shp[1][1] = b2/Jdet;
    shp[1][2] = b3/Jdet;
    shp[1][3] = b4/Jdet;

    shp[2][0] = c1/Jdet;
    shp[2][1] = c2/Jdet;
    shp[2][2] = c3/Jdet;
    shp[2][3] = c4/Jdet;

    shp[3][0] = ss[0];
    shp[3][1] = ss[1];
    shp[3][2] = ss[2];
    shp[3][3] = 1 - ss[0] - ss[1] - ss[2];


    // opserr << "ss = " << ss[0] << " "
    //   << ss[1] << " "
    //   << ss[2] << " "
    //   << shp[3][4] << endln;

    // opserr << "Jdet = " << Jdet << endln;

    // opserr << "shp = " <<  endln;
    // for(int i = 0; i <= 4; i++)
    // {
    //   for(int j = 0; j <= 4; j++)
    //   {
    //     opserr << shp[i][j] << " ";
    //   }
    //   opserr << endln;
    // }

    // Return[{Nfx,Nfy,Nfz,Jdet}]];
   return ;
}
