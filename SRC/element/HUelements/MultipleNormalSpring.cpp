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

// $Revision: 1.0 $
// $Date: 2013-06-04 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/multipleNormalSpring/MultipleNormalSpring.h,v $

// Written: Ken Ishii
// Created: June 2012
//
// Multiple Normal Spring (MNS) model
//
// Description: This file contains the implementation of the MultipleNormalSpring class.
//              This file contains the function to parse the TCL input
// element multipleNormalSpring eleTag? iNode? jNode? nDivide? -mat matTag? -shape shape? -size size? <-lambda lambda?> <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>
//

// This element works in 3-dim 6-dof model.
// The multiple normal springs are arranged parallel to local-x axis.
// Stiffness of local-y,z,rx remain zero.


#include <ID.h>
#include <Vector.h>

#include <MultipleNormalSpring.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <UniaxialMaterial.h>

#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>

//extern void printCommand(int argc, TCL_Char **argv);
extern "C" double dbesi0(double);
extern "C" double dbesi1(double);


static bool errDetected(bool ifNoError, const char *msg){
  if (ifNoError){
    opserr << "" << endln;
    opserr << "========================================" << endln;
    opserr << "multipleNormalSpring element : input error detected" << endln;
    opserr << "------------------------------" << endln;
  }
  opserr << "  " << msg << endln;
  return false;
};

void * OPS_ADD_RUNTIME_VPV(OPS_MultipleNormalSpring)
{
    // 3-dim, 6-dof
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    if (ndm != 3 || ndf != 6) {
	opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
	opserr << "WARNING multipleNormalSpring command only works when ndm is 3 and ndf is 6" << endln;
	return 0;
    }


    //arguments (necessary)
    int eleTag;
    int iNode;
    int jNode;
    int nDivide;

    //arguments (necessary, input with -???)
    int matTag;
    UniaxialMaterial *material;
    int shape;
    double size;

    //arguments (optional, input with -???)
    double lambda = -1.0;
    Vector oriX(0);
    Vector oriYp(3); oriYp(0) = 0.0; oriYp(1) = 1.0; oriYp(2) = 0.0;
    double mass = 0.0;

    // input comfirmation
    int recvMat    = 0;
    int recvShape  = 0;
    int recvSize   = 0;
    int recvLambda = 0;
    int recvOrient = 0;
    int recvMass   = 0;


    //
    Element *theElement = 0;


    // error flag
    bool ifNoError = true;



    if (OPS_GetNumRemainingInputArgs() < 4)  { //element multipleNormalSpring eleTag? iNode? jNode? nDivide?

	ifNoError = errDetected(ifNoError,"insufficient arguments");

    } else {

	//argv[2~5]
	int idata[4];
	int numdata = 4;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
	    ifNoError = errDetected(ifNoError,"invalid int inputs");
	}

	eleTag = idata[0];
	iNode = idata[1];
	jNode = idata[2];
	nDivide = idata[3];
	if (nDivide <= 0)  {
	    ifNoError = errDetected(ifNoError,"invalid nDivide");
	}

	//argv[6~]
	while (OPS_GetNumRemainingInputArgs() > 0) {
      
	    double value;
	    const char* flag = OPS_GetString();
	    if (strcmp(flag,"-mat")==0 && OPS_GetNumRemainingInputArgs()>0) { // -mat matTag?

		numdata = 1;
		if (OPS_GetIntInput(&numdata, &matTag) < 0) {
		    ifNoError = errDetected(ifNoError,"invalid matTag");
		}


		material = OPS_getUniaxialMaterial(matTag);
		if (material == 0)  {
		    ifNoError = errDetected(ifNoError,"material model not found");
		}

		recvMat++ ;

	    } else if (strcmp(flag,"-shape")==0 && OPS_GetNumRemainingInputArgs()>0) { // -shape shape?
		const char* shapeflag = OPS_GetString();
		if (strcmp(shapeflag,"round") == 0) {
		    shape = 1; //round shape
		} else if (strcmp(shapeflag,"square") == 0) {
		    shape = 2; //square
		} else {
		    ifNoError = errDetected(ifNoError,"invalid shape (\"round\" or \"square\" are available)");
		}
	
		recvShape++ ;
	
	    } else if (strcmp(flag,"-size")==0 && OPS_GetNumRemainingInputArgs()>0) { // -size size?

		numdata = 1;
		if (OPS_GetDoubleInput(&numdata, &size)<0 || size<=0) {
		    ifNoError = errDetected(ifNoError,"invalid size");
		}

		recvSize++ ;

	    } else if (strcmp(flag,"-lambda")==0 && OPS_GetNumRemainingInputArgs()>0) { // <-lambda lambda?>
		numdata = 1;
		if (OPS_GetDoubleInput(&numdata, &lambda)<0 || lambda<0) {
		    ifNoError = errDetected(ifNoError,"invalid lambda");
		}

		recvLambda++ ;

	    } else if (strcmp(flag,"-orient")==0 && OPS_GetNumRemainingInputArgs()>=6) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

		oriX.resize(3);

		for (int j=1; j<=3; j++) {
		    numdata = 1;
		    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
			ifNoError = errDetected(ifNoError,"invalid orient");
		    } else {
			oriX(j-1) = value;
		    }
		}

		for (int j=1; j<=3; j++) {
		    numdata = 1;
		    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
			ifNoError = errDetected(ifNoError,"invalid orient");
		    } else {
			oriYp(j-1) = value;
		    }
		}

		recvOrient++ ;

	    } else if (strcmp(flag,"-orient")==0 && OPS_GetNumRemainingInputArgs()>=3) { // <-orient yp1? yp2? yp3?>

		for (int j=1; j<=3; j++) {
		    numdata = 1;
		    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
			ifNoError = errDetected(ifNoError,"invalid orient");
		    } else {
			oriYp(j-1) = value;
		    }
		}

		recvOrient++ ;

	    } else if (strcmp(flag,"-mass")==0 && OPS_GetNumRemainingInputArgs()>0) { // <-mass m?> �̓ǂݍ���
		numdata = 1;
		if (OPS_GetDoubleInput(&numdata, &mass)<0 || mass<=0) {
		    ifNoError = errDetected(ifNoError,"invalid mass");
		}

		recvMass++ ;
	
	    } else { //invalid option
	
		ifNoError = errDetected(ifNoError,"invalid optional arguments");
		break;

	    }

	}
    
    } // end input
  


    // input cofirmation
    // necessary arguments
    if (recvMat != 1)  {
	char buf[100];
	sprintf(buf,"wrong number of -mat inputs (got %d inputs, but want 1 input)",recvMat);
	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvShape != 1)  {
	char buf[100];
	sprintf(buf,"wrong number of -shape inputs (got %d inputs, but want 1 input)",recvShape);
	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvSize != 1)  {
	char buf[100];
	sprintf(buf,"wrong number of -size inputs (got %d inputs, but want 1 input)",recvSize);
	ifNoError = errDetected(ifNoError,buf);
    }

    // optional arguments
    if (recvLambda >= 2)  {
	char buf[100];
	sprintf(buf,"wrong number of -lambda inputs (got %d inputs, but want 1 input)",recvLambda);
	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvOrient >= 2)  {
	char buf[100];
	sprintf(buf,"wrong number of -ori inputs (got %d inputs, but want 1 input)",recvOrient);
	ifNoError = errDetected(ifNoError,buf);
    }

    if (recvMass >= 2)  {
	char buf[100];
	sprintf(buf,"wrong number of -mass inputs (got %d inputs, but want 1 input)",recvMass);
	ifNoError = errDetected(ifNoError,buf);
    }

  
    // if error detected
    if (!ifNoError) {
	opserr << "------------------------------" << endln;
	//input:
	// printCommand(argc, argv);
	//want:
	opserr << "Want: element multipleNormalSpring eleTag? iNode? jNode? nDivide? -mat matTag? -shape shape? -size size? <-lambda lambda?> <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>\n";
	opserr << "========================================" << endln;
	opserr << "" << endln;
	return 0;
    }
 
  

    // now create the multipleNormalSpring
    theElement = new MultipleNormalSpring(eleTag, iNode, jNode, nDivide, material, shape, size, lambda, oriYp, oriX, mass);

    return theElement;
}



// initialize the class wide variables
Matrix MultipleNormalSpring::theMatrix(12,12);
Vector MultipleNormalSpring::theVector(12);
Vector MultipleNormalSpring::theLoad(12);


MultipleNormalSpring::MultipleNormalSpring(int Tag, int Nd1, int Nd2,
					   int NDivide,
					   UniaxialMaterial *Material,
					   int Shape,
					   double Size,
					   double Lambda,
					   const Vector OriYp, const Vector OriX, double Mass)
  : Element(Tag, ELE_TAG_MultipleNormalSpring),
    connectedExternalNodes(2),
    nDivide(NDivide), shape(Shape),size(Size),lambda(Lambda),oriX(OriX), oriYp(OriYp), mass(Mass),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6)
{
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "MultipleNormalSpring::setUp() - element: "
	   << this->getTag() << " failed to create an ID of size 2\n";
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;
  
  // check material input
  if (Material == 0)  {
    opserr << "MultipleNormalSpring::MultipleNormalSpring() - "
	   << "null uniaxial material pointer passed.\n";
    exit(-1);
  }

  theMaterials = new UniaxialMaterial* [nDivide*nDivide]; // materials

  // get copies of the uniaxial materials
  for (int i=0; i<nDivide*nDivide; i++)  {
    theMaterials[i] = Material->getCopy();

    if (theMaterials[i] == 0) {
      opserr << "MultipleNormalSpring::MultipleNormalSpring() - "
 	     << "failed to copy uniaxial material.\n";
      exit(-1);
    }
  }
  

  //position of centroid, distribution factor
  posLy = new double [nDivide*nDivide]; //local-y position
  posLz = new double [nDivide*nDivide]; //local-z position
  distFct = new double [nDivide*nDivide]; //distribution factor
  if (shape == 1){//round shape

    incA = (M_PI*size*size)/(4.0*nDivide*nDivide); //area of each spring

    int p = nDivide%2;
    int k = -1; //index
  
    for(int i=1; i<=((nDivide+p)/2); i++) {
      //circle
      if(i==1 && p==1) {
	k++; //k=1
	double r = 0.0;
	posLy[k] = 0.0;
	posLz[k] = 0.0;
	
	//distribution factor
	if(lambda < 0) { //uniform
	  distFct[k] = 1.0;
	} else if(lambda == 0) { //parabolic
	  distFct[k] = 2.0;
	} else {
	  distFct[k] = (1.0-1.0/dbesi0(lambda))/(1.0-2.0/lambda*dbesi1(lambda)/dbesi0(lambda));
	}

	continue;
      }

      //circular sector
      double tht = M_PI/(4.0*(2*i-1-p));
      double r1 = ((2.0*i-2-p)/(2.0*nDivide))*size;
      double r2 = ((2.0*i-p)/(2.0*nDivide))*size;
      double r = (2.0/3.0)*(sin(tht)/tht)*((r1*r1+r1*r2+r2*r2)/(r1+r2)) ;
      for(int j=1; j<=(4*(2*i-1-p)); j++){
	k++; //k=(2*i-2-p)^2+j
	double alp = (2*j-1)*tht;
	posLy[k] = r*cos(alp);
	posLz[k] = r*sin(alp);

	//distribution factor
	if(lambda < 0) {
	  distFct[k] = 1.0; //uniform
	} else if(lambda == 0) {
	  distFct[k] = 2.0*(1.0-r/(size/2)*r/(size/2)); //parabolic
	} else {
	  distFct[k] = (1.0-dbesi0(lambda*r/(size/2))/dbesi0(lambda))/(1.0-2.0/lambda*dbesi1(lambda)/dbesi0(lambda));
	}
      }
    }

  } else {//square shape
    opserr << "MultipleNormalSpring::MultipleNormalSpring() - "
	   << "square shape \n";
    exit(-1);
  }


  // initial basic stiffness matrix
  basicStiffInit.Zero();
  for (int i=0; i<(nDivide*nDivide); i++)  {

    //total stiffness
    double tmpStiff = (theMaterials[i]->getInitialTangent())*incA*distFct[i]/hgt;

    basicStiffInit(0,0) += tmpStiff ; //x,x
    basicStiffInit(0,4) += tmpStiff * posLz[i] ; //x,ry
    basicStiffInit(0,5) += tmpStiff * posLy[i] ; //x,rz
    basicStiffInit(4,0) += tmpStiff * posLz[i] ; //ry,x
    basicStiffInit(4,4) += tmpStiff * posLz[i] * posLz[i] ; //ry,ry
    basicStiffInit(4,5) += tmpStiff * posLz[i] * posLy[i] ; //ry,rz
    basicStiffInit(5,0) += tmpStiff * posLy[i] ; //rz,x
    basicStiffInit(5,4) += tmpStiff * posLy[i] * posLz[i] ; //rz,ry
    basicStiffInit(5,5) += tmpStiff * posLy[i] * posLy[i] ; //rz,rz
  }
  
  // initialize variables
  this->revertToStart();
}


MultipleNormalSpring::MultipleNormalSpring()
  : Element(0, ELE_TAG_MultipleNormalSpring),
    connectedExternalNodes(2),
    nDivide(0), shape(0), size(0), oriX(0), oriYp(0), mass(0.0),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6)
{	


  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "MultipleNormalSpring::MultipleNormalSpring() - "
	   <<  "failed to create an ID of size 2\n";
    exit(-1);
  }
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  // set material array pointers to NULL
  theMaterials = 0;

}


MultipleNormalSpring::~MultipleNormalSpring()
{
  // invoke the destructor on any objects created by the object
  // that the object still holds a pointer to

  if (theMaterials != 0) {
    for (int i=0; i<nDivide; i++)
      if (theMaterials[i] != 0)
	delete theMaterials[i];
    delete [] theMaterials;
  }

  if (posLy != 0)
    delete [] posLy;

  if (posLz != 0)
    delete [] posLz;

  if (distFct != 0)
    delete [] distFct;
}


int MultipleNormalSpring::getNumExternalNodes() const
{
  return 2;
}


const ID& MultipleNormalSpring::getExternalNodes() 
{
  return connectedExternalNodes;
}


Node** MultipleNormalSpring::getNodePtrs() 
{
  return theNodes;
}


int MultipleNormalSpring::getNumDOF() 
{
  return 12;
}


void MultipleNormalSpring::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (!theDomain)  {
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    return;
  }
  
  // first set the node pointers
  theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
  theNodes[1] = theDomain->getNode(connectedExternalNodes(1));	
  
  // if can't find both - send a warning message
  if (!theNodes[0] || !theNodes[1])  {
    if (!theNodes[0])  {
      opserr << "WARNING MultipleNormalSpring::setDomain() - Nd1: " 
	     << connectedExternalNodes(0) << " does not exist in the model for ";
    } else  {
      opserr << "WARNING MultipleNormalSpring::setDomain() - Nd2: " 
	     << connectedExternalNodes(1) << " does not exist in the model for ";
    }
    opserr << "MultipleNormalSpring ele: " << this->getTag() << endln;
    
    return;
  }
  
  // now determine the number of dof and the dimension
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  
  // if differing dof at the ends - print a warning message
  if (dofNd1 != 6)  {
    opserr << "MultipleNormalSpring::setDomain() - node 1: "
	   << connectedExternalNodes(0) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  if (dofNd2 != 6)  {
    opserr << "MultipleNormalSpring::setDomain() - node 2: "
	   << connectedExternalNodes(1) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  
  // call the base class method
  this->DomainComponent::setDomain(theDomain);
  
  // set up the transformation matrix for orientation
  this->setUp();
}   	 


int MultipleNormalSpring::commitState()
{
  int errCode = 0;
  
  // commit material models
  for (int i=0; i<nDivide; i++)
    errCode += theMaterials[i]->commitState();
  
  return errCode;
}


int MultipleNormalSpring::revertToLastCommit()
{
  int errCode = 0;
  
  // revert material models
  for (int i=0; i<nDivide; i++)
    errCode += theMaterials[i]->revertToLastCommit();
  
  return errCode;
}


int MultipleNormalSpring::revertToStart()
{   
  int errCode=0;

  // trial variables
  basicDisp.Zero();
  basicForce.Zero();

  // reset basic stiffness matrix
  basicStiff = basicStiffInit;
  
  // revert material models
  for (int i=0; i<nDivide; i++)
    errCode += theMaterials[i]->revertToStart();
  
  return errCode;
}


int MultipleNormalSpring::update()
{
  // get global trial displacements and velocities
  const Vector &dsp1 = theNodes[0]->getTrialDisp();
  const Vector &dsp2 = theNodes[1]->getTrialDisp();
  const Vector &vel1 = theNodes[0]->getTrialVel();
  const Vector &vel2 = theNodes[1]->getTrialVel();
  
  static Vector globalDisp(12), globalDispDot(12);
  for (int i=0; i<6; i++)  {
    globalDisp(i)   = dsp1(i);  globalDispDot(i)   = vel1(i);
    globalDisp(i+6) = dsp2(i);  globalDispDot(i+6) = vel2(i);
  }

  static Vector localDispDot(12), basicDispDot(6);

  
  // transform response from the global to the local system
  localDisp    = Tgl*globalDisp;
  localDispDot = Tgl*globalDispDot;
  
  // transform response from the local to the basic system
  basicDisp    = Tlb*localDisp;
  basicDispDot = Tlb*localDispDot;


  // calculate shear forces and stiffnesses in basic y- and z-direction
  // get trial shear forces of hysteretic component
  basicForce.Zero();
  basicStiff.Zero();

  for (int i=0; i<(nDivide*nDivide); i++)  {

    //strain of each spring
    double tmpStrain     = (basicDisp(0) + basicDisp(4)*posLz[i] + basicDisp(5)*posLy[i])/hgt;
    double tmpStrainRate = (basicDispDot(0) + basicDispDot(4)*posLz[i] + basicDispDot(5)*posLy[i])/hgt;

    theMaterials[i]->setTrialStrain(tmpStrain,tmpStrainRate);

    //total force
    double tmpForce  = (theMaterials[i]->getStress())*incA*distFct[i];

    basicForce(0) += tmpForce ; //x,x
    basicForce(4) += tmpForce * posLz[i] ; //ry
    basicForce(5) += tmpForce * posLy[i] ; //rz

    //total stiffness
    double tmpStiff = (theMaterials[i]->getTangent())*incA*distFct[i]/hgt;

    basicStiff(0,0) += tmpStiff ; //x,x
    basicStiff(0,4) += tmpStiff * posLz[i] ; //x,ry
    basicStiff(0,5) += tmpStiff * posLy[i] ; //x,rz
    basicStiff(4,0) += tmpStiff * posLz[i] ; //ry,x
    basicStiff(4,4) += tmpStiff * posLz[i] * posLz[i] ; //ry,ry
    basicStiff(4,5) += tmpStiff * posLz[i] * posLy[i] ; //ry,rz
    basicStiff(5,0) += tmpStiff * posLy[i] ; //rz,x
    basicStiff(5,4) += tmpStiff * posLy[i] * posLz[i] ; //rz,ry
    basicStiff(5,5) += tmpStiff * posLy[i] * posLy[i] ; //rz,rz
  }
  
  return 0;
}


const Matrix& MultipleNormalSpring::getTangentStiff()
{
  // zero the matrix
  theMatrix.Zero();
  
  // transform from basic to local system
  static Matrix localStiff(12,12);
  localStiff.addMatrixTripleProduct(0.0, Tlb, basicStiff, 1.0);
  
  // transform from local to global system
  theMatrix.addMatrixTripleProduct(0.0, Tgl, localStiff, 1.0);
  
  return theMatrix;
}


const Matrix& MultipleNormalSpring::getInitialStiff()
{
  // zero the matrix
  theMatrix.Zero();
  
  // transform from basic to local system
  static Matrix localStiff(12,12);
  localStiff.addMatrixTripleProduct(0.0, Tlb, basicStiffInit, 1.0);
  
  // transform from local to global system
  theMatrix.addMatrixTripleProduct(0.0, Tgl, localStiff, 1.0);
  
  return theMatrix;
}


const Matrix& MultipleNormalSpring::getMass()
{
  // zero the matrix
  theMatrix.Zero();
  
  // check for quick return
  if (mass == 0.0)  {
    return theMatrix;
  }    
  
  double m = 0.5*mass;
  for (int i = 0; i < 3; i++)  {
    theMatrix(i,i)     = m;
    theMatrix(i+6,i+6) = m;
  }
  
  return theMatrix; 
}


void MultipleNormalSpring::zeroLoad()
{
  theLoad.Zero();
}


int MultipleNormalSpring::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
  opserr <<"MultipleNormalSpring::addLoad() - "
	 << "load type unknown for element: "
	 << this->getTag() << endln;
  
  return -1;
}


int MultipleNormalSpring::addInertiaLoadToUnbalance(const Vector &accel)
{
  // check for quick return
  if (mass == 0.0)  {
    return 0;
  }
  
  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
    opserr << "MultipleNormalSpring::addInertiaLoadToUnbalance() - "
	   << "matrix and vector sizes are incompatible\n";
    return -1;
  }
  
  // want to add ( - fact * M R * accel ) to unbalance
  // take advantage of lumped mass matrix
  double m = 0.5*mass;
  for (int i = 0; i < 3; i++)  {
    theLoad(i)   -= m * Raccel1(i);
    theLoad(i+6) -= m * Raccel2(i);
  }
  
  return 0;
}


const Vector& MultipleNormalSpring::getResistingForce()
{
  // zero the residual
  theVector.Zero();
  
  // determine resisting forces in local system
  static Vector localForce(12);
  localForce = Tlb^basicForce;
  
  // determine resisting forces in global system
  theVector = Tgl^localForce;
  
  // subtract external load
  theVector.addVector(1.0, theLoad, -1.0);
  
  return theVector;
}


const Vector& MultipleNormalSpring::getResistingForceIncInertia()
{	
  theVector = this->getResistingForce();
  
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    theVector += this->getRayleighDampingForces();
  
  // now include the mass portion
  if (mass != 0.0)  {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();    
    
    double m = 0.5*mass;
    for (int i = 0; i < 3; i++)  {
      theVector(i)   += m * accel1(i);
      theVector(i+6) += m * accel2(i);
    }
  }
  
  return theVector;
}


int MultipleNormalSpring::sendSelf(int commitTag, Channel &sChannel)
{
  return -1;
}


int MultipleNormalSpring::recvSelf(int commitTag, Channel &rChannel,
				  FEM_ObjectBroker &theBroker)
{
  return -1;
}


int MultipleNormalSpring::displaySelf(Renderer &theViewer,
				      int displayMode, float fact, const char **modes, int numModes)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}


void MultipleNormalSpring::Print(OPS_Stream &s, int flag)
{
  
  if (flag == 0)  {
    // print everything
    s << "Element: " << this->getTag(); 
    s << "  type: MultipleNormalSpring  iNode: " << connectedExternalNodes(0);
    s << "  jNode: " << connectedExternalNodes(1) << endln;
    s << "  Material : " << theMaterials[0]->getTag() << endln;
    s << "  mass: " << mass << endln;
    // determine resisting forces in global system
    s << "  resisting force: " << this->getResistingForce() << endln;
  } else if (flag == 1)  {
    // does nothing
  }
}


Response* MultipleNormalSpring::setResponse(const char **argv, int argc,
					   OPS_Stream &output)
{
  Response *theResponse = 0;
  
  output.tag("ElementOutput");
  output.attr("eleType","MultipleNormalSpring");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  
  // global forces
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
    {
      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");
      
      theResponse = new ElementResponse(this, 1, theVector);
    }
  // local forces
  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
    {
      output.tag("ResponseType","N_ 1");
      output.tag("ResponseType","Vy_1");
      output.tag("ResponseType","Vz_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Tz_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","T_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");
      
      theResponse = new ElementResponse(this, 2, theVector);
    }
  // basic forces
  else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0)
    {
      output.tag("ResponseType","qb1");
      output.tag("ResponseType","qb2");
      output.tag("ResponseType","qb3");
      output.tag("ResponseType","qb4");
      output.tag("ResponseType","qb5");
      output.tag("ResponseType","qb6");
      
      theResponse = new ElementResponse(this, 3, Vector(6));
    }
  // local displacements
  else if (strcmp(argv[0],"localDisplacement") == 0 ||
	   strcmp(argv[0],"localDisplacements") == 0)
    {
      output.tag("ResponseType","ux_1");
      output.tag("ResponseType","uy_1");
      output.tag("ResponseType","uz_1");
      output.tag("ResponseType","rx_1");
      output.tag("ResponseType","ry_1");
      output.tag("ResponseType","rz_1");
      output.tag("ResponseType","ux_2");
      output.tag("ResponseType","uy_2");
      output.tag("ResponseType","uz_2");
      output.tag("ResponseType","rx_2");
      output.tag("ResponseType","ry_2");
      output.tag("ResponseType","rz_2");
      
      theResponse = new ElementResponse(this, 4, theVector);
    }
  // basic displacements
  else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"deformations") == 0 || 
	   strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0 ||
	   strcmp(argv[0],"basicDisplacement") == 0 || strcmp(argv[0],"basicDisplacements") == 0)
    {
      output.tag("ResponseType","ub1");
      output.tag("ResponseType","ub2");
      output.tag("ResponseType","ub3");
      output.tag("ResponseType","ub4");
      output.tag("ResponseType","ub5");
      output.tag("ResponseType","ub6");
      
      theResponse = new ElementResponse(this, 5, Vector(6));
    }
  
  output.endTag(); // ElementOutput
  
  return theResponse;
}


int MultipleNormalSpring::getResponse(int responseID, Information &eleInfo)
{
  switch (responseID)  {
  case 1:  // global forces
    return eleInfo.setVector(this->getResistingForce());
    
  case 2:  // local forces
    theVector.Zero();
    // determine resisting forces in local system
    theVector = Tlb^basicForce;
    
    return eleInfo.setVector(theVector);
    
  case 3:  // basic forces
    return eleInfo.setVector(basicForce);
    
  case 4:  // local displacements
    return eleInfo.setVector(localDisp);
    
  case 5:  // basic displacements
    return eleInfo.setVector(basicDisp);
    
  default:
    return -1;
  }
}


// set up the transformation matrix for orientation
void MultipleNormalSpring::setUp()
{ 
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	
  Vector oriXp = end2Crd - end1Crd;
  double elmLen = oriXp.Norm();
  hgt = elmLen; //initial length of each spring

  if (elmLen > DBL_EPSILON)  {
    if (oriX.Size() == 0)  {
      oriX.resize(3);
      oriX = oriXp;
    } else  {
      opserr << "WARNING MultipleNormalSpring::setUp() - " 
	     << "element: " << this->getTag() << endln
	     << "ignoring nodes and using specified "
	     << "local x vector to determine orientation\n";
    }
  }
  // check that vectors for orientation are of correct size
  if (oriX.Size() != 3 || oriYp.Size() != 3)  {
    opserr << "MultipleNormalSpring::setUp() - "
	   << "element: " << this->getTag() << endln
	   << "incorrect dimension of orientation vectors\n";
    exit(-1);
  }
    
  // establish orientation of element for the transformation matrix
  // z = x cross yp
  Vector oriZ(3);
  oriZ(0) = oriX(1)*oriYp(2) - oriX(2)*oriYp(1);
  oriZ(1) = oriX(2)*oriYp(0) - oriX(0)*oriYp(2);
  oriZ(2) = oriX(0)*oriYp(1) - oriX(1)*oriYp(0);
  
  // y = z cross x
  Vector oriY(3);
  oriY(0) = oriZ(1)*oriX(2) - oriZ(2)*oriX(1);
  oriY(1) = oriZ(2)*oriX(0) - oriZ(0)*oriX(2);
  oriY(2) = oriZ(0)*oriX(1) - oriZ(1)*oriX(0);
    
  // compute length(norm) of vectors
  double xn = oriX.Norm();
  double yn = oriY.Norm();
  double zn = oriZ.Norm();
  
  // check valid x and y vectors, i.e. not parallel and of zero length
  if (xn == 0 || yn == 0 || zn == 0)  {
    opserr << "MultipleNormalSpring::setUp() - "
	   << "element: " << this->getTag() << endln
	   << "invalid orientation vectors\n";
    exit(-1);
  }
  
  // create transformation matrix from global to local system
  Tgl.Zero();
  Tgl(0,0) = Tgl(3,3) = Tgl(6,6) = Tgl(9,9)   = oriX(0)/xn;
  Tgl(0,1) = Tgl(3,4) = Tgl(6,7) = Tgl(9,10)  = oriX(1)/xn;
  Tgl(0,2) = Tgl(3,5) = Tgl(6,8) = Tgl(9,11)  = oriX(2)/xn;
  Tgl(1,0) = Tgl(4,3) = Tgl(7,6) = Tgl(10,9)  = oriY(0)/yn;
  Tgl(1,1) = Tgl(4,4) = Tgl(7,7) = Tgl(10,10) = oriY(1)/yn;
  Tgl(1,2) = Tgl(4,5) = Tgl(7,8) = Tgl(10,11) = oriY(2)/yn;
  Tgl(2,0) = Tgl(5,3) = Tgl(8,6) = Tgl(11,9)  = oriZ(0)/zn;
  Tgl(2,1) = Tgl(5,4) = Tgl(8,7) = Tgl(11,10) = oriZ(1)/zn;
  Tgl(2,2) = Tgl(5,5) = Tgl(8,8) = Tgl(11,11) = oriZ(2)/zn;
  
  // create transformation matrix from local to basic system (linear)
  Tlb.Zero();
  Tlb(0,0) = Tlb(1,1) = Tlb(2,2) = Tlb(3,3) = Tlb(4,4) = Tlb(5,5) = -1.0;
  Tlb(0,6) = Tlb(1,7) = Tlb(2,8) = Tlb(3,9) = Tlb(4,10) = Tlb(5,11) = 1.0;
  Tlb(1,5) = Tlb(1,11) = -0.5*elmLen;
  Tlb(2,4) = Tlb(2,10) = 0.5*elmLen;
}

