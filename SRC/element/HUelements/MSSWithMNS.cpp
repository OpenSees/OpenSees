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
// $Date: 2013-06-05 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/MSSWithMNS/MSSWithMNS.h,v $

// Written: Ken Ishii
// Created: Jan 2013
//
// MSSWithMNS model
//
// Description: This file contains the implementation of the MSSWithMNS class.
//              This file contains the function to parse the TCL input
// element MSSWithMNS eleTag? iNode? jNode? -shape shape? -size size? totalRubber?
//                    -nMSS nMSS? -matMSS matMSSTag? <-limDisp limDisp?>
//                    -nMNS nMNS? -matMNS matMNSTag? <-lambda lambda?>
//                    <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>

// This element works in 3-dim 6-dof model.
// The multiple normal springs are arranged parallel to local-x axis.
// The multiple shear springs are distributed in local-yz plane.


#include <TclModelBuilder.h>
#include <ID.h>
#include <Vector.h>

#include <MSSWithMNS.h>

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

#include <DBESI0.C> // Bessel function I0, Copyright(C) 1996 Takuya OOURA
#include <DBESI1.C> // Bessel function I1


extern void printCommand(int argc, TCL_Char **argv);


bool errDetected(bool ifNoError,char *msg){
  if (ifNoError){
    opserr << "" << endln;
    opserr << "========================================" << endln;
    opserr << "MSSWithMNS element : input error detected" << endln;
    opserr << "------------------------------" << endln;
  }
  opserr << "  " << msg << endln;
  return false;
};


int TclModelBuilder_addMSSWithMNS(ClientData clientData,
    Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
    TclModelBuilder *theTclBuilder)
{

  // ensure the destructor has not been called
  if (theTclBuilder == 0)  {
    opserr << "WARNING builder has been destroyed - MSSWithMNS\n";    
    return TCL_ERROR;
  }

  //3-dim, 6dof
  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING MSSWithMNS command only works when ndm is 3 and ndf is 6" << endln;
    return TCL_ERROR;
  }

  //arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;

  //arguments (necessary, input with -???)
  int shape;
  double size;
  double totalRubber;
  int nMSS;
  int matMSSTag;
  UniaxialMaterial *matMSS;
  int nMNS;
  int matMNSTag;
  UniaxialMaterial *matMNS;

  //arguments (optional, input with -???)
  double limDisp = 0.0;
  double lambda = -1.0;
  Vector oriX(0);
  Vector oriYp(3); oriYp(0) = 0.0; oriYp(1) = 1.0; oriYp(2) = 0.0;
  double mass = 0.0;

  // input comfirmation
  int recvShape  = 0;
  int recvSize   = 0;
  int recvNMSS  = 0;
  int recvMatMSS  = 0;
  int recvLimDisp = 0;
  int recvNMNS  = 0;
  int recvMatMNS  = 0;
  int recvLambda = 0;
  int recvOrient = 0;
  int recvMass   = 0;


  //
  Element *theElement = 0;


  //error flag
  bool ifNoError = true;



  if (argc < 5)  { //element MSSWithMNS eleTag? iNode? jNode?

    ifNoError = errDetected(ifNoError,"insufficient arguments");

  } else {

    //argv[2~4]
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK)  {
      ifNoError = errDetected(ifNoError,"invalid eleTag");
    }

    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK)  {
      ifNoError = errDetected(ifNoError,"invalid iNode");
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK)  {
      ifNoError = errDetected(ifNoError,"invalid jNode");
    }

    //argv[5~]
    for (int i=5; i<=(argc-1); i++) {
      
      double value;

      if (strcmp(argv[i],"-shape")==0 && (i+1)<=(argc-1)) { // -shape shape?
	
      	if (strcmp(argv[i+1],"round") == 0) {
	  shape = 1; //round
	} else if (strcmp(argv[i+1],"square") == 0) {
	  shape = 2; //square
	} else {
	  ifNoError = errDetected(ifNoError,"invalid shape (\"round\" or \"square\" are available)");
	}
	
	recvShape++ ;
	i += 1;


      } else if (strcmp(argv[i],"-size")==0 && (i+2)<=(argc-1)) { // -size size? totalRubber?
	
	if (Tcl_GetDouble(interp, argv[i+1], &size) != TCL_OK || size <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid size");
	}

	if (Tcl_GetDouble(interp, argv[i+2], &totalRubber) != TCL_OK || totalRubber <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid totalRubber");
	}

	recvSize++ ;
	i += 2;


      } else if (strcmp(argv[i],"-nMSS")==0 && (i+1)<=(argc-1)) { // -nMSS nMSS?
	
	if (Tcl_GetInt(interp, argv[i+1], &nMSS) != TCL_OK || nMSS <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid nMSS");
	}

	recvNMSS++ ;
	i += 1;


      } else if (strcmp(argv[i],"-matMSS")==0 && (i+1)<=(argc-1)) { // -matMSS matMSSTag?
	
	if (Tcl_GetInt(interp,argv[i+1], &matMSSTag) != TCL_OK) {
	  ifNoError = errDetected(ifNoError,"invalid matMSSTag");
	}

	matMSS = OPS_getUniaxialMaterial(matMSSTag);
	if (matMSS == 0)  {
	  ifNoError = errDetected(ifNoError,"material for MSS model not found");
	}

	recvMatMSS++ ;
	i += 1;


      } else if (strcmp(argv[i],"-limDisp")==0 && (i+1)<=(argc-1)) { // <-limDisp limDisp?>
	
	if (Tcl_GetDouble(interp, argv[i+1], &limDisp) != TCL_OK || limDisp < 0)  {
	  ifNoError = errDetected(ifNoError,"invalid limDisp");
	}

	recvLimDisp++ ;
	i += 1;


      } else if (strcmp(argv[i],"-nMNS")==0 && (i+1)<=(argc-1)) { // -nMNS nMNS?
	
	if (Tcl_GetInt(interp, argv[i+1], &nMNS) != TCL_OK || nMNS <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid nMNS");
	}

	recvNMNS++ ;
	i += 1;


      } else if (strcmp(argv[i],"-matMNS")==0 && (i+1)<=(argc-1)) { // -matMNS matMNSTag?
	
	if (Tcl_GetInt(interp,argv[i+1], &matMNSTag) != TCL_OK) {
	  ifNoError = errDetected(ifNoError,"invalid matMNSTag");
	}

	matMNS = OPS_getUniaxialMaterial(matMNSTag);
	if (matMNS == 0)  {
	  ifNoError = errDetected(ifNoError,"material for MNS model not found");
	}

	recvMatMNS++ ;
	i += 1;


      } else if (strcmp(argv[i],"-lambda")==0 && (i+1)<=(argc-1)) { // <-lambda lambda?>
	
	if (Tcl_GetDouble(interp, argv[i+1], &lambda) != TCL_OK || lambda < 0)  {
	  ifNoError = errDetected(ifNoError,"invalid lambda");
	}

	recvLambda++ ;
	i += 1;


      } else if (strcmp(argv[i],"-orient")==0 && (i+6)<=(argc-1) && Tcl_GetDouble(interp,argv[i+4], &value) == TCL_OK) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

	oriX.resize(3);

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    ifNoError = errDetected(ifNoError,"invalid orient");
	  } else {
	    oriX(j-1) = value;
	  }
	}

	i += 3;

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    ifNoError = errDetected(ifNoError,"invalid orient");
	  } else {
	    oriYp(j-1) = value;
	  }
	}

	recvOrient++ ;
	i += 3;

      } else if (strcmp(argv[i],"-orient")==0 && (i+3)<=(argc-1)) { // <-orient yp1? yp2? yp3?>

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    ifNoError = errDetected(ifNoError,"invalid orient");
	  } else {
	    oriYp(j-1) = value;
	  }
	}

	recvOrient++ ;
	i += 3;

      } else if (strcmp(argv[i],"-mass")==0 && (i+1)<=(argc-1)) { // <-mass mass?>
	
	if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK || mass <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid mass");
	}

	recvMass++ ;
	i += 1;
	
      } else { //invalid option
	
	ifNoError = errDetected(ifNoError,"invalid optional arguments");
	break;

      }

    }
    
  } //end input
  

  // input cofirmation
  // necessary arguments
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

  if (recvNMSS != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -NMSS inputs (got %d inputs, but want 1 input)",recvNMSS);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvMatMSS != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -matMSS inputs (got %d inputs, but want 1 input)",recvMatMSS);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvNMNS != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -NMNS inputs (got %d inputs, but want 1 input)",recvNMNS);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvMatMNS != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -matMNS inputs (got %d inputs, but want 1 input)",recvMatMNS);
    ifNoError = errDetected(ifNoError,buf);
  }


  //optional arguments
  if (recvLimDisp >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -limDisp inputs (got %d inputs, but want 1 input)",recvLimDisp);
    ifNoError = errDetected(ifNoError,buf);
  }

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

  
  //if error detected
  if (!ifNoError) {
    opserr << "------------------------------" << endln;
    //input:
    printCommand(argc, argv);
    //want:
    opserr << "Want: element MSSWithMNS eleTag? iNode? jNode? -shape shape? -size size? totalRubber?\n";
    opserr << "                         -nMSS nMSS? -matMSS matMSSTag? <-lim limDisp?> -nMNS nMNS? -matMNS matMNSTag? <-lambda lambda?>\n";
    opserr << "                         <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>\n";
    opserr << "========================================" << endln;
    opserr << "" << endln;
    return TCL_ERROR;
  }

  // now create the MSSWithMNS
  //  theElement = new MSSWithMNS(eleTag, iNode, jNode, nDivide, material, shape, size, lambda, oriYp, oriX, mass);
  theElement = new MSSWithMNS(eleTag, iNode, jNode, shape, size, totalRubber, nMSS, matMSS, limDisp, nMNS, matMNS, lambda, oriYp, oriX, mass);

  if (theElement == 0)  {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "MSSWithMNS element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  // then add the MSSWithMNS to the domain
  if (theTclDomain->addElement(theElement) == false)  {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "MSSWithMNS element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }       
  
  // if get here we have successfully created the MSSWithMNS and added it to the domain
  return TCL_OK;
}


// initialize the class wide variables
Matrix MSSWithMNS::theMatrix(12,12);
Vector MSSWithMNS::theVector(12);
Vector MSSWithMNS::theLoad(12);

MSSWithMNS::MSSWithMNS(int Tag, int Nd1, int Nd2,
		       int Shape, double Size, double TotalRubber,
		       int NMSS, UniaxialMaterial *MatMSS, double LimDisp,
		       int NMNS, UniaxialMaterial *MatMNS, double Lambda,
		       const Vector OriYp, const Vector OriX,
		       double Mass)
  : Element(Tag, ELE_TAG_MSSWithMNS),
    connectedExternalNodes(2),
    shape(Shape),size(Size),totalRubber(TotalRubber),
    nMSS(NMSS),limDisp(LimDisp),
    nMNS(NMNS),lambda(Lambda),
    oriX(OriX), oriYp(OriYp), mass(Mass),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6)
{
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "MSSWithMNS::setUp() - element: "
	   << this->getTag() << " failed to create an ID of size 2\n";
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;
  
  // check material input
  if (MatMNS == 0)  {
    opserr << "MSSWithMNS::MSSWithMNS() - "
	   << "null uniaxial material pointer passed.\n";
    exit(-1);
  }

  theINodeMNSMaterials = new UniaxialMaterial* [nMNS*nMNS]; // material‚ðnMNS^2ŒÂ iNode‘¤

  // get copies of the uniaxial materials
  for (int i=0; i<nMNS*nMNS; i++)  {
    theINodeMNSMaterials[i] = MatMNS->getCopy();

    if (theINodeMNSMaterials[i] == 0) {
      opserr << "MSSWithMNS::MSSWithMNS() - "
 	     << "failed to copy uniaxial material.\n";
      exit(-1);
    }
  }
  

  //position of centroid, distribution factor
  posLy = new double [nMNS*nMNS]; //local-y position
  posLz = new double [nMNS*nMNS]; //local-z position
  distFct = new double [nMNS*nMNS]; //distribution factor
  if (shape == 1){//round shape

    incA = (M_PI*size*size)/(4.0*nMNS*nMNS); //area of each normal spring

    int p = nMNS%2;
    int k = -1; //index
  
    for(int i=1; i<=((nMNS+p)/2); i++) {
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
      double r1 = ((2.0*i-2-p)/(2.0*nMNS))*size;
      double r2 = ((2.0*i-p)/(2.0*nMNS))*size;
      double r = (2.0/3.0)*(sin(tht)/tht)*((r1*r1+r1*r2+r2*r2)/(r1+r2));
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
    opserr << "MSSWithMNS::MSSWithMNS() - "
	   << "square shape \n";
    exit(-1);
  }
  
  // initial basic stiffness matrix
  basicStiffInit.Zero();
  for (int i=0; i<(nMNS*nMNS); i++)  {

    //total stiffness
    double tmpStiff = (theINodeMNSMaterials[i]->getInitialTangent())*incA*distFct[i]/(totalRubber/2.0);

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


MSSWithMNS::MSSWithMNS()
  : Element(0, ELE_TAG_MSSWithMNS),
    connectedExternalNodes(2),
    shape(0),size(0.0),totalRubber(0.0),
    nMSS(0),limDisp(0.0),
    nMNS(0),lambda(0.0),
    oriX(0), oriYp(0), mass(0.0),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6)
{	


  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "MSSWithMNS::MSSWithMNS() - "
	   <<  "failed to create an ID of size 2\n";
    exit(-1);
  }
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  // set material array pointers to NULL
  theINodeMNSMaterials = 0;

}


MSSWithMNS::~MSSWithMNS()
{
  // invoke the destructor on any objects created by the object
  // that the object still holds a pointer to

  if (theINodeMNSMaterials != 0) {
    for (int i=0; i<nMNS; i++)
      if (theINodeMNSMaterials[i] != 0)
	delete theINodeMNSMaterials[i];
    delete [] theINodeMNSMaterials;
  }

  if (posLy != 0)
    delete [] posLy;

  if (posLz != 0)
    delete [] posLz;

  if (distFct != 0)
    delete [] distFct;
}


int MSSWithMNS::getNumExternalNodes() const
{
  return 2;
}


const ID& MSSWithMNS::getExternalNodes() 
{
  return connectedExternalNodes;
}


Node** MSSWithMNS::getNodePtrs() 
{
  return theNodes;
}


int MSSWithMNS::getNumDOF() 
{
  return 12;
}


void MSSWithMNS::setDomain(Domain *theDomain)
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
      opserr << "WARNING MSSWithMNS::setDomain() - Nd1: " 
	     << connectedExternalNodes(0) << " does not exist in the model for ";
    } else  {
      opserr << "WARNING MSSWithMNS::setDomain() - Nd2: " 
	     << connectedExternalNodes(1) << " does not exist in the model for ";
    }
    opserr << "MSSWithMNS ele: " << this->getTag() << endln;
    
    return;
  }
  
  // now determine the number of dof and the dimension    
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  
  // if differing dof at the ends - print a warning message
  if (dofNd1 != 6)  {
    opserr << "MSSWithMNS::setDomain() - node 1: "
	   << connectedExternalNodes(0) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  if (dofNd2 != 6)  {
    opserr << "MSSWithMNS::setDomain() - node 2: "
	   << connectedExternalNodes(1) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  
  // call the base class method
  this->DomainComponent::setDomain(theDomain);
  
  // set up the transformation matrix for orientation
  this->setUp();
}   	 


int MSSWithMNS::commitState()
{
  int errCode = 0;
  
  // commit material models
  for (int i=0; i<nMNS; i++)
    errCode += theINodeMNSMaterials[i]->commitState();
  
  return errCode;
}


int MSSWithMNS::revertToLastCommit()
{
  int errCode = 0;
  
  // revert material models
  for (int i=0; i<nMNS; i++)
    errCode += theINodeMNSMaterials[i]->revertToLastCommit();
  
  return errCode;
}


int MSSWithMNS::revertToStart()
{   
  int errCode=0;

  // trial variables
  basicDisp.Zero();
  basicForce.Zero();

  // reset basic stiffness matrix
  basicStiff = basicStiffInit;
  
  // revert material models
  for (int i=0; i<nMNS; i++)
    errCode += theINodeMNSMaterials[i]->revertToStart();
  
  return errCode;
}


int MSSWithMNS::update()
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

  for (int i=0; i<(nMNS*nMNS); i++)  {

    //strain of each spring
    double tmpStrain     = (basicDisp(0) + basicDisp(4)*posLz[i] + basicDisp(5)*posLy[i])/(totalRubber/2.0);
    double tmpStrainRate = (basicDispDot(0) + basicDispDot(4)*posLz[i] + basicDispDot(5)*posLy[i])/(totalRubber/2.0);

    theINodeMNSMaterials[i]->setTrialStrain(tmpStrain,tmpStrainRate);

    //total force
    double tmpForce  = (theINodeMNSMaterials[i]->getStress())*incA*distFct[i];

    basicForce(0) += tmpForce ; //x,x
    basicForce(4) += tmpForce * posLz[i] ; //ry
    basicForce(5) += tmpForce * posLy[i] ; //rz

    //total stiffness
    double tmpStiff = (theINodeMNSMaterials[i]->getTangent())*incA*distFct[i]/(totalRubber/2.0);

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


const Matrix& MSSWithMNS::getTangentStiff()
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


const Matrix& MSSWithMNS::getInitialStiff()
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


const Matrix& MSSWithMNS::getMass()
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
    theMatrix(i+3,i+3) = m;
  }
  
  return theMatrix; 
}


void MSSWithMNS::zeroLoad()
{
  theLoad.Zero();
}


int MSSWithMNS::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
  opserr <<"MSSWithMNS::addLoad() - "
	 << "load type unknown for element: "
	 << this->getTag() << endln;
  
  return -1;
}


int MSSWithMNS::addInertiaLoadToUnbalance(const Vector &accel)
{
  // check for quick return
  if (mass == 0.0)  {
    return 0;
  }
  
  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
    opserr << "MSSWithMNS::addInertiaLoadToUnbalance() - "
	   << "matrix and vector sizes are incompatible\n";
    return -1;
  }
  
  // want to add ( - fact * M R * accel ) to unbalance
  // take advantage of lumped mass matrix
  double m = 0.5*mass;
  for (int i = 0; i < 3; i++)  {
    theLoad(i)   -= m * Raccel1(i);
    theLoad(i+3) -= m * Raccel2(i);
  }
  
  return 0;
}


const Vector& MSSWithMNS::getResistingForce()
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


const Vector& MSSWithMNS::getResistingForceIncInertia()
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
      theVector(i+3) += m * accel2(i);
    }
  }
  
  return theVector;
}


int MSSWithMNS::sendSelf(int commitTag, Channel &sChannel)
{
  return -1;
}


int MSSWithMNS::recvSelf(int commitTag, Channel &rChannel,
				  FEM_ObjectBroker &theBroker)
{
  return -1;
}


int MSSWithMNS::displaySelf(Renderer &theViewer,
				     int displayMode, float fact)
{
  // first determine the end points of the element based on
  // the display factor (a measure of the distorted image)
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	
  
  const Vector &end1Disp = theNodes[0]->getDisp();
  const Vector &end2Disp = theNodes[1]->getDisp();
  
  static Vector v1(3);
  static Vector v2(3);
  
  for (int i = 0; i < 3; i++)  {
    v1(i) = end1Crd(i) + end1Disp(i)*fact;
    v2(i) = end2Crd(i) + end2Disp(i)*fact;    
  }
  
  return theViewer.drawLine (v1, v2, 1.0, 1.0);
}


void MSSWithMNS::Print(OPS_Stream &s, int flag)
{
  
  if (flag == 0)  {
    // print everything
    s << "Element: " << this->getTag(); 
    s << "  type: MSSWithMNS  iNode: " << connectedExternalNodes(0);
    s << "  jNode: " << connectedExternalNodes(1) << endln;
    s << "  Material : " << theINodeMNSMaterials[0]->getTag() << endln;
    s << "  mass: " << mass << endln;
    // determine resisting forces in global system
    s << "  resisting force: " << this->getResistingForce() << endln;
  } else if (flag == 1)  {
    // does nothing
  }
}


Response* MSSWithMNS::setResponse(const char **argv, int argc,
					   OPS_Stream &output)
{
  Response *theResponse = 0;
  
  output.tag("ElementOutput");
  output.attr("eleType","MSSWithMNS");
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


int MSSWithMNS::getResponse(int responseID, Information &eleInfo)
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
void MSSWithMNS::setUp()
{ 
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	
  Vector oriXp = end2Crd - end1Crd;
  //double totalHeight = oriXp.Norm();
  totalHeight = oriXp.Norm(); // total height of bearing
  //hgt = totalRubber/2.0 ; //imaginary length of axial spring

  if (totalHeight > DBL_EPSILON)  {
    if (oriX.Size() == 0)  {
      oriX.resize(3);
      oriX = oriXp;
    } else  {
      opserr << "WARNING MSSWithMNS::setUp() - " 
	     << "element: " << this->getTag() << endln
	     << "ignoring nodes and using specified "
	     << "local x vector to determine orientation\n";
    }
  }
  // check that vectors for orientation are of correct size
  if (oriX.Size() != 3 || oriYp.Size() != 3)  {
    opserr << "MSSWithMNS::setUp() - "
	   << "element: " << this->getTag() << endln
	   << "incorrect dimension of orientation vectors\n";
    exit(-1);
  }
    
  // establish orientation of element for the tranformation matrix
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
    opserr << "MSSWithMNS::setUp() - "
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
  Tlb(1,5) = Tlb(1,11) = -0.5*totalHeight;
  Tlb(2,4) = Tlb(2,10) = 0.5*totalHeight;
}

