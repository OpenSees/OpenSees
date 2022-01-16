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
// $Date: 2013-05-31 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/multipleShearSpring/MultipleShearSpring.h,v $

// Written: Ken Ishii
// Created: June 2012
//
// Multiple Shear Spring (MSS) model
//
// Description: This file contains the implementation of the MultipleShearSpring class.
//              This file contains the function to parse the TCL input
// element multipleShearSpring eleTag? iNode? jNode? nSpring? -mat matTag? <-lim limDisp?>  <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>

// This element works in 3-dim 6-dof model.
// The multiple shear springs are distributed in local-yz plane.
// Stiffness of local-x,rx,ry,rz remain zero.

#include <MultipleShearSpring.h>

#include <Domain.h>
#include <ID.h>
#include <Vector.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <UniaxialMaterial.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>


void * OPS_ADD_RUNTIME_VPV(OPS_MultipleShearSpring)
{
    // 3-dim, 6-dof
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    if (ndm != 3 || ndf != 6) {
	opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
	opserr << "WARNING multipleShearSpring command only works when ndm is 3 and ndf is 6" << endln;
	return 0;
    }

    //arguments (necessary)
    int eleTag;
    int iNode;
    int jNode;
    int nSpring;
    int matTag;

    //material
    UniaxialMaterial *material = 0;
    UniaxialMaterial **theMaterials = 0;
    int recvMat = 0;

    //arguments (optional)
    double limDisp = 0.0;
    Vector oriX(0);
    Vector oriYp(3); oriYp(0) = 0.0; oriYp(1) = 1.0; oriYp(2) = 0.0;
    double mass = 0.0;

    //
    Element *theElement = 0;


    //error flag
    bool ifNoError = true;


    if (OPS_GetNumRemainingInputArgs() < 6)  { //element multipleShearSpring eleTag? iNode? jNode? nSpring? -mat matTag?

	opserr << "WARNING insufficient arguments\n";
	ifNoError = false;

    } else {

	//argv[2~5]
	int idata[4];
	int numdata = 4;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
	    opserr << "WARNING invalid multipleShearSpring int inputs\n";
	    ifNoError = false;
	}
	eleTag = idata[0];
	iNode = idata[1];
	jNode = idata[2];
	nSpring = idata[3];
	if (nSpring <= 0) {
	    opserr << "WARNING invalid nSpring\n";
	    ifNoError = false;
	}

	//argv[6~]
	while (OPS_GetNumRemainingInputArgs() > 0) {
	    const char* flag = OPS_GetString();
      
	    double value;
      
	    if (strcmp(flag,"-mat")==0 && OPS_GetNumRemainingInputArgs()>0) { // -mat matTag?
		numdata = 1;
		if (OPS_GetIntInput(&numdata, &matTag) < 0) {
		    opserr << "WARNING invalid matTag\n";
		    ifNoError = false;
		}

		material = OPS_getUniaxialMaterial(matTag);
		if (material == 0)  {
		    opserr << "WARNING material model not found\n";
		    opserr << "uniaxialMaterial: " << matTag << endln;
		    opserr << "multipleShearSpring element: " << eleTag << endln;
		    return 0;
		}

		//opserr << "org material " << material->getClassType() << "\n";
		recvMat++ ;

	    } else if (strcmp(flag,"-nMat")==0 && OPS_GetNumRemainingInputArgs()>=nSpring) { // -mat matTag?

		theMaterials = new UniaxialMaterial *[nSpring];
		for (int j=0; j<nSpring; j++) {
		    numdata = 1;
		    if (OPS_GetIntInput(&numdata, &matTag) < 0) {
			opserr << "WARNING invalid matTag\n";
			ifNoError = false;
		    }
	  
		    theMaterials[j] = OPS_getUniaxialMaterial(matTag);
		    if (theMaterials[j] == 0)  {
			opserr << "WARNING material model not found\n";
			opserr << "uniaxialMaterial: " << matTag << endln;
			opserr << "multipleShearSpring element: " << eleTag << endln;
			return 0;
		    }
		}
		//opserr << "org material " << material->getClassType() << "\n";
		recvMat++ ;

	    } else if (strcmp(flag,"-orient")==0 && OPS_GetNumRemainingInputArgs()>=6) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

		oriX.resize(3);

		for (int j=1; j<=3; j++) {
		    numdata = 1;
		    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
			opserr << "WARNING invalid -orient value\n";
			ifNoError = false;
		    } else {
			oriX(j-1) = value;
		    }
		}

		for (int j=1; j<=3; j++) {
		    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
			opserr << "WARNING invalid -orient value\n";
			ifNoError = false;
		    } else {
			oriYp(j-1) = value;
		    }
		}

	    } else if (strcmp(flag,"-orient")==0 && OPS_GetNumRemainingInputArgs()>=3) { // <-orient yp1? yp2? yp3?> �̓ǂݍ���	  

		for (int j=1; j<=3; j++) {
		    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
			opserr << "WARNING invalid -orient value\n";
			ifNoError = false;
		    } else {
			oriYp(j-1) = value;
		    }
		}

	    } else if (strcmp(flag,"-mass")==0 && OPS_GetNumRemainingInputArgs()>0) { // <-mass m?> �̓ǂݍ���
		if (OPS_GetDoubleInput(&numdata, &mass) < 0 || mass <= 0) {
		    opserr << "WARNING invalid mass\n";
		    ifNoError = false;
		}

	    } else if (strcmp(flag,"-lim")==0 && OPS_GetNumRemainingInputArgs()>0) { // <-lim limDisp?> �̓ǂݍ���
		if (OPS_GetDoubleInput(&numdata, &limDisp) < 0 || limDisp < 0) {
		    opserr << "WARNING invalid limDisp\n";
		    ifNoError = false;
		}

	    } else { //invalid option
	
		opserr << "WARNING invalid optional arguments \n";
		ifNoError = false;
		break;

	    }

	}
    
    } //end input
  


    //confirm material
    if (recvMat != 1)  {
	opserr << "WARNING wrong number of -mat inputs\n";
	opserr << "got " << recvMat << " inputs, but want 1 input\n";
	ifNoError = false;
    }

  
    //if error detected
    if (!ifNoError) {
	//input:
	//want:
	opserr << "Want: element multipleShearSpring eleTag? iNode? jNode? nSpring? -mat matTag? <-lim dsp> <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>\n";
	return 0;
    }
  

    // now create the multipleShearSpring
    if (theMaterials == 0) {
	theElement = new MultipleShearSpring(eleTag, iNode, jNode, nSpring, material, limDisp, oriYp, oriX, mass);
    } else {
	theElement = new MultipleShearSpring(eleTag, iNode, jNode, theMaterials, nSpring, limDisp, oriYp, oriX, mass);
	delete [] theMaterials;
    }

    return theElement;
}



// initialize the class wide variables
Matrix MultipleShearSpring::theMatrix(12,12);
Vector MultipleShearSpring::theVector(12);
Vector MultipleShearSpring::theLoad(12);


MultipleShearSpring::MultipleShearSpring(int Tag, int Nd1, int Nd2,
					 int NSpring,
					 UniaxialMaterial *Material,
					 double LimDisp,
					 const Vector OriYp, const Vector OriX, double Mass)
  : Element(Tag, ELE_TAG_MultipleShearSpring),
    connectedExternalNodes(2),
    nSpring(NSpring), limDisp(LimDisp), oriX(OriX), oriYp(OriYp), mass(Mass),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6)
{
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "MultipleShearSpring::setUp() - element: "
	   << this->getTag() << " failed to create an ID of size 2\n";
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;
  
  // check material input
  if (Material == 0)  {
    opserr << "MultipleShearSpring::MultipleShearSpring() - "
	   << "null uniaxial material pointer passed.\n";
    exit(-1);
  }

  theMaterials = new UniaxialMaterial* [nSpring];

  // get copies of the uniaxial materials
  for (int i=0; i<nSpring; i++)  {
    theMaterials[i] = Material->getCopy();

    if (theMaterials[i] == 0) {
      opserr << "MultipleShearSpring::MultipleShearSpring() - "
 	     << "failed to copy uniaxial material.\n";
      exit(-1);
    }
  }

  
  //arrangement of each spring
  cosTht = new double [nSpring];
  sinTht = new double [nSpring];
  
  for (int i=0; i<nSpring; i++)  {
    cosTht[i] = cos(M_PI*i/nSpring);
    sinTht[i] = sin(M_PI*i/nSpring);
  }


  //calculate Feq and Seq
  //imaginary material to calculate Feq and Seq
  dmyMssMaterial = Material->getCopy();
  if (dmyMssMaterial == 0) {
    opserr << "MultipleShearSpring::MultipleShearSpring() - "
	   << "failed to copy uniaxial material.\n";
    exit(-1);
  }
  dmyMssMaterial->revertToStart();

  //initial Feq and Seq
  if (limDisp > 0) {
    double uRef, fRef, sRef;//u:deformation, f:force, s:stiffness
    double uCmp, fSum, sSum;

    //imaginary material (1-directional deformation)
    uRef = limDisp;
    dmyMssMaterial->setTrialStrain(uRef,0);
    fRef = dmyMssMaterial->getStress();
    sRef = dmyMssMaterial->getTangent();

    //MSS
    fSum = 0.0;
    sSum = 0.0;
    for (int i=0; i<nSpring; i++)  {
      uCmp = uRef * cosTht[i];
      dmyMssMaterial->setTrialStrain(uCmp,0);
      fSum += dmyMssMaterial->getStress()  * cosTht[i];
      sSum += dmyMssMaterial->getTangent() * cosTht[i] * cosTht[i];
    }
    mssFeq = fRef/fSum;
    mssSeq = sRef/sSum;

  } else {

    mssFeq = 1.0;
    mssSeq = 1.0;

  }


  // initial basic stiffness matrix
  basicStiffInit.Zero();
  for (int i=0; i<nSpring; i++)  {
    double tmpTangent = theMaterials[i]->getInitialTangent();

    basicStiffInit(1,1) += tmpTangent * cosTht[i] * cosTht[i];
    basicStiffInit(1,2) += tmpTangent * cosTht[i] * sinTht[i];
    basicStiffInit(2,1) += tmpTangent * sinTht[i] * cosTht[i];
    basicStiffInit(2,2) += tmpTangent * sinTht[i] * sinTht[i];
  }


  // Seq: equivalent coefficient for stiffness
  basicStiffInit *= mssSeq;


  // initialize variables
  this->revertToStart();
}




MultipleShearSpring::MultipleShearSpring(int Tag, int Nd1, int Nd2,
					 UniaxialMaterial **theMats,
					 int NSpring,
					 double LimDisp,
					 const Vector OriYp, const Vector OriX, double Mass)
  : Element(Tag, ELE_TAG_MultipleShearSpring),
    connectedExternalNodes(2),
    nSpring(NSpring), limDisp(LimDisp), oriX(OriX), oriYp(OriYp), mass(Mass),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6)
{
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "MultipleShearSpring::setUp() - element: "
	   << this->getTag() << " failed to create an ID of size 2\n";
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;
  
  // check material input
  if (theMats == 0)  {
    opserr << "MultipleShearSpring::MultipleShearSpring() - "
	   << "null uniaxial material pointer passed.\n";
    exit(-1);
  }

  theMaterials = new UniaxialMaterial* [nSpring];

  // get copies of the uniaxial materials
  for (int i=0; i<nSpring; i++)  {
    if (theMats[i] != 0) {
      theMaterials[i] = theMats[i]->getCopy();
    } else {
      theMaterials[i] = 0;
    }
	
    if (theMaterials[i] == 0) {
      opserr << "MultipleShearSpring::MultipleShearSpring() - "
 	     << "failed to copy uniaxial material.\n";
      exit(-1);
    }
  }
  
  //arrangement of each spring
  cosTht = new double [nSpring];
  sinTht = new double [nSpring];
  
  for (int i=0; i<nSpring; i++)  {
    cosTht[i] = cos(M_PI*i/nSpring);
    sinTht[i] = sin(M_PI*i/nSpring);
  }


  //calculate Feq and Seq
  //imaginary material to calculate Feq and Seq
  dmyMssMaterial = theMaterials[0]->getCopy();
  if (dmyMssMaterial == 0) {
    opserr << "MultipleShearSpring::MultipleShearSpring() - "
	   << "failed to copy uniaxial material.\n";
    exit(-1);
  }
  dmyMssMaterial->revertToStart();

  //initial Feq and Seq
  if (limDisp > 0) {
    double uRef, fRef, sRef;//u:deformation, f:force, s:stiffness
    double uCmp, fSum, sSum;

    //imaginary material (1-directional deformation)
    uRef = limDisp;
    dmyMssMaterial->setTrialStrain(uRef,0);
    fRef = dmyMssMaterial->getStress();
    sRef = dmyMssMaterial->getTangent();

    //MSS
    fSum = 0.0;
    sSum = 0.0;
    for (int i=0; i<nSpring; i++)  {
      uCmp = uRef * cosTht[i];
      dmyMssMaterial->setTrialStrain(uCmp,0);
      fSum += dmyMssMaterial->getStress()  * cosTht[i];
      sSum += dmyMssMaterial->getTangent() * cosTht[i] * cosTht[i];
    }
    mssFeq = fRef/fSum;
    mssSeq = sRef/sSum;

  } else {

    mssFeq = 1.0;
    mssSeq = 1.0;

  }


  // initial basic stiffness matrix
  basicStiffInit.Zero();
  for (int i=0; i<nSpring; i++)  {
    double tmpTangent = theMaterials[i]->getInitialTangent();

    basicStiffInit(1,1) += tmpTangent * cosTht[i] * cosTht[i];
    basicStiffInit(1,2) += tmpTangent * cosTht[i] * sinTht[i];
    basicStiffInit(2,1) += tmpTangent * sinTht[i] * cosTht[i];
    basicStiffInit(2,2) += tmpTangent * sinTht[i] * sinTht[i];
  }


  // Seq: equivalent coefficient for stiffness
  basicStiffInit *= mssSeq;


  // initialize variables
  this->revertToStart();
}


MultipleShearSpring::MultipleShearSpring()
  : Element(0, ELE_TAG_MultipleShearSpring),
    connectedExternalNodes(2),
    nSpring(0), limDisp(0), oriX(0), oriYp(0), mass(0.0),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6)
{	


  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "MultipleShearSpring::MultipleShearSpring() - "
	   <<  "failed to create an ID of size 2\n";
    exit(-1);
  }
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  // set material array pointers to NULL
  theMaterials = 0;
  dmyMssMaterial = 0;

}


MultipleShearSpring::~MultipleShearSpring()
{
  // invoke the destructor on any objects created by the object
  // that the object still holds a pointer to

  if (theMaterials != 0) {
    for (int i=0; i<nSpring; i++)
      if (theMaterials[i] != 0)
	delete theMaterials[i];
    delete [] theMaterials;
  }

  if (cosTht != 0)
    delete [] cosTht;

  if (sinTht != 0)
    delete [] sinTht;

  if (dmyMssMaterial != 0)
    delete dmyMssMaterial;


}


int MultipleShearSpring::getNumExternalNodes() const
{
  return 2;
}


const ID& MultipleShearSpring::getExternalNodes() 
{
  return connectedExternalNodes;
}


Node** MultipleShearSpring::getNodePtrs() 
{
  return theNodes;
}


int MultipleShearSpring::getNumDOF() 
{
  return 12;
}


void MultipleShearSpring::setDomain(Domain *theDomain)
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
      opserr << "WARNING MultipleShearSpring::setDomain() - Nd1: " 
	     << connectedExternalNodes(0) << " does not exist in the model for ";
    } else  {
      opserr << "WARNING MultipleShearSpring::setDomain() - Nd2: " 
	     << connectedExternalNodes(1) << " does not exist in the model for ";
    }
    opserr << "MultipleShearSpring ele: " << this->getTag() << endln;
    
    return;
  }
  
  // now determine the number of dof and the dimension    
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  
  // if differing dof at the ends - print a warning message
  if (dofNd1 != 6)  {
    opserr << "MultipleShearSpring::setDomain() - node 1: "
	   << connectedExternalNodes(0) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  if (dofNd2 != 6)  {
    opserr << "MultipleShearSpring::setDomain() - node 2: "
	   << connectedExternalNodes(1) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  
  // call the base class method
  this->DomainComponent::setDomain(theDomain);
  
  // set up the transformation matrix for orientation
  this->setUp();
}   	 


int MultipleShearSpring::commitState()
{
  int errCode = 0;
  
  // commit material models
  for (int i=0; i<nSpring; i++)
    errCode += theMaterials[i]->commitState();
  
  return errCode;
}


int MultipleShearSpring::revertToLastCommit()
{
  int errCode = 0;
  
  // revert material models
  for (int i=0; i<nSpring; i++)
    errCode += theMaterials[i]->revertToLastCommit();
  
  return errCode;
}


int MultipleShearSpring::revertToStart()
{   
  int errCode=0;

  // trial variables
  basicDisp.Zero();
  basicForce.Zero();

  // reset basic stiffness matrix
  basicStiff = basicStiffInit;
  
  // revert material models
  for (int i=0; i<nSpring; i++)
    errCode += theMaterials[i]->revertToStart();
  
  return errCode;
}


int MultipleShearSpring::update()
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

  for (int i=0; i<nSpring; i++)  {
    double tmpStrain     = basicDisp(1)*cosTht[i]    + basicDisp(2)*sinTht[i];
    double tmpStrainRate = basicDispDot(1)*cosTht[i] + basicDispDot(2)*sinTht[i];

    theMaterials[i]->setTrialStrain(tmpStrain,tmpStrainRate);

    double tmpStress  = theMaterials[i]->getStress();

    basicForce(1) += tmpStress * cosTht[i]; //basic-y
    basicForce(2) += tmpStress * sinTht[i]; //basic-z

    double tmpTangent = theMaterials[i]->getTangent();

    basicStiff(1,1) += tmpTangent * cosTht[i] * cosTht[i];
    basicStiff(1,2) += tmpTangent * cosTht[i] * sinTht[i];
    basicStiff(2,1) += tmpTangent * sinTht[i] * cosTht[i];
    basicStiff(2,2) += tmpTangent * sinTht[i] * sinTht[i];
  }

  // calculate Feq and Seq
  if (limDisp > 0) {
    double uRef, fRef, sRef;//u:deformation, f:force, s:stiffness
    double uCmp, fSum, sSum;
    double refDisp;

    //imaginary material (1-directional deformation)
    refDisp = sqrt(basicDisp(1)*basicDisp(1)+basicDisp(2)*basicDisp(2));
    uRef = (refDisp>limDisp) ? refDisp : limDisp;

    dmyMssMaterial->setTrialStrain(uRef,0);
    fRef = dmyMssMaterial->getStress();
    sRef = dmyMssMaterial->getTangent();

    //MSS
    fSum = 0.0;
    sSum = 0.0;
    for (int i=0; i<nSpring; i++)  {
      uCmp = uRef * cosTht[i];
      dmyMssMaterial->setTrialStrain(uCmp,0);
      fSum += dmyMssMaterial->getStress() * cosTht[i];
      sSum += dmyMssMaterial->getTangent() * cosTht[i] * cosTht[i];
    }
    mssFeq = fRef/fSum;
    mssSeq = sRef/sSum;

  }

  //  opserr << "forceFactor: " << mssFeq << " stiffFactor: " << mssSeq << endln;
  basicForce *= mssFeq;
  basicStiff *= mssSeq;

  
  return 0;
}


const Matrix& MultipleShearSpring::getTangentStiff()
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


const Matrix& MultipleShearSpring::getInitialStiff()
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


const Matrix& MultipleShearSpring::getMass()
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


void MultipleShearSpring::zeroLoad()
{
  theLoad.Zero();
}


int MultipleShearSpring::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
  opserr <<"MultipleShearSpring::addLoad() - "
	 << "load type unknown for element: "
	 << this->getTag() << endln;
  
  return -1;
}


int MultipleShearSpring::addInertiaLoadToUnbalance(const Vector &accel)
{
  // check for quick return
  if (mass == 0.0)  {
    return 0;
  }
  
  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
    opserr << "MultipleShearSpring::addInertiaLoadToUnbalance() - "
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


const Vector& MultipleShearSpring::getResistingForce()
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


const Vector& MultipleShearSpring::getResistingForceIncInertia()
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


int MultipleShearSpring::sendSelf(int commitTag, Channel &sChannel)
{
  return -1;
}


int MultipleShearSpring::recvSelf(int commitTag, Channel &rChannel,
				  FEM_ObjectBroker &theBroker)
{
  return -1;
}


int MultipleShearSpring::displaySelf(Renderer &theViewer,
				     int displayMode, float fact, const char **modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}


void MultipleShearSpring::Print(OPS_Stream &s, int flag)
{
  
  if (flag == 0)  {
    // print everything
    s << "Element: " << this->getTag(); 
    s << "  type: MultipleShearSpring  iNode: " << connectedExternalNodes(0);
    s << "  jNode: " << connectedExternalNodes(1) << endln;
    s << "  Material : " << theMaterials[0]->getTag() << endln;
    s << "  mass: " << mass << endln;
    // determine resisting forces in global system
    s << "  resisting force: " << this->getResistingForce() << endln;
  } else if (flag == 1)  {
    // does nothing
  }
}


Response* MultipleShearSpring::setResponse(const char **argv, int argc,
					   OPS_Stream &output)
{
  Response *theResponse = 0;
  
  output.tag("ElementOutput");
  output.attr("eleType","MultipleShearSpring");
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


int MultipleShearSpring::getResponse(int responseID, Information &eleInfo)
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
void MultipleShearSpring::setUp()
{ 
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	
  Vector oriXp = end2Crd - end1Crd;
  double elmLen = oriXp.Norm();
  
  if (elmLen > DBL_EPSILON)  {
    if (oriX.Size() == 0)  {
      oriX.resize(3);
      oriX = oriXp;
    } else  {
      opserr << "WARNING MultipleShearSpring::setUp() - " 
	     << "element: " << this->getTag() << endln
	     << "ignoring nodes and using specified "
	     << "local x vector to determine orientation\n";
    }
  }
  // check that vectors for orientation are of correct size
  if (oriX.Size() != 3 || oriYp.Size() != 3)  {
    opserr << "MultipleShearSpring::setUp() - "
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
    opserr << "MultipleShearSpring::setUp() - "
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

