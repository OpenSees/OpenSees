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
// $Date: 2012-08-24 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/yamamotoBiaxialHDR/YamamotoBiaxialHDR.h,v $

// Written: Masaru Kikuchi
// Created: August 2014
//
// Description: This file contains the implementation of the YamamotoBiaxialHDR class.
//
//

// element YamamotoBiaxialHDR eleTag? iNode? jNode? Tp? DDo? DDi? Hr?  <-coRS cr? cs?> <-orient <x1? x2? x3?> y1? y2? y3?> <-mass m?>

// This element works in 3-dim 6-dof model.
// The multiple shear springs are distributed in local-yz plane.
// Stiffness of local-x,rx,ry,rz remain zero.

#include <YamamotoBiaxialHDR.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>

#include <ID.h>
#include <Vector.h>

#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_YamamotoBiaxialHDR)
{
    // 3-dim, 6-dof
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    if (ndm != 3 || ndf != 6) {
	opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
	opserr << "WARNING YamamotoBiaxialHDR command only works when ndm is 3 and ndf is 6" << endln;
	return 0;
    }


    //arguments (necessary)
    int eleTag;
    int iNode;
    int jNode;

    int Tp = 1;
    double DDo;
    double DDi;
    double Hr;

    //arguments (optional)
    double Cr=1.0;
    double Cs=1.0;
    Vector oriX(0);
    Vector oriYp(3); oriYp(0) = 0.0; oriYp(1) = 1.0; oriYp(2) = 0.0;
    double mass = 0.0;


    //
    Element *theElement = 0;


    //error flag
    bool ifNoError = true;

    // 
    int numdata = 1;

    if (OPS_GetNumRemainingInputArgs() < 7)  {
        // element YamamotoBiaxialHDR eleTag? iNode? jNode? Tp? DDo? DDi? Hr?
	// argc =            1           2             3      4      5     6   7    8    9
	// argv =       argv[0]      argv[1]      argv[2]  argv[3] ................. argv[8]
	opserr << "WARNING insufficient arguments\n";
	ifNoError = false;
	
    } else {


	//argv[2~8]
	if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
	    opserr << "WARNING invalid YamamotoBiaxialHDR eleTag\n";
	    ifNoError = false;
	}

	// iNode
	if (OPS_GetIntInput(&numdata, &iNode) < 0) {
	    opserr << "WARNING invalid iNode\n";
	    ifNoError = false;
	}

	// jNode
	if (OPS_GetIntInput(&numdata, &jNode) < 0) {
	    opserr << "WARNING invalid jNode\n";
	    ifNoError = false;
	}

	// Tp
	const char* tparg = OPS_GetString();
	if (strcmp(tparg,"1") == 0) {
	    Tp = 1; // Bridgestone X0.6R (EESD version)
	} else {
	    opserr << "WARNING invalid YamamotoBiaxialHDR Tp" << endln;
	    ifNoError = false;
	}

	// DDo
	if (OPS_GetDoubleInput(&numdata, &DDo) < 0 || DDo <= 0.0) {
	    opserr << "WARNING invalid YamamotoBiaxialHDR DDo" << endln;
	    ifNoError = false;
	}

	// DDi
	if (OPS_GetDoubleInput(&numdata, &DDi) < 0 || DDi < 0.0) {
	    opserr << "WARNING invalid YamamotoBiaxialHDR DDi" << endln;
	    ifNoError = false;
	}

	// Hr
	if (OPS_GetDoubleInput(&numdata, &Hr) < 0 || Hr <= 0.0) {
	    opserr << "WARNING invalid YamamotoBiaxialHDR Hr" << endln;
	    ifNoError = false;
	}

	// check print--------------------------------------------/
	//  opserr << "   \n";
	//  opserr << "TclModelBuilder_addYamamotoBiaxialHDR()\n";
	//  opserr << "  tp  = " << Tp << endln;
	//  opserr << "  ddo = " << DDo << endln;
	//  opserr << "  ddi = " << DDi << endln;
	//  opserr << "  hr  = " << Hr << endln;
	//------------------------------------------------------

	// argv[9~]
	while (OPS_GetNumRemainingInputArgs() > 0) {
	    double value;
	    const char* flag = OPS_GetString();
	    if (strcmp(flag,"-orient")==0 && OPS_GetNumRemainingInputArgs() >= 6) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

		oriX.resize(3);

		// x1, x2, x3
		for (int j=1; j<=3; j++) {
		    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
			opserr << "WARNING invalid -orient value\n";
			ifNoError = false;
		    } else {
			oriX(j-1) = value;
		    }
		}
	
		// yp1, yp2, yp3
		for (int j=1; j<=3; j++) {
		    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
			opserr << "WARNING invalid -orient value\n";
			ifNoError = false;
		    } else {
			oriYp(j-1) = value;
		    }
		}
	
	    } else if (strcmp(flag,"-orient")==0 && OPS_GetNumRemainingInputArgs() >= 3) { // <-orient yp1? yp2? yp3?>
	
		for (int j=1; j<=3; j++) {
		    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
			opserr << "WARNING invalid -orient value\n";
			ifNoError = false;
		    } else {
			oriYp(j-1) = value;
		    }
		}
	
	    } else if (strcmp(flag,"-mass")==0 && OPS_GetNumRemainingInputArgs()>0) { // <-mass m?>

		if (OPS_GetDoubleInput(&numdata, &mass) < 0 || mass <= 0) {
		    opserr << "WARNING invalid mass\n";
		    ifNoError = false;
		}

	    } else if (strcmp(flag,"-coRS")==0 && OPS_GetNumRemainingInputArgs()>=2) { // <-coRS cr? cs?>
		if (OPS_GetDoubleInput(&numdata, &Cr) < 0 || Cr <= 0) {
		    opserr << "WARNING invalid cr\n";
		    ifNoError = false;
		}
		if (OPS_GetDoubleInput(&numdata, &Cs) < 0 || Cs <= 0) {
		    opserr << "WARNING invalid cs\n";
		    ifNoError = false;
		}

	    } else {
	
		opserr << "WARNING invalid optional arguments \n";
		ifNoError = false;
		break;
	
	    }
	}

    } //end input

  
    //if error detected
    if (!ifNoError) {
	//input:
	//want:
	opserr << "Want: element YamamotoBiaxialHDR eleTag? iNode? jNode? Tp? DDo? DDi? Hr?  <-coRS cr? cs?> <-orient <x1? x2? x3?> y1? y2? y3?> <-mass m?>\n";
	return 0;
    }
  

    // now create the YamamotoBiaxialHDR
    theElement = new YamamotoBiaxialHDR(eleTag, iNode, jNode, Tp, DDo, DDi, Hr, Cr, Cs, oriYp, oriX, mass);
  
    // if get here we have successfully created the YamamotoBiaxialHDR and added it to the domain
    return theElement;
}


// --- end of OPS_YamamotoBiaxialHDR





// initialize the class wide variables
Matrix YamamotoBiaxialHDR::theMatrix(12,12);
Vector YamamotoBiaxialHDR::theVector(12);
Vector YamamotoBiaxialHDR::theLoad(12);


YamamotoBiaxialHDR::YamamotoBiaxialHDR(int Tag, int Nd1, int Nd2, int Tp, double DDo, double DDi, double Hr, double Cr, double Cs,
	   const Vector OriYp, const Vector OriX, double Mass)
  : Element(Tag, ELE_TAG_YamamotoBiaxialHDR),
    tp(Tp),ddo(DDo),ddi(DDi),hr(Hr),cr(Cr),cs(Cs),
    connectedExternalNodes(2),
    oriX(OriX), oriYp(OriYp), mass(Mass),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6)
{
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "YamamotoBiaxialHDR::setUp() - element: "
	   << this->getTag() << " failed to create an ID of size 2\n";
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  // cross-section area
  ar = ( ddo*ddo - ddi*ddi )*M_PI / 4.0;

  // polar moment of inertia of area
  ip = ( pow(ddo,4.0) - pow(ddi,4.0) )*M_PI / 32.0;


  // Bridgestone X0.6R
  if ( tp == 1 ) {
    for (int i=0; i<2; i++) {
      initialStiff[i] = (0.22*cr + 1.0*cs) * 1.0e6 * ar/hr;
    }
    alpha = 0.25*hr;
    nn    = 0.7;
  }


  // initial basic stiffness matrix
  basicStiffInit.Zero();

  basicStiffInit(1,1) = this->getInitialTangent(0);
  basicStiffInit(2,2) = this->getInitialTangent(1);


  // initialize variables
  this->revertToStart();

  //check print
  opserr << "basicStiffInit:  " << basicStiff << endln;


  // check print--------------------------------------------
  //opserr << "   \n";
  //opserr << "YamamotoBiaxialHDR::YamamotoBiaxialHDR()\n";
  //opserr << "  tp = " << tp << endln;
  //opserr << "  dr = " << dr << endln;
  //opserr << "  hr = " << hr << endln;
  //opserr << "  ar = " << ar << endln;
  //opserr << "  ip = " << ip << endln;
  //opserr << "  alpha = " << alpha << endln;
  //opserr << "  nn = " << nn << endln;
  //------------------------------------------------------

}

YamamotoBiaxialHDR::YamamotoBiaxialHDR()
  : Element(0, ELE_TAG_YamamotoBiaxialHDR),
    tp(0),ddo(0.0),ddi(0.0),hr(0.0),cr(0.0),cs(0.0),
    connectedExternalNodes(2),
    oriX(0), oriYp(0), mass(0.0),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6)
{	

  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "YamamotoBiaxialHDR::YamamotoBiaxialHDR() - "
	   <<  "failed to create an ID of size 2\n";
    exit(-1);
  }
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;


  //
  for (int i=0; i<2; i++) {

    trialDeform[i]  = 0.0;
    trialForce[i]   = 0.0;
    trialStiff[i]   = 0.0;
    trialQ[i]   = 0.0;
    trialP[i]   = 0.0;

    commitDeform[i]  = 0.0;
    commitForce[i]   = 0.0;
    commitStiff[i]   = 0.0;
    commitQ[i]   = 0.0;
    commitP[i]   = 0.0;

  }


}


YamamotoBiaxialHDR::~YamamotoBiaxialHDR()
{
  // does nothing
}


int YamamotoBiaxialHDR::getNumExternalNodes() const
{
  return 2;
}


const ID& YamamotoBiaxialHDR::getExternalNodes() 
{
  return connectedExternalNodes;
}


Node** YamamotoBiaxialHDR::getNodePtrs() 
{
  return theNodes;
}


int YamamotoBiaxialHDR::getNumDOF() 
{
  return 12;
}


void YamamotoBiaxialHDR::setDomain(Domain *theDomain)
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
      opserr << "WARNING YamamotoBiaxialHDR::setDomain() - Nd1: " 
	     << connectedExternalNodes(0) << " does not exist in the model for ";
    } else  {
      opserr << "WARNING YamamotoBiaxialHDR::setDomain() - Nd2: " 
	     << connectedExternalNodes(1) << " does not exist in the model for ";
    }
    opserr << "YamamotoBiaxialHDR ele: " << this->getTag() << endln;
    
    return;
  }
  
  // now determine the number of dof and the dimension    
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  
  // if differing dof at the ends - print a warning message
  if (dofNd1 != 6)  {
    opserr << "YamamotoBiaxialHDR::setDomain() - node 1: "
	   << connectedExternalNodes(0) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  if (dofNd2 != 6)  {
    opserr << "YamamotoBiaxialHDR::setDomain() - node 2: "
	   << connectedExternalNodes(1) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  
  // call the base class method
  this->DomainComponent::setDomain(theDomain);
  
  // set up the transformation matrix for orientation
  this->setUp();
}   	 

//--------------------------------------------------------------------------------------
int YamamotoBiaxialHDR::setTrialStrain(const Vector &strain)
{

  double tstrain[2], gamma; // shear strain
  double taur,taus;

  trialDeform[0] = strain(1); // local-y
  trialDeform[1] = strain(2); // local-z

  trialP[0] = trialDeform[0];
  trialP[1] = trialDeform[1];

  tstrain[0] = trialDeform[0]/hr;
  tstrain[1] = trialDeform[1]/hr;


  // check print--------------------------------------------
  //opserr << "   \n";
  //opserr << "YamamotoBiaxialHDR::setTrialStrain()\n";
  //opserr << "  trialDeform[0] = " << trialDeform[0] << endln;
  //opserr << "  trialDeform[1] = " << trialDeform[1] << endln;
  //opserr << "  trialP[0] = " << trialP[0] << endln;
  //opserr << "  trialP[1] = " << trialP[1] << endln;
  //------------------------------------------------------


  if ( tp == 1 ) {
  //-----------------------------------------------------------------------------
  // Bridgestone X0.6R
  //-----------------------------------------------------------------------------

    double DP[2], DPl;
    double Fr, Fs;

      DP[0] = trialP[0] - commitP[0];
      DP[1] = trialP[1] - commitP[1];
      DPl = sqrt(pow(DP[0],2) + pow(DP[1],2));

      double commitQl = sqrt(pow(commitQ[0],2)+pow(commitQ[1],2));

      if (DPl < DBL_EPSILON) {
        trialQ[0] = commitQ[0];
        trialQ[1] = commitQ[1];
      }

      else if ( (DPl >= DBL_EPSILON) && (commitQl < DBL_EPSILON) ) {
        trialQ[0] = commitQ[0] + DP[0]/alpha;
        trialQ[1] = commitQ[1] + DP[1]/alpha;
      }

      else {
        trialQ[0] = commitQ[0] + DPl/alpha * (DP[0]/DPl - pow(commitQl,nn)*commitQ[0]/commitQl);
        trialQ[1] = commitQ[1] + DPl/alpha * (DP[1]/DPl - pow(commitQl,nn)*commitQ[1]/commitQl);
      }

    gamma =  sqrt(pow(tstrain[0],2) + pow(tstrain[1],2));

    // Fr
    if (gamma<1.8) {
      taur = 0.22*gamma;  // [MPa]
    } else {
      taur = 0.22*gamma + 0.20*pow((gamma-1.8),2);  // [MPa]
    }
    Fr = cr*taur*ar*1e6;  // [MPa] --> [N]

    if (Fr < DBL_EPSILON) {
      trialFr[0] = 0.0;
      trialFr[1] = 0.0;
    } else {
      trialFr[0] =   Fr    * tstrain[0]/gamma;
      trialFr[1] =   Fr    * tstrain[1]/gamma;
    }

    // Fs
    taus = 0.25 + 0.02*gamma + 0.016*pow(gamma,3);   // [MPa]
    Fs = cs*taus*ar*1e6;  // [MPa] --> [N]

    double trialQl = sqrt(pow(trialQ[0],2)+pow(trialQ[1],2));
    if (trialQl < DBL_EPSILON) {
      trialFs[0] = 0.0;
      trialFs[1] = 0.0;
    } else {
      trialFs[0] =   Fs    * trialQ[0];
      trialFs[1] =   Fs    * trialQ[1];
    }

    trialForce[0] = trialFr[0] + trialFs[0];
    trialForce[1] = trialFr[1] + trialFs[1];

    // check print--------------------------------------------
    //opserr << "  trialFr[0] = " << trialFr[0] << endln;
    //opserr << "  trialFr[1] = " << trialFr[1] << endln;
    //opserr << "  trialFs[0] = " << trialFs[0] << endln;
    //opserr << "  trialFs[1] = " << trialFs[1] << endln;
    //opserr << "  trialForce[0] = " << trialForce[0] << endln;
    //opserr << "  trialForce[1] = " << trialForce[1] << endln;
    //------------------------------------------------------

  }  //   if ( tp == 1 ) ... end


  // tangent stiffness
    for (int i=0; i<2; i++) {
      if ((trialDeform[i] - commitDeform[i]) < DBL_EPSILON) {
       //trialStiff[i] = initialStiff[i];
       //trialStiff[i] = commitStiff[i];
       trialStiff[i] = initialStiff[i];
      } else {
       // trialStiff[i] = (trialForce[i]-commitForce[i])/(trialDeform[i]-commitDeform[i]);
       trialStiff[i] = initialStiff[i];
      }

    }

    // check print--------------------------------------------
    //opserr << "  trialStiff[0] = " << trialStiff[0] << endln;
    //opserr << "  trialStiff[1] = " << trialStiff[1] << endln;
    //------------------------------------------------------


  return 0;
}

  const double& YamamotoBiaxialHDR::getStrain(int direction)
{
  return trialDeform[direction];
}

  const double& YamamotoBiaxialHDR::getStress(int direction)
{
  return trialForce[direction];
}

  const double& YamamotoBiaxialHDR::getTangent(int direction)
{
  return trialStiff[direction];
}

  const double& YamamotoBiaxialHDR::getInitialTangent(int direction)
{
  return initialStiff[direction];
}

//--------------------------------------------------------------------------------------
int YamamotoBiaxialHDR::commitState()
{
  int errCode = 0;
  

  // commit <- trial

  for (int i=0; i<2; i++) {
    commitDeform[i]  = trialDeform[i];
    commitForce[i]   = trialForce[i];
    commitStiff[i]   = trialStiff[i];
    commitQ[i] = trialQ[i];
    commitP[i] = trialP[i];
    commitFr[i] = trialFr[i];
    commitFs[i] = trialFs[i];
  }

  return errCode;
}


int YamamotoBiaxialHDR::revertToLastCommit()
{
  int errCode = 0;

  // trial <- commit
  for (int i=0; i<2; i++) {
    trialDeform[i]  = commitDeform[i];
    trialForce[i]   = commitForce[i];
    trialStiff[i]   = commitStiff[i];
    trialQ[i] = commitQ[i];
    trialP[i] = commitP[i];
    trialFr[i] = commitFr[i];
    trialFs[i] = commitFs[i];
  }

  return errCode;
}


int YamamotoBiaxialHDR::revertToStart()
{
  int errCode=0;

  // trial variables
  basicDisp.Zero();
  basicForce.Zero();

  // reset basic stiffness matrix
  basicStiff = basicStiffInit;

  for (int i=0; i<2; i++) {

    trialDeform[i]  = 0.0;
    trialForce[i]   = 0.0;
    trialStiff[i]   = initialStiff[i];
    trialQ[i] = 0.0;
    trialP[i]  = 0.0;
    trialFr[i] = 0.0;
    trialFs[i] = 0.0;

    commitDeform[i] = 0.0;
    commitForce[i] = 0.0;
    commitStiff[i] = initialStiff[i];
    commitQ[i] = 0.0;
    commitP[i]  = 0.0;
    commitFr[i] = 0.0;
    commitFs[i] = 0.0;
  }

  return errCode;
}


int YamamotoBiaxialHDR::update()
{
  // get global trial displacements and velocities
  const Vector &dsp1 = theNodes[0]->getTrialDisp();
  const Vector &dsp2 = theNodes[1]->getTrialDisp();
  
  static Vector globalDisp(12), globalDispDot(12);
  for (int i=0; i<6; i++)  {
    globalDisp(i)   = dsp1(i);
    globalDisp(i+6) = dsp2(i);
  }

  static Vector localDispDot(12);

  
  // transform response from the global to the local system
  localDisp    = Tgl*globalDisp;
  
  // transform response from the local to the basic system
  basicDisp    = Tlb*localDisp;


  // calculate shear forces and stiffnesses in basic y- and z-direction
  // get trial shear forces of hysteretic component
  basicForce.Zero();
  basicStiff.Zero();

  this->setTrialStrain(basicDisp);

  basicForce(1) = this->getStress(0);
  basicForce(2) = this->getStress(1);

  basicStiff(1,1) = this->getTangent(0);
  basicStiff(2,2) = this->getTangent(1);

  //check print
  //opserr << "basicStiff:  " << basicStiff << endln;

  return 0;
}


const Matrix& YamamotoBiaxialHDR::getTangentStiff()
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


const Matrix& YamamotoBiaxialHDR::getInitialStiff()
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


const Matrix& YamamotoBiaxialHDR::getMass()
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


void YamamotoBiaxialHDR::zeroLoad()
{
  theLoad.Zero();
}


int YamamotoBiaxialHDR::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
  opserr <<"YamamotoBiaxialHDR::addLoad() - "
	 << "load type unknown for element: "
	 << this->getTag() << endln;
  
  return -1;
}


int YamamotoBiaxialHDR::addInertiaLoadToUnbalance(const Vector &accel)
{
  // check for quick return
  if (mass == 0.0)  {
    return 0;
  }
  
  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
    opserr << "YamamotoBiaxialHDR::addInertiaLoadToUnbalance() - "
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


const Vector& YamamotoBiaxialHDR::getResistingForce()
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


const Vector& YamamotoBiaxialHDR::getResistingForceIncInertia()
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


int YamamotoBiaxialHDR::sendSelf(int commitTag, Channel &sChannel)
{
  return -1;
}


int YamamotoBiaxialHDR::recvSelf(int commitTag, Channel &rChannel,
				  FEM_ObjectBroker &theBroker)
{
  return -1;
}


int YamamotoBiaxialHDR::displaySelf(Renderer &theViewer,
				    int displayMode, float fact,
				    const char **modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}


void YamamotoBiaxialHDR::Print(OPS_Stream &s, int flag)
{
  
  if (flag == 0)  {
    // print everything
    s << "Element: " << this->getTag(); 
    s << "  type: YamamotoBiaxialHDR  iNode: " << connectedExternalNodes(0);
    s << "                            jNode: " << connectedExternalNodes(1) << endln;
    
    s << "Input parameters: " << endln;
    s << "  Tp: " << tp << endln;
    s << "  DDo: " << ddo << endln;
    s << "  DDi: " << ddi << endln;
    s << "  Hr: " << hr << endln;
    s << "  Cr: " << cr << endln;
    s << "  Cs: " << cs << endln;
    
    // determine resisting forces in global system
    //s << "  resisting force: " << this->getResistingForce() << endln;
    //s << "  Disp-x: " << this->getStrain(0) << endln;
    //s << "  Force-x: " << this->getStress(0) << endln;
    //s << "  Stiff-x: " << this->getTangent(0) << endln;
    //s << "  StiffInit-x: " << this->getInitialTangent(0) << endln;
    //s << "  Disp-y: " << this->getStrain(1) << endln;
    //s << "  Force-y: " << this->getStress(1) << endln;
    //s << "  Stiff-y: " << this->getTangent(1) << endln;
    //s << "  StiffInit-y: " << this->getInitialTangent(1) << endln;

  } else if (flag == 1)  {
    // does nothing
  }
}


Response* YamamotoBiaxialHDR::setResponse(const char **argv, int argc,
					   OPS_Stream &output)
{
  Response *theResponse = 0;
  
  output.tag("ElementOutput");
  output.attr("eleType","YamamotoBiaxialHDR");
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


int YamamotoBiaxialHDR::getResponse(int responseID, Information &eleInfo)
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
void YamamotoBiaxialHDR::setUp()
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
      opserr << "WARNING YamamotoBiaxialHDR::setUp() - " 
	     << "element: " << this->getTag() << endln
	     << "ignoring nodes and using specified "
	     << "local x vector to determine orientation\n";
    }
  }
  // check that vectors for orientation are of correct size
  if (oriX.Size() != 3 || oriYp.Size() != 3)  {
    opserr << "YamamotoBiaxialHDR::setUp() - "
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
    opserr << "YamamotoBiaxialHDR::setUp() - "
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

