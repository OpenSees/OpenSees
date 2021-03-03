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
                                                                        
// $Revision$
// $Date$
// $URL$
                                                                        
// Written: MHS 
// Created: May 2001
//
// Description: This file contains the class implementation for CorotTruss.

#include <CorotTruss.h>
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <Parameter.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <map>

#include <ElementResponse.h>

Matrix CorotTruss::M2(2,2);
Matrix CorotTruss::M4(4,4);
Matrix CorotTruss::M6(6,6);
Matrix CorotTruss::M12(12,12);

Vector CorotTruss::V2(2);
Vector CorotTruss::V4(4);
Vector CorotTruss::V6(6);
Vector CorotTruss::V12(12);


#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_CorotTrussElement()
{
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 4) {
    opserr << "Invalid Args want: element CorotTruss $tag $iNode $jNode $sectTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>";
    opserr << " or: element CorotTruss $tag $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;	
  }

  if (numRemainingArgs == 4 || numRemainingArgs == 6 || numRemainingArgs == 8 || numRemainingArgs == 10)
    return 0; // it's a CorotTrussSection

  int iData[3];
  double A = 0.0;
  double rho = 0.0;
  int matTag = 0;
  int doRayleigh = 0; // by default rayleigh not done
  int cMass = 0; // by default use lumped mass matrix
  int ndm = OPS_GetNDM();

  int numData = 3;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer (tag, iNode, jNode) in element CorotTruss " << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDouble(&numData, &A) != 0) {
    opserr << "WARNING: Invalid A: element CorotTruss " << iData[0] << 
      " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;	
  }

  numData = 1;
  if (OPS_GetInt(&numData, &matTag) != 0) {
    opserr << "WARNING: Invalid matTag: element CorotTruss " << iData[0] << 
      " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial = OPS_GetUniaxialMaterial(matTag);
    
  if (theUniaxialMaterial == 0) {
    opserr << "WARNING: Invalid material not found element CorotTruss " << iData[0] << " $iNode $jNode $A " << 
      matTag << " <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;
  }
  
  numRemainingArgs -= 5;
  while (numRemainingArgs > 1) {
    const char *argvS = OPS_GetString();
  
    if (strcmp(argvS,"-rho") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, &rho) != 0) {
	opserr << "WARNING Invalid rho in element CorotTruss " << iData[0] << 
	  " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
	return 0;
      }
    } else if (strcmp(argvS,"-cMass") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &cMass) != 0) {
	opserr << "WARNING: Invalid cMass in element CorotTruss " << iData[0] << 
	  " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
	return 0;
      }
    } else if (strcmp(argvS,"-doRayleigh") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &doRayleigh) != 0) {
	opserr << "WARNING: Invalid doRayleigh in element CorotTruss " << iData[0] << 
	  " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
	return 0;
      }
    } else {
      opserr << "WARNING: Invalid option " << argvS << "  in: element CorotTruss " << iData[0] << 
	" $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
      return 0;
    }      
    numRemainingArgs -= 2;
  }

  // now create the CorotTruss
  theElement = new CorotTruss(iData[0], ndm, iData[1], iData[2], *theUniaxialMaterial, A, rho, doRayleigh, cMass);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element CorotTruss " << iData[0] << 
      " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
  }

  return theElement;
}

void *
OPS_CorotTrussElement(const ID& info)
{
    if (info.Size() == 0) return 0;

    Element *theElement = 0;

    int iData[3];

    double A = 0.0;
    double rho = 0.0;
    int matTag = 0;
    int doRayleigh = 0; // by default rayleigh not done
    int cMass = 0; // by default use lumped mass matrix
    int ndm = OPS_GetNDM();

    static std::map<int, Vector> meshdata;
    if (info(0) == 1) {

        // save data
        int numRemainingArgs = OPS_GetNumRemainingInputArgs();

        if (numRemainingArgs < 2) {
            opserr << "Invalid Args want: element CorotTruss $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
            return 0;
        }

        int numData = 1;
        if (OPS_GetDouble(&numData, &A) != 0) {
            opserr << "WARNING: Invalid A: element CorotTruss " <<
                   " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
            return 0;
        }

        numData = 1;
        if (OPS_GetInt(&numData, &matTag) != 0) {
            opserr << "WARNING: Invalid matTag: element CorotTruss "  <<
                   " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
            return 0;
        }

        while (OPS_GetNumRemainingInputArgs() > 1) {
            const char *argvS = OPS_GetString();

            if (strcmp(argvS,"-rho") == 0) {
                numData = 1;
                if (OPS_GetDouble(&numData, &rho) != 0) {
                    opserr << "WARNING Invalid rho in element CorotTruss " <<
                           " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
                    return 0;
                }
            } else if (strcmp(argvS,"-cMass") == 0) {
                numData = 1;
                if (OPS_GetInt(&numData, &cMass) != 0) {
                    opserr << "WARNING: Invalid cMass in element CorotTruss " <<
                           " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
                    return 0;
                }
            } else if (strcmp(argvS,"-doRayleigh") == 0) {
                numData = 1;
                if (OPS_GetInt(&numData, &doRayleigh) != 0) {
                    opserr << "WARNING: Invalid doRayleigh in element CorotTruss " <<
                           " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
                    return 0;
                }
            } else {
                opserr << "WARNING: Invalid option " << argvS << "  in: element CorotTruss " <<
                       " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
                return 0;
            }

        }

        if (info.Size() < 2) {
            opserr << "WARNING: need info -- inmesh, meshtag\n";
            return 0;
        }

        // save the data for a mesh
        Vector& mdata = meshdata[info(1)];
        mdata.resize(5);
        mdata(0) = A;
        mdata(1) = rho;
        mdata(2) = (double) matTag;
        mdata(3) = (double) doRayleigh;
        mdata(4) = (double) cMass;
        return &meshdata;

    } else if (info(0) == 2) {

        if (info.Size() < 5) {
            opserr << "WARNING: need info -- inmesh, meshtag, eleTag, nd1, nd2\n";
            return 0;
        }

        // get the data for a mesh
        Vector& mdata = meshdata[info(1)];
        if (mdata.Size() < 5) return 0;

        iData[0] = info(2);
        iData[1] = info(3);
        iData[2] = info(4);

        A = mdata(0);
        rho = mdata(1);
        matTag = (int) mdata(2);
        doRayleigh = (int) mdata(3);
        cMass = (int) mdata(4);
    }

    UniaxialMaterial *theUniaxialMaterial = OPS_GetUniaxialMaterial(matTag);

    if (theUniaxialMaterial == 0) {
        opserr << "WARNING: Invalid material not found element CorotTruss " << iData[0] << " $iNode $jNode $A " <<
               matTag << " <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
        return 0;
    }

    // now create the CorotTruss
    theElement = new CorotTruss(iData[0], ndm, iData[1], iData[2], *theUniaxialMaterial, A, rho, doRayleigh, cMass);

    if (theElement == 0) {
        opserr << "WARNING: out of memory: element CorotTruss " << iData[0] <<
               " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    }

    return theElement;
}


// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the CorotTruss end nodes.
CorotTruss::CorotTruss(int tag, int dim,
			   int Nd1, int Nd2, 
			   UniaxialMaterial &theMat,
			   double a, double r,
               int damp, int cm)
  :Element(tag,ELE_TAG_CorotTruss),     
  theMaterial(0), connectedExternalNodes(2),
  numDOF(0), numDIM(dim), Lo(0.0), Ln(0.0), 
  A(a), rho(r), doRayleighDamping(damp), cMass(cm),
  R(3,3), theLoad(0), theMatrix(0), theVector(0)
{
  // get a copy of the material and check we obtained a valid copy
  theMaterial = theMat.getCopy();
  if (theMaterial == 0) {
    opserr << "FATAL CorotTruss::CorotTruss - " <<  tag <<
      "failed to get a copy of material with tag " << theMat.getTag() << endln;
    exit(-1);
  }
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2) {
    opserr << "FATAL CorotTruss::CorotTruss - " <<  tag <<
      "failed to create an ID of size 2\n";
    exit(-1);
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;        

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
CorotTruss::CorotTruss()
  :Element(0,ELE_TAG_CorotTruss),     
  theMaterial(0),connectedExternalNodes(2),
  numDOF(0), numDIM(0), Lo(0.0), Ln(0.0),
  A(0.0), rho(0.0), doRayleighDamping(0), cMass(0),
  R(3,3), theLoad(0), theMatrix(0), theVector(0)
{
  // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2) {
    opserr << "FATAL CorotTruss::CorotTruss - failed to create an ID of size 2\n";
    exit(-1);
  }

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
CorotTruss::~CorotTruss()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    if (theMaterial != 0)
        delete theMaterial;
    if (theLoad != 0)
        delete theLoad;
}

int
CorotTruss::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
CorotTruss::getExternalNodes(void) 
{
	return connectedExternalNodes;
}

Node **
CorotTruss::getNodePtrs(void) 
{
  return theNodes;
}

int
CorotTruss::getNumDOF(void) 
{
	return numDOF;
}

// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the CorotTruss element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
CorotTruss::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    Lo = 0.0;
    Ln = 0.0;
    return;
  }
  
  // first set the node pointers
  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);	
  
  // if can't find both - send a warning message
  if ((theNodes[0] == 0) || (theNodes[1] == 0)) {
    opserr << "CorotTruss::setDomain() - CorotTruss " << this->getTag() << " node " <<
      Nd1 << "does not exist in the model \n";
    
    // fill this in so don't segment fault later
    numDOF = 6;    
    
    return;
  }
  
  // now determine the number of dof and the dimesnion    
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  
  // if differing dof at the ends - print a warning message
  if (dofNd1 != dofNd2) {
    opserr << "WARNING CorotTruss::setDomain(): nodes " << Nd1 <<
      " and " << Nd2 << "have differing dof at ends for CorotTruss " << this->getTag() << endln;
    
    // fill this in so don't segment fault later
    numDOF = 6;    
    
    return;
  }	
  
  if (numDIM == 1 && dofNd1 == 1) {
    numDOF = 2;
    theMatrix = &M2;
    theVector = &V2;
  }
  else if (numDIM == 2 && dofNd1 == 2) {
    numDOF = 4;
    theMatrix = &M4;
    theVector = &V4;
  }
  else if (numDIM == 2 && dofNd1 == 3) {
    numDOF = 6;
    theMatrix = &M6;
    theVector = &V6;
  }
  else if (numDIM == 3 && dofNd1 == 3) {
    numDOF = 6;
    theMatrix = &M6;
    theVector = &V6;
  }
  else if (numDIM == 3 && dofNd1 == 6) {
    numDOF = 12;
    theMatrix = &M12;
    theVector = &V12;
  }
  else {
    opserr << " CorotTruss::setDomain -- nodal DOF " << dofNd1 << " not compatible with element\n";
    
    // fill this in so don't segment fault later
    numDOF = 6;    
    
    return;
  }

    // create the load vector
    if (theLoad == 0)
        theLoad = new Vector(numDOF);
    else if (theLoad->Size() != numDOF) {
        delete theLoad;
        theLoad = new Vector(numDOF);
    }
    if (theLoad == 0) {
        opserr << "Truss::setDomain - truss " << this->getTag()
            << "out of memory creating vector of size" << numDOF << endln;
        exit(-1);
        return;
    }          
    
	// call the base class method
	this->DomainComponent::setDomain(theDomain);

	// now determine the length, cosines and fill in the transformation
	// NOTE t = -t(every one else uses for residual calc)
	const Vector &end1Crd = theNodes[0]->getCrds();
	const Vector &end2Crd = theNodes[1]->getCrds();

	// Determine global offsets
    double cosX[3];
    cosX[0] = 0.0;  cosX[1] = 0.0;  cosX[2] = 0.0;
    int i;
    for (i = 0; i < numDIM; i++) {
        cosX[i] += end2Crd(i)-end1Crd(i);
    }

	// Set undeformed and initial length
	Lo = cosX[0]*cosX[0] + cosX[1]*cosX[1] + cosX[2]*cosX[2];
	Lo = sqrt(Lo);
	Ln = Lo;

    // Initial offsets
   	d21[0] = Lo;
	d21[1] = 0.0;
	d21[2] = 0.0;

	// Set global orientation
	cosX[0] /= Lo;
	cosX[1] /= Lo;
	cosX[2] /= Lo;

	R(0,0) = cosX[0];
	R(0,1) = cosX[1];
	R(0,2) = cosX[2];

	// Element lies outside the YZ plane
	if (fabs(cosX[0]) > 0.0) {
		R(1,0) = -cosX[1];
		R(1,1) =  cosX[0];
		R(1,2) =  0.0;

		R(2,0) = -cosX[0]*cosX[2];
		R(2,1) = -cosX[1]*cosX[2];
		R(2,2) =  cosX[0]*cosX[0] + cosX[1]*cosX[1];
	}
	// Element is in the YZ plane
	else {
		R(1,0) =  0.0;
		R(1,1) = -cosX[2];
		R(1,2) =  cosX[1];

		R(2,0) =  1.0;
		R(2,1) =  0.0;
		R(2,2) =  0.0;
	}

	// Orthonormalize last two rows of R
	double norm;
	for (i = 1; i < 3; i++) {
		norm = sqrt(R(i,0)*R(i,0) + R(i,1)*R(i,1) + R(i,2)*R(i,2));
		R(i,0) /= norm;
		R(i,1) /= norm;
		R(i,2) /= norm;
	}
}

int
CorotTruss::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "CorotTruss::commitState () - failed in base class\n";
  }    
  retVal = theMaterial->commitState();
  return retVal;
}

int
CorotTruss::revertToLastCommit()
{
	// Revert the material
	return theMaterial->revertToLastCommit();
}

int
CorotTruss::revertToStart()
{
	// Revert the material to start
	return theMaterial->revertToStart();
}

int
CorotTruss::update(void)
{
  // Nodal displacements
  const Vector &end1Disp  = theNodes[0]->getTrialDisp();
  const Vector &end2Disp  = theNodes[1]->getTrialDisp();    
  const Vector &end1Vel   = theNodes[0]->getTrialVel();
  const Vector &end2Vel   = theNodes[1]->getTrialVel();	
  
  // Initial offsets
  d21[0] = Lo; d21[1] = d21[2] = 0.0;
  v21[0] = v21[1] = v21[2] = 0.0;
  
  // Update offsets in basic system due to nodal displacements
  for (int i = 0; i < numDIM; i++) {
    double deltaDisp = end2Disp(i) - end1Disp(i);
    d21[0] += deltaDisp*R(0,i);
    d21[1] += deltaDisp*R(1,i);
    d21[2] += deltaDisp*R(2,i);
    double deltaVel = end2Vel(i) - end1Vel(i);
    v21[0] += deltaVel*R(0,i);
    v21[1] += deltaVel*R(1,i);
    v21[2] += deltaVel*R(2,i);
  }
  
  // Compute new length
  Ln = sqrt(d21[0]*d21[0] + d21[1]*d21[1] + d21[2]*d21[2]);
  
  // Compute engineering strain and strain rate
  double strain = (Ln - Lo)/Lo;
  double rate = (d21[0]*v21[0] + d21[1]*v21[1] + d21[2]*v21[2])/Ln/Lo;
  
  // Set material trial strain
  return theMaterial->setTrialStrain(strain,rate);
}

const Matrix &
CorotTruss::getTangentStiff(void)
{
    static Matrix kl(3,3);

    // Material stiffness
    //
    // Get material tangent
    double EA = A*theMaterial->getTangent();
    EA /= (Ln * Ln * Lo);

    int i,j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            kl(i,j) = EA*d21[i]*d21[j];

    // Geometric stiffness
    //
    // Get material stress
    double q = A*theMaterial->getStress();
    double SA = q/(Ln*Ln*Ln);
    double SL = q/Ln;
    
    for (i = 0; i < 3; i++) {
        kl(i,i) += SL;
        for (j = 0; j < 3; j++)
            kl(i,j) -= SA*d21[i]*d21[j];
    }
    
    // Compute R'*kl*R
    static Matrix kg(3,3);
    kg.addMatrixTripleProduct(0.0, R, kl, 1.0);

    Matrix &K = *theMatrix;
    K.Zero();

    // Copy stiffness into appropriate blocks in element stiffness
    int numDOF2 = numDOF/2;
    for (i = 0; i < numDIM; i++) {
        for (j = 0; j < numDIM; j++) {
            K(i,j)                 =  kg(i,j);
            K(i,j+numDOF2)         = -kg(i,j);
            K(i+numDOF2,j)         = -kg(i,j);
            K(i+numDOF2,j+numDOF2) =  kg(i,j);
        }
    }

    return *theMatrix;
}


const Matrix &
CorotTruss::getInitialStiff(void)
{
    static Matrix kl(3,3);

    // Material stiffness
    kl.Zero();
    kl(0,0) = A * theMaterial->getInitialTangent() / Lo;

    // Compute R'*kl*R
    static Matrix kg(3,3);
    kg.addMatrixTripleProduct(0.0, R, kl, 1.0);

    Matrix &K = *theMatrix;
    K.Zero();

    // Copy stiffness into appropriate blocks in element stiffness
    int numDOF2 = numDOF/2;
    for (int i = 0; i < numDIM; i++) {
        for (int j = 0; j < numDIM; j++) {
            K(i,j)                 =  kg(i,j);
            K(i,j+numDOF2)         = -kg(i,j);
            K(i+numDOF2,j)         = -kg(i,j);
            K(i+numDOF2,j+numDOF2) =  kg(i,j);
        }
    }

    return *theMatrix;
}


const Matrix &
CorotTruss::getDamp(void)
{
    static Matrix kl(3,3);

    Matrix a(3,1);
    a(0,0) = (Lo+d21[0])/Ln;
    a(1,0) = d21[1]/Ln;
    a(2,0) = 0.0;

    Matrix cb(1,1);
    cb(0,0) = A*theMaterial->getDampTangent()/Lo;

    kl.addMatrixTripleProduct(0.0, a, cb, 1.0);

    // Compute R'*kl*R
    static Matrix kg(3,3);
    kg.addMatrixTripleProduct(0.0, R, kl, 1.0);

    Matrix &K = *theMatrix;
    K.Zero();

    if (doRayleighDamping == 1)
        *theMatrix = this->Element::getDamp();

    // Copy stiffness into appropriate blocks in element stiffness
    int numDOF2 = numDOF/2;
    for (int i = 0; i < numDIM; i++) {
        for (int j = 0; j < numDIM; j++) {
            K(i,j)                 +=  kg(i,j);
            K(i,j+numDOF2)         += -kg(i,j);
            K(i+numDOF2,j)         += -kg(i,j);
            K(i+numDOF2,j+numDOF2) +=  kg(i,j);
        }
    }

    return *theMatrix;
}


const Matrix &
CorotTruss::getMass(void)
{
    // zero the matrix
    Matrix &mass = *theMatrix;
    mass.Zero();
    
    // check for quick return
    if (Lo == 0.0 || rho == 0.0)
        return mass;
    
    if (cMass == 0)  {
        // lumped mass matrix
        double m = 0.5*rho*Lo;
        int numDOF2 = numDOF/2;
        for (int i = 0; i < numDIM; i++) {
            mass(i,i) = m;
            mass(i+numDOF2,i+numDOF2) = m;
        }
    } else  {
        // consistent mass matrix
        double m = rho*Lo/6.0;
        int numDOF2 = numDOF/2;
        for (int i = 0; i < numDIM; i++) {
            mass(i,i) = 2.0*m;
            mass(i,i+numDOF2) = m;
            mass(i+numDOF2,i) = m;
            mass(i+numDOF2,i+numDOF2) = 2.0*m;
        }
    }
    
    return *theMatrix;
}


void 
CorotTruss::zeroLoad(void)
{
    theLoad->Zero();
}

int 
CorotTruss::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "CorotTruss::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  
  return -1;
}



int 
CorotTruss::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for quick return
    if (Lo == 0.0 || rho == 0.0)
        return 0;
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);    
    
    int nodalDOF = numDOF/2;
    
    // want to add ( - fact * M R * accel ) to unbalance
    if (cMass == 0)  {
        double m = 0.5*rho*Lo;
        for (int i=0; i<numDIM; i++) {
            (*theLoad)(i) -= m*Raccel1(i);
            (*theLoad)(i+nodalDOF) -= m*Raccel2(i);
        }
    } else  {
        double m = rho*Lo/6.0;
        for (int i=0; i<numDIM; i++) {
            (*theLoad)(i) -= 2.0*m*Raccel1(i) + m*Raccel2(i);
            (*theLoad)(i+nodalDOF) -= m*Raccel1(i) + 2.0*m*Raccel2(i);
        }
    }
    
    return 0;
}

const Vector &
CorotTruss::getResistingForce()
{
	// Get material stress
	double SA = A*theMaterial->getStress();
	SA /= Ln;

    static Vector ql(3);

	ql(0) = d21[0]*SA;
	ql(1) = d21[1]*SA;
	ql(2) = d21[2]*SA;

    static Vector qg(3);
    qg.addMatrixTransposeVector(0.0, R, ql, 1.0);

    Vector &P = *theVector;
    P.Zero();

    // Copy forces into appropriate places
    int numDOF2 = numDOF/2;
    for (int i = 0; i < numDIM; i++) {
        P(i)         = -qg(i);
        P(i+numDOF2) =  qg(i);
    }
    
    return *theVector;
}



const Vector &
CorotTruss::getResistingForceIncInertia()
{	
    Vector &P = *theVector;
    P = this->getResistingForce();
    
    // subtract external load
    P -= *theLoad;
    
    // now include the mass portion
    if (Lo != 0.0 && rho != 0.0) {
        
        // add inertia forces from element mass
        const Vector &accel1 = theNodes[0]->getTrialAccel();
        const Vector &accel2 = theNodes[1]->getTrialAccel();	
        
        int numDOF2 = numDOF/2;
        
        if (cMass == 0)  {
            // lumped mass matrix
            double m = 0.5*rho*Lo;
            for (int i=0; i<numDIM; i++) {
                P(i) += m*accel1(i);
                P(i+numDOF2) += m*accel2(i);
            }
        } else  {
            // consistent mass matrix
            double m = rho*Lo/6.0;
            for (int i=0; i<numDIM; i++) {
                (*theVector)(i) += 2.0*m*accel1(i) + m*accel2(i);
                (*theVector)(i+numDOF2) += m*accel1(i) + 2.0*m*accel2(i);
            }
        }
        
        // add the damping forces if rayleigh damping
        if (doRayleighDamping == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
            theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);
    } else  {
        
        // add the damping forces if rayleigh damping
        if (doRayleighDamping == 1 && (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
            theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);
    }
    
    return *theVector;
}

int
CorotTruss::sendSelf(int commitTag, Channel &theChannel)
{
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // truss packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(9);
  data(0) = this->getTag();
  data(1) = numDIM;
  data(2) = numDOF;
  data(3) = A;
  data(6) = rho;
  data(7) = doRayleighDamping;
  data(8) = cMass;
  
  data(4) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();

  // NOTE: we do have to ensure that the material has a database
  // tag if we are sending to a database channel.
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  data(5) = matDbTag;

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING Truss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }	      

  // truss then sends the tags of it's two end nodes
  res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr << "WARNING Truss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }

  // finally truss asks it's material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "WARNING Truss::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    return -3;
  }

  return 0;
}

int
CorotTruss::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res;
  int dataTag = this->getDbTag();

  // truss creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector

  static Vector data(9);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING Truss::recvSelf() - failed to receive Vector\n";
    return -1;
  }	      

  this->setTag((int)data(0));
  numDIM = (int)data(1);
  numDOF = (int)data(2);
  A = data(3);
  rho = data(6);
  doRayleighDamping = (int)data(7);
  cMass = (int)data(8);

  // truss now receives the tags of it's two external nodes
  res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr << "WARNING Truss::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  // finally truss creates a material object of the correct type,
  // sets its database tag and asks this new object to recveive itself.

  int matClass = (int)data(4);
  int matDb = (int)data(5);

  // check if we have a material object already & if we do if of right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

    // if old one .. delete it
    if (theMaterial != 0)
      delete theMaterial;

    // create a new material object
    theMaterial = theBroker.getNewUniaxialMaterial(matClass);
    if (theMaterial == 0) {
      opserr << "WARNING Truss::recvSelf() - " << this->getTag() << 
	"failed to get a blank Material of type: " << matClass << endln;
      return -3;
    }
  }

  theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "WARNING Truss::recvSelf() - " << this->getTag() << " failed to receive its Material\n";
    return -3;    
  }

  return 0;
}

int
CorotTruss::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// ensure setDomain() worked
	if (Ln == 0.0)
		return 0;

    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

void
CorotTruss::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "\nCorotTruss, tag: " << this->getTag() << endln;
		s << "\tConnected Nodes: " << connectedExternalNodes;
		s << "\tSection Area: " << A << endln;
		s << "\tUndeformed Length: " << Lo << endln;
		s << "\tCurrent Length: " << Ln << endln;
		s << "\tMass Density/Length: " << rho << endln;
		s << "\tConsistent Mass: " << cMass << endln;
		s << "\tRotation matrix: " << endln;
        
		if (theMaterial) {
			s << "\tAxial Force: " << A*theMaterial->getStress() << endln;
			s << "\tUniaxialMaterial, tag: " << theMaterial->getTag() << endln;
			theMaterial->Print(s, flag);
		}
	}
    
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"CorotTruss\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
		s << "\"A\": " << A << ", ";
		s << "\"massperlength\": " << rho << ", ";
		s << "\"material\": \"" << theMaterial->getTag() << "\"}";
	}
}

Response*
CorotTruss::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","Truss");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    //
    // we compare argv[0] for known response types for the CorotTruss
    //

    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0) 
        || (strcmp(argv[0],"globalForce") == 0) || (strcmp(argv[0],"globalForces") == 0)){
            char outputData[10];
            int numDOFperNode = numDOF/2;
            for (int i=0; i<numDOFperNode; i++) {
                sprintf(outputData,"P1_%d", i+1);
                output.tag("ResponseType", outputData);
            }
            for (int j=0; j<numDOFperNode; j++) {
                sprintf(outputData,"P2_%d", j+1);
                output.tag("ResponseType", outputData);
            }
            theResponse =  new ElementResponse(this, 1, Vector(numDOF));

    } else if ((strcmp(argv[0],"axialForce") == 0) ||
	       (strcmp(argv[0],"basicForce") == 0) || 
	       (strcmp(argv[0],"localForces") == 0) || 
	       (strcmp(argv[0],"basicForces") == 0)) {
            output.tag("ResponseType", "N");
            theResponse =  new ElementResponse(this, 2, 0.0);

    } else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformation") == 0 ||
        strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"basicDefo") == 0 ||
        strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0) {

            output.tag("ResponseType", "U");
            theResponse = new ElementResponse(this, 3, 0.0);

            // a material quantity
    }
    else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "-material") == 0) {
        output.tag("GaussPointOutput");
        output.attr("number", 1);
        output.attr("eta", 0.0);

        if (argc > 1) {
            // we need at least one more argument otherwise 
            // there is no need to forward this call to the material
            if (argc > 2) {
                // if we have 2 or more extra arguments, the first one 
                // could be an integer. In this case we check to see if it is the section id
                // (only 1 in this case)
                int sectionNum = atoi(argv[1]);
                if (sectionNum == 0) {
                    // if it is not a number we forward the call to the section as usual
                    theResponse = theMaterial->setResponse(&argv[1], argc - 1, output);
                }
                else {
                    // it is a number. Now we have to make sure it is within the allowed range
                    // for this element (in this case it can only be 1)
                    // If it is > 1, then we MUST return NULL, because the MPCO recorder iteratively
                    // uses this call to understand how many fibers we have in a section
                    if (sectionNum == 1) {
                        theResponse = theMaterial->setResponse(&argv[2], argc - 2, output);
                    }
                }
            }
            else {
                // otherwise forward it as usual
                theResponse = theMaterial->setResponse(&argv[1], argc - 1, output);
            }
        }
        output.endTag();
    }

    output.endTag();
    return theResponse;
}

int 
CorotTruss::getResponse(int responseID, Information &eleInfo)
{
    double strain;

    switch (responseID) {
    case 1:
        return eleInfo.setVector(this->getResistingForce());

    case 2:
        return eleInfo.setDouble(A * theMaterial->getStress());

    case 3:
        if (Lo == 0.0) {
            strain = 0.0;
        } else {
            strain = theMaterial->getStrain();
        }
        return eleInfo.setDouble(Lo * strain);

    default:
        return 0;
    }
}

int
CorotTruss::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;
  
  // Cross sectional area of the truss
  if (strcmp(argv[0],"A") == 0)
    return param.addObject(1, this);
  
  // Explicit specification of a material parameter
  if (strstr(argv[0],"material") != 0) {
    if (argc < 2)
      return -1;
    else
      return theMaterial->setParameter(&argv[1], argc-1, param);
  } 
  
  // Otherwise, send it to the material
  else
    return theMaterial->setParameter(argv, argc, param);
}

int
CorotTruss::updateParameter (int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    A = info.theDouble;
    return 0;
  default:
    return -1;
  }
}
