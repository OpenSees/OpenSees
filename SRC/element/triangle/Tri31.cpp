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
                                                                        
// $Revision: 1.00 $
// $Date: 2010/09/08 20:01:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/triangular/Tri31.h,v $
                                                                        
// Written: Roozbeh Geraili Mikola (roozbehg@berkeley.edu)
// Created: Sep 2010
// Revised: --------
//
// Description: This file contains the class definition for Tri31.

#include "Tri31.h"
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <ElementIter.h>
#include <elementAPI.h>
#include <map>

double Tri31::matrixData[36];
Matrix Tri31::K(matrixData, 6, 6);
Vector Tri31::P(6);
double Tri31::shp[3][3];
double Tri31::pts[1][2];
double Tri31::wts[1];

#define OPS_Export

static int num_Tri31 = 0;

OPS_Export void *
OPS_Tri31()
{
  if (num_Tri31 == 0) {
    num_Tri31++;
    opserr<<"Tri31 - Written by Roozbeh G. Mikola and N.Sitar, UC Berkeley\n";
  //OPS_Error("Tri31 - Written by Roozbeh G. Mikola and N.Sitar, UC Berkeley\n",1);
  }
  
  // Pointer to an element that will be returned
  Element *theElement = 0;
  
  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();
  
  if (numRemainingInputArgs < 7) {
    opserr << "Invalid #args, want: element element Tri31 eleTag? iNode? jNode? kNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
    return 0;
  }

  int iData[5];
  char *theType;
  double dData[5];
  dData[1] = 0.0;
  dData[2] = 0.0;
  dData[3] = 0.0;
  dData[4] = 0.0;
  
  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element Tri31\n";
    return 0;
  }
  
  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid thickness data: element Tri31 " << iData[0] << endln;
    return 0;
  }
  
  // if (OPS_GetStringCopy(&theType) != 0) {
  //   opserr << "WARNING invalid type, want: ""PlaneStress"" or ""PlaneStrain""  element SSPquad " << iData[0] << endln;
  //   return 0;
  // }
  theType = (char*)OPS_GetString();
  
  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
    opserr << "WARNING invalid integer data: element Tri31\n";
    return 0;
  }
  int matID = iData[4];
  
  NDMaterial *theMaterial = OPS_getNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element Tri31 " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }
  
  if (numRemainingInputArgs == 11) {
    numData = 4;
    if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
      opserr << "WARNING invalid optional data: element Tri31 " << iData[0] << endln;
      return 0;
    }
  }
  
  // parsing was successful, allocate the element
  theElement = new Tri31(iData[0], iData[1], iData[2], iData[3],
			 *theMaterial, theType, 
			 dData[0], dData[1], dData[2], dData[3], dData[4]);
  
  if (theElement == 0) {
    opserr << "WARNING could not create element of type Tri31\n";
    return 0;
  }
  
  return theElement;
}

void *
OPS_Tri31(const ID &info)
{
    if (num_Tri31 == 0) {
	num_Tri31++;
	opserr<<"Tri31 - Written by Roozbeh G. Mikola and N.Sitar, UC Berkeley\n";
	//OPS_Error("Tri31 - Written by Roozbeh G. Mikola and N.Sitar, UC Berkeley\n",1);
    }

    int iData[5];
    char *theType = (char*) "PlaneStress";
    double dData[5];
    dData[1] = 0.0;
    dData[2] = 0.0;
    dData[3] = 0.0;
    dData[4] = 0.0;
    int numData;

    // Pointer to an element that will be returned
    Element *theElement = 0;

    // regular element, not in a mesh, get tags
    if (info.Size() == 0) {
	int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

	if (numRemainingInputArgs < 4) {
	    opserr << "Invalid #args, want: element element Tri31 eleTag? iNode? jNode? kNode?\n";
	    return 0;
	}

	numData = 4;
	if (OPS_GetIntInput(&numData, iData) != 0) {
	    opserr << "WARNING invalid integer data: element Tri31\n";
	    return 0;
	}
    }

    // regular element, or in a mesh
    if (info.Size()==0 || info(0)==1) {
	if(OPS_GetNumRemainingInputArgs() < 3) {
	    opserr<<"insufficient arguments: thk? type? matTag? <pressure? rho? b1? b2?>\n";
	    return 0;
	}

	numData = 1;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
	    opserr << "WARNING invalid thickness data: element Tri31 " << endln;
	    return 0;
	}

	// if (OPS_GetStringCopy(&theType) != 0) {
	//   opserr << "WARNING invalid type, want: ""PlaneStress"" or ""PlaneStrain""  element SSPquad " << iData[0] << endln;
	//   return 0;
	// }
	theType = (char*)OPS_GetString();

	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
	    opserr << "WARNING invalid integer data: element Tri31\n";
	    return 0;
	}

	if (OPS_GetNumRemainingInputArgs() == 4) {
	    numData = 4;
	    if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
		opserr << "WARNING invalid optional data: element Tri31 " << endln;
		return 0;
	    }
	}
    }

    // store data for different mesh
    static std::map<int, Vector> meshdata;
    if (info.Size()>0 && info(0)==1) {
	if (info.Size() < 2) {
	    opserr << "WARNING: need info -- inmesh, meshtag\n";
	    return 0;
	}

	// save the data for a mesh
	Vector& mdata = meshdata[info(1)];
	mdata.resize(7);
	for (int i=0; i<5; ++i) {
	    mdata(i) = dData[i];
	}
	mdata(5) = iData[4];
	if (strcmp(theType,"PlaneStrain") == 0 ||
	    strcmp(theType,"PlaneStrain2D") == 0) {
	    mdata(6) = 1;
	} else if (strcmp(theType,"PlaneStress") == 0 ||
		   strcmp(theType,"PlaneStress2D") == 0) {
	    mdata(6) = 2;
	}

	return &meshdata;

    } else if (info.Size()>0 && info(0)==2) {
	if (info.Size() < 6) {
	    opserr << "WARNING: need info -- inmesh, meshtag, eleTag, nd1, nd2, nd3\n";
	    return 0;
	}

	// get the data for a mesh
	Vector& mdata = meshdata[info(1)];
	if (mdata.Size() < 7) return 0;

	for (int i=0; i<5; ++i) {
	    dData[i] = mdata(i);
	}
	for (int i=0; i<4; ++i) {
	    iData[i] = info(2+i);
	}
	iData[4] = mdata(5);
	if (mdata(6) == 1) {
	    theType = (char*)"PlaneStrain";
	} else if (mdata(6) == 2) {
	    theType = (char*)"PlaneStress";
	}
    }

    int matID = iData[4];
    NDMaterial *theMaterial = OPS_getNDMaterial(matID);
    if (theMaterial == 0) {
	opserr << "WARNING element Tri31 " << iData[0] << endln;
	opserr << " Material: " << matID << "not found\n";
	return 0;
    }

    // parsing was successful, allocate the element
    theElement = new Tri31(iData[0], iData[1], iData[2], iData[3],
			   *theMaterial, theType,
			   dData[0], dData[1], dData[2], dData[3], dData[4]);

    if (theElement == 0) {
	opserr << "WARNING could not create element of type Tri31\n";
	return 0;
    }

    return theElement;
}

int OPS_Tri31(Domain& theDomain, const ID& elenodes, ID& eletags)
{
    // get inputs
    int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();
    if (numRemainingInputArgs < 3) {
	opserr << "Invalid #args, want: thk? type? matTag? <pressure? rho? b1? b2?>\n";
	return -1;
    }

    int matID;
    double thk;
    char *theType;
    double dData[4];
    dData[0] = 0.0;
    dData[1] = 0.0;
    dData[2] = 0.0;
    dData[3] = 0.0;
 
    int numData = 1;
    if (OPS_GetDoubleInput(&numData, &thk) != 0) {
	opserr << "WARNING invalid thickness data: element Tri31 \n";
	return -1;
    }

    theType = (char*)OPS_GetString();
  
    numData = 1;
    if (OPS_GetIntInput(&numData, &matID) != 0) {
	opserr << "WARNING invalid integer data: element Tri31\n";
	return -1;
    }
  
    NDMaterial *theMaterial = OPS_getNDMaterial(matID);
    if (theMaterial == 0) {
	opserr << "WARNING element Tri31 \n";
	opserr << " Material: " << matID << "not found\n";
	return -1;
    }
  
    if (OPS_GetNumRemainingInputArgs() >= 4) {
	numData = 4;
	if (OPS_GetDoubleInput(&numData, &dData[0]) != 0) {
	    opserr << "WARNING invalid optional data: element Tri31\n";
	    return -1;
	}
    }

    // create elements
    ElementIter& theEles = theDomain.getElements();
    Element* theEle = theEles();
    int currTag = 0;
    if (theEle != 0) {
	currTag = theEle->getTag();
    }

    eletags.resize(elenodes.Size()/3);

    for (int i=0; i<eletags.Size(); i++) {
	theEle = new Tri31(--currTag,elenodes(3*i),elenodes(3*i+1),elenodes(3*i+2),
			   *theMaterial, theType, thk,
			   dData[0], dData[1], dData[2], dData[3]);
	if (theEle == 0) {
	    opserr<<"WARNING: run out of memory for creating element\n";
	    return -1;
	}
	if (theDomain.addElement(theEle) == false) {
	    opserr<<"WARNING: failed to add element to domain\n";
	    delete theEle;
	    return -1;
	}
	eletags(i) = currTag;
    }
  
    return 0;
}

Tri31::Tri31(int tag, int nd1, int nd2, int nd3,
	     NDMaterial &m, const char *type, double t,
	     double p, double r, double b1, double b2)
:Element (tag, ELE_TAG_Tri31), 
  theMaterial(0), connectedExternalNodes(3), 
 Q(6), pressureLoad(6), thickness(t), pressure(p), rho(r), Ki(0)
{
	pts[0][0] = 0.333333333333333;
	pts[0][1] = 0.333333333333333;

	wts[0] = 0.5;

	if (strcmp(type,"PlaneStrain") != 0 && strcmp(type,"PlaneStress") != 0
	    && strcmp(type,"PlaneStrain2D") != 0 && strcmp(type,"PlaneStress2D") != 0) {
			opserr << "Tri31::Tri31 -- improper material type: " << type << "for Tri31\n";
	        exit(-1);
	}

	// Body forces
	b[0] = b1;
	b[1] = b2;

    // Allocate arrays of pointers to NDMaterials
    theMaterial = new NDMaterial *[numgp];
    
    if (theMaterial == 0) {
      opserr << "Tri31::Tri31 - failed allocate material model pointer\n";
      exit(-1);
    }

	int i;
    for (i = 0; i < numgp; i++) {
		// Get copies of the material model for each integration point
        theMaterial[i] = m.getCopy(type);
			
        // Check allocation
        if (theMaterial[i] == 0) {
	   	    opserr << "Tri31::Tri31 -- failed to get a copy of material model\n";
	        exit(-1);
        }
    }

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    
    for (i=0; i<numnodes; i++) theNodes[i] = 0;
}

Tri31::Tri31()
:Element (0,ELE_TAG_Tri31),
  theMaterial(0), connectedExternalNodes(3), Q(6), pressureLoad(6), thickness(0.0), pressure(0.0), Ki(0)
{
	pts[0][0] = 0.333333333333333;
	pts[0][1] = 0.333333333333333;

	wts[0] = 0.5;

    for (int i=0; i<numnodes; i++) theNodes[i] = 0;
}

Tri31::~Tri31()
{    
  for (int i = 0; i < numgp; i++) {
	  if (theMaterial[i]) delete theMaterial[i];
  }

  // Delete the array of pointers to NDMaterial pointer arrays
  if (theMaterial) delete [] theMaterial;

  if (Ki != 0) delete Ki;
}

int
Tri31::getNumExternalNodes() const
{
    return numnodes;
}

const ID&
Tri31::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **
Tri31::getNodePtrs(void) 
{
  return theNodes;
}

int
Tri31::getNumDOF()
{
  int sum = 0;

  sum += theNodes[0]->getNumberDOF();
  sum += theNodes[1]->getNumberDOF();
  sum += theNodes[2]->getNumberDOF();

  return sum;
}

void
Tri31::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
		theNodes[0] = 0;
	    theNodes[1] = 0;
	    theNodes[2] = 0;
	    return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int Nd3 = connectedExternalNodes(2);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    theNodes[2] = theDomain->getNode(Nd3);

    if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0) {
		//opserr << "FATAL ERROR Tri31 (tag: %d), node not found in domain",
	    //this->getTag());
	    return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    
    if (dofNd1 != 2 && dofNd1 != 3) {
      opserr << "WARNING Tri31::setDomain() element " << this->getTag() << " does not have 2 or 3 DOF at node " << theNodes[0]->getTag() << endln;
      return;
    }
    if (dofNd2 != 2 && dofNd2 != 3) {
      opserr << "WARNING Tri31::setDomain() element " << this->getTag() << " does not have 2 or 3 DOF at node " << theNodes[1]->getTag() << endln;
      return;
    }
    if (dofNd3 != 2 && dofNd3 != 3) {
      opserr << "WARNING Tri31::setDomain() element " << this->getTag() << " does not have 2 or 3 DOF at node " << theNodes[2]->getTag() << endln;
      return;
    }

    this->DomainComponent::setDomain(theDomain);

    // Compute consistent nodal loads due to pressure
    this->setPressureLoadAtNodes();
}

int
Tri31::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
		opserr << "Tri31::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i < numgp; i++) retVal += theMaterial[i]->commitState();

    return retVal;
}

int
Tri31::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numgp; i++) retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}

int
Tri31::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numgp; i++) retVal += theMaterial[i]->revertToStart();

    return retVal;
}


int
Tri31::update()
{
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	
	static double u[2][3];

	u[0][0] = disp1(0);
	u[1][0] = disp1(1);
	u[0][1] = disp2(0);
	u[1][1] = disp2(1);
	u[0][2] = disp3(0);
	u[1][2] = disp3(1);

	static Vector eps(3);

	int ret = 0;

	// Loop over the integration points
	for (int i = 0; i < numgp; i++) {

		// Determine Jacobian for this integration point
		this->shapeFunction(pts[i][0], pts[i][1]);

		// Interpolate strains
		//eps = B*u;
		//eps.addMatrixVector(0.0, B, u, 1.0);
		eps.Zero();
		for (int beta = 0; beta < numnodes; beta++) {
			eps(0) += shp[0][beta]*u[0][beta];
			eps(1) += shp[1][beta]*u[1][beta];
			eps(2) += shp[0][beta]*u[1][beta] + shp[1][beta]*u[0][beta];
		}

		// Set the material strain
		ret += theMaterial[i]->setTrialStrain(eps);
	}

	return ret;
}

const Matrix&
Tri31::getTangentStiff()                                                                    
{

	K.Zero();

	double dvol;
	double DB[3][2];

	// Loop over the integration points
	for (int i = 0; i < numgp; i++) {
		// Determine Jacobian for this integration point
	    dvol = this->shapeFunction(pts[i][0], pts[i][1]);
	    dvol *= (thickness*wts[i]);
	  
	    // Get the material tangent
	    const Matrix &D = theMaterial[i]->getTangent();
	  
	    // Perform numerical integration
	    //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
	    //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);
	  
	    double D00 = D(0,0); double D01 = D(0,1); double D02 = D(0,2);
	    double D10 = D(1,0); double D11 = D(1,1); double D12 = D(1,2);
	    double D20 = D(2,0); double D21 = D(2,1); double D22 = D(2,2);

	    //for (int beta = 0, ib = 0, colIb =0, colIbP1 = 8; 
	    //   beta < 4; 
	    //   beta++, ib += 2, colIb += 16, colIbP1 += 16) {
	    
	    for (int alpha = 0, ia = 0; alpha < numnodes; alpha++, ia += 2) {
			for (int beta = 0, ib = 0; beta < numnodes; beta++, ib += 2) {

				DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
	            DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
	            DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
	            DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
	            DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
	            DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);
	      
	            K(ia,ib)     += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	            K(ia,ib+1)   += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	            K(ia+1,ib)   += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	            K(ia+1,ib+1) += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
	            //matrixData[colIb   +   ia] += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	            //matrixData[colIbP1 +   ia] += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	            //matrixData[colIb   + ia+1] += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	            //matrixData[colIbP1 + ia+1] += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
	      
	        }
	   }
	}
	
	return K;
}


const Matrix& 
Tri31::getInitialStiff()                                                               
{
	if (Ki != 0) return *Ki;

    K.Zero();
  
    double dvol;
    double DB[3][2];
  
    // Loop over the integration points
    for (int i = 0; i < numgp; i++) {

  	    // Determine Jacobian for this integration point
        dvol = this->shapeFunction(pts[i][0], pts[i][1]);
        dvol *= (thickness*wts[i]);
    
        // Get the materialmgp tangent
        const Matrix &D = theMaterial[i]->getInitialTangent();

        double D00 = D(0,0); double D01 = D(0,1); double D02 = D(0,2);
        double D10 = D(1,0); double D11 = D(1,1); double D12 = D(1,2);
        double D20 = D(2,0); double D21 = D(2,1); double D22 = D(2,2);
    
        // Perform numerical integration
        //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
        //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);
        //for (int beta = 0, ib = 0, colIb =0, colIbP1 = 2*numnodes; beta < numnodes; beta++, ib += 2, colIb += numnodes*numnodes, colIbP1 += numnodes*numnodes) {
	    //	  for (int alpha = 0, ia = 0; alpha < numnodes; alpha++, ia += 2) {

	    for (int alpha = 0, ia = 0; alpha < numnodes; alpha++, ia += 2) {
		    for (int beta = 0, ib = 0; beta < numnodes; beta++, ib += 2) {

		  	    DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
	            DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
	            DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
	            DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
	            DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
	            DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);

	            K(ia,ib)     += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	            K(ia,ib+1)   += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	            K(ia+1,ib)   += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	            K(ia+1,ib+1) += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
	
	            //matrixData[colIb   +   ia] += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	            //matrixData[colIbP1 +   ia] += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	            //matrixData[colIb   + ia+1] += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	            //matrixData[colIbP1 + ia+1] += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
           }
        }
    }

    Ki = new Matrix(K);
    return K;
}

const Matrix&
Tri31::getMass()
{
	K.Zero();

	int i;
	static double rhoi[1]; //numgp
	double sum = 0.0;
	for (i = 0; i < numgp; i++) {
		if (rho == 0)
			rhoi[i] = theMaterial[i]->getRho();
	    else
			rhoi[i] = rho;	    
	    sum += rhoi[i];
	}

	if (sum == 0.0) return K;

	double rhodvol, Nrho;

	// Compute a lumped mass matrix
	for (i = 0; i < numgp; i++) {
		
		// Determine Jacobian for this integration point
		rhodvol = this->shapeFunction(pts[i][0], pts[i][1]);

		// Element plus material density ... MAY WANT TO REMOVE ELEMENT DENSITY
		rhodvol *= (rhoi[i]*thickness*wts[i]);

		for (int alpha = 0, ia = 0; alpha < numnodes; alpha++, ia++) {
			Nrho = shp[2][alpha]*rhodvol;
			K(ia,ia) += Nrho;
			ia++;
			K(ia,ia) += Nrho;
		}
	}

	return K;
}

void
Tri31::zeroLoad(void)
{
    applyLoad = 0;
    appliedB[0] = 0.0;
    appliedB[1] = 0.0;

	Q.Zero();
    return;
}

int 
Tri31::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	// body forces can be applied in a load pattern
	int type;
	const Vector &data = theLoad->getData(type, loadFactor);

	if (type == LOAD_TAG_SelfWeight) {
		applyLoad = 1;
		appliedB[0] += loadFactor*data(0)*b[0];
		appliedB[1] += loadFactor*data(1)*b[1];
		return 0;
	} else {
	    opserr << "Tri31::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
		return -1;
	}

    return -1;
}

int 
Tri31::addInertiaLoadToUnbalance(const Vector &accel)
{
	int i;
	static double rhoi[1]; //numgp
	double sum = 0.0;
	for (i = 0; i < numgp; i++) {
            if(rho == 0) {
                rhoi[i] = theMaterial[i]->getRho();
            } else {
                rhoi[i] = rho;
            }
	    sum += rhoi[i];
	}

	if (sum == 0.0)
	return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = theNodes[0]->getRV(accel);
	const Vector &Raccel2 = theNodes[1]->getRV(accel);
	const Vector &Raccel3 = theNodes[2]->getRV(accel);

	if (2 != Raccel1.Size() || 2 != Raccel2.Size() || 2 != Raccel3.Size()) {
		opserr << "Tri31::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
	    return -1;
	}

	static double ra[6];

	ra[0] = Raccel1(0);
	ra[1] = Raccel1(1);
	ra[2] = Raccel2(0);
	ra[3] = Raccel2(1);
	ra[4] = Raccel3(0);
	ra[5] = Raccel3(1);

	// Compute mass matrix
	this->getMass();

	// Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	for (i = 0; i < 2*numnodes; i++) Q(i) += -K(i,i)*ra[i];

	return 0;
}

const Vector&
Tri31::getResistingForce()
{
	P.Zero();

	double dvol;

	// Loop over the integration points
	for (int i = 0; i < numgp; i++) {

		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i][0], pts[i][1]);
		dvol *= (thickness*wts[i]);

		// Get material stress response
		const Vector &sigma = theMaterial[i]->getStress();

		// Perform numerical integration on internal force
		//P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
		//P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);
		for (int alpha = 0, ia = 0; alpha < numnodes; alpha++, ia += 2) {                                    
			
			P(ia) += dvol*(shp[0][alpha]*sigma(0) + shp[1][alpha]*sigma(2));
			
			P(ia+1) += dvol*(shp[1][alpha]*sigma(1) + shp[0][alpha]*sigma(2));

			// Subtract equiv. body forces from the nodes
			//P = P - (N^ b) * intWt(i)*intWt(j) * detJ;
			//P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);
			if (applyLoad == 0) {
			  P(ia)   -= dvol*(shp[2][alpha]*b[0]);
			  P(ia+1) -= dvol*(shp[2][alpha]*b[1]);
			} else {
			  P(ia)   -= dvol*(shp[2][alpha]*appliedB[0]);
			  P(ia+1) -= dvol*(shp[2][alpha]*appliedB[1]);
			}
		}
	}

	// Subtract pressure loading from resisting force
	if (pressure != 0.0) {
		//P = P - pressureLoad;
		P.addVector(1.0, pressureLoad, -1.0);
	}
	
	// Subtract other external nodal loads ... P_res = P_int - P_ext
	//P = P - Q;
	P.addVector(1.0, Q, -1.0);

	return P;
}

const Vector&
Tri31::getResistingForceIncInertia()
{
	int i;
	static double rhoi[1]; //numgp
	double sum = 0.0;
	for (i = 0; i < numgp; i++) {
            if(rho == 0) {
		rhoi[i] = theMaterial[i]->getRho();
            } else {
                rhoi[i] = rho;
            }
	    sum += rhoi[i];
	}

	// if no mass terms .. just add damping terms
	if (sum == 0.0) {
		this->getResistingForce();

	    // add the damping forces if rayleigh damping
	    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) P += this->getRayleighDampingForces();

	    return P;
	}

	const Vector &accel1 = theNodes[0]->getTrialAccel();
	const Vector &accel2 = theNodes[1]->getTrialAccel();
	const Vector &accel3 = theNodes[2]->getTrialAccel();
	
	static double a[6];

	a[0] = accel1(0);
	a[1] = accel1(1);
	a[2] = accel2(0);
	a[3] = accel2(1);
	a[4] = accel3(0);
	a[5] = accel3(1);

	// Compute the current resisting force
	this->getResistingForce();

	// Compute the mass matrix
	this->getMass();

	// Take advantage of lumped mass matrix
	for (i = 0; i < 2*numnodes; i++) P(i) += K(i,i)*a[i];

	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) P += this->getRayleighDampingForces();

	return P;
}

int
Tri31::sendSelf(int commitTag, Channel &theChannel)                                                           
{
	int res = 0;
  
    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();
  
    // Tri31 packs its data into a Vector and sends this to theChannel
    // along with its dbTag and the commitTag passed in the arguments
    static Vector data(10);
    data(0) = this->getTag();
    data(1) = thickness;
    data(3) = b[0];
    data(4) = b[1];
    data(5) = pressure;

    data(6) = alphaM;
    data(7) = betaK;
    data(8) = betaK0;
    data(9) = betaKc;
  
    res += theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
		opserr << "WARNING Tri31::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }	      
  
    // Now Tri31 sends the ids of its materials
    int matDbTag;
	int count=0;
  
    static ID idData(2*numgp+numnodes+1);
  
    int i;
    for (i = 0; i < numgp; i++) { 
		idData(i) = theMaterial[i]->getClassTag();
        matDbTag = theMaterial[i]->getDbTag();
        // NOTE: we do have to ensure that the material has a database
        // tag if we are sending to a database channel.
        if (matDbTag == 0) {
			matDbTag = theChannel.getDbTag();
			if (matDbTag != 0) theMaterial[i]->setDbTag(matDbTag);
        }
        idData(i+numgp) = matDbTag;
    }
    count += 2*numgp;
    idData(count) = connectedExternalNodes(0); count += 1;
    idData(count) = connectedExternalNodes(1); count += 1;
    idData(count) = connectedExternalNodes(2);

    res += theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
		opserr << "WARNING Tri31::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // Finally, Tri31 asks its material objects to send themselves
    for (i = 0; i < numgp; i++) {
		res += theMaterial[i]->sendSelf(commitTag, theChannel);
        if (res < 0) {
			opserr << "WARNING Tri31::sendSelf() - " << this->getTag() << " failed to send its Material\n";
            return res;
        }
    }
  
    return res;
}

int
Tri31::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)                                 
{
	int res = 0;

	int dataTag = this->getDbTag();

	// Tri31 creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	static Vector data(10);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING Tri31::recvSelf() - failed to receive Vector\n";
	    return res;
	}

	this->setTag((int)data(0));
	thickness = data(1);
	b[0] = data(3);
	b[1] = data(4);
	pressure = data(5);

	alphaM = data(6);
	betaK  = data(7);
	betaK0 = data(8);
	betaKc = data(9);

	static ID idData(2*numgp+numnodes+1);
	// Tri31 now receives the tags of its four external nodes
	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING Tri31::recvSelf() - " << this->getTag() << " failed to receive ID\n";
	    return res;
	}

	int count = 2*numgp;

	connectedExternalNodes(0) = idData(count); count += 1;
	connectedExternalNodes(1) = idData(count); count += 1;
	connectedExternalNodes(2) = idData(count);

	if (theMaterial == 0) {
		// Allocate new materials
	    theMaterial = new NDMaterial *[numgp];
		if (theMaterial == 0) {
			opserr << "Tri31::recvSelf() - Could not allocate NDMaterial* array\n";
			return -1;
		}
		for (int i = 0; i < numgp; i++) {
			int matClassTag = idData(i);
			int matDbTag = idData(i+numgp);
			// Allocate new material with the sent class tag
			theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
			if (theMaterial[i] == 0) {
				opserr << "Tri31::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
				return -1;
			}
			// Now receive materials into the newly allocated space
			theMaterial[i]->setDbTag(matDbTag);
			res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "Tri31::recvSelf() - material " << i << " failed to recv itself\n";
				return res;
			}
		}
	}

	// materials exist , ensure materials of correct type and recvSelf on them
	else {
		for (int i = 0; i < numgp; i++) {
			int matClassTag = idData(i);
		    int matDbTag = idData(i+numgp);
		    // Check that material is of the right type; if not,
		    // delete it and create a new one of the right type
		    if (theMaterial[i]->getClassTag() != matClassTag) {
				delete theMaterial[i];
		        theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
		        if (theMaterial[i] == 0) {
					opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to create\n";
					
		            return -1;
		        }
		    }
		    // Receive the material
		    theMaterial[i]->setDbTag(matDbTag);
		    res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
		    if (res < 0) {
				opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to recv itself\n";
		        return res;
		    }
		}
	}

	return res;
}

void
Tri31::Print(OPS_Stream &s, int flag)                                                                 
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nTri31, element id:  " << this->getTag() << endln;
        s << "\tConnected external nodes:  " << connectedExternalNodes;
        s << "\tthickness:  " << thickness << endln;
        s << "\tsurface pressure:  " << pressure << endln;
        s << "\tmass density:  " << rho << endln;
        s << "\tbody forces:  " << b[0] << " " << b[1] << endln;
        theMaterial[0]->Print(s, flag);
        s << "\tStress (xx yy xy)" << endln;
        for (int i = 0; i<numgp; i++) s << "\t\tGauss point " << i + 1 << ": " << theMaterial[i]->getStress();
    }

    if (flag == 2) {

        s << "#Tri31\n";

        int i;
        const int numNodes = numnodes;
        const int nstress = numgp;

        for (i = 0; i < numNodes; i++) {
            const Vector &nodeCrd = theNodes[i]->getCrds();
            // const Vector &nodeDisp = theNodes[i]->getDisp();
            s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
        }

        // spit out the section location & invoke print on the scetion
        const int numMaterials = numgp;

        static Vector avgStress(nstress);
        static Vector avgStrain(nstress);
        avgStress.Zero();
        avgStrain.Zero();
        for (i = 0; i < numMaterials; i++) {
            avgStress += theMaterial[i]->getStress();
            avgStrain += theMaterial[i]->getStrain();
        }
        avgStress /= numMaterials;
        avgStrain /= numMaterials;

        s << "#AVERAGE_STRESS ";
        for (i = 0; i < nstress; i++) s << avgStress(i) << " ";
        s << endln;

        s << "#AVERAGE_STRAIN ";
        for (i = 0; i < nstress; i++) s << avgStrain(i) << " ";
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"Tri31\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
        s << connectedExternalNodes(1) << ", ";
        s << connectedExternalNodes(2) << "], ";
        s << "\"thickness\": " << thickness << ", ";
        s << "\"surfacePressure\": " << pressure << ", ";
        s << "\"masspervolume\": " << rho << ", ";
        s << "\"bodyForces\": [" << b[0] << ", " << b[1] << "], ";
        s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
    }
}

int
Tri31::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// get the end point display coords
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	theNodes[0]->getDisplayCrds(v1, fact, displayMode);
	theNodes[1]->getDisplayCrds(v2, fact, displayMode);
	theNodes[2]->getDisplayCrds(v3, fact, displayMode);

	// place values in coords matrix
	static Matrix coords(3, 3);
	for (int i = 0; i < 3; i++) {
		coords(0, i) = v1(i);
		coords(1, i) = v2(i);
		coords(2, i) = v3(i);
	}

	// set the quantity to be displayed at the nodes;
	// if displayMode is 1 through 3 we will plot material stresses otherwise 0.0
	static Vector values(numgp);
	if (displayMode < numgp && displayMode > 0) {
		for (int i = 0; i < numgp; i++) {
			const Vector& stress = theMaterial[i]->getStress();
			values(i) = stress(displayMode - 1);
		}
	}
	else {
		for (int i = 0; i < numgp; i++)
			values(i) = 0.0;
	}

	// draw the polygon
	return theViewer.drawPolygon(coords, values, this->getTag());
}

Response*
Tri31::setResponse(const char **argv, int argc, OPS_Stream &output)                                         
{
	Response *theResponse =0;

    output.tag("ElementOutput");
    output.attr("eleType","Tri31");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);
    output.attr("node3",connectedExternalNodes[2]);

    char dataOut[10];
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
		for (int i=1; i<=numgp; i++) {
			sprintf(dataOut,"P1_%d",i);
            output.tag("ResponseType",dataOut);
            sprintf(dataOut,"P2_%d",i);
            output.tag("ResponseType",dataOut);
        }
    
        theResponse =  new ElementResponse(this, 1, P);
    } else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {

		int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= numgp) {
			output.tag("GaussPoint");
            output.attr("number",pointNum);
            output.attr("eta",pts[pointNum-1][0]);
            output.attr("neta",pts[pointNum-1][1]);

            theResponse =  theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
            output.endTag();

       }
   } else if ((strcmp(argv[0],"stresses") ==0) || (strcmp(argv[0],"stress") ==0)) {
	   for (int i=0; i<numgp; i++) {
		   output.tag("GaussPoint");
           output.attr("number",i+1);
           output.attr("eta",pts[i][0]);
           output.attr("neta",pts[i][1]);

           output.tag("NdMaterialOutput");
           output.attr("classType", theMaterial[i]->getClassTag());
           output.attr("tag", theMaterial[i]->getTag());
      
           output.tag("ResponseType","sigma11");
           output.tag("ResponseType","sigma22");
           output.tag("ResponseType","sigma12");
      
           output.endTag(); // GaussPoint
           output.endTag(); // NdMaterialOutput
       }
       theResponse =  new ElementResponse(this, 3, Vector(3*numgp));
   }

  else if ((strcmp(argv[0],"stressesAtNodes") ==0) || (strcmp(argv[0],"stressAtNodes") ==0)) {
    for (int i=0; i<numnodes; i++) {
      output.tag("NodalPoint");
      output.attr("number",i+1);
      // output.attr("eta",pts[i][0]);
      // output.attr("neta",pts[i][1]);

      // output.tag("NdMaterialOutput");
      // output.attr("classType", theMaterial[i]->getClassTag());
      // output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");

      output.endTag(); // GaussPoint
      // output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 11, Vector(3*numnodes));
  }


   output.endTag(); // ElementOutput

   return theResponse;
}

int 
Tri31::getResponse(int responseID, Information &eleInfo)                                         
{

	if (responseID == 1) {
		return eleInfo.setVector(this->getResistingForce());
    } else if (responseID == 3) {
		// Loop over the integration points
        static Vector stresses(3*numgp);
        int cnt = 0;
        for (int i = 0; i < numgp; i++) {
			// Get material stress response
            const Vector &sigma = theMaterial[i]->getStress();
            stresses(cnt) = sigma(0);
            stresses(cnt+1) = sigma(1);
            stresses(cnt+2) = sigma(2);
            cnt += 3;
       }
       return eleInfo.setVector(stresses);

  } else if (responseID == 11) {

    // extrapolate stress from Gauss points to element nodes
    static Vector stressGP(3*numgp);
    static Vector stressAtNodes(3*numnodes); // 3*nnodes
	stressAtNodes.Zero();
    int cnt = 0;
	// first get stress components (xx, yy, xy) at Gauss points
    for (int i = 0; i < numgp; i++) {
      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      stressGP(cnt) = sigma(0);
      stressGP(cnt+1) = sigma(1);
      stressGP(cnt+2) = sigma(2);
      cnt += 3;
    }

	double We[numnodes][numgp] = {{1.0},
								  {1.0},
								  {1.0}};

	int p, l;
	for (int i = 0; i < numnodes; i++) {
	  for (int k = 0; k < 3; k++) { // number of stress components
		p = 3*i + k;
		for (int j = 0; j < numgp; j++) {
		  l = 3*j + k;
		  stressAtNodes(p) += We[i][j] * stressGP(l);
		  // opserr << "stressAtNodes(" << p << ") = We[" << i << "][" << j << "] * stressGP(" << l << ") = " << We[i][j] << " * " << stressGP(l) << " = " << stressAtNodes(p) <<  "\n";
		}
	  }
	}

    return eleInfo.setVector(stressAtNodes);

  } else

    return -1;
}

int
Tri31::setParameter(const char **argv, int argc, Parameter &param)                                     
{
	if (argc < 1) return -1;

    int res = -1;

    // Tri31 pressure loading
    if (strcmp(argv[0],"pressure") == 0) return param.addObject(2, this);

    // a material parameter
    else if (strstr(argv[0],"material") != 0) {
		if (argc < 3) return -1;

        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= numgp) return theMaterial[pointNum-1]->setParameter(&argv[2], argc-2, param);
		else return -1;
    }

    // otherwise it could be just a forall material parameter
    else {
		int matRes = res;
        for (int i=0; i<numgp; i++) {
			matRes =  theMaterial[i]->setParameter(argv, argc, param);
            if (matRes != -1) res = matRes;
        }
    }
  
    return res;
}
    
int
Tri31::updateParameter(int parameterID, Information &info)                                          
{
	switch (parameterID) {
		case -1:
			return -1;
	  
	    case 2:
			pressure = info.theDouble;
		    this->setPressureLoadAtNodes();	// update consistent nodal loads
		    return 0;
	    default: 
			/*	  
		    if (parameterID >= 100) { // material parameter
			   int pointNum = parameterID/100;
			   if (pointNum > 0 && pointNum <= 4)
			      return theMaterial[pointNum-1]->updateParameter(parameterID-100*pointNum, info);
			   else
			      return -1;
		    } else // unknown
		    */
		    return -1;
	}
}

double Tri31::shapeFunction(double xi, double eta)
{
	const Vector &nd1Crds = theNodes[0]->getCrds();
	const Vector &nd2Crds = theNodes[1]->getCrds();
	const Vector &nd3Crds = theNodes[2]->getCrds();

	shp[2][0] = xi;	        // N_1
	shp[2][1] = eta;		// N_2
	shp[2][2] = 1-xi-eta;	// N_3

	double J[2][2];

	// See p 180 "A First Course in Finite Elements" by Fish and Belytschko.
	J[0][0] = (nd1Crds(0) - nd3Crds(0));
	J[0][1] = (nd2Crds(0) - nd3Crds(0));
	J[1][0] = (nd1Crds(1) - nd3Crds(1));
	J[1][1] = (nd2Crds(1) - nd3Crds(1));

	double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
	double oneOverdetJ = 1.0/detJ;
	double L[2][2];

	// L = inv(J)
	L[0][0] =  J[1][1]*oneOverdetJ;
	L[1][0] = -J[0][1]*oneOverdetJ;
	L[0][1] = -J[1][0]*oneOverdetJ;
	L[1][1] =  J[0][0]*oneOverdetJ;

	// See Cook, Malkus, Plesha p. 169 for the derivation of these terms
    shp[0][0] = L[0][0];	            // N_1,1
    shp[0][1] = L[0][1];		        // N_2,1
    shp[0][2] = -(L[0][0] + L[0][1]);	// N_3,1
	
    shp[1][0] = L[1][0];	            // N_1,2
    shp[1][1] = L[1][1];		        // N_2,2
    shp[1][2] = -(L[1][0] + L[1][1]);	// N_3,2

    return detJ;
}

void 
Tri31::setPressureLoadAtNodes(void)
{
    pressureLoad.Zero();

	if (pressure == 0.0) return;

	const Vector &node1 = theNodes[0]->getCrds();
	const Vector &node2 = theNodes[1]->getCrds();
	const Vector &node3 = theNodes[2]->getCrds();

	double x1 = node1(0);
	double y1 = node1(1);
	double x2 = node2(0);
	double y2 = node2(1);
	double x3 = node3(0);
	double y3 = node3(1);

	double dx12 = x2-x1;
	double dy12 = y2-y1;
	double dx23 = x3-x2;
	double dy23 = y3-y2;
	double dx31 = x1-x3;
	double dy31 = y1-y3;

	double pressureOver2 = pressure/2.0;

	// Contribution from side 12
	pressureLoad(0) += pressureOver2*dy12;
	pressureLoad(2) += pressureOver2*dy12;
	pressureLoad(1) += pressureOver2*-dx12;
	pressureLoad(3) += pressureOver2*-dx12;

	// Contribution from side 23
	pressureLoad(2) += pressureOver2*dy23;
	pressureLoad(4) += pressureOver2*dy23;
	pressureLoad(3) += pressureOver2*-dx23;
	pressureLoad(5) += pressureOver2*-dx23;

	// Contribution from side 31
	pressureLoad(4) += pressureOver2*dy31;
	pressureLoad(0) += pressureOver2*dy31;
	pressureLoad(5) += pressureOver2*-dx31;
	pressureLoad(1) += pressureOver2*-dx31;

	//pressureLoad = pressureLoad*thickness;
}








