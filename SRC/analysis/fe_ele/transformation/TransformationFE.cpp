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
                                                                        
// $Revision: 1.20 $
// $Date: 2008-11-19 23:39:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/transformation/TransformationFE.cpp,v $
                                                                        
// Written: fmk 
// Created: 05/99
//
// Purpose: This file contains the code for implementing the methods
// of the TransformationFE class interface.

#include "TransformationFE.h"
#include <stdlib.h>

#include <Element.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <Integrator.h>
#include <Subdomain.h>
#include <AnalysisModel.h>
#include <Matrix.h>
#include <Vector.h>
#include <TransformationConstraintHandler.h>

#define MAX_NUM_DOF 64

// static variables initialisation
Matrix **TransformationFE::modMatrices; 
Vector **TransformationFE::modVectors;  
Matrix **TransformationFE::theTransformations; 
int TransformationFE::numTransFE(0);           
int TransformationFE::transCounter(0);           
int TransformationFE::sizeTransformations(0);          
double *TransformationFE::dataBuffer = 0;          
double *TransformationFE::localKbuffer = 0;          
int    *TransformationFE::dofData = 0;    ;          
int TransformationFE::sizeBuffer(0);            

//  TransformationFE(Element *, Integrator *theIntegrator);
//	construictor that take the corresponding model element.
TransformationFE::TransformationFE(int tag, Element *ele)
:FE_Element(tag, ele), theDOFs(0), numSPs(0), theSPs(0), modID(0), 
  modTangent(0), modResidual(0), numGroups(0), numTransformedDOF(0)
{
  // set number of original dof at ele
    numOriginalDOF = ele->getNumDOF();

    // create the array of pointers to DOF_Groups
    const ID &nodes = ele->getExternalNodes();
    Domain *theDomain = ele->getDomain();
    int numNodes = nodes.Size();
    theDOFs = new DOF_Group *[numNodes];
    if (theDOFs == 0) {
	opserr << "FATAL TransformationFE::TransformationFE() - out of memory craeting ";
	opserr << "array of size : " << numNodes << " for storage of DOF_Group\n";
	exit(-1);
    }

    numGroups = numNodes;

    // now fill the array of DOF_Group pointers
    for (int i=0; i<numNodes; i++) {
	Node *theNode = theDomain->getNode(nodes(i));
	if (theNode == 0) {
	    opserr << "FATAL TransformationFE::TransformationFE() - no Node with tag: ";
	    opserr << nodes(i) << " in the domain\n";;
	    exit(-1);
	}
	DOF_Group *theDofGroup = theNode->getDOF_GroupPtr();
	if (theDofGroup == 0) {
	    opserr << "FATAL TransformationFE::TransformationFE() - no DOF_Group : ";
	    opserr << " associated with node: " << nodes(i) << " in the domain\n";;
	    exit(-1);
	}	
	theDOFs[i] = theDofGroup;
    }

    // see if theTransformation array is big enough
    // if not delete the old and create a new one
    if (numNodes > sizeTransformations) {
	if (theTransformations != 0) 
	    delete [] theTransformations;
	
	theTransformations = new Matrix *[numNodes];
	if (theTransformations == 0) {
	    opserr << "FATAL TransformationFE::TransformationFE() - out of memory ";
	    opserr << "for array of pointers for Transformation matrices of size ";
	    opserr << numNodes;
	    exit(-1);
	}		    
	sizeTransformations = numNodes;
    }	

    // if this is the first element of this type create the arrays for 
    // modified tangent and residual matrices
    if (numTransFE == 0) {

	modMatrices = new Matrix *[MAX_NUM_DOF+1];
	modVectors  = new Vector *[MAX_NUM_DOF+1];
	dataBuffer = new double[MAX_NUM_DOF*MAX_NUM_DOF];
	localKbuffer = new double[MAX_NUM_DOF*MAX_NUM_DOF];
	dofData      = new int[MAX_NUM_DOF];
	sizeBuffer = MAX_NUM_DOF*MAX_NUM_DOF;
	
	if (modMatrices == 0 || modVectors == 0 || dataBuffer == 0 ||
	    localKbuffer == 0 || dofData == 0) {
	    opserr << "TransformationFE::TransformationFE(Element *) ";
	    opserr << " ran out of memory";	    
	}
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    modMatrices[i] = 0;
	    modVectors[i] = 0;
	}
    }

    // increment the number of transformations
    numTransFE++;
}



// ~TransformationFE();    
//	destructor.

TransformationFE::~TransformationFE()
{

    numTransFE--;
    
    if (theDOFs != 0)
	delete [] theDOFs;
    if (theSPs != 0)
	delete [] theSPs;

    int numDOF = 0;    
    if (modID != 0)
	numDOF = modID->Size();
    
    if (modID != 0)
	delete modID;

    if (numDOF > MAX_NUM_DOF) {
	// tangent and residual may have been  created specially
	if (modTangent != 0) delete modTangent;
	if (modResidual != 0) delete modResidual;
    }

    // if this is the last FE_Element, clean up the
    // storage for the matrix and vector objects
    if (numTransFE == 0) {
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    if (modVectors[i] != 0)
		delete modVectors[i];
	    if (modMatrices[i] != 0)
		delete modMatrices[i];
	}
	delete [] modMatrices;
	delete [] modVectors;
	delete [] theTransformations;
	delete [] dataBuffer;
	delete [] localKbuffer;
	delete [] dofData;
	modMatrices = 0;
	modVectors = 0;
	theTransformations = 0;
	dataBuffer = 0;
	localKbuffer = 0;
	dofData = 0;
	sizeTransformations = 0;
	sizeBuffer = 0;
	transCounter = 0;
    }
}    


const ID &
TransformationFE::getDOFtags(void) const 
{
    return this->FE_Element::getDOFtags();
}


const ID &
TransformationFE::getID(void) const
{
    // make sure that it exists
    if (modID == 0) {
	opserr << "FATAL TransformationFE::getID() called before setID()\n";
	exit(-1);
    }
    return *modID;
}


int
TransformationFE::setID(void)
{
    // determine number of DOF
    numTransformedDOF = 0;
    for (int ii=0; ii<numGroups; ii++) {
	DOF_Group *dofPtr = theDOFs[ii];
	numTransformedDOF += dofPtr->getNumDOF();
    }

    // create an ID to hold the array, cannot use existing as 
    // may be different size
    if (modID != 0)
      delete modID;
    modID = 0;

    modID = new ID(numTransformedDOF);
    if (modID == 0 || modID->Size() == 0) {
	opserr << "TransformationFE::setID() ";
	opserr << " ran out of memory for ID of size :";
	opserr << numTransformedDOF << endln;
	exit(-1);
    }

    // fill in the ID
    int current = 0;
    for (int i=0; i<numGroups; i++) {
	DOF_Group *dofPtr = theDOFs[i];
	const ID &theDOFid = dofPtr->getID();

	for (int j=0; j<theDOFid.Size(); j++)  
	    if (current < numTransformedDOF)
		(*modID)(current++) = theDOFid(j);
	    else {
		opserr << "WARNING TransformationFE::setID() - numDOF and";
		opserr << " number of dof at the DOF_Groups\n";
		return -3;
	    }		
    }
    
    // set the pointers to the modified tangent matrix and residual vector
    if (numTransformedDOF <= MAX_NUM_DOF) {
	// use class wide objects
	if (modVectors[numTransformedDOF] == 0) {
	    modVectors[numTransformedDOF] = new Vector(numTransformedDOF);
	    modMatrices[numTransformedDOF] = new Matrix(numTransformedDOF,numTransformedDOF);
	    modResidual = modVectors[numTransformedDOF];
	    modTangent = modMatrices[numTransformedDOF];
	    if (modResidual == 0 || modResidual->Size() != numTransformedDOF ||	
		modTangent == 0 || modTangent->noCols() != numTransformedDOF)	{  
		opserr << "TransformationFE::setID() ";
		opserr << " ran out of memory for vector/Matrix of size :";
		opserr << numTransformedDOF << endln;
		exit(-1);
	    }
	} else {
	    modResidual = modVectors[numTransformedDOF];
	    modTangent = modMatrices[numTransformedDOF];
	}
    } else {
	// create matrices and vectors for each object instance
	modResidual = new Vector(numTransformedDOF);
	modTangent = new Matrix(numTransformedDOF, numTransformedDOF);
	if (modResidual == 0 || modResidual->Size() ==0 ||
	    modTangent ==0 || modTangent->noRows() ==0) {
	    
	    opserr << "TransformationFE::setID() ";
	    opserr << " ran out of memory for vector/Matrix of size :";
	    opserr << numTransformedDOF << endln;
	    exit(-1);
	}
    }     

    return 0;
}

const Matrix &
TransformationFE::getTangent(Integrator *theNewIntegrator)
{
    const Matrix &theTangent = this->FE_Element::getTangent(theNewIntegrator);

    static ID numDOFs(dofData, 1);
    numDOFs.setData(dofData, numGroups);
    
    // DO THE SP STUFF TO THE TANGENT 
    
    // get the transformation matrix from each dof group & number of local dof
    // for original node.
    int numNode = numGroups;
    for (int a = 0; a<numNode; a++) {
      Matrix *theT = theDOFs[a]->getT();
      theTransformations[a] = theT;
      if (theT != 0)
	numDOFs[a] = theT->noRows(); // T^ 
      else
	numDOFs[a] = theDOFs[a]->getNumDOF();
    }

    // perform Tt K T -- as T is block diagonal do T(i)^T K(i,j) T(j)
    // where blocks are of size equal to num ele dof at a node

    int startRow = 0;
    int noRowsTransformed = 0;
    int noRowsOriginal = 0;

    static Matrix localK;

    // foreach block row, for each block col do
    for (int i=0; i<numNode; i++) {

	int startCol = 0;
	int numDOFi = numDOFs[i];	
	int noColsOriginal = 0;

	for (int j=0; j<numNode; j++) {

	    const Matrix *Ti = theTransformations[i];
	    const Matrix *Tj = theTransformations[j];
	    int numDOFj = numDOFs[j];	
	    localK.setData(localKbuffer, numDOFi, numDOFj);

	    // copy K(i,j) into localK matrix
	    // CHECK SIZE OF BUFFFER	    
	    for (int a=0; a<numDOFi; a++)
		for (int b=0; b<numDOFj; b++)
		    localK(a,b) = theTangent(noRowsOriginal+a, noColsOriginal+b);

	    // now perform the matrix computation T(i)^T localK T(j)
	    // note: if T == 0 then the Identity is assumed
	    int noColsTransformed = 0;
	    static Matrix localTtKT;
	    
	    if (Ti != 0 && Tj != 0) {
		noRowsTransformed = Ti->noCols();
		noColsTransformed = Tj->noCols();
		// CHECK SIZE OF BUFFFER
		localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
		//localTtKT = (*Ti) ^ localK * (*Tj);
		localTtKT.addMatrixTripleProduct(0.0, *Ti, localK, *Tj, 1.0);
	    } else if (Ti == 0 && Tj != 0) {
		noRowsTransformed = numDOFi;
		noColsTransformed = Tj->noCols();
		// CHECK SIZE OF BUFFFER
		localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
		// localTtKT = localK * (*Tj);	       
		localTtKT.addMatrixProduct(0.0, localK, *Tj, 1.0);
	    } else if (Ti != 0 && Tj == 0) {
		noRowsTransformed = Ti->noCols();
		noColsTransformed = numDOFj;
		localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
		//localTtKT = (*Ti) ^ localK;
		localTtKT.addMatrixTransposeProduct(0.0, *Ti, localK, 1.0);
	    } else {
		noRowsTransformed = numDOFi;
		noColsTransformed = numDOFj;
		localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
		localTtKT = localK;
	    }
	    // now copy into modTangent the T(i)^t K(i,j) T(j) product
	    for (int c=0; c<noRowsTransformed; c++) 
		for (int d=0; d<noColsTransformed; d++) 
		    (*modTangent)(startRow+c, startCol+d) = localTtKT(c,d);
	    
	    startCol += noColsTransformed;
	    noColsOriginal += numDOFj;
	}

	noRowsOriginal += numDOFi;
	startRow += noRowsTransformed;
    }

    return *modTangent;
}


const Vector &
TransformationFE::getResidual(Integrator *theNewIntegrator)

{
    const Vector &theResidual = this->FE_Element::getResidual(theNewIntegrator);
    // DO THE SP STUFF TO THE TANGENT
    
    // perform Tt R  -- as T is block diagonal do T(i)^T R(i)
    // where blocks are of size equal to num ele dof at a node

    int startRowTransformed = 0;
    int startRowOriginal = 0;
    int numNode = numGroups;

    // foreach block row, for each block col do
    for (int i=0; i<numNode; i++) {
	int noRows = 0;
	int noCols = 0;
	const Matrix *Ti = theDOFs[i]->getT();
	if (Ti != 0) {
	  noRows = Ti->noCols(); // T^
	  noCols = Ti->noRows();

	  /*
	  Vector orig(noCols);
	  Vector mod(noRows);
	  for (int k=startRowOriginal; k<startRowOriginal+noCols; k++)
	    orig(k-startRowOriginal)= theResidual(k);
	  mod = (*Ti)^orig;
	  for (int k=startRowTransformed; k<startRowTransformed+noRows; k++)
	    (*modResidual)(k) = mod (k-startRowTransformed);

	  */

	  for (int j=0; j<noRows; j++) {
	    double sum = 0.0;
	    for (int k=0; k<noCols; k++)
	      sum += (*Ti)(k,j) * theResidual(startRowOriginal + k);
	    (*modResidual)(startRowTransformed +j) = sum;
	  }

	} else {
	  noCols = theDOFs[i]->getNumDOF();
	  noRows = noCols;
	  for (int j=0; j<noRows; j++)
	    (*modResidual)(startRowTransformed +j) = theResidual(startRowOriginal + j);
	}
	startRowTransformed += noRows;
	startRowOriginal += noCols;
    }

    return *modResidual;
}




const Vector &
TransformationFE::getTangForce(const Vector &disp, double fact)
{
    opserr << "TransformationFE::getTangForce() - not yet implemented\n";
    modResidual->Zero();
    return *modResidual;
}

const Vector &
TransformationFE::getK_Force(const Vector &accel, double fact)
{
  this->FE_Element::zeroTangent();    
  this->FE_Element::addKtToTang();    
  const Matrix &theTangent = this->FE_Element::getTangent(0);

  static ID numDOFs(dofData, 1);
  numDOFs.setData(dofData, numGroups);
    
  // DO THE SP STUFF TO THE TANGENT 
  
  // get the transformation matrix from each dof group & number of local dof
  // for original node.
  int numNode = numGroups;
  for (int a = 0; a<numNode; a++) {
    Matrix *theT = theDOFs[a]->getT();
    theTransformations[a] = theT;
    if (theT != 0)
      numDOFs[a] = theT->noRows(); // T^ 
    else
      numDOFs[a] = theDOFs[a]->getNumDOF();
  }
  
  // perform Tt K T -- as T is block diagonal do T(i)^T K(i,j) T(j)
  // where blocks are of size equal to num ele dof at a node
  
  int startRow = 0;
  int noRowsTransformed = 0;
  int noRowsOriginal = 0;
  
  static Matrix localK;
  
  // foreach block row, for each block col do
  for (int i=0; i<numNode; i++) {
    
    int startCol = 0;
    int numDOFi = numDOFs[i];	
    int noColsOriginal = 0;
    
    for (int j=0; j<numNode; j++) {
      
      const Matrix *Ti = theTransformations[i];
      const Matrix *Tj = theTransformations[j];
      int numDOFj = numDOFs[j];	
      localK.setData(localKbuffer, numDOFi, numDOFj);
      
      // copy K(i,j) into localK matrix
      // CHECK SIZE OF BUFFFER	    
      for (int a=0; a<numDOFi; a++)
	for (int b=0; b<numDOFj; b++)
	  localK(a,b) = theTangent(noRowsOriginal+a, noColsOriginal+b);
      
      // now perform the matrix computation T(i)^T localK T(j)
      // note: if T == 0 then the Identity is assumed
      int noColsTransformed = 0;
      static Matrix localTtKT;
      
      if (Ti != 0 && Tj != 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK * (*Tj);
	localTtKT.addMatrixTripleProduct(0.0, *Ti, localK, *Tj, 1.0);
      } else if (Ti == 0 && Tj != 0) {
	noRowsTransformed = numDOFi;
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	// localTtKT = localK * (*Tj);	       
	localTtKT.addMatrixProduct(0.0, localK, *Tj, 1.0);
      } else if (Ti != 0 && Tj == 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = numDOFj;
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK;
	localTtKT.addMatrixTransposeProduct(0.0, *Ti, localK, 1.0);
      } else {
	noRowsTransformed = numDOFi;
	noColsTransformed = numDOFj;
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	localTtKT = localK;
      }
      // now copy into modTangent the T(i)^t K(i,j) T(j) product
      for (int c=0; c<noRowsTransformed; c++) 
	for (int d=0; d<noColsTransformed; d++) 
	  (*modTangent)(startRow+c, startCol+d) = localTtKT(c,d);
      
      startCol += noColsTransformed;
      noColsOriginal += numDOFj;
    }
    
    noRowsOriginal += numDOFi;
    startRow += noRowsTransformed;
  }
  
  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numTransformedDOF);
  for (int j=0; j<numTransformedDOF; j++) {
    int dof = (*modID)(j);
    if (dof >= 0)
      tmp(j) = accel(dof);
    else
      tmp(j) = 0.0;
  }

  modResidual->addMatrixVector(0.0, *modTangent, tmp, 1.0);

  return *modResidual;
}


const Vector &
TransformationFE::getKi_Force(const Vector &accel, double fact)
{
  this->FE_Element::zeroTangent();    
  this->FE_Element::addKiToTang();    
  const Matrix &theTangent = this->FE_Element::getTangent(0);

  static ID numDOFs(dofData, 1);
  numDOFs.setData(dofData, numGroups);
    
  // DO THE SP STUFF TO THE TANGENT 
  
  // get the transformation matrix from each dof group & number of local dof
  // for original node.
  int numNode = numGroups;
  for (int a = 0; a<numNode; a++) {
    Matrix *theT = theDOFs[a]->getT();
    theTransformations[a] = theT;
    if (theT != 0)
      numDOFs[a] = theT->noRows(); // T^ 
    else
      numDOFs[a] = theDOFs[a]->getNumDOF();
  }
  
  // perform Tt K T -- as T is block diagonal do T(i)^T K(i,j) T(j)
  // where blocks are of size equal to num ele dof at a node
  
  int startRow = 0;
  int noRowsTransformed = 0;
  int noRowsOriginal = 0;
  
  static Matrix localK;
  
  // foreach block row, for each block col do
  for (int i=0; i<numNode; i++) {
    
    int startCol = 0;
    int numDOFi = numDOFs[i];	
    int noColsOriginal = 0;
    
    for (int j=0; j<numNode; j++) {
      
      const Matrix *Ti = theTransformations[i];
      const Matrix *Tj = theTransformations[j];
      int numDOFj = numDOFs[j];	
      localK.setData(localKbuffer, numDOFi, numDOFj);
      
      // copy K(i,j) into localK matrix
      // CHECK SIZE OF BUFFFER	    
      for (int a=0; a<numDOFi; a++)
	for (int b=0; b<numDOFj; b++)
	  localK(a,b) = theTangent(noRowsOriginal+a, noColsOriginal+b);
      
      // now perform the matrix computation T(i)^T localK T(j)
      // note: if T == 0 then the Identity is assumed
      int noColsTransformed = 0;
      static Matrix localTtKT;
      
      if (Ti != 0 && Tj != 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK * (*Tj);
	localTtKT.addMatrixTripleProduct(0.0, *Ti, localK, *Tj, 1.0);
      } else if (Ti == 0 && Tj != 0) {
	noRowsTransformed = numDOFi;
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	// localTtKT = localK * (*Tj);	       
	localTtKT.addMatrixProduct(0.0, localK, *Tj, 1.0);
      } else if (Ti != 0 && Tj == 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = numDOFj;
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK;
	localTtKT.addMatrixTransposeProduct(0.0, *Ti, localK, 1.0);
      } else {
	noRowsTransformed = numDOFi;
	noColsTransformed = numDOFj;
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	localTtKT = localK;
      }
      // now copy into modTangent the T(i)^t K(i,j) T(j) product
      for (int c=0; c<noRowsTransformed; c++) 
	for (int d=0; d<noColsTransformed; d++) 
	  (*modTangent)(startRow+c, startCol+d) = localTtKT(c,d);
      
      startCol += noColsTransformed;
      noColsOriginal += numDOFj;
    }
    
    noRowsOriginal += numDOFi;
    startRow += noRowsTransformed;
  }
  
  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numTransformedDOF);
  for (int j=0; j<numTransformedDOF; j++) {
    int dof = (*modID)(j);
    if (dof >= 0)
      tmp(j) = accel(dof);
    else
      tmp(j) = 0.0;
  }

  modResidual->addMatrixVector(0.0, *modTangent, tmp, 1.0);

  return *modResidual;
}

const Vector &
TransformationFE::getM_Force(const Vector &accel, double fact)
{
  this->FE_Element::zeroTangent();    
  this->FE_Element::addMtoTang();    
  const Matrix &theTangent = this->FE_Element::getTangent(0);

  static ID numDOFs(dofData, 1);
  numDOFs.setData(dofData, numGroups);
    
  // DO THE SP STUFF TO THE TANGENT 
  
  // get the transformation matrix from each dof group & number of local dof
  // for original node.
  int numNode = numGroups;
  for (int a = 0; a<numNode; a++) {
    Matrix *theT = theDOFs[a]->getT();
    theTransformations[a] = theT;
    if (theT != 0)
      numDOFs[a] = theT->noRows(); // T^ 
    else
      numDOFs[a] = theDOFs[a]->getNumDOF();
  }
  
  // perform Tt K T -- as T is block diagonal do T(i)^T K(i,j) T(j)
  // where blocks are of size equal to num ele dof at a node
  
  int startRow = 0;
  int noRowsTransformed = 0;
  int noRowsOriginal = 0;
  
  static Matrix localK;
  
  // foreach block row, for each block col do
  for (int i=0; i<numNode; i++) {
    
    int startCol = 0;
    int numDOFi = numDOFs[i];	
    int noColsOriginal = 0;
    
    for (int j=0; j<numNode; j++) {
      
      const Matrix *Ti = theTransformations[i];
      const Matrix *Tj = theTransformations[j];
      int numDOFj = numDOFs[j];	
      localK.setData(localKbuffer, numDOFi, numDOFj);
      
      // copy K(i,j) into localK matrix
      // CHECK SIZE OF BUFFFER	    
      for (int a=0; a<numDOFi; a++)
	for (int b=0; b<numDOFj; b++)
	  localK(a,b) = theTangent(noRowsOriginal+a, noColsOriginal+b);
      
      // now perform the matrix computation T(i)^T localK T(j)
      // note: if T == 0 then the Identity is assumed
      int noColsTransformed = 0;
      static Matrix localTtKT;
      
      if (Ti != 0 && Tj != 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK * (*Tj);
	localTtKT.addMatrixTripleProduct(0.0, *Ti, localK, *Tj, 1.0);
      } else if (Ti == 0 && Tj != 0) {
	noRowsTransformed = numDOFi;
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	// localTtKT = localK * (*Tj);	       
	localTtKT.addMatrixProduct(0.0, localK, *Tj, 1.0);
      } else if (Ti != 0 && Tj == 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = numDOFj;
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK;
	localTtKT.addMatrixTransposeProduct(0.0, *Ti, localK, 1.0);
      } else {
	noRowsTransformed = numDOFi;
	noColsTransformed = numDOFj;
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	localTtKT = localK;
      }
      // now copy into modTangent the T(i)^t K(i,j) T(j) product
      for (int c=0; c<noRowsTransformed; c++) 
	for (int d=0; d<noColsTransformed; d++) 
	  (*modTangent)(startRow+c, startCol+d) = localTtKT(c,d);
      
      startCol += noColsTransformed;
      noColsOriginal += numDOFj;
    }
    
    noRowsOriginal += numDOFi;
    startRow += noRowsTransformed;
  }
  
  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numTransformedDOF);
  for (int j=0; j<numTransformedDOF; j++) {
    int dof = (*modID)(j);
    if (dof >= 0)
      tmp(j) = accel(dof);
    else
      tmp(j) = 0.0;
  }

  modResidual->addMatrixVector(0.0, *modTangent, tmp, 1.0);

  return *modResidual;
}

const Vector &
TransformationFE::getC_Force(const Vector &accel, double fact)
{
  this->FE_Element::zeroTangent();    
  this->FE_Element::addCtoTang();    
  const Matrix &theTangent = this->FE_Element::getTangent(0);

  static ID numDOFs(dofData, 1);
  numDOFs.setData(dofData, numGroups);
    
  // DO THE SP STUFF TO THE TANGENT 
  
  // get the transformation matrix from each dof group & number of local dof
  // for original node.
  int numNode = numGroups;
  for (int a = 0; a<numNode; a++) {
    Matrix *theT = theDOFs[a]->getT();
    theTransformations[a] = theT;
    if (theT != 0)
      numDOFs[a] = theT->noRows(); // T^ 
    else
      numDOFs[a] = theDOFs[a]->getNumDOF();
  }
  
  // perform Tt K T -- as T is block diagonal do T(i)^T K(i,j) T(j)
  // where blocks are of size equal to num ele dof at a node
  
  int startRow = 0;
  int noRowsTransformed = 0;
  int noRowsOriginal = 0;
  
  static Matrix localK;
  
  // foreach block row, for each block col do
  for (int i=0; i<numNode; i++) {
    
    int startCol = 0;
    int numDOFi = numDOFs[i];	
    int noColsOriginal = 0;
    
    for (int j=0; j<numNode; j++) {
      
      const Matrix *Ti = theTransformations[i];
      const Matrix *Tj = theTransformations[j];
      int numDOFj = numDOFs[j];	
      localK.setData(localKbuffer, numDOFi, numDOFj);
      
      // copy K(i,j) into localK matrix
      // CHECK SIZE OF BUFFFER	    
      for (int a=0; a<numDOFi; a++)
	for (int b=0; b<numDOFj; b++)
	  localK(a,b) = theTangent(noRowsOriginal+a, noColsOriginal+b);
      
      // now perform the matrix computation T(i)^T localK T(j)
      // note: if T == 0 then the Identity is assumed
      int noColsTransformed = 0;
      static Matrix localTtKT;
      
      if (Ti != 0 && Tj != 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK * (*Tj);
	localTtKT.addMatrixTripleProduct(0.0, *Ti, localK, *Tj, 1.0);
      } else if (Ti == 0 && Tj != 0) {
	noRowsTransformed = numDOFi;
	noColsTransformed = Tj->noCols();
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	// localTtKT = localK * (*Tj);	       
	localTtKT.addMatrixProduct(0.0, localK, *Tj, 1.0);
      } else if (Ti != 0 && Tj == 0) {
	noRowsTransformed = Ti->noCols();
	noColsTransformed = numDOFj;
	// CHECK SIZE OF BUFFFER
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	//localTtKT = (*Ti) ^ localK;
	localTtKT.addMatrixTransposeProduct(0.0, *Ti, localK, 1.0);
      } else {
	noRowsTransformed = numDOFi;
	noColsTransformed = numDOFj;
	localTtKT.setData(dataBuffer, noRowsTransformed, noColsTransformed);
	localTtKT = localK;
      }
      // now copy into modTangent the T(i)^t K(i,j) T(j) product
      for (int c=0; c<noRowsTransformed; c++) 
	for (int d=0; d<noColsTransformed; d++) 
	  (*modTangent)(startRow+c, startCol+d) = localTtKT(c,d);
      
      startCol += noColsTransformed;
      noColsOriginal += numDOFj;
    }
    
    noRowsOriginal += numDOFi;
    startRow += noRowsTransformed;
  }
  
  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numTransformedDOF);
  for (int j=0; j<numTransformedDOF; j++) {
    int dof = (*modID)(j);
    if (dof >= 0)
      tmp(j) = accel(dof);
    else
      tmp(j) = 0.0;
  }

  modResidual->addMatrixVector(0.0, *modTangent, tmp, 1.0);

  return *modResidual;
}


void  
TransformationFE::addD_Force(const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    static Vector response;
    response.setData(dataBuffer, numOriginalDOF);
		    
    for (int i=0; i<numTransformedDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalD_Force(response, fact);
}   	 

void  
TransformationFE::addM_Force(const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    static Vector response;
    response.setData(dataBuffer, numOriginalDOF);
		    
    for (int i=0; i<numTransformedDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalM_Force(response, fact);
}   	 



// CHANGE THE ID SENT
const Vector &
TransformationFE::getLastResponse(void)
{
    Integrator *theLastIntegrator = this->getLastIntegrator();
    if (theLastIntegrator != 0) {
	if (theLastIntegrator->getLastResponse(*modResidual,*modID) < 0) {
	    opserr << "WARNING TransformationFE::getLastResponse(void)";
	    opserr << " - the Integrator had problems with getLastResponse()\n";
	}
    }
    else {
	modResidual->Zero();
	opserr << "WARNING  TransformationFE::getLastResponse()";
	opserr << " No Integrator yet passed\n";
    }
    
    Vector &result = *modResidual;
    return result;
}


int 
TransformationFE::transformResponse(const Vector &modResp, 
				    Vector &unmodResp)
{
    // perform T R  -- as T is block diagonal do T(i) R(i)
    // where blocks are of size equal to num ele dof at a node

    int startRowOriginal = 0;
    int startRowTransformed = 0;
    int numNode = numGroups;
    int noRows = 0;
    int noCols = 0;

    for (int i=0; i<numNode; i++) {
	const Matrix *Ti = theDOFs[i]->getT();
	if (Ti != 0) {
	    noRows = Ti->noRows();
	    noCols = Ti->noCols();
	    for (int j=0; j<noRows; j++) {
		double sum = 0.0;
		for (int k=0; k<noCols; k++)
		    sum += (*Ti)(j,k) * modResp(startRowTransformed +k) ;
		unmodResp(startRowOriginal + j) = sum;
	    }
	} else {
	    noCols = theDOFs[i]->getNumDOF();
	    noRows = noCols;
	    for (int j=0; j<noCols; j++)
		unmodResp(startRowOriginal + j) = modResp(startRowTransformed +j);
	}
	startRowOriginal += noRows;
	startRowTransformed += noCols;
    }

    return 0;
}


// AddingSensitivity:BEGIN /////////////////////////////////
void  
TransformationFE::addD_ForceSensitivity(int gradNumber, const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    static Vector response;
    response.setData(dataBuffer, numOriginalDOF);
		    
    for (int i=0; i<numTransformedDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalD_ForceSensitivity(gradNumber, response, fact);
}   	 

void  
TransformationFE::addM_ForceSensitivity(int gradNumber, const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    static Vector response;
    response.setData(dataBuffer, numOriginalDOF);
		    
    for (int i=0; i<numTransformedDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalM_ForceSensitivity(gradNumber, response, fact);
}   	 

// AddingSensitivity:END ////////////////////////////////////
