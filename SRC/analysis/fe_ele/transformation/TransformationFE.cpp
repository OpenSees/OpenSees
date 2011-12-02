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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/transformation/TransformationFE.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/TransformationFE.C
//
// Written: fmk 

// Created: 05/99
// Revision: A
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
double *TransformationFE::dataBuffer;          
double *TransformationFE::localKbuffer;          
int TransformationFE::sizeBuffer(0);            

//  TransformationFE(Element *, Integrator *theIntegrator);
//	construictor that take the corresponding model element.
TransformationFE::TransformationFE(Element *ele,
				   TransformationConstraintHandler &theH) 
:FE_Element(ele), theDOFs(0), numSPs(0), theSPs(0), modID(0), 
 modTangent(0), modResidual(0), numGroups(0), numEleNodalDOF(0), theHandler(&theH)
{
    // create the array of pointers to DOF_Groups
    const ID &nodes = ele->getExternalNodes();
    Domain *theDomain = ele->getDomain();
    int numNodes = nodes.Size();
    theDOFs = new DOF_Group *[numNodes];
    if (theDOFs == 0) {
	cerr << "FATAL TransformationFE::TransformationFE() - out of memory craeting ";
	cerr << "array of size : " << numNodes << " for storage of DOF_Group\n";
	exit(-1);
    }

    numGroups = numNodes;
    numEleNodalDOF = ele->getNumDOF()/numNodes;

    // now fill the array of DOF_Group pointers
    for (int i=0; i<numNodes; i++) {
	Node *theNode = theDomain->getNode(nodes(i));
	if (theNode == 0) {
	    cerr << "FATAL TransformationFE::TransformationFE() - no Node with tag: ";
	    cerr << nodes(i) << " in the domain\n";;
	    exit(-1);
	}
	DOF_Group *theDofGroup = theNode->getDOF_GroupPtr();
	if (theDofGroup == 0) {
	    cerr << "FATAL TransformationFE::TransformationFE() - no DOF_Group : ";
	    cerr << " associated with node: " << nodes(i) << " in the domain\n";;
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
	    cerr << "FATAL TransformationFE::TransformationFE() - out of memory ";
	    cerr << "for array of pointers for Transformation matrices of size ";
	    cerr << numNodes;
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
	sizeBuffer = MAX_NUM_DOF*MAX_NUM_DOF;
	
	if (modMatrices == 0 || modVectors == 0) {
	    cerr << "TransformationFE::TransformationFE(Element *) ";
	    cerr << " ran out of memory";	    
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
	modMatrices = 0;
	modVectors = 0;
	theTransformations = 0;
	dataBuffer = 0;
	localKbuffer = 0;
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
	cerr << "FATAL TransformationFE::getID() called before setID()\n";
	exit(-1);
    }
    return *modID;
}


int
TransformationFE::setID(void)
{
    // get the TransformationDOF_Groups to do their ID stuff
    if (transCounter == 0) {
	transCounter++;
	theHandler->doneDOFids();
    } else if (transCounter == numTransFE)
	transCounter = 0;    
    else
	transCounter++;

    // determine number of DOF
    int numDOF = 0;
    for (int ii=0; ii<numGroups; ii++) {
	DOF_Group *dofPtr = theDOFs[ii];
	numDOF += dofPtr->getNumDOF();
    }

    // create an ID to hold the array, cannot use existing as 
    // may be different size
    modID = new ID(numDOF);
    if (modID == 0 || modID->Size() == 0) {
	cerr << "TransformationFE::setID() ";
	cerr << " ran out of memory for ID of size :";
	cerr << numDOF << endl;
	exit(-1);
    }

    // fill in the ID
    int current = 0;
    for (int i=0; i<numGroups; i++) {
	DOF_Group *dofPtr = theDOFs[i];
	const ID &theDOFid = dofPtr->getID();

	for (int j=0; j<theDOFid.Size(); j++)  
	    if (current < numDOF)
		(*modID)(current++) = theDOFid(j);
	    else {
		cerr << "WARNING TransformationFE::setID() - numDOF and";
		cerr << " number of dof at the DOF_Groups\n";
		return -3;
	    }		
    }
    
    // set the pointers to the modified tangent matrix and residual vector
    if (numDOF <= MAX_NUM_DOF) {
	// use class wide objects
	if (modVectors[numDOF] == 0) {
	    modVectors[numDOF] = new Vector(numDOF);
	    modMatrices[numDOF] = new Matrix(numDOF,numDOF);
	    modResidual = modVectors[numDOF];
	    modTangent = modMatrices[numDOF];
	    if (modResidual == 0 || modResidual->Size() != numDOF ||	
		modTangent == 0 || modTangent->noCols() != numDOF)	{  
		cerr << "TransformationFE::setID() ";
		cerr << " ran out of memory for vector/Matrix of size :";
		cerr << numDOF << endl;
		exit(-1);
	    }
	} else {
	    modResidual = modVectors[numDOF];
	    modTangent = modMatrices[numDOF];
	}
    } else {
	// create matrices and vectors for each object instance
	modResidual = new Vector(numDOF);
	modTangent = new Matrix(numDOF, numDOF);
	if (modResidual == 0 || modResidual->Size() ==0 ||
	    modTangent ==0 || modTangent->noRows() ==0) {
	    
	    cerr << "TransformationFE::setID() ";
	    cerr << " ran out of memory for vector/Matrix of size :";
	    cerr << numDOF << endl;
	    exit(-1);
	}
    }     

    return 0;
}

const Matrix &
TransformationFE::getTangent(Integrator *theNewIntegrator)
{
    // impose the SP_Constraints at the nodes
    if (transCounter == 0) {
	transCounter++;
	theHandler->enforceSPs();
    } else if (transCounter == numTransFE)
	transCounter = 0;
    else
	transCounter++;
	    
    const Matrix &theTangent = this->FE_Element::getTangent(theNewIntegrator);

    // DO THE SP STUFF TO THE TANGENT 
    
    // get the transformation matrix from each dof group
    int numNode = numGroups;
    for (int a = 0; a<numNode; a++)
	theTransformations[a] = theDOFs[a]->getT();

    // perform Tt K T -- as T is block diagonal do T(i)^T K(i,j) T(j)
    // where blocks are of size equal to num ele dof at a node
    int startRow = 0;

    Matrix *localK = new Matrix(localKbuffer, numEleNodalDOF, numEleNodalDOF);
    
    // foreach block row, for each block col do
    for (int i=0; i<numNode; i++) {
	int noRows = 0;
	int startCol = 0;
	
	for (int j=0; j<numNode; j++) {
	    
	    const Matrix *Ti = theTransformations[i];
	    const Matrix *Tj = theTransformations[j];
	    
	    // copy K(i,j) into localK matrix
	    // CHECK SIZE OF BUFFFER	    
	    for (int a=0; a<numEleNodalDOF; a++)
		for (int b=0; b<numEleNodalDOF; b++)
		    (*localK)(a,b) = theTangent(i*numEleNodalDOF+a, j*numEleNodalDOF+b);
	    
	    // now perform the matrix computation T(i)^T localK T(j)
	    // note: if T == 0 then the Identity is assumed
	    int noCols = 0;
	    Matrix *localTtKT;
	    
	    if (Ti != 0 && Tj != 0) {
		noRows = Ti->noCols();
		noCols = Tj->noCols();
		// CHECK SIZE OF BUFFFER
		localTtKT = new Matrix(dataBuffer, noRows, noCols);
		*localTtKT = (*Ti)^(*localK)*(*Tj);
	    } else if (Ti== 0 && Tj != 0) {
		noRows = numEleNodalDOF;
		noCols = Tj->noCols();
		// CHECK SIZE OF BUFFFER
		localTtKT = new Matrix(dataBuffer, noRows, noCols);
		*localTtKT = (*localK)*(*Tj);		
	    } else if (Ti != 0 && Tj == 0) {
		noRows = Ti->noCols();
		noCols = numEleNodalDOF;
		// CHECK SIZE OF BUFFFER
		localTtKT = new Matrix(dataBuffer, noRows, noCols);
		*localTtKT = (*Ti)^(*localK);
	    } else {
		noRows = numEleNodalDOF;
		noCols = numEleNodalDOF;
		localTtKT = localK;
	    }
	    // now copy into modTangent the T(i)^t K(i,j) T(j) product
	    for (int c=0; c<noRows; c++) 
		for (int d=0; d<noCols; d++) 
		    (*modTangent)(startRow+c, startCol+d) = (*localTtKT)(c,d);
	    
	    if (localTtKT  != localK)
		delete localTtKT;
	    
	    startCol += noCols;
	}
	startRow += noRows;
    }

    delete localK;
    return *modTangent;
}


const Vector &
TransformationFE::getResidual(Integrator *theNewIntegrator)
{
    // impose the SP_Constraints at the nodes
    if (transCounter == 0) {
	transCounter++;
	theHandler->enforceSPs();
    } else if (transCounter == numTransFE)
	transCounter = 0;
    else
	transCounter++;
	    
    const Vector &theResidual = this->FE_Element::getResidual(theNewIntegrator);

    // DO THE SP STUFF TO THE TANGENT
    
    // perform Tt R  -- as T is block diagonal do T(i)^T R(i)
    // where blocks are of size equal to num ele dof at a node

    int startRow = 0;
    int numNode = numGroups;
    
    // foreach block row, for each block col do
    for (int i=0; i<numNode; i++) {
	int noRows = 0;
	int noCols = 0;
	const Matrix *Ti = theDOFs[i]->getT();
	if (Ti != 0) {
	    noRows = Ti->noCols();
	    noCols = numEleNodalDOF;
	    for (int j=0; j<noRows; j++) {
		double sum = 0;
		for (int k=0; k<noCols; k++)
		    sum += (*Ti)(k,j) * theResidual(i*numEleNodalDOF + k);
		(*modResidual)(startRow +j) = sum;
	    }
	} else {
	    noRows = numEleNodalDOF;
	    for (int j=0; j<noRows; j++)
		(*modResidual)(startRow +j) = theResidual(i*numEleNodalDOF + j);
	}
	startRow += noRows;
    }
    return *modResidual;
}




const Vector &
TransformationFE::getTangForce(const Vector &disp, double fact)
{
    cerr << "TransformationFE::getTangForce() - not yet implemented\n";
    modResidual->Zero();
    return *modResidual;
}

const Vector &
TransformationFE::getKtForce(const Vector &disp, double fact)
{
    cerr << "TransformationFE::getKtForce() - not yet implemented\n";
    modResidual->Zero();
    return *modResidual;
}


const Vector &
TransformationFE::getKsForce(const Vector &disp, double fact)
{
    cerr << "TransformationFE::getKsForce() - not yet implemented\n";
    modResidual->Zero();
    return *modResidual;
}



const Vector &
TransformationFE::getD_Force(const Vector &vel, double fact)
{
    cerr << "TransformationFE::getD_Force() - not yet implemented\n";
    modResidual->Zero();
    return *modResidual;
}




const Vector &
TransformationFE::getM_Force(const Vector &accel, double fact)
{
  cerr << "TransformationFE::getM_Force() - not yet implemented\n";
  modResidual->Zero();
  return *modResidual;
}


// CHANGE THE ID SENT
const Vector &
TransformationFE::getLastResponse(void)
{
    Integrator *theLastIntegrator = this->getLastIntegrator();
    if (theLastIntegrator != 0) {
	if (theLastIntegrator->getLastResponse(*modResidual,*modID) < 0) {
	    cerr << "WARNING TransformationFE::getLastResponse(void)";
	    cerr << " - the Integrator had problems with getLastResponse()\n";
	}
    }
    else {
	modResidual->Zero();
	cerr << "WARNING  TransformationFE::getLastResponse()";
	cerr << " No Integrator yet passed\n";
    }
    
    Vector &result = *modResidual;
    return result;
}


void  
TransformationFE::addKtForce(const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    int size = numGroups*numEleNodalDOF;
    Vector response(dataBuffer, size);
		    
    int numDOF = modID->Size();    
    for (int i=0; i<numDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalKtForce(response, fact);
}   	 


void  
TransformationFE::addKsForce(const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    int size = numGroups*numEleNodalDOF;
    Vector response(dataBuffer, size);
		    
    int numDOF = modID->Size();    
    for (int i=0; i<numDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalKsForce(response, fact);
}   	 

void  
TransformationFE::addKiForce(const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    int size = numGroups*numEleNodalDOF;
    Vector response(dataBuffer, size);
		    
    int numDOF = modID->Size();    
    for (int i=0; i<numDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalKiForce(response, fact);
}   	 


void  
TransformationFE::addKcForce(const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    int size = numGroups*numEleNodalDOF;
    Vector response(dataBuffer, size);
		    
    int numDOF = modID->Size();    
    for (int i=0; i<numDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalKcForce(response, fact);
}   	 


void  
TransformationFE::addD_Force(const Vector &disp,  double fact)
{
    if (fact == 0.0)
	return;

    int size = numGroups*numEleNodalDOF;
    Vector response(dataBuffer, size);
		    
    int numDOF = modID->Size();    
    for (int i=0; i<numDOF; i++) {
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

    int size = numGroups*numEleNodalDOF;
    Vector response(dataBuffer, size);
		    
    int numDOF = modID->Size();    
    for (int i=0; i<numDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modResidual)(i) = disp(loc);
	else
	    (*modResidual)(i) = 0.0;
    }
    transformResponse(*modResidual, response);
    this->addLocalM_Force(response, fact);
}   	 

int 
TransformationFE::transformResponse(const Vector &modResp, 
				    Vector &unmodResp)
{
    // perform T R  -- as T is block diagonal do T(i) R(i)
    // where blocks are of size equal to num ele dof at a node

    int startRow = 0;
    int numNode = numGroups;

    for (int i=0; i<numNode; i++) {
	int noRows = 0;
	int noCols = 0;
	const Matrix *Ti = theDOFs[i]->getT();
	if (Ti != 0) {
	    noRows = numEleNodalDOF;
	    noCols = Ti->noCols();
	    for (int j=0; j<noRows; j++) {
		double sum = 0.0;
		for (int k=0; k<noCols; k++)
		    sum += (*Ti)(j,k) * modResp(startRow +k) ;
		unmodResp(i*numEleNodalDOF + j) = sum;
	    }
	} else {
	    noCols = numEleNodalDOF;
	    for (int j=0; j<noCols; j++)
		unmodResp(i*numEleNodalDOF + j) = modResp(startRow +j);
	}
	startRow += noCols;
    }
    
    return 0;
}
