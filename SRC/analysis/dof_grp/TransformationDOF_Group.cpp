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
                                                                        
// $Revision: 1.22 $
// $Date: 2009-03-09 22:33:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/dof_grp/TransformationDOF_Group.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 05/99
// Revision: A
//
// Purpose: This file contains the code for implementing the methods
// of the TransformationDOF_Group class interface.
//
// What: "@(#) TransformationDOF_Group.C, revA"

#include <TransformationDOF_Group.h>
#include <stdlib.h>

#include <Domain.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <TransientIntegrator.h>
#include <MP_Constraint.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <TransformationConstraintHandler.h>

#define MAX_NUM_DOF 16

// static variables initialisation
Matrix **TransformationDOF_Group::modMatrices; 
Vector **TransformationDOF_Group::modVectors;  
int TransformationDOF_Group::numTransDOFs(0);     // number of objects
TransformationConstraintHandler *TransformationDOF_Group::theHandler = 0;     // number of objects

TransformationDOF_Group::TransformationDOF_Group(int tag, Node *node, 
						 MP_Constraint *mp,
						 TransformationConstraintHandler *theTHandler)  
:DOF_Group(tag,node),
 theMP(mp),Trans(0),modTangent(0),modUnbalance(0),modID(0),theSPs(0)
{
    // determine the number of DOF 
    int numNodalDOF = node->getNumberDOF();
    const ID &retainedDOF = mp->getRetainedDOFs();
    const ID &constrainedDOF = mp->getConstrainedDOFs();    
    int numNodalDOFConstrained = constrainedDOF.Size();
    numConstrainedNodeRetainedDOF = numNodalDOF - numNodalDOFConstrained;
    needRetainedData = -1;

    int numRetainedNodeDOF = retainedDOF.Size();

    modNumDOF = numConstrainedNodeRetainedDOF + numRetainedNodeDOF;

    // create space for the SP_Constraint array
    theSPs = new SP_Constraint *[numNodalDOF];
    for (int ii=0; ii<numNodalDOF; ii++) 
	theSPs[ii] = 0;

    /***********************
    // set the SP_Constraint corresponding to the dof in modID
    Domain *theDomain=node->getDomain();
    int nodeTag = node->getTag();    
    SP_ConstraintIter &theSPIter = theDomain->getSPs();
    SP_Constraint *sp;
    while ((sp = theSPIter()) != 0) {
	if (sp->getNodeTag() == nodeTag) {
	    int dof = sp->getDOF_Number();
	    int loc = 0;
	    for (int i=0; i<dof; i++) 
		if (constrainedDOF.getLocation(i) < 0)
		    loc++;
	    theSPs[loc] = sp;
	}
    }
    *******************/
     
    // create ID and transformation matrix
    modID = new ID(modNumDOF);
    Trans = new Matrix(numNodalDOF, modNumDOF);

    if (modID == 0 || modID->Size() == 0 ||
	Trans == 0 || Trans->noRows() == 0) {
	
	opserr << "FATAL TransformationDOF_Group::TransformationDOF_Group() -";
	opserr << " ran out of memory for size: " << modNumDOF << endln;
	exit(-1);
    }

    // initially set the id values to -2 for any dof still due to constrained node
    for (int i=0; i<numConstrainedNodeRetainedDOF; i++)
	(*modID)(i) = -2;
    
    // for all the constrained dof values set to -1
    for (int j=numConstrainedNodeRetainedDOF; j<modNumDOF; j++)
	(*modID)(j) = -1;

    // for all the dof corresponding to the retained node set initially to -1
    // we don't initially assign these equation nos. - this is done in doneID()
    for (int k=numConstrainedNodeRetainedDOF; k<modNumDOF; k++)
	(*modID)(k) = -1;
    
    // if this is the first TransformationDOF_Group we now
    // create the arrays used to store pointers to class wide
    // matrix and vector objects used to return modTangent and residual
    if (numTransDOFs == 0) {
	modMatrices = new Matrix *[MAX_NUM_DOF+1];
	modVectors  = new Vector *[MAX_NUM_DOF+1];
	
	if (modMatrices == 0 || modVectors == 0) {
	    opserr << "TransformationDOF_Group::TransformationDOF_Group(Node *) ";
	    opserr << " ran out of memory";	    
	}
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    modMatrices[i] = 0;
	    modVectors[i] = 0;
	}
    }    

    // set the pointers for the modTangent and residual
    if (modNumDOF <= MAX_NUM_DOF) {
	// use class wide objects
	if (modVectors[modNumDOF] == 0) {
	    // have to create matrix and vector of size as none yet created
	    modVectors[modNumDOF] = new Vector(modNumDOF);
	    modMatrices[modNumDOF] = new Matrix(modNumDOF,modNumDOF);
	    modUnbalance = modVectors[modNumDOF];
	    modTangent = modMatrices[modNumDOF];
	    if (modUnbalance == 0 || modUnbalance->Size() != modNumDOF ||	
		modTangent == 0 || modTangent->noCols() != modNumDOF)	{  
		opserr << "TransformationDOF_Group::TransformationDOF_Group(Node *) ";
		opserr << " ran out of memory for vector/Matrix of size :";
		opserr << modNumDOF << endln;
		exit(-1);
	    }
	} else {
	    modUnbalance = modVectors[modNumDOF];
	    modTangent = modMatrices[modNumDOF];
	}
    } else {
	// create matrices and vectors for each object instance
	modUnbalance = new Vector(modNumDOF);
	modTangent = new Matrix(modNumDOF, modNumDOF);
	if (modUnbalance == 0 || modTangent ==0 ||
	    modTangent ==0 || modTangent->noRows() ==0) {
	    
	    opserr << "TransformationDOF_Group::TransformationDOF_Group(Node *) ";
	    opserr << " ran out of memory for vector/Matrix of size :";
	    opserr << modNumDOF << endln;
	    exit(-1);
	}
    }
    
    numTransDOFs++;
    theHandler = theTHandler;
}

void 
TransformationDOF_Group::setID(int dof, int value)
{
    if (theMP == 0)
	this->DOF_Group::setID(dof,value);

    else
	(*modID)(dof) = value;
}
	

TransformationDOF_Group::TransformationDOF_Group(int tag, 
						 Node *node, 
						 TransformationConstraintHandler *theTHandler)
:DOF_Group(tag,node),
 theMP(0),Trans(0),modTangent(0),modUnbalance(0),modID(0),theSPs(0) 
{
    modNumDOF = node->getNumberDOF();
    // create space for the SP_Constraint array
    theSPs = new SP_Constraint *[modNumDOF];
    for (int i=0; i<modNumDOF; i++) 
	theSPs[i] = 0;

    // set the SP_Constraint corresponding to the dof in myID
    Domain *theDomain=node->getDomain();
    int nodeTag = node->getTag();
    SP_ConstraintIter &theSPIter = theDomain->getSPs();
    SP_Constraint *sp;
    while ((sp = theSPIter()) != 0) {
	if (sp->getNodeTag() == nodeTag) {
	    int dof = sp->getDOF_Number();
	    theSPs[dof] = sp;
	}
    }    
    
    // if this is the first TransformationDOF_Group we now
    // create the arrays used to store pointers to class wide
    // matrix and vector objects used to return modTangent and residual
    if (numTransDOFs == 0) {
	modMatrices = new Matrix *[MAX_NUM_DOF+1];
	modVectors  = new Vector *[MAX_NUM_DOF+1];
	
	if (modMatrices == 0 || modVectors == 0) {
	    opserr << "TransformationDOF_Group::TransformationDOF_Group(Node *) ";
	    opserr << " ran out of memory";	    
	}
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    modMatrices[i] = 0;
	    modVectors[i] = 0;
	}
    }    

    numTransDOFs++;
    theHandler = theTHandler;
}


// ~TransformationDOF_Group();    
//	destructor.

TransformationDOF_Group::~TransformationDOF_Group()
{
    numTransDOFs--;

    // delete modTangent and residual if created specially
    if (modNumDOF > MAX_NUM_DOF) {
	if (modTangent != 0) delete modTangent;
	if (modUnbalance != 0) delete modUnbalance;
    }
    
    if (modID != 0) delete modID;
    if (Trans != 0) delete Trans;
    if (theSPs != 0) delete [] theSPs;

    // if this is the last FE_Element, clean up the
    // storage for the matrix and vector objects
    if (numTransDOFs == 0) {
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    if (modVectors[i] != 0)
		delete modVectors[i];
	    if (modMatrices[i] != 0)
		delete modMatrices[i];
	}	
	delete [] modMatrices;
	delete [] modVectors;
    }    
}    


const ID &
TransformationDOF_Group::getID(void) const
{
  if (modID != 0) 
	return *modID;
    else
	return this->DOF_Group::getID();
}

int
TransformationDOF_Group::getNumDOF(void) const
{
    return modNumDOF;
}


int
TransformationDOF_Group::getNumFreeDOF(void) const
{
    if (modID != 0) {
	int numFreeDOF = modNumDOF;
	for (int i=0; i<modNumDOF; i++)
	    if ((*modID)(i) == -1)
		numFreeDOF--;
	return numFreeDOF;	
    } else
	return this->DOF_Group::getNumFreeDOF();
}

int
TransformationDOF_Group::getNumConstrainedDOF(void) const
{   
    if (modID != 0) {    
	int numConstr = 0;
	for (int i=0; i<modNumDOF; i++)
	    if ((*modID)(i) < 0)
		numConstr++;    
	
	return numConstr;
    } else
	return this->DOF_Group::getNumConstrainedDOF();	
}    



const Matrix &
TransformationDOF_Group::getTangent(Integrator *theIntegrator) 
{


    const Matrix &unmodTangent = this->DOF_Group::getTangent(theIntegrator);
    Matrix *T = this->getT();
    if (T != 0) {
	// *modTangent = (*T) ^ unmodTangent * (*T);
	modTangent->addMatrixTripleProduct(0.0, *T, unmodTangent, 1.0);
	return *modTangent;
	
    } else 
      return unmodTangent;
}

const Vector &
TransformationDOF_Group::getUnbalance(Integrator *theIntegrator)
{
    const Vector &unmodUnbalance = 
	this->DOF_Group::getUnbalance(theIntegrator);

    Matrix *T = this->getT();
    if (T != 0) {
	// *modUnbalance = (*T) ^ unmodUnbalance;
	modUnbalance->addMatrixTransposeVector(0.0, *T, unmodUnbalance, 1.0);
	return *modUnbalance;    
    } else
	return unmodUnbalance;
}


const Vector & 
TransformationDOF_Group::getCommittedDisp(void)
{
    const Vector &responseC = myNode->getDisp();
    
    if (theMP == 0)
	return responseC;
    else {
	int retainedNode = theMP->getNodeRetained();
	Domain *theDomain = myNode->getDomain();
	Node *retainedNodePtr = theDomain->getNode(retainedNode);
	const Vector &responseR = retainedNodePtr->getDisp();
	const ID &retainedDOF = theMP->getRetainedDOFs();
	const ID &constrainedDOF = theMP->getConstrainedDOFs();    	
	int numCNodeDOF = myNode->getNumberDOF();
	int numRetainedNodeDOF = retainedDOF.Size();

	int loc = 0;
	for (int i=0; i<numCNodeDOF; i++) {
	    if (constrainedDOF.getLocation(i) < 0) {
		(*modUnbalance)(loc) = responseC(i);
		loc++;
	    } 
	}
	for (int j=0; j<numRetainedNodeDOF; j++) {
	    int dof = retainedDOF(j);
	    (*modUnbalance)(loc) = responseR(dof);
	    loc++;
	}

	return *modUnbalance;
    }
}

const Vector & 
TransformationDOF_Group::getCommittedVel(void)
{
    const Vector &responseC = myNode->getVel();
    
    if (theMP == 0)
	return responseC;
    else {
	int retainedNode = theMP->getNodeRetained();
	Domain *theDomain = myNode->getDomain();
	Node *retainedNodePtr = theDomain->getNode(retainedNode);
	const Vector &responseR = retainedNodePtr->getVel();
	const ID &retainedDOF = theMP->getRetainedDOFs();
	const ID &constrainedDOF = theMP->getConstrainedDOFs();    	
	int numCNodeDOF = myNode->getNumberDOF();
	int numRetainedNodeDOF = retainedDOF.Size();

	int loc = 0;
	for (int i=0; i<numCNodeDOF; i++) {
	    if (constrainedDOF.getLocation(i) < 0) {
		(*modUnbalance)(loc) = responseC(i);
		loc++;
	    } 
	}
	for (int j=0; j<numRetainedNodeDOF; j++) {
	    int dof = retainedDOF(j);
	    (*modUnbalance)(loc) = responseR(dof);
	    loc++;
	}
	
	return *modUnbalance;
    }
}


const Vector & 
TransformationDOF_Group::getCommittedAccel(void)
{
    const Vector &responseC = myNode->getAccel();
    
    if (theMP == 0)
	return responseC;
    else {
	int retainedNode = theMP->getNodeRetained();
	Domain *theDomain = myNode->getDomain();
	Node *retainedNodePtr = theDomain->getNode(retainedNode);
	const Vector &responseR = retainedNodePtr->getAccel();
	const ID &retainedDOF = theMP->getRetainedDOFs();
	const ID &constrainedDOF = theMP->getConstrainedDOFs();    	
	int numCNodeDOF = myNode->getNumberDOF();
	int numRetainedNodeDOF = retainedDOF.Size();

	int loc = 0;
	for (int i=0; i<numCNodeDOF; i++) {
	    if (constrainedDOF.getLocation(i) < 0) {
		(*modUnbalance)(loc) = responseC(i);
		loc++;
	    } 
	}
	for (int j=0; j<numRetainedNodeDOF; j++) {
	    int dof = retainedDOF(j);
	    (*modUnbalance)(loc) = responseR(dof);
	    loc++;
	}
	
	return *modUnbalance;
    }
}

// void setNodeDisp(const Vector &u);
//	Method to set the corresponding nodes displacements to the
//	values in u, components identified by myID;

void
TransformationDOF_Group::setNodeDisp(const Vector &u)
{
  // call base class method and return if no MP_Constraint
  if (theMP == 0) {
    this->DOF_Group::setNodeDisp(u);
    return;
  }

  const ID &theID = this->getID();
  for (int i=0; i<modNumDOF; i++) {
    int loc = theID(i);
    if (loc >= 0)
      (*modUnbalance)(i) = u(loc);
    else 
      (*modUnbalance)(i) = 0.0;
  }


  if (needRetainedData == 0) {
    int retainedNode = theMP->getNodeRetained();
    Domain *theDomain = myNode->getDomain();
    Node *retainedNodePtr = theDomain->getNode(retainedNode);
    const Vector &responseR = retainedNodePtr->getTrialDisp();
    const ID &retainedDOF = theMP->getRetainedDOFs();

    for (int i=numConstrainedNodeRetainedDOF, j=0; i<modNumDOF; i++, j++) {
      int loc = theID(i);
      if (loc < 0)
	(*modUnbalance)(i) = responseR(retainedDOF(j));
    }
  }

  Matrix *T = this->getT();
  // *unbalance = (*T) * (*modUnbalance);
  unbalance->addMatrixVector(0.0, *T, *modUnbalance, 1.0);

  const Vector &disp = myNode->getTrialDisp();

  int numDOF = myNode->getNumberDOF();
  for (int i=0; i<numDOF; i++) {
    if (theSPs[i] != 0)
      (*unbalance)(i) = disp(i);
  }
  myNode->setTrialDisp(*unbalance);
}

void
TransformationDOF_Group::setNodeVel(const Vector &u)
{
    // call base class method and return if no MP_Constraint
    if (theMP == 0) {
	this->DOF_Group::setNodeVel(u);
	return;
    }
    
   const ID &theID = this->getID();
   for (int i=0; i<modNumDOF; i++) {
     int loc = theID(i);
     if (loc >= 0)
       (*modUnbalance)(i) = u(loc);
     else   
       (*modUnbalance)(i) = 0.0;	    
   }    

  if (needRetainedData == 0) {
    int retainedNode = theMP->getNodeRetained();
    Domain *theDomain = myNode->getDomain();
    Node *retainedNodePtr = theDomain->getNode(retainedNode);
    const Vector &responseR = retainedNodePtr->getTrialVel();
    const ID &retainedDOF = theMP->getRetainedDOFs();
    
    for (int i=numConstrainedNodeRetainedDOF, j=0; i<modNumDOF; i++, j++) {
      int loc = theID(i);
      if (loc < 0)
	(*modUnbalance)(i) = responseR(retainedDOF(j));
    }
  }

  Matrix *T = this->getT();
  // *unbalance = (*T) * (*modUnbalance);
  unbalance->addMatrixVector(0.0, *T, *modUnbalance, 1.0);

  const Vector &vel = myNode->getTrialVel();
  int numDOF = myNode->getNumberDOF();
  for (int i=0; i<numDOF; i++) {
    if (theSPs[i] != 0)
      (*unbalance)(i) = vel(i);
  }
  myNode->setTrialVel(*unbalance);
}


void
TransformationDOF_Group::setNodeAccel(const Vector &u)
{
    // call base class method and return if no MP_Constraint
    if (theMP == 0) {
	this->DOF_Group::setNodeAccel(u);
	return;
    }
    
   const ID &theID = this->getID();
   for (int i=0; i<modNumDOF; i++) {
	int loc = theID(i);
	if (loc >= 0)
	    (*modUnbalance)(i) = u(loc);
	else 	// NO SP STUFF .. WHAT TO DO
	    (*modUnbalance)(i) = 0.0;	    

    }    

  if (needRetainedData == 0) {
    int retainedNode = theMP->getNodeRetained();
    Domain *theDomain = myNode->getDomain();
    Node *retainedNodePtr = theDomain->getNode(retainedNode);
    const Vector &responseR = retainedNodePtr->getTrialAccel();
    const ID &retainedDOF = theMP->getRetainedDOFs();
    
    for (int i=numConstrainedNodeRetainedDOF, j=0; i<modNumDOF; i++, j++) {
      int loc = theID(i);
      if (loc < 0)
	(*modUnbalance)(i) = responseR(retainedDOF(j));
    }
  }

    Matrix *T = this->getT();
    // *unbalance = (*T) * (*modUnbalance);
    unbalance->addMatrixVector(0.0, *T, *modUnbalance, 1.0);
    const Vector &accel = myNode->getTrialAccel();
    int numDOF = myNode->getNumberDOF();
    for (int i=0; i<numDOF; i++) {
      if (theSPs[i] != 0)
	(*unbalance)(i) = accel(i);
    }
    myNode->setTrialAccel(*unbalance);
}


// void setNodeIncrDisp(const Vector &u);
//	Method to set the corresponding nodes displacements to the
//	values in u, components identified by myID;

void
TransformationDOF_Group::incrNodeDisp(const Vector &u)
{
    // call base class method and return if no MP_Constraint
    if (theMP == 0) {
	this->DOF_Group::incrNodeDisp(u);
	return;
    }

   const ID &theID = this->getID();

   for (int i=0; i<modNumDOF; i++) {
     int loc = theID(i);
     if (loc >= 0)
       (*modUnbalance)(i) = u(loc);
     else  
       (*modUnbalance)(i) = 0.0;	    
   }    
   
   Matrix *T = this->getT();
   // *unbalance = (*T) * (*modUnbalance);
   unbalance->addMatrixVector(0.0, *T, *modUnbalance, 1.0);
   
   int numDOF = myNode->getNumberDOF();
   for (int i=0; i<numDOF; i++) {
     if (theSPs[i] != 0)
       (*unbalance)(i) = 0.0;
   }
   myNode->incrTrialDisp(*unbalance);
}


void
TransformationDOF_Group::incrNodeVel(const Vector &u)
{
  // call base class method and return if no MP_Constraint
  if (theMP == 0) {
    this->DOF_Group::incrNodeVel(u);
    return;
  }
  
  const ID &theID = this->getID();
  for (int i=0; i<modNumDOF; i++) {
    int loc = theID(i);
    if (loc >= 0)
      (*modUnbalance)(i) = u(loc);
    else   
      (*modUnbalance)(i) = 0.0;	    
  }    
  Matrix *T = this->getT();
  
  // *unbalance = (*T) * (*modUnbalance);
  unbalance->addMatrixVector(0.0, *T, *modUnbalance, 1.0);
  
  int numDOF = myNode->getNumberDOF();
  for (int i=0; i<numDOF; i++) {
    if (theSPs[i] != 0)
      (*unbalance)(i) = 0.0;
  }
  myNode->incrTrialVel(*unbalance);
}


void
TransformationDOF_Group::incrNodeAccel(const Vector &u)
{
  // call base class method and return if no MP_Constraint
  if (theMP == 0) {
    this->DOF_Group::incrNodeAccel(u);
    return;
  }
  
  const ID &theID = this->getID();
  for (int i=0; i<modNumDOF; i++) {
    int loc = theID(i);
    if (loc >= 0)
      (*modUnbalance)(i) = u(loc);
    else 	
      (*modUnbalance)(i) = 0.0;	    
  }    
  Matrix *T = this->getT();

  // *unbalance = (*T) * (*modUnbalance);
  unbalance->addMatrixVector(0.0, *T, *modUnbalance, 1.0);
  int numDOF = myNode->getNumberDOF();
  for (int i=0; i<numDOF; i++) {
    if (theSPs[i] != 0)
      (*unbalance)(i) = 0.0;
  }
  myNode->incrTrialAccel(*unbalance);
}


void
TransformationDOF_Group::setEigenvector(int mode, const Vector &u)
{
  // call base class method and return if no MP_Constraint
  if (theMP == 0) {
    this->DOF_Group::setEigenvector(mode, u);
    return;
  }
  
  const ID &theID = this->getID();
  for (int i=0; i<modNumDOF; i++) {
    int loc = theID(i);
    if (loc >= 0)
      (*modUnbalance)(i) = u(loc);
    else 	
      (*modUnbalance)(i) = 0.0;	    
  }    
  Matrix *T = this->getT();

    if (T != 0) {
      // *unbalance = (*T) * (*modUnbalance);
      unbalance->addMatrixVector(0.0, *T, *modUnbalance, 1.0);
      myNode->setEigenvector(mode, *unbalance);
    } else
      myNode->setEigenvector(mode, *modUnbalance);
}


Matrix *
TransformationDOF_Group::getT(void)
{
    if (theMP == 0)
	return 0;

    if (theMP->isTimeVarying() == false) {
	return Trans;
    }
    
    int numNodalDOF = myNode->getNumberDOF();
    const ID &retainedDOF = theMP->getRetainedDOFs();
    const ID &constrainedDOF = theMP->getConstrainedDOFs();    
    int numNodalDOFConstrained = constrainedDOF.Size();
    int numRetainedDOF = numNodalDOF - numNodalDOFConstrained;
    int numRetainedNodeDOF = retainedDOF.Size();

    Trans->Zero();
    const Matrix &Ccr = theMP->getConstraint();
    int col = 0;
    for (int i=0; i<numNodalDOF; i++) {
      int loc = constrainedDOF.getLocation(i);
      if (loc < 0) {
	(*Trans)(i,col) = 1.0;
	col++;
      } else {
	for (int j=0; j<numRetainedNodeDOF; j++)
	  (*Trans)(i,j+numRetainedDOF) = Ccr(loc,j);
      }
    }

    return Trans;    
}


int
TransformationDOF_Group::doneID(void)
{
  if (theMP == 0)
    return 0;

    // get number of DOF & verify valid
    int numNodalDOF = myNode->getNumberDOF();
    const ID &retainedDOF = theMP->getRetainedDOFs();
    const ID &constrainedDOF = theMP->getConstrainedDOFs();    
    int numNodalDOFConstrained = constrainedDOF.Size();
    int numRetainedDOF = numNodalDOF - numNodalDOFConstrained;
    int numRetainedNodeDOF = retainedDOF.Size();

    int retainedNode = theMP->getNodeRetained();
    Domain *theDomain = myNode->getDomain();
    Node *retainedNodePtr = theDomain->getNode(retainedNode);
    DOF_Group *retainedGroup = retainedNodePtr->getDOF_GroupPtr();
    const ID &otherID = retainedGroup->getID();
    
    // set the ID for those dof corresponding to dof at another node
    for (int i=0; i<numRetainedNodeDOF; i++) {
	int dof = retainedDOF(i);
	int id = otherID(dof);
	(*modID)(i+numRetainedDOF) = id;
    }

    // if constraint is not time-varying determine the transformation matrix
    if (theMP->isTimeVarying() == false) {
	Trans->Zero();
	const Matrix &Ccr = theMP->getConstraint();
	int col = 0;
	for (int i=0; i<numNodalDOF; i++) {
	    int loc = constrainedDOF.getLocation(i);
	    if (loc < 0) {
		(*Trans)(i,col) = 1.0;
		col++;
	    } else {
		for (int j=0; j<numRetainedNodeDOF; j++)
		    (*Trans)(i,j+numRetainedDOF) = Ccr(loc,j);
	    }
	}
    }
	
    // set the pointers for the tangent and residual
    if (modNumDOF <= MAX_NUM_DOF) {
	// use class wide objects
	if (modVectors[modNumDOF] == 0) {
	    // have to create matrix and vector of size as none yet created
	    modVectors[modNumDOF] = new Vector(modNumDOF);
	    modMatrices[modNumDOF] = new Matrix(modNumDOF,modNumDOF);
	    modUnbalance = modVectors[modNumDOF];
	    modTangent = modMatrices[modNumDOF];
	    if (modUnbalance == 0 || modUnbalance->Size() != modNumDOF ||	
		modTangent == 0 || modTangent->noCols() != modNumDOF)	{  
		opserr << "DOF_Group::DOF_Group(Node *) ";
		opserr << " ran out of memory for vector/Matrix of size :";
		opserr << modNumDOF << endln;
		exit(-1);
	    }
	} else {
	    modUnbalance = modVectors[modNumDOF];
	    modTangent = modMatrices[modNumDOF];
	}
    } else {
	// create matrices and vectors for each object instance
	modUnbalance = new Vector(modNumDOF);
	modTangent = new Matrix(modNumDOF, modNumDOF);
	if (modUnbalance == 0 || modUnbalance->Size() ==0 ||
	    modTangent ==0 || modTangent->noRows() ==0) {
	    
	    opserr << "DOF_Group::DOF_Group(Node *) ";
	    opserr << " ran out of memory for vector/Matrix of size :";
	    opserr << modNumDOF << endln;
	    exit(-1);
	}
    }    

    if (modID != 0) {
      for (int i=numConstrainedNodeRetainedDOF; i<modNumDOF; i++)
	if ((*modID)(i) == -1)
	  needRetainedData = 0;
    }

    return 0;
}

int 
TransformationDOF_Group::addSP_Constraint(SP_Constraint &theSP)
{
    // add the SP_Constraint
    int dof = theSP.getDOF_Number();
    theSPs[dof] = &theSP;

    // set a -1 in the correct ID location
    if (theMP == 0)
	this->setID(dof,-1);
    else {
	const ID &constrainedDOF = theMP->getConstrainedDOFs();
	int loc = 0;
	for (int i=0; i<dof; i++) 
	    if (constrainedDOF.getLocation(i) < 0)
		loc++;
	this->setID(loc,-1);
    } 
    return 0;
}

int 
TransformationDOF_Group::enforceSPs(int doMP)
{
  int numDof = myNode->getNumberDOF();
  
  if (doMP == 1) {
    for (int i=0; i<numDof; i++)
      if (theSPs[i] != 0) {
	double value = theSPs[i]->getValue();
	myNode->setTrialDisp(value, i);
      }
  } 

  else {

    if (needRetainedData == 0) {
      
      if (theMP != 0) {
	
	const ID &theID = this->getID();
	
	int retainedNode = theMP->getNodeRetained();
	Domain *theDomain = myNode->getDomain();
	Node *retainedNodePtr = theDomain->getNode(retainedNode);
	const Vector &responseR = retainedNodePtr->getTrialDisp();
	const ID &retainedDOF = theMP->getRetainedDOFs();
	
	modUnbalance->Zero();    
	for (int i=numConstrainedNodeRetainedDOF, j=0; i<modNumDOF; i++, j++) {
	  int loc = theID(i);
	  if (loc < 0)
	    (*modUnbalance)(i) = responseR(retainedDOF(j));
	}
	
	Matrix *T = this->getT();
	if (T != 0) {
	  
	  // *unbalance = (*T) * (*modUnbalance);
	  unbalance->addMatrixVector(0.0, *T, *modUnbalance, 1.0);
	  
	  const ID &constrainedDOF = theMP->getConstrainedDOFs();
	  for (int i=0; i<constrainedDOF.Size(); i++) {
	    int cDOF = constrainedDOF(i);
	    myNode->setTrialDisp((*unbalance)(cDOF), cDOF);
	  }
	}
      }
    }
  }

  return 0;
}

void
TransformationDOF_Group::addM_Force(const Vector &Udotdot, double fact)
{
    // call base class method and return if no MP_Constraint
    if (theMP == 0 || modID == 0) {
	this->DOF_Group::addM_Force(Udotdot, fact);
	return;
    }
    
   for (int i=0; i<modNumDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modUnbalance)(i) = Udotdot(loc);
	else 	// DO THE SP STUFF
	    (*modUnbalance)(i) = 0.0;	    
    }    

    Vector unmod(Trans->noRows());
    //unmod = (*Trans) * (*modUnbalance);
    unmod.addMatrixVector(0.0, *Trans, *modUnbalance, 1.0);
    this->addLocalM_Force(unmod, fact);
}


const Vector &
TransformationDOF_Group::getM_Force(const Vector &Udotdot, double fact)
{
  // call base class method and return if no MP_Constraint
  if (theMP == 0 || modID == 0) {
    return this->DOF_Group::getM_Force(Udotdot, fact);
  }

  this->DOF_Group::zeroTangent();    
  this->DOF_Group::addMtoTang();    
  const Matrix &unmodTangent = this->DOF_Group::getTangent(0);

    
  Vector data(modNumDOF);
  for (int i=0; i<modNumDOF; i++) {
    int loc = (*modID)(i);
    if (loc >= 0)
      data(i) = Udotdot(loc);
    else 	// DO THE SP STUFF
      data(i) = 0.0;	    
  }    

  Matrix *T = this->getT();
  if (T != 0) {
    // *modTangent = (*T) ^ unmodTangent * (*T);
    modTangent->addMatrixTripleProduct(0.0, *T, unmodTangent, 1.0);
    modUnbalance->addMatrixVector(0.0, *modTangent, data, 1.0);
    
    return *modUnbalance;
  } else {
      modUnbalance->addMatrixVector(0.0, unmodTangent, data, 1.0);
      return *modUnbalance;
  }
}

const Vector &
TransformationDOF_Group::getC_Force(const Vector &Udotdot, double fact)
{
  opserr << "TransformationDOF_Group::getC_Force() - not yet implemented\n";
  return *modUnbalance;
}

const Vector &
TransformationDOF_Group::getTangForce(const Vector &Udotdot, double fact)
{
  opserr << "TransformationDOF_Group::getTangForce() - not yet implemented\n";
  return *modUnbalance;
}



// AddingSensitivity:BEGIN ////////////////////////////////////////
const Vector &
TransformationDOF_Group::getDispSensitivity(int gradNumber)
{
  const Vector &result = this->DOF_Group::getDispSensitivity(gradNumber);

  Matrix *T = this->getT();
  if (T != 0) {
    // *modUnbalance = (*T) ^ unmodUnbalance;
    modUnbalance->addMatrixTransposeVector(0.0, *T, result, 1.0);
    return *modUnbalance;    
  } else
    return result;

}

const Vector &
TransformationDOF_Group::getVelSensitivity(int gradNumber)
{
  const Vector &result = this->DOF_Group::getVelSensitivity(gradNumber);

  Matrix *T = this->getT();
  if (T != 0) {
    // *modUnbalance = (*T) ^ unmodUnbalance;
    modUnbalance->addMatrixTransposeVector(0.0, *T, result, 1.0);
    return *modUnbalance;    
  } else
    return result;
}

const Vector &
TransformationDOF_Group::getAccSensitivity(int gradNumber)
{
  const Vector &result = this->DOF_Group::getAccSensitivity(gradNumber);

  Matrix *T = this->getT();
  if (T != 0) {
    // *modUnbalance = (*T) ^ unmodUnbalance;
    modUnbalance->addMatrixTransposeVector(0.0, *T, result, 1.0);
    return *modUnbalance;    
  } else
    return result;
}
	
	
int
TransformationDOF_Group::saveDispSensitivity(const Vector &u,
					     int gradNum, int numGrads)
{
  // call base class method and return if no MP_Constraint
  if (theMP == 0) {
    return this->DOF_Group::saveDispSensitivity(u, gradNum, numGrads);
  }
  
  const ID &theID = this->getID();
  for (int i=0; i<modNumDOF; i++) {
    int loc = theID(i);
    if (loc >= 0)
      (*modUnbalance)(i) = u(loc);
    // DO THE SP STUFF
  }    
  Matrix *T = this->getT();
  if (T != 0) {
    
    // *unbalance = (*T) * (*modUnbalance);
    unbalance->addMatrixVector(0.0, *T, *modUnbalance, 1.0);
    
  } else
    *unbalance = *modUnbalance;


  myNode->saveDispSensitivity(*unbalance, gradNum, numGrads);
  
  return 0;
}

int
TransformationDOF_Group::saveVelSensitivity(const Vector &u,
					    int gradNum, int numGrads)
{
  // call base class method and return if no MP_Constraint
  if (theMP == 0) {
    return this->DOF_Group::saveVelSensitivity(u, gradNum, numGrads);
  }
  
  const ID &theID = this->getID();
  for (int i=0; i<modNumDOF; i++) {
    int loc = theID(i);
    if (loc >= 0)
      (*modUnbalance)(i) = u(loc);
    // DO THE SP STUFF
  }    
  Matrix *T = this->getT();
  if (T != 0) {
    
    // *unbalance = (*T) * (*modUnbalance);
    unbalance->addMatrixVector(0.0, *T, *modUnbalance, 1.0);
    
  } else
    *unbalance = *modUnbalance;


  myNode->saveVelSensitivity(*unbalance, gradNum, numGrads);
  
  return 0;
}

int
TransformationDOF_Group::saveAccSensitivity(const Vector &u,
					    int gradNum, int numGrads)
{
  // call base class method and return if no MP_Constraint
  if (theMP == 0) {
    return this->DOF_Group::saveAccSensitivity(u, gradNum, numGrads);
  }
  
  const ID &theID = this->getID();
  for (int i=0; i<modNumDOF; i++) {
    int loc = theID(i);
    if (loc >= 0)
      (*modUnbalance)(i) = u(loc);
    // DO THE SP STUFF
  }    
  Matrix *T = this->getT();
  if (T != 0) {
    
    // *unbalance = (*T) * (*modUnbalance);
    unbalance->addMatrixVector(0.0, *T, *modUnbalance, 1.0);
    
  } else
    *unbalance = *modUnbalance;


  myNode->saveAccelSensitivity(*unbalance, gradNum, numGrads);
  
  return 0;
}

int 
TransformationDOF_Group::saveSensitivity(const Vector &u,
					 const Vector &udot,
					 const Vector &udotdot,
					 int gradNum, int numGrads)
{
  // call base class method and return if no MP_Constraint
  if (theMP == 0) {
    return this->DOF_Group::saveSensitivity(u, udot, udotdot, gradNum, numGrads);
  }
  
  int ok =0;

  ok += this->saveDispSensitivity(u, gradNum, numGrads);
  ok += this->saveVelSensitivity(udot, gradNum, numGrads);
  ok += this->saveAccSensitivity(udotdot, gradNum, numGrads);

  return ok;
}

void  
TransformationDOF_Group::addM_ForceSensitivity(const Vector &Udotdot, double fact)
{
    // call base class method and return if no MP_Constraint
    if (theMP == 0 || modID == 0) {
	this->DOF_Group::addM_ForceSensitivity(Udotdot, fact);
	return;
    }
    
   for (int i=0; i<modNumDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modUnbalance)(i) = Udotdot(loc);
	else 	// DO THE SP STUFF
	    (*modUnbalance)(i) = 0.0;	    
    }    

    Vector unmod(Trans->noRows());
    //unmod = (*Trans) * (*modUnbalance);
    unmod.addMatrixVector(0.0, *Trans, *modUnbalance, 1.0);
    this->DOF_Group::addM_ForceSensitivity(unmod, fact);
}

void
TransformationDOF_Group::addD_Force(const Vector &Udot, double fact)
{
    // call base class method and return if no MP_Constraint
    if (theMP == 0 || modID == 0) {
	this->DOF_Group::addD_Force(Udot, fact);
	return;
    }
    
   for (int i=0; i<modNumDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modUnbalance)(i) = Udot(loc);
	else 	// DO THE SP STUFF
	    (*modUnbalance)(i) = 0.0;	    
    }    

    Vector unmod(Trans->noRows());
    //unmod = (*Trans) * (*modUnbalance);
    unmod.addMatrixVector(0.0, *Trans, *modUnbalance, 1.0);
    this->DOF_Group::addD_Force(unmod, fact);
}

void
TransformationDOF_Group::addD_ForceSensitivity(const Vector &Udot, double fact)
{
    // call base class method and return if no MP_Constraint
    if (theMP == 0 || modID == 0) {
	this->DOF_Group::addD_ForceSensitivity(Udot, fact);
	return;
    }
    
   for (int i=0; i<modNumDOF; i++) {
	int loc = (*modID)(i);
	if (loc >= 0)
	    (*modUnbalance)(i) = Udot(loc);
	else 	// DO THE SP STUFF
	    (*modUnbalance)(i) = 0.0;	    
    }    

    Vector unmod(Trans->noRows());
    //unmod = (*Trans) * (*modUnbalance);
    unmod.addMatrixVector(0.0, *Trans, *modUnbalance, 1.0);
    this->DOF_Group::addD_ForceSensitivity(unmod, fact);
}

// AddingSensitivity:END //////////////////////////////////////////
