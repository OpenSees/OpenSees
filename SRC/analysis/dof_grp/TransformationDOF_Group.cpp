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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/dof_grp/TransformationDOF_Group.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/dof_grp/TransformationDOF_Group.C
//
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

#define MAX_NUM_DOF 16


// static variables initialisation
Matrix **TransformationDOF_Group::modMatrices; 
Vector **TransformationDOF_Group::modVectors;  
int TransformationDOF_Group::numTransDOFs(0);     // number of objects

TransformationDOF_Group::TransformationDOF_Group(int tag, Node *node, 
						 MP_Constraint *mp)  
:DOF_Group(tag,node),
 theMP(mp),Trans(0),modTangent(0),modUnbalance(0),modID(0),theSPs(0)
{
    // determine the number of DOF 
    int numNodalDOF = node->getNumberDOF();
    const ID &retainedDOF = mp->getRetainedDOFs();
    const ID &constrainedDOF = mp->getConstrainedDOFs();    
    int numNodalDOFConstrained = constrainedDOF.Size();
    int numConstrainedNodeRetainedDOF = numNodalDOF - numNodalDOFConstrained;
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
	
	cerr << "FATAL TransformationDOF_Group::TransformationDOF_Group() -";
	cerr << " ran out of memory for size: " << modNumDOF << endl;
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
	    cerr << "TransformationDOF_Group::TransformationDOF_Group(Node *) ";
	    cerr << " ran out of memory";	    
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
		cerr << "TransformationDOF_Group::TransformationDOF_Group(Node *) ";
		cerr << " ran out of memory for vector/Matrix of size :";
		cerr << modNumDOF << endl;
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
	    
	    cerr << "TransformationDOF_Group::TransformationDOF_Group(Node *) ";
	    cerr << " ran out of memory for vector/Matrix of size :";
	    cerr << modNumDOF << endl;
	    exit(-1);
	}
    }
    
    numTransDOFs++;
}

void 
TransformationDOF_Group::setID(int dof, int value)
{
    if (theMP == 0)
	this->DOF_Group::setID(dof,value);
    else
	(*modID)(dof) = value;
}
	

TransformationDOF_Group::TransformationDOF_Group(int tag, Node *node)
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
	    cerr << "TransformationDOF_Group::TransformationDOF_Group(Node *) ";
	    cerr << " ran out of memory";	    
	}
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    modMatrices[i] = 0;
	    modVectors[i] = 0;
	}
    }    

    numTransDOFs++;
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
	    if ((*modID)(i) < 0)
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
	*modTangent = (*T) ^ unmodTangent * (*T);
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
	*modUnbalance = (*T) ^ unmodUnbalance;
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
	// DO THE SP STUFF
    }    
    Matrix *T = this->getT();
    if (T != 0) {
	*unbalance = (*T) * (*modUnbalance);
	myNode->setTrialDisp(*unbalance);
    } else
	myNode->setTrialDisp(*modUnbalance);
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
	else 	// NO SP STUFF .. WHAT TO DO
	    (*modUnbalance)(i) = 0.0;	    

    }    
    Matrix *T = this->getT();
    if (T != 0) {
	*unbalance = (*T) * (*modUnbalance);
	myNode->setTrialVel(*unbalance);
    } else
	myNode->setTrialVel(*modUnbalance);
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
    Matrix *T = this->getT();
    if (T != 0) {
	*unbalance = (*T) * (*modUnbalance);
	myNode->setTrialAccel(*unbalance);
    } else
	myNode->setTrialAccel(*modUnbalance);
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
	else 	// DO THE SP STUFF
	    (*modUnbalance)(i) = 0.0;	    
    }    

    Matrix *T = this->getT();
    if (T != 0) {
	*unbalance = (*T) * (*modUnbalance);
	myNode->incrTrialDisp(*unbalance);
    } else 
	myNode->incrTrialDisp(*modUnbalance);
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
	else 	// DO THE SP STUFF
	    (*modUnbalance)(i) = 0.0;	    
    }    
    Matrix *T = this->getT();
    if (T != 0) {
	*unbalance = (*T) * (*modUnbalance);
	myNode->incrTrialVel(*unbalance);
    } else
	myNode->incrTrialVel(*modUnbalance);
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
	else 	// DO THE SP STUFF
	    (*modUnbalance)(i) = 0.0;	    
    }    
    Matrix *T = this->getT();
    if (T != 0) {
	*unbalance = (*T) * (*modUnbalance);
	myNode->incrTrialAccel(*unbalance);
    } else
	myNode->incrTrialAccel(*modUnbalance);
}


Matrix *
TransformationDOF_Group::getT(void)
{
    if (theMP == 0)
	return 0;

    if (theMP->isTimeVarying() == false) {
	return Trans;
    }
    
    int numNodalDOF = this->DOF_Group::getNumDOF();
    // const ID &retainedDOF = theMP->getRetainedDOFs();
    const ID &constrainedDOF = theMP->getConstrainedDOFs();    
    int numNodalDOFConstrained = constrainedDOF.Size();
    int numRetainedDOF = numNodalDOF - numNodalDOFConstrained;
    //    int numReatainedNodeDOF = retainedDOF.Size();

    Trans->Zero();
    const Matrix &Ccr = theMP->getConstraint();
    int col = 0;
    for (int i=0; i<numNodalDOF; i++) {
	int loc = constrainedDOF.getLocation(i);
	if (loc < 0) {
	    (*Trans)(i,col) = 1.0;
	    col++;
	} else {
	    for (int j=0; j<numRetainedDOF; j++)
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
    int numNodalDOF = this->getNumDOF();
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
		cerr << "DOF_Group::DOF_Group(Node *) ";
		cerr << " ran out of memory for vector/Matrix of size :";
		cerr << modNumDOF << endl;
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
	    
	    cerr << "DOF_Group::DOF_Group(Node *) ";
	    cerr << " ran out of memory for vector/Matrix of size :";
	    cerr << modNumDOF << endl;
	    exit(-1);
	}
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
TransformationDOF_Group::enforceSPs(void)
{
    Vector trialDisp(myNode->getTrialDisp());
    int numDof = myNode->getNumberDOF();
    for (int i=0; i<numDof; i++)
	if (theSPs[i] != 0) {
	    double value = theSPs[i]->getValue();
	    trialDisp(i) = value;
	}
    myNode->setTrialDisp(trialDisp);
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
    unmod = (*Trans) * (*modUnbalance);
    this->addLocalM_Force(unmod,fact);
}
