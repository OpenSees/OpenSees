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
                                                                        
// $Revision: 1.2 $
// $Date: 2000-12-13 04:25:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/dof_grp/DOF_Group.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/dof_grp/DOF_Group.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the code for implementing the methods
// of the DOF_Group class interface.
//
// What: "@(#) DOF_Group.C, revA"


#include <DOF_Group.h>
#include <stdlib.h>

#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <TransientIntegrator.h>

#define MAX_NUM_DOF 256


// static variables initialisation
Matrix DOF_Group::errMatrix(1,1);
Vector DOF_Group::errVect(1);
Matrix **DOF_Group::theMatrices; // array of pointers to class wide matrices
Vector **DOF_Group::theVectors;  // array of pointers to class widde vectors
int DOF_Group::numDOFs(0);           // number of objects


//  DOF_Group(Node *);
//	construictor that take the corresponding model node.

DOF_Group::DOF_Group(int tag, Node *node)
:unbalance(0), tangent(0), myNode(node), 
 myTag(tag), myID(node->getNumberDOF()), 
 numDOF(node->getNumberDOF())
{
    // get number of DOF & verify valid
    int numDOF = node->getNumberDOF();
    if (numDOF <= 0) {
	cerr << "DOF_Group::DOF_Group(Node *) ";
	cerr << " node must have at least 1 dof " << *node;
	exit(-1);
    }	

    // check the ID created is of appropriate size
    if (myID.Size() != numDOF) {
	cerr << "DOF_Group::DOF_Group(Node *) ";
	cerr << " ran out of memory creating ID for node " << *node;
	exit(-1);
    }

    // initially set all the IDs to be -2
    for (int i=0; i<numDOF; i++)
	myID(i) = -2;
    
    // if this is the first DOF_Group we now
    // create the arrays used to store pointers to class wide
    // matrix and vector objects used to return tangent and residual
    if (numDOFs == 0) {
	theMatrices = new Matrix *[MAX_NUM_DOF+1];
	theVectors  = new Vector *[MAX_NUM_DOF+1];
	
	if (theMatrices == 0 || theVectors == 0) {
	    cerr << "DOF_Group::DOF_Group(Node *) ";
	    cerr << " ran out of memory";	    
	}
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    theMatrices[i] = 0;
	    theVectors[i] = 0;
	}
    }    
    
    // set the pointers for the tangent and residual
    if (numDOF <= MAX_NUM_DOF) {
	// use class wide objects
	if (theVectors[numDOF] == 0) {
	    // have to create matrix and vector of size as none yet created
	    theVectors[numDOF] = new Vector(numDOF);
	    theMatrices[numDOF] = new Matrix(numDOF,numDOF);
	    unbalance = theVectors[numDOF];
	    tangent = theMatrices[numDOF];
	    if (unbalance == 0 || unbalance->Size() != numDOF ||	
		tangent == 0 || tangent->noCols() != numDOF)	{  
		cerr << "DOF_Group::DOF_Group(Node *) ";
		cerr << " ran out of memory for vector/Matrix of size :";
		cerr << numDOF << endl;
		exit(-1);
	    }
	} else {
	    unbalance = theVectors[numDOF];
	    tangent = theMatrices[numDOF];
	}
    } else {
	// create matrices and vectors for each object instance
	unbalance = new Vector(numDOF);
	tangent = new Matrix(numDOF, numDOF);
	if (unbalance == 0 || unbalance->Size() ==0 ||
	    tangent ==0 || tangent->noRows() ==0) {
	    
	    cerr << "DOF_Group::DOF_Group(Node *) ";
	    cerr << " ran out of memory for vector/Matrix of size :";
	    cerr << numDOF << endl;
	    exit(-1);
	}
    }
    
    numDOFs++;
}


DOF_Group::DOF_Group(int tag, int ndof)
:unbalance(0), tangent(0), myNode(0), 
 myTag(tag), myID(ndof), 
 numDOF(ndof)
{
    // get number of DOF & verify valid
    int numDOF = ndof;
    if (numDOF <= 0) {
	cerr << "DOF_Group::DOF_Group(int, int ndof) ";
	cerr << ndof << " ndof specified, there must be at least 1\n";
	exit(-1);
    }	

    // check the ID created is of appropriate size
    if (myID.Size() != numDOF) {
	cerr << "DOF_Group::DOF_Group(int, int ndof) ";
	cerr << " ran out of memory creating ID of size " << ndof << endl;
	exit(-1);
    }

    // initially set all the IDs to be -2
    for (int i=0; i<numDOF; i++)
	myID(i) = -2;
    
    // if this is the first DOF_Group we now
    // create the arrays used to store pointers to class wide
    // matrix and vector objects used to return tangent and residual
    if (numDOFs == 0) {
	theMatrices = new Matrix *[MAX_NUM_DOF+1];
	theVectors  = new Vector *[MAX_NUM_DOF+1];
	
	if (theMatrices == 0 || theVectors == 0) {
	    cerr << "DOF_Group::DOF_Group(int, int ndof) ";
	    cerr << " ran out of memory";	    
	}
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    theMatrices[i] = 0;
	    theVectors[i] = 0;
	}
    }    

    // set the pointers for the tangent and residual
    if (numDOF <= MAX_NUM_DOF) {
	// use class wide objects
	if (theVectors[numDOF] == 0) {
	    // have to create matrix and vector of size as none yet created
	    theVectors[numDOF] = new Vector(numDOF);
	    theMatrices[numDOF] = new Matrix(numDOF,numDOF);
	    unbalance = theVectors[numDOF];
	    tangent = theMatrices[numDOF];
	    if (unbalance == 0 || unbalance->Size() != numDOF ||	
		tangent == 0 || tangent->noCols() != numDOF)	{  
		cerr << "DOF_Group::DOF_Group(int, int ndof) ";
		cerr << " ran out of memory for vector/Matrix of size :";
		cerr << numDOF << endl;
		exit(-1);
	    }
	} else {
	    unbalance = theVectors[numDOF];
	    tangent = theMatrices[numDOF];
	}
    } else {
	// create matrices and vectors for each object instance
	unbalance = new Vector(numDOF);
	tangent = new Matrix(numDOF, numDOF);
	if (unbalance == 0 || tangent ==0 ||
	    tangent ==0 || tangent->noRows() ==0) {
	    
	    cerr << "DOF_Group::DOF_Group(int, int ndof) ";
	    cerr << " ran out of memory for vector/Matrix of size :";
	    cerr << numDOF << endl;
	    exit(-1);
	}
    }
    
    numDOFs++;
}

// ~DOF_Group();    
//	destructor.

DOF_Group::~DOF_Group()
{
    numDOFs--;

    // set the pointer in the associated Node to 0, to stop
    // segmentation fault if node tries to use this object after destroyed
    if (myNode != 0) {
	myNode->setDOF_GroupPtr(0);
    }
    
    int numDOF = unbalance->Size();
    
    // delete tangent and residual if created specially
    if (numDOF > MAX_NUM_DOF) {
	if (tangent != 0) delete tangent;
	if (unbalance != 0) delete unbalance;
    }

    // if this is the last FE_Element, clean up the
    // storage for the matrix and vector objects
    if (numDOFs == 0) {
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    if (theVectors[i] != 0)
		delete theVectors[i];
	    if (theMatrices[i] != 0)
		delete theMatrices[i];
	}	
	delete [] theMatrices;
	delete [] theVectors;
    }    
}    

// void setID(int index, int value);
//	Method to set the corresponding index of the ID to value.

void
DOF_Group::setID(int index, int value)
{
    if ((index >= 0) && (index < numDOF))
	myID(index) = value;
    else {
	cerr << "WARNING DOF_Group::setID - invalid location ";
	cerr << index << " in ID of size " << numDOF << endl;
    }	
}

// void setID(const ID &);
//	Method to set the ID to be same as that given.

void
DOF_Group::setID(const ID &copy)
{
    myID = copy;
}
 

// const ID &getID(void) const;
//	Method to return the current ID.

const ID &
DOF_Group::getID(void) const
{
    return myID;
}



int
DOF_Group::getNumDOF(void) const
{
    return numDOF;
}


int
DOF_Group::getNumFreeDOF(void) const
{
    int numFreeDOF = numDOF;
    for (int i=0; i<numDOF; i++)
	if (myID(i) < 0)
	    numFreeDOF--;
    
    return numFreeDOF;
}

int
DOF_Group::getNumConstrainedDOF(void) const
{   
    int numConstr = 0;
    for (int i=0; i<numDOF; i++)
	if (myID(i) < 0)
	    numConstr++;    

    return numConstr;
}    



const Matrix &
DOF_Group::getTangent(Integrator *theIntegrator) 
{	
    if (theIntegrator != 0)
	theIntegrator->formNodTangent(this);    
    return *tangent;
}

void  
DOF_Group::zeroTangent(void)
{
    tangent->Zero();
}


void  
DOF_Group::addMtoTang(double fact)
{
    if (myNode != 0) {
	if (tangent->addMatrix(1.0, myNode->getMass(), fact) < 0) {
	    cerr << "DOF_Group::addMtoTang(void) ";
	    cerr << " invoking addMatrix() on the tangent failed\n";	    
	}
    }
    else {
	cerr << "DOF_Group::addMtoTang(void) - no Node associated";
	cerr << " subclass should provide the method \n";	    
    }	
}



void
DOF_Group::zeroUnbalance(void) 
{
    unbalance->Zero();
}


const Vector &
DOF_Group::getUnbalance(Integrator *theIntegrator)
{
    if (theIntegrator != 0)
	theIntegrator->formNodUnbalance(this);

    return *unbalance;
}


void
DOF_Group::addPtoUnbalance(double fact)
{
    if (myNode != 0) {
	if (unbalance->addVector(1.0, myNode->getUnbalancedLoad(), fact) < 0) {
	    cerr << "DOF_Group::addPIncInertiaToUnbalance() -";
	    cerr << " invoking addVector() on the unbalance failed\n";	    
	}
    }
    else {
	cerr << "DOF_Group::addPtoUnbalance() - no Node associated";
	cerr << " subclass should provide the method \n";	    
    }		
}


void
DOF_Group::addPIncInertiaToUnbalance(double fact)
{
    if (myNode != 0) {
	if (unbalance->addVector(1.0, myNode->getUnbalancedLoadIncInertia(), 
				 fact) < 0) {

	    cerr << "DOF_Group::addPIncInertiaToUnbalance() - ";
	    cerr << " invoking addVector() on the unbalance failed\n";	    
	}
    }
     else {
	cerr << "DOF_Group::addPIncInertiaToUnbalance() - no Node associated";
	cerr << " subclass should provide the method \n";	    
    }		
}


void  
DOF_Group::addM_Force(const Vector &Udotdot, double fact)
{
    if (myNode == 0) {
	cerr << "DOF_Group::addM_Force() - no Node associated";	
	cerr << " subclass should not call this method \n";	    
	return;
    }

    Vector accel(numDOF);
    // get accel for the unconstrained dof
    for (int i=0; i<numDOF; i++) {
	int loc = myID(i);
	if (loc >= 0)
	    accel(i) = Udotdot(loc); 
	else accel(i) = 0.0;
    }
	
    if (unbalance->addMatrixVector(1.0, myNode->getMass(), accel, fact) < 0) {  
	cerr << "DOF_Group::addM_Force() ";
	cerr << " invoking addMatrixVector() on the unbalance failed\n";
    }
    else {

    }		
}


const Vector & 
DOF_Group::getCommittedDisp(void)
{
    if (myNode == 0) {
	cerr << "DOF_Group::getCommittedDisp: no associated Node ";
	cerr << " returning the error Vector\n";
	return errVect;
    }
    return myNode->getDisp();
}


const Vector & 
DOF_Group::getCommittedVel(void)
{
    if (myNode == 0) {
	cerr << "DOF_Group::getCommittedVel: no associated Node ";
	cerr << " returning the error Vector\n";
	return errVect;	
    }
    return myNode->getVel();
}


const Vector & 
DOF_Group::getCommittedAccel(void)
{
    if (myNode == 0) {
	cerr << "DOF_Group::getCommittedAccel: no associated Node ";
	cerr << " returning the error Vector\n";
	return errVect;	
    }
    return myNode->getAccel();
}

// void setNodeDisp(const Vector &u);
//	Method to set the corresponding nodes displacements to the
//	values in u, components identified by myID;

void
DOF_Group::setNodeDisp(const Vector &u)
{
    if (myNode == 0) {
	cerr << "DOF_Group::setNodeDisp: no associated Node\n";
	return;
    }
    
    Vector &disp = *unbalance;
    disp = myNode->getTrialDisp();
    int i;
    
    // get disp for my dof out of vector u
    for (i=0; i<numDOF; i++) {
	int loc = myID(i);
	if (loc >= 0)
	    disp(i) = u(loc);  
    }

    myNode->setTrialDisp(disp);
}
	
	
// void setNodeVel(const Vector &udot);
//	Method to set the corresponding nodes velocities to the
//	values in udot, components identified by myID;

void
DOF_Group::setNodeVel(const Vector &udot)
{

    if (myNode == 0) {
	cerr << "DOF_Group::setNodeVel: 0 Node Pointer\n";
	return;
    }
    
    Vector &vel = *unbalance;
    vel = myNode->getTrialVel();
    int i;
    
    // get vel for my dof out of vector udot
    for (i=0; i<numDOF; i++) {
	int loc = myID(i);	    	
	if (loc >= 0) 
	    vel(i) = udot(loc);  
    }

    myNode->setTrialVel(vel);

}



// void setNodeAccel(const Vector &udotdot);
//	Method to set the corresponding nodes accelerations to the
//	values in udotdot, components identified by myID;

void
DOF_Group::setNodeAccel(const Vector &udotdot)
{

    if (myNode == 0) {
	cerr << "DOF_Group::setNodeAccel: 0 Node Pointer\n";
	return;
    }

    Vector &accel = *unbalance;;
    accel = myNode->getTrialAccel();
    int i;
    
    // get disp for the unconstrained dof
    for (i=0; i<numDOF; i++) {
	int loc = myID(i);	    		
	if (loc >= 0)
	    accel(i) = udotdot(loc); 
    }

    myNode->setTrialAccel(accel);
}


// void setNodeIncrDisp(const Vector &u);
//	Method to set the corresponding nodes displacements to the
//	values in u, components identified by myID;

void
DOF_Group::incrNodeDisp(const Vector &u)
{
    if (myNode == 0) {
	cerr << "DOF_Group::setNodeDisp: 0 Node Pointer\n";
	exit(-1);
    }

    Vector &disp = *unbalance;;

    if (disp.Size() == 0) {
      cerr << "DOF_Group::setNodeIncrDisp - out of space\n";
      return;
    }
    int i;

    // get disp for my dof out of vector u
    for (i=0; i<numDOF; i++) {
	int loc = myID(i);	    			
	if (loc >= 0)
	    disp(i) = u(loc);  
	else disp(i) = 0.0;  
    }

    myNode->incrTrialDisp(disp);

}
	
	
// void setNodeincrVel(const Vector &udot);
//	Method to set the corresponding nodes velocities to the
//	values in udot, components identified by myID;

void
DOF_Group::incrNodeVel(const Vector &udot)
{

    if (myNode == 0) {
	cerr << "DOF_Group::setNodeVel: 0 Node Pointer\n";
	exit(-1);
    }
    
    Vector &vel = *unbalance;
    int i;
    
    // get vel for my dof out of vector udot
    for (i=0; i<numDOF; i++) {
	int loc = myID(i);
	if (loc >= 0)
	    vel(i) = udot(loc);  // -1 for dof labelled 1 through ndof
	else  vel(i) = 0.0;
    }
    myNode->incrTrialVel(vel);
}



// void setNodeIncrAccel(const Vector &udotdot);
//	Method to set the corresponding nodes accelerations to the
//	values in udotdot, components identified by myID;

void
DOF_Group::incrNodeAccel(const Vector &udotdot)
{

    if (myNode == 0) {
	cerr << "DOF_Group::setNodeAccel: 0 Node Pointer\n";
	exit(-1);
    }

    Vector &accel = *unbalance;
    int i;
    
    // get disp for the unconstrained dof
    for (i=0; i<numDOF; i++) {
	int loc = myID(i);
	if (loc >= 0)
	    accel(i) = udotdot(loc); 
	else accel(i) = 0.0;
    }    
    myNode->incrTrialAccel(accel);
}



void
DOF_Group::setEigenvector(int mode, const Vector &theVector)
{

    if (myNode == 0) {
	cerr << "DOF_Group::setNodeAccel: 0 Node Pointer\n";
	exit(-1);
    }

    Vector &eigenvector = *unbalance;
    int i;
    
    // get disp for the unconstrained dof
    for (i=0; i<numDOF; i++) {
	int loc = myID(i);
	if (loc >= 0)
	    eigenvector(i) = theVector(loc); 
	else eigenvector(i) = 0.0;
    }    
    myNode->setEigenvector(mode, eigenvector);
}



int
DOF_Group::getTag(void) const
{
    return myTag;
}

Matrix *
DOF_Group::getT(void)
{
    return 0;
}



void  
DOF_Group::addLocalM_Force(const Vector &accel, double fact)
{
    if (myNode != 0) {
	if (unbalance->addMatrixVector(1.0, myNode->getMass(), accel, fact) < 0) {  
				       
	    cerr << "DOF_Group::addM_Force() ";
	    cerr << " invoking addMatrixVector() on the unbalance failed\n"; 
	}
    }
    else {
	cerr << "DOF_Group::addM_Force() - no Node associated";
	cerr << " subclass should not call this method \n";	    
    }		
}





