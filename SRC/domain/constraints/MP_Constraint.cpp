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
// $Date: 2000-09-15 08:23:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/MP_Constraint.cpp,v $
                                                                        
                                                                        
// File: ~/domain/constraints//MP_Constraint.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementation of class MP_Constraint.
//
// The class MP_Constraint interface:
//

#include <MP_Constraint.h>

#include <stdlib.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
 
// constructor for FEM_ObjectBroker
MP_Constraint::MP_Constraint(int clasTag)
:DomainComponent(0,clasTag),
 nodeRetained(0),nodeConstrained(0),constraint(0),constrDOF(0),retainDOF(0),
 dbTag1(0), dbTag2(0)
{
    
}

// constructor for Subclass
MP_Constraint::MP_Constraint(int tag, int nodeRetain, int nodeConstr, 
			     ID &constrainedDOF, 
			     ID &retainedDOF, int clasTag)
:DomainComponent(tag, clasTag),
 nodeRetained(nodeRetain), nodeConstrained(nodeConstr), 
 constraint(0), constrDOF(0), retainDOF(0),  dbTag1(0), dbTag2(0)
{
    constrDOF = new ID(constrainedDOF);
    retainDOF = new ID(retainedDOF);    
    if (constrDOF == 0 || constrainedDOF.Size() != constrDOF->Size() ||
	retainDOF == 0 || retainedDOF.Size() != retainDOF->Size()) { 
	cerr << "MP_Constraint::MP_Constraint - ran out of memory 1\n";
	exit(-1);
    }    
}


// general constructor for ModelBuilder
MP_Constraint::MP_Constraint(int tag, int nodeRetain, int nodeConstr, Matrix &constr,
			     ID &constrainedDOF, ID &retainedDOF)
:DomainComponent(tag, CNSTRNT_TAG_MP_Constraint), 
 nodeRetained(nodeRetain), nodeConstrained(nodeConstr), 
 constraint(0), constrDOF(0), retainDOF(0), dbTag1(0), dbTag2(0)
{
    
    constrDOF = new ID(constrainedDOF);
    retainDOF = new ID(retainedDOF);    
    if (constrDOF == 0 || constrainedDOF.Size() != constrDOF->Size() ||
	retainDOF == 0 || retainedDOF.Size() != retainDOF->Size()) { 
	cerr << "MP_Constraint::MP_Constraint - ran out of memory 1\n";
	exit(-1);
    }    
    
    constraint = new Matrix(constr);
    if (constraint == 0 || constr.noCols() != constr.noCols()) { 
	cerr << "MP_Constraint::MP_Constraint - ran out of memory 2\n";
	exit(-1);
    }        
}



MP_Constraint::~MP_Constraint()
{
    // invoke the destructor on the matrix and the two ID objects
    if (constraint != 0)
	delete constraint;
    if (constrDOF != 0)
	delete constrDOF;
    if (retainDOF != 0)
	delete retainDOF;    
}


int
MP_Constraint::getNodeRetained(void) const
{
    // return id of retained node
    return nodeRetained;
}

int
MP_Constraint::getNodeConstrained(void) const
{
    // return id of constrained node    
    return nodeConstrained;
}


const ID &
MP_Constraint::getConstrainedDOFs(void) const
{
    if (constrDOF == 0) {
	cerr << "MP_Constraint::getConstrainedDOF - no ID was set, ";
	cerr << "was recvSelf() ever called? or subclass incorrect?\n";	
	exit(-1);
    }

    // return the ID corresponding to constrained DOF of Ccr
    return *constrDOF;    
}


const ID &
MP_Constraint::getRetainedDOFs(void) const
{
    if (retainDOF == 0) {
	cerr << "MP_Constraint::getRetainedDOFs - no ID was set\n ";
	cerr << "was recvSelf() ever called? or subclass incorrect?\n";		
	exit(-1);
    }

    // return the ID corresponding to retained DOF of Ccr
    return *retainDOF;    
}

int 
MP_Constraint::applyConstraint(double timeStamp)
{
    // does nothing MP_Constraint objects are time invariant
    return 0;
}

bool
MP_Constraint::isTimeVarying(void) const
{
    return false;
}


const Matrix &
MP_Constraint::getConstraint(void)
{
    if (constraint == 0) {
	cerr << "MP_Constraint::getConstraint - no Matrix was set\n";
	exit(-1);
    }    

    // return the constraint matrix Ccr
    return *constraint;    
}

int 
MP_Constraint::sendSelf(int cTag, Channel &theChannel)
{
    ID data(9);
    int dataTag = this->getDbTag();

    data(0) = this->getTag(); 
    data(1) = nodeRetained;
    data(2) = nodeConstrained;
    if (constraint == 0) data(3) = 0; else data(3) = constraint->noRows();
    if (constraint == 0) data(4) = 0; else data(4) = constraint->noCols();    
    if (constrDOF == 0) data(5) = 0; else data(5) = constrDOF->Size();    
    if (retainDOF == 0) data(6) = 0; else data(6) = retainDOF->Size();        
    
    // need two database tags for ID objects
    if (constrDOF != 0 && dbTag1 == 0) 
      dbTag1 = theChannel.getDbTag();
    if (retainDOF != 0 && dbTag2 == 0) 
      dbTag2 = theChannel.getDbTag();

    data(7) = dbTag1;
    data(8) = dbTag2;

    int result = theChannel.sendID(dataTag, cTag, data);
    if (result < 0) {
	cerr << "WARNING MP_Constraint::sendSelf - error sending ID data\n";
	return result;  
    }    
    
    if (constraint != 0 && constraint->noRows() != 0) {
	int result = theChannel.sendMatrix(dataTag, cTag, *constraint);
	if (result < 0) {
	    cerr << "WARNING MP_Constraint::sendSelf ";
	    cerr << "- error sending Matrix data\n"; 
	    return result;  
	}
    }

    if (constrDOF != 0 && constrDOF->Size() != 0) {
	int result = theChannel.sendID(dbTag1, cTag, *constrDOF);
	if (result < 0) {
	    cerr << "WARNING MP_Constraint::sendSelf ";
	    cerr << "- error sending constrained data\n"; 
	    return result;  
	}
    }

    if (retainDOF != 0 && retainDOF->Size() != 0) {
	int result = theChannel.sendID(dbTag2, cTag, *retainDOF);
	if (result < 0) {
	    cerr << "WARNING MP_Constraint::sendSelf ";
	    cerr << "- error sending retained data\n"; 
	    return result;  
	}
    }
    
    return 0;
}


int 
MP_Constraint::recvSelf(int cTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    ID data(9);
    int result = theChannel.recvID(dataTag, cTag, data);
    if (result < 0) {
	cerr << "WARNING MP_Constraint::recvSelf - error receiving ID data\n";
	return result;  
    }    

    this->setTag(data(0));
    nodeRetained = data(1);
    nodeConstrained = data(2);
    int numRows = data(3); 
    int numCols = data(4);
    dbTag1 = data(7);
    dbTag2 = data(8);
    
    if (numRows != 0 && numCols != 0) {
	constraint = new Matrix(numRows,numCols);
	
	int result = theChannel.recvMatrix(dataTag, cTag, *constraint);
	if (result < 0) {
	    cerr << "WARNING MP_Constraint::recvSelf ";
	    cerr << "- error receiving Matrix data\n"; 
	    return result;  
	}
    }    
    int size = data(5);
    if (size != 0) {
	constrDOF = new ID(size);
	int result = theChannel.recvID(dbTag1, cTag, *constrDOF);
	if (result < 0) {
	    cerr << "WARNING MP_Constraint::recvSelf ";
	    cerr << "- error receiving constrained data\n"; 
	    return result;  
	}	
    }
    
    size = data(6);
    if (size != 0) {
	retainDOF = new ID(size);
	int result = theChannel.recvID(dbTag2, cTag, *retainDOF);
	if (result < 0) {
	    cerr << "WARNING MP_Retainaint::recvSelf ";
	    cerr << "- error receiving retained data\n"; 
	    return result;  
	}	
    }    
    
    return 0;
}



void
MP_Constraint::Print(ostream &s, int flag)
{     
    s << "MP_Constraint: " << this->getTag() << "\n";
    s << "\tNode Constrained: " << nodeConstrained;
    s << " node Retained: " << nodeRetained ;
    if (constrDOF != 0)
	s << " constrained dof: " << *constrDOF;    
    if (retainDOF != 0)
	s << " retained dof: " << *retainDOF;        
    if (constraint != 0)
	s << " constraint matrix: " << *constraint << "\n";
}


