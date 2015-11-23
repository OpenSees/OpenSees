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
                                                                        
// $Revision: 1.7 $
// $Date: 2009-05-11 20:56:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSOE.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for 
// UmfpackGenLinSolver. It solves the UmfpackGenLinSOEobject by calling
// UMFPACK5.7.1 routines.
//
// What: "@(#) UmfpackGenLinSolver.h, revA"

#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <stdlib.h>


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>

UmfpackGenLinSOE::UmfpackGenLinSOE(UmfpackGenLinSolver &the_Solver)
    :LinearSOE(the_Solver, LinSOE_TAGS_UmfpackGenLinSOE), X(), B(), Ap(), Ai(), Ax()
{
    the_Solver.setLinearSOE(*this);
}


UmfpackGenLinSOE::UmfpackGenLinSOE()
    :LinearSOE(LinSOE_TAGS_UmfpackGenLinSOE), X(), B(), Ap(), Ai(), Ax()
{
}


UmfpackGenLinSOE::~UmfpackGenLinSOE()
{
}


int
UmfpackGenLinSOE::getNumEqn(void) const
{
    return X.Size();
}

int
UmfpackGenLinSOE::setSize(Graph &theGraph)
{
    int size = theGraph.getNumVertex();
    if (size < 0) {
	opserr<<"size of soe < 0\n";
	return -1;
    }

    // fist itearte through the vertices of the graph to get nnz
    Vertex *theVertex;
    int nnz = 0;
    VertexIter &theVertices = theGraph.getVertices();
    while ((theVertex = theVertices()) != 0) {
	const ID &theAdjacency = theVertex->getAdjacency();
	nnz += theAdjacency.Size() +1; // the +1 is for the diag entry
    }

    // resize A, B, X
    Ap.reserve(size+1);
    Ai.reserve(nnz);
    Ax.resize(nnz,0.0);
    B.resize(size);
    B.Zero();
    X.resize(size);
    X.Zero();

    // fill in Ai and Ap
    Ap.push_back(0);
    for (int a=0; a<size; a++) {

	theVertex = theGraph.getVertexPtr(a);
	if (theVertex == 0) {
	    opserr << "WARNING:UmfpackGenLinSOE::setSize :";
	    opserr << " vertex " << a << " not in graph! - size set to 0\n";
	    size = 0;
	    return -1;
	}

	const ID &theAdjacency = theVertex->getAdjacency();
	int idSize = theAdjacency.Size();
	ID col(0,idSize+1);

	// diagonal
	col.insert(theVertex->getTag());

	// now we have to place the entries in the ID into order in Ai
	for (int i=0; i<idSize; i++) {
	    int row = theAdjacency(i);
	    col.insert(row);
	}

	// copy to Ai
	for (int i=0; i<col.Size(); i++) {
	    Ai.push_back(col(i));
	}

	// set Ap
	Ap.push_back(Ap[a]+col.Size());
    }

    // invoke setSize() on the Solver
    LinearSOESolver *the_Solver = this->getSolver();
    int solverOK = the_Solver->setSize();
    if (solverOK < 0) {
	opserr << "WARNING:UmfpackGenLinSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }
    return 0;
}

int
UmfpackGenLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return
    if (fact == 0.0) return 0;

    int idSize = id.Size();

    // check that m and id are of similar size
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "UmfpackGenLinSOE::addA() ";
	opserr << " - Matrix and ID not of similar sizes\n";
	return -1;
    }

    int size = X.Size();
    if (fact == 1.0) { // do not need to multiply
	for (int j=0; j<idSize; j++) {
	    int col = id(j);
	    if (col<0 || col>=size) {
		continue;
	    }
	    for (int i=0; i<idSize; i++) {
		int row = id(i);
		if (row<0 || row>=size) {
		    continue;
		}

		// find place in A
		for (int k=Ap[col]; k<Ap[col+1]; k++) {
		    if (Ai[k] == row) {
			Ax[k] += m(i,j);
			break;
		    }
		}
	    }
	}
    } else {
	for (int j=0; j<idSize; j++) {
	    int col = id(j);
	    if (col<0 || col>=X.Size()) {
		continue;
	    }
	    for (int i=0; i<idSize; i++) {
		int row = id(i);
		if (row<0 || row>=X.Size()) {
		    continue;
		}

		// find place in A
		for (int k=Ap[col]; k<Ap[col+1]; k++) {
		    if (Ai[k] == row) {
			Ax[k] += fact*m(i,j);
			break;
		    }
		}
	    }
	}
    }

    return 0;
}


int
UmfpackGenLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();    
    // check that m and id are of similar size
    if (idSize != v.Size() ) {
	opserr << "UmfpackGenLinSOE::addB() ";
	opserr << " - Vector and ID not of similar sizes\n";
	return -1;
    }    

    int size = B.Size();
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0) B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0) B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0) B[pos] += v(i) * fact;
	}
    }
    
    return 0;
}


int
UmfpackGenLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  {
	B.Zero();
	return 0;
    }

    int size = B.Size();
    if (v.Size() != size) {
	opserr << "WARNING BandGenLinSOE::setB() -";
	opserr << " incomptable sizes " << size << " and " << v.Size() << endln;
	return -1;
    }

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<size; i++) {
	    B[i] = v(i);
	}
    } else if (fact == -1.0) {
	for (int i=0; i<size; i++) {
	    B[i] = -v(i);
	}
    } else {
	for (int i=0; i<size; i++) {
	    B[i] = v(i) * fact;
	}
    }
    
    return 0;
}

void
UmfpackGenLinSOE::zeroA(void)
{
    Ax.assign(Ax.size(),0.0);
}

void
UmfpackGenLinSOE::zeroB(void)
{
    B.Zero();
}

void
UmfpackGenLinSOE::setX(int loc, double value)
{
    if (loc<X.Size() && loc>=0) {
	X(loc) = value;
    }
}


void
UmfpackGenLinSOE::setX(const Vector &x)
{
    if (x.Size() == X.Size()) {
	X = x;
    }
}


const Vector &
UmfpackGenLinSOE::getX(void)
{
    return X;
}

const Vector &
UmfpackGenLinSOE::getB(void)
{
    return B;
}

double
UmfpackGenLinSOE::normRHS(void)
{
    return B.Norm();
}


int
UmfpackGenLinSOE::setUmfpackGenLinSolver(UmfpackGenLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    if (X.Size() != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:UmfpackGenLinSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return -1;
	}
    }
    return this->LinearSOE::setSolver(newSolver);
}


int
UmfpackGenLinSOE::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int
UmfpackGenLinSOE::recvSelf(int cTag, Channel &theChannel,
			   FEM_ObjectBroker &theBroker)
{
    return 0;
}
