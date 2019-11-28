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
// $Date: 2015-03-27 9:29:32 $

// Written: Minjie Zhu
// Created: March 2015
//
// Description: This file contains the class definition for PFEMGeneralLinSOE
// PFEMGeneralLinSOE is a subclass of LinearSOE. It stores the matrix equation
// Ax=b using the sparse column-compacted storage scheme for storing the
// matrix A.
//

#include <PFEMGeneralLinSOE.h>
#include <PFEMUnifiedSolver_Hybrid.h>
#include <Vertex.h>
#include <Graph.h>
#include <algorithm>
#include <Matrix.h>
#include <ID.h>

PFEMGeneralLinSOE::PFEMGeneralLinSOE()
    :PFEMLinSOE(LinSOE_TAGS_PFEMGeneralLinSOE), N(0),
     rowInd(), colPtr(), A()
{
}

PFEMGeneralLinSOE::PFEMGeneralLinSOE(PFEMUnifiedSolver_Hybrid& theSolver)
    :PFEMLinSOE(LinSOE_TAGS_PFEMGeneralLinSOE), N(0),
     rowInd(), colPtr(), A()
{
    theSolver.setLinearSOE(*this);
    LinearSOE::setSolver(theSolver);
}

PFEMGeneralLinSOE::~PFEMGeneralLinSOE()
{
}

int
PFEMGeneralLinSOE::setMatIDs(Graph& theGraph, int Ssize, int Fsize,
			     int Isize, int Psize, int Pisize)
{
    // dimension
    N = Ssize+Fsize+Isize+Psize+Pisize;
    rowInd.clear();
    colPtr.clear();
    A.clear();
    
    // dof type
    const ID& dofType = this->getDofType();
    newDofID = this->getDofID();
    int index = 0;
    for(int i=0; i<newDofID.Size(); i++) {
    	int type = dofType(i);
	if(type < 0) continue;
	newDofID(i) = index++;
    	// if(type == 2) newDofID(i) += Ssize;
    	// else if(type == 1) newDofID(i) += Ssize+Isize;
    	// else if(type == 3) newDofID(i) += Ssize+Isize+Fsize;
    	// else if(type == 4) newDofID(i) += Ssize+Isize+Fsize+Psize;
    }

    // first column
    colPtr.push_back(0);

    // for each column
    Vertex* theVertex = 0;
    for(int j=0; j<dofType.Size(); j++) {
	theVertex = theGraph.getVertexPtr(j);
	if(theVertex == 0) {
	    opserr<<"WARNING: vertex "<<j<<" not in graph -- PFEMGeneralLinSOE::setSize\n";
	    return -1;
	}

	// start and end index of current column
	int start = colPtr.back();
	int end = start;

	// column number
	int col = theVertex->getTag();

	// column type
	int coltype = dofType(col);
	if(coltype < 0) {
	    continue;
	}
	int colid = newDofID(col);

	// diagnol term
	rowInd.push_back(colid);
	end++;

	// off diagnol terms
	const ID& theAdjacency = theVertex->getAdjacency();
	for(int i=0; i<theAdjacency.Size(); i++) {
	    int row = theAdjacency(i);
	    int rowtype = dofType(row);
	    if(rowtype < 0) continue;
	    int rowid = newDofID(row);
	    rowInd.push_back(rowid);
	    end++;
	}

	// sort row indices of this column
	std::sort(rowInd.begin()+start, rowInd.begin()+end);

	// add new column
	colPtr.push_back(end);
    }

    // resize A
    A.resize(rowInd.size(),0.0);

    return 0;
}

int 
PFEMGeneralLinSOE::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int 
PFEMGeneralLinSOE::recvSelf(int cTag, Channel &theChannel, 
			    FEM_ObjectBroker &theBroker)  
{
    return 0;
}

void
PFEMGeneralLinSOE::zeroA()
{
    for(int i=0; i<(int)A.size(); i++) {
	A[i] = 0.0;
    }
}

int
PFEMGeneralLinSOE::addA(const Matrix& mat, const ID& id, double fact)
{
    // quick return
    if(fact == 0.0) return 0;

    // check that m and id are of similar size
    if(id.Size() != mat.noRows() && id.Size() != mat.noCols()) {
	opserr << "PFEMGeneralLinSOE::addA() ";
	opserr << " - Matrix and ID not of similar sizes\n";
	return -1;
    }

    const ID& dofType = this->getDofType();
    
    if(fact == 1.0) {
	for(int j=0; j<id.Size(); j++) {
	    // each column
	    int col = id(j);
	    if(col>=dofType.Size() || col<0) continue;
	    int coltype = dofType(col);
	    if(coltype < 0) continue;
	    int colid = newDofID(col);

	    for(int i=0; i<id.Size(); i++) {
		// each row
		int row = id(i);
		if(row>=dofType.Size() || row<0) continue;
		int rowtype = dofType(row);
		if(rowtype < 0) continue;
		int rowid = newDofID(row);

		// place in sparse matrix
		for(int k=colPtr[colid]; k<colPtr[colid+1]; k++) {
		    if(rowInd[k] == rowid) {
			A[k] += mat(i,j);
			break;
		    }
		}
	    }
	}
    } else {
	for(int j=0; j<id.Size(); j++) {
	    // each column
	    int col = id(j);
	    if(col>=dofType.Size() || col<0) continue;
	    int coltype = dofType(col);
	    if(coltype < 0) continue;
	    int colid = newDofID(col);

	    for(int i=0; i<id.Size(); i++) {
		// each row
		int row = id(i);
		if(row>=dofType.Size() || row<0) continue;
		int rowtype = dofType(row);
		if(rowtype < 0) continue;
		int rowid = newDofID(row);

		// place in sparse matrix
		for(int k=colPtr[colid]; k<colPtr[colid+1]; k++) {
		    if(rowInd[k] == rowid) {
			A[k] += fact*mat(i,j);
			break;
		    }
		}
	    }
	}
    }

    return 0;
}
