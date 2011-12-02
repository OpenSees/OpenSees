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
// $Date: 2006-01-10 00:42:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/SparseGenRowLinSOE.cpp,v $
                                                                        
// Written: fmk 
// Created: 04/05
//
// Description: This file contains the implementation for SparseGenRowLinSOE

#include <SparseGenRowLinSOE.h>
#include <SparseGenRowLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

SparseGenRowLinSOE::SparseGenRowLinSOE(SparseGenRowLinSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_SparseGenRowLinSOE),
 size(0), nnz(0), A(0), B(0), X(0), colA(0), rowStartA(0),
 vectX(0), vectB(0),
 Asize(0), Bsize(0),
 factored(false)
{
    the_Solver.setLinearSOE(*this);
}

/*******************************************************************************

SparseGenRowLinSOE::SparseGenRowLinSOE(int N, int NNZ, int *RowStartA, int *ColA,
				       SparseGenRowLinSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_SparseGenRowLinSOE),
 size(N), nnz(NNZ), A(0), B(0), X(0), 
 colA(ColA), rowStartA(RowStartA), 
 vectX(0), vectB(0),
 Asize(0), Bsize(0),
 factored(false)
{

    A = new double[NNZ];
	
    if (A == 0) {
	opserr << "WARNING :SparseGenRowLinSOE::SparseGenRowLinSOE :";
	opserr << " ran out of memory for A (NNZ) (";
	opserr << NNZ << ") \n";
	size = 0; nnz = 0;
    } else {
	// zero the matrix
	Asize = NNZ;
	for (int i=0; i<Asize; i++)
	    A[i] = 0;
    
	B = new double[size];
	X = new double[size];
	
	if (B == 0 || X == 0) {
	    opserr << "WARNING :SparseGenRowLinSOE::SparseGenRowLinSOE :";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << size << ") \n";
	    size = 0; Bsize = 0; NNZ = 0;
	} else {
	    Bsize = size;
	    // zero the vectors
	    for (int j=0; j<size; j++) {
		B[j] = 0;
		X[j] = 0;
	    }
	}
    }

    vectX = new Vector(X,size);
    vectB = new Vector(B,size);    

    the_Solver.setLinearSOE(*this);
    
    // invoke setSize() on the Solver        
    if (the_Solver.setSize() < 0) {
	opserr << "WARNING :SparseGenRowLinSOE::SparseGenRowLinSOE :";
	opserr << " solver failed setSize() in constructor\n";
    }    
}

*******************************************************************************/    

SparseGenRowLinSOE::~SparseGenRowLinSOE()
{
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (rowStartA != 0) delete [] rowStartA;
    if (colA != 0) delete []colA;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;        
}


int
SparseGenRowLinSOE::getNumEqn(void) const
{
    return size;
}

int 
SparseGenRowLinSOE::setSize(Graph &theGraph)
{

    int result = 0;
    int oldSize = size;
    size = theGraph.getNumVertex();

    // fist itearte through the vertices of the graph to get nnz
    Vertex *theVertex;
    int newNNZ = 0;
    VertexIter &theVertices = theGraph.getVertices();
    while ((theVertex = theVertices()) != 0) {
	const ID &theAdjacency = theVertex->getAdjacency();
	newNNZ += theAdjacency.Size() +1; // the +1 is for the diag entry
    }
    nnz = newNNZ;

    if (newNNZ > Asize) { // we have to get more space for A and colA
	if (A != 0) 
	    delete [] A;
	if (colA != 0)
	    delete [] colA;
	
	A = new double[newNNZ];
	colA = new int[newNNZ];
	
        if (A == 0 || colA == 0) {
            opserr << "WARNING SparseGenRowLinSOE::SparseGenRowLinSOE :";
	    opserr << " ran out of memory for A and colA with nnz = ";
	    opserr << newNNZ << " \n";
	    size = 0; Asize = 0; nnz = 0;
	    result =  -1;
        } 
	
	Asize = newNNZ;
    }

    // zero the matrix
    for (int i=0; i<Asize; i++)
	A[i] = 0;
	
    factored = false;
    
    if (size > Bsize) { // we have to get space for the vectors
	
	// delete the old	
	if (B != 0) delete [] B;
	if (X != 0) delete [] X;
	if (rowStartA != 0) delete [] rowStartA;

	// create the new
	B = new double[size];
	X = new double[size];
	rowStartA = new int[size+1]; 
	
        if (B == 0 || X == 0 || rowStartA == 0) {
            opserr << "WARNING SparseGenRowLinSOE::SparseGenRowLinSOE :";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << size << ") \n";
	    size = 0; Bsize = 0;
	    result =  -1;
        }
	else
	    Bsize = size;
    }

    // zero the vectors
    for (int j=0; j<size; j++) {
	B[j] = 0;
	X[j] = 0;
    }
    
    // create new Vectors objects
    if (size != oldSize) {
	if (vectX != 0)
	    delete vectX;

	if (vectB != 0)
	    delete vectB;
	
	vectX = new Vector(X,size);
	vectB = new Vector(B,size);	
    }

    // fill in rowStartA and colA
    if (size != 0) {
      rowStartA[0] = 0;
      int startLoc = 0;
      int lastLoc = 0;
      for (int a=0; a<size; a++) {

	theVertex = theGraph.getVertexPtr(a);
	if (theVertex == 0) {
	  opserr << "WARNING:SparseGenRowLinSOE::setSize :";
	  opserr << " vertex " << a << " not in graph! - size set to 0\n";
	  size = 0;
	  return -1;
	}

	colA[lastLoc++] = theVertex->getTag(); // place diag in first
	const ID &theAdjacency = theVertex->getAdjacency();
	int idSize = theAdjacency.Size();
	
	// now we have to place the entries in the ID into order in colA
	for (int i=0; i<idSize; i++) {

	  int row = theAdjacency(i);
	  bool foundPlace = false;
	  // find a place in colA for current col
	  for (int j=startLoc; j<lastLoc; j++)
	    if (colA[j] > row) { 
	      // move the entries already there one further on
	      // and place col in current location
	      for (int k=lastLoc; k>j; k--)
		
		colA[k] = colA[k-1];
	      colA[j] = row;
	      foundPlace = true;
	      j = lastLoc;
	    }
	  if (foundPlace == false) // put in at the end
	    colA[lastLoc] = row;

	  lastLoc++;
	}
	rowStartA[a+1] = lastLoc;;	    
	startLoc = lastLoc;
      }
    }
    
    // invoke setSize() on the Solver   
     LinearSOESolver *the_Solver = this->getSolver();
    int solverOK = the_Solver->setSize();
    if (solverOK < 0) {
	opserr << "WARNING:SparseGenRowLinSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }    
    return result;
}

int 
SparseGenRowLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  
	return 0;

    int idSize = id.Size();
    
    // check that m and id are of similar size
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "SparseGenRowLinSOE::addA() ";
	opserr << " - Matrix and ID not of similar sizes\n";
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply 
	for (int i=0; i<idSize; i++) {
	    int row = id(i);
	    if (row < size && row >= 0) {
		int startRowLoc = rowStartA[row];
		int endRowLoc = rowStartA[row+1];
		for (int j=0; j<idSize; j++) {
		    int col = id(j);
		    if (col <size && col >= 0) {
			// find place in A using colA
			for (int k=startRowLoc; k<endRowLoc; k++)
			    if (colA[k] == col) {
				A[k] += m(i,j);
				k = endRowLoc;
			    }
		     }
		}  // for j		
	    } 
	}  // for i
    } else {
	for (int i=0; i<idSize; i++) {
	    int row = id(i);
	    if (row < size && row >= 0) {
		int startRowLoc = rowStartA[row];
		int endRowLoc = rowStartA[row+1];
		for (int j=0; j<idSize; j++) {
		    int col = id(j);
		    if (col <size && col >= 0) {
			// find place in A using colA
			for (int k=startRowLoc; k<endRowLoc; k++)
			    if (colA[k] == col) {
				A[k] += fact * m(i,j);
				k = endRowLoc;
			    }
		     }
		}  // for j		
	    } 
	}  // for i
    }
    return 0;
}

    
int 
SparseGenRowLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();    
    // check that m and id are of similar size
    if (idSize != v.Size() ) {
	opserr << "SparseGenRowLinSOE::addB() ";
	opserr << " - Vector and ID not of similar sizes\n";
	return -1;
    }    

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] -= v(i) * fact;
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i) * fact;
	}
    }	

    return 0;
}


int
SparseGenRowLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


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
SparseGenRowLinSOE::zeroA(void)
{
    double *Aptr = A;
    for (int i=0; i<Asize; i++)
	*Aptr++ = 0;

    factored = false;
}
	
void 
SparseGenRowLinSOE::zeroB(void)
{
    double *Bptr = B;
    for (int i=0; i<size; i++)
	*Bptr++ = 0;
}

void 
SparseGenRowLinSOE::setX(int loc, double value)
{
    if (loc < size && loc >=0)
	X[loc] = value;
}

void 
SparseGenRowLinSOE::setX(const Vector &x)
{
  if (x.Size() == size && vectX != 0)
    *vectX = x;
}

const Vector &
SparseGenRowLinSOE::getX(void)
{
    if (vectX == 0) {
	opserr << "FATAL SparseGenRowLinSOE::getX - vectX == 0";
	exit(-1);
    }
    return *vectX;
}

const Vector &
SparseGenRowLinSOE::getB(void)
{
    if (vectB == 0) {
	opserr << "FATAL SparseGenRowLinSOE::getB - vectB == 0";
	exit(-1);
    }        
    return *vectB;
}

double 
SparseGenRowLinSOE::normRHS(void)
{
    double norm =0.0;
    for (int i=0; i<size; i++) {
	double Yi = B[i];
	norm += Yi*Yi;
    }
    return sqrt(norm);
    
}    


int
SparseGenRowLinSOE::setSparseGenRowSolver(SparseGenRowLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    
    if (size != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:SparseGenRowLinSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return -1;
	}
    }
    
    return this->LinearSOE::setSolver(newSolver);
}


int 
SparseGenRowLinSOE::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int 
SparseGenRowLinSOE::recvSelf(int cTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)  
{
    return 0;
}

