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
                                                                        
// $Revision: 1.6 $
// $Date: 2009-05-20 17:30:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandSPD/BandSPDLinSOE.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: February 1997

// Description: This file contains the implementation for BandSPDLinSOE

#include <stdlib.h>

#include <BandSPDLinSOE.h>
#include <BandSPDLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
//#include <f2c.h>
#include <math.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <iostream>
using std::nothrow;

BandSPDLinSOE::BandSPDLinSOE(BandSPDLinSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_BandSPDLinSOE),
 size(0), half_band(0), A(0), B(0), X(0), vectX(0), vectB(0),
 Asize(0), Bsize(0),
 factored(false)
{
    the_Solver.setLinearSOE(*this);
}


BandSPDLinSOE::BandSPDLinSOE(BandSPDLinSolver &the_Solver, int classTag)
:LinearSOE(the_Solver, classTag),
 size(0), half_band(0), A(0), B(0), X(0), vectX(0), vectB(0),
 Asize(0), Bsize(0),
 factored(false)
{

}


BandSPDLinSOE::BandSPDLinSOE(int classTag)
:LinearSOE(classTag),
 size(0), half_band(0), A(0), B(0), X(0), vectX(0), vectB(0),
 Asize(0), Bsize(0),
 factored(false)
{

}


BandSPDLinSOE::BandSPDLinSOE(int N, int numSuper,
			     BandSPDLinSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_BandSPDLinSOE),
 size(0), half_band(0), A(0), B(0), X(0), vectX(0), vectB(0),
 Asize(0), Bsize(0),
 factored(false)
{
    size = N;
    half_band = numSuper+1;

    A = new (nothrow) double[half_band*size];
	
    if (A == 0) {
	opserr << "FATAL:BandSPDLinSOE::BandSPDLinSOE :";
	opserr << " ran out of memory for A (size,ku) (";
	opserr << size <<", " << numSuper << ") \n";
	size = 0; Asize = 0;
    } else {

	// zero the matrix
	Asize = half_band*size;
	for (int j=0; j<Asize; j++)
	    A[j] = 0;
    
	B = new (nothrow) double[size];
	X = new (nothrow) double[size];
	
	if (B == 0 || X == 0 ) {
	    opserr << "FATAL:BandSPDLinSOE::BandSPDLinSOE :";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << size << ") \n";
	    size = 0; Bsize = 0;
	}
	
	// zero the vectors
	for (int i=0; i<size; i++) {
	    B[i] = 0;
	    X[i] = 0;
	}
	
    }

    vectX = new Vector(X,size);
    vectB = new Vector(B,size);    
    

    the_Solver.setLinearSOE(*this);    
    
    int solverOK = the_Solver.setSize();
    if (solverOK < 0) {
	opserr << "FATAL:BandSPDLinSOE::BandSPDLinSOE :";
	opserr << " solver failed setSize() in constructor\n";
    }   	 
}

    
BandSPDLinSOE::~BandSPDLinSOE()
{
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;        
}


int 
BandSPDLinSOE::getNumEqn(void) const
{
    return size;
}

int 
BandSPDLinSOE::setSize(Graph &theGraph)
{
    int result = 0;
    int oldSize = size;
    size = theGraph.getNumVertex();
    half_band = 0;
    
    Vertex *vertexPtr;
    VertexIter &theVertices = theGraph.getVertices();
    
    while ((vertexPtr = theVertices()) != 0) {
	int vertexNum = vertexPtr->getTag();
	const ID &theAdjacency = vertexPtr->getAdjacency();
	for (int i=0; i<theAdjacency.Size(); i++) {
	    int otherNum = theAdjacency(i);
	    int diff = vertexNum-otherNum;
	    if (half_band < diff)
		half_band = diff;
	}
    }
    half_band += 1; // include the diagonal
     
    if (half_band*size > Asize) { // we have to get another space for A

	if (A != 0) 
	    delete [] A;

	A = new (nothrow) double[half_band*size];
	
        if (A == 0) {
            opserr << "WARNING BandSPDLinSOE::setSize() :";
	    opserr << " ran out of memory for A (size,ku) (";
	    opserr << size <<", " << half_band-1 << ") \n";
	    Asize = 0; size = 0;
	    result = -1;
        }
	else
	    Asize = half_band*size;
    }

    // zero the matrix
    for (int i=0; i<half_band*size; i++)
	A[i] = 0;
	
    factored = false;
    
    if (size > Bsize) { // we have to get another space for A
	
	// delete the old	
	if (B != 0) delete [] B;
	if (X != 0) delete [] X;

	// create the new
	B = new (nothrow) double[size];
	X = new (nothrow) double[size];
	
        if (B == 0 || X == 0) {
            opserr << "WARNING BandSPDLinSOE::setSize():";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << size << ") \n";
	    Bsize = 0; size =0;
	    result = -1;
        }
    }

    // zero the vectors
    for (int j=0; j<size; j++) {
	B[j] = 0;
	X[j] = 0;
    }

    if (size != oldSize) {

	if (vectX != 0)
	    delete vectX;
	if (vectB != 0)
	    delete vectB;

	vectX = new Vector(X,size);	
	vectB = new Vector(B,size);

	if (size > Bsize)
	    Bsize = size;
    }
    
    // invoke setSize() on the Solver
    LinearSOESolver *the_Solver = this->getSolver();
    int solverOK = the_Solver->setSize();
    if (solverOK < 0) {
	opserr << "WARNING:BandSPDLinSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }    

    return result;    
}

int 
BandSPDLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;
    
    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "BandSPDLinSOE::addA()	- Matrix and ID not of similar sizes\n";
	return -1;
    }

    if (fact == 1.0) { // do not need to multiply 
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *coliiPtr = A +(col+1)*half_band -1;
		int minColRow = col - half_band + 1;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0 && 
			row <= col && row >= minColRow) { // only add upper
			 double *APtr = coliiPtr + (row-col);
			 *APtr += m(j,i);
		     }
		}  // for j
	    } 
	}  // for i
    } else {
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *coliiPtr = A +(col+1)*half_band -1;
		int minColRow = col - half_band +1;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0 && 
			row <= col && row >= minColRow) { // only add upper
			 double *APtr = coliiPtr + (row-col);
			 *APtr += m(j,i)*fact;
		     }
		}  // for j
	    } 
	}  // for i
    }
    return 0;
}


int 
BandSPDLinSOE::addColA(const Vector &colData, int col, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  return 0;
  
  if (colData.Size() != size) {
    opserr << "BanSPDLinSOE::addColA() - colData size not equal to n\n";
    return -1;
  }

  if (col > size && col < 0) {
    opserr << "BandSPDLinSOE::addColA() - col " << col << "outside range 0 to " << size << endln;
    return -1;
  }
  
  if (fact == 1.0) { // do not need to multiply 
    
    double *coliiPtr = A +(col+1)*half_band -1;
    int minColRow = col - half_band + 1;
    for (int row=0; row<size; row++) {
      if (row <size && row >= 0 && 
	  row <= col && row >= minColRow) { // only add upper
	double *APtr = coliiPtr + (row-col);
	*APtr += colData(row);
      }
    }  

  } else {

    double *coliiPtr = A +(col+1)*half_band -1;
    int minColRow = col - half_band + 1;
    for (int row=0; row<size; row++) {
      if (row <size && row >= 0 && 
	  row <= col && row >= minColRow) { // only add upper
	double *APtr = coliiPtr + (row-col);
	*APtr += colData(row) * fact;
      }
    }  
  }

  return 0;
}

    
int 
BandSPDLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    
    // check for a quick return 
    if (fact == 0.0)  return 0;

    // check that m and id are of similar size
    int idSize = id.Size();        
    if (idSize != v.Size() ) {
	opserr << "BandSPDLinSOE::addB()	- Vector and ID not of similar sizes\n";
	return -1;
    }    
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] -= v(i);
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
BandSPDLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    if (v.Size() != size) {
	opserr << "WARNING BandGenLinSOE::setB() -";
	opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
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
BandSPDLinSOE::zeroA(void)
{
    double *Aptr = A;
    int theSize = Asize; 
    for (int i=0; i<theSize; i++)
	*Aptr++ = 0;
    
    factored = false;
}
	
void 
BandSPDLinSOE::zeroB(void)
{
    double *Bptr = B;
    for (int i=0; i<size; i++)
	*Bptr++ = 0;
}


void 
BandSPDLinSOE::setX(int loc, double value)
{
    if (loc < size && loc >= 0)
	X[loc] = value;
}

void 
BandSPDLinSOE::setX(const Vector &x)
{
    if (x.Size() == size && vectX != 0)
      *vectX = x;
}




const Vector &
BandSPDLinSOE::getX(void)
{
    if (vectX == 0) {
	opserr << "FATAL BandSPDLinSOE::getX - vectX == 0";
	exit(-1);
    }
    return *vectX;
}

const Vector &
BandSPDLinSOE::getB(void)
{
    if (vectB == 0) {
	opserr << "FATAL BandSPDLinSOE::getB - vectB == 0";
	exit(-1);
    }    
    return *vectB;
}

double 
BandSPDLinSOE::normRHS(void)
{
    double norm =0.0;
    for (int i=0; i<size; i++) {
	double Bi = B[i];
	norm += Bi*Bi;
    }
    return sqrt(norm);
}    


int
BandSPDLinSOE::setBandSPDSolver(BandSPDLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    
    if (size != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:BandSPDLinSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return solverOK;
	}
    }

    return this->LinearSOE::setSolver(newSolver);
}


int 
BandSPDLinSOE::sendSelf(int tag, Channel &theChannel)
{
  opserr << "BandSPDLinSOE::sendSelf() - not implemented\n";
  return -1;
}


int 
BandSPDLinSOE::recvSelf(int tag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  opserr << "BandSPDLinSOE::recvSelf( - not implemented\n";
  return -1;
}
