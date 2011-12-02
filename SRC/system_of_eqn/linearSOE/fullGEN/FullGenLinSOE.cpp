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
// $Date: 2000-09-15 08:23:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/fullGEN/FullGenLinSOE.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/fullGEN/FullGenLinSOE.C
//
// Written: fmk 
// Created: Febuary 1997
// Revision: A
//
// Description: This file contains the implementation for FullGenLinSOE


#include <FullGenLinSOE.h>
#include <FullGenLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>


#include <Channel.h>
#include <FEM_ObjectBroker.h>

FullGenLinSOE::FullGenLinSOE(FullGenLinSolver &theSolvr)
:LinearSOE(theSolvr, LinSOE_TAGS_FullGenLinSOE),
 size(0), A(0), B(0), X(0), vectX(0), vectB(0),
 Asize(0), Bsize(0),
 factored(false)
{
    theSolvr.setLinearSOE(*this);
}


FullGenLinSOE::FullGenLinSOE(int N, FullGenLinSolver &theSolvr)
:LinearSOE(theSolvr, LinSOE_TAGS_FullGenLinSOE),
 size(0), A(0), B(0), X(0), vectX(0), vectB(0),
 Asize(0), Bsize(0),
 factored(false)
{
    size = N;

    A = new double[size*size];
	
    if (A == 0) {
	cerr << "WARNING :FullGenLinSOE::FullGenLinSOE :";
	cerr << " ran out of memory for A (size,size) (";
	cerr << size <<", " << size << ") \n";
	size = 0; 
    } else {
	// zero the matrix
	Asize = size*size;
	for (int i=0; i<Asize; i++)
	    A[i] = 0;
    
	B = new double[size];
	X = new double[size];
	
	if (B == 0 || X == 0) {
	    cerr << "WARNING :FullGenLinSOE::FullGenLinSOE :";
	    cerr << " ran out of memory for vectors (size) (";
	    cerr << size << ") \n";
	    size = 0; Bsize = 0;
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
    
    theSolvr.setLinearSOE(*this);
    
    // invoke setSize() on the Solver        
    if (theSolvr.setSize() < 0) {
	cerr << "WARNING :FullGenLinSOE::FullGenLinSOE :";
	cerr << " solver failed setSize() in constructor\n";
    }    
    
}

    
FullGenLinSOE::~FullGenLinSOE()
{
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;        
}


int
FullGenLinSOE::getNumEqn(void) const
{
    return size;
}

int 
FullGenLinSOE::setSize(Graph &theGraph)
{
    int result = 0;
    int oldSize = size;
    size = theGraph.getNumVertex();

    if (size*size > Asize) { // we have to get another space for A

	if (A != 0) 
	    delete [] A;

	A = new double[size*size];
	
        if (A == 0) {
            cerr << "WARNING FullGenLinSOE::FullGenLinSOE :";
	    cerr << " ran out of memory for A (size,size) (";
	    cerr << size <<", " << size << ") \n";
	    size = 0; Asize = 0;
	    result =  -1;
        } else
	    Asize = size*size;
    }

    // zero the matrix
    for (int i=0; i<Asize; i++)
	A[i] = 0;
	
    factored = false;
    
    if (size > Bsize) { // we have to get space for the vectors
	
	// delete the old	
	if (B != 0) delete [] B;
	if (X != 0) delete [] X;

	// create the new
	B = new double[size];
	X = new double[size];
	
        if (B == 0 || X == 0) {
            cerr << "WARNING FullGenLinSOE::FullGenLinSOE :";
	    cerr << " ran out of memory for vectors (size) (";
	    cerr << size << ") \n";
	    size = 0; Bsize = 0;
	    result =  -1;
        }
	else
	    Bsize = size;
    }

    // zero the vectors
    for (int j=0; j<Bsize; j++) {
	B[j] = 0;
	X[j] = 0;
    }

    // create new Vectors
    if (size != oldSize) {
	if (vectX != 0)
	    delete vectX;

	if (vectB != 0)
	    delete vectB;
	
	vectX = new Vector(X,Bsize);
	vectB = new Vector(B,Bsize);	

    }

    // invoke setSize() on the Solver    
    LinearSOESolver *theSolvr = this->getSolver();
    int solverOK = theSolvr->setSize();
    if (solverOK < 0) {
	cerr << "WARNING:FullGenLinSOE::setSize :";
	cerr << " solver failed setSize()\n";
	return solverOK;
    }    
    
    return result;
}

int 
FullGenLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();
    
    // check that m and id are of similar size
    if (idSize != m.noRows() && idSize != m.noCols()) {
	cerr << "FullGenLinSOE::addA()	- Matrix and ID not of similar sizes\n";
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply 
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *startColiPtr = A + col*size;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0) {
			 double *APtr = startColiPtr + row;
			 *APtr += m(j,i);
		     }
		}  // for j
	    } 
	}  // for i
    } else {
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *startColiPtr = A + col*size;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0) {
			 double *APtr = startColiPtr + row;
			 *APtr += m(j,i) * fact;
		     }
		}  // for j
	    } 
	}  // for i
    }    
    return 0;
}




int 
FullGenLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();    
    // check that m and id are of similar size
    if (idSize != v.Size() ) {
	cerr << "FullGenLinSOE::addB()	- Vector and ID not of similar sizes\n";
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
FullGenLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    if (v.Size() != size) {
	cerr << "WARNING BandGenLinSOE::setB() -";
	cerr << " incomptable sizes " << size << " and " << v.Size() << endl;
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
FullGenLinSOE::zeroA(void)
{
    double *Aptr = A;
    int theSize = size*size;
    for (int i=0; i<theSize; i++)
	*Aptr++ = 0;

    factored = false;
}
	
void 
FullGenLinSOE::zeroB(void)
{
    double *Bptr = B;
    for (int i=0; i<size; i++)
	*Bptr++ = 0;
}

void 
FullGenLinSOE::setX(int loc, double value)
{
    if (loc < size && loc >=0)
	X[loc] = value;
}

const Vector &
FullGenLinSOE::getX(void)
{
    if (vectX == 0) {
	cerr << "FATAL FullGenLinSOE::getX - vectX == 0";
	exit(-1);
    }
    return *vectX;
}

const Vector &
FullGenLinSOE::getB(void)
{
    if (vectB == 0) {
	cerr << "FATAL FullGenLinSOE::getB - vectB == 0";
	exit(-1);
    }        
    return *vectB;
}

double 
FullGenLinSOE::normRHS(void)
{
    double norm =0.0;
    for (int i=0; i<size; i++) {
	double Yi = B[i];
	norm += Yi*Yi;
    }
    return sqrt(norm);
    
}    


int
FullGenLinSOE::setFullGenSolver(FullGenLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    
    if (size != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    cerr << "WARNING:FullGenLinSOE::setSolver :";
	    cerr << "the new solver could not setSeize() - staying with old\n";
	    return -1;
	}
    }
    
    return this->LinearSOE::setSolver(newSolver);
}


int 
FullGenLinSOE::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int 
FullGenLinSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}






