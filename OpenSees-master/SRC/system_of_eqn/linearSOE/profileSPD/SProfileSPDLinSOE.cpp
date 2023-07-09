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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the implementation for ProfileSPDLinSOE

#include <stdlib.h>
#include <SProfileSPDLinSOE.h>
#include <SProfileSPDLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <iostream>
using std::nothrow;

SProfileSPDLinSOE::SProfileSPDLinSOE(SProfileSPDLinSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_SProfileSPDLinSOE),
 size(0), profileSize(0), A(0), B(0), X(0), doubleB(0), doubleX(0),
 vectX(0), vectB(0), iDiagLoc(0), Asize(0), Bsize(0), isAfactored(false), 
 isAcondensed(false), numInt(0) 
{
    the_Solver.setLinearSOE(*this);
}


SProfileSPDLinSOE::SProfileSPDLinSOE(int classTag)
:LinearSOE(classTag),
 size(0), profileSize(0), A(0), B(0), X(0), doubleB(0), doubleX(0),
 vectX(0), vectB(0), iDiagLoc(0), Asize(0), Bsize(0), isAfactored(false), 
 isAcondensed(false), numInt(0) 
{

}


SProfileSPDLinSOE::SProfileSPDLinSOE(SProfileSPDLinSolver &the_Solver, int classTag)
:LinearSOE(the_Solver, classTag),
 size(0), profileSize(0), A(0), B(0), X(0), doubleB(0), doubleX(0),
 vectX(0), vectB(0), iDiagLoc(0), Asize(0), Bsize(0), isAfactored(false), 
 isAcondensed(false), numInt(0) 
{
    the_Solver.setLinearSOE(*this);
}


SProfileSPDLinSOE::SProfileSPDLinSOE(int N, int *iLoc,
					       SProfileSPDLinSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_SProfileSPDLinSOE),
 size(0), profileSize(0), A(0), B(0), X(0), vectX(0), vectB(0),
 iDiagLoc(0), Asize(0), Bsize(0), isAfactored(false), isAcondensed(false),
 numInt(0)
{
    size = N;
    profileSize = iLoc[N-1];
    
    A = new (nothrow) float[iLoc[N-1]];
	
    if (A == 0) {
	opserr << "FATAL:BandSPDLinSOE::BandSPDLinSOE :";
	opserr << " ran out of memory for A (profileSize) (";
	opserr << size <<", " << profileSize << ") \n";
	size = 0; profileSize = 0;
    }
    else {
	// zero the matrix
	Asize = iLoc[N-1];
	for (int k=0; k<Asize; k++)
	    A[k] = 0;
    
	B = new (nothrow) float[size];
	X = new (nothrow) float[size];
	doubleB = new (nothrow) double[size];
	doubleX = new (nothrow) double[size];
	iDiagLoc = new (nothrow) int[size];
    
	if (B == 0 || X == 0 || iDiagLoc == 0 || doubleX == 0 || doubleB == 0) {
	    opserr << "WARNING SProfileSPDLinSOE::SProfileSPDLinSOE :";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << size << ") \n";
	    size = 0; Bsize = 0;
	} else
	    Bsize = size;

	// zero the vectors
	for (int i=0; i<size; i++) {
	    B[i] = 0;
	    X[i] = 0;
	    doubleB[i] = 0;
	    doubleX[i] = 0;
	    iDiagLoc[i] = iLoc[i];
	}
    }

    vectX = new Vector(doubleX,size);
    vectB = new Vector(doubleB,size);
    
    the_Solver.setLinearSOE(*this);        
    
    int solverOK = the_Solver.setSize();
    if (solverOK < 0) {
	opserr << "WARNING SProfileSPDLinSOE::SProfileSPDLinSOE :";
	opserr << " solver failed setSize() in constructor\n";
    }
}

    
SProfileSPDLinSOE::~SProfileSPDLinSOE()
{
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (B != 0) delete [] doubleB;
    if (X != 0) delete [] doubleX;
    if (iDiagLoc != 0) delete [] iDiagLoc;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;    
}


int 
SProfileSPDLinSOE::getNumEqn(void) const
{
    return size;
}

int 
SProfileSPDLinSOE::setSize(Graph &theGraph)
{
    int oldSize = size;
    int result = 0;
    size = theGraph.getNumVertex();

    // check we have enough space in iDiagLoc and iLastCol
    // if not delete old and create new
    if (size > Bsize) { 
	if (iDiagLoc != 0) delete [] iDiagLoc;
	iDiagLoc = new (nothrow) int[size];

	if (iDiagLoc == 0) {
	    opserr << "WARNING SProfileSPDLinSOE::setSize() : ";
	    opserr << " - ran out of memory for iDiagLoc\n";
	    size = 0; Asize = 0;
	    result = -1;
	}
    }

    // zero out iDiagLoc 
    for (int i=0; i<size; i++) {
	iDiagLoc[i] = 0;
    }

    // now we go through the vertices to find the height of each col and
    // width of each row from the connectivity information.
    
    Vertex *vertexPtr;
    VertexIter &theVertices = theGraph.getVertices();

    while ((vertexPtr = theVertices()) != 0) {
	int vertexNum = vertexPtr->getTag();
	const ID &theAdjacency = vertexPtr->getAdjacency();
	int iiDiagLoc = iDiagLoc[vertexNum];
	int *iiDiagLocPtr = &(iDiagLoc[vertexNum]);

	for (int i=0; i<theAdjacency.Size(); i++) {
	    int otherNum = theAdjacency(i);
	    int diff = vertexNum-otherNum;
	    if (diff > 0) {
		if (iiDiagLoc < diff) {
		    iiDiagLoc = diff;
		    *iiDiagLocPtr = diff;
		}
	    } 
	}
    }


    // now go through iDiagLoc, adding 1 for the diagonal element
    // and then adding previous entry to give current location.
    if (iDiagLoc != 0)
	iDiagLoc[0] = 1; // NOTE FORTRAN ARRAY LOCATION

    for (int j=1; j<size; j++)
	iDiagLoc[j] = iDiagLoc[j] + 1 + iDiagLoc[j-1];

    if (iDiagLoc != 0)       
    	profileSize = iDiagLoc[size-1];

    // check if we need more space to hold A
    // if so then go get it
    if (profileSize > Asize) { 

	// delete old space
	if (A != 0)
	    delete [] A;
	
	// get new space
	A = new (nothrow) float[profileSize];
	
        if (A == 0) {
            opserr << "SProfileSPDLinSOE::SProfileSPDLinSOE :";
	    opserr << " ran out of memory for A (size,Profile) (";
	    opserr << size <<", " << profileSize << ") \n";
	    size = 0;  Asize = 0;  profileSize = 0;
	    result = -1;
        }
	else 
	    Asize = profileSize;
    }

    // zero the matrix
    for (int k=0; k<profileSize; k++)
	A[k] = 0;

    isAfactored = false;
    isAcondensed = false;    

    if (size > Bsize) { // we have to get another space for A
	
	// delete the old
	if (B != 0) delete [] B;
	if (X != 0) delete [] X;
	if (doubleB != 0) delete [] doubleB;
	if (doubleX != 0) delete [] doubleX;

	// create the new	
	B = new (nothrow)  float[size];
	X = new (nothrow) float[size];
	doubleB = new (nothrow) double[size];
	doubleX = new (nothrow) double[size];
	
        if (B == 0 || X == 0 || doubleB == 0 || doubleX == 0) {
            opserr << "SProfileSPDLinSOE::SProfileSPDLinSOE :";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << size << ") \n";
	    size = 0; Bsize = 0;
	    result = -1;
        }
    }

    // zero the vectors
    for (int l=0; l<size; l++) {
	B[l] = 0;
	X[l] = 0;
	doubleB[l] = 0;
	doubleX[l] = 0;
    }
    
    if (size != oldSize) {
	
	if (vectX != 0)
	    delete vectX;
	if (vectB != 0)
	    delete vectB;
	
	vectX = new Vector(doubleX,size);
	vectB = new Vector(doubleB,size);
	
	if (size > Bsize)
	  Bsize = size;
    }
    
    // invoke setSize() on the Solver
    LinearSOESolver *the_Solver = this->getSolver();
    int solverOK = the_Solver->setSize();
    if (solverOK < 0) {
	opserr << "WARNING SProfileSPDLinSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }    

    return result;
}

int 
SProfileSPDLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;
    
    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "SProfileSPDLinSOE::addA()	- Matrix and ID not of similar sizes\n";
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply 
      for (int i=0; i<idSize; i++) {
	int col = id(i);
	if (col < size && col >= 0) {
	  float *coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 
	  int minColRow;
	  if (col == 0)
	    minColRow = 0;
	  else
	    minColRow = col - (iDiagLoc[col] - iDiagLoc[col-1]) +1;
	  for (int j=0; j<idSize; j++) {
	    int row = id(j);
	    if (row <size && row >= 0 && 
		row <= col && row >= minColRow) { 
	      
	      // we only add upper and inside profile
	      float *APtr = coliiPtr + (row-col);
	      *APtr += m(j,i);
	    }
	  }  // for j
	} 
      }  // for i
    } else {
      for (int i=0; i<idSize; i++) {
	int col = id(i);
	if (col < size && col >= 0) {
	  float *coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 		
	  int minColRow;
	  if (col == 0)
	    minColRow = 0;
	  else
	    minColRow = col - (iDiagLoc[col] - iDiagLoc[col-1]) +1;
	  for (int j=0; j<idSize; j++) {
	    int row = id(j);
	    if (row <size && row >= 0 && 
		row <= col && row >= minColRow) { 
	      
	      // we only add upper and inside profile
	      float *APtr = coliiPtr + (row-col);
	      *APtr += m(j,i) * fact;
	    }
	  }  // for j
	} 
      }  // for i
    }
    return 0;
}

    
int 
SProfileSPDLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    
    // check for a quick return 
    if (fact == 0.0)  return 0;

    // check that m and id are of similar size
    int idSize = id.Size();        
    if (idSize != v.Size() ) {
	opserr << "SProfileSPDLinSOE::addB() -";
	opserr << " Vector and ID not of similar sizes\n";
	return -1;
    }    
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i) * fact;
	}
    }	
    return 0;
}


int
SProfileSPDLinSOE::setB(const Vector &v, double fact)
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
SProfileSPDLinSOE::zeroA(void)
{
    float *Aptr = A;
    for (int i=0; i<Asize; i++)
	*Aptr++ = 0;
    
    isAfactored = false;
}
	
void 
SProfileSPDLinSOE::zeroB(void)
{
    float *Bptr = B;
    for (int i=0; i<size; i++)
	*Bptr++ = 0;
}


void 
SProfileSPDLinSOE::setX(int loc, double value)
{
    if (loc < size && loc >=0)
	X[loc] = value;
}

void 
SProfileSPDLinSOE::setX(const Vector &x)
{
  if (x.Size() == size && vectX != 0)
    *vectX = x;
}

const Vector &
SProfileSPDLinSOE::getX(void)
{
  
  if (vectX == 0) {
    opserr << "FATAL SProfileSPDLinSOE::getX - vectX == 0";
    exit(-1);
  }    
  return *vectX;
}

const Vector &
SProfileSPDLinSOE::getB(void)
{
  if (vectB == 0) {
    opserr << "FATAL SProfileSPDLinSOE::getB - vectB == 0";
    exit(-1);
  }        
  return *vectB;
}

double 
SProfileSPDLinSOE::normRHS(void)
{
    double norm =0.0;
    for (int i=0; i<size; i++) {
	double Yi = B[i];
	norm += Yi*Yi;
    }
    return sqrt(norm);
    
}    


int
SProfileSPDLinSOE::setProfileSPDSolver(SProfileSPDLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);

    if (size != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:SProfileSPDLinSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return -1;
	}
    }
    
    return this->setSolver(newSolver);
}


int 
SProfileSPDLinSOE::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}


int 
SProfileSPDLinSOE::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}
