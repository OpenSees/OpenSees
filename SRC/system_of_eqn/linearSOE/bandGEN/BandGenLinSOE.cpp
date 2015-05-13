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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandGEN/BandGenLinSOE.cpp,v $
                                                                        
                                                                        
// Written: fmk 
//
// Description: This file contains the implementation for BandGenLinSOE

#include <stdlib.h>

#include <BandGenLinSOE.h>
#include <BandGenLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <iostream>
using std::nothrow;

BandGenLinSOE::BandGenLinSOE(BandGenLinSolver &theSolvr)
:LinearSOE(theSolvr, LinSOE_TAGS_BandGenLinSOE),
 size(0), numSuperD(0), numSubD(0), A(0), B(0), X(0), 
 vectX(0), vectB(0), Asize(0), Bsize(0), factored(false)
{
    theSolvr.setLinearSOE(*this);
}

BandGenLinSOE::BandGenLinSOE()
:LinearSOE(LinSOE_TAGS_BandGenLinSOE),
 size(0), numSuperD(0), numSubD(0), A(0), B(0), X(0), 
 vectX(0), vectB(0), Asize(0), Bsize(0), factored(false)
{

}

BandGenLinSOE::BandGenLinSOE(int classTag)
:LinearSOE(classTag),
 size(0), numSuperD(0), numSubD(0), A(0), B(0), X(0), 
 vectX(0), vectB(0), Asize(0), Bsize(0), factored(false)
{

}


BandGenLinSOE::BandGenLinSOE(int N, int numSuperDiag, int numSubDiag,
			     BandGenLinSolver &theSolvr)
:LinearSOE(theSolvr, LinSOE_TAGS_BandGenLinSOE),
 size(N), numSuperD(numSuperDiag), numSubD(numSubDiag), A(0), B(0), 
 X(0), vectX(0), vectB(0), Asize(0), Bsize(0), factored(false)
{
    Asize = N * (2*numSubD + numSuperD +1);
    A = new (nothrow)double[Asize];
	
    if (A == 0) {
	opserr << "WARNING BandGenLinSOE::BandGenLinSOE :";
	opserr << " ran out of memory for A (size,super,sub) (";
	opserr << size <<", " << numSuperDiag << ", " << numSubDiag << ") \n";
	Asize = 0; size = 0; numSubD = 0; numSuperD = 0; 
    } else {
	// zero the matrix
	for (int i=0; i<Asize; i++)
	    A[i] = 0;
    
	B = new (nothrow)double[size];
	X = new (nothrow)double[size];
	
	if (B == 0 || X == 0) {
	    opserr << "WARNING BandGenLinSOE::BandGenLinSOE :";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << size << ") \n";
	    Bsize = 0; size = 0; numSubD = 0; numSuperD = 0;
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
    
    int solverOK = theSolvr.setSize();
    if (solverOK < 0) {
	opserr << "WARNING BandGenLinSOE::BandGenLinSOE :";
	opserr << " solver failed setSize() in constructor\n";
    }    
}

int
BandGenLinSOE::getNumEqn(void) const
{
    return size;
}
    
BandGenLinSOE::~BandGenLinSOE()
{
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;    
}



int 
BandGenLinSOE::setSize(Graph &theGraph)
{
    int result = 0;
    int oldSize = size;
    size = theGraph.getNumVertex();
    
    /*
     * determine the number of superdiagonals and subdiagonals
     */
    
    numSubD = 0;
    numSuperD = 0;

    Vertex *vertexPtr;
    VertexIter &theVertices = theGraph.getVertices();
    
    while ((vertexPtr = theVertices()) != 0) {
	int vertexNum = vertexPtr->getTag();
	const ID &theAdjacency = vertexPtr->getAdjacency();
	for (int i=0; i<theAdjacency.Size(); i++) {
	    int otherNum = theAdjacency(i);
	    int diff = vertexNum - otherNum;
	    if (diff > 0) {
		if (diff > numSuperD)
		    numSuperD = diff;
	    } else 
		if (diff < numSubD)
		    numSubD = diff;
	}
    }
    numSubD *= -1;

    int newSize = size * (2*numSubD + numSuperD +1);
    if (newSize > Asize) { // we have to get another space for A

	if (A != 0) 
	    delete [] A;

	A = new (nothrow)double[newSize];
	
        if (A == 0) {
            opserr << "WARNING BandGenLinSOE::BandGenLinSOE :";
	    opserr << " ran out of memory for A (size,super,sub) (";
	    opserr << size <<", " << numSuperD << ", " << numSubD << ") \n";
	    Asize = 0; size = 0; numSubD = 0; numSuperD = 0;
	    result= -1;
        }
	else  
	    Asize = newSize;
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
	B = new (nothrow) double[size];
	X = new (nothrow) double[size];
	
        if (B == 0 || X == 0) {
            opserr << "WARNING BandGenLinSOE::BandGenLinSOE :";
	    opserr << " ran out of memory for vectors (size) (";
	    opserr << size << ") \n";
	    Bsize = 0; size = 0; numSubD = 0; numSuperD = 0;
	    result = -1;
        }
	else 
	    Bsize = size;
    }

    // zero the vectors
    for (int j=0; j<size; j++) {
	B[j] = 0;
	X[j] = 0;
    }

    // get new Vector objects if size has changes
    if (oldSize != size) {
	if (vectX != 0) 
	    delete vectX;

	if (vectB != 0) 
	    delete vectB;
		
	vectX = new Vector(X,size);
	vectB = new Vector(B,size);
    }
    
    // invoke setSize() on the Solver
    LinearSOESolver *theSolvr = this->getSolver();
    int solverOK = theSolvr->setSize();
    if (solverOK < 0) {
	opserr << "WARNING:BandGenLinSOE::setSize :";
	opserr << " solver failed setSize()\n";
	return solverOK;
    }    

    return result;    
}

int 
BandGenLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
  //  opserr << "BAND- addA() " << m << id;
    
    // check for a quick return 
    if (fact == 0.0)  return 0;

    
    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "BandGenLinSOE::addA()	- Matrix and ID not of similar sizes\n";
	return -1;
    }
    
    int ldA = 2*numSubD + numSuperD + 1;


    if (fact == 1.0) { // do not need to multiply 
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *coliiPtr = A + col*ldA + numSubD + numSuperD;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0) {		    
			int diff = col - row;
			if (diff > 0) {
			    if (diff <= numSuperD) {
				double *APtr = coliiPtr - diff;
				*APtr += m(j,i);
			    }			

			} else {
			    diff *= -1;
			    if (diff <= numSubD) {
				double *APtr = coliiPtr + diff;
				*APtr += m(j,i);
			    }
			}
		    }
		}  // for j
	    } 
	}  // for i
    } else {
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *coliiPtr = A + col*ldA + numSubD + numSuperD;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0) {		    
			int diff = col - row;
			if (diff > 0) {
			    if (diff <= numSuperD) {
				double *APtr = coliiPtr - diff;
				*APtr += m(j,i) *fact;
			    }
			} else {
			    diff *= -1;
			    if (diff <= numSubD) {
				double *APtr = coliiPtr + diff;
				*APtr += m(j,i) *fact;
			    }
			}
		    }
		}  // for j
	    } 
	}  // for i
    }    

    //    for (int i=0; i<Asize; i++) opserr << A[i] << " ";
    // opserr << endln;

    return 0;
}



int 
BandGenLinSOE::addColA(const Vector &colData, int col, double fact)
{
  if (fact == 0.0)  return 0;
  
  if (colData.Size() != size) {
    opserr << "BandGenLinSOE::addColA() - colData size not equal to n\n";
    return -1;
  }

  if (col > size && col < 0) {
    opserr << "BandGenLinSOE::addColA() - col " << col << "outside range 0 to " << size << endln;
    return -1;
  }
    
  int ldA = 2*numSubD + numSuperD + 1;
  
  if (fact == 1.0) { // do not need to multiply 

    double *coliiPtr = A + col*ldA + numSubD + numSuperD;
    for (int row=0; row<size; row++) {
      if (row <size && row >= 0) {		    
	int diff = col - row;
	if (diff > 0) {
	  if (diff <= numSuperD) {
	    double *APtr = coliiPtr - diff;
	    *APtr += colData(row);
	  }			
	} else {
	  diff *= -1;
	  if (diff <= numSubD) {
	    double *APtr = coliiPtr + diff;
	    *APtr += colData(row);
	  }
	}
      }
    }  // for j
  } else {

    double *coliiPtr = A + col*ldA + numSubD + numSuperD;
    for (int row=0; row<size; row++) {
      if (row <size && row >= 0) {		    
	int diff = col - row;
	if (diff > 0) {
	  if (diff <= numSuperD) {
	    double *APtr = coliiPtr - diff;
	    *APtr += colData(row);
	  }			
	} else {
	  diff *= -1;
	  if (diff <= numSubD) {
	    double *APtr = coliiPtr + diff;
	    *APtr += colData(row) * fact;
	  }
	}
      }
    }  
  }
  return 0;
}

    
int 
BandGenLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    // check that m and id are of similar size
    int idSize = id.Size();        
    if (idSize != v.Size() ) {
	opserr << "BandGenLinSOE::addB()	- Vector and ID not of similar sizes\n";
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
BandGenLinSOE::setB(const Vector &v, double fact)
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
BandGenLinSOE::zeroA(void)
{
    double *Aptr = A;
    int theSize = Asize;
    for (int i=0; i<theSize; i++)
	*Aptr++ = 0;
    
    factored = false;
}
	
void 
BandGenLinSOE::zeroB(void)
{
    double *Bptr = B;
    for (int i=0; i<size; i++)
	*Bptr++ = 0;
}


const Vector &
BandGenLinSOE::getX(void)
{
    if (vectX == 0) {
	opserr << "FATAL BandGenLinSOE::getX - vectX == 0!";
	exit(-1);
    }    
    
    return *vectX;
}


const Vector &
BandGenLinSOE::getB(void)
{
    if (vectB == 0) {
	opserr << "FATAL BandGenLinSOE::getB - vectB == 0!";
	exit(-1);
    }    

    return *vectB;
}


double 
BandGenLinSOE::normRHS(void)
{
    double norm =0.0;
    double *Bptr = B;
    for (int i=0; i<size; i++) {
	double Yi = *Bptr++;
	norm += Yi*Yi;
    }
    return sqrt(norm);
}    


void 
BandGenLinSOE::setX(int loc, double value)
{
    if (loc < size && loc >= 0)
	X[loc] = value;
}

void 
BandGenLinSOE::setX(const Vector &x)
{
    if (x.Size() == size && vectX != 0)
      *vectX = x;
}


int
BandGenLinSOE::setBandGenSolver(BandGenLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);
    
    if (size != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:BandGenLinSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return solverOK;
	}
    }	
    
    return this->setSolver(newSolver);
}


int 
BandGenLinSOE::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}


int 
BandGenLinSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}


