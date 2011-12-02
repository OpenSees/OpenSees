// File: ~/system_of_eqn/eigenSOE/BandArpackSOE.C
//
// Written: Jun Peng 
// Created: Febuary 1999
// Revision: A
//
// Description: This file contains the class definition for BandArpackSOE
// BandArpackSOE is a subclass of EigenSOE. It uses the LAPACK storage
// scheme to store the components of the K, M matrix, which is a full matrix.
// It uses the ARPACK to do eigenvalue analysis.

#include <BandArpackSOE.h>
#include <BandArpackSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <AnalysisModel.h>

BandArpackSOE::BandArpackSOE(BandArpackSolver &theSolvr,  
			     AnalysisModel &aModel, double theShift)
:EigenSOE(theSolvr, EigenSOE_TAGS_BandArpackSOE),
 size(0), numSubD(0), numSuperD(0), A(0), 
 Asize(0), factored(false), shift(theShift), theModel(&aModel)
{
    theSolvr.setEigenSOE(*this);
}


int
BandArpackSOE::getNumEqn(void) const
{
    return size;
}
    
BandArpackSOE::~BandArpackSOE()
{
    if (A != 0) delete [] A;
}

int 
BandArpackSOE::setSize(Graph &theGraph)
{
    int result = 0;
    size = theGraph.getNumVertex();
        
    // determine the number of superdiagonals and subdiagonals
    
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

	A = new double[newSize];
	
        if (A == 0) {
            cerr << "WARNING BandGenLinSOE::BandGenLinSOE :";
	    cerr << " ran out of memory for A (size,super,sub) (";
	    cerr << size <<", " << numSuperD << ", " << numSubD << ") \n";
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

    // invoke setSize() on the Solver
    EigenSolver *theSolvr = this->getSolver();
    int solverOK = theSolvr->setSize();
    if (solverOK < 0) {
	cerr << "WARNING:BandArpackSOE::setSize :";
	cerr << " solver failed setSize()\n";
	return solverOK;
    } 
   
    return result;    
}

int 
BandArpackSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
	cerr << "BandArpackSOE::addA()	- Matrix and ID not of similar sizes\n";
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
    return 0;
}


void 
BandArpackSOE::zeroA(void)
{
    double *Aptr = A;
    int theSize = size*(2*numSubD+numSuperD+1);
    for (int i=0; i<theSize; i++)
	*Aptr++ = 0;
    
    factored = false;
}

int 
BandArpackSOE::addM(const Matrix &m, const ID &id, double fact)
{
  return this->addA(m, id, -shift);
}   
 
void 
BandArpackSOE::zeroM(void)
{
  
}


double 
BandArpackSOE::getShift(void)
{
    return shift;
}


int 
BandArpackSOE::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}
    
int 
BandArpackSOE::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
    return 0;
}

