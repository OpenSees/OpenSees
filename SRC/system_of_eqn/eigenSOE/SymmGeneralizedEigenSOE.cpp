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

#include <SymmGeneralizedEigenSOE.h>
#include <SymmGeneralizedEigenSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <AnalysisModel.h>
#include <iostream>
using std::nothrow;

SymmGeneralizedEigenSOE::SymmGeneralizedEigenSOE(SymmGeneralizedEigenSolver &theSolver,  
    AnalysisModel &aModel)
    : EigenSOE(theSolver, EigenSOE_TAGS_SymmGeneralizedEigenSOE),
    size(0), A(0), Asize(0), M(0), Msize(0),
    factored(false), theModel(&aModel)
{
    theSolver.setEigenSOE(*this);
}


SymmGeneralizedEigenSOE::~SymmGeneralizedEigenSOE()
{
    if (A != 0)
        delete [] A;
    if (M != 0)
        delete [] M;
}


int SymmGeneralizedEigenSOE::getNumEqn(void) const
{
    return size;
}


int SymmGeneralizedEigenSOE::setSize(Graph &theGraph)
{
    int result = 0;
    size = theGraph.getNumVertex();

    int newSize = size*size;

    // we have to get another space for A
    if (newSize > Asize) {
        if (A != 0) 
            delete [] A;

        A = new (nothrow) double[newSize];
        if (A == 0) {
            opserr << "WARNING SymmGeneralizedEigenSOE::setSize() - "
                << "ran out of memory for A (size,size) ("
                << size << ", " << size << ")\n";
            Asize = 0; size = 0;
            result= -1;
        }
        else
            Asize = newSize;
    }

    // zero the matrix
    for (int i=0; i<Asize; i++)
        A[i] = 0;

    // we have to get another space for M
    if (newSize > Msize) {
        if (M != 0) 
            delete [] M;

        M = new (nothrow) double[newSize];
        if (M == 0) {
            opserr << "WARNING SymmGeneralizedEigenSOE::setSize() - "
                << "ran out of memory for M (size,size) ("
                << size << ", " << size << ")\n";
            Msize = 0; size = 0;
            result= -1;
        }
        else
            Msize = newSize;
    }

    // zero the matrix
    for (int i=0; i<Msize; i++)
        M[i] = 0;

    factored = false;

    // invoke setSize() on the Solver
    EigenSolver *theSolver = this->getSolver();
    int solverOK = theSolver->setSize();
    if (solverOK < 0) {
        opserr << "WARNING SymmGeneralizedEigenSOE::setSize() - ";
        opserr << "solver failed in setSize()\n";
        return solverOK;
    } 

    return result;
}


int SymmGeneralizedEigenSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for quick return 
    if (fact == 0.0)
        return 0;

    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
        opserr << "SymmGeneralizedEigenSOE::addA() - Matrix and ID not of similar sizes\n";
        return -1;
    }

    if (fact == 1.0) { // do not need to multiply 
        for (int i=0; i<idSize; i++) {
            int col = id(i);
            if (col < size && col >= 0) {
                double *startColiPtr = A + col*size;
                for (int j=0; j<idSize; j++) {
                    int row = id(j);
                    if (row <size && row >= 0 &&
			row <= col) { // Only add upper
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
                    if (row <size && row >= 0 &&
			row <= col) { // Only add upper
                        double *APtr = startColiPtr + row;
                        *APtr += m(j,i)*fact;
                    }
                }  // for j
            } 
        }  // for i
    }

    return 0;
}


int SymmGeneralizedEigenSOE::addM(const Matrix &m, const ID &id, double fact)
{
    // check for quick return 
    if (fact == 0.0)
        return 0;

    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
        opserr << "SymmGeneralizedEigenSOE::addM() - Matrix and ID not of similar sizes\n";
        return -1;
    }

    if (fact == 1.0) { // do not need to multiply 
        for (int i=0; i<idSize; i++) {
            int col = id(i);
            if (col < size && col >= 0) {
                double *startColiPtr = M + col*size;
                for (int j=0; j<idSize; j++) {
                    int row = id(j);
                    if (row <size && row >= 0 &&
			row <= col) { // Only add upper
                        double *MPtr = startColiPtr + row;
                        *MPtr += m(j,i);
                    }
                }  // for j
            } 
        }  // for i
    } else {
        for (int i=0; i<idSize; i++) {
            int col = id(i);
            if (col < size && col >= 0) {
                double *startColiPtr = M + col*size;
                for (int j=0; j<idSize; j++) {
                    int row = id(j);
                    if (row <size && row >= 0 &&
			row <= col) { // Only add upper
                        double *MPtr = startColiPtr + row;
                        *MPtr += m(j,i)*fact;
                    }
                }  // for j
            } 
        }  // for i
    }
    
    return 0;
}


void SymmGeneralizedEigenSOE::zeroA(void)
{
    double *Aptr = A;
    for (int i = 0; i < Asize; i++)
        *Aptr++ = 0;

    factored = false;
}


void SymmGeneralizedEigenSOE::zeroM(void)
{
    double *Mptr = M;
    for (int i = 0; i < Msize; i++)
        *Mptr++ = 0;

    factored = false;
}


int SymmGeneralizedEigenSOE::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}


int SymmGeneralizedEigenSOE::recvSelf(int commitTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    return 0;
}
