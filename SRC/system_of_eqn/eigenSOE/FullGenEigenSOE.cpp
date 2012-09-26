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
// $Date: 2009-05-20 17:31:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/FullGenEigenSOE.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 11/07
// Revision: A
//
// Description: This file contains the implementation of the
// FullGenEigenSOE class.

#include <FullGenEigenSOE.h>
#include <FullGenEigenSolver.h>
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

FullGenEigenSOE::FullGenEigenSOE(FullGenEigenSolver &theSolver,  
    AnalysisModel &aModel)
    : EigenSOE(theSolver, EigenSOE_TAGS_FullGenEigenSOE),
    size(0), A(0), Asize(0), M(0), Msize(0),
    factored(false), theModel(&aModel)
{
    theSolver.setEigenSOE(*this);
}


FullGenEigenSOE::~FullGenEigenSOE()
{
    if (A != 0)
        delete [] A;
    if (M != 0)
        delete [] M;
}


int FullGenEigenSOE::getNumEqn(void) const
{
    return size;
}


int FullGenEigenSOE::setSize(Graph &theGraph)
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
            opserr << "WARNING FullGenEigenSOE::setSize() - "
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
            opserr << "WARNING FullGenEigenSOE::setSize() - "
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
        opserr << "WARNING FullGenEigenSOE::setSize() - ";
        opserr << "solver failed in setSize()\n";
        return solverOK;
    } 

    return result;
}


int FullGenEigenSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for quick return 
    if (fact == 0.0)
        return 0;

    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
        opserr << "FullGenEigenSOE::addA() - Matrix and ID not of similar sizes\n";
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
                        *APtr += m(j,i)*fact;
                    }
                }  // for j
            } 
        }  // for i
    }

    return 0;
}


int FullGenEigenSOE::addM(const Matrix &m, const ID &id, double fact)
{
    // check for quick return 
    if (fact == 0.0)
        return 0;

    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
        opserr << "FullGenEigenSOE::addM() - Matrix and ID not of similar sizes\n";
        return -1;
    }

    if (fact == 1.0) { // do not need to multiply 
        for (int i=0; i<idSize; i++) {
            int col = id(i);
            if (col < size && col >= 0) {
                double *startColiPtr = M + col*size;
                for (int j=0; j<idSize; j++) {
                    int row = id(j);
                    if (row <size && row >= 0) {
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
                    if (row <size && row >= 0) {
                        double *MPtr = startColiPtr + row;
                        *MPtr += m(j,i)*fact;
                    }
                }  // for j
            } 
        }  // for i
    }

    return 0;
}


void FullGenEigenSOE::zeroA(void)
{
    double *Aptr = A;
    for (int i = 0; i < Asize; i++)
        *Aptr++ = 0;

    factored = false;
}


void FullGenEigenSOE::zeroM(void)
{
    double *Mptr = M;
    for (int i = 0; i < Msize; i++)
        *Mptr++ = 0;

    factored = false;
}


int FullGenEigenSOE::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}


int FullGenEigenSOE::recvSelf(int commitTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    return 0;
}
