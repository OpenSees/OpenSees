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
                                                                        
// $Revision: 1.1 $
// $Date: 2005-12-06 22:20:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenRowLinSOE.h,v $
                                                                        
#ifndef DistributedSparseGenRowLinSOE_h
#define DistributedSparseGenColLinSOE_h

// Written: fmk 
// Created: 04/05
// Revision: A
//
// Description: This file contains the class definition for DistributedSparseGenColLinSOE
// DistributedSparseGenColLinSOE is a subclass of LinearSOE. It stores the matrix equation
// Ax=b using the sparse row-compacted storage scheme for storing the 
// matrix A. 
//
// What: "@(#) DistributedSparseGenColLinSOE.h, revA"

#include <LinearSOE.h>
#include <Vector.h>
#include <ID.h>

class DistributedSparseGenRowLinSolver;

class DistributedSparseGenRowLinSOE : public LinearSOE
{
  public:
    DistributedSparseGenRowLinSOE(DistributedSparseGenRowLinSolver &theSolver);        
    ~DistributedSparseGenRowLinSOE();

    int getNumEqn(void) const;
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        
    
    void zeroA(void);
    void zeroB(void);
    
    const Vector &getX(void);
    const Vector &getB(void);    
    double normRHS(void);

    void setX(int loc, double value);        
    void setX(const Vector &x);        
    int setDistributedSparseGenRowSolver(DistributedSparseGenRowLinSolver &newSolver);    

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    friend class DistributedSuperLU;        

    int setProcessID(int processTag);
    int setChannels(int numChannels, Channel **theChannels);

  protected:
    
  private:
    int processID;
    int numChannels;
    Channel **theChannels;
    ID **localCol;

    int n;                 // order of A
    int nnz;               // number of non-zeros in A
    double *A, *B, *X;     // 1d arrays containing coefficients of A, B and X
    int *colA, *rowStartA; // int arrays containing info about coeficientss in A
    Vector *vectX;
    Vector *vectB;    
    int Asize, Bsize;      // size of the 1d array holding A
    bool factored;

    int numP;                // numProcesses involved in computation
    int rank;                // rank of current process
    int startRow;            // start row number of rows assigned to current process
    int endRow;              // end row number of rows assigned to current
    int numRows;             // number of rows assigned to Process, 
                             //  first rows will be startRow through endRow; after that
                             //  come all rows not yet accounted for in rows that are in graph passed in setGraph
    ID rows;                 // the rows of A
    ID **otherProcessesRows; // the rows & row location of data that will be sent by the processes for rows
                             // not assigned to them yet for which they have data.
    ID **otherProcessesRowStart;
    int remoteDataLoc;
};


#endif

