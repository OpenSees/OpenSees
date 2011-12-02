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
// $Date: 2009-05-11 20:57:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenColLinSOE.h,v $
                                                                        
#ifndef DistributedSparseGenColLinSOE_h
#define DistributedSparseGenColLinSOE_h

// Written: fmk 
// Description: This file contains the class definition for DistributedSparseGenColLinSOE
// DistributedSparseGenColLinSOE is a subclass of LinearSOE. It uses the sparse column
// storage scheme. It collects all contributions from the matrices and assembles them onto 
// P0 (not really a distributed SOE .. but until i get a solver that will work with a distributed
// col storage scheme this is what will have to stick with).
//
// What: "@(#) DistributedSparseGenColLinSOE.h, revA"


#include <SparseGenColLinSOE.h>
#include <Vector.h>

class DistributedSparseGenColLinSolver;

class DistributedSparseGenColLinSOE : public SparseGenColLinSOE
{
  public:
    DistributedSparseGenColLinSOE(SparseGenColLinSolver &theSolver);
    DistributedSparseGenColLinSOE();
    
    ~DistributedSparseGenColLinSOE();

    // these methods need to be rewritten
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);            
    const Vector &getB(void);
    void zeroB(void);
    int solve(void);


    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    friend class SuperLU;    
    friend class ThreadedSuperLU;        
    friend class DistributedSuperLU;        

    int setProcessID(int processTag);
    int setChannels(int numChannels, Channel **theChannels);

  protected:
    
  private:
    int processID;
    int numChannels;
    Channel **theChannels;
    ID **localCol;

    double *workArea;
    int sizeWork;
    double *myB;
    Vector *myVectB;
};


#endif

