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
                                                                        
// $Revision: 1.5 $
// $Date: 2009-05-11 20:56:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.h,v $
                                                                        
#ifndef MumpsParallelSOE_h
#define MumpsParallelSOE_h

// Written: fmk 
// Description: This file contains the class definition for MumpsParallelSOE
// MumpsParallelSOE is a subclass of LinearSOE. It uses the sparse column
// storage scheme. The matrix A is kept distributed, X and B kept on all processors.
//
// matrix types (matType): 0 Unsymmetrc
//                         1 Symmetrix positive definite
//                         2 General Symmetric

// What: "@(#) MumpsParallelSOE.h, revA"


#include <mpi.h>

#include <MumpsSOE.h>
#include <Vector.h>

class MumpsParallelSolver;

class MumpsParallelSOE : public MumpsSOE
{
  public:
    MumpsParallelSOE(MumpsParallelSolver &theSolver, int matType=0);
    MumpsParallelSOE();
    
    ~MumpsParallelSOE();

    // these methods need to be rewritten
    int setSize(Graph &theGraph);

    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);            
    const Vector &getB(void);
    void zeroB(void);
    int solve(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    friend class MumpsParallelSolver;        

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

