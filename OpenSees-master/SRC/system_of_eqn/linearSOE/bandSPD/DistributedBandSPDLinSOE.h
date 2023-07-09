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
// $Date: 2010-02-25 23:21:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandSPD/DistributedBandSPDLinSOE.h,v $
                                                                        
#ifndef DistributedBandSPDLinSOE_h
#define DistributedBandSPDLinSOE_h

// Written: fmk 
// Description: This file contains the class definition for DistributedBandSPDLinSOE
// DistributedBandSPDLinSOE is a subclass of LinearSOE. It uses the LAPACK storage
// scheme to store the components of the A matrix, which is a banded 
// unsymmetric matrix.
//
// What: "@(#) DistributedBandSPDLinSOE.h, revA"


#include <BandSPDLinSOE.h>
#include <Vector.h>

class DistributedBandSPDLinSolver;

class DistributedBandSPDLinSOE : public BandSPDLinSOE
{
  public:
    DistributedBandSPDLinSOE(BandSPDLinSolver &theSolver);
    DistributedBandSPDLinSOE();
    
    ~DistributedBandSPDLinSOE();

    // these methods need to be rewritten
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);            
    void zeroB(void);
    int setSize(Graph &theGraph);
    int solve(void);
    const Vector &getB(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    friend class BandSPDLinLapackSolver;

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
    Vector *myVectB;
    double *myB;
};


#endif

