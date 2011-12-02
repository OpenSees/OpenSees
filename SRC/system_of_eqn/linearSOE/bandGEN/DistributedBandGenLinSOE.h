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
// $Date: 2009-05-11 20:57:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandGEN/DistributedBandGenLinSOE.h,v $
                                                                        
#ifndef DistributedBandGenLinSOE_h
#define DistributedBandGenLinSOE_h

// Written: fmk 
// Description: This file contains the class definition for DistributedBandGenLinSOE
// DistributedBandGenLinSOE is a subclass of LinearSOE. It uses the LAPACK storage
// scheme to store the components of the A matrix, which is a banded 
// unsymmetric matrix.
//
// What: "@(#) DistributedBandGenLinSOE.h, revA"


#include <BandGenLinSOE.h>
#include <Vector.h>

class DistributedBandGenLinSolver;

class DistributedBandGenLinSOE : public BandGenLinSOE
{
  public:
    DistributedBandGenLinSOE(BandGenLinSolver &theSolver);
    DistributedBandGenLinSOE();

    ~DistributedBandGenLinSOE();

    // these methods need to be rewritten
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);
    int setB(const Vector &, double fact = 1.0);            
    void zeroB(void);
    const Vector &getB(void);
    int solve(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    friend class BandGenLinLapackSolver;

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

