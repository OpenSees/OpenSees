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
// $Date: 2009-05-11 20:58:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/diagonal/DistributedDiagonalSOE.h,v $
                                                                        
// Written: fmk 
// Created: 05/05
//
// Description: This file contains the class definition for DistributedDiagonalSOE
// DistributedDiagonalSOE is a subclass of LinearSOE. It stores a diagonal system
// of equation, i.e. just the diagonal

// What: "@(#) DistributedDiagonalSOE.h, revA"

#ifndef DistributedDiagonalSOE_h
#define DistributedDiagonalSOE_h

#include <LinearSOE.h>
#include <Vector.h>
#include <ID.h>


class AnalysisModel;
class DistributedDiagonalSolver;

class DistributedDiagonalSOE : public LinearSOE
{
  public:
    DistributedDiagonalSOE(DistributedDiagonalSolver &theSolver);
    DistributedDiagonalSOE();
    ~DistributedDiagonalSOE();

    int getNumEqn(void) const;
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        
    
    void zeroA(void);
    void zeroB(void);

    void setX(int loc, double value);
    void setX(const Vector &x);
    
    const Vector &getX(void);
    const Vector &getB(void);
    double normRHS(void);

    int setDiagonalSolver(DistributedDiagonalSolver &newSolver);    
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    int setChannels(int nChannels, Channel **theC);

    int setAnalysisModel(AnalysisModel &theModel);

    friend class DistributedDiagonalSolver;
    
  protected:
    
  private:
    int size;
    double *A, *B, *X;
    Vector *vectX;
    Vector *vectB;
    bool isAfactored;

    int processID;
    int numProcesses;
    int numChannels;
    Channel **theChannels;
    ID **localCol;

    ID myDOFs;
    ID myDOFsShared;
    int numShared;
    double *dataShared;
    Vector *vectShared;
    AnalysisModel *theModel;
};


#endif



