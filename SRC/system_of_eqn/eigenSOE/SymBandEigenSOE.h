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
// $Date: 2004-10-05 00:17:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/SymBandEigenSOE.h,v $

// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition for
// SymBandEigenSOE, which stores a symmetric banded matrix, A, for
// standard eigenvalue computations.

#ifndef SymBandEigenSOE_h
#define SymBandEigenSOE_h

#include <EigenSOE.h>
#include <Vector.h>

class AnalysisModel;
class SymBandEigenSolver;

class SymBandEigenSOE : public EigenSOE
{
  public:
    SymBandEigenSOE(SymBandEigenSolver &theSolver, AnalysisModel &theModel);

    virtual ~SymBandEigenSOE();

    virtual int getNumEqn(void) const;
    virtual int setSize(Graph &theGraph);
    
    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual int addM(const Matrix &, const ID &, double fact = 1.0);    
   
    virtual void zeroA(void);
    virtual void zeroM(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
    friend class SymBandEigenSolver;

  protected:
    
  private:
    int size;
    int numSuperD;
    double *A;
    int Asize;
    bool factored;
    AnalysisModel *theModel;    

    double *M;  // handles diagonal mass matrix
    int Msize;
};

#endif
