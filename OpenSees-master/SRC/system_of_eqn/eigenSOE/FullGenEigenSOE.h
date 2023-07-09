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
// $Date: 2007-11-29 19:31:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/FullGenEigenSOE.h,v $


#ifndef FullGenEigenSOE_h
#define FullGenEigenSOE_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 11/07
// Revision: A
//
// Description: This file contains the class definition for
// FullGenEigenSOE, which stores full nonsymmetric matrices,
// A and M, for generalized eigenvalue computations.

#include <EigenSOE.h>
#include <Vector.h>

class AnalysisModel;
class FullGenEigenSolver;

class FullGenEigenSOE : public EigenSOE
{
public:
    FullGenEigenSOE(FullGenEigenSolver &theSolver,
        AnalysisModel &theModel);

    virtual ~FullGenEigenSOE();

    virtual int getNumEqn(void) const;
    virtual int setSize(Graph &theGraph);

    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual int addM(const Matrix &, const ID &, double fact = 1.0);    

    virtual void zeroA(void);
    virtual void zeroM(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);

    friend class FullGenEigenSolver;

protected:

private:
    int size;
    double *A;
    int Asize;
    double *M;
    int Msize;
    bool factored;
    AnalysisModel *theModel;
};

#endif
