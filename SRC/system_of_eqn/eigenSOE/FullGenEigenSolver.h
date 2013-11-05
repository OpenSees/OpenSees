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

// $Revision: 1.3 $
// $Date: 2009-05-11 21:01:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/FullGenEigenSolver.h,v $


#ifndef FullGenEigenSolver_h
#define FullGenEigenSolver_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 11/07
// Revision: A
//
// Description: This file contains the class definition for 
// FullGenEigenSolver. It computes the generalized eigenvalues
// and eigenvectors of a pair of real nonsymmetric matrices using
// the LAPACK subroutine DGGEV.

#include <EigenSolver.h>
#include <FullGenEigenSOE.h>

class FullGenEigenSolver : public EigenSolver
{
public:
    FullGenEigenSolver();    
    virtual ~FullGenEigenSolver();

    virtual int solve(int numEigen, bool generalized, bool findSmallest = true);    
    virtual int setSize(void);
    virtual int setEigenSOE(FullGenEigenSOE &theSOE);

    virtual const Vector &getEigenvector(int mode);
    virtual double getEigenvalue(int mode);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);

protected:

private:
    void sort(int length, double *x, int *id);

    FullGenEigenSOE *theSOE;
    int numEigen;

    double *eigenvalue;
    double *eigenvector;
    int *sortingID;
    Vector *eigenV;
};

#endif
