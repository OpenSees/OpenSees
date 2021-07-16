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
// $Date: 2009-05-11 21:08:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/ArpackSolver.h,v $

// Written: fmk
// Created: 05/09

// This is the solver that works on the ArpackSOE. It uses the LinearSOE
// in the SOE to perform the solve() operation if required.
// It uses the ARPACK library to perform the eigenvalue analysis.
// ARPACK is an eigen analysis package which was developed by 
// R.B.Lehoucq, D.C.Sorensen and C.Yang at Rice University. ARPACK is a
// collection of FORTRAN77 subroutines designed to solve large scale eigen
// problems. ARPACK is capable of solving large scale non-Hermitian standard 
// and generalized eigen problems. When the matrix <B>K</B> is symmetric, 
// the method is a variant of the Lanczos process called Implicitly Restarted
// Lanczos Method (IRLM).

//
// It is based on previous work of Jun Peng(Stanford)
//


#ifndef ArpackSolver_h
#define ArpackSolver_h

#include <EigenSolver.h>
#include <ArpackSOE.h>
#include <f2c.h>

class LinearSOE;

class ArpackSolver : public EigenSolver
{
  public:
    ArpackSolver();    
    ~ArpackSolver();

    int solve(int numMode, bool generalized, bool findSmallest = true);
    int setSize(void);
    int setEigenSOE(ArpackSOE &theSOE);
    
    const Vector &getEigenvector(int mode);
    double getEigenvalue(int mode);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
  protected:
    
  private:
    LinearSOE *theSOE;
    ArpackSOE *theArpackSOE;
    int numModesMax;
    int numMode;
    int size;
    double *eigenvalues;
    double *eigenvectors;
    Vector theVector;

    double shift;
    int ncv;
    double *v;
    double *workl;
    double *workd;
    double *resid;
    int iparam[11];
    int ipntr[11];
  //	long int* select;
  logical* select;
    
    void myMv(int n, double *v, double *result);
    void myCopy(int n, double *v, double *result);
    int getNCV(int n, int nev);
};

#endif


