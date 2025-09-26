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

#ifndef SymmGenEigenSolver_h
#define SymmGenEigenSolver_h

#include <EigenSolver.h>
#include <SymmGeneralizedEigenSOE.h>

class SymmGeneralizedEigenSolver : public EigenSolver
{
public:
    SymmGeneralizedEigenSolver(double msmall = 1e-10);    
    virtual ~SymmGeneralizedEigenSolver();

    virtual int solve(int numEigen, bool generalized, bool findSmallest = true);    
    virtual int setSize(void);
    virtual int setEigenSOE(SymmGeneralizedEigenSOE &theSOE);

    virtual const Vector &getEigenvector(int mode);
    virtual double getEigenvalue(int mode);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

protected:

private:
    void sort(int length, double *x, int *id);
    SymmGeneralizedEigenSOE *theSOE;
    int numEigen;

    double *eigenvalue;
    double *eigenvector;
    int *sortingID;
    Vector *eigenV;

  double msmall; // m_ii = k_ii*msmall if m_ii == 0.0
};

#endif
