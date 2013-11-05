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

// $Revision: 1.6 $
// $Date: 2009-05-19 21:54:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/FullGenEigenSolver.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 11/07
// Revision: A
//
// Description: This file contains the implementation of the
// FullGenEigenSolver class.

#include <FullGenEigenSolver.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <AnalysisModel.h>
#include <DOF_GrpIter.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <Integrator.h>


#ifdef _WIN32

extern "C" int DGGEV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                     double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                     double *BETA, double *VL, int *LDVL, double *VR,
                     int *LDVR, double *WORK, int *LWORK, int *INFO);

#else

extern "C" int dggev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                      double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                      double *BETA, double *VL, int *LDVL, double *VR,
                      int *LDVR, double *WORK, int *LWORK, int *INFO);

#endif


FullGenEigenSolver::FullGenEigenSolver()
    : EigenSolver(EigenSOLVER_TAGS_FullGenEigenSolver),
    theSOE(0), numEigen(0), eigenvalue(0),
    eigenvector(0), sortingID(0), eigenV(0)
{

}


FullGenEigenSolver::~FullGenEigenSolver()
{
    if (eigenvalue != 0)
        delete [] eigenvalue;
    if (eigenvector != 0)
        delete [] eigenvector;
    if (sortingID != 0)
        delete [] sortingID;
    if (eigenV != 0)
        delete eigenV;
}


int FullGenEigenSolver::solve(int nEigen, bool generalized, bool findSmallest)
{
  if (generalized == false) {
    opserr << "FullGenEigenSolver::solve() - only solves generalized problem\n";
    return -1;
  }
  
  if (theSOE == 0) {
    opserr << "FullGenEigenSolver::solve()- "
	   << " No EigenSOE object has been set yet\n";
    return -1;
  }

    // check for quick return
    if (nEigen < 1) {
        numEigen = 0;
        return 0;
    }

    // get the number of equations
    int n = theSOE->size;

    // set the number of eigenvalues
    numEigen = nEigen;
    if (numEigen > n)
        numEigen = n;

    // do not compute left eigenvalues and eigenvectors
    char *jobvl = "N";

    // compute right eigenvalues and eigenvectors
    char *jobvr = "V";

    // stiffness matrix data
    double *Kptr = theSOE->A;

    // leading dimension of K
    int ldK = n;

    // mass matrix data
    double *Mptr = theSOE->M;

    // leading dimension of M
    int ldM = n;

    // allocate memory for eigenvalues
    double *alphaR = new double [n];
    double *alphaI = new double [n];
    double *beta   = new double [n];

    if (eigenvalue != 0)
        delete [] eigenvalue;

    eigenvalue = new double [n];

    // allocate memory for sorting index array
    if (sortingID != 0)
        delete [] sortingID;
    sortingID = new int [n];

    // dummy left eigenvectors
    double vl[1];

    // leading dimension of dummy left eigenvectors
    int ldvl = 1;

    // allocate memory for right eigenvectors
    if (eigenvector != 0)
        delete [] eigenvector;
    eigenvector = new double [n*n];

    // leading dimension of right eigenvectors
    int ldvr = n;

    // dimension of the workspace array
    int lwork = n*(8+64);

    // allocate memory for workspace array
    double *work = new double [lwork];

    // output information
    int info = 0;

    // call the LAPACK eigenvalue subroutine
#ifdef _WIN32
    DGGEV(jobvl, jobvr, &n, Kptr, &ldK, Mptr, &ldM, alphaR, alphaI, beta,
          vl, &ldvl, eigenvector, &ldvr, work, &lwork, &info);
#else
    dggev_(jobvl, jobvr, &n, Kptr, &ldK, Mptr, &ldM, alphaR, alphaI, beta,
           vl, &ldvl, eigenvector, &ldvr, work, &lwork, &info);
#endif

    if (info < 0) {
        opserr << "FullGenEigenSolver::solve() - invalid argument number "
            << -info << " passed to LAPACK dggev routine\n";
        return info;
    }

    if (info > 0) {
        opserr << "FullGenEigenSolver::solve() - the LAPACK dggev routine "
            << "returned error code " << info << endln;
        return -info;
    }

    theSOE->factored = true;

    for (int i=0; i<n; i++) {
        double mag = sqrt(alphaR[i]*alphaR[i]+alphaI[i]*alphaI[i]);
        if (mag*DBL_EPSILON < fabs(beta[i])) {
            if (alphaI[i] == 0.0) {
                eigenvalue[i] = alphaR[i]/beta[i];
            }
            else {
                eigenvalue[i] = -mag/beta[i];
                opserr << "FullGenEigenSolver::solve() - the eigenvalue "
                    << i+1 << " is complex with magnitude "
                    << -eigenvalue[i] << endln;
            }
        }
        else {
	    eigenvalue[i] = DBL_MAX;
	}
        sortingID[i] = i;
    }


    // sort eigenvalues in ascending order and return sorting ID
    this->sort(n, eigenvalue, sortingID);

    for (int i=0; i<numEigen; i++) {
        if (eigenvalue[i] == DBL_MAX) {
	    opserr << "FullGenEigenSolver::solve() - the eigenvalue "
		    << i+1 << " is numerically undetermined or infinite\n";
        }
    }

    int lworkOpt = (int) work[0];
    if (lwork < lworkOpt) {
        opserr << "FullGenEigenSolver::solve() - optimal workspace size "
                << lworkOpt << " is larger than provided workspace size "
                << lwork << " consider increasing workspace\n";
    }

    // clean up the memory
    delete [] alphaR;
    delete [] alphaI;
    delete [] beta;
    delete [] work;

    return 0;
}


int FullGenEigenSolver::setSize()
{
    int size = theSOE->size;

    if (eigenV == 0 || eigenV->Size() != size) {
        if (eigenV != 0)
            delete eigenV;

        eigenV = new Vector(size);
        if (eigenV == 0 || eigenV->Size() != size) {
            opserr << "FullGenEigenSolver::setSize() ";
            opserr << " - ran out of memory for eigenVector of size ";
            opserr << theSOE->size << endln;
            return -2;	    
        }
    }

    return 0;
}


int FullGenEigenSolver::setEigenSOE(FullGenEigenSOE &thesoe)
{
    theSOE = &thesoe;

    return 0;
}


const Vector& FullGenEigenSolver::getEigenvector(int mode)
{
    if (mode <= 0 || mode > numEigen) {
        opserr << "FullGenEigenSolver::getEigenVector() - mode "
            << mode << " is out of range (1 - " << numEigen << ")\n";
        eigenV->Zero();
        return *eigenV;
    }

    int size = theSOE->size;
    int index = size*sortingID[mode-1];

    if (eigenvector != 0) {
        for (int i=0; i<size; i++) {
            (*eigenV)[i] = eigenvector[index++];
        }	
    }
    else {
        opserr << "FullGenEigenSolver::getEigenvector() - "
            << "eigenvectors not computed yet\n";
        eigenV->Zero();
    }      

    return *eigenV;
}


double FullGenEigenSolver::getEigenvalue(int mode)
{
    if (mode <= 0 || mode > numEigen) {
        opserr << "FullGenEigenSolver::getEigenvalue() - mode " 
            << mode << " is out of range (1 - " << numEigen << ")\n";
        return 0.0;
    }

    if (eigenvalue != 0) {
        return eigenvalue[mode-1];
    }
    else {
        opserr << "FullGenEigenSolver::getEigenvalue() - "
            << "eigenvalues not yet computed\n";
        return 0.0;
    }      
}


int FullGenEigenSolver::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}


int FullGenEigenSolver::recvSelf(int commitTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    return 0;
}


void FullGenEigenSolver::sort(int length, double *x, int *id)
{
    // this is an implementation of shell sort that
    // additionally keeps track of the sorting order
    int flag = 1;
    int d = length;
    int i, idTmp;
    double xTmp;
    
    while (flag || d>1) {
        flag = 0;
        d = (d+1)/2;
        for (i=0; i<(length-d); i++) {
            if (x[i+d] < x[i]) {
                // swap items at positions i+d and d
	            xTmp = x[i+d];  idTmp = id[i+d]; 
	            x[i+d] = x[i];  id[i+d] = id[i]; 
	            x[i] = xTmp;    id[i] = idTmp; 
	            // indicate that a swap has occurred
	            flag = 1;
            }
        }
    }

    return;
}
