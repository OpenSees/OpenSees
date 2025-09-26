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

#include <SymmGeneralizedEigenSolver.h>
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

extern "C" int DSYGVX(int *ITPYE, char *JOBZ, char *RANGE, char *UPLO,
		      int *N, double *A, int *LDA, double *B, int *LDB,
		      double *VL, double *VU, int *IL, int *IU,
		      double *ABSTOL, int *M, double *W, double *Z, int *LDZ,
		      double *WORK, int *LWORK, int *IWORK, int *IFAIL,
		      int *INFO);

#else

extern "C" int dsygvx_(int *ITPYE, char *JOBZ, char *RANGE, char *UPLO,
		      int *N, double *A, int *LDA, double *B, int *LDB,
		      double *VL, double *VU, int *IL, int *IU,
		      double *ABSTOL, int *M, double *W, double *Z, int *LDZ,
		      double *WORK, int *LWORK, int *IWORK, int *IFAIL,
		      int *INFO);

#endif


SymmGeneralizedEigenSolver::SymmGeneralizedEigenSolver(double m)
    : EigenSolver(EigenSOLVER_TAGS_SymmGeneralizedEigenSolver),
    theSOE(0), numEigen(0), eigenvalue(0),
      eigenvector(0), sortingID(0), eigenV(0), msmall(m)
{

}


SymmGeneralizedEigenSolver::~SymmGeneralizedEigenSolver()
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


int SymmGeneralizedEigenSolver::solve(int nEigen, bool generalized, bool findSmallest)
{
  if (generalized == false) {
    opserr << "SymmGeneralizedEigenSolver::solve() - only solves generalized problem\n";
    return -1;
  }
  
  if (theSOE == 0) {
    opserr << "SymmGeneralizedEigenSolver::solve()- "
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

    // Lower and upper range of eigenpairs
    int il = 1;
    int iu = numEigen;
    
    // solve K*x = lam*M*x
    int itype = 1;
    
    // compute eigenvalues and eigenvectors
    char *jobz = "V";

    // compute range of eigenvalues
    char *range = "I";

    // upper or lower triangle of A and B
    char *uplo = "U";
    
    // stiffness matrix data
    double *Kptr = theSOE->A;
      
    double *kCopy = new double[n*n];
    for (int i = 0; i < n*n; i++)
      kCopy[i] = Kptr[i];
    
    // leading dimension of K
    int ldK = n;

    // mass matrix data
    double *Mptr = theSOE->M;

    // Check for zero mass on diagonal, add some small mass
    int index = 0;
    for (int i = 0; i < n; i++) {
      if (Mptr[index] == 0.0)
	Mptr[index] = Kptr[index]*msmall;
      index += n+1;
    }
    
    double *mCopy = new double[n*n];
    for (int i=0; i<n*n; i++)
      mCopy[i] = Mptr[i];
    
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

    // Not used
    double vl = 0.0;
    double vu = 0.0;

    // allocate memory for right eigenvectors
    if (eigenvector != 0)
        delete [] eigenvector;
    eigenvector = new double [numEigen*n];

    // number of eigenvalues found
    int m = 0;
    
    // 
    int ldz = n;
    double *w = eigenvalue;
    double *z = eigenvector;
    
    double abstol = -1.0;
    
    // dimension of the workspace array
    int lwork = n*(8+64);
    int liwork = n*5;

    // allocate memory for workspace array
    double *work = new double [lwork];
    int *iwork = new int [liwork];

    // fail array
    int *ifail = new int [n];
    
    // output information
    int info = 0;

    // call the LAPACK eigenvalue subroutine
#ifdef _WIN32
    // itype=1, jobz='V', range='I', uplo='U', n=N, Kptr=A, ldK=N, Mptr=B, ldM=N,
    // vl=N/A, vu=N/A, il=1, iu=nEigen, abstol=0
    // m=out (num eigenvalues found)
    // w=out (first m eigenvalues)
    // z=out (first m eigenvectors)
    // ldz=out (leading dimension of eigenvectors)
    // work=double array of length lwork
    // lwork=8*n
    // iwork=int array of length 5*n
    // ifail=out (0 if success, if info>0 ifail has indices of eigenvectors that failed to converge)
    // info=out (0 if success, <0 arg error, >0 failed to converge)
    DSYGVX(&itype, jobz, range, uplo, &n, Kptr, &ldK, Mptr, &ldM, 
	   &vl, &vu, &il, &iu, &abstol,
	   &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);
#else
    dsygvx_(&itype, jobz, range, uplo, &n, Kptr, &ldK, Mptr, &ldM, 
	    &vl, &vu, &il, &iu, &abstol,
	    &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);
#endif
    
    if (info < 0) {
        opserr << "SymmGeneralizedEigenSolver::solve() - invalid argument number "
            << -info << " passed to LAPACK dsygvx routine\n";
        return info;
    }

    if (info > 0) {
        opserr << "SymmGeneralizedEigenSolver::solve() - the LAPACK dsygvx routine "
            << "returned error code " << info << endln;
        return -info;
    }

    theSOE->factored = true;

    for (int i=0; i<n; i++) {
      /*
        double mag = sqrt(alphaR[i]*alphaR[i]+alphaI[i]*alphaI[i]);
        if (mag*DBL_EPSILON < fabs(beta[i])) {
            if (alphaI[i] == 0.0) {
                eigenvalue[i] = alphaR[i]/beta[i];
            }
            else {
                eigenvalue[i] = -mag/beta[i];
                opserr << "SymmGeneralizedEigenSolver::solve() - the eigenvalue "
                    << i+1 << " is complex with magnitude "
                    << -eigenvalue[i] << endln;
            }
        }
        else {
	    eigenvalue[i] = DBL_MAX;
	}
      */
        sortingID[i] = i;
    }

    //
    // mass normalize the eigenvalues
    //

    Kptr = kCopy;
    Mptr = mCopy;    
    double *tmpV = new double[n];
    /*
    // mass normailze all vectors .. NOTE instead of numEigen!
    for (int k=0; k<n; k++) {

      for (int i=0; i<n; i++)
	tmpV[i]=0.0;
      
      double factor = 0.0;
      
      // tmp = M * phi
      for (int i=0; i<n; i++) { // foreach col
	double *mijPtr = &Mptr[i*n];
	for (int j=0; j<n; j++) { // foreach row
	  double mij = *mijPtr++;
	  tmpV[j] += mij*eigenvector[k*n+j];
	}
      }

      // phi^t * tmp
      for (int i=0; i<n; i++) { // foreach col
	factor += eigenvector[k*n+i]*tmpV[i];
      }

      if (factor >= 0) {
	factor=1.0/sqrt(factor);

	for (int i=0; i<n; i++)
	  eigenvector[k*n+i] = eigenvector[k*n+i]*factor;
      }
    }
    */
    delete [] kCopy;
    delete [] mCopy;
    delete [] tmpV;
    
    // sort eigenvalues based on size
    this->sort(numEigen, eigenvalue, sortingID);

    for (int i=0; i<numEigen; i++) {
        if (eigenvalue[i] == DBL_MAX) {
	    opserr << "SymmGeneralizedEigenSolver::solve() - the eigenvalue "
		    << i+1 << " is numerically undetermined or infinite\n";
        } 
    }

    int lworkOpt = (int) work[0];
    if (lwork < lworkOpt) {
        opserr << "SymmGeneralizedEigenSolver::solve() - optimal workspace size "
                << lworkOpt << " is larger than provided workspace size "
                << lwork << " consider increasing workspace\n";
    }

    // clean up the memory
    delete [] work;
    delete [] iwork;
    delete [] ifail;        

    return 0;
}


int SymmGeneralizedEigenSolver::setSize()
{
    int size = theSOE->size;

    if (eigenV == 0 || eigenV->Size() != size) {
        if (eigenV != 0)
            delete eigenV;

        eigenV = new Vector(size);
        if (eigenV == 0 || eigenV->Size() != size) {
            opserr << "SymmGeneralizedEigenSolver::setSize() ";
            opserr << " - ran out of memory for eigenVector of size ";
            opserr << theSOE->size << endln;
            return -2;	    
        }
    }

    return 0;
}


int SymmGeneralizedEigenSolver::setEigenSOE(SymmGeneralizedEigenSOE &thesoe)
{
    theSOE = &thesoe;

    return 0;
}


const Vector& SymmGeneralizedEigenSolver::getEigenvector(int mode)
{
    if (mode <= 0 || mode > numEigen) {
        opserr << "SymmGeneralizedEigenSolver::getEigenVector() - mode "
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
        opserr << "SymmGeneralizedEigenSolver::getEigenvector() - "
            << "eigenvectors not computed yet\n";
        eigenV->Zero();
    }      

    //opserr << "EIGEN VECTOR: " << *eigenV;
    
    return *eigenV;
    
}


double SymmGeneralizedEigenSolver::getEigenvalue(int mode)
{
    if (mode <= 0 || mode > numEigen) {
        opserr << "SymmGeneralizedEigenSolver::getEigenvalue() - mode " 
            << mode << " is out of range (1 - " << numEigen << ")\n";
        return 0.0;
    }

    if (eigenvalue != 0) {
        return eigenvalue[mode-1];
    }
    else {
        opserr << "SymmGeneralizedEigenSolver::getEigenvalue() - "
            << "eigenvalues not yet computed\n";
        return 0.0;
    }      
}


int SymmGeneralizedEigenSolver::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}


int SymmGeneralizedEigenSolver::recvSelf(int commitTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    return 0;
}


void SymmGeneralizedEigenSolver::sort(int length, double *x, int *id)
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
