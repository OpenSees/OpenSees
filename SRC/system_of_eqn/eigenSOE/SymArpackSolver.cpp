// File: ~/system_of_eqn/linearSOE/LawSolver/SymArpackSolver.C
//
// Written: Jun Peng
// Created: 2/1999
// Revision: A
//
// Description: This file contains the class definition for 
// SymArpackSolver. It solves the SymArpackSOE object by calling
// some "C" functions. The solver used here is generalized sparse
// solver. The user can choose three different ordering schema.
//
// What: "@(#) SymArpackSolver.C, revA"


#include "SymArpackSOE.h"
#include "SymArpackSolver.h"
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <AnalysisModel.h>
#include <DOF_GrpIter.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <Integrator.h>

extern "C" {
  #include <FeStructs.h>
  #include <globalVars.h>
  #include <tim.h>
}

#ifdef _WIN32
extern "C" int  DSAUPD(int *ido, char* bmat, 
			       int *n, char *which,
			       int *nev, 
			       double *tol, double *resid, int *ncv, double *v, 
			       int *ldv,
			       int *iparam, int *ipntr, double *workd, double *workl,
			       int *lworkl, int *info);

extern "C" int DSEUPD(bool *rvec, char *howmny,
			       logical *select, double *d, double *z,
			       int *ldz, double *sigma, char *bmat,
			       int 	*n, char *which,
			       int *nev, double *tol, double *resid, int *ncv, 
			       double *v,
			       int *ldv, int *iparam, int *ipntr, double *workd, 
			       double *workl, int *lworkl, int *info);

#else

extern "C" int dsaupd_(int *ido, char* bmat, int *n, char *which, int *nev, 
		       double *tol, double *resid, int *ncv, double *v, int *ldv,
		       int *iparam, int *ipntr, double *workd, double *workl,
		       int *lworkl, int *info);

extern "C" int dseupd_(bool *rvec, char *howmny, logical *select, double *d, double *z,
		       int *ldz, double *sigma, char *bmat, int *n, char *which,
		       int *nev, double *tol, double *resid, int *ncv, double *v,
		       int *ldv, int *iparam, int *ipntr, double *workd, 
		       double *workl, int *lworkl, int *info);

#endif

SymArpackSolver::SymArpackSolver(int numE)
:EigenSolver(EigenSOLVER_TAGS_SymArpackSolver),
 theSOE(0), factored(false), theNev(numE), value(0), eigenV(0)
{
    // nothing to do.
}


SymArpackSolver::~SymArpackSolver()
{ 
    if (value != 0)
        delete [] value;

    if (vector !=0)
        delete [] vector;

    if (eigenV !=0)
        delete eigenV;
}


extern "C" int pfsfct(int neqns, double *diag, double **penv, int nblks, int *xblk,
		      OFFDBLK **begblk, OFFDBLK *first, int *rowblks);

extern "C" void pfsslv(int neqns, double *diag, double **penv, int nblks, 
		       int *xblk, double *rhs, OFFDBLK **begblk);



int
SymArpackSolver::solve(int numModes, bool generalized, bool findSmallest)
{ 
  if (generalized == false) {
    opserr << "BandArpackSolver::solve(int numMode, bool generalized) - only solves generalized problem\n";
    return -1;
  }

  if (theSOE == 0) {
    opserr << "WARNING SymArpackSolver::solve(void)- ";
    opserr << " No EigenSOE object has been set\n";
    return -1;
  }

    int      nblks = theSOE->nblks;
    int      *xblk = theSOE->xblk;
    //    int      *invp = theSOE->invp;
    double   *diag = theSOE->diag;
    double   **penv = theSOE->penv;
    int      *rowblks = theSOE->rowblks;
    OFFDBLK  **begblk = theSOE->begblk;
    OFFDBLK  *first = theSOE->first;
    
    int n = theSOE->size;

    // check for quick return
    if (n == 0)
	return 0;

//	timer (FACTOR);
	if (factored == false) {

	   //factor the matrix
	   //call the "C" function to do the numerical factorization.
	   int factor;
	   factor = pfsfct(n, diag, penv, nblks, xblk, begblk, first, rowblks);
	   if (factor>0) {
		  opserr << "In SymArpackSolver: error in factorization.\n";
		  return -1;
	   }
	   factored = true;
	}


	int nev = numModes;

    int ncv = getNCV(n, nev);

    // set up the space for ARPACK functions.
    int ldv = n;
    
    double *v = new double[ldv * ncv];
    double *workl = new double[ncv * (ncv+8)];
    double *workd = new double[3 * n];
    double *d = new double[ncv * 2];
    double *resid = new double[n];
	
//    char *which = new char[2];
    char bmat = 'G';

    static char which[3];
    if (findSmallest == true) {
      strcpy(which, "LM");
    }  else {
    strcpy(which, "SM");
    }    

    int *iparam = new int[11];
    int *iwork = new int[n];

    int lworkl = ncv*ncv + 8*ncv;
    int maxitr, mode, nconv;

    double tol = 0.0;
    int info = 0;
    maxitr = 1000;
    mode = 3;

    iparam[2] = maxitr;
    iparam[6] = mode; 

    bool rvec = true;
    char howmy = 'A';
    logical *select = new logical[ncv];

    iparam[0] = 1;
    int ido = 0;
    int *ipntr = new int[11];
    int ierr = 0;
    unsigned int sizeWhich =2;
    unsigned int sizeBmat =1;
    unsigned int sizeHowmany =1;



    while (1) { 
 
#ifdef _WIN32
	DSAUPD(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
	       iparam, ipntr, workd, workl, &lworkl, &info);
#else
	dsaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
		iparam, ipntr, workd, workl, &lworkl, &info);
#endif
      if (ido == -1) {

	  myMv(n, &workd[ipntr[0]-1], &workd[ipntr[1]-1]); 

	  pfsslv(n, diag, penv, nblks, xblk, &workd[ipntr[1] - 1], begblk);
	  continue;
      } else if (ido == 1) {
          double ratio = 1.0;
	  myCopy(n, &workd[ipntr[2]-1], &workd[ipntr[1]-1]);

	  pfsslv(n, diag, penv, nblks, xblk, &workd[ipntr[1] - 1], begblk);

	  continue;
      } else if (ido == 2) {     
	  myMv(n, &workd[ipntr[0]-1], &workd[ipntr[1]-1]);

	  continue;
      }

      break;
    }

    if (info < 0) {
        opserr << "BandArpackSolver::Error with _saupd info = " << info <<endln;
	return info;
    } else {
        if (info == 1) {
	    opserr << "BandArpackSolver::Maximum number of iteration reached." << endln;
	} else if (info == 3) {
	    opserr << "BandArpackSolver::No Shifts could be applied during implicit," << endln;
	    opserr << "Arnoldi update, try increasing NCV." << endln;
	}
	double sigma = theSOE->shift;
	if (iparam[4] > 0) {
#ifdef _WIN32
	    DSEUPD(&rvec, &howmy, select, d, v, &ldv, &sigma, &bmat, &n, which,
		   &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd,
		   workl, &lworkl, &info);
#else


	    dseupd_(&rvec, &howmy, select, d, v, &ldv, &sigma, &bmat, &n, which,
		    &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd,
		    workl, &lworkl, &info);
#endif
	    if (info != 0) {
	        opserr << "BandArpackSolver::Error with dseupd_" << info;
		return -1;
	    }
	}
    }

    value = d;
    vector = v;

    theSOE->factored = true;

    delete [] workl;
    delete [] workd;
    delete [] resid;
    delete [] iparam;
    delete [] iwork;
#ifdef _WIN32
    // cannot invoke the destructor on select .. causes seg in windows!!!
#else
    delete [] select;
#endif    

    delete [] ipntr;

    
    return 0;
}


int
SymArpackSolver::setSize()
{
    // nothing to do
    // if iPiv not big enough, free it and get one large enough
    int size = theSOE->size;    
    
    if (eigenV == 0 || eigenV->Size() != size) {
	if (eigenV != 0)
	    delete eigenV;

	eigenV = new Vector(size);
	if (eigenV == 0 || eigenV->Size() != size) {
	    opserr << "WARNING BandGenLinLapackSolver::setSize() ";
	    opserr << " - ran out of memory for eigenVector of size ";
	    opserr << theSOE->size << endln;
	    return -2;	    
	}
    }
	    
    return 0;
}

int
SymArpackSolver::setEigenSOE(SymArpackSOE &theEigenSOE)
{
  theSOE = &theEigenSOE;
  return 0;
}


int 
SymArpackSolver::getNCV(int n, int nev)
{
    int result;
    if (2*nev > nev+8) {
        result = nev+8;
    } else {
        result = 2*nev;
    }

    if (result >= n) {
        result = n;
    }

    return result;
}


void
SymArpackSolver::myMv(int n, double *v, double *result)
{
    double *tempX = new double[n];
    int      *invp = theSOE->invp;    
    int i;
    
    for (i=0; i<n; i++) {
        tempX[i] = v[invp[i]];
    }
    for (i=0; i<n; i++) {
        v[i] = tempX[i];
    }

    Vector x(v, n);
    Vector y(result,n);

    y.Zero();
    AnalysisModel *theAnalysisModel = theSOE->theModel;

    // loop over the FE_Elements
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    while((elePtr = theEles()) != 0) {
      const Vector &a = elePtr->getM_Force(x,1.0);
      y.Assemble(a,elePtr->getID(),1.0);
    }

    // loop over the DOF_Groups
    Integrator *theIntegrator = 0;
    DOF_Group *dofPtr;
    DOF_GrpIter &theDofs = theAnalysisModel->getDOFs();
    while ((dofPtr = theDofs()) != 0) {
      const Vector &a = dofPtr->getM_Force(x,1.0);      
      y.Assemble(a,dofPtr->getID(),1.0);
    }

    for (i=0; i<n; i++) {
        tempX[invp[i]] = result[i];
    }
    for (i=0; i<n; i++) {
        result[i] = tempX[i];
    }
    delete [] tempX;
}

void
SymArpackSolver::myCopy(int n, double *v, double *result)
{
    for (int i=0; i<n; i++) {
        result[i] = v[i];
    }
}

const Vector &
SymArpackSolver::getEigenvector(int mode)
{
    int      *invp = theSOE->invp;
    
    if (mode <= 0 || mode > theNev) {
        opserr << "BandArpackSOE::mode is out of range(1 - nev)";
	exit (0);
    }
    int size = theSOE->size;

    int index = (mode - 1) * size;
    for (int i=0; i<size; i++) {
	(*eigenV)[i] = vector[index + invp[i]];
    }	

    return *eigenV;      
}


double
SymArpackSolver::getEigenvalue(int mode)
{
    if (mode <= 0 || mode > theNev) {
        opserr << "BandArpackSOE::mode is out of range(1 - nev)";
	exit (0);
    }
    return value[mode - 1];
}


int
SymArpackSolver::sendSelf(int cTAg, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
SymArpackSolver::recvSelf(int cTag,
			     Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}



