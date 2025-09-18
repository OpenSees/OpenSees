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
// $Date: 2003-02-14 23:02:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/itpack/ItpackLinSolver.cpp,v $
                                                                        
// Written: MHS
// Created: Sept 2001
//
// Description: This file contains the class definition for ItpackLinSolver.
// ItpackLinSolver is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes.  Instances of ItpackLinSolver 
// are used to solve a system of equations of type ItpackLinSOE.

#include <ItpackLinSolver.h>
#include <ItpackLinSOE.h>
#include <ID.h>
#include <elementAPI.h>

void* OPS_ItpackLinSolver()
{
  int method = 1;
  int numData = 1;
  
  int nArgs = OPS_GetNumRemainingInputArgs();
  if (nArgs == 0) {
    opserr << "WARNING Itpack -- no method specified, using JCG" << endln;
  }
  if (nArgs > 0 && OPS_GetIntInput(&numData,&method) < 0) {
    opserr << "WARNING Itpack -- error reading method\n";
    return 0;
  }
  
  int iter = 100;
  double omega = 1.0;
  bool symmetric = true;
  double zeta = 5e-6;
  while (OPS_GetNumRemainingInputArgs() > 1) {
    const char *arg = OPS_GetString();
    if (strcmp(arg,"-iter") == 0) {
      if (OPS_GetIntInput(&numData,&iter) < 0)
	return 0;
    }
    if (strcmp(arg,"-omega") == 0) {
      if (OPS_GetDoubleInput(&numData,&omega) < 0)
	return 0;
    }
    if (strcmp(arg,"-symmetric") == 0) {
      int symm;
      if (OPS_GetIntInput(&numData,&symm) < 0)
	return 0;
      symmetric = (symm != 0) ? true : false;
    }        
    if (strcmp(arg,"-zeta") == 0) {
      if (OPS_GetDoubleInput(&numData,&zeta) < 0)
	return 0;
    }
  }
  
  ItpackLinSolver *theSolver = new ItpackLinSolver(method, iter, omega, zeta);
  return new ItpackLinSOE(*theSolver, symmetric);  
}

ItpackLinSolver::ItpackLinSolver(int meth, int iter, double om, double z)
  :LinearSOESolver(SOLVER_TAGS_Itpack),
   theSOE(0), IA(0), JA(0), n(0),
   iwksp(0), wksp(0), nwksp(0), maxIter(iter),
   method(meth), omega(om), zeta(z)
{

}    

ItpackLinSolver::ItpackLinSolver()
  :LinearSOESolver(SOLVER_TAGS_Itpack),
   theSOE(0), IA(0), JA(0), n(0),
   iwksp(0), wksp(0), nwksp(0), maxIter(0),
   method(0), omega(0.0), zeta(0.0)
{

}    

ItpackLinSolver::~ItpackLinSolver()    
{
  if (IA != 0)
    delete [] IA;

  if (JA != 0)
    delete [] JA;

  if (iwksp != 0)
    delete [] iwksp;

  if (wksp != 0)
    delete [] wksp;
}    

int 
ItpackLinSolver::setLinearSOE(ItpackLinSOE &theItpackSOE)
{
  theSOE = &theItpackSOE;

  return 0;
}

int
ItpackLinSolver::setSize(void)
{
  // Get number of equations from SOE
  n = theSOE->size;

  if (n > 0) {
    if (iwksp != 0)
      delete [] iwksp;
    iwksp = new int[3*n];
  }

  // Should be 2*maxIter for symmetric storage, 4*maxIter for nonsymmetric
  int ncg = 4*maxIter;
  
  // Order of the black subsystem
  int nb = theSOE->size-1; // I think this is what it should be

  switch(method) {
  case ItpackJCG:
    nwksp = 4*n + ncg;
    break;
  case ItpackJSI: case ItpackJ:
    nwksp = 2*n;
    break;
  case ItpackSOR: case ItpackGS: case ItpackSORFixed:
    nwksp = n;
    break;
  case ItpackSSORCG:
    nwksp = 6*n + ncg;
    break;
  case ItpackSSORSI: case ItpackSSORFixed:
    nwksp = 5*n;
    break;
  case ItpackRSCG:
    nwksp = n + 3*nb + ncg;
    break;
  case ItpackRSSI: case ItpackRS:
    nwksp = n + nb;
    break;
  default:
    nwksp = 6*n + ncg;
    break;
  }

  if (nwksp > 0) {
    if (wksp != 0)
      delete [] wksp;
    wksp = new double[nwksp];
  }

  // Get number of nonzeros from the SOE
  int nnz = theSOE->nnz;
  
  if (nnz > 0) {
    if (JA != 0)
      delete [] JA;
    JA = new int [nnz];
  }

  int *jaPtr = theSOE->colA;
  int i;
  // Add one for FORTRAN indexing
  for (i = 0; i < nnz; i++)
    JA[i] = jaPtr[i] + 1;

  if (n > 0) {
    if (IA != 0)
      delete [] IA;
    IA = new int [n+1];
  }
  
  int *iaPtr = theSOE->rowStartA;
  // Add one for FORTRAN indexing
  for (i = 0; i <= n; i++) 
    IA[i] = iaPtr[i] + 1;

  //opserr << "ItpackLinSolver::setSize() -- method = " << method
  //	 << ", ITMAX = " << maxIter << endln;

  return 0;
}

int
ItpackLinSolver::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
ItpackLinSolver::recvSelf(int commitTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
  return -1;
}

#ifdef _WIN32
#define DFAULT dfault_
#define VFILL  vfill_
#define JCG    jcg_
#define JSI    jsi_
#define SOR    sor_
#define SSORCG ssorcg_
#define SSORSI ssorsi_
#define RSCG   rscg_
#define RSSI   rssi_

extern "C" int _stdcall DFAULT(int *IPARM, double *RPARM);

extern "C" int _stdcall VFILL(int *N, double *U, double *VAL);

extern "C" int _stdcall JCG(int *N, int *IA, int *JA, double *A, double *RHS,
			    double *U, int *IWKSP, int *NW, double *WKSP,
			    int *IPARM, double *RPARM, int *IER);

extern "C" int _stdcall JSI(int *N, int *IA, int *JA, double *A, double *RHS,
			    double *U, int *IWKSP, int *NW, double *WKSP,
			    int *IPARM, double *RPARM, int *IER);

extern "C" int _stdcall SOR(int *N, int *IA, int *JA, double *A, double *RHS,
			    double *U, int *IWKSP, int *NW, double *WKSP,
			    int *IPARM, double *RPARM, int *IER);

extern "C" int _stdcall SSORCG(int *N, int *IA, int *JA, double *A, double *RHS,
			       double *U, int *IWKSP, int *NW, double *WKSP,
			       int *IPARM, double *RPARM, int *IER);

extern "C" int _stdcall SSORSI(int *N, int *IA, int *JA, double *A, double *RHS,
			       double *U, int *IWKSP, int *NW, double *WKSP,
			       int *IPARM, double *RPARM, int *IER);

extern "C" int _stdcall RSCG(int *N, int *IA, int *JA, double *A, double *RHS,
			     double *U, int *IWKSP, int *NW, double *WKSP,
			     int *IPARM, double *RPARM, int *IER);

extern "C" int _stdcall RSSI(int *N, int *IA, int *JA, double *A, double *RHS,
			     double *U, int *IWKSP, int *NW, double *WKSP,
			     int *IPARM, double *RPARM, int *IER);

/* NOT HERE B..A..
#define DFAULT dfault_
#define VFILL  vfill_
#define JCG    jcg_
#define JSI    jsi_
#define SOR    sor_
#define SSORCG ssorcg_
#define SSORSI ssorsi_
#define RSCG   rscg_
#define RSSI   rssi_
*/

#else

extern "C" int dfault_(int *iparm, double *rparm);

extern "C" int vfill_(int *n, double *u, double *val);

extern "C" int jcg_(int *n, int *ia, int *ja, double *a, double *rhs,
		    double *u, int *iwksp, int *nw, double *wksp,
		    int *iparm, double *rparm, int *ier);

extern "C" int jsi_(int *n, int *ia, int *ja, double *a, double *rhs,
		    double *u, int *iwksp, int *nw, double *wksp,
		    int *iparm, double *rparm, int *ier);

extern "C" int sor_(int *n, int *ia, int *ja, double *a, double *rhs,
		    double *u, int *iwksp, int *nw, double *wksp,
		    int *iparm, double *rparm, int *ier);

extern "C" int ssorcg_(int *n, int *ia, int *ja, double *a, double *rhs,
		       double *u, int *iwksp, int *nw, double *wksp,
		       int *iparm, double *rparm, int *ier);

extern "C" int ssorsi_(int *n, int *ia, int *ja, double *a, double *rhs,
		       double *u, int *iwksp, int *nw, double *wksp,
		       int *iparm, double *rparm, int *ier);

extern "C" int rscg_(int *n, int *ia, int *ja, double *a, double *rhs,
		     double *u, int *iwksp, int *nw, double *wksp,
		     int *iparm, double *rparm, int *ier);

extern "C" int rssi_(int *n, int *ia, int *ja, double *a, double *rhs,
		     double *u, int *iwksp, int *nw, double *wksp,
		     int *iparm, double *rparm, int *ier);

#endif

int
ItpackLinSolver::solve(void)
{
  // Let ITPACK fill in default parameter values
  dfault_(iparm, rparm);

  // Override defaults for "textbook" methods
  switch (method) {
  case ItpackJ:
    iparm[5] = 0; iparm[6] = 2; break;
  case ItpackGS:
    iparm[5] = 0; break;
  case ItpackSORFixed:
    iparm[5] = 0; rparm[4] = omega; break;
  case ItpackSSORFixed:
    iparm[5] = 0; rparm[4] = omega; break;
  case ItpackRS:
    iparm[5] = 0; break;
  default:
    break;
  }

  // Overwrite default max number of iterations
  iparm[0] = maxIter;

  /* Set convergence tolerance for inexact Newton methods.
   * ITPACK stopping criterion: ||Ax-b|| / ||x|| * C < rparm[0]
   * For (1+p)-order convergence: rparm[0] < ||b||^(1+p)
   * Reference: Dembo et al. (1982) Inexact Newton Methods. 
   *            SIAM Journal on Numerical Analysis, 19(2), 400–408. 
   *            https://doi.org/10.1137/0719025
   */
  double norm_b = theSOE->normRHS();
  double norm_b_squared = norm_b * norm_b;
  rparm[0] = (zeta < norm_b_squared) ? zeta : norm_b_squared;

  // Print flag
  iparm[1] = 0;
  
  // Sparse matrix storage scheme (0 = symmetric, 1 = nonsymmetric)
  iparm[4] = theSOE->symmetric ? 0 : 1;

  double *aPtr = theSOE->A;
  double *xPtr = theSOE->X;
  double *bPtr = theSOE->B;

  int *iaPtr = IA;
  int *jaPtr = JA;

  // Number of non-zero entries in matrix
  int nnz = iaPtr[n]-1;

  // Copy original ordering of column indices from the SOE
  // because ITPACK will reorder the sparse matrix representation
  if (theSOE->Aformed == false) {
    int *soeColA = theSOE->colA;
    for (int i = 0; i < nnz; i++) {
      jaPtr[i] = soeColA[i] + 1;  // Add one for FORTRAN indexing
    }
  }

  int ier = 0;

  // Order of the black subsystem
  int nb = theSOE->size-1; // I think this is what it should be

  // Fill the x vector with zeros as initial guess to solution of Ax=b
  double val = 0.0;
  vfill_(&n, xPtr, &val);

  switch (method) {
  case ItpackJCG:
    jcg_(&n, iaPtr, jaPtr, aPtr, bPtr, xPtr,
	 iwksp, &nwksp, wksp, iparm, rparm, &ier);
    break;
  case ItpackJSI: case ItpackJ:
    jsi_(&n, iaPtr, jaPtr, aPtr, bPtr, xPtr,
	 iwksp, &nwksp, wksp, iparm, rparm, &ier);
    break;
  case ItpackSOR: case ItpackGS: case ItpackSORFixed:
    sor_(&n, iaPtr, jaPtr, aPtr, bPtr, xPtr,
	 iwksp, &nwksp, wksp, iparm, rparm, &ier);
    break;
  case ItpackSSORCG:
    ssorcg_(&n, iaPtr, jaPtr, aPtr, bPtr, xPtr,
	    iwksp, &nwksp, wksp, iparm, rparm, &ier);
    break;
  case ItpackSSORSI: case ItpackSSORFixed:
    ssorsi_(&n, iaPtr, jaPtr, aPtr, bPtr, xPtr,
	    iwksp, &nwksp, wksp, iparm, rparm, &ier);
    break;
  case ItpackRSCG:
    iparm[8] = nb;
    rscg_(&n, iaPtr, jaPtr, aPtr, bPtr, xPtr,
	  iwksp, &nwksp, wksp, iparm, rparm, &ier);
    break;
  case ItpackRSSI: case ItpackRS:
    iparm[8] = nb;
    rssi_(&n, iaPtr, jaPtr, aPtr, bPtr, xPtr,
	  iwksp, &nwksp, wksp, iparm, rparm, &ier);
    break;
  default:
    opserr << method << " -- unknown method type in ItpackLinSolver::solve()" << endln;
    break;    
  }

  // Tell the SOE that solve() has been called
  theSOE->Aformed = true;

  if (ier > 0) {
    opserr << "ItpackLinSolver::solve() -- returned ier = " << ier << endln;
    opserr << "Maximum iterations: " << iparm[0] << ", rparm = " << rparm[0] << endln;
    return -ier;
  }
  else {
    //opserr << "Converged in " << iparm[0] << " iterations" << endln;
    return 0;
  }
}
