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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/SuperLU.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/sparseGEN/SuperLU.h
//
// Written: fmk 
// Created: Tue 11/96
// Revision: A
//
// Description: This file contains the implementation of SuperLU
//
// What: "@(#) SuperLU.h, revA"

#include <SuperLU.h>
#include <SparseGenColLinSOE.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


SuperLU::SuperLU(int perm, double thrsh, int relx, int panel)
:SparseGenColLinSolver(SOLVER_TAGS_SuperLU),perm_r(0),perm_c(0),
 etree(0), sizePerm(0),
 relax(relx), permSpec(perm),  panelSize(panel), thresh(thrsh)
{
  refact[0] = 'N';    
}


SuperLU::~SuperLU()
{
  if (perm_r != 0)
	delete [] perm_r;
  if (perm_c != 0)
	delete [] perm_c;
  if (etree != 0)
	delete [] etree;
}

extern "C" void  dgstrf (char *refact, SuperMatrix* A, double pivThresh, 
			 double dropTol, int relax, 
			 int panelSize, int *etree,
			 void *work, int lwork, int *perm_r, int *perm_c, 
			 SuperMatrix *L, SuperMatrix *U, int *info);

extern "C" void  dgstrs (char *trans, SuperMatrix *L, SuperMatrix *U, 
			 int *perm_r, int *perm_c, 
			 SuperMatrix *B, int *info);    


extern "C" void    StatInit(int, int);

extern "C" void    dCreate_CompCol_Matrix(SuperMatrix *, int, int, int, double *,
					  int *, int *, Stype_t, Dtype_t, Mtype_t);


extern "C" void    get_perm_c(int, SuperMatrix *, int *);

extern "C" void    sp_preorder (char*, SuperMatrix*, int*, int*, SuperMatrix*);

extern "C" void    dCreate_Dense_Matrix(SuperMatrix *, int, int, double *, int,
	       	                     Stype_t, Dtype_t, Mtype_t);


int
SuperLU::solve(void)
{
    if (theSOE == 0) {
	cerr << "WARNING SuperLU::solve(void)- ";
	cerr << " No LinearSOE object has been set\n";
	return -1;
    }
    
    int n = theSOE->size;
    
    // check for quick return
    if (n == 0)
	return 0;

    if (sizePerm == 0) {
	cerr << "WARNING SuperLU::solve(void)- ";
	cerr << " size for row and col permutations 0 - has setSize() been called?\n";
	return -1;
    }

    // first copy B into X
    double *Xptr = theSOE->X;
    double *Bptr = theSOE->B;
    for (int i=0; i<n; i++)
	*(Xptr++) = *(Bptr++);

    if (theSOE->factored == false) {
	// factor the matrix
	int info;
	dgstrf(refact, &AC, thresh, 0.0, relax, panelSize, etree, NULL, 0,
	       perm_r, perm_c, &L, &U, &info);

	if (info != 0) {	
	   cerr << "WARNING SuperLU::solve(void)- ";
	   cerr << " Error " << info << " returned in factorization dgstrf()\n";
	   return -info;
	}
	refact[0] = 'Y';
	theSOE->factored = true;
    }	

    // do forward and backward substitution
    char trans[1]; trans[0] = 'N';
    int info;
    dgstrs (trans, &L, &U, perm_r, perm_c, &B, &info);    

    if (info != 0) {	
       cerr << "WARNING SuperLU::solve(void)- ";
       cerr << " Error " << info << " returned in substitution dgstrs()\n";
       return -info;
    }

    return 0;
}




int
SuperLU::setSize()
{
    int n = theSOE->size;
    if (n > 0) {

      // create space for the permutation vectors 
      // and the elimination tree
      if (sizePerm < n) {

	if (perm_r != 0)
	  delete [] perm_r;
	perm_r = new int[n];		

	if (perm_c != 0)
	  delete [] perm_c;
	perm_c = new int[n];		

	if (etree != 0)
	  delete [] etree;
	etree = new int[n];		

	if (perm_r == 0 || perm_c == 0 || etree == 0) {
	  cerr << "WARNING SuperLU::setSize()";
	  cerr << " - ran out of memory\n";
	  sizePerm = 0;
	  return -1;
	}		
	sizePerm = n;
      }

      // initialisation
      StatInit(panelSize, relax);

      // create the SuperMatrix A	
      dCreate_CompCol_Matrix(&A, n, n, theSOE->nnz, theSOE->A, 
			     theSOE->rowA, theSOE->colStartA, 
			     NC, _D, GE);
      // obtain and apply column permutation to give SuperMatrix AC
      get_perm_c(permSpec, &A, perm_c);

      sp_preorder(refact, &A, perm_c, etree, &AC);

      // create the rhs SuperMatrix B 
      dCreate_Dense_Matrix(&B, n, 1, theSOE->X, n, DN, _D, GE);
	
      // set the refact variable to 'N' after first factorization with new size 
      // can set to 'Y'.
      refact[0] = 'N';

    } else if (n == 0)
	return 0;
    else {
	cerr << "WARNING SuperLU::setSize()";
	cerr << " - order of system <  0\n";
	return -1;	
    }
	
    return 0;
}

int
SuperLU::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
SuperLU::recvSelf(int ctag,
		  Channel &theChannel, 
		  FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}















