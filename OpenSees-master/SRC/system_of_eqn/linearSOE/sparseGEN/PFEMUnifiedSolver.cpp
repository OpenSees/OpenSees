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
          
// $Revision: 1.0 $
// $Date: 2014-09-06 13:53:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMUnifiedSolver.h,v $

// Written: Minjie Zhu

#include <PFEMUnifiedSolver.h>
#include <PFEMUnifiedLinSOE.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <Timer.h>
#include <Matrix.h>

#ifdef _WIN32
extern "C" int  DGESV(int *N, int *NRHS, double *A, int *LDA,
			      int *iPiv, double *B, int *LDB, int *INFO);
#else
extern "C" int dgesv_(int *N, int *NRHS, double *A, int *LDA, int *iPiv,
		      double *B, int *LDB, int *INFO);
#endif

PFEMUnifiedSolver::PFEMUnifiedSolver()
    :LinearSOESolver(SOLVER_TAGS_PFEMUnifiedSolver), theSOE(0),
     Msym(0), Mnum(0)
{
}

PFEMUnifiedSolver::~PFEMUnifiedSolver()
{
    if(Msym != 0) {
        cs_sfree(Msym);
    }
    if(Mnum != 0) {
        cs_nfree(Mnum);
    }
}

int
PFEMUnifiedSolver::solve()
{
    cs* M = theSOE->M;
    cs* G = theSOE->G;
    cs* Gt = theSOE->Gt;
    cs* Mp = theSOE->Mp;
    const Vector& B = theSOE->getB();
    const ID& dofType = theSOE->getDofType();
    const ID& dofID = theSOE->getDofID();

    int Vsize = M->n;
    int Psize = Mp->n;
    int size = B.Size();

    if(Vsize <= 0) {
        opserr<<"WARNING: Vsize<=0 -- ";
        opserr<<"PFEMUnifiedSolver::solve\n";
        return -1;
    }

    Timer timer;

    timer.start();

    // numeric LU factorization of K
    if(Msym == 0) {
        opserr<<"WARNING: setSize has not been called";
        opserr<<" -- PFEMSolver::solve\n";
        return -1;
    }
    if(Mnum != 0) {
        cs_nfree(Mnum);
        Mnum = 0;
    }
    Mnum = cs_lu(M, Msym, 1e-6);
    if(Mnum == 0) {
        opserr<<"WARNING: failed to do LU factorization of K";
        opserr<<" -- PFEMSolver::solve\n";
        return -1;
    }
    timer.pause();
    opserr<<"LU M time = "<<timer.getReal()<<"\n";

    timer.start();

    // unified predictor : deltaV1 = M^{-1} * rv
    Vector dV1(Vsize);

    for(int i=0; i<size; i++) {        // row
        int rowtype = dofType(i);      // row type
        int rowid = dofID(i);          // row id
        if(rowtype < 3 && rowtype>=0) {
            dV1(rowid) = B(i);         // rv
        }
    }

    // M^{-1}*rv
    Vector x(Vsize);
    double* dV1_ptr = &dV1(0);
    double* x_ptr = &x(0);
    cs_ipvec(Mnum->pinv, dV1_ptr, x_ptr, Vsize);
    cs_lsolve(Mnum->L, x_ptr);
    cs_usolve(Mnum->U, x_ptr);
    cs_ipvec(Msym->q, x_ptr, dV1_ptr, Vsize);
    timer.pause();
    opserr<<"dV1 = "<<dV1.Norm()<<"\n";
    opserr<<"predictor time = "<<timer.getReal()<<"\n";
    timer.start();
    
    // no pressure
    if(Psize <= 0) {
        Vector X(size);
        for(int i=0; i<size; i++) {            // row
            int rowtype = dofType(i);          // row type
            int rowid = dofID(i); 
            if(rowtype>=0 && rowtype<3) {
                X(i) = dV1(rowid);
            }
        }
        theSOE->setX(X);
        return 0;
    }

    // L = Mp-Gt*M^{-1}*G
    Matrix L(Psize,Psize);

    // each column of L
    int index = 0;
    css* msym = Msym;
    csn* mnum = Mnum;
    omp_set_num_threads(omp_get_num_procs());

    timer.start();
#pragma omp parallel for default(none),shared(Mp,mnum,msym,Gt,G,L,Psize,Vsize),private(index),schedule(static)
    for(index=0; index<Psize; index++) {
    	int j = index;

    	// -M^{-1}*Gj, Mpj
    	Vector Gj(Vsize), Mpj(Psize);
    	for(int k=G->p[j]; k<G->p[j+1]; k++) {
	    Gj(G->i[k]) = -(G->x[k]);
    	}
    	for(int k=Mp->p[j]; k<Mp->p[j+1]; k++) {
	    Mpj(Mp->i[k]) = Mp->x[k];
    	}
    	Vector x(Vsize);
    	double* x_ptr = &x(0);
    	double* Gj_ptr = &Gj(0);
    	cs_ipvec(mnum->pinv, Gj_ptr, x_ptr, Vsize);
    	cs_lsolve(mnum->L, x_ptr);
    	cs_usolve(mnum->U, x_ptr);
    	cs_ipvec(msym->q, x_ptr, Gj_ptr, Vsize);

    	// Mpj-Gt*M^{-1}*Gj
	double* Mpj_ptr = &Mpj(0);
    	cs_gaxpy(Gt,Gj_ptr,Mpj_ptr);

    	// add to L
    	for(int i=0; i<Psize; i++) {
	    L(j,i) = Mpj(i);
    	}
    }

    timer.pause();
    opserr<<"compute L time = "<<timer.getReal()<<"\n";

    // pressure : 
    timer.start();
    Vector dP(Psize);
    double* dP_ptr = &dP(0);

    // Gt*dV1
    dV1_ptr = &dV1(0);
    cs_gaxpy(Gt,dV1_ptr,dP_ptr);

    // rp-Gt*dV1
    for(int i=0; i<size; i++) {        // row

        int rowtype = dofType(i);      // row type
        int rowid = dofID(i);          // row id
        if(rowtype == 3) {
            dP(rowid) = B(i)-dP(rowid);          
        }
    }
    
    // solve
    double* L_ptr = &L(0,0);
    int nrhs = 1;
    int ldL = Psize;
    int ldP = Psize;
    ID ipiv(Psize);
    int info;

#ifdef _WIN32
    DGESV(&n,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
#else
    dgesv_(&Psize,&nrhs,L_ptr,&ldL,&ipiv(0),dP_ptr,&ldP,&info);
#endif

    timer.pause();
    opserr<<"dP = "<<dP.Norm()<<"\n";
    opserr<<"pressure time = "<<timer.getReal()<<"\n";

    // unified corrector
    timer.start();
    Vector dV(Vsize);
    double* dV_ptr = &dV(0);

    // G*dP
    dP_ptr = &dP(0);
    cs_gaxpy(G,dP_ptr,dV_ptr);

    // M^{-1}*G*dP
    x.resize(Vsize); x.Zero();
    x_ptr = &x(0);
    cs_ipvec(Mnum->pinv, dV_ptr, x_ptr, Vsize);
    cs_lsolve(Mnum->L, x_ptr);
    cs_usolve(Mnum->U, x_ptr);
    cs_ipvec(Msym->q, x_ptr, dV_ptr, Vsize);

    // dV1-M^{-1}*G*dP
    dV = dV1-dV;
    timer.pause();
    opserr<<"corrector time = "<<timer.getReal()<<"\n";
    
    // copy to X
    Vector X(size);
    for(int i=0; i<size; i++) {            // row
        int rowtype = dofType(i);          // row type
        int rowid = dofID(i); 
        if(rowtype>=0 && rowtype<3) {
            X(i) = dV(rowid);
        } else if(rowtype == 3) {
            X(i) = dP(rowid);
        }

    }
    theSOE->setX(X);
    
    return 0;
}

int PFEMUnifiedSolver::setSize()
{
    Timer timer;
    timer.start();
    if(Msym != 0) {
        cs_sfree(Msym);
        Msym = 0;
    }

    if(Mnum != 0) {
        cs_nfree(Mnum);
        Mnum = 0;
    }

    cs* M = theSOE->M;
    if(M->n > 0) {
        Msym = cs_sqr(3, M, 0);
        if(Msym == 0) {
            opserr<<"WARNING: failed to do symbolic analysis of M";
            opserr<<" -- PFEMSolver::setSize\n";
            return -1;
        }
    }

    timer.pause();
    opserr<<"analysis time = "<<timer.getReal()<<"\n";

    return 0;
}

int
PFEMUnifiedSolver::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMUnifiedSolver::recvSelf(int ctag,
		  Channel &theChannel, 
		  FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int 
PFEMUnifiedSolver::setLinearSOE(PFEMUnifiedLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}
