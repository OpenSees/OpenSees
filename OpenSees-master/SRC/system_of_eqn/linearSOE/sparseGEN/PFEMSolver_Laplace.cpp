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

// $Revision$
// $Date$
// $URL$

//
// Written: Minjie
// Created: Sep 17 2012
//

#include "PFEMSolver_Laplace.h"
#include "PFEMLinSOE.h"
#include <iostream>
#include <cmath>
#include <Timer.h>
#include <elementAPI.h>

void* OPS_PFEMSolver_Laplace()
{
    bool once = false;
    if (OPS_GetNumRemainingInputArgs() > 0) {
	const char* opt = OPS_GetString();
	if (strcmp(opt, "-once") == 0) {
	    once = true;
	}
    }
    PFEMSolver_Laplace* theSolver = new PFEMSolver_Laplace(once);
    return new PFEMLinSOE(*theSolver);
}

PFEMSolver_Laplace::PFEMSolver_Laplace(bool once)
    :PFEMSolver(), MSym(0), MNum(0), LSym(0), LNum(0), theSOE(0),
     numonce(once), factored(false)
{
}

PFEMSolver_Laplace::~PFEMSolver_Laplace()
{
    if (MSym != 0) {
	umfpack_di_free_symbolic(&MSym);
    }
    if (LSym != 0) {
	umfpack_di_free_symbolic(&LSym);
    }
    if (MNum != 0) {
	umfpack_di_free_numeric(&MNum);
    }
    if (LNum != 0) {
	umfpack_di_free_numeric(&LNum);
    }
}

int PFEMSolver_Laplace::setSize()
{
    // reorder rows
    cs* M = theSOE->M;
    cs* Gft = theSOE->Gft;
    cs* Git = theSOE->Git;
    cs* L = theSOE->L;
    cs* mats[4] = {M,Gft,Git,L};
    for (int i=0; i<4; i++) {
	cs* mat = mats[i];
	for (int j=0; j<mat->n; j++) {
	    ID col(0, mat->p[j+1]-mat->p[j]);
	    for (int k=mat->p[j]; k<mat->p[j+1]; k++) {
		col.insert(mat->i[k]);
	    }
	    int index = 0;
	    for (int k=mat->p[j]; k<mat->p[j+1]; k++) {
		mat->i[k] = col[index++];
	    }
	}
    }

    // set default control parameters
    umfpack_di_defaults(Control);
    Control[UMFPACK_PIVOT_TOLERANCE] = 1.0;

    Timer timer;
    timer.start();

    // symbolic analysis of M
    int n = M->n;
    int nnz = M->nzmax;
    if (n>0 && nnz>0) {

	int* Mp = M->p;
	int* Mi = M->i;
	double* Mx = M->x;

	if (MSym != 0) {
	    umfpack_di_free_symbolic(&MSym);
	    MSym = 0;
	}
	int status = umfpack_di_symbolic(n,n,Mp,Mi,Mx,&MSym,Control,Info);

	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: symbolic analysis of M returns ";
	    opserr<<status<<" -- PFEMSolver_Laplace::setsize\n";
	    return -1;
	}
    }

    // symbolic analysis of L
    n = L->n;
    nnz = L->nzmax;
    if (n>0 && nnz>0) {

	int* Lp = L->p;
	int* Li = L->i;
	double* Lx = L->x;

	if (LSym != 0) {
	    umfpack_di_free_symbolic(&LSym);
	    LSym = 0;
	}
	int status = umfpack_di_symbolic(n,n,Lp,Li,Lx,&LSym,Control,Info);

	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: symbolic analysis of L returns ";
	    opserr<<status<<" -- PFEMSolver_Laplace::setsize\n";
	    return -1;
	}
    }

    timer.pause();
    opserr<<"analysis time = "<<timer.getReal()<<"\n";
    return 0;
}

int
PFEMSolver_Laplace::solve()
{
    Timer timer;
    timer.start();
    cs* M = theSOE->M;
    cs* Gft = theSOE->Gft;
    cs* Git = theSOE->Git;
    cs* L = theSOE->L;
    Vector& Mf = theSOE->Mf;
    Vector& X = theSOE->X;
    Vector& B = theSOE->B;
    ID& dofType = theSOE->dofType;
    ID& dofID = theSOE->dofID;

    int Msize = M->n;
    int Isize = Git->n;
    int Ssize = Msize-Isize;
    int Fsize = Mf.Size();
    int Psize = L->n;
    int size = X.Size();

    // check MNum
    if (!numonce && MNum != 0) {
	umfpack_di_free_numeric(&MNum);
	MNum = 0;
    }

    // check LNum
    if (!numonce && LNum != 0) {
	umfpack_di_free_numeric(&LNum);
	LNum = 0;
    }

    // Numerical factorization
    if (MNum==0 && Msize>0) {

	// numeric LU factorization of M
	if(MSym == 0) {
	    opserr<<"WARNING: setSize has not been called";
	    opserr<<" -- PFEMSolver_Laplace::solve\n";
	    return -1;
	}
	
	// numerical analysis
	int* Mp = M->p;
	int* Mi = M->i;
	double* Mx = M->x;
	int status = umfpack_di_numeric(Mp,Mi,Mx,MSym,&MNum,Control,Info);
	
	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: numeric analysis of M returns ";
	    opserr<<status<<" -- PFEMSolver_Laplace::solve\n";
	    return -1;
	}
    }

    if (LNum==0 && Psize>0) {
	
	// numeric LU factorization of L
	if(LSym == 0) {
	    opserr<<"WARNING: setSize has not been called";
	    opserr<<" -- PFEMSolver_Laplace::solve\n";
	    return -1;
	}
	
	// numerical analysis
	int* Lp = L->p;
	int* Li = L->i;
	double* Lx = L->x;
	int status = umfpack_di_numeric(Lp,Li,Lx,LSym,&LNum,Control,Info);
	
	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: numeric analysis of L returns ";
	    opserr<<status<<" -- PFEMSolver_Laplace::solve\n";
	    return -1;
	}
    }

    timer.pause();
    opserr<<"factorization time = "<<timer.getReal()<<"\n";
    timer.start();
    // structure and interface velocity predictor : deltaV1 = M^{-1} * rsi
    Vector deltaV1(Msize);
    if(Msize > 0) {

        // rsi
        for(int i=0; i<size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if(rowtype == 2) {
                deltaV1(rowid+Ssize) = B(i);   // rsi
            } else if(rowtype == 0) {
                deltaV1(rowid) = B(i);         // rsi
            }
        }
	// opserr<<"ri = "<<ri.Norm()<<"\n";
	// opserr<<"rs = "<<rs.Norm()<<"\n";

        // M^{-1}*rsi
        Vector x(Msize);
        double* deltaV1_ptr = &deltaV1(0);
        double* x_ptr = &x(0);
	int* Ap = M->p;
	int* Ai = M->i;
	double* Ax = M->x;
	int status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,x_ptr,deltaV1_ptr,MNum,Control,Info);
	deltaV1 = x;

	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: solving M returns "<<status<<" -- PFEMSolver_Laplace::solve\n";
	    return -1;
	}
    }

    // fluid velocity predictor: deltaVf1 = Mf^{-1} * rf
    Vector deltaVf1(Fsize);
    if(Fsize > 0) {
        // rf
        for(int i=0; i<size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if(rowtype == 1) {
                if(Mf(rowid) == 0) {
                    opserr<<"WANING: Zero Mf at location "<<rowid<<" ";
                    opserr<<" - PFEMLinSOE::solve()\n";
                    return -1;
                }
                deltaVf1(rowid) = B(i)/Mf(rowid);         // rf
            }
        }

    }

    timer.pause();
    // opserr<<"dV1 = "<<deltaV1.Norm()<<"\n";
    // opserr<<"dVf1 = "<<deltaVf1.Norm()<<"\n";
    opserr<<"predictor  time = "<<timer.getReal()<<"\n";
    timer.start();

    // fluid pressure
    Vector deltaP(Psize);
    if (Psize > 0) {
	double* deltaP_ptr = &deltaP(0);

	// Gft*deltaVf1
        if(Fsize > 0) {
            double* deltaVf1_ptr = &deltaVf1(0);
            cs_gaxpy(Gft, deltaVf1_ptr, deltaP_ptr);
        }

        // Git*deltaVi1
        if(Isize > 0) {
            double* deltaVi1_ptr = &deltaV1(0) + Ssize;
            cs_gaxpy(Git, deltaVi1_ptr, deltaP_ptr);
        }

	// rp-Git*deltaVi1-Gft*deltaVf1
        for(int i=0; i<size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if(rowtype == 3) {             // pressure
		B(i) += deltaP(rowid);
                deltaP(rowid) = B(i);
            }
        }

	// L^{-1}*rp
        Vector x(Psize);
        double* x_ptr = &x(0);
	int* Ap = L->p;
	int* Ai = L->i;
	double* Ax = L->x;
	int status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,x_ptr,deltaP_ptr,LNum,Control,Info);
	deltaP = x;

	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: solving L returns "<<status<<" -- PFEMSolver_Laplace::solve\n";
	    return -1;
	}
    }

    timer.pause();
    // opserr<<"deltaP = "<<deltaP.Norm()<<"\n";
    opserr<<"pressure  time = "<<timer.getReal()<<"\n";

    // Gi, Gf
    cs* Gi = cs_transpose(Git, 1);
    cs* Gf = cs_transpose(Gft, 1);

    // structural and interface velocity corrector
    Vector deltaV(Msize);
    if (Isize > 0 && Psize > 0) {
	double* deltaV_ptr = &deltaV(0);
	double* deltaP_ptr = &deltaP(0);

	// Gi*deltaP
	cs_gaxpy(Gi, deltaP_ptr, deltaV_ptr+Ssize);

	// solve
	Vector x(Msize);
	double* x_ptr = &x(0);
	int* Ap = M->p;
	int* Ai = M->i;
	double* Ax = M->x;
	int status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,x_ptr,deltaV_ptr,MNum,Control,Info);
	deltaV = x;

	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: solving M returns "<<status<<" -- PFEMSolver_Laplace::solve\n";
	    return -1;
	}
    }
    deltaV *= -1.0;
    deltaV += deltaV1;

    // fluid velocity corrector
    Vector deltaVf(Fsize);
    if (Fsize > 0 && Psize > 0) {
	double* deltaVf_ptr = &deltaVf(0);
	double* deltaP_ptr = &deltaP(0);

	// Gf*deltaP
	cs_gaxpy(Gf, deltaP_ptr, deltaVf_ptr);

	// solve
	for(int i=0; i<deltaVf.Size(); i++) {
	    deltaVf(i) /= Mf(i);
        }
    }
    deltaVf *= -1.0;
    deltaVf += deltaVf1;

    timer.pause();
    opserr<<"corrector  time = "<<timer.getReal()<<"\n";
    // timer.start();
    // copy to X
    X.Zero();
    for(int i=0; i<size; i++) {            // row
        int rowtype = dofType(i);          // row type
        int rowid = dofID(i);
        if(rowtype == 0) {
            X(i) = deltaV(rowid);
        } else if(rowtype == 2) {
            X(i) = deltaV(rowid+Ssize);
        } else if(rowtype == 1) {
            X(i) = deltaVf(rowid);
        } else if(rowtype == 3) {
            X(i) = deltaP(rowid);
        }

    }
    // opserr<<"dvi = "<<dvi.Norm()<<"\n";
    // timer.pause();
    // opserr<<"solving time for PFEMSolver_Laplace = "<<timer.getReal()<<"\n";

    cs_spfree(Gi);
    cs_spfree(Gf);

    return 0;
}



int
PFEMSolver_Laplace::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMSolver_Laplace::recvSelf(int ctag,
		  Channel &theChannel,
		  FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int
PFEMSolver_Laplace::setLinearSOE(PFEMLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}
