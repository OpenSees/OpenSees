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

// Written: Minjie
// Created: May 22 2018

#include "PFEMSolver_LumpM.h"
#include "PFEMLinSOE.h"
#include <iostream>
#include <cmath>
#include <Timer.h>
#include <elementAPI.h>

void* OPS_PFEMSolver_LumpM()
{
    bool once = false;
    if (OPS_GetNumRemainingInputArgs() > 0) {
	const char* opt = OPS_GetString();
	if (strcmp(opt, "-once") == 0) {
	    once = true;
	}
    }
    PFEMSolver_LumpM* theSolver = new PFEMSolver_LumpM(once);
    return new PFEMLinSOE(*theSolver);
}

PFEMSolver_LumpM::PFEMSolver_LumpM(bool once)
    :PFEMSolver(), MSym(0), MNum(0), LSym(0), LNum(0), theSOE(0),
     numonce(once), factored(false)
{
}

PFEMSolver_LumpM::~PFEMSolver_LumpM()
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

int PFEMSolver_LumpM::setSize()
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

    // clear
    if (MSym != 0) {
	umfpack_di_free_symbolic(&MSym);
	MSym = 0;
    }

    if (LSym != 0) {
	umfpack_di_free_symbolic(&LSym);
	LSym = 0;
    }

    return 0;
}

int
PFEMSolver_LumpM::solve()
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
    if (LSym != 0) {
	umfpack_di_free_symbolic(&LSym);
	LSym = 0;
    }
    if (!numonce && MNum!=0) {
	umfpack_di_free_numeric(&MNum);
	MNum = 0;
    }

    // check LNum
    if (!numonce && LNum!=0) {
	umfpack_di_free_numeric(&LNum);
	LNum = 0;
    }

    // symbolic analysis of M
    if (MSym==0 && Msize>0) {

	int* Mp = M->p;
	int* Mi = M->i;
	double* Mx = M->x;

	int status = umfpack_di_symbolic(Msize,Msize,Mp,Mi,Mx,&MSym,Control,Info);

	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: symbolic analysis of M returns ";
	    opserr<<status<<" -- PFEMSolver_LumpM::solve\n";
	    return -1;
	}
    }

    // Numerical factorization of M
    if (MNum==0 && Msize>0) {

	// numerical analysis
	int* Mp = M->p;
	int* Mi = M->i;
	double* Mx = M->x;
	int status = umfpack_di_numeric(Mp,Mi,Mx,MSym,&MNum,Control,Info);
	
	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: numeric analysis of M returns ";
	    opserr<<status<<" -- PFEMSolver_LumpM::solve\n";
	    return -1;
	}
    }

    timer.pause();
    opserr<<"factorization time of M = "<<timer.getReal()<<"\n";
    timer.start();
    
    // structure and interface velocity predictor : deltaV1 = M^{-1} * rsi
    Vector deltaV1(Msize);
    if(Msize > 0) {

        // rsi
        for(int i=0; i<size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if(rowtype == 2) {
                deltaV1(rowid+Ssize) = B(i);   // ri
            } else if(rowtype == 0) {
                deltaV1(rowid) = B(i);         // rs
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
	    opserr<<"WARNING: solving M returns "<<status<<" -- PFEMSolver_LumpM::solve\n";
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
    opserr<<"predictor time = "<<timer.getReal()<<"\n";
    timer.start();

    // Mi
    Vector Mi(Isize);
    for(int j=Ssize; j<Msize; j++) {
	for(int k=M->p[j]; k<M->p[j+1]; k++) {
	    Mi(j-Ssize) += M->x[k];
	}
    }

    // Gi, Gf
    cs* Gi = cs_transpose(Git, 1);
    cs* Gf = cs_transpose(Gft, 1);
    if(Fsize > 0) {
        for(int j=0; j<Psize; j++) {
            for(int k=Gf->p[j]; k<Gf->p[j+1]; k++) {
                int i = Gf->i[k];
                double& x = Gf->x[k];
                x /= Mf(i);
            }
        }
    }
    if(Isize > 0) {
        for(int j=0; j<Psize; j++) {
            for(int k=Gi->p[j]; k<Gi->p[j+1]; k++) {
                int i = Gi->i[k];
                double& x = Gi->x[k];
                x /= Mi(i);
            }
        }
    }

    timer.pause();
    opserr<<"Mi,Gi,Gf time = "<<timer.getReal()<<"\n";
    timer.start();

    // fluid pressure
    Vector deltaP(Psize);
    if (Psize > 0) {
	double* deltaP_ptr = &deltaP(0);

	// GfT*Mf^{-1}*Gf+GiT*Mi^{-1}*Gi
	cs* S = 0;
	if (Fsize > 0) {
	    cs* Lf = cs_multiply(Gft,Gf);
	    S = cs_add(Lf,L,1.0,1.0);
	    cs_spfree(Lf);
	}
	if (Isize > 0) {
	    cs* Li = cs_multiply(Git,Gi);
	    if (S == 0) {
		S = cs_add(Li,L,1.0,1.0);
	    } else {
		cs* temp2 = cs_add(Li,S,1.0,1.0);
		cs_spfree(S);
		S = temp2;
	    }
	    cs_spfree(Li);
	}

	if (S == 0) {
	    opserr<<"WARNING: S=L -- PFEMSolver_LumpM::solve\n";
	    return -1;
	}

	int* Sp = S->p;
	int* Si = S->i;
	double* Sx = S->x;
	for (int j=0; j<Psize; j++) {
	    ID col(0, Sp[j+1]-Sp[j]);
	    Vector colval(Sp[j+1]-Sp[j]);
	    ID col0(colval.Size());
	    int index = 0;
	    for (int k=Sp[j]; k<Sp[j+1]; k++) {
		col.insert(Si[k]);
		col0(index) = Si[k];
		colval(index++) = Sx[k];
	    }
	    index = 0;
	    for (int k=Sp[j]; k<Sp[j+1]; k++) {
		Si[k] = col[index++];
		Sx[k] = colval(col0.getLocation(Si[k]));
	    }
	}

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
                deltaP(rowid) = B(i) - deltaP(rowid);
            }
        }

	timer.pause();
	opserr<<"matrix time for L = "<<timer.getReal()<<"\n";
	timer.start();

	// symbolic analysis of L
	if (LSym == 0) {

	    int status = umfpack_di_symbolic(Psize,Psize,Sp,Si,Sx,&LSym,Control,Info);

	    // check error
	    if (status!=UMFPACK_OK) {
		opserr<<"WARNING: symbolic analysis of L returns ";
		opserr<<status<<" -- PFEMSolver_LumpM::solve\n";
		opserr<<"UMFPACK_ERROR_n_nonpositive = "<<UMFPACK_ERROR_n_nonpositive<<"\n";
		opserr<<"UMFPACK_ERROR_invalid_matrix = "<<UMFPACK_ERROR_invalid_matrix<<"\n";
		opserr<<"UMFPACK_ERROR_out_of_memory = "<<UMFPACK_ERROR_out_of_memory<<"\n";
		opserr<<"UMFPACK_ERROR_argument_missing = "<<UMFPACK_ERROR_argument_missing<<"\n";
		opserr<<"UMFPACK_ERROR_internal_error = "<<UMFPACK_ERROR_internal_error<<"\n";
		return -1;
	    }
	}

	timer.pause();
	opserr<<"Symbolic time for L = "<<timer.getReal()<<"\n";
	timer.start();
	
	
	// numerical analysis of L
	if (LNum == 0) {
	    int status = umfpack_di_numeric(Sp,Si,Sx,LSym,&LNum,Control,Info);
	
	    // check error
	    if (status!=UMFPACK_OK) {
		opserr<<"WARNING: numeric analysis of L returns ";
		opserr<<status<<" -- PFEMSolver_LumpM::solve\n";
		return -1;
	    }
	}

	timer.pause();
	opserr<<"Numerical time for L = "<<timer.getReal()<<"\n";
	timer.start();

	// L^{-1}*rp
        Vector x(Psize);
        double* x_ptr = &x(0);
	int* Ap = S->p;
	int* Ai = S->i;
	double* Ax = S->x;
	int status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,x_ptr,deltaP_ptr,LNum,Control,Info);
	deltaP = x;

	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: solving L returns "<<status<<" -- PFEMSolver_LumpM::solve\n";
	    return -1;
	}

	cs_spfree(S);
    }

    timer.pause();
    // opserr<<"deltaP = "<<deltaP.Norm()<<"\n";
    opserr<<"pressure  time = "<<timer.getReal()<<"\n";

    

    // structural and interface velocity corrector
    cs_spfree(Gi);
    Gi = cs_transpose(Git, 1);
    
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
	    opserr<<"WARNING: solving M returns "<<status<<" -- PFEMSolver_LumpM::solve\n";
	    return -1;
	}
    }
    deltaV += deltaV1;

    // fluid velocity corrector
    Vector deltaVf = deltaVf1;
    if (Fsize > 0 && Psize > 0) {
	double* deltaVf_ptr = &deltaVf(0);
	double* deltaP_ptr = &deltaP(0);

	// Vf1+Gf*deltaP
	cs_gaxpy(Gf, deltaP_ptr, deltaVf_ptr);
    }

    timer.pause();
    opserr<<"corrector  time = "<<timer.getReal()<<"\n";
    timer.start();
    
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
    // opserr<<"solving time for PFEMSolver_LumpM = "<<timer.getReal()<<"\n";

    cs_spfree(Gi);
    cs_spfree(Gf);

    return 0;
}



int
PFEMSolver_LumpM::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMSolver_LumpM::recvSelf(int ctag,
		  Channel &theChannel,
		  FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int
PFEMSolver_LumpM::setLinearSOE(PFEMLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}
