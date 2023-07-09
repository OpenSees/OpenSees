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
// $Date: 2012-09-17 10:51:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMSolver.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMSolver.h
//
// Written: Minjie 
// Created: Sep 17 2012
//

#include <PFEMSolver.h>
#include <PFEMLinSOE.h>
#include <iostream>
#include <cmath>

void* OPS_PFEMSolver()
{
    PFEMSolver* theSolver = new PFEMSolver();
    return new PFEMLinSOE(*theSolver);
}

PFEMSolver::PFEMSolver()
    :LinearSOESolver(SOLVER_TAGS_PFEMSolver), theSOE(0), Msym(0), Mnum(0)
{
}

PFEMSolver::~PFEMSolver()
{
    if(Msym != 0) {
        cs_sfree(Msym);
    }
    if(Mnum != 0) {
        cs_nfree(Mnum);
    }
}

int
PFEMSolver::solve()
{
    cs* M = theSOE->M;
    cs* Gft = theSOE->Gft;
    cs* Git = theSOE->Git;
    cs* L = theSOE->L;
    cs* Qt = theSOE->Qt;
    Vector& Mhat = theSOE->Mhat;
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
    int Pisize = Mhat.Size();
    int size = X.Size();

    // numeric LU factorization of M
    if(Msize > 0) {
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
            opserr<<"WARNING: failed to do LU factorization of M";
            opserr<<" -- PFEMSolver::solve\n";
            return -1;
        }
    }

    // structure and interface predictor : deltaV1 = M^{-1} * rsi
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

        // M^{-1}*rsi
        Vector x(Msize);
        double* deltaV1_ptr = &deltaV1(0);
        double* x_ptr = &x(0);
        cs_ipvec(Mnum->pinv, deltaV1_ptr, x_ptr, Msize);
        cs_lsolve(Mnum->L, x_ptr);
        cs_usolve(Mnum->U, x_ptr);
        cs_ipvec(Msym->q, x_ptr, deltaV1_ptr, Msize);
    }
    
    // fluid predictor: deltaVf1 = Mf^{-1} * rf
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

    // Mi^{-1}, Msi^{-1}
    cs* invMi1 = cs_spalloc(Isize, Isize, 1, 1, 1);
    cs* invMsi1 = cs_spalloc(Ssize, Isize, 1, 1, 1);
    if(Msize > 0) {
        Vector eyes(Msize);
        Vector x(Msize);
        double* eyes_ptr = &eyes(0);
        double* x_ptr = &x(0);
        for(int j=0; j<Isize; j++) {

            // rhs
            eyes.Zero();
            eyes(j+Ssize) = 1.0;
            x.Zero();

            // M^{-1}*eyes
            cs_ipvec(Mnum->pinv, eyes_ptr, x_ptr, Msize);
            cs_lsolve(Mnum->L, x_ptr);
            cs_usolve(Mnum->U, x_ptr);
            cs_ipvec(Msym->q, x_ptr, eyes_ptr, Msize);

            // copy
            for(int i=0; i<Msize; i++) {
                if(eyes(i) != 0.0) {
                    if(i >= Ssize) {
                        cs_entry(invMi1, i-Ssize, j, eyes(i));
                    } else {
                        cs_entry(invMsi1, i, j, eyes(i));
                    }
                }
            }
        }
    }
    cs* invMi = cs_compress(invMi1);
    cs* invMsi = cs_compress(invMsi1);
    cs_spfree(invMi1);
    cs_spfree(invMsi1);
    if(Mnum != 0) {
        cs_nfree(Mnum);
        Mnum = 0;
    }
    
    // Gi, Mf^{-1}*Gf
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

    // solve for pressure
    Vector deltaP(Psize);
    if(Psize>0) {
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
                deltaP(rowid) = B(i)-deltaP(rowid);   // rp-Git*deltaVi1-Gft*deltaVf1
            }
        }

        // S = L + Git*Mi{-1}*Gi + Gft*Mf{-1}*Gf
        cs* S = 0;
        if(Isize > 0) {
            cs* S1 = cs_multiply(Git, invMi);
            S = cs_multiply(S1, Gi);
            cs_spfree(S1);
        }
        if(Fsize > 0) {
            cs* S1 = cs_multiply(Gft, Gf);
            if(S == 0) {
                S = S1;
            } else {
                cs* S2 = cs_add(S, S1, 1.0, 1.0);
                cs_spfree(S);
                cs_spfree(S1);
                S = S2;
            }
        }
        if(S == 0) {
            S = L;

            // solve
            cs_lusol(3, S, deltaP_ptr, 1e-6);  

        } else {
            cs* S1 = cs_add(S, L, 1.0, 1.0);
            cs_spfree(S);
            S = S1;

            // solve
            cs_lusol(3, S, deltaP_ptr, 1e-6);  
            cs_spfree(S);
        }
    }

    // structure and interface corrector : deltaV = deltaV1 + M^{-1}*G*deltaP
    Vector deltaV(Msize);
    if(Isize > 0) {

        // Gi*deltaP
        Vector Gip(Isize);
        double* Gip_ptr = &Gip(0);
        if(Psize > 0) {
            double* deltaP_ptr = &deltaP(0);
            cs_gaxpy(Gi, deltaP_ptr, Gip_ptr);

            // Msi^{-1}*Gi*deltaP
            if(Ssize > 0) {
                Vector vs(Ssize);
                double* vs_ptr = &vs(0);
                cs_gaxpy(invMsi, Gip_ptr, vs_ptr);
                for(int i=0; i<Ssize; i++) {
                    deltaV(i) += vs(i);
                }
            }

            // Mi^{-1}*Gi*deltaP
            Vector vi(Isize);
            double* vi_ptr = &vi(0);
            cs_gaxpy(invMi, Gip_ptr, vi_ptr);
            for(int i=0; i<Isize; i++) {
                deltaV(i+Ssize) = vi(i);
            }

        }

    }
    deltaV += deltaV1;
    cs_spfree(Gi);
    cs_spfree(invMi);
    cs_spfree(invMsi);

    // fluid corrector: deltaVf = deltaVf1 + Mf^{-1}*Gf*deltaP
    Vector deltaVf(Fsize);
    if(Fsize > 0) {
        if(Psize > 0) {
            double* deltaVf_ptr = &deltaVf(0);
            double* deltaP_ptr = &deltaP(0);
            cs_gaxpy(Gf, deltaP_ptr, deltaVf_ptr);
        }

        deltaVf += deltaVf1;
    }
    cs_spfree(Gf);

    // deltaPi = Mhatd^{-1}*(Qt*deltaP - rpi)
    Vector deltaPi(Pisize);
    if(Pisize > 0) {

        // Qt*deltaP
        if(Psize > 0) {
            double* deltaPi_ptr = &deltaPi(0);
            double* deltaP_ptr = &deltaP(0);
            cs_gaxpy(Qt, deltaP_ptr, deltaPi_ptr);   // Qt*deltaP
        }

        // rpi
        for(int i=0; i<size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if(rowtype != 4) continue;     // rpi
            if(Mhat(rowid) == 0.0) {
                opserr<<"Zero Mhat at location "<<rowid<<" ";
                opserr<<" - PFEMLinSOE::solve()\n";
                return -1;
            }
            deltaPi(rowid) = (B(i)-deltaPi(rowid))/Mhat(rowid); // rpi
        }
    }

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
        } else if(rowtype == 4) {
            X(i) = deltaPi(rowid);
        }

    }

    return 0;
}

int PFEMSolver::setSize()
{
    cs* M = theSOE->M;
    if(M->n > 0) {
        if(Msym != 0) {
            cs_sfree(Msym);
            Msym = 0;
        }
        Msym = cs_sqr(3, M, 0);
        if(Msym == 0) {
            opserr<<"WARNING: failed to do symbolic analysis of M";
            opserr<<" -- PFEMSolver::setSize\n";
            return -1;
        }
    }
    return 0;
}

int
PFEMSolver::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMSolver::recvSelf(int ctag,
		  Channel &theChannel, 
		  FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int 
PFEMSolver::setLinearSOE(PFEMLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}
