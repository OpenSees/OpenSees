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


// Written: Minjie 

#include <PFEMSolver_Mumps.h>
#include <PFEMLinSOE.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <Timer.h>
#ifdef _MUMPS
#include <mpi.h>
#endif
#include <elementAPI.h>

#ifdef _AMGCL
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/solver/fgmres.hpp>
#include <amgcl/solver/lgmres.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#endif

void* OPS_PFEMSolver_Mumps()
{
    int numdata = 1;
    int relax = 20;
    int err = 0;
    int sym = 0;
    int add = 0;
    int print = 0;
    double ptol = 1e-4;
    double Bitol = 1e-16;
    int maxiter = 100;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();
        if (strcmp(opt, "-relax") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numdata, &relax) < 0) {
                    opserr << "WARNING: failed to get relax\n";
                    return 0;
                }
            }

        } else if (strcmp(opt, "-err") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numdata, &err) < 0) {
                    opserr << "WARNING: failed to get err\n";
                    return 0;
                }
            }

        } else if (strcmp(opt, "-sym") == 0) {

            sym = 1;

        } else if (strcmp(opt, "-print") == 0) {

            print = 1;

        } else if (strcmp(opt, "-added-mass") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numdata, &add) < 0) {
                    opserr << "WARNING: failed to get add\n";
                    return 0;
                }
            }

        } else if (strcmp(opt, "-ptol") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &ptol) < 0) {
                    opserr << "WARNING: failed to get ptol\n";
                    return 0;
                }
            }

        } else if (strcmp(opt, "-Bitol") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &Bitol) < 0) {
                    opserr << "WARNING: failed to get Bitol\n";
                    return 0;
                }
            }

        } else if (strcmp(opt, "-pmaxiter") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numdata, &maxiter) < 0) {
                    opserr << "WARNING: failed to get err\n";
                    return 0;
                }
            }
        }
    }

    PFEMSolver_Mumps* theSolver = new PFEMSolver_Mumps(relax,err,add,sym,print,
                                                       ptol,maxiter,Bitol);
    return new PFEMLinSOE(*theSolver);
}

PFEMSolver_Mumps::PFEMSolver_Mumps(int r, int e, int a, int s, int p,
                                   double tol, int niter, double bitol)
        :PFEMSolver(), theSOE(0),
#ifdef _MUMPS
         sid(),
#endif
         relax(r), err(e), add(a), sym(s), print(p),
         ptol(tol), Bitol(bitol), pmaxiter(niter)
{
#ifdef _MUMPS
    // mumps id
    sid.job = JOB_INIT;
    sid.par = 1;
    sid.sym = sym;
    sid.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&sid);

    sid.irn = 0;
    sid.jcn = 0;
#endif
}

PFEMSolver_Mumps::~PFEMSolver_Mumps()
{
#ifdef _MUMPS
    sid.job = JOB_END;
    dmumps_c(&sid);

    if(sid.irn != 0) delete [] sid.irn;
    if(sid.jcn != 0) delete [] sid.jcn;
#endif
}

int
PFEMSolver_Mumps::solve()
{
#ifdef _MUMPS
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

    // numeric LU factorization of M
    if(Msize > 0) {
        sid.job = JOB_FACTORIZATION;
        dmumps_c(&sid);

        if(sid.info[0] != 0) {
            opserr<<"WARNING: failed to factorize -- PFEMSolver_Mumps::solve\n";
            return -1;
        }
    }

    // structure and interface predictor : deltaV1 = M^{-1} * rsi
    std::vector<double> deltaV1;
    if(Msize > 0) {
        deltaV1.assign(Msize, 0.0);
        // rsi
        for (int i = 0; i < size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if (rowtype == 2) {
                deltaV1[rowid + Ssize] = B(i);   // rsi
            } else if (rowtype == 0) {
                deltaV1[rowid] = B(i);         // rsi
            }
        }
        sid.rhs = &deltaV1[0];
        sid.nrhs = 1;
        sid.job = JOB_SOLUTION;
        dmumps_c(&sid);
        if (sid.info[0] != 0) {
            opserr << "WARNING: failed to solve predictor -- PFEMSolver_Mumps::solve\n";
            return -1;
        }
    }

    // fluid predictor: deltaVf1 = Mf^{-1} * rf
    std::vector<double> deltaVf1;
    if (Fsize > 0) {
        deltaVf1.assign(Fsize,0.0);

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
                deltaVf1[rowid] = B(i)/Mf(rowid);         // rf
            }
        }
    }

    // Gi, Gf
    cs* Gi = 0;
    cs* Gf = 0;
    if(Fsize>0) {
        Gf = cs_transpose(Gft, 1);
    }
    if(Isize>0) {
        Gi = cs_transpose(Git, 1);
    }

    // solve for pressure
    std::vector<double> deltaP, rhsP;
    if(Psize>0) {
        deltaP.assign(Psize, 0.0);
        rhsP.assign(Psize, 0.0);
    }
    if(Psize>0) {

        // S = L + Git*Mi{-1}*Gi + Gft*Mf{-1}*Gf
        cs* S = 0;

        // Gft*deltaVf1
        if(Fsize > 0) {
            cs_gaxpy(Gft, &deltaVf1[0], &rhsP[0]);
        }

        // Git*deltaVi1
        if(Isize > 0) {
            cs_gaxpy(Git, &deltaV1[0]+Ssize, &rhsP[0]);
        }

        // rp-Git*deltaVi1-Gft*deltaVf1
        for(int i=0; i<size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if(rowtype == 3) {             // pressure
                rhsP[rowid] = B(i)-rhsP[rowid];   // rp-Git*deltaVi1-Gft*deltaVf1
            }
        }

        // Git*Mi{-1}*Gi
        if (Isize > 0) {
            std::vector<int> irhs_ptr(Msize+1,1);
            std::vector<int> irhs_row(Isize*Isize);
            std::vector<double> rhs_val(Isize*Isize);

            for (int j=0; j<Isize; ++j) {
                for (int i=0; i<Isize; ++i) {
                    irhs_row[Isize*j+i] = i+Ssize+1;
                }
                irhs_ptr[j+1+Ssize] = (j+1)*Isize+1;
            }

            ICNTL(sid,20,3);
            ICNTL(sid,30,1);
            sid.nz_rhs = (int)rhs_val.size();
            sid.nrhs = Msize;
            sid.irhs_ptr = &irhs_ptr[0];
            sid.irhs_sparse = &irhs_row[0];
            sid.rhs_sparse = &rhs_val[0];

            sid.job = JOB_SOLUTION;
            dmumps_c(&sid);
            if(sid.info[0] != 0) {
                opserr<<"info[1] = "<<sid.info[0]<<", info[2] ="<<sid.info[1]<<"\n";
                opserr<<"WARNING: failed to solve Bi -- PFEMSolver_Mumps::solve\n";
                return -1;
            }
            ICNTL(sid,20,0);
            ICNTL(sid,30,0);

            for (int j=Isize; j<Msize; ++j) {
                irhs_ptr[j+1] -= 1;
            }
            for (int k=0; k<(int)irhs_row.size(); ++k) {
                irhs_row[k] -= 1;
                irhs_row[k] -= Ssize;
            }

            cs* Bi1 = cs_spalloc(Isize,Isize,1,1,1);
            int numignore = 0;
            for (int j=0; j<Isize; ++j) {
                for (int i=0; i<Isize; ++i) {
                    if (fabs(rhs_val[Isize*j+i]) > Bitol) {
                        cs_entry(Bi1, i, j,rhs_val[Isize*j+i]);
                    } else {
                        numignore++;
                    }
                }
            }
            cs* Bi = cs_compress(Bi1);
            cs_spfree(Bi1);

            cs* S1 = cs_multiply(Git, Bi);
            S = cs_multiply(S1, Gi);
            cs_spfree(S1);
            cs_spfree(Bi);
        }

        // Gft*Mf{-1}*Gf
        if(Fsize > 0) {
            for(int j=0; j<Psize; j++) {
                for(int k=Gf->p[j]; k<Gf->p[j+1]; k++) {
                    Gf->x[k] /= Mf(Gf->i[k]);
                }
            }
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

        // S
        if(S == 0) {
            S = L;
        } else {
            cs* S1 = cs_add(S, L, 1.0, 1.0);
            cs_spfree(S);
            S = S1;
        }

        // solve
        if (S->nzmax > 0) {
#ifdef _AMGCL
            // solve
            amgcl::profiler<> prof;
            typedef
            amgcl::make_solver<
                    amgcl::amg<
                            amgcl::backend::builtin<double>,
                            amgcl::coarsening::smoothed_aggregation,
                            amgcl::relaxation::spai0
                    >,
                    amgcl::solver::lgmres<amgcl::backend::builtin<double> >
            > Solver;


            // parameter
            Solver::params prm;
            prm.solver.tol = ptol;
            prm.solver.maxiter = pmaxiter;

            // setup
            prof.tic("setup");
            std::vector<std::ptrdiff_t> ptr(Psize+1), num(S->nzmax);
            for (int i=0; i<S->nzmax; ++i) {
                num[i] = S->i[i];
            }
            for (int i = 0; i < Psize+1; ++i) {
                ptr[i] = S->p[i];
            }
            double* val = &(S->x[0]);
            Solver solve(amgcl::adapter::zero_copy(Psize,&ptr[0],&num[0],val),prm);
            prof.toc("setup");

            // solve
            int iters;
            double error;
            prof.tic("solve");
            std::tie(iters, error) = solve(rhsP, deltaP);
            prof.toc("solve");

            if (print) {
                std::cout << solve << std::endl;
                std::cout << "iters: " << iters << std::endl
                          << "error: " << error << std::endl
                          << prof << std::endl;
            }

            if (iters>=pmaxiter && error>ptol) {
                opserr<<"WARNING: failed to solve pressure\n";
                return -1;
            }
#endif
        }

        // release
        if(S != L) cs_spfree(S);
    }

    // structure and interface corrector : deltaV = deltaV1 + M^{-1}*G*deltaP
    std::vector<double> deltaV;
    if (Msize > 0) {
        deltaV.assign(Msize, 0.0);
    }
    if(Isize > 0) {
        // Gi*deltaP
        if(Psize > 0) {
            cs_gaxpy(Gi, &deltaP[0], &deltaV[0]+Ssize);
        }

        sid.rhs = &deltaV[0];
        sid.nrhs = 1;
        sid.job = JOB_SOLUTION;
        dmumps_c(&sid);
        if(sid.info[0] != 0) {
            opserr<<"WARNING: failed to solve corrector -- PFEMSolver_Mumps::solve\n";
            return -1;
        }
    }
    for (int i = 0; i < Msize; ++i) {
        deltaV[i] += deltaV1[i];
    }

    // fluid corrector: deltaVf = deltaVf1 + Mf^{-1}*Gf*deltaP
    std::vector<double> deltaVf;
    if(Fsize > 0) {
        deltaVf.assign(Fsize, 0.0);
    }
    if(Fsize > 0) {
        if(Psize > 0) {
            cs_gaxpy(Gf, &deltaP[0], &deltaVf[0]);
        }
    }
    for (int i=0; i<Fsize; ++i) {
        deltaVf[i] = deltaVf[i] + deltaVf1[i];
    }

    // delete Gi, Gf
    if(Fsize>0) {
        cs_spfree(Gf);
    }
    if(Isize>0) {
        cs_spfree(Gi);
    }

    // copy to X
    X.Zero();

    for(int i=0; i<size; i++) {            // row
        int rowtype = dofType(i);          // row type
        int rowid = dofID(i);
        if(rowtype == 0) {
            X(i) = deltaV[rowid];
        } else if(rowtype == 2) {
            X(i) = deltaV[rowid+Ssize];
        } else if(rowtype == 1) {
            X(i) = deltaVf[rowid];
        } else if(rowtype == 3) {
            X(i) = deltaP[rowid];
        }
    }
#endif
    return 0;
}

int PFEMSolver_Mumps::setSize()
{
#ifdef _MUMPS
    // assembled format
    ICNTL(sid,5,0);

    // input matrix is centralized on the host
    ICNTL(sid,18,0);

    // workspace relaxation: 20%
    if(relax <= 0) {
        relax = 20;
    }
    ICNTL(sid,14,relax);

    // dense right hand side
    ICNTL(sid,20,0);

    // centralized right hand side
    ICNTL(sid,21,0);

    // No error messages
    if(err < 0) err = 0;
    ICNTL(sid,1,err);
    ICNTL(sid,2,err);
    ICNTL(sid,3,err);
    ICNTL(sid,4,err);

    // check M size
    cs* M = theSOE->M;
    int Msize = M->n;
    if(Msize <= 0) return 0;

    // analysis
    sid.n = M->n;
    sid.nz = M->nzmax;
    sid.a = M->x;
    if(sid.irn != 0) delete [] sid.irn;
    if(sid.jcn != 0) delete [] sid.jcn;
    sid.irn = new int[sid.nz];
    sid.jcn = new int[sid.nz];

    for(int j=0; j<sid.n; j++) {
        for(int k=M->p[j]; k<M->p[j+1]; k++) {
            sid.irn[k] = M->i[k]+1;
            sid.jcn[k] = j+1;
        }
    }

    // call mumps
    sid.job = JOB_ANALYSIS;
    dmumps_c(&sid);
    if(sid.info[0] != 0) {
        opserr<<"WARNING: failed to analyze -- PFEMSolver_Mumps::setSize\n";
        return -1;
    }
#endif
    return 0;
}

int
PFEMSolver_Mumps::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMSolver_Mumps::recvSelf(int ctag,
                           Channel &theChannel,
                           FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int
PFEMSolver_Mumps::setLinearSOE(PFEMLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}

#ifdef _MUMPS
void
PFEMSolver_Mumps::ICNTL(DMUMPS_STRUC_C& id, int I, int val)
{
    id.icntl[I-1] = val;
}
#endif
