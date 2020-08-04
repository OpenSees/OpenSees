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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Umfpack.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Umfpack.h
//
// Written: Minjie 
// Created: Sep 17 2012
//

#include <PFEMSolver_Umfpack.h>
#include <PFEMLinSOE.h>
#include <iostream>
#include <cmath>
// #include <Timer.h>
#include <elementAPI.h>
#include <vector>

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


void* OPS_PFEMSolver_Umfpack()
{
    int numdata = 1;
    int print = 0;
    double ptol = 1e-4;
    int maxiter = 100;

    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* opt = OPS_GetString();
	if (strcmp(opt, "-print") == 0) {

	    print = 1;

	} else if (strcmp(opt, "-ptol") == 0) {

	    if (OPS_GetNumRemainingInputArgs() > 0) {
		if (OPS_GetDoubleInput(&numdata, &ptol) < 0) {
		    opserr << "WARNING: failed to get ptol\n";
		    return 0;
		}
	    }

	    
	} else if (strcmp(opt, "-pmaxiter") == 0) {

	    if (OPS_GetNumRemainingInputArgs() > 0) {
		if (OPS_GetIntInput(&numdata, &maxiter) < 0) {
		    opserr << "WARNING: failed to get max iteration for pressure\n";
		    return 0;
		}
	    }
	}
    }
    
    PFEMSolver_Umfpack* theSolver = new PFEMSolver_Umfpack(ptol,maxiter,print);
    return new PFEMLinSOE(*theSolver);
}

PFEMSolver_Umfpack::PFEMSolver_Umfpack(double tol, int niter, int p)
    :PFEMSolver(), Symbolic(0), theSOE(0), print(p),
     ptol(tol), pmaxiter(niter)
{
}

PFEMSolver_Umfpack::~PFEMSolver_Umfpack()
{
    if (Symbolic != 0) {
	umfpack_di_free_symbolic(&Symbolic);
    }
}

int
PFEMSolver_Umfpack::solve()
{
    // Timer timer;
    // timer.start();
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
    void* Numeric = 0;
    if(Msize > 0) {
	if(Symbolic == 0) {
            opserr<<"WARNING: setSize has not been called";
            opserr<<" -- PFEMSolver_Umfpack::solve\n";
            return -1;
        }
	// numerical analysis
	int* Ap = M->p;
	int* Ai = M->i;
	double* Ax = M->x;
	int status = umfpack_di_numeric(Ap,Ai,Ax,Symbolic,&Numeric,Control,Info);

	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: numeric analysis returns "<<status<<" -- PFEMSolver_Umfpack::solve\n";
	    return -1;
	}

    }
    // timer.pause();
    //opserr<<"factorization time = "<<timer.getReal()<<"\n";
    // timer.start();
    // structure and interface predictor : deltaV1 = M^{-1} * rsi
    Vector deltaV1(Msize);
    if(Msize > 0) {

        // rsi
	Vector rs,ri;
        for(int i=0; i<size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if(rowtype == 2) {
                deltaV1(rowid+Ssize) = B(i);   // rsi
		ri[ri.Size()] = B(i);
            } else if(rowtype == 0) {
                deltaV1(rowid) = B(i);         // rsi
		rs[rs.Size()] = B(i);
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
	int status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,x_ptr,deltaV1_ptr,Numeric,Control,Info);
	deltaV1 = x;
	
	// check error
	if (status!=UMFPACK_OK) {
	    opserr<<"WARNING: solving returns "<<status<<" -- PFEMSolver_Umfpack::solve\n";
	    if (Numeric != 0) {
		umfpack_di_free_numeric(&Numeric);
	    }
	    return -1;
	}	
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

    // timer.pause();
    // opserr<<"dV1 = "<<deltaV1.Norm()<<"\n";
    // opserr<<"dVf1 = "<<deltaVf1.Norm()<<"\n";
    // opserr<<"predictor  time = "<<timer.getReal()<<"\n";
    // timer.start();
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
	    int* Ap = M->p;
	    int* Ai = M->i;
	    double* Ax = M->x;
	    int status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,x_ptr,eyes_ptr,Numeric,Control,Info);
	    eyes = x;

	    // check error
	    if (status!=UMFPACK_OK) {
		opserr<<"WARNING: solving interface returns "<<status<<" -- PFEMSolver_Umfpack::solve\n";
		if (Numeric != 0) {
		    umfpack_di_free_numeric(&Numeric);
		}
	    }
	    
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
    // for (int j=0; j<Isize; j++) {
    // 	opserr<<"col "<<j<<"\n";
    // 	for (int k=invMi->p[j]; k<invMi->p[j+1]; k++) {
    // 	    opserr<<" i = "<<invMi->i[k]<<" ";
    // 	}
    // 	opserr<<"\n";
    // }

    if (Numeric != 0) {
	umfpack_di_free_numeric(&Numeric);
    }
    //opserr<<"invMi = "<<cs_norm(invMi)<<"\n";
    
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
    std::vector<double> deltaP, rhsP;
    if(Psize>0) {
	deltaP.assign(Psize, 0.0);
	rhsP.assign(Psize, 0.0);

        // Gft*deltaVf1
        if(Fsize > 0) {
            double* deltaVf1_ptr = &deltaVf1(0);
            cs_gaxpy(Gft, deltaVf1_ptr, &rhsP[0]);
        }

        // Git*deltaVi1
        if(Isize > 0) {
            double* deltaVi1_ptr = &deltaV1(0) + Ssize;
            cs_gaxpy(Git, deltaVi1_ptr, &rhsP[0]);
        }

        // rp-Git*deltaVi1-Gft*deltaVf1
        for(int i=0; i<size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if(rowtype == 3) {             // pressure
                rhsP[rowid] = B(i)-rhsP[rowid];   // rp-Git*deltaVi1-Gft*deltaVf1
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
            //S = L;
	    opserr<<"WARNING: S=L -- PFEMSolver_Umfpack\n";
	    return -1;

        } 
	cs* S1 = cs_add(S, L, 1.0, 1.0);
	cs_spfree(S);
	S = S1;
	// timer.pause();
	//opserr<<"pressure setup  time = "<<timer.getReal()<<"\n";
	// timer.start();
	if (S->nzmax > 0) {
#ifdef _AMGCL
        try {
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
        // prm.precond.coarse_enough = 1;
        // prm.precond.max_levels = maxlev;
        // prm.precond.direct_coarse = false;
		
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
        } catch (...) {
            opserr << "Pressure: AMGCL solver failed -- Fall back to Umfpack solver\n";
#endif
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

            // symbolic analysis
            void* Ssymbolic = 0;
	        int status = umfpack_di_symbolic(Psize,Psize,Sp,Si,Sx,&Ssymbolic,Control,Info);
	        // check error
	        if (status!=UMFPACK_OK) {
	            opserr<<"WARNING: pressure symbolic analysis returns "<<status<<" -- PFEMSolver_Umfpack::setsize\n";
	            return -1;
	        }

            // numerical analysis
            void* SNumeric = 0;
	        status = umfpack_di_numeric(Sp,Si,Sx,Ssymbolic,&SNumeric,Control,Info);
	        umfpack_di_free_symbolic(&Ssymbolic);

            // check error
	        if (status!=UMFPACK_OK) {
	            opserr<<"WARNING: pressure numeric analysis returns "<<status<<" -- PFEMSolver_Umfpack::solve\n";
	            return -1;
	        }

            // solve
	        std::vector<double> soln(Psize);
	        double* soln_ptr = &soln[0];
            double* rhsP_ptr = &rhsP[0];
	        status = umfpack_di_solve(UMFPACK_A,Sp,Si,Sx,soln_ptr,rhsP_ptr,SNumeric,Control,Info);

            // delete Numeric
	        if (SNumeric != 0) {
	            umfpack_di_free_numeric(&SNumeric);
	        }

            // check error
	        if (status!=UMFPACK_OK) {
	            opserr<<"WARNING: pressure solving returns "<<status<<" -- PFEMSolver_Umfpack::solve\n";
	            return -1;
            }

            deltaP = soln;
#ifdef _AMGCL
        }
#endif
	}

	cs_spfree(S);
    }
    // timer.pause();
    // opserr<<"deltaP = "<<deltaP.Norm()<<"\n";
    // opserr<<"pressure  time = "<<timer.getReal()<<"\n";

    // timer.start();
    // structure and interface corrector : deltaV = deltaV1 + M^{-1}*G*deltaP
    Vector deltaV(Msize);
    if(Isize > 0) {

        // Gi*deltaP
        Vector Gip(Isize);
        double* Gip_ptr = &Gip(0);
        if(Psize > 0) {
            double* deltaP_ptr = &deltaP[0];
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
    //opserr<<"dV = "<<deltaV.Norm()<<"\n";

    // fluid corrector: deltaVf = deltaVf1 + Mf^{-1}*Gf*deltaP
    Vector deltaVf(Fsize);
    if(Fsize > 0) {
        if(Psize > 0) {
            double* deltaVf_ptr = &deltaVf(0);
            double* deltaP_ptr = &deltaP[0];
            cs_gaxpy(Gf, deltaP_ptr, deltaVf_ptr);
        }

        deltaVf += deltaVf1;
    }
    cs_spfree(Gf);
    //opserr<<"dVf = "<<deltaVf.Norm()<<"\n";

    // deltaPi = Mhatd^{-1}*(Qt*deltaP - rpi)
    Vector deltaPi(Pisize);
    if(Pisize > 0) {

        // Qt*deltaP
        if(Psize > 0) {
            double* deltaPi_ptr = &deltaPi(0);
            double* deltaP_ptr = &deltaP[0];
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

    // timer.pause();
    //opserr<<"corrector  time = "<<timer.getReal()<<"\n";
    // timer.start();
    // copy to X
    X.Zero();
    Vector dvi;
    for(int i=0; i<size; i++) {            // row
        int rowtype = dofType(i);          // row type
        int rowid = dofID(i); 
        if(rowtype == 0) {
            X(i) = deltaV(rowid);            
        } else if(rowtype == 2) {
            X(i) = deltaV(rowid+Ssize);
	    dvi[dvi.Size()] = X(i);
        } else if(rowtype == 1) {
            X(i) = deltaVf(rowid);
        } else if(rowtype == 3) {
            X(i) = deltaP[rowid];
        } else if(rowtype == 4) {
            X(i) = deltaPi(rowid);
        }

    }
    // opserr<<"dvi = "<<dvi.Norm()<<"\n";
    // timer.pause();
    // opserr<<"solving time for PFEMSolver_Umfpack = "<<timer.getReal()<<"\n";

    return 0;
}

int PFEMSolver_Umfpack::setSize()
{
    // reorder rows
    cs* M = theSOE->M;
    cs* Gft = theSOE->Gft;
    cs* Git = theSOE->Git;
    cs* L = theSOE->L;
    cs* Qt = theSOE->Qt;
    cs* mats[5] = {M,Gft,Git,L,Qt};
    for (int i=0; i<5; i++) {
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
    
    // Timer timer;
    // timer.start();

    int n = M->n;
    int nnz = M->nzmax;
    if (n == 0 || nnz==0) return 0;

    int* Ap = M->p;
    int* Ai = M->i;
    double* Ax = M->x;

    // symbolic analysis
    if (Symbolic != 0) {
	umfpack_di_free_symbolic(&Symbolic);
    }
    int status = umfpack_di_symbolic(n,n,Ap,Ai,Ax,&Symbolic,Control,Info);

    // check error
    if (status!=UMFPACK_OK) {
	opserr<<"WARNING: symbolic analysis returns "<<status<<" -- PFEMSolver_Umfpack::setsize\n";
	Symbolic = 0;
	return -1;
    }

    // timer.pause();
    //opserr<<"analysis time = "<<timer.getReal()<<"\n";
    return 0;
}

int
PFEMSolver_Umfpack::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMSolver_Umfpack::recvSelf(int ctag,
		  Channel &theChannel, 
		  FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int 
PFEMSolver_Umfpack::setLinearSOE(PFEMLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}
